#include "aclass.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <ctime>
#include <chrono>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/CholmodSupport>
#include <gurobi_c++.h>
#include "mymesh/readers.h"
typedef std::chrono::high_resolution_clock Clock;


bool Graphs::readGraphsFromCrossections(CrossSections &CS){



    vpos = CS.V;
    nV = vpos.size()/3;



    plane2vers = CS.planeV2globalV;
    nP = plane2vers.size();
    ver2planes.clear();
    ver2planes.resize(nV);
    for(int i=0;i<nP;i++){
        for(auto a:plane2vers[i])ver2planes[a].push_back(i);
    }


    isectIndices = CS.isectIndices;


    lines = CS.lineIndices;
    lines2Edges = CS.line2Edges;
    nL = lines2Edges.size();
    for(int i=0;i<nL;i++){
        int ne = lines2Edges[i].size()/2;
        for(int j=0;j<ne*2;++j){
            linesEdges.push_back(lines2Edges[i][j]);
        }
    }



    lineVertices2nv.clear();
    lineVertices2nv.resize(nV,vector<int>(0));



    for(int j=0,k=linesEdges.size()/2;j<k;++j){
        int v1 = linesEdges[j*2],v2=linesEdges[j*2+1];
        lineVertices2nv[v1].push_back(v2);
        lineVertices2nv[v2].push_back(v1);
    }



    plane2lines = CS.plane2lines;


    planeinfo.resize(nP);
    for(int i=0;i<nP;++i)planeinfo[i] = CS.planes_info[i].planepara;


    //readValues();
    oldval = CS.val;
    nMat = CS.n_materials;


    mapplaneV2globalV = CS.planeV2globalV;



    nodeMarkers = CS.nodeMarkers;

    E = CS.newE;



    vector<uint>f2v;
    vector<int>fmarker;
    for(int i=0;i<nP;i++){
        auto mapplaneV2globalV_local = mapplaneV2globalV[i];
        int nt = CS.T[i].size()/3;
        for(auto &a:CS.T[i]){
            f2v.push_back(mapplaneV2globalV_local[a]);
        }
        for(int j=0;j<nt;++j)fmarker.push_back(i);
    }


    map<vector<int>,int>edgemarkers;
    for(int i=0;i<nL;++i){
        int ne = lines2Edges[i].size()/2;

        for(int j=0;j<ne;++j){
            vector<int>eee;
            eee.push_back(lines2Edges[i][j*2]);eee.push_back(lines2Edges[i][j*2+1]);
            sort(eee.begin(),eee.end());
            edgemarkers[eee] = NONMANIFOLDBEGIN + i;
        }

    }

    nmGraph.init(vpos,f2v,fmarker,edgemarkers);





    return true;


}
int Graphs::labelonPlane(vector<vector<double>>&val,int iPnt,int iplane){

    auto p_dbe = val[iPnt].begin() + iplane*nMat;
    auto p_ded = val[iPnt].begin() + (1+iplane)*nMat;

    return max_element(p_dbe,p_ded)-p_dbe;


}
int maxInd(double *pval,int n){

    int ind = -1;
    double maxv = -INFTY_DOUBLE;
    for(int i=0;i<n;++i)if(pval[i]>maxv){
        maxv = pval[i];
        ind = i;
    }
    return ind;

}
int Graphs::averageAcrossPlanes(vector<vector<double>>&val,int iPnt, vector<double>&ave){


    ave.clear();
    ave.resize(nMat,0);

    auto p_d = val[iPnt].data();

    for(int i=0;i<ver2planes[iPnt].size();++i){
        int iP = ver2planes[iPnt][i];
        for(int j=0;j<nMat;++j)ave[j] += p_d[iP*nMat+j];

    }

    for(auto &a:ave)a/=ver2planes[iPnt].size();
    return max_element(ave.begin(),ave.end())-ave.begin();


}

void Graphs::solveLine(){


    int totaleffectiveLine = 0,totalLine = 0;
    vector<bool>isIsectWithotherLines(nV,false);
    for(auto a:isectIndices)isIsectWithotherLines[a] = true;

    for(int iplane = 0;iplane < nP;++iplane){


        for(int i = 0;i<plane2lines[iplane].size();++i){

            bool ischanged = false;

            int iline = plane2lines[iplane][i];

            vector<int >vind(nV,-1);

            vector<nodeinfo>nodes;
            vector<int>ischangerecord;

            for(int j = 0;j<lines[iline].size();++j){

                vector<double>ave;

                int iPnt = lines[iline][j];
                int avelabel = averageAcrossPlanes(oldval,iPnt,ave);

                if(isIsectWithotherLines[iPnt]>0){
                    //nodes.push_back(nodeinfo(oldval[iPnt].data()+iplane*nMat,nMat,avelabel,true));
                    if(labelonPlane(oldval,iPnt,iplane)!=avelabel)ischanged = true;
                    ischangerecord.push_back(avelabel);
                    nodes.push_back(nodeinfo(oldval[iPnt].data()+iplane*nMat,iPnt,avelabel,true));
                }else{
                    nodes.push_back(nodeinfo(oldval[iPnt].data()+iplane*nMat,iPnt,avelabel,false));
                    //nodes.push_back(nodeinfo(oldval[iPnt].data()+iplane*nMat,nMat,avelabel,false));
                    ischangerecord.push_back(9);
                }
                vind[iPnt] = j;
            }

            ++totalLine;
            if(!ischanged)continue;
            ++totaleffectiveLine;
            cout<<"iline: "<<iline<<endl;


            //vector<vector<double> >graphconnection(nodes.size(),vector<double>(nodes.size(),0));
            vector< triple >graphconnection;
            vector<double>central(nodes.size(),0);

            for(int j=0;j<lines2Edges[iline].size()/2;++j){

                auto p_ev = lines2Edges[iline].data()+j*2;
                int v0 = vind[p_ev[0]];
                int v1 = vind[p_ev[1]];

                double edge_val = -1;//use cottan weigh here

                central[v0] -= edge_val;central[v1] -= edge_val;

                graphconnection.push_back(triple(v0,v1,edge_val));
                graphconnection.push_back(triple(v1,v0,edge_val));

                //graphconnection[v0][v1] = graphconnection[v1][v0] = 1;

            }
            for(int j=0;j<graphconnection.size();++j){
                graphconnection[j].val/=sqrt(central[graphconnection[j].i]*central[graphconnection[j].j]);
                //graphconnection[j].val/=central[graphconnection[j].i];
            }

            for(int j=0;j<nodes.size();++j){
                graphconnection.push_back(triple(j,j,1));
            }

            vector<double>solveval;

            buildMatrixandSolveLine(nodes,graphconnection,solveval);


            for(int k=0;k<nodes.size();++k){
                double *pval = newval[nodes[k].iPnt].data()+iplane*nMat;
                for(int j=0;j<nMat;++j){
                    pval[j] = solveval[k*nMat+j];
                }
            }

            if(0){
                debugline(iline,iplane);
                for(int j = 0;j<lines[iline].size();++j){
                    int iPnt = lines[iline][j];
                    cout<<iPnt<<' ';
                }cout<<endl;
                for(int j = 0;j<lines[iline].size();++j){
                    cout<<ischangerecord[j]<<' ';
                }
                cout<<endl;
                for(int j = 0;j<lines[iline].size();++j){
                    int iPnt = lines[iline][j];
                    cout<<maxInd(oldval[iPnt].data()+iplane*nMat,nMat)<<' ';
                }
                cout<<endl;
                for(int j = 0;j<lines[iline].size();++j){
                    int iPnt = lines[iline][j];
                    cout<<maxInd(newval[iPnt].data()+iplane*nMat,nMat)<<' ';
                }
                cout<<endl;
                exit(0);
            }

        }
    }

    cout<<"totalLine: "<<totalLine<<' '<<totaleffectiveLine<<endl;



}

void Graphs::debugline(int iline,int iplane){

    string p1("put an address");

    vector<double>vposition;
    vector<int>ind(nV,-1);
    auto p_vstart = vpos.data()+lines[iline][0]*3;
    vector<double>kpos(3,0);
    for(int j = 0;j<lines[iline].size();++j){
        int iPnt = lines[iline][j];
        auto p_v = vpos.data()+iPnt*3;
        kpos[0] = _VerticesDistance(p_v,p_vstart);
        for(int i=0;i<3;++i)vposition.push_back(kpos[i]);
        ind[iPnt] = j;
    }
    vector<uint>edges;

    for(auto a:lines2Edges[iline]){
        edges.push_back(ind[a]);
    }

    for(int i=0;i<3;++i){
        vposition.clear();
        for(int j = 0;j<lines[iline].size();++j){
            int iPnt = lines[iline][j];
            auto p_v = vpos.data()+iPnt*3;
            kpos[0] = _VerticesDistance(p_v,p_vstart);
            kpos[2] = oldval[iPnt][iplane*nMat+i];
            for(int i=0;i<3;++i)vposition.push_back(kpos[i]);

        }

        writeObjFile_line(p1+to_string(i),vposition,edges);

        vposition.clear();
        for(int j = 0;j<lines[iline].size();++j){
            int iPnt = lines[iline][j];
            auto p_v = vpos.data()+iPnt*3;
            kpos[0] = _VerticesDistance(p_v,p_vstart);
            kpos[2] = newval[iPnt][iplane*nMat+i];
            for(int i=0;i<3;++i)vposition.push_back(kpos[i]);

        }

        writeObjFile_line(p1+to_string(10+i),vposition,edges);
    }

}

void Graphs::debugline2(vector<vector<double> >&val){

    int n_line = lines.size();
    vector<vector<int>>line2plane(n_line);

    for(int j=0;j<plane2lines.size();++j){
        for(int i=0;i<plane2lines[j].size();++i){
            line2plane[plane2lines[j][i]].push_back(j);
        }
    }

    cout<<"debugline2"<<endl;


    for(int i = 0;i<n_line;++i){
        if(i!=1)continue;
        cout<<"iline: "<<i<<endl;
        for(int j=0; j<line2plane[i].size();++j){
            int iplane = line2plane[i][j];
            for(auto ipnt:lines[i]){
                cout<<maxInd(val[ipnt].data()+iplane*nMat,nMat)<<' ';
            }cout<<endl;
        }
        for(auto ipnt:lines[i]){
            auto &pval = newval[ipnt];
            vector<double>ave(nMat,0);
            for(int j=0;j<nP;++j){
                auto pv = pval.begin()+nMat*j;
                bool isNull = true;
                for(int k=0;k<nMat;++k)isNull  = isNull&&(pv[k]==-1);
                if(isNull)continue;
                for(int k=0;k<nMat;++k)ave[k]+=pv[k];
            }
            int label = max_element(ave.begin(),ave.end())-ave.begin();
            if(label ==6)
                cout<<ipnt<<endl;
            cout<< label <<' ';
        }cout<<endl;

    }


    exit(9);





}


void Subset(vector<int>&inputs,int s_subset,vector<vector<int> >& outputs){


    std::function<void(int *,int,int,int,vector<int>&)>  func;
    func = [&outputs,&func](int *parr,int size,int left, int index,vector<int>&container){

        if(left==0){
            outputs.push_back(container);
            return;
        }
        for(int i=index;i<size;++i){
            container.push_back(parr[i]);
            func(parr,size,left-1,i+1,container);
            container.pop_back();
        }

    };

    outputs.clear();

    vector<int>container;

    func(inputs.data(),inputs.size(),s_subset,0,container);



}


void Subset(int n,int s_subset,vector<vector<int> >& outputs){
    vector<int>mats(n);for(int i=0;i<n;++i)mats[i] = i;
    Subset(mats,s_subset,outputs);
}



typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> SpMat;
void Graphs::buildMatrixandSolveLine(vector<nodeinfo>&nodes, vector<triple> &graphconnection_triple, vector<double>&solveval){


    double alpha = 1e-3;
    int n_node = nodes.size();
    int n_var = n_node*nMat;


    vector<triple>cAb;

    vector<triple>Mmatrix;



    int nc = 0;

    for(int i=0;i<n_node;++i)if(nodes[i].isConstrained){
        double *pval = nodes[i].pval;
        int iPnt = nodes[i].iPnt;
        int label = nodes[i].label;
        for(int j=0;j<nMat;++j)if(j!=label){
            cAb.push_back(triple(i*nMat+label,i*nMat+j,pval[j]-pval[label]+10*THRES));
        }
        ++nc;
    }

    if(nc==0){
        solveval.clear();
        solveval.resize(n_var,0.);
    }else{


        vector<T>Ts;
        for(int i=0;i<graphconnection_triple.size();++i){
            auto &pp = graphconnection_triple[i];
            Ts.push_back(T(pp.i,pp.j,pp.val));

        }

        SpMat L(n_node,n_node);
        L.setFromTriplets(Ts.begin(),Ts.end());
		L = L*L;

        SpMat M(n_var,n_var);
        M.setZero();


        auto constructM = [&M,&L,n_node,n_var,this](int m1,int m2){
            vector<T>Ts;
            SpMat A(n_node,n_var);
            for(int i=0;i<n_node;++i){
                Ts.push_back(T(i,m1+i*nMat,1));
                Ts.push_back(T(i,m2+i*nMat,-1));
            }
            A.setFromTriplets(Ts.begin(),Ts.end());
            M += A.transpose() * L * A;
        };


        vector<vector<int> >matset;

        Subset(nMat,2,matset);

        for(auto &a:matset)constructM(a[0],a[1]);



        Ts.clear();
        for(int i=0;i<n_var;++i)Ts.push_back(T(i,i,alpha));
        SpMat Eye(n_var,n_var);
        Eye.setFromTriplets(Ts.begin(),Ts.end());

        M += Eye;

        cout<<"M.size(): "<<M.nonZeros()<<endl;

        for (int k=0; k<M.outerSize(); ++k)
            for (SpMat::InnerIterator it(M,k); it; ++it)
            {
                Mmatrix.push_back(triple(it.row(),it.col(),it.value()));
            }

        solveQuadraticProgramming(Mmatrix,cAb,n_var,solveval);


    }


    for(int i=0;i<n_node;++i){
        double *pval = nodes[i].pval;
        for(int j=0;j<nMat;++j){
            solveval[i*nMat+j]+=pval[j];
        }
    }



}

void Graphs::solveQuadraticProgramming(vector<triple> &M, vector<triple> &Ab,int n, vector<double>&solveval){

    vector<GRBVar>vars(n);

    try{

        GRBEnv env = GRBEnv();
        //cout << "optimize 0" << endl;
        GRBModel model = GRBModel(env);

        model.getEnv().set(GRB_IntParam_OutputFlag, 0);

        conVars.resize(nP);


        for(int i=0;i<n;++i)vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,to_string(i));




        //cout << "optimize 1" << endl;
        GRBQuadExpr obj = 0;

        for(int i=0;i<M.size();++i){

            obj += vars[M[i].i] * vars[M[i].j] * M[i].val;


        }

        // cout << "optimize 2 " <<codeer<< endl;
        model.setObjective(obj, GRB_MINIMIZE);


        //model.update();

        for(int i=0;i<Ab.size();++i){
            model.addConstr(vars[Ab[i].i] - vars[Ab[i].j] >= Ab[i].val);
        }



        auto t1 = Clock::now();

        model.optimize();

        auto t2 = Clock::now();

        cout << "Total opt time: " << std::chrono::nanoseconds(t2 - t1).count()/1e9<< endl<< endl;
        //cout << "Status: " << model.get(GRB_IntAttr_Status) << endl;
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            solveval.resize(n);
            cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
            for(int i=0;i<n;++i)solveval[i] = vars[i].get(GRB_DoubleAttr_X) ;
        }


    }catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }



}


void Graphs::solvePlane(){

    plane2nActiveVars.resize(nP);
    plane2nUnActiveVars.resize(nP);

    nplaneActiveVertices.resize(nP);
    planeActiveVertices.resize(nP);

    nplaneVars.resize(nP);
    varInvOrder.resize(nP);
    varOrder.resize(nP);
    pEn.resize(nP);
    pTran.resize(nP);

    plane2nActiveVertices.resize(nP);
    plane2nUnActiveVertices.resize(nP);



    for(int iplane = 0; iplane < nP; ++iplane){



        vector<int>vind(nV,-1);
        int ind = 0;
        vector<nodeinfo>nodes;
        vector<int>nodeorder;
        for(int i=0;i< plane2vers[iplane].size();++i){
            int iPnt = plane2vers[iplane][i];
            if(nodeMarkers[iPnt] != 0 ){
                vector<double>ave;
                int avelabel = averageAcrossPlanes(newval,iPnt,ave);
                nodes.push_back(nodeinfo(oldval[iPnt].data()+iplane*nMat,iPnt,avelabel,true));
                vind[iPnt] = ind++;
                nodeorder.push_back(i);
            }

        }
        int cutInd = nodes.size();
        planeActiveVertices[iplane] = nodeorder;


        for(int i=0;i< plane2vers[iplane].size();++i){
            int iPnt = plane2vers[iplane][i];
            if(nodeMarkers[iPnt] == 0 ){
                int avelabel = labelonPlane(newval,iPnt,iplane);
                nodes.push_back(nodeinfo(oldval[iPnt].data()+iplane*nMat,iPnt,avelabel,false));
                vind[iPnt] = ind++;
                nodeorder.push_back(i);
            }
        }

        vector< triple >graphconnection;
        vector<double>central(nodes.size(),0);
        auto &EE = E[iplane];
        int ccc = 0;
        for(int j=0;j<EE.size()/2;++j){

            auto p_ev = EE.data()+j*2;
            int v0 = vind[p_ev[0]];
            int v1 = vind[p_ev[1]];

            double edge_val = -1;//use cottan weigh here

            central[v0] -= edge_val;central[v1] -= edge_val;

            graphconnection.push_back(triple(v0,v1,edge_val));
            graphconnection.push_back(triple(v1,v0,edge_val));

            //graphconnection[v0][v1] = graphconnection[v1][v0] = 1;

        }
        for(int j=0;j<graphconnection.size();++j){
            graphconnection[j].val/=sqrt(central[graphconnection[j].i]*central[graphconnection[j].j]);
        }

        for(int j=0;j<nodes.size();++j){
            graphconnection.push_back(triple(j,j,1));
        }

        vector<double>solveval;


        BuildMatrixPlane(nodes,graphconnection,nodeorder,cutInd,iplane);
    }



}



void Graphs::BuildMatrixPlane(vector<nodeinfo>&nodes, vector< triple >&graphconnection_triple,vector<int>&nodeorder,int cutind,int iplane){


    double alpha = lambda;

    int n_node = nodes.size();
    int n_var = n_node*nMat;
    int n_ConstrainedVar = cutind*nMat;
    int n_UnConstrainedVar = n_var - n_ConstrainedVar;

    vector<triple>cAb;

    vector<triple>Mmatrix;

    vector<T>Ts;
    for(int i=0;i<graphconnection_triple.size();++i){
        auto &pp = graphconnection_triple[i];
        Ts.push_back(T(pp.i,pp.j,pp.val));

    }

    SpMat L(n_node,n_node);
    L.setFromTriplets(Ts.begin(),Ts.end());
	L =L*L;
    cout<<"L: "<<L.nonZeros()<<endl;

    SpMat M(n_var,n_var);
    M.setZero();

    auto constructM = [&M,&L,n_node,n_var,this](int m1,int m2){
        vector<T>Ts;
        SpMat A(n_node,n_var);
        for(int i=0;i<n_node;++i){
            Ts.push_back(T(i,m1+i*nMat,1));
            Ts.push_back(T(i,m2+i*nMat,-1));
        }
        A.setFromTriplets(Ts.begin(),Ts.end());
        M += A.transpose() * L * A;
    };


    vector<vector<int> >matset;

    Subset(nMat,2,matset);

    for(auto &a:matset)constructM(a[0],a[1]);



    Ts.clear();
    for(int i=0;i<n_var;++i)Ts.push_back(T(i,i,alpha));
    SpMat Eye(n_var,n_var);
    Eye.setFromTriplets(Ts.begin(),Ts.end());

    M += Eye;

    //cout<<"M.size(): "<<M.nonZeros()<<endl;

    SpMat dcomA = M.block(0,0,n_ConstrainedVar,n_ConstrainedVar);
    SpMat dcomC = M.block(n_ConstrainedVar,0,n_UnConstrainedVar,n_ConstrainedVar);
    SpMat dcomD = M.block(n_ConstrainedVar,n_ConstrainedVar,n_UnConstrainedVar,n_UnConstrainedVar);

    //cout<<  "N_val: "<<n_var<<endl;
    //cout<<"dcomA.size(): "<<dcomA.nonZeros()<<endl;
    //cout<<"dcomC.size(): "<<dcomC.nonZeros()<<endl;
    //cout<<"dcomD.size(): "<<dcomD.nonZeros()<<endl;

    // DC = D\C

    Eigen::CholmodSupernodalLLT< SpMat > solver;

    cout<<"Cholmod solve begin"<<endl;


    auto t1 = Clock::now();
    solver.compute(dcomD);
    SpMat DC = solver.solve(dcomC);
    auto t2 = Clock::now();

    cout << "Total solve time: " << std::chrono::nanoseconds(t2 - t1).count()/1e9<< endl;


    cout<<"Cholmod solve end"<<endl;

    SpMat newM = dcomA - dcomC.transpose()*DC;
    //cout << newM.nonZeros() <<" "<<n_ConstrainedVar*n_ConstrainedVar<<endl;



    plane2nActiveVars[iplane] = n_ConstrainedVar;
    plane2nUnActiveVars[iplane] = n_UnConstrainedVar;

    plane2nActiveVertices[iplane] = cutind;
    plane2nUnActiveVertices[iplane] = n_node-cutind;

    nplaneActiveVertices[iplane] = cutind;
    nplaneVars[iplane] = n_var;

    varOrder[iplane].resize(nplaneVars[iplane]);
    varInvOrder[iplane].resize(nplaneVars[iplane]);
    for(int j=0;j<n_node;++j){
        int ipnt_plane = nodeorder[j];
        for(int k=0;k<nMat;++k)varOrder[iplane][j*nMat+k] = ipnt_plane*nMat+k;
    }
    for(int j=0;j<nplaneVars[iplane];++j)varInvOrder[iplane][varOrder[iplane][j]] = j;

    auto &En = pEn[iplane];
    En.resize(n_ConstrainedVar*n_ConstrainedVar);


    for (int k=0; k<newM.outerSize(); ++k)
        for (SpMat::InnerIterator it(newM,k); it; ++it)
        {
            En[it.row() + it.col() * n_ConstrainedVar ] = it.value();
        }

    auto &Tran = pTran[iplane];
    Tran.resize(n_UnConstrainedVar*n_ConstrainedVar);


    for (int k=0; k<DC.outerSize(); ++k)
        for (SpMat::InnerIterator it(DC,k); it; ++it)
        {
            Tran[it.row() + it.col() * n_UnConstrainedVar] = it.value();
        }


}



bool Graphs::CrossectionsPipeline(CrossSections &CS, string outpath, double lambda, int maxiter){




    readGraphsFromCrossections(CS);

    newval = oldval;

    this->lambda = lambda;

    solveLine();

    solvePlane();

    //GetContour2(solval,outpath);

    doubleIter_opt(maxiter);

    GetContour2(solval,outpath);

}

