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



void Graphs::doubleIter_opt(int maxIter){

    isInit=false;

    if(1){

        double curobj;
        int iter;
        typedef std::chrono::high_resolution_clock Clock;
        auto t1 = Clock::now();
        

        vector<double>objs;

        r_iters.clear();

        solval = newval;

        IterInner_prebuild(true);

        InitialobjModel();

        InitializeMarkers();

        double preobj = solValues_modifyingConstraints_flippinglabel_QPonly();
        objs.push_back(preobj);


        newval = solval;

        for(iter=0;iter<maxIter;++iter){

            gradientvector = gradientvector_tmp;


            curobj = EpsilonIterwithMarkersandSelectiveLabels_Inner(preobj);


            newval = solval;

            objs.push_back(preobj);
            if(curobj==-1){++iter;curobj=preobj;break;};
            preobj = curobj;

        }

        auto t2 = Clock::now();
        auto totaltime = std::chrono::nanoseconds(t2 - t1).count()/1e9;

        cout<<"Total Outer Iteration: "<<iter<<endl;
        cout<<"Total time:      " << totaltime << endl;
        cout<<"Obj:             "<<curobj<<endl;

        for(auto a:r_iters)cout<<a<<' ';cout<<endl;
        int totaliter = 0;
        for(auto a:r_iters)totaliter+=a;
        cout<<"Total Inner Iteration: "<<totaliter<<endl;
        cout<<"Ave inner time: "<<totaltime/totaliter<<endl;
        cout<<setprecision(10);
        for(auto a:objs)cout<<a<<' ';cout<<endl;

    }


}


void Graphs::IterInner_prebuild(bool isrefresh){

    if(isrefresh){
        nTotalActiveV = 0;
        mapV2ActiveV.clear();
        mapV2ActiveV.resize(nV,-1);



        for(int i=0;i<nP;++i){
            auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
            for(auto a:planeActiveVertices[i])mapV2ActiveV[mapplaneV2globalV_local[a]] = 0;
        }
        for(auto &a:mapV2ActiveV)if(a!=-1)a=nTotalActiveV++;
        for(int i=0;i<nV;++i)if(mapV2ActiveV[i]!=-1)globalactiveV.push_back(i);



        planev2oldval.resize(nP);
        for(int i=0;i<nP;++i){
            auto &pvov = planev2oldval[i];
            auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
            int npv = mapplaneV2globalV_local.size();
            pvov.resize(npv*nMat);
            for(int j=0;j<npv;++j){
                auto pv = pvov.data()+j*nMat;
                auto pvo = oldval[mapplaneV2globalV_local[j]].data()+i*nMat;
                for(int k=0;k<nMat;++k){
                    pv[k] = pvo[k];
                }
            }
        }
    }

    planev2newval.resize(nP);
    for(int i=0;i<nP;++i){
        auto &pvov = planev2newval[i];
        auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
        int npv = mapplaneV2globalV_local.size();
        pvov.resize(npv*nMat);
        for(int j=0;j<npv;++j){
            auto pv = pvov.data()+j*nMat;
            auto pvo = newval[mapplaneV2globalV_local[j]].data()+i*nMat;
            for(int k=0;k<nMat;++k){
                pv[k] = pvo[k];
            }
        }
    }

    FirstlabelsOnNewVal.resize(nV,-1);
    SecondlabelsOnNewVal.resize(nV,-1);
    for(int i=0;i<nV;++i){
        auto &pval = newval[i];
        vector<double>ave(nMat,0);
        for(int j=0;j<nP;++j){
            auto pv = pval.begin()+nMat*j;
            bool isNull = true;
            for(int k=0;k<nMat;++k)isNull  = isNull&&(pv[k]==-1);
            if(isNull)continue;
            for(int k=0;k<nMat;++k)ave[k]+=pv[k];
        }

        FirstlabelsOnNewVal[i] = max_element(ave.begin(),ave.end())-ave.begin();
        ave[FirstlabelsOnNewVal[i]] = -1e20;
        SecondlabelsOnNewVal[i] = max_element(ave.begin(),ave.end())-ave.begin();

    }





    VertexMarker.clear();



    VertexMarker.resize(nTotalActiveV,-1);

    nFlipping = 0;




    for(int i=0;i<nV;++i){
        int vbinary = mapV2ActiveV[i];
        if(vbinary==-1)continue;
        //VertexMarker[vbinary] = nodeMarkers[i]==2?0:-1;
        VertexMarker[vbinary] = nodeMarkers[i]!=0?0:-1;
    }

    if(isrefresh){
        nFlipping = 0;
        for(auto &a:VertexMarker)if(a==1)a = nFlipping++;else a = -1;
    }
}

void Graphs::InitialobjModel(){

    try{
        conVars.resize(nP);

        int codeer = 0;
        for(int i=0;i<nP;++i){
            auto &conVar = conVars[i];
            conVar.resize(plane2nActiveVars[i]);
            auto &pvov = planev2newval[i];
            auto &pvovold = planev2oldval[i];
            auto &pvarorder = varOrder[i];

            for(int j=0;j<plane2nActiveVertices[i];++j){
                for(int k=0;k<nMat;++k){
                    int pos = pvarorder[j*nMat+k];
                    double initval = pvov[pos]-pvovold[pos];
                    conVar[j*nMat+k] = objmodel.addVar(-GRB_INFINITY, GRB_INFINITY, isInit?initval:0, GRB_CONTINUOUS,to_string(codeer++));
                }
            }
        }


        //cout << "optimize 1" << endl;
        GRBQuadExpr obj = 0;

        for(int k=0;k<nP;++k){

            auto &en = pEn[k];
            auto &variables = conVars[k];
            int nCon = plane2nActiveVars[k];
            for(int i=0;i<nCon;++i){
                for(int j=0;j<nCon;++j){
                    double a = en[i*nCon+j];
                    obj+=a*variables[i]*variables[j];
                }
            }

        }


        //cout << "optimize 2 " <<codeer<< endl;
        objmodel.setObjective(obj, GRB_MINIMIZE);


        objmodel.update();



    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }

}

void Graphs::InitializeMarkers(){


    previousFlipGradient.clear();
    previousEnergy.clear();

    //previousFlipGradient.resize(nV,1e10);
    previousEnergy.resize(nV,-1);

    flipGradient.clear();
    flipGradient.resize(nV,1e10);

}


double Graphs::solValues_modifyingConstraints_flippinglabel_QPonly(){

    double theta = 1e-4;

    //solval =  oldval;
    //for(auto &a:solval)for(auto &b:a)b=-1;

    planev2solveval.resize(nP);
    isChanged.clear();
    isChanged.resize(nV,0);


    try{

        //cout << "optimize 0" << endl;

        typedef std::chrono::high_resolution_clock Clock;
        auto t01 = Clock::now();

        GRBModel model(objmodel);
        //cout << "optimize 0" << endl;
        // Quiet output from Gurobi.

        GRBVar *pVar = model.getVars();
        int codeer = 0;
        for(int i=0;i<nP;++i){
            auto &conVar = conVars[i];
            conVar.resize(plane2nActiveVars[i]);
            for(int j=0;j<plane2nActiveVertices[i];++j){
                for(int k=0;k<nMat;++k){
                    conVar[j*nMat+k] = pVar[codeer++];
                }
            }
        }



        auto t012 = Clock::now();
        cout<<"Setting up Obj time:      " << std::chrono::nanoseconds(t012 - t01).count()/1e9<< endl;



        model.update();

        {

            cout<<"setting up constraint"<<endl;

            for(int k=0;k<nP;++k){

                int bc = 0,lc = 0;
                int nConvertice = plane2nActiveVertices[k];
                auto &mapplaneV2globalV_local = mapplaneV2globalV[k];
                auto &oldvalue = planev2oldval[k];
                auto &pvarOrder = varOrder[k];
                auto &variables = conVars[k];
                auto &planeActiveVertices_local = planeActiveVertices[k];
                for(int v = 0;v<nConvertice;++v){
                    int bias = nMat*v;
                    int globalvind = mapplaneV2globalV_local[planeActiveVertices_local[v]];
                    int vbinary = mapV2ActiveV[ globalvind ];

                    //flip or not
                    if(VertexMarker[vbinary]!=-1){

                        int picki = SecondlabelsOnNewVal[globalvind];
                        //cout<<pvarOrder[bias]<<endl;
                        auto olvaluei = oldvalue[pvarOrder[bias+picki]];
                        for(int j=0;j<nMat;++j)if(j!=picki){
                            auto olvaluej = oldvalue[pvarOrder[bias+j]];
                            model.addConstr((variables[bias+picki]+olvaluei)-(variables[bias+j]+olvaluej)>=theta);
                            //cout<<bias+j<<endl;
                        }

                        bc++;
                    }else{

                        lc++;
                        //continue;
                        int picki = FirstlabelsOnNewVal[globalvind];
                        //cout<<pvarOrder[bias]<<endl;
                        auto olvaluei = oldvalue[pvarOrder[bias+picki]];
                        for(int j=0;j<nMat;++j)if(j!=picki){
                            auto olvaluej = oldvalue[pvarOrder[bias+j]];
                            model.addConstr((variables[bias+picki]+olvaluei)-(variables[bias+j]+olvaluej)>=theta);
                            //cout<<bias+j<<endl;
                        }


                    }
                }
                //cout<<bc<<' '<<lc<<endl;
            }


        }




        cout << "QP optimize begin" << endl;
        cout<< "Total continuous variables: "<< nTotalActiveV*nMat <<endl;
        cout << "nActiveV, nMat, nFlipping: " << nTotalActiveV << ' '<<nMat <<' '<<nFlipping<<endl;

        auto t1 = Clock::now();
        model.optimize();
        auto t2 = Clock::now();

        cout << "optimize end" << endl;
        cout << "Total opt time: " << std::chrono::nanoseconds(t2 - t1).count()/1e9<< endl<< endl;

        // Extract Results.

        cout << "Status: " << model.get(GRB_IntAttr_Status) << endl;
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            cout << "Optimal"<<endl;
        }


        return extractResultfromGurobi(model);



    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }

}



double Graphs::IterInnerwithMarkersandSelectiveLabels_mainIter(double beginningEn,vector<vector<int>>&comps){

    if(comps.size()==0){
        r_iters.push_back(comps.size());
        return -1;
    }
    auto batchOperation = [this](){
        nFlipping = 0;
        for(auto &a:VertexMarker)if(a==1)a = nFlipping++;else a = -1;
    };
    double preobj = beginningEn, curobj = 0;

    auto backup_VertexMarker = VertexMarker;

    cout<<"comps size: "<<comps.size()<<endl;
    for(int i=0;i<comps.size();++i){

        cout<<"comps[i][j]: "<<endl;
        for(int j=0;j<comps[i].size();++j){
            int vi = comps[i][j];

            VertexMarker[mapV2ActiveV[vi]] = 1;

            cout<<vi<<": ";

        }

        batchOperation();
        //curobj = solValues(2,isInit);

        curobj = solValues_modifyingConstraints_flippinglabel_QPonly();

        previousEnergy[comps[i][0]] = curobj-preobj;
        if(preobj-1e-8>curobj){
            //newval = solval;
            r_iters.push_back(i+1);
            return curobj;

        }else{
            VertexMarker = backup_VertexMarker;
        }


    }
    r_iters.push_back(comps.size());
    return -1;


}



double Graphs::EpsilonIterwithMarkersandSelectiveLabels_Inner(double beginningEn){

    static int counter = 0;

    IterInner_prebuild(false);


    vector<int>markers(nV,-1);


    vector<bool>epsilon(nV,false);

    for(int i=0;i<nP;++i){
        auto &pvov = planev2newval[i];
        auto &pvovold = planev2oldval[i];
        auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
        int npv = mapplaneV2globalV_local.size();
        pvov.resize(npv*nMat);
        //int nnnn=0;
        for(int j=0;j<npv;++j){
            int vglobalind = mapplaneV2globalV_local[j];
            int vbinary = mapV2ActiveV[mapplaneV2globalV_local[j]];
            if(vbinary==-1)continue;
            if(VertexMarker[vbinary]!=0)continue;

            vector<double>valuevec(nMat);
            auto pv = pvov.begin()+j*nMat;
            auto pvold = pvovold.begin()+j*nMat;

            for(int nn =0;nn<nMat;++nn)valuevec[nn] = pv[nn];

            sort(valuevec.begin(),valuevec.end(),greater<double>());
            //for(int k=0;k<nMat;++k)if(fabs(pvold[k])<1e-3)ko=true;
            //cout<<nnnn++<<' ';for(int k=0;k<nMat;++k)cout<<valuevec[k]<<' ';cout<<endl;

            if(fabs(valuevec[0]-valuevec[1])<1.01e-4){
                epsilon[vglobalind] = true;
            }


        }
    }




    vector<vector<int>>comps;
    for(int i=0;i<nV;++i)if(epsilon[i])comps.push_back(vector<int>({i}));

    comps = orderingComps_markers(comps);

    ++counter;
    return IterInnerwithMarkersandSelectiveLabels_mainIter(beginningEn,comps);


}


double Graphs::extractResultfromGurobi(GRBModel &model){


    vector<int>planev2numchangeL(nP,0);
    vector<int>sollabel(nV,-1);
    vector<bool>vialateV(nV,false);
    int nviolation = 0;

    for(int i=0;i<nP;++i){
        vector<double>gurobisolvec(plane2nActiveVars[i]);
        auto &pvov = planev2solveval[i];
        auto &pvovold = planev2oldval[i];
        auto &pvovnew = planev2newval[i];
        auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
        int npv = mapplaneV2globalV_local.size();
        pvov.resize(npv*nMat);
        auto &conVar = conVars[i];

        auto &pvarOrder = varOrder[i];
        for(int j=0;j<plane2nActiveVertices[i];++j){
            auto psolv = conVar.data()+j*nMat;

            //auto pvo = solval[mapplaneV2globalV_local[localind]].data()+i*nMat;
            for(int k=0;k<nMat;++k){
                auto pv = pvov.data()+pvarOrder[j*nMat+k];
                auto pvold = pvovold.data()+pvarOrder[j*nMat+k];
                gurobisolvec[j*nMat+k] =  psolv[k].get(GRB_DoubleAttr_X);
                pv[0] =  psolv[k].get(GRB_DoubleAttr_X) + pvold[0];
            }
        }

        auto &pTra = pTran[i];
        int w = plane2nUnActiveVars[i],h = plane2nActiveVars[i];
        vector<double>matrixMul(w);
        for(int j=0;j<w;++j){
            double aaa = 0;
            for(int k=0;k<h;++k)aaa+=pTra[k*w+j]*gurobisolvec[k];
            matrixMul[j] = -aaa;
        }
        int offset = plane2nActiveVertices[i]*nMat;
        for(int j=0;j<plane2nUnActiveVertices[i];++j){
            int uoffset = offset+j*nMat;
            auto pmat = matrixMul.data()+j*nMat;

            for(int k=0;k<nMat;++k){
                auto pv = pvov.data()+pvarOrder[uoffset+k];
                auto pvold = pvovold.data()+pvarOrder[uoffset+k];
                pv[0] = pmat[k]+pvold[0];
                //cout<<pvnew[k] - pv[k]<<' ';
            }
        }



        for(int j=0;j<mapplaneV2globalV_local.size();++j){
            auto pglobalsol = solval[mapplaneV2globalV_local[j]].data()+i*nMat;
            auto pvsol = pvov.begin()+j*nMat;
            for(int k=0;k<nMat;++k){
                pglobalsol[k] = pvsol[k];
            }

            int label = max_element(pglobalsol,pglobalsol+nMat)-pglobalsol;
            if(sollabel[mapplaneV2globalV_local[j]]!=-1)if(sollabel[mapplaneV2globalV_local[j]]!=label){
                ++nviolation;
                vialateV[mapplaneV2globalV_local[j]] = true;
            }
            sollabel[mapplaneV2globalV_local[j]]=label;

        }


        if(1)for(int j=0;j<plane2nActiveVertices[i]+plane2nUnActiveVertices[i];++j){
            auto pvsol = pvov.begin()+j*nMat;
            auto pvnew = pvovnew.begin()+j*nMat;

            int a1 = max_element(pvsol,pvsol+nMat)-pvsol;
            int a2 = max_element(pvnew,pvnew+nMat)-pvnew;
            if(a1!=a2){
                planev2numchangeL[i]++;
                if(j<plane2nActiveVertices[i])isChanged[mapplaneV2globalV_local[j]] = 2;
                else isChanged[mapplaneV2globalV_local[j]] = 1;
            }

        }

    }


    if(0)for(int i=0;i<nV;++i)if(vialateV[i]){
        printinfo(i,solval);
    }


    gradientvector_tmp.resize(nP);
    for(int k=0;k<nP;++k){
        int nCon = plane2nActiveVars[k];
        auto &conVar = conVars[k];
        vector<double>gurobisolvec(nCon);
        for(int j=0;j<nCon;++j){
            gurobisolvec[j] =  conVar[j].get(GRB_DoubleAttr_X);
        }
        auto &en = pEn[k];
        vector<double>tmp(nCon,0);

        for(int i=0;i<nCon;++i){
            for(int j=0;j<nCon;++j){
                double a = en[i*nCon+j];
                tmp[i]+=a*gurobisolvec[j];
            }
        }
        gradientvector_tmp[k] = tmp;
    }


    cout<<"obj:         "<<model.get(GRB_DoubleAttr_ObjVal)<<endl;
    cout<<"nviolation:  "<<nviolation<<endl;
    //if(nviolation != 0)cout<<"The solution contains inconsistent result, please consider enlarge the active region."<<endl;
    cout<<"nplaneV:     ";for(int i=0;i<nP;++i)cout<<plane2nActiveVertices[i]+plane2nUnActiveVertices[i]<<'\t';cout<<endl;
    cout<<"lable change:";for(int i=0;i<nP;++i)cout<<planev2numchangeL[i]<<'\t';cout<<endl;
    cout<<"Portion:     ";for(int i=0;i<nP;++i)cout<<double(planev2numchangeL[i])/(plane2nActiveVertices[i]+plane2nUnActiveVertices[i])<<'\t';cout<<endl;

    return model.get(GRB_DoubleAttr_ObjVal);


}


vector<vector<int>> Graphs::orderingComps_markers(vector<vector<int>>&in_comps){


    //double flipingstep = 1e-3;


    computeFlipGradient();

    flipEnergy.clear();
    flipEnergy.resize(nV,0);



    vector<pair<double,int>>c2info(in_comps.size());
    for(int i=0;i<in_comps.size();++i){
        c2info[i].second = i;


        double val = 1e10;
        for(auto a:in_comps[i])val = fmin(val,flipGradient[a]);

        c2info[i].first = val;
    }

    sort(c2info.begin(), c2info.end(), [](pair<double,int> &left, pair<double,int> &right) {
        return left.first < right.first;
    });

    vector<vector<int>>sortcomp;

    for(int i=0;i<in_comps.size();++i){

        int index = in_comps[c2info[i].second][0];
        if(fabs(previousFlipGradient[index]-flipGradient[index])<1e-3)if(previousEnergy[index]>0)continue;
        sortcomp.push_back( in_comps[c2info[i].second] );
    }

    for(auto &a:sortcomp)for(auto b:a)cout<<b<<' ';cout<<endl;

    return sortcomp;


}


void Graphs::computeFlipGradient(){


    swap(previousFlipGradient,flipGradient);
    flipGradient.clear();
    flipGradient.resize(nV,0);

    for(int k=0;k<nP;++k){

        int bc = 0,lc = 0;
        int nConvertice = plane2nActiveVertices[k];

        auto &mapplaneV2globalV_local = mapplaneV2globalV[k];

        auto &newvalue = planev2newval[k];
        auto &pvarOrder = varOrder[k];

        auto &planeActiveVertices_local = planeActiveVertices[k];
        auto &gradient = gradientvector[k];
        for(int v = 0;v<nConvertice;++v){
            int bias = nMat*v;
            int globalvind = mapplaneV2globalV_local[planeActiveVertices_local[v]];

            int firstlabel = FirstlabelsOnNewVal[globalvind];
            int secondlabel = SecondlabelsOnNewVal[globalvind];

            double firstgradient = gradient[bias+firstlabel];
            double secondgradient = gradient[bias+secondlabel];

            if(newvalue[pvarOrder[bias+firstlabel]]-newvalue[pvarOrder[bias+secondlabel]]<1.01e-4)
                cout<<globalvind<<' '<<firstgradient<<' '<<secondgradient<<' '<<flipGradient[globalvind]<<' '<<firstlabel<<' '<<secondlabel<<endl;

            //if(fabs(flipEnergy[globalvind])>1e-6)continue;
            flipGradient[globalvind] += -firstgradient+secondgradient;
            //flipGradient[globalvind] = min(min(-firstgradient,secondgradient),flipGradient[globalvind]);

            //cout<<v<<endl;

        }
        //cout<<bc<<' '<<lc<<endl;
    }



}


void Graphs::printinfo(int globalv, vector<vector<double> >&val){

    cout<<"v index: "<<globalv<<endl;
    cout<<"nodeMarkers: "<<nodeMarkers[globalv]<<endl;
    cout<<"v planes: "; for(auto a:ver2planes[globalv])cout<<a<<" ";cout<<endl;
    for(int j=0;j<nP;++j){
        auto pglobalsol = val[globalv].data()+j*nMat;
        int label = max_element(pglobalsol,pglobalsol+nMat)-pglobalsol;
        cout<<label<<' ';
        for(int k=0;k<nMat;++k)cout<<pglobalsol[k]<<' ';cout<<endl;
    }

}
