#include "contour.h"
#include<iostream>
#include<fstream>
#include<set>
#include<assert.h>
#include<map>
#include"UnionFind.h"
#include <iomanip>

void Contour::ReadContour(ifstream &in){

    for(int i=0;i<4;++i)in>>plane_para[i];
    in>>n_vertices;
    in>>n_edges;

    vertices.resize(n_vertices*3);
    edges.resize(n_edges*2);
    m1.resize(n_edges);
    m2.resize(n_edges);

    for(int i=0;i<vertices.size();++i)in>>vertices[i];
    for(int i=0;i<n_edges;++i){
        in>>edges[i*2];
        in>>edges[i*2+1];
        in>>m1[i];
        in>>m2[i];
    }

    set<int>mmm;
    for(auto a:m1)mmm.insert(a);
    for(auto a:m2)mmm.insert(a);

    if(0)n_materials = mmm.size();
    else n_materials = *max_element(mmm.begin(),mmm.end())+1;
    //for(auto )
    //cin>>n_materials;

}
void Contour::splitContour(vector<Contour> &outC){

    if(n_vertices==0 || n_edges==0 || n_materials<=0)return;

    outC.resize(n_materials);
    for(int i=0;i<n_materials;++i){

        vector<int>vpick(n_vertices,-2);
        vector<bool>epick(n_edges,false);
        auto &pvertices = outC[i].vertices;
        auto &pedges = outC[i].edges;
        auto &pm1 = outC[i].m1;
        auto &pm2 = outC[i].m2;

        for(int j = 0; j<n_edges;++j){
            if(m1[j]==i || m2[j] == i)epick[j] = true;
            else epick[j] = false;
        }

        for(int j = 0; j<n_edges;++j){
            if(epick[j]){
                auto p_ev = ev_begin(j);
                vpick[p_ev[0]] = vpick[p_ev[1]] = -1;
            }
        }

        int vind = 0;
        for(int j = 0; j<n_vertices;++j){
            if(vpick[j]==-1){vpick[j] = vind;++vind;}
        }

        int eind = 0;
        for(int j = 0; j<n_edges;++j){
            if(m1[j]==i || m2[j] == i){epick[j] = true;++eind;}
            else epick[j] = false;
        }

        outC[i].n_vertices = vind;
        outC[i].n_edges = eind;
        outC[i].n_materials = 2;
        pvertices.resize(0);
        pedges.resize(0);

        pm1.resize(0);
        pm2.resize(0);

        for(int j=0;j<4;++j)outC[i].plane_para[j] = plane_para[j];
        for(int j = 0; j<n_edges;++j)if(epick[j]){
            auto p_ev = ev_begin(j);
            pedges.push_back(vpick[p_ev[0]]);
            pedges.push_back(vpick[p_ev[1]]);

            if(i==0){
                pm1.push_back(m1[j]==0?0:1);
                pm2.push_back(m2[j]==0?0:1);
            }else{
                pm1.push_back(m1[j]==i?1:0);
                pm2.push_back(m2[j]==i?1:0);
            }
        }


        for(int j = 0; j<n_vertices;++j){
            if(vpick[j]>-1){
                auto p_v = v_begin(j);
                for(int k=0;k<3;++k)pvertices.push_back(p_v[k]);


            }
        }

    }


}
void Contour::WriteContour(ofstream &out){

    for(int i=0;i<4;++i)out<<plane_para[i]<<' ';out<<endl;
    out<<n_vertices<<' ';
    out<<n_edges<<endl;



    for(int i=0;i<n_vertices;++i){
        int ind = i*3;
        for(int j=0;j<3;++j)out<<vertices[ind+j]<<' ';
        out<<endl;
    }
    for(int i=0;i<n_edges;++i){
        out<<edges[i*2]<<' ';
        out<<edges[i*2+1]<<' ';
        out<<m1[i]<<' ';
        out<<m2[i]<<endl;
    }

}
void Contour::WriteContourCtrFormat(ofstream &out){

    vector<vector<double>>verticesSeqs;
    vector<int>loops2mats;
    vector<vector<int>>verticesloops;

    ComputeMat2VerticesLoop(verticesloops,verticesSeqs,loops2mats);


    vector<double>tmp(3,1);
    for(int j=0;j<3;++j)if(fabs(plane_para[j])>1e-1){
        int v1,v2;
        if(j==0){v1=1;v2=2;}else if(j==1){v1=0;v2=2;}else if(j==2){v1=0;v2=1;}
        tmp[j] = (0-plane_para[v1]*tmp[v1]-plane_para[v2]*tmp[v2])/plane_para[j];
        break;
    }
    normalize(tmp.data());

    for(int i=0;i<loops2mats.size()/2;++i){

        for(int i=0;i<4;++i)out<<plane_para[i]<<' ';
        for(int j=0;j<3;++j)out<<verticesSeqs[i][j]<<' ';
        for(int j=0;j<2;++j)out<<tmp[j]<<' ';out<<tmp[2]<<endl;

        auto p_v = verticesSeqs[i].data();
        int nv = verticesSeqs[i].size()/3;

        out<<nv<<' ';
        out<<nv<<endl;

        for(int k=0;k<nv;++k){
            int ind = k*3;
            for(int j=0;j<2;++j)out<<p_v[ind+j]<<' ';
            out<<p_v[ind+2]<<endl;
        }
        for(int k=0;k<nv-1;++k){
            out<<k<<" "<<k+1<<" "<<0<<" "<<1<<endl;
        }
        out<<nv-1<<" "<<0<<" "<<0<<" "<<1<<endl;
    }

}

void Contour::WriteContourCSLFormat(ofstream &out){

    vector<vector<double>>verticesSeqs;
    vector<int>loops2mats;
    vector<vector<int>>verticesloops;

    ComputeMat2VerticesLoop(verticesloops,verticesSeqs,loops2mats);


    out<<n_vertices<<' '<<verticesloops.size()-1<<' ';
    for(int i=0;i<4;++i)out<<plane_para[i]<<' ';
    out<<endl;
    out<<endl;

    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;

        out<<p_v[0]<<' '<<p_v[1]<<' '<<p_v[2]<<endl;

    }
    out<<endl;
    out<<endl;


    for(int i=0;i<verticesloops.size();++i)if(loops2mats[i]!=0){
        out<<verticesloops[i].size()<<' '<<loops2mats[i]<<' ';
        for(int j=verticesloops[i].size()-1;j>=0;--j)out<<verticesloops[i][j]<<' ';
        out<<endl;
        out<<endl;
    }

}




void Contour::CalSelfProperty(){

    leftcorner.clear();
    rightcorner.clear();
    midpoint.clear();
    leftcorner.resize(3,1000000);
    rightcorner.resize(3,-1000000);
    midpoint.resize(3,0);

    n_edges = edges.size()/2;
    n_vertices = vertices.size()/3;
    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;
        for(int j=0;j<3;++j){
            leftcorner[j] = min(leftcorner[j],p_v[j]);
            rightcorner[j] = max(rightcorner[j],p_v[j]);
            midpoint[j] += p_v[j];
        }
    }

    for(int j=0;j<3;++j){

        midpoint[j] /= n_vertices;
    }

    plane_info.build(midpoint,plane_para);
}
bool Contour::ValidationCheck(){
    for(auto a:edges)assert(a<n_vertices);
    double baa=0;
    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        if((baa = VerticesDistance(p_ev[0],p_ev[1]))<1e-6){
            cout<<"Contour degenerate edge: "<<baa<<' '<<p_ev[0]<<' '<<p_ev[1]<<endl;
        }
    }
    cout<<"Contour validation check"<<endl;

}

void Contour::BuildContourFromImage(vector<vector<int>>&im, double z){

    /*
    if(im.size()==0)return;

    plane_para[0] = 0;plane_para[1] = 0;plane_para[2] = 1;plane_para[3] = -z;

    int height = im.size();
    int width = im[0].size();

    int bw = width*2-1;
    int bh = height*2-1;

    vector<vector<int>>fbuf(bh);
    for(int i=0;i<bh;++i){
        fbuf[i].resize(bw,-1);
    }


    for(int i = 0; i<height-1;++i ){
        auto &buf1 = fbuf[2*i];
        auto &im1 = im[i];
        for(int j = 0; j<width-1;++j ){
            if(im1[j]!=im1[j+1])buf1[j*2+1] = 0;
        }
    }
    for(int i = 0; i<height-1;++i ){
        auto &buf1 = fbuf[2*i+1];
        auto &im1 = im[i];
        auto &im2 = im[i+1];
        for(int j = 0; j<width-1;++j ){
            if(im1[j]!=im2[j])buf1[j*2] = 0;
        }
    }
    for(int i = 0; i<height-1;++i ){
        auto &buf1 = fbuf[2*i];
        auto &buf2 = fbuf[2*i+1];
        auto &buf3 = fbuf[2*i+2];
        for(int j = 1; j<bw;j+=2 ){
            if(buf1[j]==0 || buf3[j]==0)buf2[j] = 0;
        }
        for(int j = 0; j<bw-2;j+=2 ){
            if(buf2[j]==0 || buf2[j+2]==0)buf2[j+1] = 0;
        }
    }
    int vind = 0;
    vertices.clear();
    for(int i=1;i<bh;i+=2){
        auto &buf1 = fbuf[i];
        for(int j = 1; j<bw;j+=2 )if(buf1[j]==0){
            vertices.push_back(j);
            vertices.push_back(i);
            vertices.push_back(z);
            buf1[j] = vind++;
        }
    }

    edges.clear();m1.clear();m2.clear();
    for(int i = 0; i<height-1;++i ){
        auto &buf1 = fbuf[2*i+1];
        auto &im1 = im[i];
        auto &im2 = im[i+1];
        for(int j = 1; j<width-1;++j ){
            if(im1[j]!=im2[j]){
                edges.push_back(buf1[j*2-1]);
                edges.push_back(buf1[j*2+1]);
                m1.push_back(im2[j]);
                m2.push_back(im1[j]);

            }
        }
    }
    for(int i = 1; i<height-1;++i ){
        auto &buf1 = fbuf[2*i-1];
        auto &buf2 = fbuf[2*i+1];
        auto &im1 = im[i];
        for(int j = 0; j<width-1;++j ){
            if(im1[j]!=im1[j+1]){
                edges.push_back(buf1[j*2+1]);
                edges.push_back(buf2[j*2+1]);
                m1.push_back(im1[j]);
                m2.push_back(im1[j+1]);

            }
        }
    }


    n_vertices = vertices.size()/3;
    n_edges = edges.size()/2;
    if(1){
        vector<uint>edgesMat;
        vector<int>edges2Cell(n_edges*2,0);
        for(int i=0;i<n_edges;++i)edges2Cell[i*2]=1;
        for(int i=0;i<n_edges;++i){
            edgesMat.push_back(m1[i]);
            edgesMat.push_back(m2[i]);
        }
        Curve simCurve;
        simCurve.ImportCurve(vertices,edges);
        simCurve.CurveNetAnalayze(edgesMat,edges2Cell,2);
        simCurve.CurveNetSimplification(6,vertices,edges,edgesMat);//6
        n_vertices = vertices.size()/3;
        n_edges = edges.size()/2;
        m1.resize(n_edges);m2.resize(n_edges);
        for(int i=0;i<n_edges;++i){
            m1[i] = edgesMat[2*i];
            m2[i] = edgesMat[2*i+1];
        }
    }
    CurveSmoothing(800);

    */
}

void Contour::ImportContour(double *para,vector<double>&inCtrV,vector<uint>&inCtrE,vector<int>&inCtrEMat){
    for(int i=0;i<4;++i)plane_para[i] = para[i];
    vertices=inCtrV;
    edges = inCtrE;
    n_vertices = vertices.size()/3;
    n_edges = edges.size()/2;
    m1.resize(n_edges);m2.resize(n_edges);
    for(int i=0;i<n_edges;++i){
        m1[i] = inCtrEMat[2*i];
        m2[i] = inCtrEMat[2*i+1];
    }

}

void Contour::OutputContour(double *para, vector<double>&outCtrV, vector<uint>&outCtrE, vector<int>&outCtrEMat){
    for(int i=0;i<4;++i)para[i] = plane_para[i];
    outCtrV = vertices;
    outCtrE = edges;
    outCtrEMat.clear();
    for(int i=0;i<n_edges;++i){
        outCtrEMat.push_back(m1[i]);
        outCtrEMat.push_back(m2[i]);
    }


}


void Contour::CurveSmoothing(int maxiter){


    vector<vector<int>>vnn(n_vertices);


    for(int i=0;i<n_edges;++i){
        int v1 = edges[i*2],v2 = edges[i*2+1];
        vnn[v1].push_back(v2);
        vnn[v2].push_back(v1);
    }



    vector<double>buffer1(vertices.size());
    vector<double>buffer2(vertices.size());
    vector<double>*pre = &buffer1,*cur = &buffer2;
    *pre = vertices;

    double kaipa = 0.1,lamda = 0.5;
    auto gama = 1 / (kaipa - 1 / lamda);

    auto oneiter = [&vnn,this](double coef,vector<double>*pre,vector<double>*cur){
        double ocoef = 1- coef;

        vector<double>aaa(3);
        for(int i =0;i<n_vertices;++i){

            for(int j=0;j<3;++j)aaa[j] = 0;
            auto &vv = vnn[i];
            if(vv.size()!=2){
                //cout<<"Adas"<<endl;
                for(int j=0;j<3;++j)cur->at(i*3+j) = pre->at(i*3+j);
                continue;
            }

            for(int k=0;k<vv.size();++k){
                for(int j=0;j<3;++j){aaa[j]+=pre->at((vv[k]*3)+j);}
            }
            for(int j=0;j<3;++j)aaa[j]/=(vv.size());
            //for(int j=0;j<3;++j)cout<<aaa[j]<<" ";cout<<endl;
            for(int j=0;j<3;++j)cur->at(i*3+j) = ocoef * pre->at(i*3+j) + coef * aaa[j];
        }

    };

    for(int iter=0;iter<maxiter;++iter){

        //cout<<iter<<endl;
        oneiter(lamda,pre,cur);
        swap(pre,cur);
        oneiter(gama,pre,cur);
        swap(pre,cur);

    }
    vertices=(*pre);




}

void Contour::MergeBoundaryBox(vector<double>&inCtrV, vector<uint>&inCtrE, vector<int>&inCtrEMat){


    if(inCtrV.size()==0 || inCtrE.size()==0)return;
    vector<vector<uint>>v2v(n_vertices);
    vector<vector<uint>>v2e(n_vertices);

    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        v2v[p_ev[0]].push_back(p_ev[1]);
        v2v[p_ev[1]].push_back(p_ev[0]);

        v2e[p_ev[0]].push_back(i);
        v2e[p_ev[1]].push_back(i);
    }

    int nv = inCtrV.size()/3;
    int ne = inCtrE.size()/2;


    vector<vector<uint>>in_v2v(nv);
    vector<vector<uint>>in_v2e(nv);

    for(int i=0;i<ne;++i){
        auto p_ev = inCtrE.data()+i*2;
        in_v2v[p_ev[0]].push_back(p_ev[1]);
        in_v2v[p_ev[1]].push_back(p_ev[0]);

        in_v2e[p_ev[0]].push_back(i);
        in_v2e[p_ev[1]].push_back(i);
    }

    for(auto &a:in_v2v)assert(a.size()==2);






    auto closestPoint = [this,&nv](double *pv, vector<double>&inV){


        vector<double>dist(nv);
        auto p_dv = inV.data();
        for(int i=0;i<nv;++i){
            dist[i] = vecSquareDist(pv,p_dv+i*3);
        }

        return min_element(dist.begin(),dist.end())-dist.begin();

    };



    vector<int>newAddInd(nv,-1);
    vector<int>attachV;
    vector<int>attachinV;
    vector<bool>modifiedV(nv,false);
    vector<bool>modifiedE(ne,false);
    vector<uint>addE;
    vector<int>addEMat;
    for(int i=0;i<n_vertices;++i){
        assert(v2v[i].size()>0);
        if(v2v[i].size()>1)continue;

        int cInd = closestPoint(v_begin(i),inCtrV);

        assert(newAddInd[cInd]==-1);//extreme cases may appear
        modifiedV[cInd] = true;
        newAddInd[cInd] = i;
        attachV.push_back(i);
        attachinV.push_back(cInd);

    }

    cout<<"attachV: "<<attachV.size()<<endl;

    int outsideM = numeric_limits<int>::max();

    if(attachV.size()==0){
        addE = inCtrE;
        addEMat = inCtrEMat;
        for(auto &a:addEMat)if(a==1)a = outsideM;
    }else {


        for(int i=0;i<attachV.size();++i){

            int oriV = attachV[i],inV = attachinV[i];
            int ml,mr;

            int e = v2e[oriV][0];
            if(ev_begin(e)[1]==oriV){ml = m1[e];mr = m2[e];}
            else {ml = m2[e];mr = m1[e];}

            auto &p_ve = in_v2e[inV];
            auto &p_vv = in_v2v[inV];
            for(int j=0;j<2;++j)if(!modifiedV[p_vv[j]]){
                int ine = p_ve[j];
                int nextv = p_vv[j];
                int currv = inV;
                int inml,inmr;
                auto p_ev = inCtrE.data()+ine*2;
                auto p_em = inCtrEMat.data()+ine*2;
                if(p_ev[0]==inV){inml = p_em[0];inmr = p_em[1];}
                else {inml = p_em[1];inmr = p_em[0];}
                bool isleftturn = (inml == 0);

                while(true){
                    addE.push_back(currv);addE.push_back(nextv);
                    if(isleftturn){addEMat.push_back(ml);addEMat.push_back(outsideM);}
                    else {addEMat.push_back(outsideM);addEMat.push_back(mr);}
                    if(modifiedV[nextv])break;
                    modifiedV[nextv] = true;
                    int tmpv = nextv;
                    auto &pvvnext = in_v2v[nextv];
                    nextv = pvvnext[0]==currv?pvvnext[1]:pvvnext[0];
                    currv = tmpv;

                }

            }
        }

    }
    int addInd = n_vertices;
    for(int i=0;i<nv;++i)if(newAddInd[i]==-1){
        newAddInd[i] = addInd++;
        auto p_v = inCtrV.data()+i*3;
        for(int j=0;j<3;++j)vertices.push_back(p_v[j]);
    }
    for(auto &a:addE)a = newAddInd[a];
    edges.insert(edges.end(),addE.begin(),addE.end());
    for(int i=0;i<addEMat.size()/2;++i){
        m1.push_back(addEMat[i*2]);
        m2.push_back(addEMat[i*2+1]);
    }
    CalSelfProperty();
    cout<<"aa"<<endl;


}
void Contour::LableTrianglesWithWindingNumber(vector<double>&testFaceCenter, vector<int>&outfaceMat){


    if(n_vertices == 0 || n_edges == 0){
        outfaceMat.clear();
        outfaceMat.resize(testFaceCenter.size()/3,numeric_limits<int>::max());
        return;
    }
    vector<vector<int>>edgesloops;
    vector<vector<int>>verticesloops;
    vector<int>loops2mat;

    vector<vector<int>>lable2edges;

    vector<vector<uint>>v2v(n_vertices);
    vector<vector<uint>>v2e(n_vertices);
    vector<vector<uint>>e2e(n_edges);

    //map<pair<int,int>,int>edgesmap;
    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        v2v[p_ev[0]].push_back(p_ev[1]);
        v2v[p_ev[1]].push_back(p_ev[0]);

        v2e[p_ev[0]].push_back(i);
        v2e[p_ev[1]].push_back(i);

        //edgesmap[make_pair(p_ev[0],p_ev[1])] = i;
        //edgesmap[make_pair(p_ev[1],p_ev[0])] = i;
    }
    for(auto &v2ee : v2e){
        for(auto a:v2ee)for(auto b:v2ee)if(a!=b)e2e[a].push_back(b);
    }

    set<int>labelset;vector<int>labelvec;
    map<int,int>labelmap;
    for(int i=0;i<n_edges;++i){
        labelset.insert(m1[i]);
        labelset.insert(m2[i]);
    }
    for(auto a:labelset)labelvec.push_back(a);
    for(int i=0;i<labelvec.size();++i)labelmap[labelvec[i]] = i;
    lable2edges.resize(labelvec.size());
    for(int i=0;i<n_edges;++i){
        lable2edges[labelmap[m1[i]]].push_back(i);
        lable2edges[labelmap[m2[i]]].push_back(i);
    }


    for(int i=0;i<lable2edges.size();++i){

        int mat = i;
        vector<int>&ele  = lable2edges[i];

        UnionFind unifind;
        unifind.SetElements(ele);
        vector<bool>chosenEdges( n_edges, false);

        for(auto b:ele)chosenEdges[b] = true;

        for(auto b:ele)for(auto c:e2e[b])if(chosenEdges[c])unifind.Union(b,c);

        vector<vector<int>>tmp_edgesloops;
        unifind.ExtractComponents(tmp_edgesloops);
        for(auto &b:tmp_edgesloops){
            edgesloops.push_back(b);
            loops2mat.push_back(mat);
        }

    }

    verticesloops.resize(loops2mat.size());
    vector<int>countVn(n_vertices,0);
    for(int i=0;i<loops2mat.size();++i){
        vector<int>&v2l  = verticesloops[i];
        vector<int>&ele  = edgesloops[i];
        vector<bool>chosenEdges( n_edges, false);
        vector<bool>chosenVs( n_vertices, false);
        for(auto b:ele)chosenEdges[b] = true;
        for(auto b:ele)for(int j=0;j<2;++j)chosenVs[edges[b*2+j]] = true;

        int pickE = 0;
        int curV,preV;
        if(m2[ele[pickE]] == labelvec[loops2mat[i]]){curV = edges[ele[pickE]*2+1];preV = edges[ele[pickE]*2];}
        else {curV = edges[ele[pickE]*2];preV = edges[ele[pickE]*2+1];}

        int ne = ele.size();
        for(int j=0;j<ne;++j){
            countVn[curV]++;
            v2l.push_back(curV);
            for(auto b:v2v[curV])if(chosenVs[b] && b!=preV){
                preV = curV;
                curV = b;
                break;
            }
        }
        //        ele.clear();
        //        int curE = ele[0],preE = ele[0];
        //        for(int i=0;i<ne;++i){
        //            ele.push_back(curE);
        //            for(auto b:e2e[curE])if(chosenEdges[b] && b!=preE){
        //                preE = curE;
        //                curE = b;
        //                break;
        //            }
        //        }

    }
    //for(int i=0;i<n_vertices;++i)assert(countVn[i]==v2e[i].size());

    //    for(int i=0;i<verticesloops.size();++i){
    //        auto &a=verticesloops[i];
    //        int nv = a.size();
    //        int curV = a[1],preV = a[0];
    //        int gmat = labelvec[loops2mat[i]];
    //        int mat,matinvese;
    //        int kkk = 0;
    //        for(int j=1;j<nv;++j){
    //            int e = edgesmap[make_pair(curV,preV)];

    //            if(edges[e*2]==preV){mat = m2[e];matinvese = m1[e];}
    //            else {mat = m1[e];matinvese = m2[e];}
    //            //assert(gmat==mat);
    //            if(gmat!=mat && gmat==matinvese){
    //                ++kkk;
    //            }
    //            preV = curV;
    //            curV = a[j+1];

    //        }
    //        curV = a[0],preV = a[nv-1];
    //        int e = edgesmap[make_pair(curV,preV)];
    //        if(edges[e*2]==preV){mat = m2[e];matinvese = m1[e];}
    //        else {mat = m1[e];matinvese = m2[e];}
    //        //assert(gmat==mat);
    //        if(gmat!=mat && gmat==matinvese){
    //            ++kkk;
    //        }
    //        cout<<"kkk: "<<kkk<<endl;

    //    }

    //unordered_map<int,double>windingnumber;
    auto calculateWindingnumber = [this](double *p_vpos,vector<int>&vseq){
        int nv = vseq.size();
        vector<double>rays(nv*3);
        auto p_rayd = rays.data();
        for(int i=0;i<nv;++i){
            auto p_rays = p_rayd+i*3;
            minusVec(v_begin(vseq[i]),p_vpos,p_rays);
            normalize(p_rays);
        }
        vector<double>angles(nv);
        double cs[3];
        double isreverse = -1;

        for(int i=0;i<nv-1;++i){
            angles[i] = angleNor(p_rayd+i*3,p_rayd+i*3+3);
            cross(plane_para,p_rayd+i*3,cs);
            if(dot(cs,p_rayd+i*3+3)*isreverse<0)angles[i] = -angles[i];
        }
        angles[nv-1] = angleNor(p_rayd+(nv-1)*3,p_rayd);
        cross(plane_para,p_rayd+(nv-1)*3,cs);
        if(dot(cs,p_rayd)*isreverse<0)angles[nv-1] = -angles[nv-1];

        double windingnumber = 0;
        for(auto a:angles)windingnumber+=a;
        windingnumber/=(2*my_PI);
        return windingnumber;


    };
    auto calculateWindingnumberBatch = [&calculateWindingnumber,this,&loops2mat,&verticesloops](double *p_vpos,vector<double>&outwindingnumber){

        outwindingnumber.resize(loops2mat.size());
        for(int i=0;i<loops2mat.size();++i){
            outwindingnumber[i] = calculateWindingnumber(p_vpos,verticesloops[i]);
        }

    };

    int ntcv = testFaceCenter.size()/3;
    outfaceMat.resize(ntcv);
    for(int i=0;i<ntcv;++i){
        auto p_tcv = testFaceCenter.data()+i*3;
        vector<double>outwindingnumber;
        calculateWindingnumberBatch(p_tcv,outwindingnumber);
        //cout<<"ori: ";for(auto a:outwindingnumber)cout<<a<<' ';cout<<endl;
        vector<double>windingnumberEveryMat(labelvec.size(),0);for(int j=0;j<labelvec.size();++j)windingnumberEveryMat[j] = 0;
        for(int j=0;j<loops2mat.size();++j)windingnumberEveryMat[loops2mat[j]] += outwindingnumber[j];


        windingnumberEveryMat[labelvec.size()-1]+=1;
        //for(auto a:windingnumberEveryMat)cout<<a<<' ';cout<<endl;

        for(auto &a:windingnumberEveryMat)a = fabs(1-a);

        outfaceMat[i] = labelvec[min_element(windingnumberEveryMat.begin(),windingnumberEveryMat.end())-windingnumberEveryMat.begin()];

    }

    cout<<endl;






}

int Contour::getnumberofloops(){

    vector<vector<double>>verticesSeqs;
    vector<int>loops2mats;
    vector<vector<int>>verticesloops;

    ComputeMat2VerticesLoop(verticesloops,verticesSeqs,loops2mats);
    return loops2mats.size();
}

void Contour::ComputeMat2VerticesLoop(vector<vector<int>>&verticesloops,vector<vector<double>>&verticesSeqs, vector<int>&loops2mats){


    if(n_vertices == 0 || n_edges == 0){
        verticesSeqs.clear();
        loops2mats.clear();
        return;
    }
    vector<vector<int>>edgesloops;

    vector<int>loops2mat;

    vector<vector<int>>lable2edges;

    vector<vector<uint>>v2v(n_vertices);
    vector<vector<uint>>v2e(n_vertices);
    vector<vector<uint>>e2e(n_edges);

    //map<pair<int,int>,int>edgesmap;
    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        v2v[p_ev[0]].push_back(p_ev[1]);
        v2v[p_ev[1]].push_back(p_ev[0]);

        v2e[p_ev[0]].push_back(i);
        v2e[p_ev[1]].push_back(i);

        //edgesmap[make_pair(p_ev[0],p_ev[1])] = i;
        //edgesmap[make_pair(p_ev[1],p_ev[0])] = i;
    }
    for(auto &v2ee : v2e){
        for(auto a:v2ee)for(auto b:v2ee)if(a!=b)e2e[a].push_back(b);
    }

    set<int>labelset;vector<int>labelvec;
    map<int,int>labelmap;
    for(int i=0;i<n_edges;++i){
        labelset.insert(m1[i]);
        labelset.insert(m2[i]);
    }
    for(auto a:labelset)labelvec.push_back(a);
    for(int i=0;i<labelvec.size();++i)labelmap[labelvec[i]] = i;
    lable2edges.resize(labelvec.size());
    for(int i=0;i<n_edges;++i){
        lable2edges[labelmap[m1[i]]].push_back(i);
        lable2edges[labelmap[m2[i]]].push_back(i);
    }


    for(int i=0;i<lable2edges.size();++i){

        int mat = i;
        vector<int>&ele  = lable2edges[i];

        UnionFind unifind;
        unifind.SetElements(ele);
        vector<bool>chosenEdges( n_edges, false);

        for(auto b:ele)chosenEdges[b] = true;

        for(auto b:ele)for(auto c:e2e[b])if(chosenEdges[c])unifind.Union(b,c);

        vector<vector<int>>tmp_edgesloops;
        unifind.ExtractComponents(tmp_edgesloops);
        for(auto &b:tmp_edgesloops){
            edgesloops.push_back(b);
            loops2mat.push_back(mat);
        }

    }

    verticesloops.resize(loops2mat.size());
    vector<int>countVn(n_vertices,0);
    for(int i=0;i<loops2mat.size();++i){
        vector<int>&v2l  = verticesloops[i];
        vector<int>&ele  = edgesloops[i];
        vector<bool>chosenEdges( n_edges, false);
        vector<bool>chosenVs( n_vertices, false);
        for(auto b:ele)chosenEdges[b] = true;
        for(auto b:ele)for(int j=0;j<2;++j)chosenVs[edges[b*2+j]] = true;

        int pickE = 0;
        int curV,preV;
        if(m2[ele[pickE]] == labelvec[loops2mat[i]]){curV = edges[ele[pickE]*2+1];preV = edges[ele[pickE]*2];}
        else {curV = edges[ele[pickE]*2];preV = edges[ele[pickE]*2+1];}

        int ne = ele.size();
        for(int j=0;j<ne;++j){
            countVn[curV]++;
            v2l.push_back(curV);
            for(auto b:v2v[curV])if(chosenVs[b] && b!=preV){
                preV = curV;
                curV = b;
                break;
            }
        }
    }


    loops2mats = loops2mat;
    for(auto &a:loops2mats)a = labelvec[a];
    verticesSeqs.clear();
    verticesSeqs.resize(loops2mats.size());
    for(int i=0;i<loops2mats.size();++i){
        auto &vsq = verticesSeqs[i];
        vsq.clear();
        auto &vlp = verticesloops[i];
        for(int j=0;j<vlp.size();++j){
            auto p_v = v_begin(vlp[j]);
            for(int k=0;k<3;++k)vsq.push_back(p_v[k]);
        }
    }


}



void Contour::MergeMaterial(vector<int>&mappingMat){

    if(1){
        cout<<n_materials<<endl;
        vector<set<int>>contactness(n_materials);
        for(int i=0;i<n_edges;++i){
            contactness[m1[i]].insert(m2[i]);
            contactness[m2[i]].insert(m1[i]);
        }

        for(int i=0;i<n_materials;++i){
            cout<<i<<": ";
            for(auto a:contactness[i])cout<<a<<' ';cout<<endl;
        }

        //exit(11);
    }

    for(auto &a:m1)a = mappingMat[a];
    for(auto &a:m2)a = mappingMat[a];

    vector<bool>restEdge(n_edges,false);

    for(int i=0;i<n_edges;++i){
        if(m1[i]!=m2[i])restEdge[i]=true;
    }

    vector<int>newVer(n_vertices,-1);
    for(int i=0;i<n_edges;++i)if(restEdge[i]){
        auto p_ev= ev_begin(i);
        for(int j=0;j<2;++j)newVer[p_ev[j]] = 0;

    }
    int neInd = 0;
    for(auto &a:newVer)if(a==0)a=neInd++;

    vector<double>verticesPos;
    vector<uint>edges2vertices;

    for(int i=0;i<n_vertices;++i)if(newVer[i]!=-1){
        auto p_v = v_begin(i);
        for(int j=0;j<3;++j)verticesPos.push_back(p_v[j]);
    }
    for(int i=0;i<n_edges;++i)if(restEdge[i]){
        auto p_ev= ev_begin(i);
        for(int j=0;j<2;++j)edges2vertices.push_back(newVer[p_ev[j]]);
    }

    vector<int>_groupN1,_groupN2;
    for(int i=0;i<n_edges;++i)if(restEdge[i]){
        _groupN1.push_back(m1[i]);
        _groupN2.push_back(m2[i]);
    }


    edges = edges2vertices;
    vertices = verticesPos;

    m1 = _groupN1;
    m2 = _groupN2;

    n_vertices = vertices.size()/3;
    n_edges = edges.size()/2;

    CurveSmoothing(100);



}
void Contour::AddContour(Contour &mc){

    int offset = vertices.size()/3;
    for(auto a:mc.vertices)vertices.push_back(a);
    for(auto a:mc.edges)edges.push_back(a+offset);
    for(auto a:mc.m1)m1.push_back(a);
    for(auto a:mc.m2)m2.push_back(a);
    CalSelfProperty();

}
void Contour::InversePlaneNormal(){

    for(int i=0;i<4;++i)plane_para[i] = -plane_para[i];
    swap(m1,m2);

}
void Contour::RotatingandShifting(double r_angle, vector<double>&shiftbeforeR,vector<double>&raxis, vector<double>&shiftv){
/*
    vector<Eigen::Vector3d>ori(n_vertices);
    vector<Eigen::Vector3d>tranformed(n_vertices);

    for(int j = 0;j<n_vertices;++j){
        auto p_spv = v_begin(j);
        for(int k=0;k<3;++k)ori[j](k) = p_spv[k] + shiftbeforeR[k];
    }

    Eigen::Vector3d ra;
    for(int k=0;k<3;++k)ra(k) = raxis[k];


    Eigen::AngleAxisd t(r_angle,ra);


    for(int j = 0;j<n_vertices;++j){
        tranformed[j] = t*ori[j];
    }
    for(int j = 0;j<n_vertices;++j){
        auto p_spv = v_begin(j);
        for(int k=0;k<3;++k)p_spv[k] = tranformed[j](k) - shiftbeforeR[k] + shiftv[k];
    }
    Eigen::Vector3d vnormal;
    for(int i=0;i<3;++i)vnormal(i) = plane_para[i];
    vnormal = t*vnormal;
    for(int i=0;i<3;++i)plane_para[i] = vnormal(i);


    //plane_para[3]-=dot(shiftv.data(),plane_para);
    plane_para[3]=-dot(vertices.data(),plane_para);
*/

}

int CrossSections::SpecialMerging(string infilename, vector<double>&points, vector<uint> &edges){
    ifstream reader(infilename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Contour file " << infilename << endl;
        return -1;
    }
    int n_c;
    reader>>n_c;

    CSs.resize(n_c);


    for(int i=0;i<n_c;++i)CSs[i].ReadContour(reader);

    vector<vector<int>>mappingMat(n_c);
    if(0){
        for(int i=0;i<n_c;++i){
            mappingMat[i].resize(18,0);
            //for(int j=10;j<17;++j)mappingMat[i][j] = 1;
            for(int j=0;j<18;++j)mappingMat[i][j] = j;
        }
        {
            for(int i=4;i<18;++i)mappingMat[0][i] = 3;
            mappingMat[0][15] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[1][i] = 3;
            mappingMat[1][8] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[2][i] = 3;
            mappingMat[2][8] = 8;
            mappingMat[2][7] = 7;
        }

        {
            //for(int i=4;i<18;++i)mappingMat[2][i] = 3;
            mappingMat[3][8] = 0;
            mappingMat[3][10] = 0;
            mappingMat[3][11] = 0;
            mappingMat[3][16] = 6;
        }
        {

            mappingMat[4][12] = 11;
            mappingMat[4][17] = 11;
            mappingMat[4][2] = 6;
            mappingMat[4][10] = 1;

        }
    }

    if(1){
        for(int i=0;i<n_c;++i){
            mappingMat[i].resize(18,0);
            //for(int j=10;j<17;++j)mappingMat[i][j] = 1;
            for(int j=0;j<18;++j)mappingMat[i][j] = j;
        }
        {
            for(int i=4;i<18;++i)mappingMat[0][i] = 1;
            mappingMat[0][13] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[1][i] = 3;
            mappingMat[1][15] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[2][i] = 3;
            mappingMat[2][8] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[3][i] = 3;
            mappingMat[3][8] = 8;
            mappingMat[3][7] = 7;
        }

        {
            //for(int i=4;i<18;++i)mappingMat[2][i] = 3;
            mappingMat[4][8] = 0;
            mappingMat[4][10] = 0;
            mappingMat[4][11] = 0;
            mappingMat[4][16] = 6;
        }
        {

            mappingMat[5][12] = 11;
            mappingMat[5][17] = 11;
            mappingMat[5][2] = 6;
            mappingMat[5][10] = 1;

        }
        {

            mappingMat[6][12] = 11;
            mappingMat[6][17] = 11;
            mappingMat[6][2] = 6;
            mappingMat[6][10] = 1;

        }
    }

    int testCC = 6;
    for(int i=0;i<n_c;++i){
        //if(i!=testCC)continue;
        CSs[i].MergeMaterial(mappingMat[i]);
    }
    int n_materials = CSs[0].n_materials;

    if(0){
        points = CSs[testCC].vertices;
        edges = CSs[testCC].edges;
    }else Stackup(points,edges);
    return 0;


}
void CrossSections::ReadCrossSectionsfromImages(vector<vector<vector<int>>>&im,double zdis){

    n_CrossSec = im.size();

    CSs.resize(n_CrossSec);
    for(int i=0;i<n_CrossSec;++i){

        CSs[i].BuildContourFromImage(im[i],(i+1)*zdis);
    }




}

bool CrossSections::ReadCrossSections(string filename){
    ifstream fin(filename.data());
    if(fin.fail()){
        cout<<"Fail to open input file: "<<filename<<endl;
        return false;
    }

    fin>>n_CrossSec;

    //n_CrossSec = 5;


    CSs.resize(n_CrossSec);
    for(int i=0;i<n_CrossSec;++i){
        CSs[i].ReadContour(fin);
    }
    CalSelfProperty();

    return true;
}

bool CrossSections::ReadCrossSections(vector<vector<double>>&planeinfo, vector<vector<double>>&plane2vertices,
                            vector<vector<uint> > &plane2edges, vector<vector<int>>&plane2edgesMat){



    n_CrossSec = planeinfo.size();

    CSs.resize(n_CrossSec);
    for(int i=0;i<n_CrossSec;++i){
        CSs[i].ImportContour(planeinfo[i].data(),plane2vertices[i],plane2edges[i],plane2edgesMat[i]);
    }
    CalSelfProperty();

    return true;


}



void CrossSections::ShiftCrossSections(vector<vector<double>>&trlvec){
    vector<double>rota({0,0,1});
    vector<double>shiftbeforeR({0,0,0});

    for(int i=0;i<CSs.size();++i){

        CSs[i].RotatingandShifting(0,shiftbeforeR,rota,trlvec[i]);
    }


}
void CrossSections::ManipulateCrossSections(){

    //RotatingandShifting

    vector<double>rota({0,0,1});
    vector<double>shiftbeforeR({-500,-500,0});
    vector<double>shiftv({0,0,0});


    if(0){
        n_CrossSec = 5;
        CSs.resize(n_CrossSec);
        {
            CSs[2] = CSs[1];
            shiftv[2] = 100;shiftv[0] = 0;shiftv[1] = 0;
            CSs[2].RotatingandShifting(20/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        {
            CSs[3] = CSs[1];
            shiftv[2] = 200;shiftv[0] = 30;shiftv[1] = -60;
            CSs[3].RotatingandShifting(20/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        {
            CSs[4] = CSs[0];
            shiftv[2] = 400;shiftv[0] = 60;shiftv[1] = -30;
            CSs[4].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv);
        }
    }

    if(1){
        n_CrossSec = 5;
        CSs.resize(n_CrossSec);
        {
            CSs[1] = CSs[0];
            shiftv[2] = 100;shiftv[0] = 75;shiftv[1] = 70;
            CSs[1].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv);
        }
        {
            CSs[2] = CSs[0];
            shiftv[2] = 200;shiftv[0] = 00;shiftv[1] = 00;
            CSs[2].RotatingandShifting(20/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        {
            CSs[3] = CSs[0];
            shiftv[2] = 300;shiftv[0] = -80;shiftv[1] = -90;
            CSs[3].RotatingandShifting(-20/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        {
            CSs[4] = CSs[0];
            shiftv[2] = 400;shiftv[0] = -10;shiftv[1] = -10;
            CSs[4].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        //        {
        //            CSs[5] = CSs[0];
        //            shiftv[2] = 500;shiftv[0] = 0;shiftv[1] = -0;
        //            CSs[5].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv);
        //        }
    }

    if(0){
        n_CrossSec = 4;
        CSs.resize(n_CrossSec);
        vector<double>shiftbeforeRBB({(-800-1000)/2,-00,-00});
        vector<double>rotaRR({0,1,0});
        vector<double>shiftv11({0,0,-100});
        CSs[0].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv11);

        {
            Contour tmpcc = CSs[0];
            shiftv[2] = 0;shiftv[0] =800;shiftv[1] = 0;
            tmpcc.RotatingandShifting(-180/180.*3.1415,shiftbeforeR,rotaRR,shiftv);
            tmpcc.InversePlaneNormal();
            CSs[0].AddContour(tmpcc);
        }

        if(n_CrossSec>1){
            CSs[1] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[1].RotatingandShifting(90/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }

        if(n_CrossSec>2){
            CSs[2] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[2].RotatingandShifting(45/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }

        if(n_CrossSec>3){
            CSs[3] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[3].RotatingandShifting(135/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }


    }
    if(0){
        n_CrossSec = 3;
        CSs.resize(n_CrossSec);
        int middistance = 1500;
        vector<double>shiftbeforeR({-500,-500,0});
        vector<double>shiftbeforeRBB({0,-00,-00});
        shiftbeforeRBB[0] = (-middistance-1000)/2;
        vector<double>rotaRR({0,1,0});
        vector<double>shiftv11({0,0,-200});
        CSs[0].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv11);

        {
            Contour tmpcc = CSs[0];
            shiftv[2] = 0;shiftv[0] =middistance;shiftv[1] = 0;
            tmpcc.RotatingandShifting(-180/180.*3.1415,shiftbeforeR,rotaRR,shiftv);
            tmpcc.InversePlaneNormal();
            CSs[0].AddContour(tmpcc);
        }

        if(n_CrossSec>1){
            CSs[1] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[1].RotatingandShifting(60/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }

        if(n_CrossSec>2){
            CSs[2] = CSs[0];
            shiftv[2] = 50;shiftv[0] = 50;shiftv[1] = 0;
            CSs[2].RotatingandShifting(120/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }

        if(n_CrossSec>3){
            CSs[3] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[3].RotatingandShifting(135/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }


    }

    if(0){
        vector<double>shiftbeforeR({-0,-0,0});
        vector<double>shiftv11({10,0,0});
        for(int i=0;i<CSs.size();++i){
            shiftv11[0] = i*20;
            CSs[i].RotatingandShifting(0,shiftbeforeR,rota,shiftv11);
        }

    }
}


void CrossSections::CalSelfProperty(){
    for(auto &a:CSs)a.CalSelfProperty();
    width = height = -1;

    n_CrossSec = CSs.size();
    leftcorner.clear();
    rightcorner.clear();
    midpoint.clear();
    leftcorner.resize(3,1000000);
    rightcorner.resize(3,-1000000);
    midpoint.resize(3,0);
    for(int i=0;i<n_CrossSec;i++)for(int j=0;j<3;++j)leftcorner[j] = min(leftcorner[j],CSs[i].leftcorner[j]);
    for(int i=0;i<n_CrossSec;i++)for(int j=0;j<3;++j)rightcorner[j] = max(rightcorner[j],CSs[i].rightcorner[j]);
    for(int i=0;i<n_CrossSec;i++)for(int j=0;j<3;++j)midpoint[j] += CSs[i].midpoint[j];
    for(int j=0;j<3;++j)midpoint[j] /= n_CrossSec;


    total_nver = 0;
    total_nedges = 0;
    for(auto &a:CSs)total_nedges+=a.n_edges;
    acc_nver.resize(n_CrossSec+1,0);
    for(int i=0;i<n_CrossSec;i++)acc_nver[i+1] = CSs[i].n_vertices+acc_nver[i];
    total_nver = acc_nver[n_CrossSec];

    set<int>matts;
    for(auto &a:CSs)for(auto b:a.m1)matts.insert(b);
    for(auto &a:CSs)for(auto b:a.m2)matts.insert(b);
    n_materials = matts.size();

    planes_info.clear();
    for(auto &a:CSs)planes_info.push_back(a.plane_info);


}
void CrossSections::Stackup(vector<double>&points, vector<uint> &edges){

    CalSelfProperty();

    points.clear();

    for(int i=0;i<n_CrossSec;i++)for(int j=0;j<CSs[i].n_vertices*3;++j)points.push_back( CSs[i].vertices[j]);
    edges.clear();
    for(int i=0;i<n_CrossSec;i++)for(int j=0;j<CSs[i].n_edges*2;++j){
        edges.push_back( CSs[i].edges[j]+acc_nver[i]);
    }

}

void ReadCrossSection(string filename){





}

void CrossSections::MapMaterials(){

    map<int,int>mapping;
    mapping[11]=0;
    mapping[15] = 2;
    mapping[17] = 1;

    for(auto &a:CSs){
        for(auto &b:a.m1)b = mapping[b];
        for(auto &b:a.m2)b = mapping[b];
    }


}
int CrossSections::WriteCrossSections(string filename){
    string outCname = filename + string("mapping.contour");
    ofstream outer(outCname.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }

    outer<<std::setprecision(12);
    outer<<CSs.size()<<endl;
    for(auto &b:CSs)b.WriteContour(outer);
    outer.close();
    if(0)for(int i=0;i<CSs.size();++i){
        string outCname = filename + string("mapping")+to_string(i)+ string(".contour");
        ofstream outer(outCname.data(), ofstream::out);
        if (!outer.good()) {
            cout << "Can not open the Contour file " << outCname << endl;
            return -1;
        }

        outer<<1<<endl;
        CSs[i].WriteContour(outer);
        outer.close();
    }

    return 0;

}

int CrossSections::WriteCrossSections(string filename,vector<int>&pickedCs){


    string outCname = filename;
    ofstream outer(outCname.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }

    outer<<std::setprecision(12);outer<<std::fixed;
    outer<<pickedCs.size()<<endl;
    for(auto b:pickedCs)CSs[b].WriteContour(outer);
    outer.close();

}


int CrossSections::WriteCrossSectionsCtrFormat(string filename){
    string outCname = filename + string("mapping.ctr");
    ofstream outer(outCname.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }

    outer<<std::setprecision(12);
    int totalloop = 0;
    for(auto &b:CSs)totalloop+=b.getnumberofloops();
    outer<<totalloop/2<<endl;
    for(auto &b:CSs)b.WriteContourCtrFormat(outer);
    outer.close();

    return 0;



}


int CrossSections::WriteCrossSectionsCSLFormat(string filename){
    string outCname = filename + string(".csl");
    ofstream outer(outCname.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }

    CalSelfProperty();

    int ncs = CSs.size();
    //int ncs = 5;

    outer<<std::setprecision(12);outer<<std::fixed;

    outer<<"CSLC"<<endl;
    outer<<ncs<<' '<<n_materials<<endl;
    outer<<endl;
    for(int i=0;i<ncs;++i){
        outer<<i+1<<' ';
        CSs[i].WriteContourCSLFormat(outer);
        outer<<endl;
    }


    return 0;




}


void Contour::WriteContourToFCM(ofstream &out){


    vector<vector<double>>verticesSeqs;
    vector<int>loops2mats;
    vector<vector<int>>verticesloops;

    ComputeMat2VerticesLoop(verticesloops,verticesSeqs,loops2mats);

    //    out<<"{0,0,0}";
    //    return;
    // {{v, e},{v, e}}
    //loop1 = {loop1, eloop1, 1};
    //{loop1, loop2, loop3}
    //
    out<<"{";
    for(int i=0;i<loops2mats.size();++i){
        out<<"{";
        //v
        //out<<"{{0,0,0}},";
        out<<"{";
        auto p_v = verticesSeqs[i].data();
        int nv = verticesSeqs[i].size()/3;
        int j=0;
        for(j=0;j<nv-1;++j){
            out<<"{";
            for(int k=0;k<2;++k)out<<p_v[j*3+k]<<", ";
            out<<p_v[j*3+2];
            out<<"}, ";
        }
        out<<"{";
        for(int k=0;k<2;++k)out<<p_v[j*3+k]<<", ";
        out<<p_v[j*3+2];
        out<<"}";
        out<<"}, ";

        //e
        //out<<"{{0,1}},";

        out<<"{";
        for(int k=1;k<nv;++k){
            out<<"{";
            out<<k<<", "<<k+1;
            out<<"}, ";
        }
        out<<"{";out<<nv<<", "<<1;out<<"}";
        out<<"},";

        out<<loops2mats[i]+1;

        if(i!=loops2mats.size()-1)out<<"}, ";
        else out<<"}";
    }
    out<<"}";


}

int CrossSections::WriteCrossSectionsToFCM(string filename,int nMat){

    string outCname = filename;
    ofstream out(outCname.data(), ofstream::out);
    if (!out.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }
    cout << "write:  " << outCname << endl;

    out<<std::setprecision(12);
    out<<std::fixed;
    out<<"{";
    //loops
    out<<"{";
    for(int i=0;i<CSs.size();++i){
        CSs[i].WriteContourToFCM(out);
        if(i!=CSs.size()-1)out<<",";

    }

    //out<<"},";
    out<<"},";

    //plane paras
    out<<"{";
    for(int i=0;i<CSs.size();++i){
        out<<"{";
        out<<"{0,0,0},";
        //          auto p_cv = CSs[i].v_begin(0);
        //          out<<"{";
        //          out<<p_cv[0]<<","<<p_cv[1]<<","<<p_cv[2];
        //          out<<"},";
        out<<"{";
        for(int j=0;j<3;++j)out<<CSs[i].plane_para[j]<<",";
        out<<CSs[i].plane_para[3];
        out<<"},";
        vector<double>tmp(3,1);
        for(int j=0;j<3;++j)if(fabs(CSs[i].plane_para[j])>1e-1){
            int v1,v2;
            if(j==0){v1=1;v2=2;}else if(j==1){v1=0;v2=2;}else if(j==2){v1=0;v2=1;}
            tmp[j] = (0-CSs[i].plane_para[v1]*tmp[v1]-CSs[i].plane_para[v2]*tmp[v2])/CSs[i].plane_para[j];
            break;
        }
        normalize(tmp.data());
        out<<"{"<<tmp[0]<<","<<tmp[1]<<","<<tmp[2]<<"}";

        if(i!=CSs.size()-1) out<<"},";
        else out<<"}";
    }
    out<<"},";

    out<<nMat<<"}";
    out.close();


    cout<<"write fin"<<endl;
    return 0;

}

void Contour::WriteContourToFCMtCtr(ofstream &out){





    int nV = vertices.size()/3;
    out<<"{";
    out<<"{";

    //v
    out<<"{";
    for(int i=0;i<nV;++i){
        auto p_v = v_begin(i);
        out<<"{";
        out<<p_v[0]<<","<<p_v[1]<<","<<p_v[2];
        if(i!=nV-1)out<<"},";else out<<"}";
    }
    out<<"},";

    //e
    int ne = edges.size()/2;
    out<<"{";
    for(int i=0;i<ne;++i){
        auto p_e = ev_begin(i);
        out<<"{";
        out<<p_e[0]+1<<","<<p_e[1]+1;
        if(i!=ne-1)out<<"},";else out<<"}";
    }
    out<<"}";


    out<<"}";
    out<<"}";

}

int CrossSections::WriteCrossSectionsToFCMtCtr(string filename){

    string outCname = filename;
    ofstream out(outCname.data(), ofstream::out);
    if (!out.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }
    cout << "write:  " << outCname << endl;
    out<<std::setprecision(12);
    out<<std::fixed;
    out<<"{";
    //loops
    out<<"{";
    for(int i=0;i<CSs.size();++i){
        CSs[i].WriteContourToFCMtCtr(out);
        if(i!=CSs.size()-1)out<<",";

    }

    //out<<"},";
    out<<"},";

    //plane paras
    out<<"{";
    for(int i=0;i<CSs.size();++i){
        out<<"{";
        //          auto p_cv = CSs[i].v_begin(0);
        //          out<<"{";
        //          out<<p_cv[0]<<","<<p_cv[1]<<","<<p_cv[2];
        //          out<<"},";
        out<<"{0,0,0},";

        out<<"{";
        for(int j=0;j<3;++j)out<<CSs[i].plane_para[j]<<",";
        out<<CSs[i].plane_para[3];
        out<<"},";
        vector<double>tmp(3,1);
        for(int j=0;j<3;++j)if(fabs(CSs[i].plane_para[j])>1e-1){
            int v1,v2;
            if(j==0){v1=1;v2=2;}else if(j==1){v1=0;v2=2;}else if(j==2){v1=0;v2=1;}
            tmp[j] = (0-CSs[i].plane_para[v1]*tmp[v1]-CSs[i].plane_para[v2]*tmp[v2])/CSs[i].plane_para[j];
            break;
        }
        normalize(tmp.data());
        out<<"{"<<tmp[0]<<","<<tmp[1]<<","<<tmp[2]<<"}";

        if(i!=CSs.size()-1) out<<"},";
        else out<<"}";
    }
    out<<"}";

    out<<"}";

    out.close();


    cout<<"write fin"<<endl;
    return 0;

}


int CrossSections::WriteCrossSectionsObjFormat(string filename){


    filename = filename + ".obj";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output Obj file " << filename << endl;
        //return false;
    }


    vector<double>vertices;vector<uint>edge2vertices;
    Stackup(vertices,edge2vertices);

    int n_vertices = vertices.size()/3;

    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }


    {
        int n_edges = edge2vertices.size()/2;
        for(int i=0;i<n_edges;++i){
            auto p_ev = edge2vertices.data()+i*2;
            outer << "l " << p_ev[0]+1<< " "<< p_ev[1]+1 << endl;
        }
    }


    outer.close();
    cout<<"saving finish: "<<filename<<endl;


}




void CrossSections::SubSetPlanes(vector<int>&reserveplanes){



    vector<Contour>newCSs;

    for(auto a:reserveplanes)newCSs.push_back(CSs[a]);


    swap(CSs , newCSs);

    CalSelfProperty();


}







void CrossSections::naiveConstruction(){
/*
    cout<<"naiveConstruction"<<endl;
    //--------------------- Construct the frame ---------------------
    // <1> frameVs.
    // Two addtional layers are added on the top and bottom. The spaceing of each
    // addtional layer to it adjacent layer is set to be the minimum spacing among
    // all the spacings.
    // |spacings| = _nCrossSection-1; |extSpacings| = _nCrossSection+2 = nLayer.
    vector<float> spacings(n_CrossSec-1, 100.0f);

    vector<vector<double> > _frameVs;
    vector<vector<int> > _frameFs;
    vector<vector<int>> _frameCells;
    vector<vector<double> > _frameF2PlaneParas;


    CalSelfProperty();


    float bdSpacing = *min_element(spacings.begin(), spacings.end());
    vector<float> extSpacings = spacings;
    // For the second layer (first true layer).
    extSpacings.insert(extSpacings.begin(), bdSpacing);
    // For the first layer.
    extSpacings.insert(extSpacings.begin(), 0.0);
    // For the last layer.
    extSpacings.push_back(bdSpacing);
    int nLayer = n_CrossSec + 2;
    int _nCell = n_CrossSec + 1;
    for (int i = 1; i < extSpacings.size(); ++i) {
        extSpacings[i] += extSpacings[i - 1];
    }
    _frameVs.reserve(nLayer * 4);
    for (int i = 0; i < nLayer; ++i) {
        _frameVs.push_back({0, 0, extSpacings[i]});
        _frameVs.push_back({1.0 * width, 0., extSpacings[i]});
        _frameVs.push_back({1.0 * width, 1.0 * height, extSpacings[i]});
        _frameVs.push_back({0, 1.0 * height, extSpacings[i]});
    }
    // <2> frameFs.
    // frameFs = shared Fs and side Fs.
    // Shared Fs.
    for(auto a:extSpacings)cout<<a<<' ';cout<<endl;
    for (int i = 0; i < nLayer; ++i) {
        _frameFs.push_back({i * 4, i * 4 + 1, i * 4 + 2, i * 4 + 3});
    }
    // Side Fs.
    for (int i = 0; i < _nCell; ++i) {
        int pre = i * 4;
        _frameFs.push_back({pre + 0, pre + 1, pre + 5, pre + 4});
        _frameFs.push_back({pre + 1, pre + 2, pre + 6, pre + 5});
        _frameFs.push_back({pre + 2, pre + 3, pre + 7, pre + 6});
        _frameFs.push_back({pre + 3, pre + 0, pre + 4, pre + 7});
    }
    // <3> frameCells.
    for (int i = 0; i < _nCell; ++i) {
        int pre = i * 4 + nLayer;
        _frameCells.push_back({i, i + 1, pre, pre + 1, pre + 2, pre + 3});
    }
    // <4> frameF2PlaneParas.
    // Four parameters used to define each plane {a,b,c,d} in the form of
    // ax+by+cz+d=0.
    for (int i = 0; i < nLayer; ++i) {
        _frameF2PlaneParas.push_back({0, 0, 1, -1.0f * extSpacings[i]});
    }

    //return;

    // ----------------------- Generate tets ------------------------
    // Find the centroid of each cell for classifying tets in to cells.
    vector<vector<double>> cellCenter(_nCell, vector<double>(3));
    for (int i = 0; i < _nCell; ++i) {
        auto& res = cellCenter[i];
        for (int p = i * 4; p < i * 4 + 8; ++p) {
            std::transform(res.begin(), res.end(), _frameVs[p].begin(), res.begin(),
                           std::plus<double>());
        }
        for(auto &a:res)a/=8;
        //Utility::MultiplyVectorsByValue(0.125f, res);
    }
    // Prepare the data structure for tetgen.
    tetgenio in, out;
    tetgenio::facet* f;
    tetgenio::polygon* p;
    // PLC: pointlist.
    cout<<"PLC: pointlist. "<<_frameVs.size()<<' '<<total_nver<<endl;
    in.numberofpoints = _frameVs.size()+total_nver;
    //in.numberofpoints = _frameVs.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    in.pointmarkerlist = new int[in.numberofpoints];
    for (int i = 0; i < _frameVs.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            in.pointlist[i * 3 + j] = _frameVs[i][j];
        }
        in.pointmarkerlist[i] = -1;
    }
    for(int i=0;i<n_CrossSec;++i){
        int offset = (_frameVs.size()+acc_nver[i])*3;
        for(int j=0;j<CSs[i].vertices.size();++j)in.pointlist[offset+j] = CSs[i].vertices[j];
        offset = offset/3;
        for(int j=0;j<CSs[i].n_vertices;++j)in.pointmarkerlist[offset+j] = -10-i;
    }


    // PLC: facetlist.
    cout<<"PLC: facetlist."<<endl;
    in.numberoffacets = _frameFs.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (int i = 0; i < _frameFs.size(); ++i) {
        f = &in.facetlist[i];
        f->numberofpolygons = 1;
        if(i<1 || i>n_CrossSec)f->numberofpolygons = 1;
        else f->numberofpolygons=1+CSs[i-1].n_edges;

        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = _frameFs[i].size();
        p->vertexlist = new int[p->numberofvertices];
        for (int j = 0; j < _frameFs[i].size(); ++j) {
            p->vertexlist[j] = _frameFs[i][j];
        }
        in.facetmarkerlist[i] = i;

        if(1)if(!(i<1 || i>n_CrossSec)){
            int offset = _frameVs.size()+acc_nver[i-1];
            for(int j=0;j<CSs[i-1].n_edges;++j){
                p = &f->polygonlist[j+1];
                p->numberofvertices = 2;p->vertexlist = new int[p->numberofvertices];
                auto p_ev = CSs[i-1].ev_begin(j);
                for (int k = 0; k < 2; ++k)p->vertexlist[k] = p_ev[k]+ offset;
            }
        }
    }

    // PLC: region.
    cout<<"PLC: region."<<endl;
    in.numberofregions = _frameCells.size();
    in.regionlist = new REAL[in.numberofregions * 5];
    for (int i = 0; i < _frameCells.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            in.regionlist[i * 5 + j] = cellCenter[i][j];
        }
        // Region attribute (marker).
        // Starting from -1. {-1, -2, -3, ...}.
        in.regionlist[i * 5 + 3] = -i - 1;
        // Region volume constraint (not used if specify overall constraint after -a
        // switch).
        in.regionlist[i * 5 + 4] = 100000;
    }


    char tetin[] = "../tetgen1.4.3/c++tetin";
    // in.save_nodes(tetin);
    // in.save_poly(tetin);

    // Tetrahedralize the PLC.
    // Switches are chosen to read PLC (p), do quality mesh generation (q) with a
    // specified quality bound (1.414), apply a maximum volume constraint
    // (a10000) and omit terminal output except errors (Q).
    // -Qpq1.4a10000
    float tetVolLimit = 5000.0f;
    string switchesStr = "Qpq1.414Aa" + to_string(tetVolLimit);
    // char switches[switchesS.size()] = switchesS.c_str();
    // char switches[] = "Qpq1.414Aa10000";
    std::vector<char> tmp(switchesStr.begin(), switchesStr.end());
    tmp.push_back('\0');
    char* switches = &tmp[0];
    // char switches[] = "Qpq1.414Aa10";
    tetrahedralize(switches, &in, &out);


    char tetout[] = "/Users/Research/Geometry/MM/ConvertFolder/tetout";
    out.save_nodes(tetout);
    out.save_elements(tetout);
    out.save_faces(tetout);
*/


}



