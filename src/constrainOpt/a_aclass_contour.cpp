#include "aclass.h"
#include <iostream>
#include <fstream>
#include<math.h>
#include <iomanip>
#include <assert.h>
#include <chrono>
#include "../mymesh/contour.h"
inline void weightedAddVec(const double weight1,const double weight2, const double *orivec1,const double *orivec2,double *desvec){
    for(int i =0;i<3;++i)desvec[i] = weight1*orivec1[i] + weight2*orivec2[i];
}


void Graphs::UnifiedValue(vector<vector<double>>&inval,vector<int>&pickv,int pickplane,vector<vector<double>>&outvalue){

    int npickv = pickv.size();
    outvalue.clear();
    outvalue.resize(pickv.size());

    for(int i=0;i<npickv;++i){
        vector<double>&invalpv = inval[pickv[i]];
        vector<double>vbuf(nMat,0);
        int counc = 0;
        for(int j=0;j<nP;++j){
            bool isvoid = true;
            for(int k=0;k<nMat;++k)if(invalpv[j*nMat+k]!=-1)isvoid=false;

            if(j==pickplane)assert(!isvoid);
            if(!isvoid){
                for(int k=0;k<nMat;++k)vbuf[k]+=invalpv[j*nMat+k];
                counc++;
            }

        }

        for(int k=0;k<nMat;++k)vbuf[k]/=counc;

        outvalue[i] = vbuf;
    }



}


void Graphs::UnifiedValue(vector<vector<double>>&inval,vector<int>&pickv,vector<vector<double>>&outvalue){

    int npickv = pickv.size();
    outvalue.clear();
    outvalue.resize(pickv.size());

    for(int i=0;i<npickv;++i){
        vector<double>&invalpv = inval[pickv[i]];
        vector<double>vbuf(nMat,0);
        int counc = 0;
        for(int j=0;j<nP;++j){
            bool isvoid = true;
            for(int k=0;k<nMat;++k)if(invalpv[j*nMat+k]!=-1)isvoid=false;

            //if(j==pickplane)assert(!isvoid);
            if(!isvoid){
                for(int k=0;k<nMat;++k)vbuf[k]+=invalpv[j*nMat+k];
                counc++;
            }

        }

        for(int k=0;k<nMat;++k)vbuf[k]/=counc;

        outvalue[i] = vbuf;
    }



}

void Graphs::GetContour(vector<vector<double> > &val,string fpath){



    vector<DynamicMesh::onePlane>planeMeshes;
    nmGraph.outputToMarker(planeMeshes);

    assert(planeMeshes.size()==nP);

    CrossSections CrossSs;
    vector<Contour> &CSs = CrossSs.CSs;



    for(int i=0;i<nP;++i){


//        struct onePlane{
//            int marker;
//            vector<double>vertices;
//            vector<unsigned int>faces2vertices;
//            vector<double>vattr;
//            vector<int>inverseVmapping;
//            vector<int>inverseFmapping;

//        };

        DynamicMesh singleplane;
        singleplane.init(planeMeshes[i].vertices,planeMeshes[i].faces2vertices);

        vector<vector<double>>curvalue(0);
        UnifiedValue(val,planeMeshes[i].inverseVmapping,i,curvalue);

        vector<double> outCtrV;
        vector<uint> outCtrE;
        vector<int>outCtrEM;
        singleplane.Contouring_EdgeMidpoint(curvalue,planeinfo[i].data(),outCtrV,outCtrE,outCtrEM);


        CSs.resize(CSs.size()+1);
        Contour &ccCs = CSs[CSs.size()-1];
        ccCs.ImportContour(planeinfo[i].data(),outCtrV,outCtrE,outCtrEM);

        //ccCs.CurveSmoothing(10);
       // cout<<"CurveSmoothing"<<endl;



    }


    CrossSs.WriteCrossSectionsToFCM(fpath+string("tctrloop.txt"),nMat);

    CrossSs.WriteCrossSectionsToFCMtCtr(fpath+string("tctr.txt"));


    CrossSs.WriteCrossSections(fpath+string("tctrContour"));

    CrossSs.WriteCrossSectionsCtrFormat(fpath+string("outctr"));






}


void Graphs::GetContour2(vector<vector<double> > &val, string fpath){


    cout<<"Contouring"<<endl;
    DynamicMesh::ContourStruct ccStruct;

    vector<int>pickV(nmGraph.nV);
    for(int i=0;i<nmGraph.nV;++i)pickV[i] = i;

    vector<vector<double>>curvalue(0);
    UnifiedValue(val,pickV,curvalue);
    nmGraph.ContouringBatch_EdgeMidpoint(curvalue,planeinfo,ccStruct);


    ccStruct.speicalmarkersmoothing(10);
    //ccStruct.CurveNetReUniformSampling(2);
//    Contour container;
//    container.ImportContour(planeinfo[0].data(),ccStruct.outCtrV,ccStruct.outCtrE,ccStruct.outCtrEM);

//    container.CurveSmoothing(500);

//    ccStruct.outCtrV = container.vertices;
    vector<DynamicMesh::ContourStruct>ccSSSSs;


    nmGraph.separateContourbyMarker(ccStruct,ccSSSSs);

    assert(ccSSSSs.size()==nP);

    CrossSections CrossSs;
    vector<Contour> &CSs = CrossSs.CSs;

    for(int i=0;i<nP;++i){


        DynamicMesh::ContourStruct &ppCC = ccSSSSs[i];
        CSs.resize(CSs.size()+1);
        Contour &ccCs = CSs[CSs.size()-1];
        ccCs.ImportContour(planeinfo[i].data(),ppCC.outCtrV,ppCC.outCtrE,ppCC.outCtrEM);

        //ccCs.CurveSmoothing(10);
       // cout<<"CurveSmoothing"<<endl;

    }

    //CrossSs.WriteCrossSectionsToFCM(fpath+string("tctrloop.txt"),nMat);

    //CrossSs.WriteCrossSectionsToFCMtCtr(fpath+string("tctr.txt"));


    CrossSs.WriteCrossSections(fpath+string("outputContour"));

    //CrossSs.WriteCrossSectionsCtrFormat(fpath+string("outctr"));



    //CrossSs.WriteCrossSectionsCSLFormat(fpath+string("Csl"));

    CrossSs.WriteCrossSectionsObjFormat(fpath+string("outputView"));





}
