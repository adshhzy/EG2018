#include "contour.h"
#include<iostream>
#include<fstream>
#include<set>
#include<assert.h>
#include<map>
#include <iomanip>


#include <eigen3/Eigen/Geometry>

#define REAL double
#define VOID int
#define __cplusplus
#include "triangle.h"

void planeinfo::build(vector<double>&info){

    midpoint.resize(3);
    planepara.resize(4);
    for(int i=0;i<3;++i){

        midpoint[i] = info[i];
    }

    for(int i=0;i<4;++i){

        planepara[i] = info[i+3];
    }




}

void planeinfo::build(double a1,double a2,double a3,double a4,double a5,double a6,double a7){

    midpoint.resize(3);
    planepara.resize(4);

    midpoint[0] = a1;midpoint[1] = a2;midpoint[2] = a3;

    planepara[0] = a4;planepara[1] = a5;planepara[2] = a6;planepara[3] = a7;




}

void planeinfo::build(vector<double>&midpoint,double *para){

    this->midpoint = midpoint;

    planepara.resize(4);
    for(int i=0;i<4;++i)planepara[i] = para[i];

}


void lineinfo::build(double a1,double a2,double a3,double a4,double a5,double a6){

    shootpoint.resize(3);
    dir.resize(3);

    shootpoint[0] = a1;shootpoint[1] = a2;shootpoint[2] = a3;

    dir[0] = a4;dir[1] = a5;dir[2] = a6;


}

void lineinfo::build(vector<double>&p,vector<double>&vec){

    shootpoint = p;
    dir = vec;




}


void seg::build(double a1,double a2,double a3,double a4,double a5,double a6){


    v1.resize(3);
    v2.resize(3);

    v1[0] = a1;v1[1] = a2;v1[2] = a3;

    v2[0] = a4;v2[1] = a5;v2[2] = a6;

}

void seg::build(vector<double>&vv1,vector<double>&vv2){

    v1 = vv1;
    v2 = vv2;



}


void seg::build(double *vv1,double *vv2){

    v1.resize(3);
    v2.resize(3);
    for(int i=0;i<3;++i)v1[i] = vv1[i];
    for(int i=0;i<3;++i)v2[i] = vv2[i];



}


bool isPointOnPlane(planeinfo &p1,vector<double>&v){

    return fabs(dot(p1.planepara.data(),(const double*)(v.data()))+p1.planepara[3])<THRES;


}

bool isPointOnLine(lineinfo &l,vector<double>&v){

    double vec[3],cs[3];
    minusVec(l.shootpoint.data(),v.data(),vec);
    cross(l.dir.data(),vec,cs);
    return normVec(cs)<THRES;


}
bool isPointOnSeg(seg &s,double *v){

    double vec1[3],vec2[3];

    minusVec(v,s.v1.data(),vec1);
    minusVec(s.v2.data(),v,vec2);

    if( normVec(vec1)<THRES ||  normVec(vec2)<THRES) return true;

    return fabs(cosine(vec1,vec2)-1)<THRES;
}

bool isPointOnSeg(seg &s,vector<double>&v){


    return isPointOnSeg(s,v.data());

}

double planeIntersectLine(planeinfo &p,lineinfo &l,vector<double>&intersectionV){

    /*
     planeIntersectLine[plane_, line_] :=

        Module[{pnt = {}, D, dT = 0, p},

        p = {plane[[2, 1]], plane[[2, 2]], plane[[2, 3]]};
        D = p.line[[2]];

        If[approxEqual[D, 0, 0.0001],
            If[isPointOnPlane[plane, line[[1]]], pnt = line[[1]]];,
            (
                    dT = -(p.line[[1]] + plane[[2, 4]])/D;
                    pnt = line[[1]] + dT*line[[2]];
            )
        ];

       {pnt, dT}
    ]
    */


    double D = dot(p.planepara.data(),l.dir.data());

    intersectionV.clear();intersectionV.resize(3,INFTY_DOUBLE);

    if(fabs(D)<1e-5){
        if(isPointOnPlane(p,l.shootpoint)){
            intersectionV = l.shootpoint;
            return 0;
        }else{

            return INFTY_DOUBLE;
        }
    }else{

        double dT = -(dot(p.planepara.data(),l.shootpoint.data())+p.planepara[3])/D;
        weightedAddVec(1.,dT,l.shootpoint.data(),l.dir.data(),intersectionV.data());
        return dT;
    }




}


bool planeIntersectPlane(planeinfo &p1,planeinfo &p2,lineinfo &l_re){


    /*
        planeIntersectPlane[plane1_, plane2_] :=

            Module[{line = {}, pnt1, norm1, norm2, vecToLine, vecLine, pntLine,
                     dT},
                pnt1 = plane1[[1]];
                norm1 = Normalize[{plane1[[2, 1]], plane1[[2, 2]], plane1[[2, 3]]}];
                norm2 = Normalize[{plane2[[2, 1]], plane2[[2, 2]], plane2[[2, 3]]}];
                If[! approxEqual[Abs[norm1.norm2], 1, 0.0001], (
                    vecLine = Normalize[Cross[norm1, norm2]];
                    vecToLine = Normalize[Cross[vecLine, norm1] ];
                    {pntLine, dT} = planeIntersectLine[plane2, {pnt1, vecToLine}];
                    line = {pntLine, vecLine};
                )];

                line
         ]

    */

    double dd = fabs(cosine(p1.planepara.data(),p2.planepara.data()))-1;
    if(fabs(dd)<THRES){
        return false;

    }else{

        l_re.dir.resize(3);
        vector<double>vecToLine(3);
        cross(p1.planepara.data(),p2.planepara.data(),l_re.dir.data());
        normalize(l_re.dir.data());

        cross(l_re.dir.data(),p1.planepara.data(),vecToLine.data());
        normalize(vecToLine.data());

        lineinfo ll;
        ll.build(p1.midpoint,vecToLine);
        planeIntersectLine(p2,ll,l_re.shootpoint);

        return true;


    }








}

bool lineIntersectLine(lineinfo &l1,lineinfo &l2,vector<double>&intersectV,double &dT1,double &dT2){


    auto pv1 = l1.shootpoint.data();
    auto pv2 = l2.shootpoint.data();

    auto pdir1 = l1.dir.data();
    auto pdir2 = l2.dir.data();


    double vec12[3],vec21[3];
    minusVec(pv2,pv1,vec12);
    minusVec(pv1,pv2,vec21);

    double cross212[3],cross21[3],cross121[3];
    cross(pdir2,vec12,cross212);
    double a = normVec(cross212);
    cross(pdir2,pdir1,cross21);
    double b = normVec(cross21);
    cross(pdir1,vec21,cross121);
    double c = normVec(cross121);




    if(fabs(b)<THRES){

        //cout<<"1"<<endl;
        return false;

    }else if(normVec(cross212)<THRES || normVec(cross121)<THRES || fabs(fabs(cosine(cross212,cross121))-1)<THRES){
        //cout<<"2"<<endl;
        dT1 = a/b;
        intersectV.resize(3);
        if(dot(cross212,cross21)<0)dT1*=-1;
        weightedAddVec(1.,dT1,pv1,pdir1,intersectV.data());

        dT2 = c/b;
        if(dot(cross121,cross21)>0)dT2*=-1;

        return true;

    }else{

        //cout<<"3"<<endl;
        return false;
    }


}

bool segIntersectSeg(seg &s1,seg &s2,vector<double>&intersectionV,double &dT1,double &dT2){

        vector<double>dir1(3),dir2(3);
        minusVec(s1.v2.data(),s1.v1.data(),dir1.data());
        minusVec(s2.v2.data(),s2.v1.data(),dir2.data());

        lineinfo l1(s1.v1,dir1);
        lineinfo l2(s2.v1,dir2);
        if( lineIntersectLine(l1,l2,intersectionV,dT1,dT2) ){

            //cout<<"adsad"<<endl;
            if( (dT1<0 || dT1>1) && !(fabs(dT1)<THRES || fabs(dT1-1)< THRES ) )return false;
            if( (dT2<0 || dT2>1) && !(fabs(dT2)<THRES || fabs(dT2-1)< THRES ) )return false;

            return true;

        }
        return false;


}



bool isLinesEqual(lineinfo &l1,lineinfo &l2){

    double d = fabs(fabs(cosine(l1.dir.data(),l2.dir.data()))-1);
    if(d>THRES)return false;

    double vec12[3];
    add(l2.shootpoint.data(),l2.dir.data(),vec12);
    minusVec(l1.shootpoint.data(),vec12,vec12);


    d = fabs(fabs(cosine(l1.dir.data(),vec12))-1);

    if(d<THRES)return true;
    else return false;


}




double findLineT(lineinfo &l,vector<double>&p){

    double vec[3];
    minusVec(p.data(),l.shootpoint.data(),vec);
    double dT = normVec(vec)/normVec(l.dir.data());
    if(dot(vec,l.dir.data())<0)dT*=-1;
    return dT;


}

double findLineT_seg(seg &s,vector<double>&p){

    double dir[3];
    minusVec(s.v2.data(),s.v1.data(),dir);
    double vec[3];
    minusVec(p.data(),s.v1.data(),vec);
    double dT = normVec(vec)/normVec(dir);
    if(dot(vec,dir)<0)dT*=-1;
    return dT;


}

double point2segDistance(double *p1,double *p2, double *pnt, double &dT){

    double vec[3],pp1[3];
    minusVec(p2,p1,vec);
    minusVec(pnt,p1,pp1);

    double d1 = dot(vec,pp1);
    double d2 = dot(pp1,pp1);
    if(d1<=0)dT = 0;
    else if(d2<d1)dT = 1;
    else dT = d1/d2;

    weightedAddVec(1,dT,p1,vec,pp1);
    double dist = _VerticesDistance(pnt,pp1);
    return dist;


}

double point2segDistance(double *p1,double *p2, double *pnt){
    double dT;
    return point2segDistance(p1,p2,pnt,dT);

}

double calculateWindingnumber(double *p_vpos,double *plane_para,vector<double>&vseq,double *mindist){
    int nv = vseq.size()/3;
    vector<double>rays(nv*3);
    auto p_rayd = rays.data();
    vector<double>dist(nv);

    for(int i=0;i<nv;++i){
        auto p_rays = p_rayd+i*3;
        minusVec(vseq.data()+i*3,p_vpos,p_rays);
        dist[i] = normalize(p_rays);

    }
    int minind = min_element(dist.begin(),dist.end())-dist.begin();
    if(dist[minind]<THRES){
        if(mindist)*mindist = 0;
        return 0;
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

    if(mindist){

        int m1 =  (minind-1+nv)%nv, m2 =  (minind+1)%nv;
        double d1 = point2segDistance(vseq.data()+m1*3,vseq.data()+minind*3,p_vpos);
        double d2 = point2segDistance(vseq.data()+m2*3,vseq.data()+minind*3,p_vpos);
        *mindist = min(d1,d2);
    }


    return windingnumber;


}

/******************************************************************/
/******************************************************************/

void ComputeAngleAxisMatrix(Eigen::AngleAxisd &m_aa,double *oriVec,double *newVec){
    Eigen::Vector3d oriVecE,newVecE,nc;
    for(int j=0;j<3;++j)oriVecE(j) = oriVec[j];
    for(int j=0;j<3;++j)newVecE(j) = newVec[j];
    nc = oriVecE.cross(newVecE);
    if(len(nc.data())<1e-6){nc(0) = 1.0;nc(1) = 0.0;nc(2) = 0.0;}
    nc.normalize();
    double angle = acos(cosine(oriVec,newVec));
    m_aa = Eigen::AngleAxisd(angle,nc);
}

void RotateVertices(vector<double>&inv,double *oriVec,double *newVec,vector<double>&outv){

    int nvs = inv.size()/3;
    Eigen::AngleAxisd t;
    vector<Eigen::Vector3d>ori(nvs);
    vector<Eigen::Vector3d>tranformed(nvs);

    ComputeAngleAxisMatrix(t,oriVec,newVec);
    for(int j = 0;j<nvs;++j){
        auto p_spv = inv.data()+j*3;
        for(int k=0;k<3;++k)ori[j](k) = p_spv[k];
    }
    for(int j = 0;j<nvs;++j)tranformed[j] = t*ori[j];

    outv.resize(inv.size());
    for(int j = 0;j<nvs;++j){
        auto p_spv = outv.data()+j*3;
        for(int k=0;k<3;++k)p_spv[k] = tranformed[j](k);
    }

}

void ProjectXYcoordinate(vector<double>&Vpos3D,vector<double>&Vpos2D){
    int n_v = Vpos3D.size()/3;
    Vpos2D.clear();
    for(int i=0;i<n_v;++i){
        auto p_v = Vpos3D.data()+i*3;
        Vpos2D.push_back(p_v[0]);Vpos2D.push_back(p_v[1]);
    }

}


void PatchZcoordinate(double z,vector<double>&Vpos2D,vector<double>&Vpos3D){
    int n_v = Vpos2D.size()/2;
    Vpos3D.clear();
    for(int i=0;i<n_v;++i){
        auto p_v = Vpos2D.data()+i*2;
        Vpos3D.push_back(p_v[0]);Vpos3D.push_back(p_v[1]);Vpos3D.push_back(z);
    }

}

void ThresholdColoring(double degree,uchar *p_rgb);



void CallTriangleForSinglePlane(vector<double>&Vpos2D,vector<uint>&E2V, string &cmline,vector<double>&out_Vpos2D,vector<uint>&out_Triangles){


    triangulateio in, out;
    int nv = Vpos2D.size()/2;
    int ne = E2V.size()/2;


    in.numberofpoints = nv;
    in.numberofpointattributes = 0;
    in.pointlist = new REAL[in.numberofpoints*2];

    for(int i=0;i<nv*2;++i)in.pointlist[i] = Vpos2D[i];

    in.pointattributelist = NULL;
    in.pointmarkerlist = new int[in.numberofpoints];
    for(int i=0;i<in.numberofpoints;++i)in.pointmarkerlist[i] = i+10000;

    in.numberofsegments = ne;
    in.segmentlist = new int[in.numberofsegments*2];

    for(int i=0;i<ne*2;++i)in.segmentlist[i] = E2V[i];


    in.numberofholes = 0;
    in.numberofregions = 0;
    in.regionlist = NULL;


    out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
    /* Not needed if -N switch used or number of point attributes is zero: */
    out.pointattributelist = (REAL *) NULL;
    out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
    out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
    /* Not needed if -E switch used or number of triangle attributes is zero: */
    out.triangleattributelist = (REAL *) NULL;
    out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
    /* Needed only if segments are output (-p or -c) and -P not used: */
    out.segmentlist = (int *) NULL;
    /* Needed only if segments are output (-p or -c) and -P and -B not used: */
    out.segmentmarkerlist = (int *) NULL;
    out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
    out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

    //string cmline("pczQYY");//pczPNEO  pczAevn pczenQYY
    //cout<<"11"<<endl;
    try{
        triangulate((char *)(cmline.data()), &in, &out, NULL);
    }catch(...){
        cerr << "error" << endl;
    }

    //triangulate("pczQYY", &in, &out, NULL);
    //cout<<"22"<<endl;
    cout<<"triangulate "<<out.numberofpoints<<' '<<in.numberofpoints<<endl;
    out_Triangles.clear();
    for (int i = 0; i < out.numberoftriangles; i++) {
        for (int j = 0; j < out.numberofcorners; j++) {
            out_Triangles.push_back( uint(out.trianglelist[i * out.numberofcorners + j]) );
        }
    }

    out_Vpos2D.clear();
    for (int i = 0; i < out.numberofpoints; i++) {
        for (int j = 0; j < 2; j++) {
            out_Vpos2D.push_back(out.pointlist[i * 2 + j]);
        }
    }

    // for (int i = 0; i < out.numberofpoints; i++) cout<<out.pointmarkerlist[i]<<' ';cout<<endl;



    //        for (int i = 0; i < out.numberoftriangles; i++) {
    //            for (int j = 0; j < out.numberofcorners; j++) {
    //                cout<< out.trianglelist[i * out.numberofcorners + j]<<' ';
    //            }
    //            cout<<endl;
    //        }

    //        cout<<endl;

    delete []in.pointlist;
    delete []in.segmentlist;
    delete []in.pointmarkerlist;
    free(out.pointlist);
    free(out.pointattributelist);
    free(out.trianglelist);
    free(out.triangleattributelist);


}




void CallTriangleForSinglePlane_Voronoi(vector<double>&Vpos2D,vector<uint>&E2V, string &cmline, vector<double>&out_Vpos2D,vector<uint>&out_Triangles){


    triangulateio in, out;
    int nv = Vpos2D.size()/2;
    int ne = E2V.size()/2;


    in.numberofpoints = nv;
    in.numberofpointattributes = 0;
    in.pointlist = new REAL[in.numberofpoints*2];

    for(int i=0;i<nv*2;++i)in.pointlist[i] = Vpos2D[i];

    in.pointattributelist = NULL;
    in.pointmarkerlist = new int[in.numberofpoints];
    for(int i=0;i<in.numberofpoints;++i)in.pointmarkerlist[i] = i+10000;

    in.numberofsegments = ne;
    in.segmentlist = new int[in.numberofsegments*2];

    for(int i=0;i<ne*2;++i)in.segmentlist[i] = E2V[i];


    in.numberofholes = 0;
    in.numberofregions = 0;
    in.regionlist = NULL;


    out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
    /* Not needed if -N switch used or number of point attributes is zero: */
    out.pointattributelist = (REAL *) NULL;
    out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
    out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
    /* Not needed if -E switch used or number of triangle attributes is zero: */
    out.triangleattributelist = (REAL *) NULL;
    out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
    /* Needed only if segments are output (-p or -c) and -P not used: */
    out.segmentlist = (int *) NULL;
    /* Needed only if segments are output (-p or -c) and -P and -B not used: */
    out.segmentmarkerlist = (int *) NULL;
    out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
    out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

    //string cmline("pczQYY");//pczPNEO  pczAevn pczenQYY

    triangulate((char *)(cmline.data()), &in, &out, NULL);

    cout<<"triangulate "<<out.numberofpoints<<' '<<in.numberofpoints<<endl;
    out_Triangles.clear();
    for (int i = 0; i < out.numberoftriangles; i++) {
        for (int j = 0; j < out.numberofcorners; j++) {
            out_Triangles.push_back( uint(out.trianglelist[i * out.numberofcorners + j]) );
        }
    }

    out_Vpos2D.clear();
    for (int i = 0; i < out.numberofpoints; i++) {
        for (int j = 0; j < 2; j++) {
            out_Vpos2D.push_back(out.pointlist[i * 2 + j]);
        }
    }

    // for (int i = 0; i < out.numberofpoints; i++) cout<<out.pointmarkerlist[i]<<' ';cout<<endl;



    //        for (int i = 0; i < out.numberoftriangles; i++) {
    //            for (int j = 0; j < out.numberofcorners; j++) {
    //                cout<< out.trianglelist[i * out.numberofcorners + j]<<' ';
    //            }
    //            cout<<endl;
    //        }

    //        cout<<endl;

    delete []in.pointlist;
    delete []in.segmentlist;
    delete []in.pointmarkerlist;
    free(out.pointlist);
    free(out.pointattributelist);
    free(out.trianglelist);
    free(out.triangleattributelist);


}





void RatateToZaxis3Dto2D(vector<double>&Vpos3D,vector<double>&plane_para,vector<double>&Vpos2D,double &zpos){


    double zaxis[3] = {0,0,1};

    vector<double>rotateVpos;

    RotateVertices(Vpos3D,plane_para.data(),zaxis,rotateVpos);

    ProjectXYcoordinate(rotateVpos,Vpos2D);

    if(rotateVpos.size()>2)zpos = rotateVpos[2];

}
void RatateToNormalAxis2Dto3D(vector<double>&Vpos2D,vector<double>&plane_para,double &zpos,vector<double>&Vpos3D){


    double zaxis[3] = {0,0,1};
    vector<double>tmpVpos3D;
    PatchZcoordinate(zpos,Vpos2D,tmpVpos3D);
    RotateVertices(tmpVpos3D,zaxis,plane_para.data(),Vpos3D);

}

void GetFaceCenter(vector<double>&vpos,vector<uint>&faces,vector<double>&faceCenters){
    int nf = faces.size()/3;
    faceCenters.clear();
    faceCenters.resize(nf*3,0);
    for(int i=0;i<nf;++i){
        auto p_fv = faces.data()+i*3;
        auto p_fc = faceCenters.data()+i*3;
        for(int j = 0;j<3;++j){
            auto p_v = vpos.data()+p_fv[j]*3;
            for(int k=0;k<3;++k)p_fc[k]+=p_v[k];
        }
        for(int k=0;k<3;++k)p_fc[k]/=3;

    }

}

void CallTriangleForSinglePlane_3D(vector<double>&Vpos3D,vector<uint>&E2V, vector<double>&plane_para, string &cmline, vector<double>&out_Vpos3D,vector<uint>&out_Triangles){

    vector<double>outVpos2D,rotateVpos2D;
    double zpos;

    RatateToZaxis3Dto2D(Vpos3D,plane_para,rotateVpos2D,zpos);

    //cout<<"44"<<endl;
    CallTriangleForSinglePlane(rotateVpos2D,E2V,cmline,outVpos2D,out_Triangles);
    //cout<<"55"<<endl;
    RatateToNormalAxis2Dto3D(outVpos2D,plane_para,zpos,out_Vpos3D);


}



