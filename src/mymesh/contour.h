#ifndef CONTOUR_H
#define CONTOUR_H



#include<vector>
#include<math.h>
#include"utility.h"

using namespace std;
using namespace MyUtility;

#define INFTY_DOUBLE std::numeric_limits<double>::max()
#define THRES 1.e-5


class planeinfo{
public:
    vector<double>midpoint;
    vector<double>planepara;

    planeinfo(){}

    void build(vector<double>&info);

    void build(double a1,double a2,double a3,double a4,double a5,double a6,double a7);
    void build(vector<double>&midpoint,double *para);

};

class lineinfo{
public:
    vector<double>shootpoint;
    vector<double>dir;

    lineinfo(){}

    lineinfo(vector<double>&p,vector<double>&vec):shootpoint(p),dir(vec){}

    //void build(vector<double>&info);

    void build(double a1,double a2,double a3,double a4,double a5,double a6);

    void build(vector<double>&p,vector<double>&vec);
};

class seg{
public:
    vector<double>v1;
    vector<double>v2;

    seg(){}

    seg(vector<double>&vv1,vector<double>&vv2):v1(vv1),v2(vv2){}

    seg(double *vv1,double *vv2){build(vv1,vv2);}

    //void build(vector<double>&info);

    void build(double a1,double a2,double a3,double a4,double a5,double a6);

    void build(vector<double>&vv1,vector<double>&vv2);
    void build(double *vv1,double *vv2);
};





class Contour{
public:
    int n_vertices,n_edges,n_materials;
    double plane_para[4];
    vector<double>vertices;
    vector<uint>edges;
    vector<int>m1;
    vector<int>m2;
public:
    vector<double>leftcorner;
    vector<double>midpoint;
    vector<double>rightcorner;
    planeinfo plane_info;

public:
    Contour():n_vertices(0),n_edges(0),n_materials(0){}
    void ReadContour(ifstream &in);
    void splitContour(vector<Contour> &outC);
    void WriteContour(ofstream &out);
    void WriteContourCtrFormat(ofstream &out);
    void WriteContourCSLFormat(ofstream &out);
    void CalSelfProperty();
    bool ValidationCheck();
    void BuildContourFromImage(vector<vector<int>>&im,double z);
    void CurveSmoothing(int maxiter);
    void ImportContour(double *para, vector<double>&inCtrV, vector<uint>&inCtrE, vector<int>&inCtrEMat);
    void MergeMaterial(vector<int>&mappingMat);
    void RotatingandShifting(double r_angle, vector<double> &shiftbeforeR, vector<double> &raxis, vector<double> &shiftv);
    void AddContour(Contour &mc);
    void InversePlaneNormal();
    void MergeBoundaryBox(vector<double>&inCtrV, vector<uint>&inCtrE, vector<int>&inCtrEMat);
    void OutputContour(double *para, vector<double>&outCtrV, vector<uint>&outCtrE, vector<int>&outCtrEMat);
    void LableTrianglesWithWindingNumber(vector<double> &testFaceCenter, vector<int>&outfaceMat);
    void ComputeMat2VerticesLoop(vector<vector<int> > &verticesloops, vector<vector<double>>&verticesSeqs, vector<int>&loops2mats);

    int getnumberofloops();
    void WriteContourToFCM(ofstream &out);
    void WriteContourToFCMtCtr(ofstream &out);

    inline uint* ev_begin(uint e_ind){ return &(edges[e_ind * 2]); }
    inline uint* ev_end(uint e_ind){ return &(edges[e_ind * 2 + 2]); }
    inline double* v_begin(int v_ind){ return &(vertices[v_ind * 3]); }
    inline double* v_end(int v_ind){ return &(vertices[(v_ind+1) * 3]); }
    double VerticesDistance(int first, int second)const{
        double dist = 0.0;
        auto p_v1 = &(vertices[first*3]);
        auto p_v2 = &(vertices[second*3]);
        for (int i = 0; i < 3; ++i)
            dist += pow((p_v1[i] - p_v2[i] ), 2);
        return sqrt(dist);
    }
};
class CrossSections{
public:
    int n_materials,n_CrossSec,total_nver,total_nedges;
    vector<Contour>CSs;
    vector<int>acc_nver;
public:
    double width,height;

public:
    vector<double>leftcorner;
    vector<double>midpoint;
    vector<double>rightcorner;
    vector<planeinfo>planes_info;

public:
    CrossSections():n_CrossSec(0){}
    CrossSections(vector<Contour>&a):CSs(a),n_CrossSec(a.size()){
        n_materials=-1;
        for(auto &b:CSs)n_materials = max(n_materials,b.n_materials);
    }
    void naiveConstruction();
    void CalSelfProperty();
    void Stackup(vector<double>&points,vector<uint>&edges);
    void MapMaterials();
    int WriteCrossSections(string filename);
    int WriteCrossSections(string filename,vector<int>&pickedCs);
    void ReadCrossSectionsfromImages(vector<vector<vector<int>>>&im, double zdis);
    int SpecialMerging(string infilename, vector<double>&points, vector<uint> &edges);
    void ManipulateCrossSections();
    bool ReadCrossSections(string filename);
    bool ReadCrossSections(    vector<vector<double>>&planeinfo, vector<vector<double>>&plane2vertices,
                                vector<vector<uint>>&plane2edges, vector<vector<int>>&plane2edgesMat);
    void ShiftCrossSections(vector<vector<double>>&trlvec);

    int WriteCrossSectionsToFCM(string filename, int nMat);
    int WriteCrossSectionsToFCMtCtr(string filename);
    int WriteCrossSectionsCtrFormat(string filename);

    int WriteCrossSectionsCSLFormat(string filename);

    int WriteCrossSectionsObjFormat(string filename);

    void SubSetPlanes(vector<int> &reserveplanes);



public:

    vector<double>V;

    vector<vector<uint>>E;

    vector<vector<uint>>T;

    vector<int>isectIndices;
    vector<vector<int> >lineIndices;
    vector<vector<int> >line2plane;

    vector<vector<int> >plane2lines;

    vector<vector<int>>line2Edges;

    vector<vector<int>>ver2planes;
    vector<vector<int> >planeV2globalV;
    vector<int>planeVcutoffInd;


    vector<vector<double> >newV;
    vector<vector<uint>>newE;

    //vector<>


    vector<vector<double> >val;

    vector<int>nodeMarkers;


public:

    //io interface
    void writePlanesVal(string filename);
    void writeSegs(string filename, vector<seg>&segs);
    void writeSegs2(string filename, vector<seg>&segs,vector<double>&newv);

    //main pipeline interface
    void FCMGeometryProcessPipeline(string infilename, string outpath, double iline_step, double portion, string triangle_cmd);

public:

    //find the bounding box of the input
    void FindBoundingBox(vector<planeinfo > &boundingbox, double coef);

    //find the intersection line between cross-sectional planes and the bounding box
    void findPlaneBoundLines(planeinfo &p, vector<planeinfo > &boundingbox, vector<seg>&segments);

    //find the intersection line among cross-sectional planes
    void findPlaneLines(vector<planeinfo > &planes, vector<planeinfo > &boundingbox, vector<seg>&p2psegs, vector<vector<int>>&planelines, vector<vector<seg> > &p2bsegs);

    //triangluate all cross-sectional planes with common vertices on intersection lines
    void findPlaneGraphMM(vector<seg>&p2psegs, vector<vector<int>>&planelines, vector<vector<seg> > &p2bsegs, double sampW, double boxSampCoef);


    //pipeline for calling Trianlge for each plane
    void CallTriangle(string &cmline, string outprefix);


public:

    //compute the sign distance function on each vertices for each plane and each label
    void findImplicitValue(string outprefix);

public:
    //determine the active vertices
    //(whose value of the implicit function will be explicitly consider as
    //continous variables for quadratic programming)
    void LineSafeVertices(double portion);






};

















#endif // CONTOUR_H
