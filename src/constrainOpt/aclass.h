#ifndef ACLASS_H
#define ACLASS_H


#include <vector>
#include "gurobi_c++.h"
#include "mymesh/contour.h"
#include "DynamicMesh.h"
using namespace std;

class Graphs{
public:

    bool isInit;
    double theta = 1e-4;

    double lambda;

public:

    DynamicMesh nmGraph;

public:

    //deprecated
    void GetContour(vector<vector<double>>&val, string fpath);

    //global contouring, do contouring for all planes at the same time
    void GetContour2(vector<vector<double>>&val, string fpath);

    void UnifiedValue(vector<vector<double>>&inval,vector<int>&pickv,int pickplane,vector<vector<double>>&outvalue);
    void UnifiedValue(vector<vector<double>>&inval,vector<int>&pickv,vector<vector<double>>&outvalue);



public:
    int nV;
    int nL;
    int nP;
    int nMat;

    vector<double>vpos;

    //0: ordinary vertices;  1: inconsistent vertices;
    //2: consistent but small difference between first and second label(portion)
    vector<int>nodeMarkers;


    vector<vector<uint> >E;

    vector<vector<double>>planeinfo;
    vector<vector<int>>lines;

    vector<vector<int>>plane2vers;
    vector<vector<int>>ver2planes;

    vector<vector<int>>plane2lines;

    //original value from sign distance function on each plane
    vector<vector<double>>oldval;

    //swapping buffers
    vector<vector<double>>newval;
    vector<vector<double>>solval;


    //number of continuous variables = nactive vertices X n material
    vector<int>plane2nActiveVars;

     //number of continuous variables that can be inferred from active variables
    vector<int>plane2nUnActiveVars;

    //number of active vertices (nodeMarker != 0)
    vector<int>nplaneActiveVertices;
    vector<vector<int>>planeActiveVertices;

    vector<int>nplaneVars;

    //reorder the variables by [active vertices, inactive vertices]
    vector<vector<int>>varInvOrder;
    vector<vector<int>>varOrder;

    //matrix representating the quadratic energy
    vector<vector<double>>pEn;

    //matrix for infering inactive variables from active variables
    vector<vector<double>>pTran;

    vector<int>plane2nActiveVertices;
    vector<int>plane2nUnActiveVertices;

private:

    vector<int>isectIndices;
    vector<int>linesEdges;
    vector<vector<int>>lines2Edges;
    vector<vector<int>>lineVertices2nv;

private:

    vector<vector<double>>planev2oldval;
    vector<vector<double>>planev2newval;
    vector<vector<double>>planev2solveval;

    vector<vector<int>>globalv2planeActiveV;

    vector<vector<int>>mapplaneV2globalV;
    vector<int>mapV2ActiveV;

    vector<int>globalactiveV;


    vector<int>FirstlabelsOnNewVal;
    vector<int>SecondlabelsOnNewVal;


    //buffer for indicating whether need to do flipping on active vertices
    vector<int>VertexMarker;



    int nFlipping;
    int nTotalActiveV;


    vector<int>isChanged;



private:
    vector<int>r_iters;

private:

    vector<vector<double>>gradientvector;
    vector<vector<double>>gradientvector_tmp;

    vector<double>flipEnergy;
    vector<double>flipGradient;
    vector<double>rowSum;

    vector<double>previousFlipGradient;
    vector<double>previousEnergy;


private:

    //gurobi model for solving quadratic programming
    GRBEnv objenv;
    GRBModel objmodel;

    vector<vector<GRBVar>> conVars;
    vector<GRBVar>binaryVars;

public:
    Graphs():objenv(GRBEnv()),objmodel(GRBModel(objenv)){

        objmodel.getEnv().set(GRB_IntParam_OutputFlag, 0);

    }



public:

    struct nodeinfo{
        //vector<double>val;
        double *pval;
        int iPnt;
        int label;
        bool isConstrained;
        nodeinfo(double *pval,int iPnt,int label,bool isc):pval(pval),iPnt(iPnt),label(label),isConstrained(isc){}
        nodeinfo(int iPnt,int label,bool isc):iPnt(iPnt),label(label),isConstrained(isc){}

    };

    struct triple{
        int i;
        int j;
        double val;
        triple(int i,int j,double val):i(i),j(j),val(val){}
    };


public:
    bool CrossectionsPipeline(CrossSections &CS, string outpath, double lambda, int maxiter);

private:
    bool readGraphsFromCrossections(CrossSections &CS);

private:

    int averageAcrossPlanes(vector<vector<double>>&val,int iPnt, vector<double>&ave);
    int labelonPlane(vector<vector<double>>&val,int iPnt,int iplane);

    void solveQuadraticProgramming(vector<triple>&M, vector<triple>&Ab, int n, vector<double> &solveval);


private:

    //perform quadratic programming on each line to get the initial label
    void solveLine();
    void buildMatrixandSolveLine(vector<nodeinfo>&nodes, vector< triple >&graphconnection_triple, vector<double>&solveval);

    //compute the matrix for quadratic energy
    void solvePlane();
    void BuildMatrixPlane(vector<nodeinfo>&nodes, vector< triple >&graphconnection_triple, vector<int>&nodeorder, int cutind, int iplane);

private:

    //main iteration
    void doubleIter_opt(int maxIter = 30);

private:

    //model set up
    void InitialobjModel();
    void InitializeMarkers();
    void IterInner_prebuild(bool isrefresh = false);

private:

    //solving quadratic programming using gurobi
    double solValues_modifyingConstraints_flippinglabel_QPonly();

    //extract result
    double extractResultfromGurobi(GRBModel &model);


private:

    //perfrom computation on flipping label
    double IterInnerwithMarkersandSelectiveLabels_mainIter(double beginningEn,vector<vector<int>>&comps);

    //compute the vertices and associated labels that are going to be flip
    double EpsilonIterwithMarkersandSelectiveLabels_Inner(double beginningEn);

private:

    //ordering candidate of flipping labels by gradient
    vector<vector<int>> orderingComps_markers(vector<vector<int>>&in_comps);
    void computeFlipGradient();


private:
    void debugline(int iline, int iplane);

    void debugline2(vector<vector<double> >&val);

    void printinfo(int globalv, vector<vector<double> >&val);


};




















#endif // ACLASS_H
