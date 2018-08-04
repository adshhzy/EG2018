#ifndef DYNAMICMESH_H
#define DYNAMICMESH_H

#include"mymesh/my_mesh.h"
#include <vector>
#include <map>
#include <set>

#define NONMANIFOLDBEGIN -100


class DynamicMesh{

public:

    class Vertices{
    public:
        double pos[3];
        vector<int>ve;
        vector<int>vf;
        vector<int>vv;

        double attr;
        int marker;
        Vertices():marker(-1){}
    };
    class Edges{
      public:
        int ev[2];
        vector<int>ef;
        double attr;
        int marker;
        Edges():marker(-1){}
    };

    class Faces{
      public:
        int fv[3];
        vector<int>fe;
        vector<int>ff;
        double attr;
        int marker;
        Faces():marker(-1){}

    };

    int nV,nE,nF,nMarker,nEdgemarker;

    vector<Vertices>vertices;
    vector<Edges>edges;
    vector<Faces>faces;

    set<int>markerset;
    set<int>edge_markerset;

    void init(vector<double>&in_vertices, vector<uint>&in_faces);
    void init(vector<double>&in_vertices, vector<uint>&in_faces, vector<int>&fmarkers, map<vector<int>, int> &edgemarkers);
    void setparameters();
    void splitEdges(int iedges, int &newvertex, vector<int>&newedges);
    void propagateFaceMarkertoEdges();
    void propagateFaceMarkertoEdges(int iedges);

    void checkNeighborhoodtable();
    void writeObjfile(string filename);
    void outputToMathematica(string fpath);
    void outputToInside(vector<double>&out_vertices,vector<uint>&out_faces2v);


    void triangleMidPoint(int tind,vector<double>&midpoint);


    bool triangleInequality(int iface);

    void findBadTriangles(vector<int>&badtri);
    void findBadEdges(vector<int>&badedges,double attr);


    void reinit();

    void Contouring_EdgeMidpoint(vector<vector<double> > &val, double *planepara, vector<double> &outCtrV, vector<uint> &outCtrE, vector<int> &outCtrEM);




public:
    struct onePlane{
        int marker;
        vector<double>vertices;
        vector<unsigned int>faces2vertices;
        vector<double>vattr;
        vector<int>inverseVmapping;
        vector<int>inverseFmapping;

    };
    void outputToMarker(vector<onePlane>&planeMeshes);

    class ContourStruct{
    public:
        vector<double>outCtrV;
        vector<uint>outCtrE;
        vector<int>outCtrEM;
        vector<int>outCtrEMarker;

    private:


    public:
        ContourStruct(){}
        void speicalmarkersmoothing(int maxiter);

        void CurveNetReUniformSampling(int s);

    };


    void ContouringBatch_EdgeMidpoint(vector<vector<double>>&val,vector<vector<double>>&planepara,ContourStruct &ccstruct);


    void separateContourbyMarker(ContourStruct &ccstruct,vector<ContourStruct>&cccSs);
};







#endif // DYNAMICMESH_H
