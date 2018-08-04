#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include<mymesh/readers.h>
#include<mymesh/contour.h>
using namespace std;

void TransformCtrIntoMathematica(){



    //string filename("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/FerretBrain/ferret1");
    string filename("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/Liver/liver1");
    //string filename("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/Toy/blob");



    string infilename = filename+string(".ctr");
    ifstream reader(infilename.data(), ofstream::in);
    if (!reader.good()) {
        cout << "Can not open file: " << infilename << endl;
    }else {
        cout << "Reading: "<<infilename<<endl;
    }


    int nP;
    reader>>nP;
    vector<vector<double>>planeinfo(nP);
    vector<int>plane2nV(nP);
    vector<int>plane2nE(nP);
    vector<vector<double>>plane2vertices(nP);
    vector<vector<unsigned int>>plane2edges(nP);
    vector<vector<int>>plane2edgesMat(nP);


    for(int i=0;i<nP;++i){
        auto &p_planeinfo = planeinfo[i];
        p_planeinfo.resize(11);
        for(int j=0;j<11;++j){
            reader>>p_planeinfo[j];
        }

        reader>>plane2nV[i];
        reader>>plane2nE[i];

        auto &p_plane2vertices = plane2vertices[i];
        p_plane2vertices.resize(plane2nV[i]*3);
        for(int nv = plane2nV[i]*3,j=0;j<nv;++j){
            reader>>p_plane2vertices[j];
        }

        auto &p_plane2edges = plane2edges[i];
        auto &p_plane2edgesMat = plane2edgesMat[i];
        p_plane2edges.resize(plane2nE[i]*2);
        p_plane2edgesMat.resize(plane2nE[i]*2);
        for(int j=0;j<plane2nE[i];++j){
            int i1 =j*2, i2=i1+1;
            reader>>p_plane2edges[i1];
            reader>>p_plane2edges[i2];
            reader>>p_plane2edgesMat[i1];
            reader>>p_plane2edgesMat[i2];
        }


    }



    /*******************************************************/
    /*******************************************************/


    CrossSections CS;
    CS.ReadCrossSections(planeinfo,plane2vertices,plane2edges,plane2edgesMat);

    //CS.ManipulateCrossSections();

    CS.WriteCrossSectionsToFCM(filename+string("_tctr.txt"),CS.n_materials);
    CS.WriteCrossSectionsToFCMtCtr(filename+string("_loop.txt"));


    return;





    /*******************************************************/
    /*******************************************************/





    string outfilename = filename+string("_loop.txt");
    ofstream outer(outfilename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not open file: " << outfilename << endl;
    }else {
        cout << "Writing: "<<outfilename<<endl;
    }

    outer<<setprecision(12);
    outer<<fixed;
    auto oneloop = [](ofstream &out,vector<double>&vpos,vector<unsigned int>&edges,int mat){
        out<<"{";

            out<<"{";
            auto p_v = vpos.data();
            int nv = vpos.size()/3;
            int j=0;
            for(j=0;j<nv-1;++j){
                out<<"{";
                for(int k=0;k<2;++k)out<<p_v[j*3+k]<<",";
                out<<p_v[j*3+2];
                out<<"}, ";
            }
            out<<"{";
            for(int k=0;k<2;++k)out<<p_v[j*3+k]<<",";
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
            out<<"{";out<<nv<<","<<0;out<<"}";
            out<<"},";

            out<<mat+1;

        out<<"}";
    };

    outer<<"{";
    outer<<"{";

        for(int i=0;i<nP;++i){
            outer<<"{";

            for(int j=0;j<2;++j){

                oneloop(outer,plane2vertices[i],plane2edges[i],j);
                if(j!=1)outer<<",";

            }


            if(i!=nP-1)outer<<"},";
            else outer<<"}";
        }

    outer<<"},";

    outer<<"{";

    for(int i=0;i<nP;++i){
        outer<<"{";
        auto &p_planeinfo = planeinfo[i];
        outer<<"{"<<p_planeinfo[4]<<','<<p_planeinfo[5]<<','<<p_planeinfo[6]<<"},";
        outer<<"{"<<p_planeinfo[0]<<','<<p_planeinfo[1]<<','<<p_planeinfo[2]<<','<<p_planeinfo[3]<<"},";
        outer<<"{"<<p_planeinfo[8]<<','<<p_planeinfo[9]<<','<<p_planeinfo[10]<<"}";
        if(i!=nP-1)outer<<"},";
        else outer<<"}";
    }

    outer<<"},";


    outer<<2;

    outer<<"}";



    outer.close();
}


void TransformMingToSuf(){

    vector<double>vertices;
    vector<unsigned int>faces2vertices;
    vector<double>vertices_normal;


    string filepath("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/visualization2/ferretbrain/");
    string fobjname = filepath + string("Result_afterSmooth.obj");

    readObjFile(fobjname,vertices,faces2vertices,vertices_normal);

    string filename = filepath + string("ContuorEdges.txt");
    ifstream reader(filename.data(), ofstream::in);
    if (!reader.good()) {
        cout << "Can not open the OFF file " << filename << endl;
        return;
    }else {
        cout << "Reading: "<<filename<<endl;
    }

    int n_cedge;
    vector<int>CtrEdges;

    reader>>n_cedge;
    CtrEdges.resize(n_cedge*2);

    for(int i=0;i<CtrEdges.size();++i)reader>>CtrEdges[i];

    int n_faces = faces2vertices.size()/3;
    vector<int>faceMat(n_faces*2,0);
    for(int i=0;i<n_faces;++i){
        faceMat[i*2] = 1;
    }

    string outfilename = filepath + string("outsuf");



    writeSufFile(outfilename,vertices,faces2vertices,faceMat,CtrEdges);











}





void TransformCslIntoMathematica(){




    string filename("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/Bermano/");




    string infilename = filename+string("contour_timepoint1.csl");
    ifstream reader(infilename.data(), ofstream::in);
    if (!reader.good()) {
        cout << "Can not open file: " << infilename << endl;
    }else {
        cout << "Reading: "<<infilename<<endl;
    }


    string header;
    reader>>header;
    int nP, nMat;
    reader>>nP; reader>>nMat;
    vector<vector<double>>planeinfo(nP);
    vector<int>plane2nV(nP);
    vector<int>plane2nE(nP);
    vector<vector<double>>plane2vertices(nP);
    vector<vector<unsigned int>>plane2edges(nP);
    vector<vector<int>>plane2edgesMat(nP);


    for(int i=0;i<nP;++i){
        int val;
        int ncomp;
        reader>>val;
        reader>>plane2nV[i];
        reader>>ncomp;
        auto &p_planeinfo = planeinfo[i];
        p_planeinfo.resize(4);
        for(int j=0;j<4;++j){
            reader>>p_planeinfo[j];
        }

//        reader>>plane2nV[i];
//        reader>>plane2nE[i];

        auto &p_plane2vertices = plane2vertices[i];
        p_plane2vertices.resize(plane2nV[i]*3);
        for(int nv = plane2nV[i]*3,j=0;j<nv;++j){
            reader>>p_plane2vertices[j];
        }


        reader>>plane2nE[i];
        reader>>val;
        auto &p_plane2edges = plane2edges[i];
        auto &p_plane2edgesMat = plane2edgesMat[i];
//        p_plane2edges.resize(plane2nE[i]*2);
//        p_plane2edgesMat.resize(plane2nE[i]*2);


        int val1,val2;
        reader>>val1;
        for(int j=1;j<plane2nE[i];++j){
            reader>>val2;
            p_plane2edges.push_back(val1);p_plane2edges.push_back(val2);
            p_plane2edgesMat.push_back(0);p_plane2edgesMat.push_back(1);

            val1 = val2;
        }

        p_plane2edges.push_back(val1);p_plane2edges.push_back(p_plane2edges[0]);
        p_plane2edgesMat.push_back(0);p_plane2edgesMat.push_back(1);

    }



    /*******************************************************/
    /*******************************************************/


    CrossSections CS;
    CS.ReadCrossSections(planeinfo,plane2vertices,plane2edges,plane2edgesMat);

    //CS.ManipulateCrossSections();

    CS.WriteCrossSectionsToFCM(filename+string("tctrloop.txt"),CS.n_materials);
    CS.WriteCrossSectionsToFCMtCtr(filename+string("tctr.txt"));


    return;





    /*******************************************************/
    /*******************************************************/



}


void TransformerSelectivePlanes(){



    //string filename("/Users/Research/Geometry/MM/ConvertFolder/chickenheart_complete6mapping.contour");
    //string fileout("/Users/Research/Geometry/MM/ConvertFolder/chickenheart_complete6mapping_subset.contour");

    string filename("/Users/Research/Geometry/MM/ConvertFolder/ctrmapping.contour");
    string fileout("/Users/Research/Geometry/MM/ConvertFolder/ctrmapping_2.contour");

    vector<int>reserveCSs({1});



    CrossSections CS;



    CS.ReadCrossSections(filename);

    CS.CalSelfProperty();

    CS.WriteCrossSections(fileout,reserveCSs);




}
