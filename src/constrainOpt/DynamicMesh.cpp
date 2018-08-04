#include <iostream>
#include <iomanip>
#include "DynamicMesh.h"

#include "mymesh/readers.h"
#include "mymesh/contour.h"
#include "mymesh/geo_curv.h"
#include <set>
#include <fstream>



void DynamicMesh::setparameters(){


    nV = vertices.size();
    nF = faces.size();
    nE = edges.size();



}


bool DynamicMesh::triangleInequality(int iface)
{
    auto &p_fe = faces[iface].fe;
    double l0 = edges[p_fe[0]].attr;
    double l1 = edges[p_fe[1]].attr;
    double l2 = edges[p_fe[2]].attr;
    if((l0+l1)>=l2 && (fabs(l0-l1)<=l2)){return true;}
    else{return false;}
}



void DynamicMesh::checkNeighborhoodtable(){


    vector<bool>checkpoint(nV,false);

    for(int i=0;i<nE;++i){
        for(auto f:edges[i].ef){
            int ll = 0;
            for(auto e:faces[f].fe)if(e==i){ll++;}
            assert(ll==1);
            for(int j=0;j<3;++j)checkpoint[faces[f].fv[j]]=true;
            for(int j=0;j<2;++j)assert(checkpoint[edges[i].ev[j]]);
            //            for(int j=0;j<2;++j)if(!checkpoint[edges[i].ev[j]]){
            //                auto ee = edges[i];
            //                auto ff = faces[f];
            //                cout<<endl;
            //            }



            for(int j=0;j<3;++j)checkpoint[faces[f].fv[j]]=false;


        }
    }

    for(int i=0;i<nF;++i){

        for(int j=0;j<3;++j)checkpoint[faces[i].fv[j]]=true;
        for(auto e:faces[i].fe){
            int ll = 0;
            for(auto f:edges[e].ef)if(f==i){ll++;}
            assert(ll==1);
            for(int j=0;j<2;++j)assert(checkpoint[edges[e].ev[j]]);
        }
        for(int j=0;j<3;++j)checkpoint[faces[i].fv[j]]=false;
    }


}

void DynamicMesh::triangleMidPoint(int tind,vector<double>&midpoint){

    auto p_fv = faces[tind].fv;
    midpoint.resize(3);
    _TriangleMidpoint(vertices[p_fv[0]].pos,vertices[p_fv[1]].pos,vertices[p_fv[2]].pos,midpoint.data());



}

void DynamicMesh::findBadTriangles(vector<int>&badtri){

    badtri.clear();
    for(int i=0;i<nF;++i)if(!triangleInequality(i))badtri.push_back(i);

}
void DynamicMesh::findBadEdges(vector<int>&badedges, double attr){

    badedges.clear();
    for(int i=0;i<nE;++i)if(edges[i].attr<attr)badedges.push_back(i);

}

void DynamicMesh::init(vector<double>&in_vertices, vector<uint>&in_faces){

    Mesh meshconverter;

    meshconverter.vertices = in_vertices;

    meshconverter.faces2vertices = in_faces;

    meshconverter.setparameters();
    meshconverter.CheckValidation();
    meshconverter.BuildNeighborTable_nonmanifold();


    nV = in_vertices.size()/3;
    nF = in_faces.size()/3;
    nE = meshconverter.n_edges;

    vertices.resize(nV);

    for(int i=0;i<nV;++i){
        Vertices &v = vertices[i];
        for(int j=0;j<3;++j)v.pos[j] = in_vertices[i*3+j];

    }

    faces.resize(nF);

    for(int i=0;i<nF;++i){
        Faces &f = faces[i];
        for(int j=0;j<3;++j)f.fv[j] = in_faces[i*3+j];
        //int nf

        //f.marker = fmarkers[i];

        f.fe.resize(3);
        auto p_fe = meshconverter.fe_begin(i);
        for(int j=0;j<3;++j)f.fe[j] = p_fe[j];

    }


    edges.resize(nE);
    for(int i=0;i<nE;++i){
        Edges &e = edges[i];
        auto p_ev = meshconverter.ev_begin(i);
        for(int j=0;j<2;++j)e.ev[j] = p_ev[j];

        int n_ef = meshconverter.efnon_num(i);
        auto p_ef = meshconverter.efnon_begin(i);

        e.ef.resize(n_ef);
        for(int j=0;j<n_ef;++j)e.ef[j] = p_ef[j];
    }




}

void DynamicMesh::init(vector<double> &in_vertices, vector<uint> &in_faces, vector<int> &fmarkers, map<vector<int>, int> &edgemarkers){


    Mesh meshconverter;

    meshconverter.vertices = in_vertices;

    meshconverter.faces2vertices = in_faces;

    meshconverter.setparameters();
    meshconverter.CheckValidation();
    meshconverter.BuildNeighborTable_nonmanifold();


    nV = in_vertices.size()/3;
    nF = in_faces.size()/3;
    nE = meshconverter.n_edges;

    vertices.resize(nV);

    for(int i=0;i<nV;++i){
        Vertices &v = vertices[i];
        for(int j=0;j<3;++j)v.pos[j] = in_vertices[i*3+j];

    }

    faces.resize(nF);

    for(int i=0;i<nF;++i){
        Faces &f = faces[i];
        for(int j=0;j<3;++j)f.fv[j] = in_faces[i*3+j];
        //int nf

        f.marker = fmarkers[i];

        f.fe.resize(3);
        auto p_fe = meshconverter.fe_begin(i);
        for(int j=0;j<3;++j)f.fe[j] = p_fe[j];

    }


    edges.resize(nE);
    for(int i=0;i<nE;++i){
        Edges &e = edges[i];
        auto p_ev = meshconverter.ev_begin(i);
        for(int j=0;j<2;++j)e.ev[j] = p_ev[j];

        int n_ef = meshconverter.efnon_num(i);
        auto p_ef = meshconverter.efnon_begin(i);

        e.ef.resize(n_ef);
        for(int j=0;j<n_ef;++j)e.ef[j] = p_ef[j];
    }


    propagateFaceMarkertoEdges();


    edge_markerset.clear();
    vector<int>nonedges;
    for(int i=0;i<nE;++i)if(edges[i].ef.size()>2){
        vector<int>ind(2);
        ind[0] = edges[i].ev[0];ind[1] = edges[i].ev[1];
        sort(ind.begin(),ind.end());

        nonedges.push_back(ind[0]);nonedges.push_back(ind[1]);
        auto iter = edgemarkers.find(ind);

        if(iter != edgemarkers.end()){
            edges[i].marker = iter->second;
            edge_markerset.insert(iter->second);
        }else{
            // auto &ees = edges[i];

            //            int key1 = 128,key2 = 129;
            //            vector<Faces>ffss;
            //            for(int k=0;k<nF;++k){
            //                for(int j=0;j<3;++j)if(faces[k].fv[j]==key1 ||faces[k].fv[j]==key2){
            //                    ffss.push_back(faces[k]);
            //                }
            //            }

            cout<<"Do not Find"<<endl;
        }


    }

    nEdgemarker = edge_markerset.size();
    markerset.clear();
    for(auto a:fmarkers)markerset.insert(a);
    nMarker = markerset.size();
    cout<<edgemarkers.size()<<endl;
    cout<<endl;
}

void DynamicMesh::reinit(){

    vector<double>vpos(nV*3);
    for(int i=0;i<nV;++i){
        auto p_v = vertices[i].pos;
        for(int j=0;j<3;++j)vpos[i*3+j] = p_v[j];

    }

    vector<int>fmarker(nF);
    vector<unsigned int>faces2vertices(nF*3);
    for(int i=0;i<nF;++i){
        auto p_fv = faces[i].fv;
        for(int j=0;j<3;++j)faces2vertices[i*3+j] = p_fv[j];
        fmarker[i] = faces[i].marker;
    }
    //void init(vector<double>&in_vertices, vector<uint>&in_faces, vector<int>&fmarkers, map<vector<int>, int> &edgemarkers);
    map<vector<int>, int> edgemarkers;
    for(int i=0;i<nE;++i)if(edges[i].ef.size()>2){

        auto p_ev = edges[i].ev;
        vector<int>eee;
        eee.push_back(p_ev[0]);eee.push_back(p_ev[1]);
        sort(eee.begin(),eee.end());
        edgemarkers[eee] = edges[i].marker;


    }
    init(vpos,faces2vertices,fmarker,edgemarkers);
}


void DynamicMesh::outputToInside(vector<double>&out_vertices,vector<uint>&out_faces2v){


    out_vertices.resize(nV*3);
    for(int i=0;i<nV;++i){
        auto p_v = vertices[i].pos;
        for(int j=0;j<3;++j)out_vertices[i*3+j] = p_v[j];

    }


    out_faces2v.resize(nF*3);
    for(int i=0;i<nF;++i){
        auto p_fv = faces[i].fv;
        for(int j=0;j<3;++j)out_faces2v[i*3+j] = p_fv[j];
    }


}


void DynamicMesh::splitEdges(int iedges, int &newvertex, vector<int> &newedges){


    newedges.clear();
    //important
    edges.reserve(edges.size()+edges[iedges].ef.size()*2);
    faces.reserve(faces.size()+edges[iedges].ef.size()*2);

    auto Estrutcture = edges[iedges];
    int nef = edges[iedges].ef.size();
    vector<int>p_ev(2);
    for(int i=0;i<2;++i)p_ev[i] = edges[iedges].ev[i];
    vector<int>topvertice(nef);
    vector<vector<int>>newface;


    newvertex = vertices.size();

    vertices.resize(newvertex+1);
    MyUtility::_VerticesMidpoint(vertices[p_ev[0]].pos,vertices[p_ev[1]].pos,vertices[newvertex].pos);

    vector<int>newsplitee(2);
    newsplitee[0] = iedges;newsplitee[1]  = edges.size();
    edges.resize(edges.size()+1);
    newedges.push_back(newsplitee[0]);newedges.push_back(newsplitee[1]);

    edges[newsplitee[0]].ev[0] = p_ev[0];edges[newsplitee[0]].ev[1] = newvertex;
    edges[newsplitee[1]].ev[0] = p_ev[1];edges[newsplitee[1]].ev[1] = newvertex;

    edges[newsplitee[0]].marker = edges[newsplitee[1]].marker = edges[iedges].marker;


    //edges[newsplitee[0]].ef = Estrutcture.ef;

    edges[newsplitee[1]].ef.resize(nef);
    for(int i=0;i<nef;++i){
        edges[newsplitee[1]].ef[i] = nF+i;
    }

    for(int i=0;i<nef;++i){
        int iface = Estrutcture.ef[i];

        auto p_fv = faces[iface].fv;
        auto &p_fe = faces[iface].fe;
        int e1 = -1,e2 = -1;
        for(int j=0;j<3;++j)if(p_fe[j]!=iedges){
            auto evev = edges[p_fe[j]].ev;
            bool istop = false;
            for(int k=0;k<2;++k)if(evev[k]==p_ev[0])istop=true;
            if(istop){assert(e1==-1);e1 = p_fe[j];}
            else{assert(e2==-1);e2 = p_fe[j];}

        }
        assert(e1!=-1);assert(e2!=-1);
        bool test = false;
        for(int j=0;j<edges[e2].ef.size();++j)if(edges[e2].ef[j]==iface){edges[e2].ef[j] = faces.size();test=true;break;}
        if(!test){
            auto eeeee = edges[e2];
            assert(test);
        }
        assert(test);

        for(int j=0;j<3;++j)if(p_fv[j]!=p_ev[0] && p_fv[j]!=p_ev[1]){
            topvertice[i] = p_fv[j];
            break;
        }
        int inewedge = edges.size();
        edges.resize(inewedge+1);
        edges[inewedge].ev[0] = topvertice[i];
        edges[inewedge].ev[1] = newvertex;

        edges[inewedge].ef.resize(2);
        edges[inewedge].ef[0] = iface;
        edges[inewedge].ef[1] = faces.size();

        newedges.push_back(inewedge);

        int newfaceind = faces.size();
        faces.resize(newfaceind+1);
        faces[newfaceind] = faces[iface];

        bool test1 = false;
        for(int j=0;j<3;++j)if(faces[iface].fv[j]==p_ev[1]){
            faces[iface].fv[j] = newvertex;test1=true;
            break;
        }
        assert(test1);

        bool test2 = false;
        for(int j=0;j<3;++j)if(faces[iface].fe[j]==e2){
            faces[iface].fe[j] = inewedge;test2=true;
            break;
        }
        assert(test2);






        bool test3 = false;
        for(int j=0;j<3;++j)if(faces[newfaceind].fv[j]==p_ev[0]){
            faces[newfaceind].fv[j] = newvertex;test3=true;
            break;
        }
        assert(test3);
        //faces[newfaceind].fv[j] = newvertex;


        faces[newfaceind].fe.resize(3);
        faces[newfaceind].fe[0] = e2;faces[newfaceind].fe[1] = newsplitee[1];faces[newfaceind].fe[2] =inewedge;



        faces[newfaceind].marker  =  faces[iface].marker;


        propagateFaceMarkertoEdges(inewedge);
        //checkNeighborhoodtable();
    }


    setparameters();
    //checkNeighborhoodtable();







}


void DynamicMesh::propagateFaceMarkertoEdges(int iedges){

    auto &ef = edges[iedges].ef;
    if(ef.size()==2){
        if(faces[ef[0]].marker != faces[ef[1]].marker){
            cout<<"edge marker error"<<endl;
        }else{

            edges[iedges].marker = faces[ef[0]].marker;
        }
    }
    else if(ef.size()==1){
        edges[iedges].marker = faces[ef[0]].marker;

    }else if(ef.size()>2){
        edges[iedges].marker = NONMANIFOLDBEGIN;
    }else{
        cout<<"isolate edge error"<<endl;

    }

}

void DynamicMesh::propagateFaceMarkertoEdges(){

    for(int i=0;i<nE;++i){
        propagateFaceMarkertoEdges(i);
    }

}

void DynamicMesh::writeObjfile(string filename){


    //bool writeObjFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices);
    vector<double>vpos(nV*3);
    for(int i=0;i<nV;++i){
        auto p_v = vertices[i].pos;
        for(int j=0;j<3;++j)vpos[i*3+j] = p_v[j];

    }

    vector<unsigned int>faces2vertices(nF*3);
    for(int i=0;i<nF;++i){
        auto p_fv = faces[i].fv;
        for(int j=0;j<3;++j)faces2vertices[i*3+j] = p_fv[j];

    }

    //    Mesh meshconverter;

    //    meshconverter.vertices = vpos;

    //    meshconverter.faces2vertices = faces2vertices;

    //    meshconverter.setparameters();
    //    meshconverter.BuildNeighborTable_nonmanifold();
    //    meshconverter.CheckValidation();
    writeObjFile( filename,vpos,faces2vertices);

}
void DynamicMesh::outputToMathematica(string fpath){

    string vfname = fpath + string("toMathematicaV.txt");
    ofstream ofs(vfname, std::ofstream::out);
    ofs << std::fixed;
    ofs << setprecision(13);
    ofs << "{";
    for(int i=0;i<nV;++i){
        ofs << "{";
        for(int j=0;j<2;++j)ofs<<vertices[i].pos[j]<<",";
        ofs<<vertices[i].pos[2];
        if(i!=nV-1)ofs<<"},";
        else ofs<<"}";
    }
    ofs << "}";
    ofs.close();


    string efname = fpath + string("toMathematicaE.txt");
    ofs.open(efname, std::ofstream::out);
    int nP = nMarker;
    ofs<<"{";
    for(int i=0;i<nP;++i){
        vector<bool>pickE(nE,false);
        vector<int>pickEbuffer(0);
        for(auto &f:faces)if(f.marker==i)for(auto e:f.fe)pickE[e] = true;
        for(int j=0;j<nE;++j)if(pickE[j])pickEbuffer.push_back(j);
        ofs<<"{";
            int nwE = pickEbuffer.size();
            for(int j=0;j<nwE;++j){
                auto &eStruc = edges[pickEbuffer[j]];
                ofs<<"{";
                ofs<<eStruc.ev[0]+1<<","<<eStruc.ev[1]+1;
                if(j!=nwE-1)ofs<<"},";
                else ofs<<"}";
            }


        if(i!=nP-1)ofs<<"},";
        else ofs<<"}";
    }
    ofs<<"}";
    ofs.close();


    vector< vector<int> >vertices2Plane(nV);
    for(int i=0;i<nP;++i){
        vector<bool>pickV(nV,false);
        for(auto &f:faces)if(f.marker==i)for(auto v:f.fv)pickV[v] = true;
        for(int j=0;j<nV;++j)if(pickV[j])vertices2Plane[j].push_back(i+1);
    }
    string v2pfname = fpath + string("toMathematicaV2P.txt");
    ofs.open(v2pfname, std::ofstream::out);
    ofs<<"{";
        for(int i=0;i<nV;++i){
            ofs<<"{";
            for(int j=0;j<vertices2Plane[i].size();++j){
                ofs<<vertices2Plane[i][j];
                if(j!=vertices2Plane[i].size()-1)ofs<<",";
            }
            if(i!=nV-1)ofs<<"},";
            else ofs<<"}";

        }

    ofs<<"}";
    ofs.close();


    string lfname = fpath + string("toMathematicaL.txt");
    ofs.open(lfname, std::ofstream::out);
    int nL = nEdgemarker;
    ofs<<"{";
    for(int i=0;i<nL;++i){
        vector<bool>pickV(nV,false);
        vector<int>pickVbuffer(0);
        for(auto &e:edges)if(e.marker==NONMANIFOLDBEGIN+i)for(auto v:e.ev)pickV[v] = true;
        for(int j=0;j<nV;++j)if(pickV[j])pickVbuffer.push_back(j+1);
        ofs<<"{";
            int nwV = pickVbuffer.size();
            for(int j=0;j<nwV;++j){
                ofs<<pickVbuffer[j];
                if(j!=nwV-1)ofs<<",";
            }


        if(i!=nL-1)ofs<<"},";
        else ofs<<"}";
    }
    ofs<<"}";
    ofs.close();



    string lefname = fpath + string("toMathematicaLE.txt");
    ofs.open(lefname, std::ofstream::out);
    ofs<<"{";
    for(int i=0;i<nL;++i){
        vector<bool>pickE(nE,false);
        vector<int>pickEbuffer(0);
        for(int j=0;j<nE;++j)if(edges[j].marker==NONMANIFOLDBEGIN+i)pickE[j] = true;
        for(int j=0;j<nE;++j)if(pickE[j])pickEbuffer.push_back(j);
        ofs<<"{";
            int nwE = pickEbuffer.size();
            for(int j=0;j<nwE;++j){
                auto &eStruc = edges[pickEbuffer[j]];
                ofs<<"{";
                ofs<<eStruc.ev[0]+1<<","<<eStruc.ev[1]+1;
                if(j!=nwE-1)ofs<<"},";
                else ofs<<"}";
            }


        if(i!=nL-1)ofs<<"},";
        else ofs<<"}";
    }
    ofs<<"}";
    ofs.close();

    string tfname = fpath + string("toMathematicaT.txt");
    ofs.open(tfname, std::ofstream::out);
    ofs<<"{";
    for(int i=0;i<nP;++i){
        vector<bool>pickF(nF,false);
        vector<int>pickFbuffer(0);
        for(int j=0;j<nF;++j)if(faces[j].marker==i)pickF[j] = true;
        for(int j=0;j<nF;++j)if(pickF[j])pickFbuffer.push_back(j);
        ofs<<"{";
            int nwF = pickFbuffer.size();
            for(int j=0;j<nwF;++j){
                auto &fStruc = faces[pickFbuffer[j]];
                ofs<<"{";
                ofs<<fStruc.fv[0]+1<<","<<fStruc.fv[1]+1<<","<<fStruc.fv[2]+1;
                if(j!=nwF-1)ofs<<"},";
                else ofs<<"}";
            }


        if(i!=nP-1)ofs<<"},";
        else ofs<<"}";
    }
    ofs<<"}";
    ofs.close();



}

void DynamicMesh::outputToMarker(vector<onePlane>&planeMeshes){


    planeMeshes.clear();

    planeMeshes.resize(nMarker);

    vector<int>mvec;
    for(auto a:markerset)mvec.push_back(a);
    for(int imarker=0;imarker<nMarker;++imarker){

        int marker = mvec[imarker];
        onePlane &iplane = planeMeshes[imarker];
        iplane.marker = marker;

        vector<int>pickV(nV,-1);
        vector<bool>pickT(nF,false);

        int nnewF = 0;
        for(int i=0;i<nF;++i)if(faces[i].marker == marker){
            pickT[i]=true;
            for(int j=0;j<3;++j)pickV[faces[i].fv[j]] = 0;
            nnewF++;
        }
        int nnewV = 0;
        for(auto &a:pickV)if(a!=-1)a=nnewV++;

        auto &ver = iplane.vertices;
        auto &vattr = iplane.vattr;
        auto &invver = iplane.inverseVmapping;
        for(int i = 0;i<nV;++i)if(pickV[i]!=-1){
            auto p_v = vertices[i].pos;
            for(int j=0;j<3;++j)ver.push_back(p_v[j]);
            vattr.push_back(vertices[i].attr);
            invver.push_back(i);
        }

        auto &fac = iplane.faces2vertices;
        auto &invfac = iplane.inverseFmapping;
        for(int i = 0;i<nF;++i)if(pickT[i]){
            auto &p_fv = faces[i].fv;
            for(int j=0;j<3;++j)fac.push_back(pickV[p_fv[j]]);
            invfac.push_back(i);
        }


    }


}




void DynamicMesh::Contouring_EdgeMidpoint(vector<vector<double>>&val,double *planepara,vector<double>&outCtrV,vector<uint>&outCtrE,vector<int>&outCtrEM){




    cout<<"Contouring_EdgeMidpoint"<<endl;
    assert(val.size()==nV);


    outCtrV.clear();outCtrE.clear();outCtrEM.clear();


    vector<int>vlabel(nV);
    for(int i=0;i<nV;++i)vlabel[i] = max_element(val[i].begin(),val[i].end()) - val[i].begin();

    vector<int>activeEdge(nE,-1);
    vector<int>activeEdgebuffer;
    vector<int>activeFace(nF,0);
    int ncountAE = 0;
    for(int i=0;i<nE;++i)if(vlabel[edges[i].ev[0]]!=vlabel[edges[i].ev[1]]){
        activeEdge[i]=ncountAE++;activeEdgebuffer.push_back(i);
        for(auto a:edges[i].ef){
            activeFace[a]++;
        }
    }

    vector<int>activeFacebuffer;
    for(int i=0;i<nF;++i)if(activeFace[i]!=0)activeFacebuffer.push_back(i);


    vector<vector<double>>edgewmidpoint(ncountAE);
    for(int i=0;i<activeEdgebuffer.size();++i){
        int eind = activeEdgebuffer[i];
        int v0 = edges[eind].ev[0],v1=edges[eind].ev[1];

        vector<double>v0value = val[v0];
        vector<double>v1value = val[v1];

        double weighted0,weighted1;


        int l0 = vlabel[v0],l1 = vlabel[v1];
        weighted0 =v0value[l0]-v0value[l1];
        weighted1 =v1value[l1]-v1value[l0];

        double sumweight = (weighted0+weighted1);
        vector<double>nvpos(3);

        weightedAddVec(weighted1/sumweight,weighted0/sumweight,vertices[v0].pos,vertices[v1].pos,nvpos.data());

        edgewmidpoint[i] = nvpos;

        for(auto a:nvpos)outCtrV.push_back(a);

    }


    vector<bool>vflagbuf(nV,false);
    for(int iff=0;iff<activeFacebuffer.size();++iff){

        int find = activeFacebuffer[iff];
        assert(activeFace[find]>=2);assert(activeFace[find]<=3);
        for(int j=0;j<3;++j)vflagbuf[faces[find].fv[j]] = true;
        if(activeFace[find]==2){

            int topv,bottomv[2],oute = -1;
            vector<int>acee;
            for(auto e:faces[find].fe)if(activeEdge[e]!=-1){
                acee.push_back(activeEdge[e]);
            }else{
                oute = e;
                for(int i=0;i<2;++i)bottomv[i] = edges[oute].ev[i];
            }
            assert(oute!=-1);
            assert(vlabel[bottomv[0]]==vlabel[bottomv[1]]);
            for(int i=0;i<2;++i)vflagbuf[bottomv[i] ]= false;
            for(int j=0;j<3;++j)if(vflagbuf[faces[find].fv[j]]){topv = faces[find].fv[j];break;}

            int topvlabel = vlabel[topv],bvlabel = vlabel[bottomv[0]];
            for(auto a:acee)outCtrE.push_back(a);
            double edgevec[3],ev12topvvec[3];
            minusVec(edgewmidpoint[acee[1]].data(),edgewmidpoint[acee[0]].data(),edgevec);
            minusVec(vertices[topv].pos,edgewmidpoint[acee[0]].data(),ev12topvvec);
            double crsVec[3];
            cross(edgevec,ev12topvvec,crsVec);
            if(dot(crsVec,planepara)>0){
                outCtrEM.push_back(topvlabel);
                outCtrEM.push_back(bvlabel);
            }else{
                outCtrEM.push_back(bvlabel);
                outCtrEM.push_back(topvlabel);
            }


        }else if(activeFace[find]==3){

            int  midpointInd = outCtrV.size()/3;

            vector<double>trimid(3);

            triangleMidPoint(find,trimid);
            for(auto a:trimid)outCtrV.push_back(a);
            auto p_fe = faces[find].fe;
            for(int i=0;i<3;++i){
                int eind = p_fe[i];
                outCtrE.push_back(activeEdge[eind]);outCtrE.push_back(midpointInd);
                int ev0 = edges[eind].ev[0],ev1 = edges[eind].ev[1];
                double edgevec[3],ev12topvvec[3];
                minusVec(trimid.data(),edgewmidpoint[activeEdge[eind]].data(),edgevec);
                minusVec(vertices[ev0].pos,edgewmidpoint[activeEdge[eind]].data(),ev12topvvec);

                double crsVec[3];
                cross(edgevec,ev12topvvec,crsVec);
                if(dot(crsVec,planepara)>0){
                    outCtrEM.push_back(vlabel[ev0]);
                    outCtrEM.push_back(vlabel[ev1]);
                }else{
                    outCtrEM.push_back(vlabel[ev1]);
                    outCtrEM.push_back(vlabel[ev0]);
                }

            }




        }
        for(int j=0;j<3;++j)vflagbuf[faces[find].fv[j]] = false;


    }
    cout<<"Contouring_EdgeMidpoint end"<<endl;

}


void DynamicMesh::ContouringBatch_EdgeMidpoint(vector<vector<double>>&val,vector<vector<double>>&planepara,ContourStruct &ccstruct){




    cout<<"ContouringBatch_EdgeMidpoint"<<endl;
    assert(val.size()==nV);


    auto &outCtrV = ccstruct.outCtrV;
    auto &outCtrE = ccstruct.outCtrE;
    auto &outCtrEM = ccstruct.outCtrEM;
    auto &outCtrEMarker = ccstruct.outCtrEMarker;
    outCtrV.clear();outCtrE.clear();outCtrEM.clear();outCtrEMarker.clear();


    vector<int>vlabel(nV);
    for(int i=0;i<nV;++i)vlabel[i] = max_element(val[i].begin(),val[i].end()) - val[i].begin();

    vector<int>activeEdge(nE,-1);
    vector<int>activeEdgebuffer;
    vector<int>activeFace(nF,0);
    int ncountAE = 0;
    for(int i=0;i<nE;++i)if(vlabel[edges[i].ev[0]]!=vlabel[edges[i].ev[1]]){
        activeEdge[i]=ncountAE++;activeEdgebuffer.push_back(i);
        for(auto a:edges[i].ef){
            activeFace[a]++;
        }
    }

    vector<int>activeFacebuffer;
    for(int i=0;i<nF;++i)if(activeFace[i]!=0)activeFacebuffer.push_back(i);


    vector<vector<double>>edgewmidpoint(ncountAE);
    vector<double>edgew(ncountAE);
    for(int i=0;i<activeEdgebuffer.size();++i){
        int eind = activeEdgebuffer[i];
        int v0 = edges[eind].ev[0],v1=edges[eind].ev[1];

        vector<double>v0value = val[v0];
        vector<double>v1value = val[v1];

        double weighted0,weighted1;


        int l0 = vlabel[v0],l1 = vlabel[v1];
        weighted0 =v0value[l0]-v0value[l1];
        weighted1 =v1value[l1]-v1value[l0];

        double sumweight = (weighted0+weighted1);
        edgew[i] = weighted1/min(1e-8,weighted0);
        vector<double>nvpos(3);

        weightedAddVec(weighted1/sumweight,weighted0/sumweight,vertices[v0].pos,vertices[v1].pos,nvpos.data());

        edgewmidpoint[i] = nvpos;

        for(auto a:nvpos)outCtrV.push_back(a);

    }


    vector<bool>vflagbuf(nV,false);
    for(int iff=0;iff<activeFacebuffer.size();++iff){

        int find = activeFacebuffer[iff];
        assert(activeFace[find]>=2);assert(activeFace[find]<=3);
        for(int j=0;j<3;++j)vflagbuf[faces[find].fv[j]] = true;
        if(activeFace[find]==2){

            int topv,bottomv[2],oute = -1;
            vector<int>acee;
            for(auto e:faces[find].fe)if(activeEdge[e]!=-1){
                acee.push_back(activeEdge[e]);
            }else{
                oute = e;
                for(int i=0;i<2;++i)bottomv[i] = edges[oute].ev[i];
            }
            assert(oute!=-1);
            assert(vlabel[bottomv[0]]==vlabel[bottomv[1]]);
            for(int i=0;i<2;++i)vflagbuf[bottomv[i] ]= false;
            for(int j=0;j<3;++j)if(vflagbuf[faces[find].fv[j]]){topv = faces[find].fv[j];break;}

            int topvlabel = vlabel[topv],bvlabel = vlabel[bottomv[0]];
            for(auto a:acee)outCtrE.push_back(a);
            outCtrEMarker.push_back(faces[find].marker);
            double edgevec[3],ev12topvvec[3];
            minusVec(edgewmidpoint[acee[1]].data(),edgewmidpoint[acee[0]].data(),edgevec);
            minusVec(vertices[topv].pos,edgewmidpoint[acee[0]].data(),ev12topvvec);
            double crsVec[3];
            cross(edgevec,ev12topvvec,crsVec);
            if(dot(crsVec,planepara[faces[find].marker].data())>0){
                outCtrEM.push_back(topvlabel);
                outCtrEM.push_back(bvlabel);
            }else{
                outCtrEM.push_back(bvlabel);
                outCtrEM.push_back(topvlabel);
            }


        }else if(activeFace[find]==3){

            int  midpointInd = outCtrV.size()/3;

            vector<double>trimid(3,0);

            triangleMidPoint(find,trimid);

//            for(auto e:faces[find].fe){
//                auto p_mide  = edgewmidpoint[activeEdge[e]].data();
//                for(int kk=0;kk<3;++kk)trimid[kk]+=p_mide[kk];
//            }
//            for(int kk=0;kk<3;++kk)trimid[kk]/=3;

//            vector<vector<double>>weightlist(3,vector<double>(3));
//            auto p_fv = faces[find].fv;
//            for(auto eind:faces[find].fe){
//                 int v0 = edges[eind].ev[0],v1=edges[eind].ev[1];
//                 int ind0,ind1;
//                 double ww = edgew[activeEdge[eind]];
//                 for(int i=0;i<3;++i)if(p_fv[i]==v0){ind0=i;break;}
//                 for(int i=0;i<3;++i)if(p_fv[i]==v1){ind1=i;break;}
//                 weightlist[ind0][ind1] = min(1e-8,ww);
//                 weightlist[ind1][ind0] = 1/min(1e-8,ww);
//            }

//            double y = weightlist[0][2];double yx = y/weightlist[0][1];
//            double c = 1/(1+y+yx), a = y*c, b = yx*c;
//            vector<double*>p_fvpos(3);
//            cout<<a<<' '<<b<<' '<<c<<' '<<endl;
//            for(int i=0;i<3;++i)p_fvpos[i] = vertices[p_fv[i]].pos;
//            for(int i=0;i<3;++i)trimid[i] = a * p_fvpos[0][i] +  b * p_fvpos[1][i] +  c * p_fvpos[2][i];


//            auto p_fv = faces[find].fv;
//            vector<int>vlabelset;
//            for(int i=0;i<3;++i)vlabelset.push_back(vlabel[p_fv[i]]);
//            vector<double>vweight(3,0);
//            for(int i=0;i<3;++i){
//                  for(int j=0;j<3;++j)vweight[i] += fabs(val[p_fv[i]][vlabelset[j]]);
//                  for(int j=0;j<3;++j)cout<<val[p_fv[i]][vlabelset[j]]<<' ';cout<<endl;
//            }
//            cout<<endl;
//            //for(int i=0;i<3;++i)cout<<vweight[i]<<' ';cout<<endl;
//            double sumvw = 0;
//            for(auto &a:vweight){a=1/a;sumvw+=a;}
//            for(auto &a:vweight)a/=sumvw;
//            //for(int i=0;i<3;++i)cout<<vweight[i]<<' ';cout<<endl;
//            vector<double*>p_fvpos(3);
//            for(int i=0;i<3;++i)p_fvpos[i] = vertices[p_fv[i]].pos;
//            for(int i=0;i<3;++i)trimid[i] = vweight[0] * p_fvpos[0][i] +  vweight[1] * p_fvpos[1][i] + vweight[2]* p_fvpos[2][i];








            for(auto a:trimid)outCtrV.push_back(a);
            auto p_fe = faces[find].fe;
            for(int i=0;i<3;++i){
                int eind = p_fe[i];
                outCtrE.push_back(activeEdge[eind]);outCtrE.push_back(midpointInd);
                outCtrEMarker.push_back(faces[find].marker);
                int ev0 = edges[eind].ev[0],ev1 = edges[eind].ev[1];
                double edgevec[3],ev12topvvec[3];
                minusVec(trimid.data(),edgewmidpoint[activeEdge[eind]].data(),edgevec);
                minusVec(vertices[ev0].pos,edgewmidpoint[activeEdge[eind]].data(),ev12topvvec);

                double crsVec[3];
                cross(edgevec,ev12topvvec,crsVec);
                if(dot(crsVec,planepara[faces[find].marker].data())>0){
                    outCtrEM.push_back(vlabel[ev0]);
                    outCtrEM.push_back(vlabel[ev1]);
                }else{
                    outCtrEM.push_back(vlabel[ev1]);
                    outCtrEM.push_back(vlabel[ev0]);
                }

            }




        }
        for(int j=0;j<3;++j)vflagbuf[faces[find].fv[j]] = false;


    }
    cout<<"ContouringBatch_EdgeMidpoint end"<<endl;

}


void DynamicMesh::separateContourbyMarker(ContourStruct &ccstruct,vector<ContourStruct>&cccSs){



    cccSs.clear();

    cccSs.resize(nMarker);

    int nCE = ccstruct.outCtrE.size()/2;
    int nCV = ccstruct.outCtrV.size()/3;

    vector<int>mvec;
    for(auto a:markerset)mvec.push_back(a);
    for(int imarker=0;imarker<nMarker;++imarker){

        int marker = mvec[imarker];
        ContourStruct &iplane = cccSs[imarker];
        //iplane.marker = marker;

        vector<int>pickV(nCV,-1);
        vector<bool>pickE(nCE,false);

        int nnewE = 0;
        for(int i=0;i<nCE;++i)if(ccstruct.outCtrEMarker[i] == marker){
            pickE[i]=true;
            for(int j=0;j<2;++j)pickV[ccstruct.outCtrE[i*2+j]] = 0;
            nnewE++;
        }
        int nnewV = 0;
        for(auto &a:pickV)if(a!=-1)a=nnewV++;

        auto &ver = iplane.outCtrV;
        //auto &vattr = iplane.vattr;
        //auto &invver = iplane.inverseVmapping;
        for(int i = 0;i<nCV;++i)if(pickV[i]!=-1){
            auto p_v = ccstruct.outCtrV.data()+i*3;
            for(int j=0;j<3;++j)ver.push_back(p_v[j]);
            //vattr.push_back(vertices[i].attr);
            //invver.push_back(i);
        }

        auto &edg = iplane.outCtrE;
        auto &edgM = iplane.outCtrEM;
        //auto &invfac = iplane.inverseFmapping;
        for(int i = 0;i<nCE;++i)if(pickE[i]){
            auto p_ev = ccstruct.outCtrE.data()+i*2;
            for(int j=0;j<2;++j)edg.push_back(pickV[p_ev[j]]);
            for(int j=0;j<2;++j)edgM.push_back(ccstruct.outCtrEM[i*2+j]);
            //invfac.push_back(i);
        }


    }


}


void DynamicMesh::ContourStruct::speicalmarkersmoothing(int maxiter){

    int (n_vertices) = outCtrV.size()/3;
    int n_edges = outCtrEMarker.size();
    assert(outCtrEMarker.size()==outCtrE.size()/2 && outCtrEM.size()/2==outCtrEMarker.size());
    vector<vector<int>>vnn(n_vertices);


    set<int>markset;
    for(auto a:outCtrEMarker)markset.insert(a);
    int nnp = markset.size();
    vector<int>vmarker(n_vertices);


    for(int i=0;i<n_edges;++i){
        int v1 = outCtrE[i*2],v2 = outCtrE[i*2+1];
        vnn[v1].push_back(v2);
        vnn[v2].push_back(v1);
        vmarker[v1] = outCtrEMarker[i];vmarker[v2] = outCtrEMarker[i];
    }



    vector<double>buffer1(outCtrV.size());
    vector<double>buffer2(outCtrV.size());
    vector<double>*pre = &buffer1,*cur = &buffer2;
    *pre = outCtrV;

    double kaipa = 0.1,lamda = 0.5;
    auto gama = 1 / (kaipa - 1 / lamda);

    auto oneiter = [&vnn,n_vertices,this,&vmarker,nnp](double coef,vector<double>*pre,vector<double>*cur){
        double ocoef = 1- coef;

        vector<bool>isspecial(n_vertices,false);
        vector<double>normalupdateVec(n_vertices*3);
        vector<double>specialupdateVec(n_vertices*3,0);
        vector<double>aaa(3);
        for(int i =0;i<n_vertices;++i){

            for(int j=0;j<3;++j)aaa[j] = 0;
            auto &vv = vnn[i];
            if(vv.size()!=2){
                //cout<<"Adas"<<endl;
                for(int j=0;j<3;++j)cur->at(i*3+j) = pre->at(i*3+j);
                for(int j=0;j<3;++j)normalupdateVec[i*3+j]=0;
                isspecial[i]=true;
                continue;
                if(vv.size()!=4)continue;
                vector<vector<int>>mapper(nnp);
                for(auto a:vv)mapper[vmarker[a]].push_back(a);

                for(auto &a:mapper)if(a.size()!=0){
                    assert(a.size()==2);
                    vector<double>mp(3,0);
                    for(int k=0;k<3;++k)mp[k] = (pre->at((a[0]*3)+k)+pre->at((a[1]*3)+k))/2;
                    for(int k=0;k<3;++k)mp[k] = pre->at((i*3)+k)-mp[k];

                    for(int k=0;k<3;++k)specialupdateVec[a[0]*3+k] = 0.3*mp[k];
                    for(int k=0;k<3;++k)specialupdateVec[a[1]*3+k] = 0.3*mp[k];
                }




                continue;
            }

            for(int k=0;k<vv.size();++k){
                for(int j=0;j<3;++j){aaa[j]+=pre->at((vv[k]*3)+j);}
            }
            for(int j=0;j<3;++j)aaa[j]/=(vv.size());
            for(int j=0;j<3;++j)normalupdateVec[i*3+j] = coef * aaa[j];
            //for(int j=0;j<3;++j)cout<<aaa[j]<<" ";cout<<endl;
            //for(int j=0;j<3;++j)cur->at(i*3+j) = ocoef * pre->at(i*3+j) + coef * aaa[j];
        }
        for(int i =0;i<n_vertices;++i)if(!isspecial[i])
            for(int j=0;j<3;++j)cur->at(i*3+j) = ocoef * pre->at(i*3+j) + normalupdateVec[i*3+j]+specialupdateVec[i*3+j];

    };

    for(int iter=0;iter<maxiter;++iter){

        //cout<<iter<<endl;
        oneiter(lamda,pre,cur);
        swap(pre,cur);
        oneiter(gama,pre,cur);
        swap(pre,cur);

    }
    outCtrV=(*pre);

}



void DynamicMesh::ContourStruct::CurveNetReUniformSampling(int s){


    n_rf::Curve curveContainer;

    curveContainer.ImportCurve(this->outCtrV,this->outCtrE);

    vector<int>doubleCrossInd2FakedCellInd(this->outCtrEM.size());
    set<int>Markerset;
    for(int i=0;i<outCtrEMarker.size();++i){
        int fakeInd = outCtrEMarker[i];
        Markerset.insert(fakeInd);
        doubleCrossInd2FakedCellInd[i*2] = fakeInd;
    }
    int n_markers = Markerset.size();
    for(int i=0;i<outCtrEMarker.size();++i){
        doubleCrossInd2FakedCellInd[i*2+1] = outCtrEMarker[i] + n_markers;
    }


    vector<int>returnCrossInd2FakedCellInd;

    curveContainer.CurveNetAnalayze(this->outCtrEM,doubleCrossInd2FakedCellInd,n_markers*2);

    curveContainer.CurveNetUniformSampling(0.5,outCtrV,outCtrE,outCtrEM,returnCrossInd2FakedCellInd);

    outCtrEMarker.resize(outCtrE.size()/2);
    for(int i=0;i<outCtrEM.size()/2;++i){
        if(returnCrossInd2FakedCellInd[i*2]<n_markers)outCtrEMarker[i] = returnCrossInd2FakedCellInd[i*2];
        else outCtrEMarker[i] = returnCrossInd2FakedCellInd[i*2+1];

    }








}



void CtrManipulation(){


    CrossSections Cs;

    string filename("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/visualization2/tctrContourmapping.contour");
    string ofilename("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/visualization2/tctrContourmapping2.contour");

    Cs.ReadCrossSections(filename);


    vector<int>pickedCSs({4,5});
    Cs.WriteCrossSections(ofilename,pickedCSs);



}
