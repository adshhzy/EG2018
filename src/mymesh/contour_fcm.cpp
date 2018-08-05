#include "contour.h"
#include<iostream>
#include<fstream>
#include<set>
#include<assert.h>
#include<map>
#include <iomanip>
#include "readers.h"
#include "my_mesh.h"
#include <unistd.h>



double planeIntersectLine(planeinfo &p,lineinfo &l,vector<double>&intersectionV);
bool isPointOnPlane(planeinfo &p1,vector<double>&v);
bool planeIntersectPlane(planeinfo &p1,planeinfo &p2,lineinfo &l_re);
bool isLinesEqual(lineinfo &l1,lineinfo &l2);

bool segIntersectSeg(seg &s1,seg &s2,vector<double>&intersectionV,double &dT1,double &dT2);
double findLineT_seg(seg &s,vector<double>&p);
bool isPointOnLine(lineinfo &l,vector<double>&v);
bool isPointOnSeg(seg &s,vector<double>&v);


void CallTriangleForSinglePlane_3D(vector<double>&Vpos3D,vector<uint>&E2V, vector<double>&plane_para, string &cmline,  vector<double>&out_Vpos3D,vector<uint>&out_Triangles);
double calculateWindingnumber(double *p_vpos,double *plane_para,vector<double>&vseq,double *mindist);
void RatateToZaxis3Dto2D(vector<double>&Vpos3D,vector<double>&plane_para,vector<double>&Vpos2D,double &zpos);

void CrossSections::writeSegs(string filename, vector<seg>&segs){

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

    for(int i=0;i<segs.size();++i){
        auto p_v = segs[i].v1.data();
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
        p_v = segs[i].v2.data();
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;

    }


    {
        int n_edges = edge2vertices.size()/2;
        for(int i=0;i<n_edges;++i){
            auto p_ev = edge2vertices.data()+i*2;
            outer << "l " << p_ev[0]+1<< " "<< p_ev[1]+1 << endl;
        }
    }

    for(int i=0;i<segs.size();++i){
        outer << "l " << n_vertices+i*2 + 1<< " "<< n_vertices+i*2+2<< endl;

    }



    outer.close();
    cout<<"saving finish: "<<filename<<endl;

}


void CrossSections::writeSegs2(string filename, vector<seg>&segs,vector<double>&newv){

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

    for(int i=0;i<segs.size();++i){
        auto p_v = segs[i].v1.data();
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
        p_v = segs[i].v2.data();
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;

    }


    {
        int n_edges = edge2vertices.size()/2;
        for(int i=0;i<n_edges;++i){
            auto p_ev = edge2vertices.data()+i*2;
            outer << "l " << p_ev[0]+1<< " "<< p_ev[1]+1 << endl;
        }
    }

    for(int i=0;i<segs.size();++i){
        outer << "l " << n_vertices+i*2 + 1<< " "<< n_vertices+i*2+2<< endl;

    }


    int n_newv = newv.size()/3;
    for(int i=0;i<n_newv;++i){
        auto p_v = newv.data()+i*3;
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }


    outer.close();
    cout<<"saving finish: "<<filename<<endl;




}



void CrossSections::FCMGeometryProcessPipeline(string infilename, string outpath, double iline_step, double portion, string triangle_cmd, bool iswrite_intermedium_result){


    iswritemediaresult = iswrite_intermedium_result;
    ReadCrossSections(infilename);

    CalSelfProperty();

    if(iswritemediaresult)WriteCrossSectionsObjFormat(outpath + ("mediumresult_input"));

    vector<planeinfo>boundingbox;
    FindBoundingBox(boundingbox,8);

    vector<seg>segs;
    vector<vector<seg> >bsegs;
    vector<vector<int> >planelines;

    findPlaneLines(planes_info,boundingbox,segs,planelines,bsegs);


//    vector<seg>testsegs(segs);
//    for(auto &a:bsegs)for(auto &b:a)testsegs.push_back(b);
//    writeSegs(outpath + ("mediumresult"),testsegs);

    findPlaneGraphMM(segs,planelines,bsegs,iline_step, 4);

    //cout<<"Ads"<<endl;
    CallTriangle(triangle_cmd, outpath + ("mediumresult_tri"));
    //cout<<"Ads"<<endl;

    findImplicitValue(outpath + ("mediumresult_value"));

    LineSafeVertices(portion);



}

void CrossSections::FindBoundingBox(vector<planeinfo>&boundingbox,double coef){

    double scale[3];
    for(int i=0;i<3;++i)scale[i] = rightcorner[i] - leftcorner[i];
    for(int i=0;i<3;++i)scale[i]/=coef;
    boundingbox.resize(6);

    double al[3],ar[3];
    for(int i=0;i<3;++i)al[i] = leftcorner[i]-scale[i];
    for(int i=0;i<3;++i)ar[i] = rightcorner[i]+scale[i];



    boundingbox[0].build( al[0],  midpoint[1], midpoint[2], 1,0,0,-al[0] );

    boundingbox[1].build( ar[0],  midpoint[1], midpoint[2], -1,0,0,ar[0] );

    boundingbox[2].build( midpoint[0],  al[1], midpoint[2], 0,1,0,-al[1] );

    boundingbox[3].build( midpoint[0],  ar[1], midpoint[2], 0,-1,0,ar[1] );

    boundingbox[4].build( midpoint[0],  midpoint[1], al[2], 0,0,1,-al[2] );

    boundingbox[5].build( midpoint[0],  midpoint[1], ar[2], 0,0,-1,ar[2] );

}



void CrossSections::findPlaneBoundLines(planeinfo &p, vector<planeinfo > &boundingbox, vector<seg>&segments){


    segments.clear();
    for(int ibb = 0; ibb<boundingbox.size();++ibb){
        lineinfo line;

        if(!planeIntersectPlane(p,boundingbox[ibb],line))continue;
        vector<vector<double> >pnts;
        for(int jbb = 0;jbb<boundingbox.size();++jbb){

            if(jbb==ibb)continue;
            vector<double>intV;
            double dT = planeIntersectLine(boundingbox[jbb],line,intV);
            if(dT!=INFTY_DOUBLE)pnts.push_back(intV);

        }

        vector<double>v1,v2;
        //cout<<"pnts_size(): "<<pnts.size()<<endl;
        v1.clear();v2.clear();
        //int k=0;
        for(int i=0;i<pnts.size();++i){


            bool isvalid  = true;
            for(int j=0;j<boundingbox.size();++j){

                double val = dot(pnts[i].data(),boundingbox[j].planepara.data())+boundingbox[j].planepara[3];
                if(val < 0 && fabs(val)>THRES )isvalid = false;
                //cout<<val<<' ';
            }//cout<<endl;

            if(isvalid){
                if(v1.size()==0)v1 = pnts[i];
                else v2 = pnts[i];
                //++k;
            }
            if(v2.size()>0){
                segments.push_back( seg(v1,v2) );
                break;
            }

            //cout<<k<<endl;


        }



    }


}


void CrossSections::findPlaneLines(vector<planeinfo > &planes, vector<planeinfo > &boundingbox, vector<seg>&p2psegs,vector<vector<int>>&planelines,vector<vector<seg> >&p2bsegs){


    p2psegs.clear();
    planelines.clear();
    p2bsegs.clear();

    vector<lineinfo>lines;

    int nplane = planes.size();

    vector<set<int>>planelines_set(nplane);
    p2bsegs.resize(nplane);
    for(int ipp = 0; ipp<planes.size();++ipp){


        findPlaneBoundLines(planes[ipp],boundingbox,p2bsegs[ipp]);

        for(int jpp = ipp; jpp<planes.size();++jpp){

            lineinfo line;

            if(!planeIntersectPlane(planes[ipp],planes[jpp],line))continue;

            int iline = -1;

            for(int i=0;i<lines.size();++i)if(isLinesEqual(lines[i],line)){iline = i;break;}

            if(iline!=-1){
                planelines_set[ipp].insert(iline);
                planelines_set[jpp].insert(iline);

            }else{

                double dT1 = INFTY_DOUBLE,dT2 = INFTY_DOUBLE;

                vector<double>pnt1,pnt2;

                for(int ibb = 0;ibb<boundingbox.size();++ibb){

                    vector<double>pnt;
                    double dT = planeIntersectLine(boundingbox[ibb],line,pnt);
                    if(dT==INFTY_DOUBLE)continue;
                    bool isvalid = true;
                    for(int jbb=0;jbb<boundingbox.size();++jbb){
                        double val = dot(pnt.data(),boundingbox[jbb].planepara.data())+boundingbox[jbb].planepara[3];
                        if(val < 0 && fabs(val)>THRES )isvalid = false;
                    }
                    if(isvalid){
                        if(fabs(dT)<fabs(dT1)){dT2=dT1;pnt2=pnt1;dT1=dT;pnt1=pnt;}
                        else if(fabs(dT)<fabs(dT2)){dT2 = dT;pnt2=pnt;}
                    }

                }

                if(pnt1.size()==3 &&pnt2.size()==3){

                    if(dT1<dT2)p2psegs.push_back(seg(pnt1,pnt2));else p2psegs.push_back(seg(pnt2,pnt1));

                    iline = lines.size();
                    lines.push_back(line);

                    planelines_set[ipp].insert(iline);
                    planelines_set[jpp].insert(iline);

                }


            }

        }


    }


    planelines.resize(nplane);
    for(int i=0;i<nplane;++i){
        for(auto a:planelines_set[i])planelines[i].push_back(a);
    }








}



void CrossSections::findPlaneGraphMM(vector<seg>&p2psegs, vector<vector<int>>&planelines, vector<vector<seg> > &p2bsegs, double sampW, double boxSampCoef){


    double pointThres = THRES*10;

    auto listContain_ori = [](vector<double>&l,vector<double>&v,double thres){
        double dif[3];
        int nv = l.size()/3;
        auto p_d = l.data();
        auto p_v = v.data();
        for(int i=0; i<nv; ++i){
            minusVec(p_d+i*3, p_v, dif);
            if(normVec(dif)<thres)return i;
        }
        return -1;
    };
    auto listContain = [&pointThres,listContain_ori](vector<double>&l,vector<double>&v){
        return listContain_ori(l,v,pointThres);
    };




    int n_lines = p2psegs.size();
    int n_planes = planelines.size();

    V.clear();
    isectIndices.clear();
    vector<set<int>>lineIsects(n_lines);

    E.clear();
    E.resize(n_planes);

    auto PushbackV_ori = [this](double *v){
        for(int i=0;i<3;++i)V.push_back(v[i]);
        return V.size()/3-1;
    };
    auto PushbackV = [PushbackV_ori](vector<double>&v){
        return PushbackV_ori(v.data());

    };


    /*(*---- Find all intersections of lines ----*)*/


    line2plane.clear();
    line2plane.resize(n_lines);
    for(int i=0;i<n_planes;++i)for(auto a:planelines[i])line2plane[a].push_back(i);

    vector<vector<pair<int,double> > >lineTs(n_lines);

    for(int i=0;i<n_lines;++i){
        for(int j=i;j<n_lines;++j){
            vector<double>pnt;
            double dT1,dT2;

            if(!segIntersectSeg(p2psegs[i],p2psegs[j],pnt,dT1,dT2))continue;

            int ind = listContain(V,pnt);
            if(ind==-1){
                ind = PushbackV(pnt);
                isectIndices.push_back(ind);
            }

            if(lineIsects[i].find(ind)==lineIsects[i].end()){
                lineIsects[i].insert(ind);
                lineTs[i].push_back( make_pair( ind, findLineT_seg(p2psegs[i],pnt)) );
            }
            if(lineIsects[j].find(ind)==lineIsects[j].end()){
                lineIsects[j].insert(ind);
                lineTs[j].push_back( make_pair( ind, findLineT_seg(p2psegs[j],pnt)) );
            }
        }
    }


    /*(*---- Find all intersections of lines with contours ----*)*/

    vector<vector<int> >ctrIndices(n_planes);

    for(int i=0;i<n_planes;++i){
        Contour &Cs = CSs[i];

        int nedges = Cs.n_edges;
        ctrIndices[i].resize(Cs.n_vertices,-1);

        for(int j = 0; j<planelines[i].size();++j){
            int iline = planelines[i][j];
            for(int k=0; k<nedges;++k){

                auto p_ev = Cs.ev_begin(k);
                auto p_v1 = Cs.v_begin(p_ev[0]);
                auto p_v2 = Cs.v_begin(p_ev[1]);

                vector<double>pnt;
                double dT1,dT2;
                seg ctrseg(p_v1,p_v2);

                if(!segIntersectSeg(p2psegs[iline],ctrseg,pnt,dT1,dT2))continue;

                int ind = listContain(V,pnt);
                if(ind==-1){
                    ind = PushbackV(pnt);
                    //isectIndices.push_back(ind);
                }

                if(lineIsects[iline].find(ind)==lineIsects[iline].end()){
                    lineIsects[iline].insert(ind);
                    lineTs[iline].push_back(make_pair( ind, findLineT_seg(p2psegs[iline],pnt)) );
                }

                if(dT2<0.5)ctrIndices[i][p_ev[0]] = ind;
                else ctrIndices[i][p_ev[1]] = ind;

            }

        }

    }


    /*(*---- Find all points on lines ----*)*/
    lineIndices.clear();
    lineIndices.resize(n_lines);
    //vector<vector<int> >lineEdges(n_lines);
    line2Edges.clear();line2Edges.resize(n_lines);
    vector<vector<int> >lineEndPoint(n_lines);
    for(int i=0;i<n_lines;++i){
        sort(lineTs[i].begin(),lineTs[i].end(),[](pair<int,double>&l,pair<int,double>&r){return l.second<r.second;});
        double linedir[3];

        seg &s = p2psegs[i];
        minusVec(s.v2.data(),s.v1.data(),linedir);

        double normlinedir = normVec(linedir);
        int n_samp = floor(normlinedir/sampW);
        //if(normVec(linedir)-n_samp*sampW > sampW*0.25)++n_samp;

        normalize(linedir);
        int indexLast = -1;
        bool pickbe = true;
        for(int j = 0;j<n_samp+1;++j){
            double pnt[3];
            if(j==n_samp)copyVec(s.v2.data(),pnt);
            else weightedAddVec(1,j*sampW,s.v1.data(),linedir,pnt);

            vector<int>indices;
            for(int k=0;k<lineTs[i].size();++k){
                if(fabs(lineTs[i][k].second*normlinedir-j*sampW)<sampW/2.)
                    indices.push_back(lineTs[i][k].first);
            }

            if(indices.size()==0){

                indices.push_back(PushbackV_ori(pnt));


            }//else cout<<"bb"<<endl;
            if(pickbe){pickbe=false;lineEndPoint[i].push_back(indices[0]);}
            //if(indices.size()>1)cout<<"xx"<<endl;
            for(int k=0;k<indices.size();++k){
                lineIndices[i].push_back(indices[k]);
                if(indexLast>-1){line2Edges[i].push_back(indexLast);line2Edges[i].push_back(indices[k]);}
                indexLast = indices[k];
                //if((j==0 && k==1) || (j == n_samp-1 && k) )
            }



        }
        lineEndPoint[i].push_back(indexLast);

    }

    plane2lines.clear();plane2lines.resize(n_planes);
    for(int i=0;i<n_lines;++i){
        for(auto b:line2plane[i])plane2lines[b].push_back(i);
    }
    for(int i=0;i<n_lines;++i){
        for(auto b:line2plane[i])
            for(auto a:line2Edges[i])E[b].push_back(a);
    }



    /*(*---- Add verts and edges from contours ----*)*/
    for(int i=0;i<n_planes;++i){
        Contour &Cs = CSs[i];
        int n_vertices = Cs.n_vertices;
        vector<int>&ctrIndex = ctrIndices[i];
        for(int j=0;j<n_vertices;++j){
            if(ctrIndex[j]==-1){
                ctrIndex[j] = PushbackV_ori(Cs.v_begin(j));
            }
        }
        int n_edges = Cs.n_edges;
        auto &EE = E[i];
        for(int j=0;j<n_edges;++j){
            auto p_e = Cs.ev_begin(j);
            for(int k=0;k<2;++k){
                EE.push_back(ctrIndex[p_e[k]]);
            }
        }
    }


    /*(*---- Add verts and edges from bounding box ----*)*/

    vector<double>boundPnts;
    vector<int>boundIndices;
    for(int i=0;i<n_planes;++i){
        //continue;
        //if(i!=3)continue;
        vector<seg>&boundingbox = p2bsegs[i];
        for(int j=0;j<boundingbox.size();++j){
            //if(j!=1 &&j!=2 && j!=3)continue;
            seg &s = boundingbox[j];

            int startind = listContain(boundPnts,s.v1);
            if(startind==-1){
                startind = PushbackV(s.v1);
                boundIndices.push_back(startind);
                for(auto a:s.v1)boundPnts.push_back(a);
            }else{
                startind = boundIndices[startind];
            }
            int endind = listContain(boundPnts,s.v2);
            if(endind==-1){
                endind = PushbackV(s.v2);
                boundIndices.push_back(endind);
                for(auto a:s.v2)boundPnts.push_back(a);
            }else{
                endind = boundIndices[endind];
            }


            vector<pair<int,double> >segTs;
            segTs.push_back(make_pair(startind,0));
            segTs.push_back(make_pair(endind,1));



            for(int k=0;k<n_lines;++k){
                vector<double>v1,v2;
                auto p_v1 = V.data() + 3*lineEndPoint[k][0];
                auto p_v2 = V.data() + 3*lineEndPoint[k][1];
                for(int lll = 0;lll<3;++lll)v1.push_back(p_v1[lll]);
                for(int lll = 0;lll<3;++lll)v2.push_back(p_v2[lll]);
                if(isPointOnSeg(s,v1)){
                    segTs.push_back(make_pair(lineEndPoint[k][0],findLineT_seg(s,v1)));
                }
                if(isPointOnSeg(s,v2)){
                    segTs.push_back(make_pair(lineEndPoint[k][1],findLineT_seg(s,v2)));
                }
            }

            sort(segTs.begin(),segTs.end(),[](pair<int,double>&l,pair<int,double>&r){return l.second<r.second;});
            //cout<<segTs.size()<<endl;
            for(int k = 0; k < segTs.size()-1; ++k){

                //E[i].push_back(segTs[k].first);E[i].push_back(segTs[k+1].first);
                //continue;


                double linedir[3];
                double *v1 = V.data()+3*segTs[k].first;
                double *v2 = V.data()+3*segTs[k+1].first;
                minusVec(v2,v1,linedir);
                double sampWB = sampW*boxSampCoef;
                int n_samp = floor(normVec(linedir)/sampWB);
                //if( normVec(linedir)-n_samp*sampWB > sampWB*0.25 )++n_samp;

                normalize(linedir);
                int indexLast = segTs[k].first;
                int ind;
                for(int m = 1;m<n_samp+1;++m){
                    double pnt[3];
                    if(m==n_samp){
                        copyVec(v2,pnt);
                        ind = segTs[k+1].first;
                    }
                    else{
                        weightedAddVec(1,m*sampWB,v1,linedir,pnt);
                        ind = PushbackV_ori(pnt);
                    }

                    E[i].push_back(indexLast);E[i].push_back(ind);
                    indexLast = ind;

                }
                if(n_samp == 0 ){E[i].push_back(indexLast);E[i].push_back(segTs[k+1].first);}
            }
        }
    }

    //vector<uint>stackE;
    //for(auto &a:E)for(auto b:a)stackE.push_back(b);
    //for(auto &a:lineEdges)for(auto b:a)stackE.push_back(b);
    //writeSegs2(filename_test,p2psegs,V);
    //writeObjFile_line(filename_test,V,stackE);


}




void CrossSections::CallTriangle(string &cmline, string outprefix){


    int n_v = V.size()/3;
    int n_plane = E.size();
    newV.resize(n_plane);
    newE.resize(n_plane);
    T.resize(n_plane);
    planeV2globalV.resize(n_plane);
    planeVcutoffInd.resize(n_plane);
    for(int i=0;i<n_plane;++i){
        //void CallTriangleForSinglePlane_3D(vector<double>&Vpos3D,vector<uint>&E2V, vector<double>&plane_para, vector<double>&out_Vpos3D,vector<uint>&out_Triangles);
        planeV2globalV[i].clear();
        auto &planeV2globalV_i = planeV2globalV[i];
        vector<int>v_inplane(n_v,-1);
        int ind = 0;
        vector<double>Vpos3D;
        vector<uint>E2V;

        for(auto a:E[i])v_inplane[a] = 0;
        for(int j=0;j<n_v;++j)if(v_inplane[j]!=-1){
            v_inplane[j] = ind++;
            auto p_v = V.data()+j*3;
            for(int k=0;k<3;++k)Vpos3D.push_back(p_v[k]);
            planeV2globalV[i].push_back(j);
        }

        for(auto a:E[i])E2V.push_back(v_inplane[a]);
        planeVcutoffInd[i] = ind;

        //cout<<"00"<<endl;
        CallTriangleForSinglePlane_3D(Vpos3D,E2V, planes_info[i].planepara, cmline, newV[i],T[i]);
        //cout<<i<<endl;
        int newttv = newV[i].size()/3;
        int offset = V.size()/3;
        for(int j=ind;j<newttv;++j){
            auto pv = newV[i].data()+j*3;
            for(int k=0;k<3;++k)V.push_back(pv[k]);
            planeV2globalV_i.push_back(j+offset-ind);
        }

        Mesh my_mesh;
        my_mesh.ImportMesh(newV[i],T[i],true,true);
        newE[i] = my_mesh.edges;
        for(auto &a:newE[i])a = planeV2globalV_i[a];

    }



    vector<uint>AllT;
    for(int i=0;i<n_plane;++i){
        auto &newind = planeV2globalV[i];
        for(auto a:T[i])AllT.push_back(newind[a]);

    }


    //writeSegs2(outprefix,p2psegs,V);
    if(iswritemediaresult)writeObjFile(outprefix,V,AllT);

}


void CrossSections::findImplicitValue(string outprefix){


    int n_vertices = V.size()/3;
    int n_plane = planeV2globalV.size();
    val.clear();
    val.resize(n_vertices,vector<double>(n_plane*n_materials,-1));


    //void ComputeMat2VerticesLoop(vector<vector<int> > &verticesloops, vector<vector<double>>&verticesSeqs, vector<int>&loops2mats);

    for(int iplane = 0;iplane<n_plane;++iplane){
        auto p_planepara = planes_info[iplane].planepara.data();
        int npv = planeV2globalV[iplane].size();

        vector<vector<int> > verticesloops;
        vector<vector<double> >verticesSeqs;
        vector<int>loops2mats;


        CSs[iplane].ComputeMat2VerticesLoop(verticesloops, verticesSeqs, loops2mats);


        for(int i=0;i<npv;++i){

            vector<double>mindist(n_materials,INFTY_DOUBLE);
            vector<unsigned char>vwindnum(n_materials,false);
            int globalv = planeV2globalV[iplane][i];

            for(int j=0;j<loops2mats.size();++j){
                int mat = loops2mats[j];

                auto p_v = V.data()+globalv*3;
                double cmindist;
                double wnum = calculateWindingnumber(p_v,p_planepara,verticesSeqs[j],&cmindist);
                vwindnum[mat]^=(unsigned char)(round(wnum)!=0);

                mindist[mat] = min(mindist[mat],cmindist);
            }

            for(int j=0;j<n_materials;++j){
                if(mindist[j]==INFTY_DOUBLE)mindist[j] = -INFTY_DOUBLE;
                else {
                    if(j!=0)mindist[j]*=(vwindnum[j]?1:-1);
                    else mindist[j]*=(vwindnum[j]?-1:1);
                }
            }

            double maxvalue = *max_element(mindist.begin(),mindist.end());

            for(int j=0;j<n_materials;++j)if(mindist[j]==-INFTY_DOUBLE)mindist[j] = 0 - THRES * 1e10;

            auto p_vval = val[globalv].data();
            for(int j=0;j<n_materials;++j)p_vval[iplane*n_materials+j] = mindist[j];

        }

        //cout<<"adasdas"<<endl;

    }



    //if(iswritemediaresult)writePlanesVal(filename_test);




}


void CrossSections::LineSafeVertices(double portion){

    int n_vertices = V.size()/3;
    int n_plane = planeV2globalV.size();
    int n_line = lineIndices.size();

    nodeMarkers.clear();
    nodeMarkers.resize(n_vertices,0);

    vector<double>minGap(n_vertices,INFTY_DOUBLE);
    vector<bool>lineV(n_vertices,false);

    auto LineSingleSpecial = [this](double *vals,double &gap){

        double maxval = -INFTY_DOUBLE,secondmaxval = -INFTY_DOUBLE;
        int maxind = -1,secondmaxind = -1;
        for(int i=0;i<n_materials;++i){
            if(vals[i]>maxval){
                secondmaxval = maxval;secondmaxind = maxind;
                maxval = vals[i];maxind = i;
            }else if(vals[i]>secondmaxval){
                secondmaxval = vals[i];secondmaxind = i;
            }
        }

        gap = maxval-secondmaxval;
        if(gap<THRES)maxind = -1;
        return maxind;


    };




    for(int i=0;i<n_line;++i){

        for(auto ind: lineIndices[i]){
            auto p_val = val[ind].data();
            int ll = -1;
            bool bp=false;
            lineV[ind] = true;
            for(auto iplane: line2plane[i]){


                double gap;

                int lltmp = LineSingleSpecial(p_val+iplane*n_materials,gap);

                minGap[ind] = min(minGap[ind],gap);


                if(lltmp==-1){bp=true;continue;}
                if(ll==-1)ll = lltmp;
                else if(ll!=lltmp){bp=true;}

            }
            if(bp)nodeMarkers[ind] = 2;


        }

    }

    vector<pair<double,int> >consistentLineV;

    for(int i=0;i<n_vertices;++i){
        if(lineV[i]&&nodeMarkers[i]!=2)consistentLineV.push_back(make_pair(minGap[i],i));
    }

    int nconsist = consistentLineV.size();
    sort(consistentLineV.begin(),consistentLineV.end(),[](pair<double,int>&l,pair<double,int>&r){return l.first<r.first;});


    int endInd = nconsist * portion;

    for(int i=0;i<endInd;++i)nodeMarkers[consistentLineV[i].second] = 1;




}





void CrossSections::writePlanesVal(string filename){


    int pickplane = 3;

    //int pickmaterial = 2;


    vector<double>outVpos2D,rotateVpos2D,Vpos3D;
    double zpos;

    RatateToZaxis3Dto2D(newV[pickplane],planes_info[pickplane].planepara,rotateVpos2D,zpos);




    int n_v = planeV2globalV[pickplane].size();
    auto p_vind = planeV2globalV[pickplane].data();
    Vpos3D.clear();
    for(int i=0;i<n_v;++i){
        auto p_v = rotateVpos2D.data()+i*2;
        Vpos3D.push_back(p_v[0]);Vpos3D.push_back(p_v[1]);Vpos3D.push_back(0);
    }
    writeObjFile(filename+string("_p"),Vpos3D,T[pickplane]);
    for(int pickmaterial=0;pickmaterial<n_materials;++pickmaterial){
        for(int i=0;i<n_v;++i){
            double z = val[p_vind[i]][pickplane*n_materials+pickmaterial];
            Vpos3D[i*3+2] = z;
        }

        writeObjFile(filename+string("_")+to_string(pickmaterial),Vpos3D,T[pickplane]);
    }
}

