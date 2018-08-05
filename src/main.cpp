#include <iostream>
#include<time.h>
#include<math.h>
#include"constrainOpt/aclass.h"
#include"mymesh/contour.h"
#include<unistd.h>
using namespace std;

int main(int argc,char** argv)
{
    //cout << "Hello World!" << endl;


    string infilename("data/mousebrain/mousebrain.contour");
    string outpath("data/mousebrain/");
    string cmline("pzQYYa0.001");
    double iline_step = 0.02;
    double portion = 0.1;
    double lambda = 1e-3;
    int maxiter = 30;
    bool iswrite_intermedium_result = false;

    int c;
    optind=1;
    while ((c = getopt(argc, argv, "i:o:l:s:p:m:t:w")) != -1) {
        switch (c) {
        case 'i':
            infilename = optarg;
            break;
        case 'o':
            outpath = string(optarg);
            break;
        case 's':
            iline_step = atof(optarg);
            break;
        case 'p':
            portion = atof(optarg);
            break;
        case 'm':
            maxiter = atoi(optarg);
            break;
        case 't':
            cmline = string(optarg);
            break;
        case 'l':
            lambda = atof(optarg);
            break;
        case 'w':
            iswrite_intermedium_result = true;
            break;
        case '?':
            cout << "Bad argument setting!" << endl;
            break;
        }
    }

    cout<<"input file: "<<infilename<<endl;
    cout<<"output path: "<<outpath<<endl;
    cout<<"triangle cmd: "<<cmline<<endl;

    cout<<"line sample step: "<<iline_step<<endl;
    cout<<"active portion: "<<portion<<endl;
    cout<<"lambda: "<<lambda<<endl;
    cout<<"iswrite intermedium result"<<iswrite_intermedium_result<<endl;


    Graphs a;

    CrossSections css;
    css.FCMGeometryProcessPipeline(infilename, outpath, iline_step, portion, cmline, iswrite_intermedium_result);
    a.CrossectionsPipeline(css, outpath, lambda, maxiter);


    return 0;
}


