TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11
CONFIG += m64


HEADERS += \
    constrainOpt/aclass.h \
    mymesh/InfoStruct.h \
    mymesh/my_mesh.h \
    mymesh/readers.h \
    mymesh/utility.h \
    mymesh/contour.h \
    mymesh/UnionFind.h \
    mymesh/geo_curv.h \
    constrainOpt/DynamicMesh.h

SOURCES += \
    constrainOpt/DynamicMesh.cpp\
    constrainOpt/transformer.cpp \
    mymesh/my_mesh.cpp \
    mymesh/readers.cpp \
    main.cpp \
    constrainOpt/a_aclass_contour.cpp \
    mymesh/contour.cpp \
    mymesh/geo_curv.cpp \
    mymesh/contour_fcm.cpp \
    mymesh/contour_geo.cpp \
    constrainOpt/a_aclass_crosssection.cpp \
    constrainOpt/a_class_iteropt.cpp \
    mymesh/UnionFind.cpp



INCLUDEPATH += $$/opt/local/include
INCLUDEPATH += $$/opt/local/include/eigen3


unix|win32: LIBS += -L$$/Library/gurobi702/mac64/lib/ -lgurobi_c++ -lgurobi70

INCLUDEPATH += $$/Library/gurobi702/mac64/include
DEPENDPATH += $$/Library/gurobi702/mac64/include

macx: LIBS += -L$$PWD/mymesh/ -ltriangle

LIBS += -L/usr/local/Cellar/suite-sparse/5.2.0/lib/ -lcholmod -lcolamd
