#ifndef _CONFIG_H_
#define _CONFIG_H_

// this header fixes a number of platform specific problems
// make sure you include it at the top of every file you write

#ifdef _MSC_VER
// turn off symbol length warnings
#pragma warning (disable: 4786)
#pragma warning (disable: 4503)
//fix non-compliant 'for' scoping
#define for if(false) {} else for 
// conversion double -> double
#pragma warning (disable: 4244)
// truncation double -> double
#pragma warning (disable: 4305)
#endif

// define usefull constants

#ifndef M_PI
#define M_PI               3.1415926535897932384626433832795
#endif

#include <cassert>
#include <iostream>
/*
#include <FL/gl.h>
#include <fl/math.h>
*/
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

#include "mymath.h"
#include <vector>
#include <set>
#include <map>
#include <stack>

typedef vector<int> intvector;
typedef vector<double> doublevector;
typedef set<int> intset;
const double TOLERANCE_HALF = 0.5;
const double TOLERANCE_POS_ONE = 1;
const double TOLERANCE_SIX = 0.000001;
const double TOLERANCE_THREE = 0.001;
const double TOLERANCE_FOUR = 0.0001;
const double TOLERANCE_TWO = 0.01;
const double TOLERANCE_ONE = 0.1;
const double TOLERANCE_ANGLE = 0.01;
const double POINT_ON_PLANE_TOLERANCE = 0.000001;
const double PLANE_PARAM_PARALLEL_TOLERANCE = 0.000001;

const int CTRTYPE_SUBVER = 1;
const int CTRTYPE_SUBEDGE = 2;
const int CTRTYPE_SUBFACE = 3;
const int CTRTYPE_MAJPT = 1;
const int CTRTYPE_MASEAM = 2;
const int CTRTYEP_MASHEET = 3;
const int NOT_EXIST = -1;

//expected length for contour edge when preprocessing
const double EXPECTLEN = 20;
const double TOLERANCE_SAME_VER = 0.01;

#ifdef LINUX
using std::min;
using std::max;
#endif

#define debug
const double PROCSIZE = 1000;

extern void RedrawWindow();

//global variables
extern double DIM;
extern bool _global_IsCtrGraph;
extern vector<set<int> > _gloabl_radialPlaneSet;

#endif /* _CONFIG_H_ */
