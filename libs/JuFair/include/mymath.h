#ifndef _MYMATH_H_
#define _MYMATH_H_

#include "config.h"

class MyMath{
public:
	static int getSign( double val );
	static bool isEqualInToler( double num1, double num2, double toler);
	//represented in x, y, z
	static double vectorLen(double x, double y , double z);
	//represented in array
	static double dotProduct( double vec1[3], double vec2[3]);
	//compute the cos of angle pt1 pt0 pt2. angle between vector pt0pt1 and pt0pt2
	static double getCosOfAngle(double pt1[3], double pt2[3], double pt0[3]);
	//get vector between two poitns, and save the vector into vec
	static void getVec(double pt1[3], double pt2[3], double vec[3]);
	static void getNeg(double vec[3]);
	static void getPtOnSeg(double pt1[3], double pt2[2], double ratio, double pt[3]);
	static void getEndPtOfSeg(double ept[ 3 ], double vec[ 3 ], double resultpt[ 3 ]);
	static void getPtOnRay(double pt[ 3 ], double dir[ 3 ], double t, double endpt[ 3 ]);
	//compute the length of the vector between the two points
	static double vectorlen( double pt1[3], double pt2[3]);
	//compute the vector's length
	static double vectorlen(double vec[3]);
	//normalize the svec, and keep it unchanged, save the normalized value into dvec
	static void normalize(double svec[3], double dvec[3]);
	//normalize vec, and save the normalized value into vec
	static void normalize(double vec[3] );
	//compute center of two points p1, and p2, save it into center
	static void center(double p1[3], double p2[3], double center[3]);
	//compute the relative position to center, and unit lenght is len.
	static void getrelativepos(double pt[3], double center[3], double len);
	/**
	* return normalized crossproduct of vec1 and vec2  and save it into cp
	*/	
	static void crossProduct(double vec1[3], double vec2[3], double cp[3]);
	/**
	* return non-normalized cross product of vec1, vec2 and save it in cp
	*/
	static void crossProductNotNorm(double vec1[3], double vec2[3], double cp[3]);
	//a naive algorithm to compute the difference for vertex sver to all the vertices in dver.
	static double computeDiff( double* sver,  double* dver,int dvernum);
	
	static void stretchVec( double vec[ 3 ], double ratio[ 3 ] );
	static void stretchVec( double vec[ 3 ], double ratio);
	//another naive algorithm to compute the difference for vertex sver to another mesh, dvers, dtrians

	//dist: save the resulting distance from the point to the line segment lies on
	//dist = the distance from the point to the line the segment lies on , if the perpendicular line from sver crosses the seg
	//dist = the minimal distance from the point to one of the two endpoints of the segment, otherwise
	static void computeDistPt2Seg( double* sver, double* dver, int ind[ 2 ], double& dist);
	//compute distance from point to a triangle 
	//dist = distance from the point to the plane the triangle is on, if the perpendicular line crosses trian
	//dist = nearest distance from the point to the three segments
	static void computeDistPt2Trian(double* sver, double* dver, int ind[ 3 ], double& dist);
	
	static double computeDiff( double* sver, double* dver, int* dtrians, int triannum);
};

#endif