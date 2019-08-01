#include "mymath.h"

int MyMath::getSign( double val )
{
	if( val > 0 )
		return 1;
	if( val == 0) 
		return 0;
	return -1;
}
void MyMath::getPtOnRay(double pt[ 3 ], double dir[ 3 ], double t, double endpt[ 3 ])
{
	for( int i = 0 ;i < 3; i++)
	{
		endpt[ i ] = pt[ i ] + dir[ i ] * t;
	}

}
bool MyMath::isEqualInToler( double num1, double num2, double toler)
{
	if( abs( num1 - num2 ) < toler )
		return true;
	return false;
}

double MyMath::vectorLen(double x, double y , double z)
{
	return sqrt(x*x + y*y + z* z);
}
//represented in array
double MyMath::dotProduct( double vec1[3], double vec2[3])
{
	double dotp = 0;
	for( int i = 0; i < 3; i ++)
		dotp += (vec1[ i ] * vec2[ i ]);
	return dotp;
}
//compute the cos of angle pt1 pt0 pt2. angle between vector pt0pt1 and pt0pt2
double MyMath::getCosOfAngle(double pt1[3], double pt2[3], double pt0[3])
{	
	double vec1[ 3 ], vec2[3] ;
	getVec( pt0, pt1, vec1);
	getVec( pt0, pt2, vec2);
	//////////////////////////////////////////////////////////////////////////
	//	cout<<"length:"<<vectorlen(vec1)<<"  "<<vectorlen(vec2)<<endl;
	return (dotProduct( vec1, vec2)/(vectorlen(vec1)*vectorlen(vec2)));
}
//get vector between two poitns, and save the vector into vec
void MyMath::getVec(double pt1[3], double pt2[3], double vec[3])
{
	for( int i = 0; i < 3; i ++)
		vec[ i ] = pt2[ i ] - pt1[ i ];
}
void MyMath::getNeg(double vec[3])
{
	for( int i = 0;i < 3; i++)
		vec[ i ] = -vec[ i ];
}
void MyMath::getPtOnSeg(double pt1[3], double pt2[2], double ratio, double pt[3])
{
	for(int i = 0; i < 3; i ++)
	{
		pt[ i ] = (1-ratio)*pt1[ i ] + ratio*pt2[ i ];
	}
}
void MyMath::getEndPtOfSeg(double ept[ 3 ], double vec[ 3 ], double resultpt[ 3 ])
{
	for( int i = 0; i< 3; i ++)
		resultpt[ i ] = ept[ i ] + vec[ i ];
}
//compute the length of the vector between the two points
double MyMath::vectorlen( double pt1[3], double pt2[3])
{
	double len = 0; 
	for( int i = 0; i < 3; i ++)
		len += pow((pt2[i] - pt1[i]),2 );
	return sqrt( len );
}
//compute the vector's length
 double MyMath::vectorlen(double vec[3])
{
	double len = 0;
	for( int i = 0; i < 3; i ++)
		len = len + vec[i] * vec[i];
	len = sqrt( len );
	return len;
}
//normalize the svec, and keep it unchanged, save the normalized value into dvec
void MyMath::normalize(double svec[3], double dvec[3])
{
	double len = vectorlen( svec );
	for( int i = 0; i < 3; i ++)
		dvec[ i ] = svec[ i ]/len;
}
//normalize vec, and save the normalized value into vec
void MyMath::normalize(double vec[3] )
{
	double len = vectorlen( vec) ;
	for( int i = 0; i < 3; i ++)
		vec[ i ] = vec[ i ]/len;
}
//compute center of two points p1, and p2, save it into center
void MyMath::center(double p1[3], double p2[3], double center[3])
{
	for( int i = 0; i < 3; i ++)
		center[ i ] = (p1[ i ] + p2[ i ])/2;
}
//compute the relative position to center, and unit lenght is len.
void MyMath::getrelativepos(double pt[3], double center[3], double len)
{
	for( int i = 0; i < 3; i ++)
	{
		pt[ i ] = (pt[ i ] - center[ i ])/len;
	}
}

/**
* return normalized crossproduct of vec1 and vec2  and save it into cp
*/	
void MyMath::crossProduct(double vec1[3], double vec2[3], double cp[3])
{
	cp[ 0 ] = vec1[ 1 ] * vec2[ 2 ] - vec1[ 2 ] * vec2[ 1 ];
	cp[ 1 ] = vec1[ 2 ] * vec2[ 0 ] - vec1[ 0 ] * vec2[ 2 ];
	cp[ 2 ] = vec1[ 0 ] * vec2[ 1 ] - vec1[ 1 ] * vec2[ 0 ];
	normalize( cp );
}
/**
* return non-normalized cross product of vec1, vec2 and save it in cp
*/
void MyMath::crossProductNotNorm(double vec1[3], double vec2[3], double cp[3])
{
	cp[ 0 ] = vec1[ 1 ] * vec2[ 2 ] - vec1[ 2 ] * vec2[ 1 ];
	cp[ 1 ] = vec1[ 2 ] * vec2[ 0 ] - vec1[ 0 ] * vec2[ 2 ];
	cp[ 2 ] = vec1[ 0 ] * vec2[ 1 ] - vec1[ 1 ] * vec2[ 0 ];
}

//a naive algorithm to compute the difference for vertex sver to all the vertices in dver.
double MyMath::computeDiff( double* sver,  double* dver,int dvernum)
{
	double minval = vectorlen( sver, dver );
	double val;
	for( int i = 1; i < dvernum; i ++)
	{
		val = vectorlen( sver, dver + 3*i);
		if( val < minval )
			minval = val;
	}
	return minval;
}

void MyMath::stretchVec( double vec[ 3 ], double ratio[ 3 ] )
{
	for( int i = 0; i < 3; i++)
		vec[ i ] *= ratio[ 3 ];
}

void MyMath::stretchVec( double vec[ 3 ], double ratio)
{
	for( int i = 0; i < 3; i++)
		vec[ i ] *= ratio;
}

//another naive algorithm to compute the difference for vertex sver to another mesh, dvers, dtrians

//dist: save the resulting distance from the point to the line segment lies on
//dist = the distance from the point to the line the segment lies on , if the perpendicular line from sver crosses the seg
//dist = the minimal distance from the point to one of the two endpoints of the segment, otherwise
void MyMath::computeDistPt2Seg( double* sver, double* dver, int ind[ 2 ], double& dist)
{
	double vec[ 3 ];
	getVec( &dver[ 3 * ind[ 0 ]], &dver[ 3*ind[ 1 ]], vec );
	double veclen = vectorlen( vec );

	//degenerate case
	if( veclen < 0.001 )
	{
		getVec( sver, &dver[ 3*ind[ 0 ]], vec);
		dist = vectorlen( vec );
		return;
	}

	stretchVec( vec, 1/veclen);
	double vec2[ 3 ];
	getVec(  &dver[ 3 * ind[ 0 ]], sver, vec2);
	double projlen = dotProduct( vec, vec2 );
	bool result = true;
	if( projlen < 0 || projlen > veclen )
		result = false;

	if( result )
	{
		double vec2len = vectorlen( vec2 );
		dist = sqrt( vec2len * vec2len - projlen*projlen);
		return;
	}

	dist = vectorlen( sver, &dver[ind[0]*3]);
	double tdist = vectorlen( sver, &dver[ind[ 1 ]*3]);
	if( tdist < dist )
		dist = tdist;
}

//compute distance from point to a triangle 
//dist = distance from the point to the plane the triangle is on, if the perpendicular line crosses trian
//dist = nearest distance from the point to the three segments
void MyMath::computeDistPt2Trian(double* sver, double* dver, int ind[ 3 ], double& dist)
{
	bool result = true;

	//decide if perpendicular line crosses triangle or not
	double tvec1[ 3 ], tvec2[ 3 ];
	double vecnorm[ 3 ][ 3 ];
	for( int i = 0; i < 3; i ++)
	{
		getVec( sver, &dver[ 3*ind[ i ]], tvec1);
		getVec( sver, &dver[ 3*ind[ (i+1)%3 ]], tvec2);
		crossProductNotNorm( tvec1, tvec2, vecnorm[ i ]);
	}
	double dotp;
	for( int i = 0; i < 3; i ++)
	{
		dotp = dotProduct( vecnorm[ i ], vecnorm[( i + 1)%3 ]);
		if( dotp < 0 )
		{
			result = false;
			break;
		}
	}

	//compute the distance from the point to the plane which triangle lies in.
	double norm[ 3 ];
	getVec( &dver[ 3*ind[ 0 ]], &dver[ 3*ind[ 1 ]], tvec1);
	getVec( &dver[ 3* ind[ 0 ]], &dver[ 3*ind[ 2 ]], tvec2);
	crossProductNotNorm( tvec1, tvec2, norm);
	double len = vectorlen( norm );
	//degenerate case - two segments lie on the same line
	if( len < 0.000001 )
	{
		result = false;
	}
	stretchVec( norm, 1/len);

	if( result )	//compute the distance from the point to the intersection point and return
	{
		getVec( sver, &dver[ 3*ind[ 0 ]], tvec2);
		dist = dotProduct( tvec2, norm);
		if( dist < 0 ) dist = -dist;
		return;
	}

	//compute the distance from the point to the three segments and select the minimal
	double tdist;
	int tind[ 2 ];
	for( int i = 0; i < 3; i ++)
	{
		tind[ 0 ] = ind[ i ];
		tind[ 1 ] = ind[ (i+1)%3 ];
		computeDistPt2Seg( sver, dver, tind, tdist);
		if( i == 0 )
			dist = tdist;
		else
		{
			if( tdist < dist )
				dist = tdist;
		}
	}
}


double MyMath::computeDiff( double* sver, double* dver, int* dtrians, int triannum)
{
	//go through all the triangles
	bool succ = true;
	double dist, mdist;
	for( int i = 0; i < triannum; i ++)
	{
		computeDistPt2Trian( sver, dver, &dtrians[ 3*i ], dist);
		if( i == 0 )
		{
			mdist = dist;
		}
		else
		{
			if( dist < mdist )
				mdist = dist;
		}
	}
	return mdist;
}