#include "meshFair.h"
#include "HashMap.h"
#include "mymath.h"

LuMesh::LuMesh()
{
}
LuMesh::~LuMesh()
{
}
void LuMesh::InputData( const std::vector<double> &mver,const intvector &mface, 
	const intvector &ctrmedge, const bool &doswap)
{
	interpolate = true;
	doFirstSwap = doswap;
	sufvernum = mver.size()/3;
	suffacenum = mface.size()/3;

	sufver = new double[ sufvernum * 3];
	sufface = new int[ suffacenum * 3 ];
	suffacenorm = new double[ suffacenum * 3];
	sufmat = new int[ suffacenum * 2 ];

	///////// check if all faces' direction is consistent;
	std::vector< std::vector< std::pair<int,std::vector< std::pair<int,int> > > > >verInfo(sufvernum);
	//   vertex    adjvers           adjverID       adjfaces    faceID& dierection    

	for(int i=0;i<suffacenum;i++){
		std::vector< std::pair<int,int> > edges;
		edges.push_back(std::pair<int,int>(mface[i*3],mface[i*3+1]));
		edges.push_back(std::pair<int,int>(mface[i*3+1],mface[i*3+2]));
		edges.push_back(std::pair<int,int>(mface[i*3+2],mface[i*3]));

		for(int c=0;c<3;c++){
			int v1,v2,dir; 
			if(edges[c].first<=edges[c].second){ v1=edges[c].first;v2=edges[c].second; dir=1;}
			else{v1=edges[c].second;v2=edges[c].first; dir=-1;}

			bool newEdge=true;
			for(int j=0;j<verInfo[v1].size();j++){
				if(verInfo[v1][j].first==v2){ //find edge;
					verInfo[v1][j].second.push_back(std::pair<int,int>(i,dir));
					newEdge=false;
				}
			}
			if(newEdge==true){
				std::vector< std::pair<int,int> > temp;
				temp.push_back(std::pair<int,int>(i,dir));
				verInfo[v1].push_back(std::pair<int,std::vector< std::pair<int,int> > >(v2,temp));
			}
		}
	}

	std::vector<bool> usedFace(suffacenum,false);
	std::vector<bool> flipFaceNormal(suffacenum,false);
	std::vector<int> floodPool;
	floodPool.push_back(0); usedFace[0]=true;
	while (floodPool.empty()==false){
		int i=floodPool.back(); floodPool.pop_back();
		//floodPool.erase(floodPool.begin());

		std::vector< std::pair<int,int> > edges;
		edges.push_back(std::pair<int,int>(mface[i*3],mface[i*3+1]));
		edges.push_back(std::pair<int,int>(mface[i*3+1],mface[i*3+2]));
		edges.push_back(std::pair<int,int>(mface[i*3+2],mface[i*3]));

		for(int c=0;c<3;c++){
			int v1,v2; 
			if(edges[c].first<=edges[c].second){ v1=edges[c].first;v2=edges[c].second;}
			else{v1=edges[c].second;v2=edges[c].first;}
			
			for(int j=0;j<verInfo[v1].size();j++){
				if(verInfo[v1][j].first==v2){ //find edge;
					int dir;
					for(int k=0;k<verInfo[v1][j].second.size();k++){
						if(verInfo[v1][j].second[k].first==i) dir=verInfo[v1][j].second[k].second; 
					}

					for(int k=0;k<verInfo[v1][j].second.size();k++){
						int faceID=verInfo[v1][j].second[k].first;
						if(usedFace[faceID]==false){
							usedFace[faceID]=true;
							flipFaceNormal[faceID] = dir==verInfo[v1][j].second[k].second?
								!flipFaceNormal[i]:flipFaceNormal[i];
							floodPool.push_back(faceID);
						}
					}
				}
			}
		}
	}

	//////////////////
	
	int matMark[18];
	for( int i = 0; i < 18; i ++)
		matMark[ i ] = 0;

	int ind = 0;
	for( int i = 0; i < sufvernum; i++)
	{
		for( int j= 0; j < 3; j ++ )
		{
            sufver[ ind ] = mver[ ind ];
			ind ++;
		}		
	}
	
	//read face in
	double vec1[3];
	double vec2[3];
	ind = 0;
	int ind2 = 0;
	int ind3 = 0;
	for( int i = 0; i< suffacenum; i++)
	{
		for( int j= 0; j < 3; j++ )
		{
			sufface[ ind2 ] = mface[ ind ++ ];
			ind2 ++;
		}
		for( int j = 0; j < 2; j ++)
		{
			if(!flipFaceNormal[i])
				sufmat[ ind3 ] = j+1/*mface[ ind ++ ]*/;
			else
				sufmat[ ind3 ] = 2-j/*mface[ ind ++ ]*/;

			ind3 ++;
		}

		matMark[sufmat[i * 2]] = 1;
		matMark[sufmat[i * 2 + 1]] = 1;
		MyMath::getVec(&sufver[ sufface[ i * 3] * 3 ], &sufver[ sufface[ i * 3 + 1] * 3 ], vec1 );
		MyMath::getVec(&sufver[ sufface[ i * 3] * 3 ], &sufver[ sufface[ i * 3 + 2] * 3 ], vec2 );
		MyMath::crossProduct( vec1, vec2, &suffacenorm[ i * 3 ]);
	}

/*
	for( int i = 0; i < 18;i ++ )
	{
		if( matMark[ i ] == 1)
		{
			matlist.push_back( i );
		}
	}
*/

	//read contour edge in
	sufctredgenum = ctrmedge.size()/2;
	if( sufctredgenum != 0 )
		sufctredge = new int[ sufctredgenum * 2 ];
	ind = 0;
	for( int i = 0; i < sufctredgenum * 2 ; i ++)
	{
		sufctredge[ i ] = ctrmedge[ i ];
	}
}

void LuMesh::OutputData( std::vector<double> &mver, intvector &mface)
{
	mver.clear();
	mface.clear();
	int ind = 0;
	for( int i = 0; i < sufvernum; i++)
	{
		for( int j= 0; j < 3; j ++ )
		{
			mver.push_back(sufver[ ind ]);
			ind ++;
		}		
	}

	ind = 0;
	int ind2 = 0;
	for( int i = 0; i< suffacenum; i++)
	{
		for( int j= 0; j < 3; j++ )
		{
			mface.push_back(sufface[ ind2 ]);
			ind2 ++;
		}
	}
}

void LuMesh::JUFair(double ratio, int times)
{
	cout<<"In JuFair! GOOD!"<<"ratio:"<<ratio<<"times:"<<times<<endl;
	cout<<"smoothing....\t";
	//initialization
	double* oldver;
	double* newver;
	double* oldfirstdiffer;
	double* newfirstdiffer;
	double* seconddiffer;
	oldver = new double[ sufvernum*3 ];
	newver = new double[ sufvernum*3 ];
	oldfirstdiffer = new double[sufvernum*3];
	newfirstdiffer = new double[sufvernum*3];
	seconddiffer = new double[sufvernum*3];

	/* Added by tao */
	double* mags = new double[ sufvernum ] ;
	double* norms = new double[ sufvernum * 3 ] ;
	/* end adding */


	for( int i = 0; i < sufvernum*3; i++)
	{	
		newver[ i ] = sufver[ i ];
	}

	//mark all the types of the vertices	0 - normal vertex 1 - contour vertex 2 - manifold vertex 3-both
	//find all the neighbors for each vertex
	vector<int> vermark;
	vector<intvector> verneighbr;
	vector<int> edgelist;
	HashMap ver2edgehash;
	gatherInfoForFair(vermark, verneighbr, ver2edgehash, edgelist);
	vector<double>wlist;
	wlist.resize( edgelist.size()/2);
	for( int i = 0; i < times; i ++)
	{	
		double* temp;
		temp = oldver;
		oldver = newver;
		newver = temp;
		temp = NULL;
		JUFairCenter(ratio, oldver, newver, oldfirstdiffer, newfirstdiffer,seconddiffer, 
			vermark, verneighbr, ver2edgehash, edgelist,wlist,norms,mags);
	}	
	cout<<"done!"<<endl;

/*
	debugpts.clear();
	if( selectVerList.size() > 0 ){
		cout<<"current ver:"<<selectVerList[ 0 ]<<"neighbors number:"<<verneighbr[ selectVerList[ 0 ]].size()<<endl;
		for( int i = 0;  i < verneighbr[selectVerList[0]].size(); i++)
		{
			debugpts.push_back( verneighbr[ selectVerList[ 0 ]][ i ]);
			cout<<verneighbr[ selectVerList[ 0 ]][ i ]<<"  " ;
		}
		cout<<endl;
	}
*/

	//deallocate
	for( int i = 0; i < sufvernum*3; i ++)
		sufver[ i ] = newver[ i ];
	delete []oldver;
	delete []newver;
	delete []oldfirstdiffer;
	delete []newfirstdiffer;
	delete []seconddiffer;
	delete []mags;
	delete []norms;
	vermark.clear();
	for( int i = 0; i < verneighbr.size(); i++)
		verneighbr[ i ].clear();
	verneighbr.clear();
	edgelist.clear();
	wlist.clear();

	//reset face normal
	double vec1[3];
	double vec2[3];
	for( int i = 0; i< suffacenum; i++)
	{
		MyMath::getVec(&sufver[ sufface[ i * 3] * 3 ], &sufver[ sufface[ i * 3 + 1] * 3 ], vec1 );
		MyMath::getVec(&sufver[ sufface[ i * 3] * 3 ], &sufver[ sufface[ i * 3 + 2] * 3 ], vec2 );
		MyMath::crossProduct( vec1, vec2, &suffacenorm[ i * 3 ]);
	}
}



void LuMesh::gatherInfoForFair(vector<int>& vermark, vector<intvector>& verneighbr, HashMap& ver2edgehash, 
							 vector<int>&edgelist)
{
	//vector<int> edgelist;
	vector<intvector> edge2facelist;
	//edgelist and edge2facelist and ver2edgehash
	int edgelen = 0;
	int edgeindex;
	for( int i = 0;  i < suffacenum; i++)
	{
		for( int j = 0; j < 3; j ++)
		{
			edgeindex = ver2edgehash.findInsertSort(sufface[ i*3 + j ], sufface[ i*3 + (j+1)%3 ], edgelen);
			if(edgeindex == edgelen)
			{
				edgelen ++;
				edgelist.push_back(sufface[ i*3 + j ]);
				edgelist.push_back(sufface[ i*3 + (j+1)%3 ]);
				intvector facelist;
				facelist.push_back( i );	//i face
				edge2facelist.push_back( facelist );	//in c++, copy and push back
				facelist.clear();		//clear temp variable, avoid memory leakage
			}
			else
			{
				edge2facelist[ edgeindex ].push_back( i );	//face i
			}
		}
	}

	//edgemark
	int* edgemark = new int[ edgelen ];		// 0 - normal edge 1 - contour edge 2 - nonmanifold edge
	for( int i= 0; i < edgelen; i++)
		edgemark[ i ] = 0;
	for( int i = 0; i < sufctredgenum; i++)
	{
		edgeindex = ver2edgehash.findInsertSort( sufctredge[ i * 2], sufctredge[i * 2 +1], edgelen);
		if( edgeindex == edgelen )
		{
			cout<<sufctredge[ i * 2]<<"\t"<<sufctredge[i * 2 +1]<<endl;
			cout<<"not valid input!"<<endl;
			continue;
		}
		edgemark[edgeindex] = 1;
	}

	for( int i = 0; i < edgelen; i ++)
	{
		if( edge2facelist[ i ].size() > 2)
		{
			edgemark[ i ] = 2;
		}
	}
	vermark.resize( sufvernum );
	for( int i = 0; i < sufvernum; i ++)	//0-normal vertex 1 - contour vertex 2 - non-manifold vertex 3 - contour and non-manifold
	{
		vermark[ i ] = 0;
	}
	for( int i = 0; i < edgelen; i ++)
	{
		if(edgemark[ i ] == 1)	//contour edge
		{
			if( vermark[edgelist[ i * 2 ]] == 2 )
				vermark[edgelist[ i * 2 ]] = 3;
			else
				vermark[edgelist[ i * 2 ]] = 1;
			if(vermark[ edgelist[i*2 + 1]] == 2)
				vermark[edgelist[ i * 2 + 1 ]] = 3;
			else
				vermark[ edgelist[i*2 + 1]] = 1;
		}
		else if(edgemark[ i ] == 2)	//non-manifold edge
		{
			if( vermark[edgelist[ i * 2 ]] == 1 )
				vermark[edgelist[ i * 2 ]] = 3;
			else
				vermark[edgelist[ i * 2 ]] = 2;
			if(vermark[ edgelist[i*2 + 1]] == 1)
				vermark[edgelist[ i * 2 + 1 ]] = 3;
			else
				vermark[ edgelist[i*2 + 1]] = 2;
		}
	}

	//verneighbr
	verneighbr.resize( sufvernum );
	int veri[2];
	for( int i = 0; i < edgelen; i ++)
	{
		veri[ 0 ] = edgelist[ 2*i ];
		veri[ 1 ] = edgelist[ 2*i + 1];
		for( int j = 0; j < 2; j ++)
		{
			if( vermark[veri[ j ]] == 2 || vermark[ veri[j]] == 3) //non-manifold vertex
			{
				if ( edgemark[ i ] == 2 ) //non-manifold edge
					verneighbr[ veri[ j ]].push_back( veri[ 1 - j ]);
			}
			else	//either contour vertex or normal vertex
			{
				verneighbr[ veri[ j ]].push_back( veri[ 1 - j ]);
			}
			//if( vermark[veri[ j ]] == 2 ) //non-manifold vertex
			//{
			//	if ( vermark[veri[ 1 - j ]] == 2 ) //non-manifold vertex too!
			//		verneighbr[ veri[ j ]].push_back( veri[ 1 - j ]);
			//}
			//else	if( vermark[veri[ j ]] == 0 ) //either contour vertex or normal vertex
			//{
			//	if ((vermark[veri[ 1 - j ]] == 2) || (vermark[ veri[1-j]] == 0)) //non-manifold vertex or normal vertex
			//		verneighbr[ veri[ j ]].push_back( veri[ 1 - j ]);
			//}
		}		
	}
	delete []edgemark;
	//edgelist.clear();
	for( int i = 0; i < edge2facelist.size(); i ++)
		edge2facelist[ i ].clear();
	edge2facelist.clear();
}


void inline LuMesh::JUFairCenter(double ratio, double* oldver, double* newver,  double* oldfirstdiffer, double* newfirstdiffer, 
							   double* seconddiffer, vector<int>&vermark, vector<intvector>& verneighbr,
							   HashMap& ver2edgehash, vector<int>&edgelist, vector<double>&wlist,
							   double* norms, double* mags)
{
	double* dis = new double[ 3 * sufvernum ] ;
	double* vecs = new double[ 3 * sufvernum ] ;

	//	vector<double> wlist;
	computeWeightList(edgelist, wlist);
	//compute oldfirstdiffer
	computeDiffer(oldver, oldfirstdiffer, wlist, verneighbr, ver2edgehash);

	//compute local coordinates
	computeNorm( oldver, norms, mags );

	for ( int i = 0 ; i < sufvernum ; i ++ )
	{
		if ( vermark[i] < 2 )
		{
			// manifold points
			double s = mags[ i ] ;
			double m = s * s ;
			dis[3 * i] = MyMath::dotProduct( &(oldfirstdiffer[ 3 * i ]), &(norms[ 3 * i ]) ) / ( m * s ) ;
			for ( int k = 0 ; k < 3 ; k ++ )
			{
				norms[ 3 * i + k ] *= ( s / m ) ;
				vecs[ 3 * i + k ] = oldfirstdiffer[ 3 * i + k ] - dis[3 * i] * norms[ 3 * i + k ] ;
			}
		}
		else
		{
			// non-manifold
			//dis[3*i] = 0 ;
			dis[ 3*i ] = oldfirstdiffer[ 3* i];
			dis[ 3*i + 1 ] = oldfirstdiffer[ 3*i + 1];
			dis[ 3*i + 2 ] = oldfirstdiffer[ 3*i + 2];
		}
	}

	//averaging orthogonal distances
	computeDiffer(dis, seconddiffer, wlist, verneighbr, ver2edgehash);

	vector<double> w2list;
	w2list.resize(sufvernum);
	computeRefreshWeightList(verneighbr, w2list);

	// finally, compute orthogonal and tangential movements
//	int tmove=0;
	for ( int i = 0 ; i < sufvernum ; i ++ )
	{
		if( interpolate )
		{
			if( vermark[ i ] == 1 || vermark[ i ] == 3)
			{
				// Contour points: don't move
				for( int j = 0; j< 3; j++)
					newver[ 3*i + j] = oldver[ 3*i + j];
				continue;
			}
		}
		if ( vermark[ i ] == 2 || vermark[ i ] == 3 )
		{
			// Non-manifold points: Laplacian
			for ( int k = 0 ; k < 3 ; k ++ )
			{
				newver[ 3*i + k] = oldver[ 3*i + k] - seconddiffer[ 3*i + k ]*w2list[ i ];
				//	newver[ 3 * i + k ] = oldver[ 3 * i + k ] + ratio * oldfirstdiffer[ 3 * i + k ] ;
			}
		}
		else
		{
			// Manifold points: new averaging
			for ( int k = 0 ; k < 3 ; k ++ )
			{
				newver[ 3 * i + k ] = oldver[ 3 * i + k ] - ratio * seconddiffer[ 3 * i ] * norms[ 3 * i + k ] + ratio * vecs[ 3 * i + k ] ;
			}
/*
			double tjj=(newver[3*i+0]-oldver[3*i+0])*(newver[3*i+0]-oldver[3*i+0]);
			tjj+=(newver[3*i+1]-oldver[3*i+1])*(newver[3*i+1]-oldver[3*i+1]);
			tjj+=(newver[3*i+2]-oldver[3*i+2])*(newver[3*i+2]-oldver[3*i+2]);
			if(tjj<0.0001) tmove++;
*/
		}
	}
/*
	cout<<"tmove:"<<tmove<<endl;
*/
	delete vecs ;
	delete dis ;
	w2list.clear();
}


void inline LuMesh::computeWeightList(vector<int>&edgelist, vector<double>& wlist)
{
	//inverseedgelenslist
	int edgelen = edgelist.size()/2;
	////wlist.resize(edgelen);
	//double temp;
	//for( int i = 0; i < edgelen; i ++)
	//{
	//	temp = MyMath::vectorlen(&sufver[ 3 * edgelist[ i * 2 ]], &sufver[ 3* edgelist[ i * 2 + 1]]);;
	//	if( temp == 0)
	//	{
	//		wlist[ i ] = 1000000;
	//		cout<<"Edge length is 0 when computing the weight!"<<endl;
	//	}
	//	else
	//		wlist[ i ] = 1/temp;
	//}
	//////////////////////////////////////////////////////////////////////////
	for( int i = 0; i < edgelen; i ++)
	{
		wlist[ i ] = 1;
	}
}

void inline LuMesh::computeDiffer(double* oldval, double* differ, vector<double>& wlist, vector<int>&vermark,
								vector<intvector>&verneighbr, HashMap& ver2edgehash)
{
	memset( differ, 0, sizeof(double)*3*sufvernum);
	//only compute difference for non-manifold and normal vertex
	for( int i = 0 ; i < sufvernum ; i++)
	{
		double wsum = 0;
		for( int j = 0; j < verneighbr[i].size(); j++)	//all neighbors
		{
			int ver2 = verneighbr[ i ][ j ];
			if( (vermark[ver2] == 1) || (vermark[ver2] == 3))	//contour vertex 
				continue;
			int edgei = ver2edgehash.findInsertSort( i, ver2, -1);
			if( edgei == -1 )
			{
				cout<<"Error in SDUmbraFair when finding current edge!"<<endl;
				continue;
			}
			double curw = wlist[ edgei ];
			wsum += curw;
			for( int k = 0; k < 3; k ++)	//for three coord
			{
				differ[ 3 * i + k ] += (curw * oldval[ 3 * ver2+ k]);
			}
		}
		for( int k = 0; k < 3; k ++)
			differ[ 3 * i + k ] = differ[ 3 * i + k ]/wsum - oldval[ 3 * i + k ];
	}
}
void inline LuMesh::computeDiffer(double* oldval, double* differ, vector<double>& wlist, 
								vector<intvector>&verneighbr, HashMap& ver2edgehash)
{
	memset( differ, 0, sizeof(double)*3*sufvernum);
	for( int i = 0 ; i < sufvernum ; i++)
	{
		//	for( int k = 0; k < 3; k ++)
		//		differ[ 3*i + k ] = 0;
		double wsum = 0;
		for( int j = 0; j < verneighbr[i].size(); j++)	//all neighbors
		{
			int ver2 = verneighbr[ i ][ j ];
			int edgei = ver2edgehash.findInsertSort( i, ver2, -1);
			if( edgei == -1 )
			{
				cout<<"Error in SDUmbraFair when finding current edge!"<<endl;
				continue;
			}
			double curw = wlist[ edgei ];
			wsum += curw;
			for( int k = 0; k < 3; k ++)	//for three coord
			{
				differ[ 3 * i + k ] += (curw * oldval[ 3 * ver2+ k]);
			}
		}
		for( int k = 0; k < 3; k ++)
			differ[ 3 * i + k ] = differ[ 3 * i + k ]/wsum - oldval[ 3 * i + k ];
	}
}
void inline LuMesh::computeDiffer(vector<double>& oldval, vector<double>& differ, vector<double>& wlist, 
								vector<intvector>&verneighbr, HashMap& ver2edgehash)
{
	for( int i = 0 ; i < sufvernum ; i++)
	{
		for( int k = 0; k < 3; k ++)
			differ[ 3*i + k ] = 0;
		double wsum = 0;
		for( int j = 0; j < verneighbr.size(); j++)	//all neighbors
		{
			int ver2 = verneighbr[ i ][ j ];
			int edgei = ver2edgehash.findInsertSort( i, ver2, -1);
			if( edgei == -1 )
			{
				cout<<"Error in SDUmbraFair when finding current edge!"<<endl;
				continue;
			}
			double curw = wlist[ edgei ];
			wsum += curw;
			for( int k = 0; k < 3; k ++)	//for three coord
			{
				differ[ 3 * i + k ] += (curw * oldval[ 3 * ver2+ k]);
			}
		}
		for( int k = 0; k < 3; k ++)
			differ[ 3 * i + k ] = differ[ 3 * i + k ]/wsum - oldval[ 3 * i + k ];

	}
}

void inline LuMesh::computeNorm(double* oldval, double* norms, double* mags )
{
	memset( norms, 0, sizeof(double)*3*sufvernum);

	for( int i = 0 ; i < suffacenum ; i++)
	{
		double vec1[3], vec2[3], nm[3];
		MyMath::getVec(&(oldval[ sufface[ i * 3] * 3 ]), &(oldval[ sufface[ i * 3 + 1] * 3 ]), vec1 );
		MyMath::getVec(&(oldval[ sufface[ i * 3] * 3 ]), &(oldval[ sufface[ i * 3 + 2] * 3 ]), vec2 );
		MyMath::crossProductNotNorm( vec1, vec2, nm);

		if ( sufmat[i*2] > sufmat[i*2 + 1] )
		{
			nm[0] = -nm[0] ;
			nm[1] = -nm[1] ;
			nm[2] = -nm[2] ;
		}

		for ( int j = 0 ; j < 3 ; j ++ )
		{
			int ind = sufface[ i * 3 + j ] ;
			// for three vertices of the triangle
			// if ( vermark[ ind ] < 2 )
			{
				for( int k = 0; k < 3; k ++)	// for three coord
				{
					norms[ 3 * ind + k ] += nm[ k ] ;
				}
			}
		}
		// printf("%d %d\n", sufmat[ i * 2 ], sufmat[ i * 2 + 1 ] );
	}
	for ( int i = 0 ; i < sufvernum ; i ++ )
	{
		mags[ i ] = sqrt( MyMath::vectorlen( &(norms[ 3 * i ]) ) );
		// printf("%f\n", mags[i]) ;
	}

}


void inline LuMesh::computeRefreshWeightList(vector<intvector>&verneighbr, vector<double>& w2list)
{
	// 1/(1 + 1/ni(signma /1nij))
	int* valence = new int[ sufvernum ];
	for( int i = 0; i < sufvernum; i ++)
		valence[ i ] = verneighbr[ i ].size();
	//////////////////////////////////////////////////////////////////////////
	//	for( int i = 0; i < sufvernum; i++)
	//		cout<<valence[i]<<" ";
	//	cout<<endl;
	//////////////////////////////////////////////////////////////////////////
	for( int i = 0; i < sufvernum; i ++)
	{
		w2list[ i ] = 0;
		for( int j = 0; j < verneighbr[ i ].size(); j++)
		{
			w2list[ i ] += (1/(double)valence[verneighbr[i][j]]);
		}
		w2list[ i ] /= valence[ i ];
		w2list[ i ] += 1;
		w2list[ i ] = 1/w2list[ i ];
		if( w2list[ i ] > 0.4)
			w2list[ i ] = 0.4;
		//	cout<<w2list[i]<<" ";
	}
	//	cout<<endl;
}





void LuMesh::LiepaRefine(double alpha)
{
	/*edge*/
	vector<int> edgelist;
	//	edgelist.clear();
	//	edge2facelist.clear();
	vector<intvector> edge2facelist;
	vector<int> ctredgelist;
	vector<int>	nmedgelist;
	//	nmedgelist.clear();
	vector<int> normedgelist;
	vector<double> edgelenlist;
	cout<<"gather edge information...."<<"\t";
	HashMap ver2edgehash;
	gatherEdgeInfo(edgelist, edge2facelist,ver2edgehash,ctredgelist,nmedgelist,normedgelist,edgelenlist);
	cout<<"done!"<<endl;
/*
	cout<<"edgelist len:"<<edgelist.size()/2<<endl;
	cout<<"ctredgelist len:"<<ctredgelist.size()<<endl;
	cout<<"nmedgelist len:"<<nmedgelist.size()<<endl;
	cout<<"normedgelist len:"<<normedgelist.size()<<endl;
*/

	/*ver*/
	cout<<"gather vertex information......"<<"\t";
	vector<double> verlist;
	vector<double> verattrlist;
	vector<int> vermark;
	gatherVerInfo( verlist,verattrlist,vermark, edgelist, edgelenlist,nmedgelist, normedgelist,ctredgelist );
	int oldedgenum = edgelist.size()/2;		//before splitting the non-manifold edge, for the splitting those non-contour edges with two contour vertices
	//////////////////////////////////////////////////////////////////////////
	cout<<"done!"<<endl;

	cout<<"verlist len:"<<verlist.size()/3<<endl;
	cout<<"edgelist len:"<<edgelist.size()/2<<endl;
	cout<<"ctredgelist len:"<<ctredgelist.size()<<endl;
	cout<<"nmedgelist len:"<<nmedgelist.size()<<endl;
	cout<<"normedgelist len:"<<normedgelist.size()<<endl;

	/*face*/
	vector<int> facelist;
	gatherFaceInfo(facelist);
	cout<<"facelist len:"<<facelist.size()/5<<endl;

	/*---debug--*/
	//dbedgelist.clear();
	//dbnmedgelist.clear();
	//dbedgelist.resize(edgelist.size());
	//for( int i =0; i < edgelist.size(); i ++)
	//{
	//	dbedgelist[ i ] = edgelist[ i ];
	//}
	//dbnmedgelist.resize( nmedgelist.size());
	//for( int i = 0; i < nmedgelist.size(); i ++)
	//{
	//	dbnmedgelist[ i ] = nmedgelist[ i ];
	//}
	////find the edge that connects the third vertex and nmvertex together
	//int dbedgelen = edgelist.size()/2;
	//for( int i = 0; i < nmedgelist.size(); i ++)
	//{
	//	int edgei = nmedgelist[ i ];
	//	int verpos[ 2 ] = {edgelist[ edgei * 2 ], edgelist[ edgei * 2  + 1 ]};
	//	intvector& curfacelist = edge2facelist[ edgei ];
	//	for( int j = 0; j < curfacelist.size()/2; j ++)
	//	{
	//		int facei = curfacelist[ j * 2  ];
	//		int edgepos = curfacelist[ j * 2 + 1 ];
	//		int thirdver = facelist[ 5*facei + (edgepos + 2)%3 ];
	//		for( int k = 0; k < 2; k ++)
	//		{                
	//			int edgeindex = ver2edgehash.findInsertSort( thirdver, verpos[ k ], dbedgelen );
	//			if( edgeindex == dbedgelen )
	//			{
	//				cout<<"HASH IS NOT CORRECTLY BUILT!!!!"<<endl;
	//			}
	//			else
	//				dbnmedgelist.push_back( edgeindex);
	//		}
	//	}
	//}
	/*------*/
	/*-----------debug----------*/
	/*dbnmedgelist.resize(nmedgelist.size()*2);
	for( int i = 0; i < nmedgelist.size(); i++)
	{
	dbnmedgelist[ i * 2] = edgelist[nmedgelist[ i ]*2];
	dbnmedgelist[ i * 2 + 1] = edgelist[nmedgelist[ i ]*2 + 1];
	}*/
	/*------------------------------*/



	/*split non-manifold edge*/
	cout<<"splitting non-manifold edge.....\t";
	splitNMEdge(verlist,verattrlist,edgelist,edge2facelist,nmedgelist,normedgelist,facelist,ver2edgehash);
	cout<<"done!"<<endl;

	/*split edges with two convex vertices*/
	//	splitNCtrEdgeTwoCtrVer( verlist, verattrlist,vermark, edgelist, edge2facelist,nmedgelist, normedgelist,facelist,ver2edgehash, oldedgenum);
	//
	//	//-------------------///
	//	/*for( int i  = 10; i < 25; i ++)
	//	{
	//		cout<<"edge "<<normedgelist[i]<<"\n";
	//		for( int j = 0; j < 2; j ++)
	//		{
	//			int facei = edge2facelist[ normedgelist[ i ] ][ 2*j ];
	//			cout<<facei<<"\t"<<edge2facelist[ normedgelist[ i ] ][ 2*j + 1]<<"\t";
	//			for( int k = 0; k < 5; k ++)
	//			{
	//				cout<<facelist[ 5*facei + k ]<<"\t";
	//			}
	//			cout<<endl;
	//		}
	//	}*/
	//
	/*swap edges*/
	HashMap ver2edgehash2;
	for( int i = 0; i < edgelist.size()/2; i ++)
	{
		if( ver2edgehash2.findInsertSort( edgelist[ i * 2 ], edgelist[ i * 2 + 1], i ) != i)
			cout<<"Error occurs when building hash map for edges!"<<endl;
	}

	if(doFirstSwap){
		cout<<"swap edges ....."<<"\t";
		swapEdge(verlist,edgelist,edge2facelist,normedgelist,facelist,ver2edgehash2);
		cout<<"done!"<<endl;
	}


	/*edge type*/
	//	getEdgeTypeList(normedgelist,nmedgelist,ctredgelist,edgetypelist);

	/*insert new vertex and swap*/
	cout<<"split triangles ....."<<"\t";
	int count = 0;
	while( count < 1024 )
	{
		count++;
		if(splitTriangle(verlist, verattrlist,ver2edgehash2,edgelist, normedgelist,edge2facelist,facelist, alpha))
		{
			swapEdge(verlist,edgelist,edge2facelist,normedgelist,facelist,ver2edgehash2);
		}
		else
			break;
	}
	cout<<"done!"<<endl;

	//	//////////////////////////////////////////////////////////////////////////
	//	/*dbnmedgelist.clear();
	//	for( int i =  0;  i < normedgelist.size(); i ++)
	//	{
	//		if( edge2facelist[ normedgelist[ i ]].size() < 4 )
	//		{
	//			dbnmedgelist.push_back( edgelist[normedgelist[ i ]*2]);
	//			dbnmedgelist.push_back( edgelist[normedgelist[ i ]*2 + 1]);
	//		}
	//	}*/
	//
	/*put the new mesh into array*/
	//vertex ctredge face facemat
	cout<<"refresh mesh...."<<"\t";
	sufvernum = verlist.size()/3;
	delete []sufver;
	sufver = new double[ sufvernum*3];
	for( int i = 0; i < sufvernum*3; i++)
		sufver[ i ] = verlist[i];
	suffacenum = facelist.size()/5;
	delete []sufface;
	delete []sufmat;
	delete []suffacenorm;
	sufface = new int[ suffacenum * 3];
	suffacenorm = new double[suffacenum*3];
	sufmat = new int[ suffacenum * 2];
	double vec1[3],vec2[3];
	for( int i = 0; i < suffacenum; i ++)
	{
		for( int j = 0; j < 3; j ++)
		{
			sufface[ i*3 + j] = facelist[ i * 5 + j];
		}
		for( int j = 0; j < 2; j ++)
		{
			sufmat[ i * 2 + j ] = facelist[i * 5 + 3 + j];
		}
		MyMath::getVec(&sufver[ sufface[ i * 3] * 3 ], &sufver[ sufface[ i * 3 + 1] * 3 ], vec1 );
		MyMath::getVec(&sufver[ sufface[ i * 3] * 3 ], &sufver[ sufface[ i * 3 + 2] * 3 ], vec2 );
		MyMath::crossProduct( vec1, vec2, &suffacenorm[ i * 3 ]);
	}

	/*	sufctredgenum = ctredgelist.size();
	delete []sufctredge;
	sufctredge = new int[ sufctredgenum*2 ];
	for( int i = 0; i < sufctredgenum; i ++)
	{
	sufctredge[ i * 2 ] = edgelist[ ctredgelist[ i ] * 2 ];
	sufctredge[ i * 2 + 1] = edgelist[ ctredgelist[ i ] * 2 + 1];
	} */
	cout<<"done!"<<endl;  

	//debug
	//	dbnmedgelist.clear();
	//	dbnmedgelist.resize( nmedgelist.size());
	//for( int i = 0; i < nmedgelist.size(); i++)
	//{
	//	dbnmedgelist.push_back(edgelist[nmedgelist[i]*2]);
	//	dbnmedgelist.push_back(edgelist[nmedgelist[i]*2+1]);
	//}

	/*release the memory*/
	for( int i = 0; i < edge2facelist.size(); i ++)
		edge2facelist[i].clear();
	edge2facelist.clear();
	edgelist.clear();
	ctredgelist.clear();
	nmedgelist.clear();
	normedgelist.clear();
	edgelenlist.clear();
	verlist.clear();
	verattrlist.clear();
	vermark.clear();

/*
	cout<<"sufvernum len:"<<sufvernum<<endl;
	cout<<"suffacenum len:"<<suffacenum<<endl;
*/
}
void LuMesh::gatherVerInfo(vector<double>& verlist,vector<double>& verattrlist, vector<int>&edgelist, vector<double>&edgelenlist,
						 vector<int>&nmedgelist, vector<int>&normedgelist, vector<int>&ctredgelist )
{
	//verlist
	verlist.resize(sufvernum*3, 0);
	for( int i = 0; i < sufvernum * 3; i ++)
		verlist[ i ] = sufver[ i ];

	//verattrilist: the expected lengh at that vertex
	//	cout<<"sufvernum:"<<sufvernum;
	verattrlist.resize(sufvernum, 0);
	int* vermark = new int[ sufvernum ];	//0 - normal ver 1 - ver on contour edge 2 - ver on non manifold edge not on contour edge
	int* vercount = new int[ sufvernum ];
	for( int i = 0; i < sufvernum; i ++)
	{
		vermark[ i ] = 0;
		vercount[ i ] = 0;
		verattrlist[ i ] = 0;
	}
	//cout<<"ctredgelist.size:"<<ctredgelist.size()<<endl;
	for( int i = 0;i < ctredgelist.size(); i ++)
	{
		vermark[ edgelist[ ctredgelist[ i ] * 2 ] ] = 1;
		vermark[ edgelist[ ctredgelist[ i ] * 2 + 1]] = 1;
	}
	for( int i = 0; i < nmedgelist.size(); i ++)
	{
		if( vermark[ edgelist[ nmedgelist[ i ]* 2 ] ] == 0 )
			vermark[ edgelist[ nmedgelist[ i ]* 2 ] ] = 2;
		if( vermark[ edgelist[ nmedgelist[ i ]* 2  + 1] ] == 0 )
			vermark[ edgelist[ nmedgelist[ i ]* 2  + 1] ] = 2;
	}

	int veri = 0;
	//	cout<<"ctredgelist.size:"<<ctredgelist.size()<<endl;
	for( int i = 0; i < ctredgelist.size(); i ++)
	{
		for( int j = 0;  j < 2; j ++)
		{
			veri = edgelist[ctredgelist[ i ]*2 + j ];
			//	cout<<"veri:"<<veri<<endl;
			//	cout<<"ctredgelist[ i ]:"<<ctredgelist[ i ]<<endl;
			//	cout<<"edgelenlist[ ctredgelist[ i ] ]"<<edgelenlist[ ctredgelist[ i ] ]<<endl;
			verattrlist[ veri ] = verattrlist[ veri ] + edgelenlist[ ctredgelist[ i ] ];
			//	cout<<edgelenlist[ ctredgelist[ i ] ]<<"\t";
			vercount[ veri ] ++;
		}
	}
	//	cout<<endl;
	for( int i = 0; i < sufvernum; i ++)
	{
		if( vermark[ i ] == 1)			//contour vertex
		{
			verattrlist[ i ] /= vercount[ i ];
			vercount[ i ] = 0;
		}
		//cout<<verattrlist[ i ]<<"\t";
	}
	//attribute for other vertices: ring by ring propagating the attributes
	int veri1, veri2;
	int edgenum = edgelist.size()/2;
	bool again = true;
	while( again )
	{
		//	memset(vercount, 0, sizeof(int)*sufvernum);
		again = false;
		for( int i = 0; i < edgenum; i ++)
		{
			veri1 = edgelist[ i * 2 ];
			veri2 = edgelist[ i * 2 + 1 ];
			if( ((vercount[ veri1 ]== 0) && (verattrlist[ veri1 ] != 0) ) && ((verattrlist[ veri2 ] == 0) ||(vercount[veri2] != 0 ) ) )		//veri1 is one that has attri, while veri2 not
			{
				verattrlist[ veri2 ] += verattrlist[ veri1 ];
				vercount[ veri2 ] ++;
			}
			else if( ((vercount[ veri2 ]== 0) && (verattrlist[ veri2 ] != 0) ) && ((verattrlist[ veri1 ] == 0) ||(vercount[veri1] != 0 ) )  )	//veri2 has while veri1 not
			{
				verattrlist[ veri1 ] += verattrlist[ veri2 ];
				vercount[ veri1 ] ++;
			}
		}

		for( int i = 0; i < sufvernum; i ++)
		{
			if( vercount[ i ] != 0 )
			{
				again = true;
				verattrlist[ i ] /= vercount[ i ];
				vercount[ i ] = 0;
			}
		}
	}
	//attribute for other vertices: old implementation
	//	int veri2 = 0;
	////	cout<<"normedgelist.size:"<<normedgelist.size()<<endl;
	//	for( int i = 0; i < normedgelist.size(); i ++)
	//	{
	//		veri = edgelist[ normedgelist[i] * 2];
	//		veri2 = edgelist[ normedgelist[i] * 2 + 1];
	//		if( vermark[ veri ] == 1)
	//		{
	//			verattrlist[ veri2 ] += verattrlist[ veri ];
	//			vercount[ veri2 ] ++;
	//		}
	//		if( vermark[ veri2 ] == 1)
	//		{
	//			verattrlist[ veri ] += verattrlist[ veri2 ];
	//			vercount[ veri ] ++;
	//		}
	//	}
	////	cout<<"nmedgelist.size:"<<nmedgelist.size()<<endl;
	//	for( int i = 0; i < nmedgelist.size(); i++)
	//	{
	//		veri = edgelist[ nmedgelist[i] * 2];
	//		veri2 = edgelist[ nmedgelist[i] * 2 + 1];
	//		if( vermark[ veri ] == 1)
	//		{
	//			verattrlist[ veri2 ] += verattrlist[ veri ];
	//			vercount[ veri2 ] ++;
	//		}
	//		if( vermark[ veri2 ] == 1)
	//		{
	//			verattrlist[ veri ] += verattrlist[ veri2 ];
	//			vercount[ veri ] ++;
	//		}
	//	}
	//	for( int i = 0; i < sufvernum; i ++)
	//	{
	//		if( vermark[ i ] != 1 && vercount[ i ] != 0 )
	//			verattrlist[ i ] /= vercount[ i ];
	//
	//	//	cout<<verattrlist[i]<<"\t";
	//	}
	//	double maxattr,minattr;
	//	maxattr = minattr = verattrlist[ edgelist[ctredgelist[ 0 ] * 2 ]];
	//	for( int i = 0; i < sufvernum; i ++)
	//	{
	//		if( vercount[ i ]!=0)
	//		{
	//			if( verattrlist[ i ] < minattr)
	//				minattr = verattrlist[ i ];
	//			else if( verattrlist[ i ] > maxattr )
	//				maxattr = verattrlist[ i ];
	//		}
	//	}
	//	//maxattr = (maxattr + minattr)/2;
	////	cout<<"maxattr"<<maxattr<<endl;
	//	for( int i = 0; i < sufvernum; i ++)
	//	{
	//		if( vercount[ i ] == 0)
	//			verattrlist[ i ] = maxattr;
	//	}
	//	cout<<endl;	
	delete []vermark;
}
void LuMesh::gatherVerInfo(vector<double>& verlist,vector<double>& verattrlist,vector<int>&vermark,  vector<int>&edgelist, vector<double>&edgelenlist,
						 vector<int>&nmedgelist, vector<int>&normedgelist, vector<int>&ctredgelist )
{
	//verlist
	verlist.resize(sufvernum*3, 0);
	for( int i = 0; i < sufvernum * 3; i ++)
		verlist[ i ] = sufver[ i ];

	//verattrilist: the expected lengh at that vertex
	//	cout<<"sufvernum:"<<sufvernum;
	verattrlist.resize(sufvernum, 0);
	//int* vermark = new int[ sufvernum ];	//0 - normal ver 1 - ver on contour edge 2 - ver on non manifold edge not on contour edge
	vermark.resize(sufvernum);
	int* vercount = new int[ sufvernum ];
	for( int i = 0; i < sufvernum; i ++)
	{
		vermark[ i ] = 0;
		vercount[ i ] = 0;
		verattrlist[ i ] = 0;
	}
	//cout<<"ctredgelist.size:"<<ctredgelist.size()<<endl;
	for( int i = 0;i < ctredgelist.size(); i ++)
	{
		vermark[ edgelist[ ctredgelist[ i ] * 2 ] ] = 1;
		vermark[ edgelist[ ctredgelist[ i ] * 2 + 1]] = 1;
	}
	for( int i = 0; i < nmedgelist.size(); i ++)
	{
		if( vermark[ edgelist[ nmedgelist[ i ]* 2 ] ] == 0 )
			vermark[ edgelist[ nmedgelist[ i ]* 2 ] ] = 2;
		if( vermark[ edgelist[ nmedgelist[ i ]* 2  + 1] ] == 0 )
			vermark[ edgelist[ nmedgelist[ i ]* 2  + 1] ] = 2;
	}

	int veri = 0;
	//	cout<<"ctredgelist.size:"<<ctredgelist.size()<<endl;
	for( int i = 0; i < ctredgelist.size(); i ++)
	{
		for( int j = 0;  j < 2; j ++)
		{
			veri = edgelist[ctredgelist[ i ]*2 + j ];
			//	cout<<"veri:"<<veri<<endl;
			//	cout<<"ctredgelist[ i ]:"<<ctredgelist[ i ]<<endl;
			//	cout<<"edgelenlist[ ctredgelist[ i ] ]"<<edgelenlist[ ctredgelist[ i ] ]<<endl;
			verattrlist[ veri ] = verattrlist[ veri ] + edgelenlist[ ctredgelist[ i ] ];
			//	cout<<edgelenlist[ ctredgelist[ i ] ]<<"\t";
			vercount[ veri ] ++;
		}
	}
	//	cout<<endl;
	for( int i = 0; i < sufvernum; i ++)
	{
		if( vermark[ i ] == 1)			//contour vertex
		{
			verattrlist[ i ] /= vercount[ i ];
			vercount[ i ] = 0;
		}
		//cout<<verattrlist[ i ]<<"\t";
	}
	//attribute for other vertices: ring by ring propagating the attributes
	int veri1, veri2;
	int edgenum = edgelist.size()/2;
	bool again = true;
	while( again )
	{
		//	memset(vercount, 0, sizeof(int)*sufvernum);
		again = false;
		for( int i = 0; i < edgenum; i ++)
		{
			veri1 = edgelist[ i * 2 ];
			veri2 = edgelist[ i * 2 + 1 ];
			if( ((vercount[ veri1 ]== 0) && (verattrlist[ veri1 ] != 0) ) && ((verattrlist[ veri2 ] == 0) ||(vercount[veri2] != 0 ) ) )		//veri1 is one that has attri, while veri2 not
			{
				verattrlist[ veri2 ] += verattrlist[ veri1 ];
				vercount[ veri2 ] ++;
			}
			else if( ((vercount[ veri2 ]== 0) && (verattrlist[ veri2 ] != 0) ) && ((verattrlist[ veri1 ] == 0) ||(vercount[veri1] != 0 ) )  )	//veri2 has while veri1 not
			{
				verattrlist[ veri1 ] += verattrlist[ veri2 ];
				vercount[ veri1 ] ++;
			}
		}

		for( int i = 0; i < sufvernum; i ++)
		{
			if( vercount[ i ] != 0 )
			{
				again = true;
				verattrlist[ i ] /= vercount[ i ];
				vercount[ i ] = 0;
			}
		}
	}
	//	delete []vermark;
}

void LuMesh::gatherEdgeInfo(vector<int>& edgelist,vector<intvector>& edge2facelist,
						  HashMap& ver2edgehash,vector<int>& ctredgelist,vector<int>&	nmedgelist,
						  vector<int>& normedgelist,vector<double>& edgelenlist)
{
	/*---------------------*/
	/*edgelist.clear();
	for(int i = 0; i < edge2facelist.size(); i ++)
	edge2facelist[ i ].clear();
	edge2facelist.clear();
	ctredgelist.clear();
	nmedgelist.clear();
	normedgelist.clear();
	edgelenlist.clear();*/
	/*---------------------*/

	//edgelist and edge2facelist and ver2edgehash
	int edgelen = 0;
	int edgeindex;
	for( int i = 0;  i < suffacenum; i++)
	{
		for( int j = 0; j < 3; j ++)
		{
			edgeindex = ver2edgehash.findInsertSort(sufface[ i*3 + j ], sufface[ i*3 + (j+1)%3 ], edgelen);
			if(edgeindex == edgelen)
			{
				edgelen ++;
				edgelist.push_back(sufface[ i*3 + j ]);
				edgelist.push_back(sufface[ i*3 + (j+1)%3 ]);
				intvector facelist;
				facelist.push_back( i );	//i face
				facelist.push_back( j );	//j edge in i face
				edge2facelist.push_back( facelist );	//in c++, copy and push back
				facelist.clear();		//clear temp variable, avoid memory leakage
			}
			else
			{
				edge2facelist[ edgeindex ].push_back( i );	//face i
				edge2facelist[ edgeindex ].push_back( j );	//edge j in face i
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//	cout << "edgelen:"<<edgelen<<"edgelist size:"<<edgelist.size()/2<<endl;
	//////////////////////////////////////////////////////////////////////////
	//ctredgelist nmedgelist normedgelist
	int* edgemark = new int[ edgelen ];		// 0 - normal edge 1 - contour edge 2 - nonmanifold edge
	for( int i= 0; i < edgelen; i++)
		edgemark[ i ] = 0;
	//////////////////////////////////////////////////////////////////////////
	//	cout<<"sufctredgenum:"<<sufctredgenum<<endl;
	for( int i = 0; i < sufctredgenum; i++)
	{
		edgeindex = ver2edgehash.findInsertSort( sufctredge[ i * 2], sufctredge[i * 2 +1], edgelen);
		//cout<<"edgeindex:"<<edgeindex<<endl;
		if( edgeindex == edgelen )
		{
			cout<<sufctredge[ i * 2]<<"\t"<<sufctredge[i * 2 +1]<<endl;
			cout<<"not valid input!"<<endl;
			continue;
		}
		edgemark[edgeindex] = 1;
	}
	for( int i = 0; i < edgelen; i ++)
	{
		if( edge2facelist[ i ].size() > 4)
		{
			edgemark[ i ] = 2;
		}
	}
	for(int i = 0; i < edgelen; i ++)
	{
		switch( edgemark[i] )
		{
		case 0:
			normedgelist.push_back( i );
			break;
		case 1:
			ctredgelist.push_back( i );
			break;
		case 2:
			nmedgelist.push_back( i );
			break;
		default:
			break;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//cout<<"contour edge number:"<<ctredgelist.size()<<"\t";
	//cout<<"non manifold edge number:"<<nmedgelist.size()<<endl;
	delete []edgemark;

	//edgelenlist
	edgelenlist.resize(edgelen,0);
	//	for( int i = 0; i < edgelen; i ++)
	//		edgelenlist[i] = 0;
	for(int i = 0; i < ctredgelist.size(); i++)
	{
		edgelenlist[ ctredgelist[i]] = MyMath::vectorlen(&sufver[edgelist[ctredgelist[i]*2]*3],
			&sufver[edgelist[ctredgelist[i]*2 + 1]*3]);
	}	
}
void LuMesh::gatherFaceInfo(vector<int>& facelist)
{
	facelist.resize(suffacenum*5, 0);
	for( int i = 0; i < suffacenum; i ++)
	{
		for( int j = 0; j < 3; j ++)
		{
			facelist[ i * 5 + j ] = sufface[ i * 3 + j];
		}
		for( int j = 0; j < 2; j ++)
			facelist[ i * 5 + j + 3 ] = sufmat[ i * 2 + j];
	}
}

void LuMesh::splitOneNMEdge(int segnum, int i, vector<double>& verlist,vector<double>& verattrlist,vector<int>& edgelist,vector<intvector>& edge2facelist,
						  vector<int>&nmedgelist, vector<int>& normedgelist,vector<int>&facelist, HashMap& ver2edgehash)
{
	int verpos[2] = {edgelist[nmedgelist[i]*2],edgelist[nmedgelist[i]*2 + 1]};
	//	cout<<"verpos are:"<<verpos[0]<<"   "<<verpos[1]<<endl;
	double newver[3];
	int verlen = verlist.size()/3;
	for( int j = 1; j < segnum; j++ )
	{
		MyMath::getPtOnSeg(&verlist[verpos[0]*3], &verlist[verpos[1]*3], j/(double)segnum, newver);
		verlist.push_back( newver[ 0 ]) ;
		verlist.push_back( newver[ 1 ]);
		verlist.push_back( newver[ 2 ]);
		verattrlist.push_back( (1 - j/(double)segnum)*verattrlist[verpos[0]] +  j/(double)segnum*verattrlist[verpos[1]] );
	}
	//------------//
	//	cout<<"vertex added!"<<endl;
	//new edge on the old non-manifold edge
	int edgelen = edgelist.size()/2;
	int oldnmedgelen = nmedgelist.size();
	nmedgelist.resize(oldnmedgelen + segnum - 1);
	for( int j = 1; j < segnum - 1; j ++)
	{
		edgelist.push_back( verlen + j - 1);
		edgelist.push_back( verlen + j);
		nmedgelist[ oldnmedgelen + j - 1] = edgelen;
		ver2edgehash.findInsertSort( verlen + j -1, verlen + j, edgelen);
		edgelen ++;
	}
	nmedgelist[ oldnmedgelen + segnum - 2] = edgelen;
	edgelist.push_back( verlen + segnum - 2);	//the last segment
	edgelist.push_back( verpos[1] );
	ver2edgehash.findInsertSort( verlen + segnum - 2, verpos[1], edgelen);
	edgelen++;
	edgelist[ nmedgelist[ i ]*2 + 1 ] = verlen;	//second point of the first segment

	int facenum = facelist.size()/5;	//old face number
	//edgelen = edgelist.size()/2;		//edge number after adding new edges on non-manifold edge
	//------------//
	//	cout<<"new edge added!"<<endl;
	//every face
	int curfacenum = edge2facelist[ nmedgelist[ i ]].size()/2;	//incident faces number
	//intvector& curfacelist = edge2facelist[ nmedgelist[i]];		//incident faces list	/*-------*/

	int oldnormedgelen = normedgelist.size();
	normedgelist.resize( oldnormedgelen + curfacenum*(segnum - 1));
	for( int j = 0; j < (segnum - 1)*curfacenum; j++)
		normedgelist[ oldnormedgelen + j] = edgelen + j;
	edge2facelist.resize(edgelen + curfacenum * (segnum - 1));

	intvector& curfacelist = edge2facelist[ nmedgelist[i]];		//incident faces list	/*-------*/
	//------------//
	//		cout<<"incident face number:"<<curfacenum<<endl;
	for( int j = 0; j < curfacenum; j ++)
	{
		//	cout<<"curface:"<<curfacelist[2*j]<<endl;
		int edgepos = curfacelist[ 2*j+1 ];
		int thirdverpos = facelist[5 * curfacelist[ 2*j ]+(edgepos + 2)%3];
		//	cout<<"current face:"<<facelist[5 * curfacelist[ 2*j ]]<<" "<<facelist[5 * curfacelist[ 2*j ] + 1]<<
		//			" "<<facelist[5 * curfacelist[ 2*j ] + 2]<<" "<<facelist[5 * curfacelist[ 2*j ]+ 3]<<" "
		//			<<facelist[5 * curfacelist[ 2*j ]+4]<<endl;
		//	cout<<"edgepos"<<edgepos<<"\t";
		//	cout<<"thirdver:"<<thirdverpos<<"\n";
		bool iscw = (verpos[ 0 ] == facelist[ 5 * curfacelist[ 2*j ] + edgepos ]);
		int mat[2] = {facelist[ 5 * curfacelist[ 2*j ] + 3], facelist[ 5 * curfacelist[ 2*j ] + 4]};

		//add new edge, add new face, refresh edge2facelist
		//add new edge
		for( int k = 0; k < segnum - 1; k++ )
		{
			edgelist.push_back( thirdverpos );
			edgelist.push_back( verlen + k );
			ver2edgehash.findInsertSort( thirdverpos, verlen + k,  edgelen + k + j * (segnum-1) );
		}
		//------------//
		//	cout<<"new edge added for current face!"<<endl;
		//new face
		if( iscw )
		{
			for( int k = 0; k < segnum - 2; k ++ )
			{
				facelist.push_back( thirdverpos );
				facelist.push_back( verlen + k );
				facelist.push_back( verlen + k + 1);
				facelist.push_back( mat[ 0 ]);
				facelist.push_back( mat[ 1 ]);
			}
			for( int k = 0; k < 5; k ++)
				facelist.push_back( facelist[ 5 * curfacelist[ 2*j ] + k ]);
			facelist[ ( facenum + segnum - 2 ) * 5 + edgepos ] = verlen + segnum - 2;
			facelist[ 5 * curfacelist[ 2*j ] + (edgepos + 1)%3 ] = verlen;
		}
		else
		{
			for( int k = 0; k < segnum - 2; k ++)
			{
				facelist.push_back( thirdverpos );
				facelist.push_back( verlen + k + 1);
				facelist.push_back( verlen + k );					
				facelist.push_back( mat[ 0 ]);
				facelist.push_back( mat[ 1 ]);
			}
			for( int k = 0; k < 5; k ++)
				facelist.push_back( facelist[ 5 * curfacelist[ 2*j ] + k ]);
			facelist[ ( facenum + segnum - 2 ) * 5 + (edgepos + 1)%3] = verlen + segnum - 2;
			facelist[ 5 * curfacelist[ 2*j ] + edgepos ] = verlen;
		}
		//---------------------//
		/*cout<<"faces:"<<facelist[5 * curfacelist[ 2*j ]]<<" "<<facelist[5 * curfacelist[ 2*j ] + 1]<<
		" "<<facelist[5 * curfacelist[ 2*j ] + 2]<<" "<<facelist[5 * curfacelist[ 2*j ]+ 3]<<" "
		<<facelist[5 * curfacelist[ 2*j ]+4]<<endl;
		for( int k = 0; k < segnum - 1; k++ )
		{
		cout<<"faces:"<<facelist[5 * (facenum+k)]<<" "<<facelist[5 *(facenum+k) + 1]<<
		" "<<facelist[5 *(facenum+k) + 2]<<" "<<facelist[5 * (facenum+k)+ 3]<<" "
		<<facelist[5 * (facenum+k)+4]<<endl;
		}*/
		facenum += (segnum - 1);

		//edge2facelist
		//new edge on the nmedge
		int edgelen2 = edgelist.size()/2;
		edge2facelist[ edgelen - 1].push_back( facenum - 1);	//the last new edge
		edge2facelist[ edgelen - 1].push_back( edgepos );
		for( int k = 2; k < segnum; k ++)
		{
			edge2facelist[ edgelen - k ].push_back( facenum - k);
			edge2facelist[ edgelen - k].push_back( 1 );
		}
		//old edge in the old face
		int oldedge = ver2edgehash.findInsertSort( thirdverpos, verpos[ 1 ], edgelen2);
		//if( iscw)
		//{
		//	oldedge = ver2edgehash.findInsertSort( thirdverpos, verpos[ 1 ], edgelen2);
		//	//////////////////////////////////////////////////////////////////////////
		//	if( oldedge == 25 )
		//	{
		//		cout<<"cw oldedge is :"<<oldedge<<" "<<thirdverpos<<" "<<verpos[1]<<endl;
		//	}
		//}
		//else
		//{
		//	oldedge = ver2edgehash.findInsertSort( thirdverpos, verpos[ 0 ], edgelen2);
		//	//////////////////////////////////////////////////////////////////////////
		//	if( oldedge == 25 )
		//	{
		//		cout<<"ccw oldedge is :"<<oldedge<<" "<<thirdverpos<<" "<<verpos[0]<<endl;
		//	}
		//}

		if(oldedge == edgelen2 )
		{
			cout<<"trying to find the old edge fails!!!"<<endl;
		}
		else
		{
			for( int k = 0; k < edge2facelist[ oldedge ].size(); k++)
			{
				if( edge2facelist[ oldedge ][ k * 2 ] == curfacelist[ j * 2 ])
				{
					edge2facelist[ oldedge ][ k * 2 ] = facenum - 1;
					break;
				}
			}
		}
		//new edge on the face
		int pos = edgelen + (segnum - 1) * j;
		for( int k = 0; k < segnum - 2; k++)
		{
			edge2facelist[ pos + k ].push_back( facenum - segnum + k + 1);
			if( iscw )
				edge2facelist[ pos + k ].push_back( 0 );
			else
				edge2facelist[ pos + k  ].push_back( 2 );
		}
		//last new edge
		edge2facelist[ pos + segnum - 2].push_back( facenum - 1 );
		if( iscw )
			edge2facelist[ pos + segnum - 2].push_back( (edgepos + 2)%3 );
		else
			edge2facelist[ pos + segnum - 2].push_back( (edgepos + 1)%3 );
		for( int k = segnum - 2; k >= 1; k --)
		{
			edge2facelist[ pos + k ].push_back( facenum - segnum + k);
			if( iscw )
				edge2facelist[ pos + k ].push_back( 2 );
			else
				edge2facelist[ pos + k ].push_back( 0 );
		}
		edge2facelist[ pos ].push_back( curfacelist[ j * 2 ] );
		if( iscw )
			edge2facelist[ pos ].push_back( (edgepos + 1)%3 );
		else
			edge2facelist[ pos ].push_back( (edgepos + 2)%3 );
		//	cout<<"current face is done!"<<endl;
	}
}

void LuMesh::splitNMEdge(vector<double>& verlist,vector<double>& verattrlist,vector<int>& edgelist,vector<intvector>& edge2facelist,
					   vector<int>&nmedgelist, vector<int>& normedgelist,vector<int>&facelist, HashMap& ver2edgehash)
{
	//go through each non-manifold edge, if the length of it is longer than parameter, split it
	int nmedgelen = nmedgelist.size();
	double curlen;	//current non-manifold edge length
	double param;	//the bigger parameter of the two endpoints of current non-manifold edge
	int segnum;		//how many segments it's going to be cut into
	//for debug
	//	for( int i = 0; i < 4; i ++)
	for( int i = 0; i < nmedgelen; i++)
	{
		curlen = MyMath::vectorlen(&verlist[edgelist[nmedgelist[i]*2]*3],	&verlist[edgelist[nmedgelist[i]*2 + 1]*3]);
		param = verattrlist[ edgelist[ nmedgelist[ i ]* 2]];
		if( verattrlist[edgelist[nmedgelist[i] * 2  + 1]] > param )
			param = verattrlist[ edgelist[nmedgelist[i]*2 + 1]];
		//for debug
		//segnum = 4;
		//	cout<<"param"<<param<<"\t"<<curlen<<endl;
		//	segnum = (int)(curlen / param + 0.6);
		segnum = curlen/param;
		if( segnum <= 1)continue;
		//split current edge into segnum segments.
		//add new vertex: refresh verlist, and verattrilist
		//add new edge: refresh edgelist, edge2facelist, nmedgelist, normedgelist, and ctredgelist
		//add new face: refresh facelist
		//ver
		//------------//
		//	cout<<"segnum"<<segnum<<endl;

		splitOneNMEdge(segnum,  i,  verlist, verattrlist,edgelist, edge2facelist,nmedgelist, normedgelist,facelist,  ver2edgehash);
	}
	//!!!old implementation below!!!
	//		int verpos[2] = {edgelist[nmedgelist[i]*2],edgelist[nmedgelist[i]*2 + 1]};
	//	//	cout<<"verpos are:"<<verpos[0]<<"   "<<verpos[1]<<endl;
	//		double newver[3];
	//		int verlen = verlist.size()/3;
	//		for( int j = 1; j < segnum; j++ )
	//		{
	//			MyMath::getPtOnSeg(&verlist[verpos[0]*3], &verlist[verpos[1]*3], j/(double)segnum, newver);
	//			verlist.push_back( newver[ 0 ]) ;
	//			verlist.push_back( newver[ 1 ]);
	//			verlist.push_back( newver[ 2 ]);
	//			verattrlist.push_back( (1 - j/(double)segnum)*verattrlist[verpos[0]] +  j/(double)segnum*verattrlist[verpos[1]] );
	//		}
	//		//------------//
	//	//	cout<<"vertex added!"<<endl;
	//		//new edge on the old non-manifold edge
	//		int edgelen = edgelist.size()/2;
	//		int oldnmedgelen = nmedgelist.size();
	//		nmedgelist.resize(oldnmedgelen + segnum - 1);
	//		for( int j = 1; j < segnum - 1; j ++)
	//		{
	//			edgelist.push_back( verlen + j - 1);
	//			edgelist.push_back( verlen + j);
	//			nmedgelist[ oldnmedgelen + j - 1] = edgelen + j - 1;
	//		}
	//		nmedgelist[ oldnmedgelen + segnum - 2] = edgelen + segnum - 2;
	//		edgelist.push_back( verlen + segnum - 2);	//the last segment
	//		edgelist.push_back( verpos[1] );
	//		edgelist[ nmedgelist[ i ]*2 + 1 ] = verlen;	//second point of the first segment
	//        
	//		int facenum = facelist.size()/5;	//old face number
	//		edgelen = edgelist.size()/2;		//edge number after adding new edges on non-manifold edge
	//		//------------//
	//	//	cout<<"new edge added!"<<endl;
	//		//every face
	//		int curfacenum = edge2facelist[ nmedgelist[ i ]].size()/2;	//incident faces number
	//		intvector& curfacelist = edge2facelist[ nmedgelist[i]];		//incident faces list
	//
	//		int oldnormedgelen = normedgelist.size();
	//		normedgelist.resize( oldnormedgelen + curfacenum*(segnum - 1));
	//		for( int j = 0; j < (segnum - 1)*curfacenum; j++)
	//			normedgelist[ oldnormedgelen + j] = edgelen + j;
	//		edge2facelist.resize(edgelen + curfacenum * (segnum - 1));
	//		//------------//
	////		cout<<"incident face number:"<<curfacenum<<endl;
	//		for( int j = 0; j < curfacenum; j ++)
	//		{
	//		//	cout<<"curface:"<<curfacelist[2*j]<<endl;
	//			int edgepos = curfacelist[ 2*j+1 ];
	//			int thirdverpos = facelist[5 * curfacelist[ 2*j ]+(edgepos + 2)%3];
	//		//	cout<<"current face:"<<facelist[5 * curfacelist[ 2*j ]]<<" "<<facelist[5 * curfacelist[ 2*j ] + 1]<<
	//	//			" "<<facelist[5 * curfacelist[ 2*j ] + 2]<<" "<<facelist[5 * curfacelist[ 2*j ]+ 3]<<" "
	//	//			<<facelist[5 * curfacelist[ 2*j ]+4]<<endl;
	//		//	cout<<"edgepos"<<edgepos<<"\t";
	//		//	cout<<"thirdver:"<<thirdverpos<<"\n";
	//			bool iscw = (verpos[ 0 ] == facelist[ 5 * curfacelist[ 2*j ] + edgepos ]);
	//			int mat[2] = {facelist[ 5 * curfacelist[ 2*j ] + 3], facelist[ 5 * curfacelist[ 2*j ] + 4]};
	//
	//			//add new edge, add new face, refresh edge2facelist
	//			//add new edge
	//            for( int k = 0; k < segnum - 1; k++ )
	//			{
	//				edgelist.push_back( thirdverpos );
	//				edgelist.push_back( verlen + k );
	//			}
	//			//------------//
	//		//	cout<<"new edge added for current face!"<<endl;
	//			//new face
	//			if( iscw )
	//			{
	//				for( int k = 0; k < segnum - 2; k ++ )
	//				{
	//					facelist.push_back( thirdverpos );
	//					facelist.push_back( verlen + k );
	//					facelist.push_back( verlen + k + 1);
	//					facelist.push_back( mat[ 0 ]);
	//					facelist.push_back( mat[ 1 ]);
	//				}
	//				for( int k = 0; k < 5; k ++)
	//					facelist.push_back( facelist[ 5 * curfacelist[ 2*j ] + k ]);
	//				facelist[ ( facenum + segnum - 2 ) * 5 + edgepos ] = verlen + segnum - 2;
	//				facelist[ 5 * curfacelist[ 2*j ] + (edgepos + 1)%3 ] = verlen;
	//			}
	//			else
	//			{
	//				for( int k = 0; k < segnum - 2; k ++)
	//				{
	//					facelist.push_back( thirdverpos );
	//					facelist.push_back( verlen + k + 1);
	//					facelist.push_back( verlen + k );					
	//					facelist.push_back( mat[ 0 ]);
	//					facelist.push_back( mat[ 1 ]);
	//				}
	//				for( int k = 0; k < 5; k ++)
	//					facelist.push_back( facelist[ 5 * curfacelist[ 2*j ] + k ]);
	//				facelist[ ( facenum + segnum - 2 ) * 5 + (edgepos + 1)%3] = verlen + segnum - 2;
	//				facelist[ 5 * curfacelist[ 2*j ] + edgepos ] = verlen;
	//			}
	//			//---------------------//
	//			/*cout<<"faces:"<<facelist[5 * curfacelist[ 2*j ]]<<" "<<facelist[5 * curfacelist[ 2*j ] + 1]<<
	//				" "<<facelist[5 * curfacelist[ 2*j ] + 2]<<" "<<facelist[5 * curfacelist[ 2*j ]+ 3]<<" "
	//				<<facelist[5 * curfacelist[ 2*j ]+4]<<endl;
	//			for( int k = 0; k < segnum - 1; k++ )
	//			{
	//				cout<<"faces:"<<facelist[5 * (facenum+k)]<<" "<<facelist[5 *(facenum+k) + 1]<<
	//					" "<<facelist[5 *(facenum+k) + 2]<<" "<<facelist[5 * (facenum+k)+ 3]<<" "
	//					<<facelist[5 * (facenum+k)+4]<<endl;
	//			}*/
	//			facenum += (segnum - 1);
	//		
	//			//edge2facelist
	//			//new edge on the nmedge
	//			int edgelen2 = edgelist.size()/2;
	//			edge2facelist[ edgelen - 1].push_back( facenum - 1);	//the last new edge
	//			edge2facelist[ edgelen - 1].push_back( edgepos );
	//			for( int k = 2; k < segnum; k ++)
	//			{
	//				edge2facelist[ edgelen - k ].push_back( facenum - k);
	//				edge2facelist[ edgelen - k].push_back( 1 );
	//			}
	//			//old edge in the old face
	//			int oldedge = ver2edgehash.findInsertSort( thirdverpos, verpos[ 1 ], edgelen2);
	//			//if( iscw)
	//			//{
	//			//	oldedge = ver2edgehash.findInsertSort( thirdverpos, verpos[ 1 ], edgelen2);
	//			//	//////////////////////////////////////////////////////////////////////////
	//			//	if( oldedge == 25 )
	//			//	{
	//			//		cout<<"cw oldedge is :"<<oldedge<<" "<<thirdverpos<<" "<<verpos[1]<<endl;
	//			//	}
	//			//}
	//			//else
	//			//{
	//			//	oldedge = ver2edgehash.findInsertSort( thirdverpos, verpos[ 0 ], edgelen2);
	//			//	//////////////////////////////////////////////////////////////////////////
	//			//	if( oldedge == 25 )
	//			//	{
	//			//		cout<<"ccw oldedge is :"<<oldedge<<" "<<thirdverpos<<" "<<verpos[0]<<endl;
	//			//	}
	//			//}
	//
	//			if(oldedge == edgelen2 )
	//			{
	//				cout<<"trying to find the old edge fails!!!"<<endl;
	//			}
	//			else
	//			{
	//				for( int k = 0; k < edge2facelist[ oldedge ].size(); k++)
	//				{
	//					if( edge2facelist[ oldedge ][ k * 2 ] == curfacelist[ j * 2 ])
	//					{
	//						edge2facelist[ oldedge ][ k * 2 ] = facenum - 1;
	//						break;
	//					}
	//				}
	//			}
	//			//new edge on the face
	//			int pos = edgelen + (segnum - 1) * j;
	//			for( int k = 0; k < segnum - 2; k++)
	//			{
	//				edge2facelist[ pos + k ].push_back( facenum - segnum + k + 1);
	//				if( iscw )
	//					edge2facelist[ pos + k ].push_back( 0 );
	//				else
	//					edge2facelist[ pos + k  ].push_back( 2 );
	//			}
	//			//last new edge
	//			edge2facelist[ pos + segnum - 2].push_back( facenum - 1 );
	//			if( iscw )
	//				edge2facelist[ pos + segnum - 2].push_back( (edgepos + 2)%3 );
	//			else
	//				edge2facelist[ pos + segnum - 2].push_back( (edgepos + 1)%3 );
	//			for( int k = segnum - 2; k >= 1; k --)
	//			{
	//				edge2facelist[ pos + k ].push_back( facenum - segnum + k);
	//				if( iscw )
	//					edge2facelist[ pos + k ].push_back( 2 );
	//				else
	//					edge2facelist[ pos + k ].push_back( 0 );
	//			}
	//			edge2facelist[ pos ].push_back( curfacelist[ j * 2 ] );
	//			if( iscw )
	//				edge2facelist[ pos ].push_back( (edgepos + 1)%3 );
	//			else
	//				edge2facelist[ pos ].push_back( (edgepos + 2)%3 );
	//		//	cout<<"current face is done!"<<endl;
	//		}
	//	}
}

void LuMesh::swapEdge(vector<double>&verlist,vector<int>& edgelist,vector<intvector>& edge2facelist,vector<int>&normedgelist, vector<int>&facelist, HashMap& ver2edgehash)
{
	int edgelen = edgelist.size()/2;
	//	HashMap ver2edgehash;
//	int temp;
	/*	for( int i = 0; i < edgelen; i ++)
	{
	if( ver2edgehash.findInsertSort( edgelist[ i * 2 ], edgelist[ i * 2 + 1], i ) != i)
	cout<<"Error occurs when building hash map for edges!"<<endl;
	}*/
	bool again = true;
	int count = 0;
	edgelen = normedgelist.size();

	//////////////////////////////////////////////////////////////////////////
	/*	FILE* fout = fopen("out.txt","w");
	for( int i = 0; i < edgelen; i ++)
	{
	fprintf( fout, "%d\t%d\t%d\n",i,edgelist[ 2*i ], edgelist[ 2*i+1 ]);
	//all the faces
	for( int j = 0; j < edge2facelist[ i ].size()/2; j++)
	{
	int facei = edge2facelist[ i ][ 2*j ];
	fprintf(fout, "%d\t%d\t%d\t%d\n", facelist[5*facei], facelist[5*facei+1],facelist[5*facei+2],
	edge2facelist[i][j*2+1]);
	}
	}
	fclose(fout);*/
	//////////////////////////////////////////////////////////////////////////
	while( again && count < 1024)
	{
		again = false;
		count ++;

		/////////////////////////check/////////////////////////////////////////////////
		/*for( int i = 0;  i< edgelen ; i ++)
		{
		for( int j = 0; j < edge2facelist[i].size()/2; j++)
		{
		int facei = edge2facelist[i][ j*2 ];
		int pos = edge2facelist[i][j*2+1];
		bool issame = false;
		issame = ((edgelist[ 2*i ] == facelist[ 5*facei + pos ]) && (edgelist[2*i+1]==facelist[5*facei+(pos+1)%3]));
		if(issame)continue;
		issame = ((edgelist[ 2*i ] == facelist[ 5*facei +(pos+1)%3 ]) && (edgelist[2*i+1]==facelist[5*facei+ pos]));
		if( issame )continue;
		cout<<"edge :"<<i<<" "<<edgelist[ 2*i ]<<edgelist[ 2*i+1 ]<<"pos:"<<pos<<endl;
		cout<<"face:"<<facelist[5*facei]<<" "<<facelist[5*facei+1]<<" "<<facelist[5*facei+2]<<endl;
		}
		}*/
		//go through each interior edge, swap it if needed
		//for( int i = 0; i < edgelen; i ++)
		//	cout<<edgelen;
		for( int i = 0; i < edgelen; i ++)
		{			
			intvector& facepair = edge2facelist[ normedgelist[ i ]];
			if( facepair.size() != 4)
			{
				/*-------*/
				if( facepair.size() > 4)
					cout<<"incident faces number is over 4, error!"<<endl;
				else 
					cout<<"incident faces number is smaller than 4, error!"<<endl;
				continue;	
			}

			int verpos[2] = {edgelist[ normedgelist[ i ]*2], edgelist[ normedgelist[ i ]*2+1]};
			int thirdver[ 2 ] = { facelist[ 5 * facepair[ 0 ] + (facepair[ 1 ] + 2)%3 ], 
				facelist[ 5 * facepair[ 2 ] + (facepair[ 3 ] + 2)%3 ]};
			//	cout<<verpos[0]<<" "<<verpos[1]<<" "<<thirdver[0]<<" "<<thirdver[1]<<endl;
			//	cout<<"face:"<<facepair[0]<<"edge:"<<facepair[1]<<" "<<facelist[ 5 * facepair[ 0 ]]<<" "<<facelist[ 5 * facepair[ 0 ]+1]<<" "<<facelist[ 5 * facepair[ 0 ]+2]<<endl;
			//	cout<<"face:"<<facepair[2]<<"edge:"<<facepair[3]<<" "<<facelist[ 5 * facepair[ 2 ]]<<" "<<facelist[ 5 * facepair[ 2 ]+1]<<" "<<facelist[ 5 * facepair[ 2 ]+2]<<endl;

			//compute the angle of the third angles
			double cosab[2] ;
			cosab[ 0 ] = MyMath::getCosOfAngle(&verlist[3*verpos[0]], &verlist[3*verpos[1]], &verlist[ 3 * thirdver[0]]);
			cosab[ 1 ] = MyMath::getCosOfAngle(&verlist[3*verpos[0]], &verlist[3*verpos[1]], &verlist[ 3 * thirdver[1]]);
			//////////////////////////////////////////////////////////////////////////
			/*if( cosab[0]<0)
			cout<<"cos1:"<<cosab[0];
			if( cosab[1]<0)
			cout<<"\tcos2:"<<cosab[1]<<endl;*/
			//////////////////////////////////////////////////////////////////////////
			/*if( cosab[0] >= 1 || cosab[ 0 ] <= -1)
			{
			cout<<"pt1:"<<verlist[3*verpos[0]]<<", "<<verlist[3*verpos[0]+1]<<", "<<verlist[3*verpos[0]+2]<<endl;
			cout<<"pt2:"<<verlist[3*verpos[1]]<<", "<<verlist[3*verpos[1]+1]<<", "<<verlist[3*verpos[1]+2]<<endl;
			cout<<"pt0:"<<verlist[ 3 * thirdver[0]]<<", "<<verlist[3 * thirdver[0]+1]<<", "<<verlist[ 3 * thirdver[0]+2]<<endl;
			cout<<MyMath::getCosOfAngle(&verlist[3*verpos[0]], &verlist[3*verpos[1]], &verlist[ 3 * thirdver[0]])<<endl;
			}
			if( cosab[1] >= 1 && cosab[ 1 ] <= -1)
			{
			cout<<"pt1:"<<verlist[3*verpos[0]]<<", "<<verlist[3*verpos[0]+1]<<", "<<verlist[3*verpos[0]+2]<<endl;
			cout<<"pt2:"<<verlist[3*verpos[1]]<<", "<<verlist[3*verpos[1]+1]<<", "<<verlist[3*verpos[1]+2]<<endl;
			cout<<"pt0:"<<verlist[ 3 * thirdver[1]]<<", "<<verlist[3 * thirdver[1]+1]<<", "<<verlist[ 3 * thirdver[1]+2]<<endl;
			cout<<MyMath::getCosOfAngle(&verlist[3*verpos[0]], &verlist[3*verpos[1]], &verlist[ 3 * thirdver[1]])<<endl;
			}*/
			double sinab[2];
			sinab[ 0 ] = sqrt(1 - cosab[0]*cosab[ 0 ]);
			sinab[  1 ] = sqrt(1 -cosab[ 1 ]*cosab[ 1 ]);
			//	cout<<cosab[0]*cosab[ 0 ]<<"  "<<sinab[0]<<cosab[ 1 ]*cosab[ 1 ]<<" "<<sinab[  1 ]<<endl;
			//cout<<(cosab[ 0 ]* sinab[ 1 ] + cosab[ 1 ] * sinab[ 0 ])<<" ";
			if( (cosab[ 0 ]* sinab[ 1 ] + cosab[ 1 ] * sinab[ 0 ]) < -0.000001)	//>180
			{
				//see if the new edge going to add is in the mesh or not
				if(ver2edgehash.findInsertSort(thirdver[0], thirdver[1], -1) != -1)
					continue;
				again = true;
				//flip edge,change two faces, change the affected edge's edge2facelist
				//edge list and hash
				edgelist[ normedgelist[i] * 2 ] = thirdver[ 0 ];
				edgelist[ normedgelist[i] * 2 + 1] = thirdver[ 1 ];
				//-----------------//
				//	cout<<endl;
				//	cout<<"the two faces are:"<<endl;
				/*if( i == 85)
				for( int j = 0; j < 2; j ++)
				{
				cout<<facelist[ facepair[2*j]*5]<<"  "<<facelist[ facepair[2*j]*5+1]<<"  "<<facelist[ facepair[2*j]*5+2]<<"  "<<
				facelist[ facepair[2*j]*5+3]<<"  "<<facelist[ facepair[2*j]*5+4]<<endl;
				}
				cout<<"thirdver:"<<thirdver[0]<<"\t"<<thirdver[1]<<endl;*/
				//replace the old edge if it exists
				ver2edgehash.findInsertSortReplace(verpos[0], verpos[1], -1);	//edge is not in the mesh now
				ver2edgehash.findInsertSortReplace( thirdver[0], thirdver[1], normedgelist[ i ]);
				//facelist
				int edgei;
				//////////////////////////////////////////////////////////////////////////
				//	cout<<"new edge:"<<thirdver[0]<<"  "<<thirdver[1]<<endl;
				for( int j = 0; j < 2; j ++)
				{
					//for facepair[ 2*j ], change its verpos[ j ] to thirdver[ 1 - j ]
					int facei = facepair[ j * 2 ];
					int edgepos = facepair[ j*2 + 1 ];
					if( facelist[ 5*facei + facepair[ j*2 + 1 ] ] == verpos[ j ])
					{
						//cout<<"type1"<<endl;
						facelist[ 5*facei + facepair[ j*2 + 1 ] ] = thirdver[ 1 - j ];
						/*if( i == 85)
						cout<<"facepair[ j*2 + 1 ]"<<facepair[ j*2 + 1 ]<<endl;*/
						/*if( i == 85)
						cout<<"trying to find old edge:"<< thirdver[1-j]<<" "<<facelist[ 5*facei + (facepair[ j*2 + 1 ]+1)%3]<<endl;*/
						edgei = ver2edgehash.findInsertSort( thirdver[1-j], facelist[ 5*facei + (facepair[ j*2 + 1 ]+1)%3], -1);
						//////////////////////////////////////////////////////////////////////////
						//	cout<<"vers:"<<thirdver[1-j]<<" "<<facelist[ 5*facei + (facepair[ j*2 + 1 ]+1)%3]<<" edge:"<<edgei<<"edgepos:"<<edgepos<<endl;
						//	cout<<"current triangle:"<<facelist[5*facei]<<" "<<facelist[5*facei+1]<<" "<<facelist[5*facei+2]<<endl;
						facepair[ j * 2 + 1 ] = (facepair[ j * 2 + 1] + 2)%3;
						////////////////////////////////check//////////////////////////////////////////
						/*if(!(((facelist[5*facei + facepair[j*2+1]] == thirdver[0])&&(facelist[5*facei + (facepair[j*2+1]+1)%3] == thirdver[1]))
						||((facelist[5*facei + facepair[j*2+1]] == thirdver[1])&&(facelist[5*facei + (facepair[j*2+1]+1)%3] == thirdver[0]))))
						{							
						cout<<"thirdver:"<<thirdver[0]<<"\t"<<thirdver[1]<<endl;
						cout<<"face:"<<facelist[ 5*facei ]<<" "<<facelist[5*facei+1]<<" "<<facelist[5*facei+2]<<" edge:"<<facepair[2*j+1]<<endl;

						cout<<"new edge error! count:"<<count<<"edge "<<i<<"type1"<<endl;
						}*/
						//	cout<<"face:"<<facelist[5*facei]<<" "<<facelist[facei*5+1]<<"  "<<facelist[facei*5+2]<<" pos:"<<facepair[j*2+1]<<endl;
						if( edgei == -1)
						{
							cout<<"error occurs while finding the edge in the reconfigured face!"<<endl;
						}
						else
						{
							intvector& tfacelist = edge2facelist[ edgei ];	
							for( int k = 0; k < tfacelist.size()/2; k++)
							{
								if( tfacelist[ k * 2 ] == facepair[ (1-j)*2 ])
								{
									tfacelist[ k * 2 ] = facei;
									tfacelist[ k*2 + 1] = edgepos ;
									//////////////////////////////////////////////////////////////////////////
									//	 cout<<thirdver[1-j]<<" "<<
									/*if( i == 85)
									cout<<"found!"<<tfacelist[ k * 2 ]<<" "<<tfacelist[ k * 2 + 1]<<endl;*/
									break;
								}
							}
						}
					}
					else
					{
						//	cout<<"type2"<<endl;
						facelist[ 5*facei + (facepair[ j*2 + 1 ]+1)%3 ] = thirdver[ 1 - j ];
						/*if( i == 85)
						cout<<"trying to find old edge:"<<" "<< thirdver[1-j]<<facelist[ 5*facei + facepair[ j*2 + 1 ]]<<endl;*/
						edgei = ver2edgehash.findInsertSort( thirdver[1-j], facelist[ 5*facei + facepair[ j*2 + 1 ]], -1);
						//	cout<<"vers:"<<thirdver[1-j]<<" "<<facelist[ 5*facei + (facepair[ j*2 + 1 ]+1)%3]<<" edge:"<<edgei<<"edgepos:"<<edgepos<<endl;
						//	cout<<"current triangle:"<<facelist[5*facei]<<" "<<facelist[5*facei+1]<<" "<<facelist[5*facei+2]<<endl;
						facepair[ j*2 + 1] = (facepair[ j*2 + 1 ]+1)%3;
						//	cout<<"face:"<<facelist[5*facei]<<" "<<facelist[facei*5+1]<<"  "<<facelist[facei*5+2]<<" pos:"<<facepair[j*2+1]<<endl;
						////////////////////////////////check//////////////////////////////////////////
						/*if(!( ((facelist[5*facei + facepair[j*2+1]] == thirdver[0])&&(facelist[5*facei + (facepair[j*2+1]+1)%3] == thirdver[1]))
						||((facelist[5*facei + facepair[j*2+1]] == thirdver[1])&&(facelist[5*facei + (facepair[j*2+1]+1)%3] == thirdver[0]))))
						{
						cout<<"thirdver:"<<thirdver[0]<<"\t"<<thirdver[1]<<endl;
						cout<<"face:"<<facelist[ 5*facei ]<<" "<<facelist[5*facei+1]<<" "<<facelist[5*facei+2]<<" edge:"<<facepair[2*j+1]<<endl;

						cout<<"new edge error! count:"<<count<<"edge "<<i<<"type2"<<endl;
						}*/
						if( edgei == -1)
						{
							cout<<"error occurs while finding the edge in the reconfigured face!"<<endl;
						}
						else
						{
							intvector& tfacelist = edge2facelist[ edgei ];	
							for( int k = 0; k < tfacelist.size()/2; k++)
							{
								if( tfacelist[ k * 2 ] == facepair[ (1-j)*2 ])
								{
									tfacelist[ k * 2 ] = facei;
									tfacelist[ k*2 + 1] = edgepos ;
									/*if( i == 85)
									cout<<"found!"<<tfacelist[ k * 2 ]<<" "<<tfacelist[ k * 2 + 1]<<endl;*/
									break;
								}
							}
						}
					}
				}
				//	cout<<"the two faces are:"<<endl;
				/*if( i == 85)
				for( int j = 0; j < 2; j ++)
				{
				cout<<facelist[ facepair[2*j]*5]<<"  "<<facelist[ facepair[2*j]*5+1]<<"  "<<facelist[ facepair[2*j]*5+2]<<"  "<<
				facelist[ facepair[2*j]*5+3]<<"  "<<facelist[ facepair[2*j]*5+4]<<endl;
				}*/
				//	cout<<"one edge is flipped!"<<endl;
				//		break;
			}			
		}
		//	cout<<"again is "<<again<<endl;
		//	if(count == 2)
		//		break;
		//	break;
	}
	cout<<"called "<<count<<" times!";
}
bool LuMesh::splitTriangle(vector<double>&verlist, vector<double>&verattrilist, HashMap& ver2edgehash2,
						 vector<int>&edgelist, vector<int>&normedgelist, vector<intvector>&edge2facelist,
						 vector<int>&facelist, double alpha)
{
	bool splitExist = false;
	//	bool again = true;
	//	int count = 0;
	int edgenum = edgelist.size()/2;
	//	while( again && count < MAXTIMES)
	//	{
	//		again = false;
	//		count ++;
	int oldfacenum , facenum;
	oldfacenum = facenum = facelist.size()/5;
	for(int curtri = 0; curtri < oldfacenum ; curtri++)
	{
		//centroid
		double centerpos[3] = {0,0,0};
		for( int i = 0; i < 3; i ++)
		{
			for( int j = 0; j < 3; j++)
			{
				centerpos[ i ] += verlist[ 3 * facelist[ curtri * 5 + j ] + i];
			}
			centerpos[ i ]/=3;
		}
		//distance
		double distance[3] = {0,0,0};
		bool split = true;
		double centerattr = 0;
		for( int i = 0; i < 3; i ++)
		{
			distance[ i ] = MyMath::vectorlen( centerpos, &verlist[ 3 * facelist[ curtri*5 + i ]]);
			if( alpha * distance[ i ] < verattrilist[ facelist[curtri*5+i]] )
			{
				split = false;
				break;
			}
			centerattr += verattrilist[ facelist[curtri*5+i]];
		}
		if( !split )continue;

		//////////////////////////////////////////////////////////////////////////
		/*cout<<"faces vertex:";
		for( int i = 0; i < 5; i ++)
		{
		cout<<facelist[5*curtri + i]<<" "; 
		}
		cout<<endl;
		cout<<"edges : ";
		for( int i = 0; i< 3; i ++)
		{
		cout<<ver2edgehash2.findInsertSort(facelist[5*curtri + i], facelist[5*curtri + (i+1)%3], -1)<<"  ";
		}
		cout<<endl;*/
		//////////////////////////////////////////////////////////////////////////
		splitExist = true;
		//split
		//add new vertex
		verlist.push_back( centerpos[0] );
		verlist.push_back( centerpos[1]);
		verlist.push_back( centerpos[2]);
		centerattr/=3;
		verattrilist.push_back( centerattr );
		int verlen = verlist.size()/3;
		//new edge
		for( int i = 0; i< 3; i ++)
		{
			edgelist.push_back(verlen - 1);
			edgelist.push_back(facelist[curtri*5 + i]);
			int tedge = 	ver2edgehash2.findInsertSort( verlen - 1, facelist[ curtri*5 + i ], edgenum + i);
			if( tedge != edgenum + i)
			{
				cout<<"tedge is :"<<tedge<<"should be:"<<edgenum+i<<endl;
			}
			normedgelist.push_back( edgenum + i);
		}			
		edgenum += 3;			

		//face
		//add new
		for( int i = 1; i< 3; i ++)
		{
			facelist.push_back(facelist[curtri*5+ i ]);
			facelist.push_back( facelist[ curtri*5 + (i+1)%3]);
			facelist.push_back( verlen - 1);
			facelist.push_back( facelist[ curtri * 5 + 3]);
			facelist.push_back( facelist[ curtri*5 + 4]);				
		}
		facelist[ 5* curtri + 2] = verlen - 1;
		facenum += 2;

		//edge2facelist
		int ver1;
		int ver2;
		int edgei;
		//for the two old face edges
		for( int i = 1; i <=2; i++)
		{
			ver1 = facelist[ 5* (facenum + i - 3)];
			ver2 = facelist[ 5* (facenum + i - 3) + 1];
			edgei = ver2edgehash2.findInsertSort(ver1, ver2, -1);
			if( edgei == -1)
			{
				cout<<"Error occurs when trying to find the old edge!"<<endl;
			}
			else
			{
				//////////////////////////////////////////////////////////////////////////
				//	cout<<"ver1:"<<ver1<<" ver2:"<<ver2<<" edgei:"<<edgei<<endl;
				//////////////////////////////////////////////////////////////////////////
				vector<int>& tfacelist = edge2facelist[ edgei ];
				for( int j = 0; j < tfacelist.size()/2; j ++)
				{
					if( tfacelist[ j * 2 ] == curtri )
					{
						//////////////////////////////////////////////////////////////////////////
						//			cout<<"curtri is found in the list!"<<endl;
						//////////////////////////////////////////////////////////////////////////
						tfacelist[ j * 2 ] = facenum + i - 3;
						tfacelist[ j*2 + 1] = 0;
						break;
					}
				}
			}
		}
		//for the three new face edges
		edge2facelist.resize( edgenum );
		intvector* tfacelist = &edge2facelist[ edgenum - 3];
		tfacelist->push_back( curtri );
		tfacelist->push_back( 2 );
		tfacelist->push_back( facenum - 1);
		tfacelist->push_back(1);
		tfacelist = &edge2facelist[ edgenum - 2];
		tfacelist->push_back( curtri );
		tfacelist->push_back( 1 );
		tfacelist->push_back( facenum - 2 );
		tfacelist->push_back( 2 );
		tfacelist = &edge2facelist[ edgenum - 1];
		tfacelist->push_back( facenum - 2);
		tfacelist->push_back( 1 ) ;
		tfacelist->push_back( facenum - 1);
		tfacelist->push_back( 2 );

		//////////////////////////////////////////////////////////////////////////
		/*	cout<<"edgenum:"<<edgenum<<endl;
		for( int i = 0; i < normedgelist.size(); i ++)
		{
		//cout<<edge2facelist[ normedgelist[i] ].size()<<" " ;
		if(edge2facelist[ normedgelist[i] ].size()!=4)
		cout<<"edge "<<normedgelist[i]<<edge2facelist[ normedgelist[i] ].size()<<" ";
		}
		cout<<endl;*/
		//////////////////////////////////////////////////////////////////////////
		//	break;	////
		//propagate the swap

	}

	//	}
	return splitExist;
}
