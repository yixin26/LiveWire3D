/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 **************************************************************/

#include "PathUtils.h"
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <float.h>
#include <map>
#include <time.h>

using namespace std;

PathUtils::PathUtils()
{
	init();
}
PathUtils::~PathUtils()
{
}

clock_t myclock()
{
	unsigned __int64 current_tics;
	FILETIME ct;

	GetSystemTimeAsFileTime( &ct );

	current_tics = (unsigned __int64)ct.dwLowDateTime +
		(((unsigned __int64)ct.dwHighDateTime) << 32);

	/* calculate the elapsed number of 100 nanosecond units */
	//current_tics -= start_tics;

	/* return number of elapsed milliseconds */
	//cout<<endl<<"time called"<<current_tics / 10<<endl;
	return (clock_t)(current_tics / 10);

}

double PathUtils::pointCloseness = 1e-8;

void PathUtils::init()
{
	//init and reset 
	points.clear();	faces.clear(); myMesh = NULL;

	triMesh = NULL; strcpy(curvatureScale, "0.");	strcpy(curvatureSmooth, "0.1");

	m_anisAlgo = OUR_MIN;
	euclideanGeodesy = ourGeodesy1 = ourGeodesy2 = pottmannGeodesy = kovacsGeodesy1 = kovacsGeodesy2 = campenGeodesy1 = campenGeodesy2 = NULL;
	euclideanVisual = ourVisual1 = ourVisual2 = pottmannVisual = kovacsVisual1 = kovacsVisual2 = campenVisual1 = campenVisual2 = NULL;

	m_curves.clear(); m_curveVE.clear(); m_curveVs.clear();
	m_Path.clear(); m_PathIndex.clear(); m_VE.clear();

	showCurveNetwork = true; nodeSize = 6; lineWidth = 3;
	showMesh = true;	meshType = 1;

	showFieldDir1 = showFieldDir2 = showFieldEllipse = showField8Shape = false;
	showCurvatureMin = showCurvatureMax = showAnistropy = false;

	//shortest path on mesh and subgraph
	//user selection
	m_snapItem.pos = Vec3(0, 0, 0); 
	m_snapItem.dis = 1000;
	isDoEuclideanGeodesy = false;
	isDoOursGeodesy1=true;	isDoOursGeodesy2=true;
	//isDoOursGeodesy1 = false;	isDoOursGeodesy2 = false;
	isDoPottmannsGeodesy = false;
	isDoKovacsGeodesy1 = false;	isDoKovacsGeodesy2 = false;
	isDoCampensGeodesy1 = false;	isDoCampensGeodesy2 = false;

	//editing on curve
	selectedCurve.first = selectedCurve.second = 0;

	meshColor = Vec3(.5, 0.539, 0.527);
	geoColor = Vec3(1, 0, 0);
	ourColor1 = Vec3(0, 0, 1);
	ourColor2 = Vec3(0, 0, 1);
	pottColor = Vec3(0, 1, 0);
	kovColor1 = Vec3(0, 1, 0);
	kovColor2 = Vec3(0, 1, 0);
	camColor1 = Vec3(1, 0, 1);
	camColor2 = Vec3(1, 0, 1);
	curveColor = Vec3(0.188235, 0.188235, 0.188235);

	meshSegments = NULL;
	meshPatches.clear();
	showPatches = false;

	m_interaction = NULL;

	edgeWidth = 1.0;
	edgeColor = Vec3(0.5, 0.5, 0.5);
}

void PathUtils::computeCurvature()
{
	cout << "\n*compute curvature*" << endl;
	TriMesh *tmesh= new TriMesh(*triMesh);
	unsigned vNum=tmesh->vertices.size();
	tmesh->colors.resize(triMesh->vertices.size());
	colorbycurv(tmesh,curvatureScale,curvatureSmooth);

	unsigned i=0; double maxAnis=0; 
	for(auto v_it=myMesh->getVertices().begin(); v_it!=myMesh->getVertices().end();v_it++)
	{
		v_it->normal()=Vec3(tmesh->normals[i][0],tmesh->normals[i][1],tmesh->normals[i][2]);
		v_it->color()=Vec3(tmesh->colors[i][0],tmesh->colors[i][1],tmesh->colors[i][2]);
		v_it->direction(0)=Vec3(tmesh->pdir1[i][0],tmesh->pdir1[i][1],tmesh->pdir1[i][2]);
		v_it->direction(1)=Vec3(tmesh->pdir2[i][0],tmesh->pdir2[i][1],tmesh->pdir2[i][2]);
		v_it->magnitude(0)=fabs(tmesh->curv1[i]);
		v_it->magnitude(1)=fabs(tmesh->curv2[i]);

		double& k1=v_it->magnitude(0);
		double& k2=v_it->magnitude(1);

		double diff = fabs(k1-k2);

		if(maxAnis<diff) 
		{
			maxAnis=diff;
		}
		i++;
	}
	AnisGeodesic::maxAnis = maxAnis;
	delete tmesh;
}
void PathUtils::computeAnisField()
{
	cout << "\n*compute anisotropic filed*" << endl;
	if(isDoEuclideanGeodesy)	{ euclideanGeodesy->anisInit_Euclidean();euclideanGeodesy->meshEmbedding_Euclidean();}
	if(isDoOursGeodesy1)		{ ourGeodesy1->anisInit_OurMax(ourGeodesy1->metricVariation);ourGeodesy1->meshEmbedding_Tensor();}
	if(isDoOursGeodesy2)		{ ourGeodesy2->anisInit_OurMin(ourGeodesy2->metricVariation);ourGeodesy2->meshEmbedding_Tensor();}
	if(isDoPottmannsGeodesy)	{ pottmannGeodesy->anisInit_Pott();pottmannGeodesy->meshEmbedding_Pott();}
	if(isDoKovacsGeodesy1)		{ kovacsGeodesy1->anisInit_KovsMax();kovacsGeodesy1->meshEmbedding_Tensor();}
	if(isDoKovacsGeodesy2)		{ kovacsGeodesy2->anisInit_KovsMin();kovacsGeodesy2->meshEmbedding_Tensor();}
	if(isDoCampensGeodesy1)		{ campenGeodesy1->anisInit_CamMax();campenGeodesy1->meshEmbedding_Tensor();}
	if(isDoCampensGeodesy2)		{ campenGeodesy2->anisInit_CamMin();campenGeodesy2->meshEmbedding_Tensor();}
}
void PathUtils::anisGeodesyInit()
{
	cout << "\n*initializing geodesy*" << endl;
	if(isDoEuclideanGeodesy)	euclideanGeodesy->geodesyInit();
	if(isDoOursGeodesy1)		ourGeodesy1->geodesyInit();
	if(isDoOursGeodesy2)		ourGeodesy2->geodesyInit();
	if(isDoPottmannsGeodesy)	pottmannGeodesy->geodesyInit();
	if(isDoKovacsGeodesy1)		kovacsGeodesy1->geodesyInit();
	if(isDoKovacsGeodesy2)		kovacsGeodesy2->geodesyInit();
	if(isDoCampensGeodesy1)		campenGeodesy1->geodesyInit();
	if(isDoCampensGeodesy2)		campenGeodesy2->geodesyInit();
}
void PathUtils::anisGeodesyPropation(SnapItem& s)
{
	if(isDoEuclideanGeodesy)	euclideanGeodesy->geodesyPropagation(s.type,s.id,s.pos);
	if(isDoOursGeodesy1)		ourGeodesy1->geodesyPropagation(s.type,s.id,s.pos);
	if(isDoOursGeodesy2)		ourGeodesy2->geodesyPropagation(s.type,s.id,s.pos);
	if(isDoPottmannsGeodesy)	pottmannGeodesy->geodesyPropagation(s.type,s.id,s.pos);
	if(isDoKovacsGeodesy1)		kovacsGeodesy1->geodesyPropagation(s.type,s.id,s.pos);
	if(isDoKovacsGeodesy2)		kovacsGeodesy2->geodesyPropagation(s.type,s.id,s.pos);
	if(isDoCampensGeodesy1)		campenGeodesy1->geodesyPropagation(s.type,s.id,s.pos);
	if(isDoCampensGeodesy2)		campenGeodesy2->geodesyPropagation(s.type,s.id,s.pos);
}
void PathUtils::anisGeodesyTracing(SnapItem& t)
{
	if(isDoEuclideanGeodesy)	euclideanGeodesy->geodesyTraceBack(t.type,t.id,t.pos);
	if(isDoOursGeodesy1)		ourGeodesy1->geodesyTraceBack(t.type,t.id,t.pos);
	if(isDoOursGeodesy2)		ourGeodesy2->geodesyTraceBack(t.type,t.id,t.pos);
	if(isDoPottmannsGeodesy)	pottmannGeodesy->geodesyTraceBack(t.type,t.id,t.pos);
	if(isDoKovacsGeodesy1)		kovacsGeodesy1->geodesyTraceBack(t.type,t.id,t.pos);
	if(isDoKovacsGeodesy2)		kovacsGeodesy2->geodesyTraceBack(t.type,t.id,t.pos);
	if(isDoCampensGeodesy1)		campenGeodesy1->geodesyTraceBack(t.type,t.id,t.pos);
	if(isDoCampensGeodesy2)		campenGeodesy2->geodesyTraceBack(t.type,t.id,t.pos);
}
void PathUtils::srcToAllGeodesy()
{
	AnisVisual* iVisual=NULL; 
	AnisGeodesic* iGeodesic=NULL; 
	bool isInit = false;
	if(euclideanVisual->showMesh)	{iVisual = euclideanVisual; iGeodesic=euclideanGeodesy; isInit=isDoEuclideanGeodesy;}
	if(ourVisual1->showMesh)		{iVisual = ourVisual1; iGeodesic=ourGeodesy1; isInit=isDoOursGeodesy1;}
	if(ourVisual2->showMesh)		{iVisual = ourVisual2; iGeodesic=ourGeodesy2; isInit=isDoOursGeodesy2;}
	if(pottmannVisual->showMesh)	{iVisual = pottmannVisual; iGeodesic=pottmannGeodesy; isInit=isDoPottmannsGeodesy;}
	if(kovacsVisual1->showMesh)		{iVisual = kovacsVisual1; iGeodesic=kovacsGeodesy1; isInit=isDoKovacsGeodesy1;}
	if(kovacsVisual2->showMesh)		{iVisual = kovacsVisual2; iGeodesic=kovacsGeodesy2; isInit=isDoKovacsGeodesy2;}
	if(campenVisual1->showMesh)		{iVisual = campenVisual1; iGeodesic=campenGeodesy1; isInit=isDoCampensGeodesy1;}
	if(campenVisual2->showMesh)		{iVisual = campenVisual2; iGeodesic=campenGeodesy2; isInit=isDoCampensGeodesy2;}

	if(!isInit) return;

	iVisual->srcToAllPaths.clear(); 
	iVisual->srcToAllPathsLength.clear();
	for(unsigned j=0;j<points.size();j++)
	{
		iGeodesic->geodesyTraceBack(AnisGeodesic::VERTEX, j);
		iVisual->srcToAllPaths.push_back(iGeodesic->geoPath);
		iVisual->srcToAllPathsLength.push_back(iGeodesic->geoPathCost);
	}

	iVisual->genDistanceFieldColor();

	iVisual->genIsoLine();
	iVisual->genDisFieldRank();
}
void PathUtils::updateCurvature()
{
	computeCurvature();
	updateMetrics();
}
void PathUtils::updateMetrics()
{
	computeAnisField();
	anisGeodesyInit();
	if(!m_PathIndex.empty())
	{
		anisGeodesyPropation(m_PathIndex.back());
	}
}
void PathUtils::globalParameter()
{
	if(isDoEuclideanGeodesy)	euclideanGeodesy->computeParameter();
// 	if(isDoOursGeodesy1)		ourGeodesy1->computeParameter();
// 	if(isDoOursGeodesy2)		ourGeodesy2->computeParameter();
// 	if(isDoPottmannsGeodesy)	pottmannGeodesy->computeParameter();
// 	if(isDoKovacsGeodesy1)		kovacsGeodesy1->computeParameter();
// 	if(isDoKovacsGeodesy2)		kovacsGeodesy2->computeParameter();
// 	if(isDoCampensGeodesy1)		campenGeodesy1->computeParameter();
// 	if(isDoCampensGeodesy2)		campenGeodesy2->computeParameter();
	if (isDoOursGeodesy1 || isDoOursGeodesy2)	
	{
		ourGeodesy1->computeParameter();
		ourGeodesy2->computeParameter();
	}
	if (isDoPottmannsGeodesy)	pottmannGeodesy->computeParameter();
	if (isDoKovacsGeodesy1 || isDoKovacsGeodesy2)	
	{
		kovacsGeodesy1->computeParameter();
		kovacsGeodesy2->computeParameter();
	}
	if (isDoCampensGeodesy1 || isDoCampensGeodesy2)	
	{
		campenGeodesy1->computeParameter();
		campenGeodesy2->computeParameter();
	}
}

void PathUtils::modifyField()
{
	std::vector<unsigned> edgeSequence;//edge sequence along unfolded triangle strip of scribbles;
	std::vector<MyInteraction::SketchFace> ToAdd; //sequence of triangles;
	std::vector<unsigned> sketchFaces;
	std::vector<Vec3> sketchCurves;
	m_interaction->processSketchForSeg(sketchFaces, sketchCurves);
	m_interaction->userSketchSamplingFromSketchFaces(sketchFaces, sketchCurves, ToAdd, edgeSequence);//find sketch curve from unordered sketch faces;
	std::vector<std::pair<unsigned, Vec3> > newDirection(ToAdd.size());
	unsigned sfid = 0;
	for (auto sf = ToAdd.begin(); sf != ToAdd.end(); sf++, sfid++)
	{
		auto& newD = newDirection[sfid];
		newD.first = sf->fid;
		newD.second = sf->ep1 - sf->ep2;
	}

	if (m_anisAlgo == OUR_MAX)			{ ourGeodesy1->modifyAnisField(newDirection); }
	else if (m_anisAlgo == OUR_MIN)		{ ourGeodesy2->modifyAnisField(newDirection); }
	else if (m_anisAlgo == KOVACS_MAX)	{ kovacsGeodesy1->modifyAnisField(newDirection); }
	else if (m_anisAlgo == KOVACS_MIN)	{ kovacsGeodesy2->modifyAnisField(newDirection); }
	else if (m_anisAlgo == CAMPEN_MAX)	{ campenGeodesy1->modifyAnisField(newDirection); }
	else if (m_anisAlgo == CAMPEN_MIN)	{ campenGeodesy2->modifyAnisField(newDirection); }
}
void PathUtils::resetField(MetricType algoType)
{
	if(algoType==EUCLIDEAN || algoType==POTTMANN) return;
// 	AnisGeodesic* tPathUtils;
// 	if(algoType==OUR_MAX)			{tPathUtils=ourGeodesy1;}
// 	else if(algoType==OUR_MIN)		{tPathUtils=ourGeodesy2;}
// 	else if(algoType==KOVACS_MAX)	{tPathUtils=kovacsGeodesy1;}
// 	else if(algoType==KOVACS_MIN)	{tPathUtils=kovacsGeodesy2;}
// 	else if(algoType==CAMPEN_MAX)	{tPathUtils=campenGeodesy1;}
// 	else if(algoType==CAMPEN_MIN)	{tPathUtils=campenGeodesy2;}

	updateMetrics();
}


//interaction;
int PathUtils::pickMeshVertex(int mouseX,int mouseY,unsigned seletionMethod,bool isStore)
{
	bool hitted;
	if(seletionMethod==0) 
		hitted = m_interaction->glSelection1(mouseX, mouseY, m_interaction->hit.hittedTriangle);
	else if(seletionMethod==1) 
		hitted = m_interaction->glSelection2(mouseX, mouseY, m_interaction->hit.hittedTriangle);

	bool isNewNode=false;
	if (hitted)
	{ //a hit on a triangle, then looks for the closest node of the triangle;
		unsigned nodeID = m_interaction->getClosestVertex(mouseX, mouseY, m_interaction->hit.hittedTriangle);
		if (m_interaction->hit.hittedNode != nodeID)
		{
			m_interaction->hit.hittedNode = nodeID;
			isNewNode=true;
		}
	}

	SnapItem hitPoint; ////if the hit is very close to a drawn curve, then snaps to it.
	auto &point = points[m_interaction->hit.hittedNode];
	m_snapItem.dis=DBL_MAX;
	double dist;
	for(unsigned i=0;i<m_curves.size();i++)
	{
		for(unsigned j=0;j<m_curves[i].size();j++)
		{
			auto vec=point-m_curves[i][j];
			dist = vec.length();
			if(m_snapItem.dis>dist)
			{
				m_snapItem.dis=dist; m_snapItem.pos=m_curves[i][j];
				m_snapItem.type=m_curveVE[i][j].first; m_snapItem.id=m_curveVE[i][j].second;
				m_snapItem.cid=i; m_snapItem.vid=j;
			}
		}
	}
	if(m_snapItem.dis<0.01)
	{
		hitPoint = m_snapItem;
	}
	else
	{
		hitPoint.type = AnisGeodesic::VERTEX;
		hitPoint.id = m_interaction->hit.hittedNode;
		hitPoint.pos = points[m_interaction->hit.hittedNode];
	}

	// tracing a path and storing
	if(isNewNode && !m_PathIndex.empty())
	{// trace a free geodesy
		anisGeodesyHighlighted(hitPoint,false);
	}


	if(isStore)
	{
		breakCurvebySnapPoint(); //snap to a curve, and break the curve from snap point
		m_PathIndex.push_back(hitPoint);
		if(m_PathIndex.size()==2) //store a free geodesy
		{
			anisGeodesyHighlighted(m_PathIndex.back(), true);
		}

		anisGeodesyPropation(m_PathIndex.back()); //propogate from the new source

		if(AnisVisual::isDoSrcToAll)
		{
			srcToAllGeodesy();
		}
	}

	if (hitted)
		return m_interaction->hit.hittedNode;
	else
		return -1;
}
void PathUtils::pickCurve(int mouseX,int mouseY,bool isDel)
{
	if (m_curves.empty()) return;
	double minDistance = m_interaction->getClosestVertex(mouseX, mouseY, m_curves, selectedCurve);

	if(minDistance<0.01 && isDel) deleteCurrentCurve();
}
void PathUtils::deleteCurrentCurve()
{
	if(m_curves.empty() || m_curves.size()<selectedCurve.first ) return;
	m_curves.erase(m_curves.begin()+selectedCurve.first);
	m_curveVE.erase(m_curveVE.begin()+selectedCurve.first);
	selectedCurve.first=selectedCurve.second=0;
}
void PathUtils::deleteCurrentVertex()
{
	if (m_interaction->isControlPress)
	{
		if(!m_PathIndex.empty())
		{
			m_PathIndex.pop_back();
		}

		euclideanGeodesy->geoPath.clear(); ourGeodesy1->geoPath.clear(); ourGeodesy2->geoPath.clear(); pottmannGeodesy->geoPath.clear();
		kovacsGeodesy1->geoPath.clear(); kovacsGeodesy2->geoPath.clear(); campenGeodesy1->geoPath.clear(); campenGeodesy2->geoPath.clear();

		if(!m_PathIndex.empty())
		{
			anisGeodesyPropation(m_PathIndex.back());
		}
	}
	else if (m_interaction->isShiftPress)
	{
		m_interaction->scribbleTriangleStrip.clear();
		m_interaction->scribbleCurve.clear();
		if (!m_interaction->scribbleCurves.empty())
		{
			m_interaction->scribbleCurves.pop_back();
			m_interaction->scribbleTriangleStrips.pop_back();
		}
	}
}
void PathUtils::clearCurve()
{
	m_Path.clear();
	m_PathIndex.clear();
	m_VE.clear();

	euclideanGeodesy->geoPath.clear(); ourGeodesy1->geoPath.clear(); ourGeodesy2->geoPath.clear(); pottmannGeodesy->geoPath.clear();
	kovacsGeodesy1->geoPath.clear(); kovacsGeodesy2->geoPath.clear(); campenGeodesy1->geoPath.clear(); campenGeodesy2->geoPath.clear();

	m_interaction->scribbleTriangleStrip.clear();
	m_interaction->scribbleCurve.clear();
	m_interaction->scribbleCurves.clear();
	m_interaction->scribbleTriangleStrips.clear();
	m_snapItem.dis=1000;
}
void PathUtils::breakCurvebySnapPoint()
{
	if(m_snapItem.dis<0.01 && m_snapItem.cid<m_curves.size()
		&&m_snapItem.vid!=0 && m_snapItem.vid<m_curves[m_snapItem.cid].size()-1 )
	{//break the edge
		std::vector<Vec3> &dcur=m_curves[m_snapItem.cid];
		std::vector<Vec3> tcur(dcur.begin()+m_snapItem.vid,dcur.end());
		m_curves[m_snapItem.cid].resize(m_snapItem.vid+1);
		m_curves.push_back(tcur);

		std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> > &dVE = m_curveVE[m_snapItem.cid];
		std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> > tVE(dVE.begin() + m_snapItem.vid, dVE.end());
		m_curveVE[m_snapItem.cid].resize(m_snapItem.vid+1);
		m_curveVE.push_back(tVE);
	}
}
bool PathUtils::storeSelectedGeodesy()
{
	if(m_PathIndex.empty()) return false;

	if(m_anisAlgo == EUCLIDEAN)				{m_Path=euclideanGeodesy->geoPath; m_VE=euclideanGeodesy->geoPathVE;}
	else if(m_anisAlgo == OUR_MAX)			{m_Path=ourGeodesy1->geoPath; m_VE=ourGeodesy1->geoPathVE;}
	else if(m_anisAlgo == OUR_MIN)			{m_Path=ourGeodesy2->geoPath; m_VE=ourGeodesy2->geoPathVE;}
	else if(m_anisAlgo == POTTMANN)			{m_Path=pottmannGeodesy->geoPath; m_VE=pottmannGeodesy->geoPathVE;}
	else if(m_anisAlgo == KOVACS_MAX)		{m_Path=kovacsGeodesy1->geoPath; m_VE=kovacsGeodesy1->geoPathVE;}
	else if(m_anisAlgo == KOVACS_MIN)		{m_Path=kovacsGeodesy2->geoPath; m_VE=kovacsGeodesy2->geoPathVE;}
	else if(m_anisAlgo == CAMPEN_MAX)		{m_Path=campenGeodesy1->geoPath; m_VE=campenGeodesy1->geoPathVE;}
	else if(m_anisAlgo == CAMPEN_MIN)		{m_Path=campenGeodesy2->geoPath; m_VE=campenGeodesy2->geoPathVE;}

	if(m_Path.size()<=1) return true;

	m_curves.push_back(m_Path);
	m_curveVE.push_back(m_VE);

	m_PathIndex.erase(m_PathIndex.begin());
	m_Path.clear();
	m_VE.clear();
	euclideanGeodesy->geoPath.clear(); ourGeodesy1->geoPath.clear(); ourGeodesy2->geoPath.clear(); pottmannGeodesy->geoPath.clear();
	kovacsGeodesy1->geoPath.clear(); kovacsGeodesy2->geoPath.clear(); campenGeodesy1->geoPath.clear(); campenGeodesy2->geoPath.clear();

	return true;
}
void PathUtils::anisGeodesyHighlighted(SnapItem &tar, bool isStore)
{
	anisGeodesyTracing(tar);
	if(!isStore)
	{ //lowest cost path is highlighted
		std::vector<double> dis(8,DBL_MAX); 
		if(isDoEuclideanGeodesy)	dis[0] = euclideanGeodesy->geoPathCost;
		if(isDoOursGeodesy1)		dis[1] = ourGeodesy1->geoPathCost;
		if(isDoOursGeodesy2)		dis[2] = ourGeodesy2->geoPathCost;
		if(isDoPottmannsGeodesy)	dis[3] = pottmannGeodesy->geoPathCost;
		if(isDoKovacsGeodesy1)		dis[4] = kovacsGeodesy1->geoPathCost;
		if(isDoKovacsGeodesy2)		dis[5] = kovacsGeodesy2->geoPathCost;
		if(isDoCampensGeodesy1)		dis[6] = campenGeodesy1->geoPathCost;
		if(isDoCampensGeodesy2)		dis[7] = campenGeodesy2->geoPathCost;
		unsigned ind= std::min_element(dis.begin(),dis.end())-dis.begin();
		m_anisAlgo=MetricType(ind);
	}
	else
	{
		storeSelectedGeodesy();
	}
}

void PathUtils::metricActived(MetricType mt, bool isDo)
{
	switch(mt)
	{
	case EUCLIDEAN:	
		isDoEuclideanGeodesy=isDo;
		if(isDo)
		{
			if(euclideanGeodesy->anisMesh==NULL)
			{ //initialize
				euclideanGeodesy->anisInit_Euclidean();
				euclideanGeodesy->meshEmbedding_Euclidean();
				euclideanGeodesy->geodesyInit();
			}
			if(!m_PathIndex.empty())
			{ //if source exists, then propagate 
				SnapItem& s = m_PathIndex.back();
				euclideanGeodesy->geodesyPropagation(s.type,s.id,s.pos);
			}
		}
		break;

	case OUR_MAX:	
		isDoOursGeodesy1=isDo;
		if(isDo)
		{
			if(ourGeodesy1->anisMesh==NULL)
			{ //initialize
				ourGeodesy1->anisInit_OurMax(ourGeodesy1->metricVariation);
				ourGeodesy1->meshEmbedding_Tensor();
				ourGeodesy1->geodesyInit();
			}
			if(!m_PathIndex.empty()){ //if source exists, then propagate 
				SnapItem& s = m_PathIndex.back();
				ourGeodesy1->geodesyPropagation(s.type,s.id,s.pos);
			}
		}
		break;

	case OUR_MIN:	
		isDoOursGeodesy2=isDo;
		if(isDo)
		{
			if(ourGeodesy2->anisMesh==NULL)
			{ //initialize
				ourGeodesy2->anisInit_OurMin(ourGeodesy2->metricVariation);
				ourGeodesy2->meshEmbedding_Tensor();
				ourGeodesy2->geodesyInit();
			}
			if(!m_PathIndex.empty())
			{ //if source exists, then propagate 
				SnapItem& s = m_PathIndex.back();
				ourGeodesy2->geodesyPropagation(s.type,s.id,s.pos);
			}
		}
		break;

	case POTTMANN:	
		isDoPottmannsGeodesy=isDo;
		if(isDo)
		{
			if(pottmannGeodesy->anisMesh==NULL)
			{ //initialize
				pottmannGeodesy->anisInit_Pott();
				pottmannGeodesy->meshEmbedding_Pott();
				pottmannGeodesy->geodesyInit();
			}
			if(!m_PathIndex.empty()){ //if source exists, then propagate 
				SnapItem& s = m_PathIndex.back();
				pottmannGeodesy->geodesyPropagation(s.type,s.id,s.pos);
			}
		}
		break;

	case KOVACS_MAX:
		isDoKovacsGeodesy1=isDo;
		if(isDo)
		{
			if(kovacsGeodesy1->anisMesh==NULL)
			{ //initialize
				kovacsGeodesy1->anisInit_KovsMax();
				kovacsGeodesy1->meshEmbedding_Tensor();
				kovacsGeodesy1->geodesyInit();
			}
			if(!m_PathIndex.empty())
			{ //if source exists, then propagate 
				SnapItem& s = m_PathIndex.back();
				kovacsGeodesy1->geodesyPropagation(s.type,s.id,s.pos);
			}
		}
		break;

	case KOVACS_MIN:
		isDoKovacsGeodesy2=isDo;
		if(isDo)
		{
			if(kovacsGeodesy2->anisMesh==NULL)
			{ //initialize
				kovacsGeodesy2->anisInit_KovsMin();
				kovacsGeodesy2->meshEmbedding_Tensor();
				kovacsGeodesy2->geodesyInit();
			}
			if(!m_PathIndex.empty())
			{ //if source exists, then propagate 
				SnapItem& s = m_PathIndex.back();
				kovacsGeodesy2->geodesyPropagation(s.type,s.id,s.pos);
			}
		} 
		break;

	case CAMPEN_MAX:
		isDoCampensGeodesy1=isDo;
		if(isDo)
		{
			if(campenGeodesy1->anisMesh==NULL)
			{ //initialize
				campenGeodesy1->anisInit_CamMax();
				campenGeodesy1->meshEmbedding_Tensor();
				campenGeodesy1->geodesyInit();
			}
			if(!m_PathIndex.empty())
			{ //if source exists, then propagate 
				SnapItem& s = m_PathIndex.back();
				campenGeodesy1->geodesyPropagation(s.type,s.id,s.pos);
			}
		} 
		break;

	case CAMPEN_MIN:		
		isDoCampensGeodesy2=isDo;
		if(isDo)
		{
			if(campenGeodesy2->anisMesh==NULL)
			{ //initialize
				campenGeodesy2->anisInit_CamMin();
				campenGeodesy2->meshEmbedding_Tensor();
				campenGeodesy2->geodesyInit();
			}
			if(!m_PathIndex.empty())
			{ //if source exists, then propagate 
				SnapItem& s = m_PathIndex.back();
				campenGeodesy2->geodesyPropagation(s.type,s.id,s.pos);
			}
		}
		break;

	default: break;
	}
}

//visualization
void PathUtils::genDisplayList()
{
	AnisVisual::count = glGenLists(66);

	euclideanVisual->genDisplayList();
	ourVisual1->genDisplayList();
	ourVisual2->genDisplayList();
	pottmannVisual->genDisplayList();
	kovacsVisual1->genDisplayList();
	kovacsVisual2->genDisplayList();
	campenVisual1->genDisplayList();
	campenVisual2->genDisplayList();
}
void PathUtils::setMeshSelected(MetricType mt)
{
	euclideanVisual->showMesh=ourVisual1->showMesh=ourVisual2->showMesh=pottmannVisual->showMesh=
		kovacsVisual1->showMesh=kovacsVisual2->showMesh=campenVisual1->showMesh=campenVisual2->showMesh=false;
	if(mt==EUCLIDEAN)			euclideanVisual->showMesh=true;
	else if(mt==OUR_MAX)		ourVisual1->showMesh=true;
	else if(mt==OUR_MIN)		ourVisual2->showMesh=true;
	else if(mt==POTTMANN)		pottmannVisual->showMesh=true;
	else if(mt==KOVACS_MAX)		kovacsVisual1->showMesh=true;
	else if(mt==KOVACS_MIN)		kovacsVisual2->showMesh=true;
	else if(mt==CAMPEN_MAX)		campenVisual1->showMesh=true;
	else if(mt==CAMPEN_MIN)		campenVisual2->showMesh=true;
}
void PathUtils::setFieldLineSize(int val)
{ 
	if(val>0)
	{
		if(AnisVisual::shapeSize<5) 
			AnisVisual::shapeSize++;
		else		
			AnisVisual::shapeSize+=5; 
	}
	else if(val<0)
	{
		if(AnisVisual::shapeSize>5) 
			AnisVisual::shapeSize-=5;
		else if(AnisVisual::shapeSize>0)
			AnisVisual::shapeSize--;
	}

	if(euclideanVisual->showMesh)	{euclideanVisual->tensorLineList(0);euclideanVisual->tensorLineList(1);euclideanVisual->tensorEllipseList();}
	if(ourVisual1->showMesh)		{ourVisual1->tensorLineList(0);ourVisual1->tensorLineList(1);ourVisual1->tensorEllipseList();}
	if(ourVisual2->showMesh)		{ourVisual2->tensorLineList(0);ourVisual2->tensorLineList(1);ourVisual2->tensorEllipseList();}
	if(pottmannVisual->showMesh)	{pottmannVisual->tensorLineList(0);pottmannVisual->tensorLineList(1);pottmannVisual->tensorEllipseList();}
	if(kovacsVisual1->showMesh)		{kovacsVisual1->tensorLineList(0);kovacsVisual1->tensorLineList(1);kovacsVisual1->tensorEllipseList();}
	if(kovacsVisual2->showMesh)		{kovacsVisual2->tensorLineList(0);kovacsVisual2->tensorLineList(1);kovacsVisual2->tensorEllipseList();}
	if(campenVisual1->showMesh)		{campenVisual1->tensorLineList(0);campenVisual1->tensorLineList(1);campenVisual1->tensorEllipseList();}
	if(campenVisual2->showMesh)		{campenVisual2->tensorLineList(0);campenVisual2->tensorLineList(1);campenVisual2->tensorEllipseList();}
}
void PathUtils::setCurveSize(bool isIncrease)
{
	if(isIncrease)
	{
		if(lineWidth<50) 
			lineWidth++;
	}
	else
	{
		if(lineWidth>0)
			lineWidth--;
	}
	AnisVisual::lineWidth = (unsigned)(ceil(((double)lineWidth)/3.0));

	if(euclideanVisual->showMesh)	{euclideanVisual->isoLineList();euclideanVisual->allPathsList();}
	if(ourVisual1->showMesh)		{ourVisual1->isoLineList();ourVisual1->allPathsList();}
	if(ourVisual2->showMesh)		{ourVisual2->isoLineList();ourVisual2->allPathsList();}
	if(pottmannVisual->showMesh)	{pottmannVisual->isoLineList();pottmannVisual->allPathsList();}
	if(kovacsVisual1->showMesh)		{kovacsVisual1->isoLineList();kovacsVisual1->allPathsList();}
	if(kovacsVisual2->showMesh)		{kovacsVisual2->isoLineList();kovacsVisual2->allPathsList();}
	if(campenVisual1->showMesh)		{campenVisual1->isoLineList();campenVisual1->allPathsList();}
	if(campenVisual2->showMesh)		{campenVisual2->isoLineList();campenVisual2->allPathsList();}
}
void PathUtils::setShowFieldEllipse(bool val)
{
	showFieldEllipse=val;
	if(showFieldEllipse)
	{
		AnisVisual::tensorShapeType = 0;
		setFieldLineSize(0);
	}
}
void PathUtils::setShowFieldby8Shape(bool val)
{
	showField8Shape=val;
	if(showField8Shape)
	{
		AnisVisual::tensorShapeType = 1;
		setFieldLineSize(0);
	}
}

void PathUtils::genSegmentation()
{
	makePatches(myMesh,	m_curves,m_curveVs,	m_curveVE,	meshSegments,faceLabels,faceLabelsNum,cycleConnected,curveOfPatch);

	setPatchColors(false);
}

void preMakePatches(MyMesh* inMesh,
	std::vector<std::vector<Vec3> >& m_curves,
	std::vector<std::vector<unsigned> >& m_curveVs,
	std::vector<std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> > >& m_curveVE,
	MyMesh*& meshSegments,
	std::vector<std::vector<unsigned> >& intersetOnEdge)
{
	double pointCloseness = PathUtils::pointCloseness;
	cout << "pointCloseness:" << pointCloseness << endl;
	// for each edge, how many intersections lie on?
	meshSegments = new MyMesh(inMesh);
	intersetOnEdge.clear();intersetOnEdge.resize(meshSegments->getEdges().size());
	m_curveVs.clear();
	//if vertex on curves is too close to vertex on mesh, then snap to mesh;
	for (unsigned i = 0; i < m_curveVE.size(); i++)
	{
		for (unsigned j = 0; j < m_curveVE[i].size(); j++)
		{
			if (m_curveVE[i][j].first == AnisGeodesic::EDGE)
			{
				auto eit = meshSegments->getEIter()[m_curveVE[i][j].second];
				if ((m_curves[i][j] - eit->vertex_iter(0)->coordinate()).length() < pointCloseness)
				{
					m_curves[i][j] = eit->vertex_iter(0)->coordinate();
					m_curveVE[i][j].first = AnisGeodesic::VERTEX;
					m_curveVE[i][j].second = eit->vertex_iter(0)->id();
				}
				else if ((m_curves[i][j] - eit->vertex_iter(1)->coordinate()).length() < pointCloseness)
				{
					m_curves[i][j] = eit->vertex_iter(1)->coordinate();
					m_curveVE[i][j].first = AnisGeodesic::VERTEX;
					m_curveVE[i][j].second = eit->vertex_iter(1)->id();
				}
			}
		}
	}

	for (unsigned i = 0; i < m_curves.size(); i++)
	{
		std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> >  tVE(1, m_curveVE[i].front());
		std::vector<Vec3> tCV(1,m_curves[i].front());
		for (unsigned j = 1; j < m_curves[i].size(); j++)
		{
			if ((m_curves[i][j - 1] - m_curves[i][j]).length() > pointCloseness)
			{
				tCV.push_back(m_curves[i][j]);
				tVE.push_back(m_curveVE[i][j]);
			}
		}
		m_curves[i] = tCV;
		m_curveVE[i] = tVE;
	}

	// add vertex to new mesh;
	for (unsigned i = 0; i < m_curveVE.size(); i++)
	{
		std::vector<unsigned> curveV;
		for (unsigned j = 0; j < m_curveVE[i].size(); j++)
		{
			if (m_curveVE[i][j].first == AnisGeodesic::EDGE)
			{
				// add a new vertex;
				//but first check if the vertex already exists;
				bool isExist = false;
				std::pair<unsigned, unsigned> vpair;
				for (unsigned k = 0; k < i; k++)//compare with other curves
				{
					if (isExist == true)break;
					for (unsigned m = 0; m < m_curves[k].size(); m++)
					{
						if ((m_curves[k][m] - m_curves[i][j]).length() < pointCloseness)
						{
							isExist = true;
							vpair.first = k;
							vpair.second = m;
							break;
						}
					}
				}
				for (unsigned m = 0; m < curveV.size(); m++)//compare with current curve
				{
					if ((m_curves[i][m] - m_curves[i][j]).length() < pointCloseness)
					{
						vpair.first = i;
						vpair.second = m;
						isExist = true;
						break;
					}
				}

				if (isExist)
				{
					if (vpair.first == i)
					{
						curveV.push_back(curveV[vpair.second]);
					}
					else
					{
						curveV.push_back(m_curveVs[vpair.first][vpair.second]);
					}
				}
				else
				{	//a new vertex, add to meshSegments
					MyMesh::Vertex newVer;
					newVer.coordinate() = m_curves[i][j];
					newVer.id() = meshSegments->getVertices().size();
					meshSegments->getVertices().push_back(newVer);

					curveV.push_back(newVer.id());
					intersetOnEdge[m_curveVE[i][j].second].push_back(newVer.id());
				}
			}
			else
			{
				curveV.push_back(m_curveVE[i][j].second);
			}
		}
		m_curveVs.push_back(curveV);
	}//m_curveVs stores vertex index of geodesics on the subdivided mesh(meshSegments)
}

void makePatches(MyMesh* inMesh,
	std::vector<std::vector<Vec3> >& m_curves,
	std::vector<std::vector<unsigned> >& m_curveVs,
	std::vector<std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> > >& m_curveVE,
	MyMesh*& meshSegments,
	std::vector<unsigned >& faceLabels,
	unsigned& faceLabelsNum,
	std::vector<std::vector<bool> >& cycleConnected,
	std::vector<std::vector<unsigned> >& curveOfPatch)
{
	cout << "clean curves before making patches..." << endl;
	std::vector<std::vector<unsigned> > intersetOnEdge;
	preMakePatches(inMesh, m_curves, m_curveVs, m_curveVE, meshSegments, intersetOnEdge);

	cout << "subdivide triangles according to geodesics" << endl;
	//////////////////////////////////////////////////////////////////////////
	// for each face, how many line segments lie on?
	std::vector<std::vector<std::pair<unsigned,unsigned> > > faceSegs(meshSegments->getFaces().size());
	for(unsigned i=0;i<m_curveVE.size();i++)
	{
		for(unsigned j=1;j<m_curveVE[i].size();j++)
		{
			unsigned m = j-1; 
			if (m_curveVE[i][m].first == AnisGeodesic::VERTEX) m += 1;
			if (m_curveVE[i][m].first == AnisGeodesic::VERTEX) m += 1;
			if(m>j) continue; //on the edge;

			unsigned n = m==j? j-1:j;

			auto e_it = meshSegments->getEIter()[m_curveVE[i][m].second];
			auto f1_it = e_it->face_iter(0);
			auto f2_it = e_it->face_iter(1);

			try
			{
				if (m_curveVE[i][n].first == AnisGeodesic::VERTEX){//

					auto v_it = meshSegments->getVIter()[m_curveVE[i][n].second];
					if (f1_it->opposite_vertex(e_it)->id() == v_it->id())
					{
						faceSegs[f1_it->id()].push_back(std::pair<unsigned, unsigned>(m_curveVs[i][m], m_curveVs[i][n]));
					}
					else if (f2_it->opposite_vertex(e_it)->id() == v_it->id())
					{
						faceSegs[f2_it->id()].push_back(std::pair<unsigned, unsigned>(m_curveVs[i][m], m_curveVs[i][n]));
					}
					else
					{
						throw 1.0;
					}
				}
				else
				{
					auto e2_it = meshSegments->getEIter()[m_curveVE[i][n].second];
					auto v1_it = e_it->vertex_iter(0);
					auto v2_it = e_it->vertex_iter(1);
					if (f1_it->opposite_edge(v1_it)->id() == e2_it->id() || f1_it->opposite_edge(v2_it)->id() == e2_it->id())
					{
						faceSegs[f1_it->id()].push_back(std::pair<unsigned, unsigned>(m_curveVs[i][m], m_curveVs[i][n]));
					}
					else if (f2_it->opposite_edge(v1_it)->id() == e2_it->id() || f2_it->opposite_edge(v2_it)->id() == e2_it->id())
					{
						faceSegs[f2_it->id()].push_back(std::pair<unsigned, unsigned>(m_curveVs[i][m], m_curveVs[i][n]));
					}
					else
					{
						throw 2.0;
					}
				}
			}
			catch (double e)
			{
				if (e == 1.0)
					cout << "segment(1) doesn't lie on the triangle" << endl;
				else if (e == 2.0)
					cout << "segment(2) doesn't lie on the triangle" << endl;
			}
		}
	}

	std::vector<Vec3> newpoints;
	std::vector<std::vector<unsigned> > newfaces;
	for(auto iv_it=meshSegments->getVertices().begin();iv_it!=meshSegments->getVertices().end();iv_it++)
	{
		newpoints.push_back(iv_it->coordinate());
	}
	for(auto if_it=meshSegments->getFaces().begin();if_it!=meshSegments->getFaces().end();if_it++)
	{
		std::vector<unsigned> tvers(3);
		for(unsigned k=0;k<3;k++)	
		{
			tvers[k] = if_it->vertex_iter(k)->id();
		}
		newfaces.push_back(tvers);
	}

	// for each face , subdivide...	
	std::vector<unsigned> subTris;
	unsigned i=0;
	for (auto f_it = meshSegments->getFaces().begin();f_it!=meshSegments->getFaces().end();f_it++,i++)
	{
		if(faceSegs[i].empty())continue;

		subTris.push_back(i);
		//the order of the segment is very important, direct graph in triangle.
		MyMesh::EdgeIter e_it[]={f_it->edge_iter(0),f_it->edge_iter(1),f_it->edge_iter(2)};
		auto vs_it = e_it[0]->vertex_iter(0);
		auto ve_it = e_it[0]->vertex_iter(1);
		std::vector<std::pair<unsigned,unsigned> > edgeSegs;
		for(unsigned j=0;j<3;j++)
		{
			unsigned eid=e_it[j]->id();
			if(j>0)
			{
				unsigned strid=0;
				unsigned endid=1;
				if(vs_it->id()== e_it[j]->vertex_iter(0)->id() || ve_it->id()== e_it[j]->vertex_iter(1)->id())
				{
					strid=1;endid=0;
				}
				vs_it = e_it[j]->vertex_iter(strid);
				ve_it = e_it[j]->vertex_iter(endid);
			}
			std::vector<unsigned> tadjver = intersetOnEdge[eid];
			std::vector<double> tdis;
			for(unsigned k=0;k<tadjver.size();k++)
			{
				unsigned vid=tadjver[k];
				Vec3 tver = meshSegments->getVIter()[vid]->coordinate() - vs_it->coordinate();
				tdis.push_back(tver.length());
			}
			std::pair<unsigned,unsigned> tseg; tseg.first = vs_it->id();
			while(!tadjver.empty())
			{
				unsigned tind = std::min_element(tdis.begin(),tdis.end())-tdis.begin();
				tseg.second = tadjver[tind];
				edgeSegs.push_back(tseg);
				tseg.first = tseg.second;
				tdis.erase(tdis.begin()+tind);
				tadjver.erase(tadjver.begin()+tind);
			}
			tseg.second = ve_it->id();
			edgeSegs.push_back(tseg);
		}
		bool sameOrientation=false;
		for(unsigned j=0;j<3;j++)
		{
			if(e_it[0]->vertex_iter(0)->id() == f_it->vertex_iter(j)->id())
			{
				if(e_it[0]->vertex_iter(1)->id()==f_it->vertex_iter((j+1)%3)->id())
				{
					sameOrientation = true;
					break;
				}
			}
		}
		std::vector<std::vector<unsigned> > verOut;
		std::vector<unsigned > verInd;
		std::map<unsigned ,unsigned > mapPointToIndex;
		std::map<unsigned ,unsigned >::iterator itMap;
		for(unsigned j=0;j<faceSegs[i].size();j++)
		{
			unsigned &n1 = faceSegs[i][j].first;
			unsigned &n2 = faceSegs[i][j].second;

			itMap = mapPointToIndex.find(n1);
			if (itMap != mapPointToIndex.end())
			{
				verOut[itMap->second].push_back(n2);
			}
			else 
			{
				verOut.push_back(std::vector<unsigned>(1,n2));
				mapPointToIndex.insert(std::pair<unsigned,unsigned>(n1,verOut.size()-1));
				verInd.push_back(n1);
			}

			itMap = mapPointToIndex.find(n2);
			if (itMap != mapPointToIndex.end())
			{
				verOut[itMap->second].push_back(n1);
			}
			else 
			{
				verOut.push_back(std::vector<unsigned>(1,n1));
				mapPointToIndex.insert(std::pair<unsigned,unsigned>(n2,verOut.size()-1));
				verInd.push_back(n2);
			}
		}
		for(unsigned j=0;j<edgeSegs.size();j++)
		{
			unsigned &n1 = edgeSegs[j].first;
			unsigned &n2 = edgeSegs[j].second;

			itMap = mapPointToIndex.find(n1);
			if (itMap != mapPointToIndex.end())
			{
				verOut[itMap->second].push_back(n2);
			}
			else 
			{
				verOut.push_back(std::vector<unsigned>(1,n2));
				mapPointToIndex.insert(std::pair<unsigned,unsigned>(n1,verOut.size()-1));
				verInd.push_back(n1);
			}
		}
		//check if there are redundant edge;
		for(int j=0;j<verOut.size();j++)
		{
			for(int v1=0;v1<verOut[j].size()-1;v1++)
			{
				for(int v2=v1+1; v2<verOut[j].size();v2++)
				{
					if(verOut[j][v1]==verOut[j][v2])
					{
						verOut[j].erase(verOut[j].begin()+v2);
						itMap = mapPointToIndex.find(verOut[j][v1]);
						unsigned rev = itMap->second;
						for(int k=0;k<verOut[rev].size();k++)
						{
							if(verOut[rev][k]==verInd[j])
							{
								verOut[rev].erase(verOut[rev].begin()+k);
								break;
							}
						}
						v1=1000;break;
					}
				}
			}
		}

		unsigned firstind=0;
		for(unsigned j=1;j<verOut.size();j++)
		{
			if(verOut[j].size()>verOut[firstind].size())
				firstind = j;
		}
		std::vector<std::vector<unsigned> > verLoops;
		unsigned traceStr=verInd[firstind];
		bool visitAll=false;
		while(!visitAll)
		{
			std::vector<unsigned> verloop;
			verloop.push_back(traceStr);
			unsigned traceCurrent = traceStr;

			unsigned traceNext;
			itMap = mapPointToIndex.find(traceStr);
			traceNext=verOut[itMap->second].front();verOut[itMap->second].erase(verOut[itMap->second].begin());
			while (traceStr!=traceNext)
			{
				verloop.push_back(traceNext);
				unsigned traceTemp = traceNext;
				itMap = mapPointToIndex.find(traceNext);
				traceNext=verOut[itMap->second].front(); 
				if(traceCurrent==traceNext)
				{
					traceNext=verOut[itMap->second][1]; 
					verOut[itMap->second].erase(verOut[itMap->second].begin()+1);
				}
				else
				{
					verOut[itMap->second].erase(verOut[itMap->second].begin());
				}
				traceCurrent=traceTemp;
			}
			verLoops.push_back(verloop);

			visitAll=true;
			firstind=0;
			for(unsigned j=1;j<verOut.size();j++)
			{
				if(verOut[j].size()>verOut[firstind].size())
					firstind = j;
			}
			if(!verOut[firstind].empty())
			{
				traceStr=verOut[firstind].front();
				visitAll=false;
			}
		}

		//add new vertices and faces to `...
		if(sameOrientation==false)
		{
			for(unsigned j=0;j<verLoops.size();j++)
			{
				reverse(verLoops[j].begin(),verLoops[j].end());
			}
		}
		for(unsigned j=0;j<verLoops.size();j++)
		{
			if(verLoops[j].size()==3)
			{
				std::vector<unsigned> tface(3);
				for(unsigned k=0;k<3;k++) tface[k]=verLoops[j][k];
				newfaces.push_back(tface);
			}
			else
			{
				Vec3 newVertex(0,0,0);
				for(unsigned k=0;k<verLoops[j].size();k++){
					newVertex+= newpoints[verLoops[j][k]];
				}
				newVertex /= verLoops[j].size();
				newpoints.push_back(newVertex);
				unsigned n1 = newpoints.size()-1;
				verLoops[j].push_back(verLoops[j].front());
				for(unsigned k=1;k<verLoops[j].size();k++)
				{
					unsigned n2 = verLoops[j][k-1];
					unsigned n3 = verLoops[j][k];
					std::vector<unsigned> tface(3);
					tface[0]=n1;tface[1]=n2;tface[2]=n3;
					newfaces.push_back(tface);
				}
			}
		}
	}
	for(int j=subTris.size()-1;j>=0;j--)
	{
		newfaces.erase(newfaces.begin()+subTris[j]);
	}
	delete meshSegments;

	cout << "trace out segmented boundaries..." ;
	extractBoundary(meshSegments, newpoints, newfaces, m_curveVs, faceLabels, faceLabelsNum, cycleConnected, curveOfPatch);
	cout << "done" << endl;
}

void extractBoundary(MyMesh*& meshSegments,
	std::vector<Vec3> newpoints,
	std::vector<std::vector<unsigned> > newfaces,
	std::vector<std::vector<unsigned> >& m_curveVs,
	std::vector<unsigned >& faceLabels,
	unsigned& faceLabelsNum,
	std::vector<std::vector<bool> >& cycleConnected,
	std::vector<std::vector<unsigned> >& curveOfPatch)
{
	meshSegments = new MyMesh;
	meshSegments->initialize(newpoints, newfaces); /*	BL_Patch();return;*/

// 	faceLabelsNum=1;
// 	faceLabels.resize(meshSegments->getFaces().size(),0);
// 	return;
// 	savePatchesFile();
	//read this mesh file in MeshLab, and check if there exist non manifold edge or face, then you can find where the bad path lies.

	//ok, now we have subdivide mesh, and patch boundary, let's cut out the patches from mesh
	std::vector<std::vector<unsigned> > boundaryEdges;
	int i;
	for (i = 0; i<m_curveVs.size(); i++)
	{
		for (unsigned j = 1; j<m_curveVs[i].size(); j++)
		{
			std::vector<unsigned> teg(3);
			teg[0] = min(m_curveVs[i][j - 1], m_curveVs[i][j]);
			teg[1] = max(m_curveVs[i][j - 1], m_curveVs[i][j]);
			teg[2] = i; //boundary ancestor curve id;
			if (teg[0] == teg[1])continue;
			boundaryEdges.push_back(teg);
		}
	}
	std::sort(boundaryEdges.begin(), boundaryEdges.end());
	for (int i = 0; i<boundaryEdges.size() - 1; i++)
	{
		if (boundaryEdges[i] == boundaryEdges[i + 1])
		{
			boundaryEdges.erase(boundaryEdges.begin() + i);	i--;
		}
	}
	std::vector<bool> isboundaryEdges(meshSegments->getEdges().size(), false);
	std::vector<unsigned > boundaryAncestor(meshSegments->getEdges().size(), 0);
	i = 0;
	while (i<boundaryEdges.size())
	{
		for (auto ie_it = meshSegments->getEdges().begin(); ie_it != meshSegments->getEdges().end(); ie_it++)
		{
			unsigned id1 = ie_it->vertex_iter(0)->id();
			unsigned id2 = ie_it->vertex_iter(1)->id();
			if (id1 == boundaryEdges[i][0] && id2 == boundaryEdges[i][1])
			{
				isboundaryEdges[ie_it->id()] = true;
				boundaryAncestor[ie_it->id()] = boundaryEdges[i][2];
				i++;
				if (i >= boundaryEdges.size()) break;
			}
		}
		i++;
	}

	auto& vertices_iter = meshSegments->getVIter();
	auto& edges_iter = meshSegments->getEIter();
	auto& faces_iter = meshSegments->getFIter();

	unsigned faceNum, fcount;
	faceNum = meshSegments->getFaces().size();
	fcount = 0;
	std::vector<bool> visitTriangles(faceNum, false);
	faceLabels.resize(faceNum);
	unsigned labelId = 0;
	std::vector<std::vector<unsigned> > edgeShared(meshSegments->getEdges().size());
	i = 0;
	while (fcount<faceNum)
	{
		std::list<unsigned> faceQueue;
		for (unsigned j = 0; j<visitTriangles.size(); j++)
		{
			if (visitTriangles[j] == false)
			{
				faceQueue.push_back(j); visitTriangles[j] = true; break;
			}
		}
		//std::vector<unsigned> meshPatch;
		while (!faceQueue.empty())
		{
			unsigned fid = faceQueue.front(); faceQueue.erase(faceQueue.begin());
			//meshPatch.push_back(fid);
			faceLabels[fid] = labelId;
			fcount++;

			auto f_it = faces_iter[fid];
			for (unsigned k = 0; k<3; k++)
			{
				auto e_it = f_it->edge_iter(k);
				if (e_it->id() >= meshSegments->getEdges().size())
				{
					cout << "edge index over" << endl;
				}

				if (isboundaryEdges[e_it->id()] == false)
				{
					unsigned newfid = e_it->opposite_face(f_it)->id();
					if (visitTriangles[newfid] == false)
					{
						faceQueue.push_back(newfid);
						visitTriangles[newfid] = true;
					}
				}
				else
				{
					edgeShared[e_it->id()].push_back(labelId); i++;
				}
			}
		}
		labelId++;
	}
	faceLabelsNum = labelId;

	TriMesh* tmesh = new TriMesh;
	tmesh->vertices.resize(newpoints.size());
	for (int i = 0; i<newpoints.size(); i++)
	{
		for (int j = 0; j<3; j++)
		{
			tmesh->vertices[i][j] = newpoints[i][j];
		}
	}
	tmesh->faces.resize(newfaces.size());
	for (int i = 0; i<newfaces.size(); i++)
	{
		for (int j = 0; j<3; j++)
		{
			tmesh->faces[i][j] = newfaces[i][j];
		}
	}
	tmesh->need_normals();

	i = 0;
	for (auto v_it = meshSegments->getVertices().begin(); v_it != meshSegments->getVertices().end(); v_it++)
	{
		v_it->normal() = Vec3(tmesh->normals[i][0], tmesh->normals[i][1], tmesh->normals[i][2]);
		i++;
	}
	delete tmesh;

	cycleConnected.resize(faceLabelsNum, std::vector<bool>(faceLabelsNum, false));
	curveOfPatch.clear(); curveOfPatch.resize(faceLabelsNum);
	for (unsigned i = 0; i<edgeShared.size(); i++)
	{
		if (edgeShared[i].size() != 2) continue;
		unsigned &pat1 = edgeShared[i][0];
		unsigned &pat2 = edgeShared[i][1];
		cycleConnected[pat1][pat2] = true;
		cycleConnected[pat2][pat1] = true;

		unsigned &curId = boundaryAncestor[i];
		if (curveOfPatch[pat1].empty())
		{
			curveOfPatch[pat1].push_back(curId);
		}
		else if (std::find(curveOfPatch[pat1].begin(), curveOfPatch[pat1].end(), curId) == curveOfPatch[pat1].end())
		{
				curveOfPatch[pat1].push_back(curId);
		}
		if (curveOfPatch[pat2].empty())
		{
			curveOfPatch[pat2].push_back(curId);
		}
		else if (std::find(curveOfPatch[pat2].begin(), curveOfPatch[pat2].end(), curId) == curveOfPatch[pat2].end())
		{
			curveOfPatch[pat2].push_back(curId);
		}
	}
}
