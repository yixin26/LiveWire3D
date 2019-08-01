/***************************************************************
* Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
**************************************************************/

#include "core/PathUtils.h"
#include <fstream>
#include <iomanip>
#include "TriMesh.h"
#include "TriMesh_algo.h"

using namespace std;
void PathUtils::rescaleMesh()
{
	cout << "mesh normalized to unit cube" << endl;

	MyMesh::BBox m_bbox;

	m_bbox.bBoxMin = Vec3(1e10, 1e10, 1e10);
	m_bbox.bBoxMax = Vec3(-1e10, -1e10, -1e10);
	m_bbox.centroid = Vec3(0, 0, 0);

	unsigned int numPoints = 0;

	for (unsigned i = 0; i<points.size(); i++)
	{
		auto &p = points[i];
		for (unsigned int d = 0; d<3; d++)
		{
			if (p[d]>m_bbox.bBoxMax[d])
				m_bbox.bBoxMax[d] = p[d];
			if (p[d]<m_bbox.bBoxMin[d])
				m_bbox.bBoxMin[d] = p[d];
		}
		m_bbox.centroid += p;
		numPoints++;
	}

	m_bbox.centroid *= 1.0 / numPoints;

	//find the sides of bbox
	Vec3 bBox;
	for (unsigned int d = 0; d<3; d++)
	{
		bBox[d] = abs(m_bbox.bBoxMin[d] - m_bbox.bBoxMax[d]);
	}

	//find the largest side
	double largestSide = (bBox[0]>bBox[1]) ? bBox[0] : bBox[1];
	largestSide = (largestSide > bBox[2]) ? largestSide : bBox[2];

	//translate the points so centroid is at origin
	//scale the points so bbox fits in canonical cube
	for (unsigned i = 0; i<points.size(); i++)
	{
		auto &p = points[i];
		p = p - m_bbox.centroid;
		//scale the points so bbox fits in canonical cube
		p *= 1.0 / largestSide;
	}

	m_bbox.isInit = true;
}
void PathUtils::readTrimesh(const char* fileName)
{
	triMesh = TriMesh::read(fileName);

	points.clear();	faces.clear();
	int vNum = triMesh->vertices.size();
	points.reserve(vNum);
	for (int i = 0; i<vNum; i++)
	{
		points.push_back(Vec3(triMesh->vertices[i][0], triMesh->vertices[i][1], triMesh->vertices[i][2]));
	}
	int fNum = triMesh->faces.size();
	faces.reserve(fNum);
	std::vector<unsigned> tf(3);
	for (int i = 0; i<fNum; i++)
	{
		tf[0] = triMesh->faces[i][0]; tf[1] = triMesh->faces[i][1]; tf[2] = triMesh->faces[i][2];
		faces.push_back(tf);
	}

	rescaleMesh();

	for (int i = 0; i<vNum; i++)
	{
		for (int j = 0; j<3; j++)
		{
			triMesh->vertices[i][j] = points[i][j];
		}
	}

	triMesh->need_normals();

	myMesh = new MyMesh;
	myMesh->initialize(points, faces);
	cout << "number of edges:" << myMesh->getEdges().size() << endl;

	computeCurvature();

	euclideanGeodesy = new AnisGeodesic;
	euclideanGeodesy->anisMetric = EUCLIDEAN;
	euclideanGeodesy->orgMesh = myMesh;

	ourGeodesy1 = new AnisGeodesic;
	ourGeodesy1->anisMetric = OUR_MAX;
	ourGeodesy1->scaleParameter = 0.05;
	ourGeodesy1->orgMesh = myMesh;

	ourGeodesy2 = new AnisGeodesic;
	ourGeodesy1->anisMetric = OUR_MIN;
	ourGeodesy2->scaleParameter = 0.05;
	ourGeodesy2->orgMesh = myMesh;


	pottmannGeodesy = new AnisGeodesic;
	pottmannGeodesy->anisMetric = POTTMANN;
	pottmannGeodesy->orgMesh = myMesh;

	kovacsGeodesy1 = new AnisGeodesic;
	kovacsGeodesy1->anisMetric = KOVACS_MAX;
	kovacsGeodesy1->orgMesh = myMesh;

	kovacsGeodesy2 = new AnisGeodesic;
	kovacsGeodesy2->anisMetric = KOVACS_MIN;
	kovacsGeodesy2->orgMesh = myMesh;

	campenGeodesy1 = new AnisGeodesic; campenGeodesy1->scaleParameter = 10;
	campenGeodesy1->anisMetric = CAMPEN_MAX;
	campenGeodesy1->orgMesh = myMesh;

	campenGeodesy2 = new AnisGeodesic; campenGeodesy2->scaleParameter = 10;
	campenGeodesy2->anisMetric = CAMPEN_MIN;
	campenGeodesy2->orgMesh = myMesh;

	euclideanVisual = new AnisVisual;
	euclideanVisual->anisGeodesy = euclideanGeodesy; euclideanVisual->showMesh = true;
	ourVisual1 = new AnisVisual;
	ourVisual1->anisGeodesy = ourGeodesy1;
	ourVisual2 = new AnisVisual;
	ourVisual2->anisGeodesy = ourGeodesy2;
	pottmannVisual = new AnisVisual;
	pottmannVisual->anisGeodesy = pottmannGeodesy;
	kovacsVisual1 = new AnisVisual;
	kovacsVisual1->anisGeodesy = kovacsGeodesy1;
	kovacsVisual2 = new AnisVisual;
	kovacsVisual2->anisGeodesy = kovacsGeodesy2;
	campenVisual1 = new AnisVisual;
	campenVisual1->anisGeodesy = campenGeodesy1;
	campenVisual2 = new AnisVisual;
	campenVisual2->anisGeodesy = campenGeodesy2;

	computeAnisField();
	anisGeodesyInit();
}

bool PathUtils::ReadCurveData(const char* fileName)
{
	if (!fileName)
		return false;

	std::ifstream reader(fileName);

	if (!reader.good())
		return false;

	m_curves.clear();

	unsigned int numCurves = 0;
	int numPoints, j, k;

	reader >> numCurves;
	if (numCurves <= 0) return false;

	//read point set data
	Vec3 point;

	for (k = 0; k < numCurves; k++)
	{
		reader >> numPoints >> j >> j;
		std::vector<Vec3> curve;
		for (int i = 0; i < numPoints; i++)
		{
			reader >> point.x >> point.y >> point.z;
			curve.push_back(point);
		}
		m_curves.push_back(curve);
	}

	reader.close();
	return true;
}
bool PathUtils::ReadPathData(const char* fileName)
{
	if (!fileName)
		return false;

	std::ifstream reader(fileName);

	if (!reader.good())
		return false;

	m_curves.clear(); m_curveVE.clear();
	int isHaveCur, isHaveMtr, isHaveNod, isHaveEdge;
	reader >> isHaveCur >> isHaveMtr >> isHaveNod >> isHaveEdge;
	if (isHaveMtr != 1 || isHaveNod != 1) isHaveMtr = isHaveNod = 0;

	unsigned int numCurves = 0;
	unsigned numPoints;

	if (isHaveMtr == 0)
	{
		reader.clear();
		reader.seekg(0, ios::beg);
	}
	reader >> numCurves;
	if (numCurves <= 0) return false;

	//read point set data
	Vec3 point; 
	unsigned mtr; 
	AML::uint2 nod; 
	unsigned i1,i2;
	for (unsigned k = 0; k < numCurves; k++)
	{
		reader >> numPoints;
		if (isHaveMtr == 1) reader >> mtr;
		if (isHaveNod == 1) reader >> nod.x >> nod.y;
		std::vector<Vec3> curve;
		std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> > cVE;
		for (int i = 0; i < numPoints; i++) 
		{
			reader >> point.x >> point.y >> point.z;
			reader >> i1 >> i2;
			curve.push_back(point);
			cVE.push_back(std::pair<AnisGeodesic::IntersectionType, unsigned>(AnisGeodesic::IntersectionType(i1), i2));
		}
		m_curves.push_back(curve);
		m_curveVE.push_back(cVE);
	}
	reader.close();
	return true;
}
bool PathUtils::SavePath()
{
	QString FileName = m_filename;
	FileName.replace(FileName.lastIndexOf(".") + 1, 3, "path");

	std::ofstream File;
	File.open(FileName.toStdString(), std::ios::out);
	if (!File.good())
		return false;

	File << m_curves.size() << endl;
	for (int i = 0; i < m_curves.size(); i++)
	{
		File << m_curves[i].size() << endl;
		for (int j = 0; j < m_curves[i].size(); j++)
		{
			Vec3 &point = m_curves[i][j];
			File << point[0] << ' ' << point[1] << ' ' << point[2] << ' ' << m_curveVE[i][j].first << ' ' << m_curveVE[i][j].second << endl;
		}
	}
	File.close();
	return true;
}
bool PathUtils::saveCurvesFile()
{
	QString FileName = m_filename;
	FileName.replace(FileName.lastIndexOf(".") + 1, 3, "curve");

	std::ofstream File;
	File.open(FileName.toStdString(), std::ios::out);
	if (!File.good())
		return false;

	File << m_curves.size() << endl;
	for (int i = 0; i < m_curves.size(); i++)
	{
		File << m_curves[i].size() << endl;
		for (int j = 0; j < m_curves[i].size(); j++)
		{
			Vec3 &point = m_curves[i][j];
			File << point[0] << ' ' << point[1] << ' ' << point[2] << ' ' << endl;
		}
	}
	File.close();
	return true;
}
bool PathUtils::savePatchCycleFile()
{
	QString FileName = m_filename;
	FileName.replace(FileName.lastIndexOf(".") + 1, 3, "cycle");

	std::ofstream File;
	File.open(FileName.toStdString(), std::ios::out);
	if (!File.good())
		return false;

	if (!patchLoops.empty())
	{
		File << patchLoops.size() << endl;
		for (unsigned i = 0; i < patchLoops.size(); i++)
		{
			File << patchLoops[i].size() << endl;
			for (unsigned j = 0; j < patchLoops[i].size(); j++)
			{
				File << patchLoops[i][j][0] << ' ' << patchLoops[i][j][1] << ' ' << patchLoops[i][j][2] << ' ';
				File << patchLoopsNormal[i][j][0] << ' ' << patchLoopsNormal[i][j][1] << ' ' << patchLoopsNormal[i][j][2] << endl;
			}
		}
	}

	File.close();
}
bool PathUtils::savePatchesFile()
{
	QString FileName = m_filename;
	FileName.replace(FileName.lastIndexOf(".") + 1, 3, "obj");

	std::ofstream File;
	File.open(FileName.toStdString(), std::ios::out);
	if (!File.good())
		return false;

	for (MyMesh::VertexIter v_it = meshSegments->getVertices().begin(); v_it != meshSegments->getVertices().end(); v_it++)
	{
		File << "v " << v_it->coordinate().x << ' ' << v_it->coordinate().y << ' ' << v_it->coordinate().z << endl;
	}

	auto& faces_iter=meshSegments->getFIter();
// 	for (auto if_it = meshSegments->getFaces().begin(); if_it != meshSegments->getFaces().end(); if_it++)
// 	{
// 		faces_iter.push_back(if_it);
// 	}

	meshPatches.resize(faceLabelsNum);
	for (unsigned i = 0; i < faceLabels.size(); i++)
	{
		meshPatches[faceLabels[i]].push_back(i);
	}
	for (unsigned i = 0; i < meshPatches.size(); i++)
	{
		File << "g " << "group" << i << endl;
		for (unsigned j = 0; j < meshPatches[i].size(); j++)
		{
			auto f_it = faces_iter[meshPatches[i][j]];
			File << "f " << f_it->vertex_iter(0)->id() + 1 << ' ' << f_it->vertex_iter(1)->id() + 1 << ' ' << f_it->vertex_iter(2)->id() + 1 << endl;
		}
	}

	File.close();
}

