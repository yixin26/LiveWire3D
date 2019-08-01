/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 **************************************************************/

#ifndef PATH_UTILS
#define PATH_UTILS

#include <vector>
#include <iostream>
#include <algorithm>

#include "AnisGeodesic.h"
#include "AnisVisual.h"
#include "Interaction.h"

#include <QtCore/QString>
#include <GL/glu.h>
#include <glut.h>
#include <set>

class TriMesh;

class PathUtils{
public:

	struct SnapItem
	{
		AnisGeodesic::IntersectionType	type;
		unsigned						id;
		double							dis;
		Vec3							pos;
		unsigned						cid;
		unsigned						vid;
	};

	PathUtils();
	~PathUtils();

	void init();

	//input&output
	TriMesh * getTrimesh(){ return triMesh; }
	MyMesh** getMyMesh(){return &myMesh;}
	void readTrimesh(const char* fileName);
	void rescaleMesh();
	bool ReadCurveData(const char* fileName);
	bool ReadPathData(const char* fileName);

	bool SavePath();
	bool saveCurvesFile();
	bool savePatchCycleFile();
	bool savePatchesFile();

	//visualization;
	void DrawCurve();
	void DrawGraph();

	MyInteraction** getInteractionUtils(){ return &m_interaction; }
	int pickMeshVertex(int mouseX, int mouseY, unsigned seletionMethod = 1, bool isStore = false);
	void pickCurve(int mouseX,int mouseY,bool isDel);

	//compute anis field and anis geodesy
	void computeCurvature();
	void computeAnisField();
	void updateCurvature();
	void updateMetrics();

	void globalParameter();

	//visualization;
	void setShowFieldEllipse(bool val);
	void setShowFieldby8Shape(bool val);

	void setFieldLineSize(int val);
	void setCurveSize(bool isIncrease);
	void setMeshSelected(MetricType sel);

	void metricActived(MetricType mt, bool isDo);

	//interaction
	void modifyField();
	void resetField(MetricType algoType);

	void deleteCurrentVertex();
	void deleteCurrentCurve();
	void clearCurve();
	void breakCurvebySnapPoint();

	//result;
	void genSegmentation();
	void setPatchColors(bool randomColor);

private:

	void anisGeodesyInit();
	void anisGeodesyPropation(SnapItem& s);

	void anisGeodesyTracing(SnapItem& t);
	void anisGeodesyHighlighted(SnapItem &tar, bool isStore);
	bool storeSelectedGeodesy();

	void drawPaths();
	void drawSelectNodesAndCurve();
	void drawCurveNetwork();
	void drawAnistropy(unsigned sel);
	void drawPatches(bool meshtype);

	void genDisplayList();

	void srcToAllGeodesy();

private:

	//initial
	std::vector<Vec3>					points;
	std::vector<std::vector<unsigned> > faces;
	MyMesh*								myMesh; 

	//compute curvatures;
	TriMesh*							triMesh;

	AnisVisual*							euclideanVisual;
	AnisVisual*							ourVisual1;
	AnisVisual*							ourVisual2;
	AnisVisual*							pottmannVisual;
	AnisVisual*							kovacsVisual1;
	AnisVisual*							kovacsVisual2;
	AnisVisual*							campenVisual1;
	AnisVisual*							campenVisual2;

	//geodesic path;
	std::vector<std::vector<Vec3> >		m_curves;
	std::vector<std::vector<unsigned> > m_curveVs;
	std::vector<std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> > > m_curveVE;

	std::vector<Vec3 >					m_Path;
	std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> > m_VE;
	std::vector<SnapItem>				m_PathIndex;

	//segmented patches
	MyMesh*								meshSegments;
	std::vector<std::vector<unsigned> > meshPatches;
	std::vector<unsigned >				faceLabels;
	unsigned							faceLabelsNum;
	std::vector<std::vector<bool> >		cycleConnected;
	std::vector<std::vector<unsigned> > curveOfPatch;
	std::vector<std::vector<Vec3> >		patchLoops;
	std::vector<std::vector<Vec3> >		patchLoopsScaled;
	std::vector<std::vector<Vec3> >		patchLoopsNormal;
	std::vector<ColorEngine::color>		patchColors;
	std::vector<unsigned>				patchColorsInd;

	//others
	SnapItem							m_snapItem; // a free geodesy will be snapped to a existed point(of curve network) if too close ;
	std::pair<unsigned,unsigned>		selectedCurve; //curve edit;

	Vec3								curveColor;
	Vec3 								meshColor;
	Vec3 								geoColor;
	Vec3 								ourColor1;
	Vec3 								ourColor2;
	Vec3 								pottColor;
	Vec3 								camColor1;
	Vec3 								camColor2;
	Vec3 								kovColor1;
	Vec3 								kovColor2;

public:

	QString							m_filename;

	MetricType							m_anisAlgo;
	AnisVisual*							anisVisual;
	AnisGeodesic*						anisGeodesy;

	AnisGeodesic*						euclideanGeodesy;
	AnisGeodesic*						ourGeodesy1;
	AnisGeodesic*						ourGeodesy2;
	AnisGeodesic*						pottmannGeodesy;
	AnisGeodesic*						kovacsGeodesy1;
	AnisGeodesic*						kovacsGeodesy2;
	AnisGeodesic*						campenGeodesy1;
	AnisGeodesic*						campenGeodesy2;

	MyInteraction*						m_interaction;

	char								curvatureScale[100];
	char								curvatureSmooth[100];

	bool			showCurveNetwork;
	bool			showMesh;
	unsigned		meshType,nodeSize,lineWidth;
	bool			showFieldDir1,showFieldDir2,showFieldEllipse,showField8Shape,showAnistropy,showCurvatureMax,showCurvatureMin; //fields

	bool			isDoEuclideanGeodesy;
	bool			isDoOursGeodesy1, isDoOursGeodesy2;
	bool			isDoPottmannsGeodesy;
	bool			isDoKovacsGeodesy1, isDoKovacsGeodesy2;
	bool			isDoCampensGeodesy1, isDoCampensGeodesy2;

	bool			showPatches;

	static double	pointCloseness;

	//debug;
	double edgeWidth;
	Vec3 edgeColor;
};

static void preMakePatches(MyMesh* inMesh,
	std::vector<std::vector<Vec3> >& m_curves,
	std::vector<std::vector<unsigned> >& m_curveVs,
	std::vector<std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> > >& m_curveVE,
	MyMesh*& meshSegments
	);
static void makePatches(MyMesh* inMesh, 
	std::vector<std::vector<Vec3> >& m_curves,
	std::vector<std::vector<unsigned> >& m_curveVs,
	std::vector<std::vector<std::pair<AnisGeodesic::IntersectionType, unsigned> > >& m_curveVE,
	MyMesh*& meshSegments,
	std::vector<unsigned >& faceLabels,
	unsigned& faceLabelsNum,
	std::vector<std::vector<bool> >& cycleConnected,
	std::vector<std::vector<unsigned> >& curveOfPatch
	);
static void extractBoundary(MyMesh*& meshSegments,
	std::vector<Vec3> newpoints,
	std::vector<std::vector<unsigned> > newfaces,
	std::vector<std::vector<unsigned> >& m_curveVs,
	std::vector<unsigned >& faceLabels,
	unsigned& faceLabelsNum,
	std::vector<std::vector<bool> >& cycleConnected,
	std::vector<std::vector<unsigned> >& curveOfPatch
	);

#endif