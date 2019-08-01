/***************************************************************
* Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
**************************************************************/

#include "core/PathUtils.h"
#include <algorithm>

void PathUtils::drawSelectNodesAndCurve()
{
	if (m_interaction->isAltPress && !m_curves.empty())
	{
		glColor3f(1, 1, 0);
		glLineWidth(double(lineWidth) / 2.);
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j<m_curves[selectedCurve.first].size(); ++j)
		{
			Vec3 p = m_curves[selectedCurve.first][j];
			glVertex3d(p[0], p[1], p[2]);
		}
		glEnd();//GL_LINE_STRIP
		glLineWidth(1.);
	}

	if (m_interaction->isControlPress)
	{
		if (m_snapItem.dis<0.01)
		{//draw snap item if a current selected vertex is too close to exsiting curve node.
			glPushMatrix();
			glTranslated(m_snapItem.pos.x, m_snapItem.pos.y, m_snapItem.pos.z);
			glColor3d(0., 1., 1.);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glutSolidSphere(.002*nodeSize, 50, 50);
			glPopMatrix();
		}
		if (m_snapItem.dis >= 0.01)
		{//draw current selected vertex;
			glPushMatrix();
			glTranslated(points[m_interaction->hit.hittedNode].x, points[m_interaction->hit.hittedNode].y, points[m_interaction->hit.hittedNode].z);
			glColor3d(0., 0., 1.);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glutSolidSphere(.001*nodeSize, 50, 50);
			glPopMatrix();
		}
		if (!m_PathIndex.empty())
		{
			SnapItem& index = m_PathIndex.back();
			glPushMatrix();
			glTranslated(index.pos.x, index.pos.y, index.pos.z);
			glColor3d(1., 0., 0.);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glutSolidSphere(.001*nodeSize, 50, 50);
			glPopMatrix();
		}
	}

	if (m_interaction->isShiftPress && false)
	{
		//draw sketch curves on mesh
		glColor3d(1, 0, 0);
		glLineWidth(1.0);

		glBegin(GL_LINE_STRIP);
		for (unsigned i = 0; i<m_interaction->scribbleCurve.size(); i++)
		{
			glVertex3d(m_interaction->scribbleCurve[i].x, m_interaction->scribbleCurve[i].y, m_interaction->scribbleCurve[i].z);
		}
		glEnd();

		for (unsigned c = 0; c<m_interaction->scribbleCurves.size(); c++)
		{
			glBegin(GL_LINE_STRIP);
			for (unsigned i = 0; i<m_interaction->scribbleCurves[c].size(); i++)
			{
				glVertex3d(m_interaction->scribbleCurves[c][i].x, m_interaction->scribbleCurves[c][i].y, m_interaction->scribbleCurves[c][i].z);
			}
			glEnd();
		}


		//draw cross sign, indicating the mouse location
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();

		glTranslated(m_interaction->hit.hittedNearVertex.x, m_interaction->hit.hittedNearVertex.y, m_interaction->hit.hittedNearVertex.z);
		float modelview_matrix[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, modelview_matrix);
		for (int i = 0; i<12; i++)
		{
			if (i % 5 == 0) modelview_matrix[i] = 1;
			else modelview_matrix[i] = 0;
		}
		modelview_matrix[15] = 1;
		glLoadIdentity();
		glMultMatrixf(modelview_matrix);
		double eyeChange = modelview_matrix[14] / -2.0;
		glScalef(eyeChange, eyeChange, eyeChange);

		glBegin(GL_LINES);
		glVertex2f(-0.015, 0);
		glVertex2f(+0.015, 0);
		glVertex2f(0, -0.015);
		glVertex2f(0, +0.015);
		glEnd();

		glPopMatrix();
	}

}
void PathUtils::drawPaths()
{
	if (true/*isControlPress*/)
	{
		double lsize = double(lineWidth) / 2;
		double l1 = lsize + 1.; double l21 = lsize; double l22 = lsize; double l30 = lsize;
		double l31 = lsize;	double l32 = lsize; double l41 = lsize; double l42 = lsize;

		bool s1, s21, s22, s30, s31, s32, s41, s42; s1 = s30 = false; s21 = true; s22 = false; s31 = true; s32 = false; s41 = true; s42 = false;
		if (m_anisAlgo == OUR_MAX)	{ l1 = lsize; l21 = lsize + 2.; s21 = false; s22 = true; }
		else if (m_anisAlgo == OUR_MIN){ l1 = lsize; l22 = lsize + 2.; }
		else if (m_anisAlgo == POTTMANN){ l1 = lsize; l30 = lsize + 2.; }
		else if (m_anisAlgo == KOVACS_MAX){ l1 = lsize; l31 = lsize + 2.; s31 = false; s32 = true; }
		else if (m_anisAlgo == KOVACS_MIN){ l1 = lsize; l32 = lsize + 2.; }
		else if (m_anisAlgo == CAMPEN_MAX){ l1 = lsize; l41 = lsize + 2.; s41 = false; s42 = true; }
		else if (m_anisAlgo == CAMPEN_MIN){ l1 = lsize; l42 = lsize + 2.; }

		Vec3 color5 = Vec3(.9, .1, .9); double l5 = 4.0;

		glLineStipple(2, 0x0F0F);
		if (isDoEuclideanGeodesy)
		{
			if (s1) 	glEnable(GL_LINE_STIPPLE);
			glLineWidth(l1);
			glBegin(GL_LINE_STRIP);
			glColor3d(geoColor.x, geoColor.y, geoColor.z);
			std::vector<Vec3>& m_geodesicPathForView = euclideanGeodesy->geoPath;
			for (unsigned i = 0; i<m_geodesicPathForView.size(); i++)
			{
				Vec3 position(m_geodesicPathForView[i][0], m_geodesicPathForView[i][1], m_geodesicPathForView[i][2]);
				glVertex3d(position.x, position.y, position.z);
			}
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		}
		if (isDoOursGeodesy1)
		{
			if (s21) 	glEnable(GL_LINE_STIPPLE);
			glLineWidth(l21);
			glBegin(GL_LINE_STRIP);
			glColor3d(ourColor1.x, ourColor1.y, ourColor1.z);
			std::vector<Vec3>& m_curvaturePathForView1 = ourGeodesy1->geoPath;
			for (unsigned i = 0; i<m_curvaturePathForView1.size(); i++)
			{
				Vec3 position(m_curvaturePathForView1[i][0], m_curvaturePathForView1[i][1], m_curvaturePathForView1[i][2]);
				glVertex3d(position.x, position.y, position.z);
			}
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		}
		if (isDoOursGeodesy2)
		{
			if (s22) 	glEnable(GL_LINE_STIPPLE);
			glLineWidth(l22);
			glBegin(GL_LINE_STRIP);
			glColor3d(ourColor2.x, ourColor2.y, ourColor2.z);
			std::vector<Vec3>& m_curvaturePathForView2 = ourGeodesy2->geoPath;
			for (unsigned i = 0; i<m_curvaturePathForView2.size(); i++)
			{
				Vec3 position(m_curvaturePathForView2[i][0], m_curvaturePathForView2[i][1], m_curvaturePathForView2[i][2]);
				glVertex3d(position.x, position.y, position.z);
			}
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		}
		if (isDoPottmannsGeodesy)
		{
			if (s30) 	glEnable(GL_LINE_STIPPLE);
			glLineWidth(l30);
			glBegin(GL_LINE_STRIP);
			glColor3d(pottColor.x, pottColor.y, pottColor.z);
			std::vector<Vec3>& m_normalPathForView = pottmannGeodesy->geoPath;
			for (unsigned i = 0; i<m_normalPathForView.size(); i++)
			{
				Vec3 position(m_normalPathForView[i][0], m_normalPathForView[i][1], m_normalPathForView[i][2]);
				glVertex3d(position.x, position.y, position.z);
			}
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		}
		if (isDoKovacsGeodesy1)
		{
			if (s31) 	glEnable(GL_LINE_STIPPLE);
			glLineStipple(2, 0x0F0F);
			glLineWidth(l31);
			glBegin(GL_LINE_STRIP);
			glColor3d(kovColor1.x, kovColor1.y, kovColor1.z);
			std::vector<Vec3>& m_kovacsPathForView1 = kovacsGeodesy1->geoPath;
			for (unsigned i = 0; i<m_kovacsPathForView1.size(); i++)
			{
				Vec3 position(m_kovacsPathForView1[i][0], m_kovacsPathForView1[i][1], m_kovacsPathForView1[i][2]);
				glVertex3d(position.x, position.y, position.z);
			}
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		}
		if (isDoKovacsGeodesy2)
		{
			if (s32) 	glEnable(GL_LINE_STIPPLE);
			glLineWidth(l32);
			glBegin(GL_LINE_STRIP);
			glColor3d(kovColor2.x, kovColor2.y, kovColor2.z);
			std::vector<Vec3>& m_kovacsPathForView2 = kovacsGeodesy2->geoPath;
			for (unsigned i = 0; i<m_kovacsPathForView2.size(); i++)
			{
				Vec3 position(m_kovacsPathForView2[i][0], m_kovacsPathForView2[i][1], m_kovacsPathForView2[i][2]);
				glVertex3d(position.x, position.y, position.z);
			}
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		}
		if (isDoCampensGeodesy1)
		{
			if (s41) 	glEnable(GL_LINE_STIPPLE);
			glLineStipple(2, 0x0F0F);
			glLineWidth(l41);
			glBegin(GL_LINE_STRIP);
			glColor3d(camColor1.x, camColor1.y, camColor1.z);
			std::vector<Vec3>& m_campenPathForView1 = campenGeodesy1->geoPath;
			for (unsigned i = 0; i<m_campenPathForView1.size(); i++)
			{
				Vec3 position(m_campenPathForView1[i][0], m_campenPathForView1[i][1], m_campenPathForView1[i][2]);
				glVertex3d(position.x, position.y, position.z);
			}
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		}
		if (isDoCampensGeodesy2)
		{
			if (s42) 	glEnable(GL_LINE_STIPPLE);
			glLineWidth(l42);
			glBegin(GL_LINE_STRIP);
			glColor3d(camColor2.x, camColor2.y, camColor2.z);
			std::vector<Vec3>& m_campenPathForView2 = campenGeodesy2->geoPath;
			for (unsigned i = 0; i<m_campenPathForView2.size(); i++)
			{
				Vec3 position(m_campenPathForView2[i][0], m_campenPathForView2[i][1], m_campenPathForView2[i][2]);
				glVertex3d(position.x, position.y, position.z);
			}
			glEnd();
			glDisable(GL_LINE_STIPPLE);
		}
		glLineWidth(1.0);
	}
}
void PathUtils::drawCurveNetwork()
{
	glLineWidth(double(lineWidth) / 2);
	glColor3d(curveColor.x, curveColor.y, curveColor.z);
	for (int i = 0; i<m_curves.size(); i++)
	{
		std::vector<Vec3>& points = m_curves[i];
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j<points.size(); j++)
		{
			Vec3 p = points[j];
			glVertex3d(p[0], p[1], p[2]);
		}
		glEnd();
	}
	glLineWidth(1.);

	if (nodeSize != 0)
	{
		glPointSize((nodeSize)* 2);
		glBegin(GL_POINTS);
		glColor3d(0., 1., 1.);
		for (int i = 0; i<m_curves.size(); i++)
		{
			if (m_curves[i].size() <= 1)continue;
			Vec3& strp = m_curves[i].front();
			Vec3& endp = m_curves[i].back();

			glVertex3d(strp.x, strp.y, strp.z);
			glVertex3d(endp.x, endp.y, endp.z);
		}
		glEnd();
	}

}
void PathUtils::drawPatches(bool meshtype)
{
	if (meshSegments == NULL) return;

	if (meshtype == 2)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	GLfloat specular_light0[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular_light0);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular_light0);
	GLfloat As[4] = { 0.05f, 0.05f, 0.05f, 1.0f };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, As);

	glBegin(GL_TRIANGLES);
	for (auto f_it = meshSegments->getFaces().begin(); f_it != meshSegments->getFaces().end(); f_it++)
	{
		MyMesh::VertexIter v[] = { f_it->vertex_iter(0), f_it->vertex_iter(1), f_it->vertex_iter(2) };
		ColorEngine::color &tcolor = patchColors[patchColorsInd[faceLabels[f_it->id()]]];
		glColor3d(tcolor.r, tcolor.b, tcolor.g);
		for (unsigned j = 0; j<3; j++)
		{
			glNormal3f(v[j]->normal().x, v[j]->normal().y, v[j]->normal().z);
			glVertex3f(v[j]->coordinate().x, v[j]->coordinate().y, v[j]->coordinate().z);
		}
	}
	glEnd();

	GLfloat specular_light[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular_light);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular_light);
	GLfloat As2[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, As2);

}
void PathUtils::drawAnistropy(unsigned sel)
{
	if (myMesh == NULL) return;

	std::vector<double> anis(myMesh->getVertices().size());
	if (sel == 0)
	{
		auto v_it = myMesh->getVertices().begin();
		for (unsigned i = 0; i<anis.size(); i++)
		{
			anis[i] = fabs(v_it->magnitude(0) - v_it->magnitude(1)); v_it++;
		}
	}
	else if (sel == 1)
	{
		auto v_it = myMesh->getVertices().begin();
		for (unsigned i = 0; i<anis.size(); i++)
		{
			anis[i] = v_it->magnitude(0); v_it++;
		}
	}
	else if (sel == 2)
	{
		auto v_it = myMesh->getVertices().begin();
		for (unsigned i = 0; i<anis.size(); i++)
		{
			anis[i] = v_it->magnitude(1); v_it++;
		}
	}

	double maxL = *std::max_element(anis.begin(), anis.end());
	double minL = *std::min_element(anis.begin(), anis.end());

	maxL -= minL; maxL = 0.66 / maxL;
	std::vector<double> mags;
	for (unsigned i = 0; i<anis.size(); i++)
	{
		mags.push_back(0.66 - (anis[i] - minL)*maxL);//reverse color
	}

	std::vector<ColorEngine::color> colorsList;
	ColorEngine::HslToRgb(mags, colorsList);

	glDisable(GL_LIGHTING);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_TRIANGLES);
	for (auto f_it = myMesh->getFaces().begin(); f_it != myMesh->getFaces().end(); f_it++)
	{
		MyMesh::VertexIter v[] = { f_it->vertex_iter(0), f_it->vertex_iter(1), f_it->vertex_iter(2) };
		for (unsigned j = 0; j<3; j++)
		{
			glColor3f(colorsList[v[j]->id()].r, colorsList[v[j]->id()].g, colorsList[v[j]->id()].b);
			glNormal3f(v[j]->normal().x, v[j]->normal().y, v[j]->normal().z);
			glVertex3f(v[j]->coordinate().x, v[j]->coordinate().y, v[j]->coordinate().z);
		}
	}
	glEnd();
}
void PathUtils::DrawCurve()
{
	if (points.empty())	return;

	if (showMesh)
	{
		glDisable(GL_LIGHTING);
		if (euclideanVisual->showMesh)
		{
			if (meshType == 0 || meshType == 2)	euclideanVisual->drawWireMesh(edgeWidth,edgeColor);
			if (showFieldDir1) euclideanVisual->drawDirection(0);
			if (showFieldDir2) euclideanVisual->drawDirection(1);
			if (AnisVisual::showAllPaths)	euclideanVisual->drawAllPaths();
			if (AnisVisual::showIsoline)		euclideanVisual->drawIsoline();
		}
		else if (ourVisual1->showMesh)
		{
			if (meshType == 0 || meshType == 2)	ourVisual1->drawWireMesh();
			if (showFieldDir1) ourVisual1->drawDirection(0);
			if (showFieldDir2) ourVisual1->drawDirection(1);
			if (AnisVisual::showAllPaths)	ourVisual1->drawAllPaths();
			if (AnisVisual::showIsoline)		ourVisual1->drawIsoline();
		}
		else if (ourVisual2->showMesh)
		{
			if (meshType == 0 || meshType == 2)	ourVisual2->drawWireMesh();
			if (showFieldDir1) ourVisual2->drawDirection(0);
			if (showFieldDir2) ourVisual2->drawDirection(1);
			if (AnisVisual::showAllPaths)	ourVisual2->drawAllPaths();
			if (AnisVisual::showIsoline)		ourVisual2->drawIsoline();
		}
		else if (pottmannVisual->showMesh)
		{
			if (meshType == 0 || meshType == 2)	pottmannVisual->drawWireMesh();
			if (showFieldDir1) pottmannVisual->drawDirection(0);
			if (showFieldDir2) pottmannVisual->drawDirection(1);
			if (AnisVisual::showAllPaths)	pottmannVisual->drawAllPaths();
			if (AnisVisual::showIsoline)		pottmannVisual->drawIsoline();
		}
		else if (kovacsVisual1->showMesh)
		{
			if (meshType == 0 || meshType == 2)	kovacsVisual1->drawWireMesh();
			if (showFieldDir1) kovacsVisual1->drawDirection(0);
			if (showFieldDir2) kovacsVisual1->drawDirection(1);
			if (AnisVisual::showAllPaths)	kovacsVisual1->drawAllPaths();
			if (AnisVisual::showIsoline)		kovacsVisual1->drawIsoline();
		}
		else if (kovacsVisual2->showMesh)
		{
			if (meshType == 0 || meshType == 2)	kovacsVisual2->drawWireMesh();
			if (showFieldDir1) kovacsVisual2->drawDirection(0);
			if (showFieldDir2) kovacsVisual2->drawDirection(1);
			if (AnisVisual::showAllPaths)	kovacsVisual2->drawAllPaths();
			if (AnisVisual::showIsoline)		kovacsVisual2->drawIsoline();
		}
		else if (campenVisual1->showMesh)
		{
			if (meshType == 0 || meshType == 2)	campenVisual1->drawWireMesh();
			if (showFieldDir1) campenVisual1->drawDirection(0);
			if (showFieldDir2) campenVisual1->drawDirection(1);
			if (AnisVisual::showAllPaths)	campenVisual1->drawAllPaths();
			if (AnisVisual::showIsoline)		campenVisual1->drawIsoline();
		}
		else if (campenVisual2->showMesh)
		{
			if (meshType == 0 || meshType == 2)	campenVisual2->drawWireMesh();
			if (showFieldDir1) campenVisual2->drawDirection(0);
			if (showFieldDir2) campenVisual2->drawDirection(1);
			if (AnisVisual::showAllPaths)	campenVisual2->drawAllPaths();
			if (AnisVisual::showIsoline)		campenVisual2->drawIsoline();
		}
	}
	if (showCurveNetwork)
	{
		drawCurveNetwork();
		drawPaths();
		drawSelectNodesAndCurve();
	}
}
void PathUtils::DrawGraph()
{
	if (points.empty())	return;

	if (showMesh)
	{
		glEnable(GL_LIGHTING);
		if (showAnistropy)
			drawAnistropy(0);
		else if (showCurvatureMax)
			drawAnistropy(1);
		else if (showCurvatureMin)
			drawAnistropy(2);
		else if (showPatches)
		{
			GLfloat mat_emit[] = { 0.1f, 0.1f, .1f, 1.0f };
			glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emit);

			if (meshType == 1 || meshType == 2)
				drawPatches(2);
			if (meshType == 0 || meshType == 2)
				drawPatches(1);
		}
		else 
		{
			bool fieldInited = false; // if draw distance field, and field is not yet initialized, then draw solid mesh;
			if (euclideanVisual->showMesh)
			{
				if (AnisVisual::showDistanceField)	fieldInited = euclideanVisual->drawDistanceField();
				if ((meshType == 1 || meshType == 2) && !fieldInited)	euclideanVisual->drawSolidMesh();
				if (showFieldEllipse || showField8Shape) euclideanVisual->drawEllipse();
			}
			else if (ourVisual1->showMesh)
			{
				if (AnisVisual::showDistanceField)	fieldInited = ourVisual1->drawDistanceField();
				if ((meshType == 1 || meshType == 2) && !fieldInited)	ourVisual1->drawSolidMesh();
				if (showFieldEllipse || showField8Shape) ourVisual1->drawEllipse();
			}
			else if (ourVisual2->showMesh)
			{
				if (AnisVisual::showDistanceField)	fieldInited = ourVisual2->drawDistanceField();
				if ((meshType == 1 || meshType == 2) && !fieldInited)	ourVisual2->drawSolidMesh();
				if (showFieldEllipse || showField8Shape) ourVisual2->drawEllipse();
			}
			else if (pottmannVisual->showMesh)
			{
				if (AnisVisual::showDistanceField)	pottmannVisual->drawDistanceField();
				else if (meshType == 1 || meshType == 2)	pottmannVisual->drawSolidMesh();
				if (showFieldEllipse || showField8Shape) pottmannVisual->drawEllipse();
			}
			else if (kovacsVisual1->showMesh)
			{
				if (AnisVisual::showDistanceField)	kovacsVisual1->drawDistanceField();
				else if (meshType == 1 || meshType == 2)	kovacsVisual1->drawSolidMesh();
				if (showFieldEllipse || showField8Shape) kovacsVisual1->drawEllipse();
			}
			else if (kovacsVisual2->showMesh)
			{
				if (AnisVisual::showDistanceField)	kovacsVisual2->drawDistanceField();
				else if (meshType == 1 || meshType == 2)	kovacsVisual2->drawSolidMesh();
				if (showFieldEllipse || showField8Shape) kovacsVisual2->drawEllipse();
			}
			else if (campenVisual1->showMesh)
			{
				if (AnisVisual::showDistanceField)	campenVisual1->drawDistanceField();
				else if (meshType == 1 || meshType == 2)	campenVisual1->drawSolidMesh();
				if (showFieldEllipse || showField8Shape) campenVisual1->drawEllipse();
			}
			else if (campenVisual2->showMesh)
			{
				if (AnisVisual::showDistanceField)	campenVisual2->drawDistanceField();
				else if (meshType == 1 || meshType == 2)	campenVisual2->drawSolidMesh();
				if (showFieldEllipse || showField8Shape) campenVisual2->drawEllipse();
			}
		}
	}
}
void PathUtils::setPatchColors(bool randomColor)
{
// 	std::vector<unsigned> colorMap(faceLabelsNum);
// 	std::vector<int> colorInd;
// 	int newColor = 0;
// 	colorInd.push_back(newColor);
// 	colorMap[0] = colorInd[0];
// 	for (unsigned i = 1; i<faceLabelsNum; i++){
// 		std::vector<int> tcolors = colorInd;
// 		for (unsigned j = 0; j<i; j++){
// 			if (cycleConnected[i][j]){
// 				tcolors[colorMap[j]] = -1;
// 			}
// 		}
// 		if (*std::max_element(tcolors.begin(), tcolors.end())<0 || i<4){
// 			newColor++;
// 			colorInd.push_back(newColor);
// 			colorMap[i] = newColor;
// 		}
// 		else{
// 			int randInd = -1;
// 			while (randInd<0){
// 				randInd = tcolors[rand() % tcolors.size()];
// 			}
// 			colorMap[i] = tcolors[randInd];
// 		}
// 	}
	std::vector<unsigned> colorMap(faceLabelsNum);
	for (int i = 0; i < colorMap.size(); i++)colorMap[i] = i;

	std::vector<ColorEngine::color> tcolors;
	unsigned renum = colorMap.size();
	while (renum > 0)
	{
		renum--;

		ColorEngine::color newcolor;
		newcolor.r = 1.0 - double(rand() % 1000) / 2000.0; //range from 0.5-1.0
		newcolor.g = 1.0 - double(rand() % 1000) / 2000.0;
		newcolor.b = 1.0 - double(rand() % 1000) / 2000.0;

		tcolors.push_back(newcolor);
	}

	if (randomColor)
	{
// 		std::vector<ColorEngine::color> colors;
// 		ColorEngine::colorHSV(colorInd.size() * 10, colors);
// 
// 		std::vector<ColorEngine::color> tcolors;
// 		for (unsigned j = 0; j<colorInd.size(); j++){
// 			unsigned trand = rand() % colors.size();
// 			tcolors.push_back(colors[trand]);
// 			colors.erase(colors.begin() + trand);
// 		}

		patchColorsInd = colorMap;
		patchColors = tcolors;
	}
	else
	{
// 		std::vector<ColorEngine::color> tcolors;
// 		ColorEngine::colorHSV(colorInd.size(), tcolors);

		patchColorsInd = colorMap;
		patchColors = tcolors;
	}
	//BL_Patch();
}