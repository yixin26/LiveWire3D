#include "MeshFairing.h"
//#include "SmoothPatch/SmoothPatch.h"
#include <set>
#include <iostream>
#include "meshFair.h"

void MeshFairing::run()
{
	//move boundary
	smoothBoundary(boundarySmoothTimes);

	std::vector<unsigned> stats(patchNum, 0);
	for (auto i = m_mesh->getFaces().begin(); i != m_mesh->getFaces().end(); i++)
	{
		stats[faceLabels[i->id()]]++;
	}
	std::list<std::vector<MyMesh::FaceIter> > patches(patchNum);
	std::vector<std::list<std::vector<MyMesh::FaceIter> >::iterator> iPatches(patchNum);
	unsigned ind = 0;
	for (auto i = patches.begin(); i != patches.end(); i++, ind++)
	{
		i->reserve(stats[ind]);
		iPatches[ind] = i;
	}
	for (auto i = m_mesh->getFaces().begin(); i != m_mesh->getFaces().end(); i++)
	{
		iPatches[faceLabels[i->id()]]->push_back(i);
	}

	int tid = 0;
	for (auto i = patches.begin(); i != patches.end(); i++, tid++)
	{
		std::vector<MyMesh::VertexIter> patchesVers;
		std::vector<MyMesh::EdgeIter>  patchesBoundaries;
		std::map<unsigned, unsigned> globalToLocalInd;

		std::set<unsigned> vs;
		std::map<unsigned,unsigned> es;
		for (auto j = i->begin(); j != i->end(); j++)
		{
			for (unsigned k = 0; k < 3; k++)
			{
				vs.insert((*j)->vertex_iter(k)->id());
				if (es.find((*j)->edge_iter(k)->id()) == es.end())
					es[(*j)->edge_iter(k)->id()] = 1;
				else
					es[(*j)->edge_iter(k)->id()] = 2;
			}
		}

		
		for (auto j = vs.begin(); j != vs.end(); j++)
		{
			globalToLocalInd[*j] = patchesVers.size();
			patchesVers.push_back(m_mesh->getVIter()[*j]);
		}

		for (auto j = es.begin(); j != es.end(); j++)
		{
			if (j->second == 1)
				patchesBoundaries.push_back(m_mesh->getEIter()[j->first]);
		}


		//move in-patches with boundary constrainted.
		//NathanSmoothPatch(patchesVers, *i, globalToLocalInd, natParam);
 		JuSmoothPatch(patchesVers, patchesBoundaries, *i, globalToLocalInd);
	}
}

void MeshFairing::smoothBoundary(unsigned t)
{
	std::vector<bool> isBoundaryEdge(m_mesh->getEdges().size(), false);
	for (auto i = m_mesh->getEdges().begin(); i != m_mesh->getEdges().end(); i++)
	{
		if (i->manifold())
		{
			isBoundaryEdge[i->id()] = faceLabels[i->face_iter(0)->id()] != faceLabels[i->face_iter(1)->id()];
		}
	}

	std::vector<unsigned> verJoints(m_mesh->getVertices().size(), 0);
	for (auto i = m_mesh->getEdges().begin(); i != m_mesh->getEdges().end(); i++)
	{
		if (isBoundaryEdge[i->id()])
		{
			verJoints[i->vertex_iter(0)->id()]++;
			verJoints[i->vertex_iter(1)->id()]++;
		}
	}

	std::map<unsigned, std::pair<unsigned,unsigned> > verAdjs;
	for (auto i = m_mesh->getEdges().begin(); i != m_mesh->getEdges().end(); i++)
	{
		if (isBoundaryEdge[i->id()])
		{
			unsigned v1 = i->vertex_iter(0)->id();
			unsigned v2 = i->vertex_iter(1)->id();
			if (verJoints[v1] == 2)
			{
				if (verAdjs.find(v1) == verAdjs.end())
				{
					verAdjs[v1] = std::pair < unsigned, unsigned >(v2, 0);
				}
				else
				{
					verAdjs[v1] = std::pair < unsigned, unsigned >(verAdjs[v1].first, v2);
				}
			}
			if (verJoints[v2] == 2)
			{
				if (verAdjs.find(v2) == verAdjs.end())
				{
					verAdjs[v2] = std::pair < unsigned, unsigned >(v1, 0);
				}
				else
				{
					verAdjs[v2] = std::pair < unsigned, unsigned >(verAdjs[v2].first, v1);
				}
			}
		}
	}

	std::map<unsigned, Vec3 > verJointsCoordinate;
	for (unsigned i = 0; i < verJoints.size(); i++)
	{
		if (verJoints[i] > 0)
		{
			verJointsCoordinate[i] = m_mesh->getVIter()[i]->coordinate();
		}
	}

	for (unsigned i = 0; i < t; i++)
	{
		auto prevCoords = verJointsCoordinate;
		for (auto j = verAdjs.begin(); j != verAdjs.end(); j++)
		{
			auto tv = (prevCoords[j->second.first] + prevCoords[j->second.second]) / 2.0;
			verJointsCoordinate[j->first] = (prevCoords[j->first] + tv) / 2.0;
		}
	}

	for (auto i = verJointsCoordinate.begin(); i != verJointsCoordinate.end(); i++)
	{
		m_mesh->getVIter()[i->first]->coordinate() = i->second;
	}
}
/*
void MeshFairing::NathanSmoothPatch(std::vector<MyMesh::VertexIter>& vers, std::vector<MyMesh::FaceIter>& faces, std::map<unsigned, unsigned>& globalToLocalInd, NatParam param)
{
	unsigned numofpoints = vers.size();
	unsigned numoftilingtris = faces.size();

	SP::Mesh mesh, outputMesh;
	mesh.allocateVertices(numofpoints);
	mesh.allocateFaces(numoftilingtris);

	float *pPositions(mesh.getPositions()),
		*pNormals(mesh.getNormals());

	for (unsigned i = 0; i < numofpoints; i++)
	{
		pPositions[3 * i] = (float)vers[i]->coordinate().x;
		pPositions[3 * i + 1] = (float)vers[i]->coordinate().y;
		pPositions[3 * i + 2] = (float)vers[i]->coordinate().z;
// 		pNormals[3 * i] = (float)vers[i]->normal().x;
// 		pNormals[3 * i + 1] = (float)vers[i]->normal().y;
// 		pNormals[3 * i + 2] = (float)vers[i]->normal().z;
	}

	int *pFaceIndices = mesh.getFaceIndices();
	int nTempIndex(0);
	for (unsigned i = 0; i < numoftilingtris; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			pFaceIndices[nTempIndex] = globalToLocalInd[faces[i]->vertex_iter(j)->id()];
			++nTempIndex;
		}
	}

	SP::SmoothPatchSettings settings;
	settings.mbConstrainNormals = false;
	settings.mbRemesh = false;
	settings.mfTension = param.Tension;
	settings.mnNumSubdivisions = 0;
	settings.mnNumLaplacianSmooths = param.NumLaplacianSmooths;
	settings.mnLaplacianOrder = param.LaplacianOrder;
	SP::SmoothPatchBuilder smoothPatchBuilder;

	cout << "tension:" << settings.mfTension << " lapTimes:" << settings.mnNumLaplacianSmooths << " lapOrder:" << settings.mnLaplacianOrder << endl;
	smoothPatchBuilder.buildSmoothPatch(settings, mesh, outputMesh);

	double** outFaces;
	int *outNum;
	double**outPoints;
	float** outNorms;
	int *outPn;
	double *outTris;

	float *pPoints = outputMesh.getPositions();
// 	*outNum = outputMesh.getFaceCount();
// 	outTris = new double[(*outNum) * 9];
// 	int *pInd = outputMesh.getFaceIndices();
// 	for (int i = 0; i < (*outNum); i++)
// 	{
// 		for (int j = 0; j < 3; j++)
// 		{
// 			int pointID = pInd[i * 3 + j];
// 			for (int k = 0; k < 3; k++)
// 			{
// 				outTris[i * 9 + j * 3 + k] = (double)pPoints[pointID * 3 + k];
// 			}
// 		}
// 	}

// 	float *pPointNormals = outputMesh.getNormals();
// 	if (withNormal)
// 	{
// 		float * outNs = new float[numofpoints * 3];
// 		for (int i = 0; i < numofpoints * 3; i++)
// 		{
// 			outNs[i] = pPointNormals[i];
// 		}
// 		*outNorms = outNs;
// 	}

	unsigned ind = 0;
	for (unsigned i = 0; i < vers.size(); i++)
	{
		vers[i]->coordinate() = Vec3(pPoints[ind], pPoints[ind + 1], pPoints[ind + 2]);
		ind += 3;
	}
}
*/

void MeshFairing::JuSmoothPatch(std::vector<MyMesh::VertexIter>& vers, std::vector<MyMesh::EdgeIter>& edges, std::vector<MyMesh::FaceIter>& faces, std::map<unsigned, unsigned>& globalToLocalInd, int laps)
{
	unsigned numofpoints = vers.size();
	unsigned numofctrlEdges = edges.size();
	unsigned numoftilingtris = faces.size();

	std::vector<double> m_verLu(numofpoints * 3);
	for (unsigned i = 0; i < numofpoints; i += 3)
	{
		m_verLu[i] = vers[i]->coordinate().x;
		m_verLu[i + 1] = vers[i]->coordinate().y;
		m_verLu[i + 2] = vers[i]->coordinate().z;
	}

	intvector m_faceLu(numoftilingtris * 3);
	for (unsigned i = 0; i < numoftilingtris; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			m_faceLu[i * 3 + j] = globalToLocalInd[faces[i]->vertex_iter(j)->id()];
		}
	}

	intvector m_ctrmedgeLu(numofctrlEdges * 2);
	for (unsigned i = 0; i < numofctrlEdges; i++)
	{
		m_ctrmedgeLu[i * 2] = globalToLocalInd[edges[i]->vertex_iter(0)->id()];
		m_ctrmedgeLu[i * 2 + 1] = globalToLocalInd[edges[i]->vertex_iter(1)->id()];
	}

	double times = laps;
	double alphan = 1.414;
	double alpha0 = 0;

	double delta = (alphan - alpha0) / (times == 0 ? 1 : times);
	double alpha = alpha0;

	double ratio = 0.5;

	//entry, include vertices, faces, and controalEdges. doSwap is true or false;
	double doSwap = false;
	LuMesh m_meshFair;
	m_meshFair.InputData(m_verLu, m_faceLu, m_ctrmedgeLu, doSwap);

	for (int i = 0; i < times; i++)
	{
		cout << endl << "smooth:" << i << endl;
		m_meshFair.LiepaRefine(alpha);
		m_meshFair.JUFair(ratio, laps);
		alpha += delta;


	}
	m_verLu.clear();
	m_faceLu.clear();
	m_ctrmedgeLu.clear();
	m_meshFair.OutputData(m_verLu, m_faceLu);
}