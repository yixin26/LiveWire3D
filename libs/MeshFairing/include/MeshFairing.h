/***************************************************************
* Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
**************************************************************/

#ifndef MY_MESHFAIRING
#define MY_MESHFAIRING

#include "myMesh.h"
#include <map>

struct NatParam
{
	NatParam() :Tension(0.0), NumLaplacianSmooths(3), LaplacianOrder(2){}
	float Tension = 0.0;
	unsigned NumLaplacianSmooths = 3;
	unsigned LaplacianOrder = 2;
};
class MeshFairing
{
public:
	MeshFairing(){ m_mesh = NULL; boundarySmoothTimes = 3; }

	~MeshFairing(){ if (m_mesh != NULL) delete m_mesh; }

	void run();

private:

	void smoothBoundary(unsigned t = 3);
	//void NathanSmoothPatch(std::vector<MyMesh::VertexIter>& vers, std::vector<MyMesh::FaceIter>& faces, std::map<unsigned, unsigned>&, NatParam param);
	void JuSmoothPatch(std::vector<MyMesh::VertexIter>& vers, std::vector<MyMesh::EdgeIter>& edges, std::vector<MyMesh::FaceIter>& faces, std::map<unsigned, unsigned>&, int laps = 3);

public:

	MyMesh* m_mesh;

	unsigned patchNum;
	std::vector<unsigned> faceLabels;//corresponding to patches;

	unsigned boundarySmoothTimes;
	NatParam natParam;
};

#endif
