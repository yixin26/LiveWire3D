#pragma once



namespace SP
{

/* Basic mesh data structure which contains a vertex array */
class Mesh
{
public:
	Mesh();
	~Mesh();
	// clears the mesh
	void clear();
	// allocated space for nNumVertices
	void allocateVertices(int nNumVertices);
	// allocated space for nNumFaces
	void allocateFaces(int nNumFaces);
	// boolean flags the size of nNumVertices true if constrained
	void allocateVertexConstraints();
	// edge constraints (i,j) index pairs
	void allocateEdgeConstraints(int nNumConstraints);

	// gets the number of vertices
	int getVertexCount()const;
	// gets the number of faces
	int getFaceCount()const;
	// get Edge constraints count
	int getEdgeConstraintCount()const;

	// pointer to an array of packed "float3" values containing vertex positions
	float *getPositions()const;

	// pointer to an array of packed "float3" values containing vertex normals
	bool  *getConstraints()const;
	int   *getEdgeConstraints()const; // 2* size of numEdgeConstraints, i,j pairs
	float *getNormals()const;
	float *getTangents()const;
	float *getBinormals()const;

	// pointer to an array of packed "int3" values containing indices into the positions and normals
	int* getFaceIndices()const;

	// copy operator
	Mesh& operator=(const Mesh& other);
private:
	int	  mnVertices;
	float *mpPositions;
	bool  *mpConstraints; // true if vertex is constrained (size of bool)
	int   *mpEdgeConstraints; // true if vertex is constrained (size of bool)
	float *mpNormals,*mpTangents,*mpBinormals;
	int   mnFaces;
	int   mnEdgeConstaints;
	int   *mpFaceIndices;
};

/* define a basic curve loops */
class Loop
{
public:
	enum LoopType{ DEFAULT=0x0,HOLE=0x1 };
	Loop();
	// clears the loop
	void clear();
	// allocated space for nNumVertices
	void allocateVertices(int nNumVertices,bool bNormals=false);
	// returns true if has normals
	bool hasNormals();
	// returns the number of vertices
	int getVertexCount()const;
	// pointer to an array of packed "float3" values containing vertex positions
	float *getPositions()const;
	// pointer to an array of packed "float3" values containing vertex normals
	float *getNormals()const;
	// gets the type of the loop
	Loop::LoopType getLoopType()const;
	// sets the loop type
	void setLoopType(const LoopType loopType);
	// gets the loop group
	int getLoopGroup();
	// loop group index
	void setLoopGroup(int nIndex);
	//
	void computeNormalFromVertices(float *pNormal)const;
	//
	void computeAverageNormal(float *pNormal)const;
private:
	int	  mnVertices;
	bool  mbNormals;
	float *mpPositions;
	float *mpNormals;
	LoopType mLoopType;
	int      mnLoopGroup;
};


class SurfaceSmoothSettings
{
public:
	SurfaceSmoothSettings():mnNumSubdivisions(0), // don't add triangle to the mesh
                            mnLaplacianOrder(2),  // BiLaplacian
                            mnLaplaceIterations(1){}
	
	int		mnNumSubdivisions;		///< number of subdivisions after remeshing - (note beyond 2 the matrix can go unstable)
	int		mnLaplacianOrder;		///< 1 - Laplacian, 2 - BiLaplacian, 3 - TriLaplacian, ...
    int     mnLaplaceIterations;    ///< 1 is usually fine 
};

/* smooths a surface given by an input mesh with constrained points
 NOTE: surface should be closed and manifold */
class SurfaceSmoother
{
public:
	static void smoothSurface(const SurfaceSmoothSettings& settings,
                              const Mesh& inputMesh,Mesh& outputMesh);
};


/* settings for building a smooth patch surface from an input patch */
class SmoothPatchSettings
{
public:
	SmoothPatchSettings():mbConstrainNormals(true),mbRemesh(true),mfTension(0.0),mnNumSubdivisions(2),mnNumLaplacianSmooths(10),mnLaplacianOrder(2){}
	
	bool	mbConstrainNormals;	///< true if boundary normals are provided otherwise normals are computed from the input faces
	bool	mbRemesh;				///< true for remeshing (otherwise the same mesh is returned
	float	mfTension;			///< tension parameter for adjusting the surface tension (default 1.0)
	int		mnNumSubdivisions;		///< number of subdivisions after remeshing - (note beyond 2 the matrix can go unstable)
	int		mnNumLaplacianSmooths; ///< number of smoothing iterations between subdivisions (smoothing between subdivisions to better space the triangles)
	int		mnLaplacianOrder;		///< 1 - Laplacian, 2 - BiLaplacian, 3 - TriLaplacian, ...
};

/* builds the smooth patch */
class SmoothPatchBuilder
{
public:
	// build the patch given settings (see above) and an input mesh patch.  Output mesh is then given
	static void buildSmoothPatch(const SmoothPatchSettings& settings, const Mesh& inputMesh,Mesh& outputMesh);
	static void generateCurves(const Mesh& inputMesh,const char *sCurveFileName,const char *sCurveNormalFileName);
};


}; // end of SP