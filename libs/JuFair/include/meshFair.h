#include "HashMap.h"
#include "mymath.h"

class LuMesh{

public:
	LuMesh();
	~LuMesh();

private:
	int sufvernum;
	int suffacenum;
	double* sufver;
	int* sufface;
	double* suffacenorm;
	int sufctredgenum;
	int* sufctredge;
	int* sufmat;
	vector<int> selectVerList;
	bool interpolate;
	bool doFirstSwap;

	typedef vector<int> intvector;


	void gatherInfoForFair(vector<int>& vermark, vector<intvector>& verneighbr,HashMap& ver2edgehash, vector<int>&edgelist);
	void inline JUFairCenter(double ratio, double* oldver, double* newver, double* oldfirstdiffer, double* newfirstdiffer, 
		double* seconddiffer, vector<int>&vermark, vector<intvector>& verneighbr,HashMap& ver2edgehash, vector<int>&edgelist, vector<double>&wlist, double* norms, double* mags);
	void inline computeWeightList(vector<int>&edgelist, vector<double>& wlist);
	void inline computeDiffer(vector<double>& oldval, vector<double>& differ, vector<double>& wlist, vector<intvector>&verneighbr, HashMap& ver2edgehash);
	void inline computeDiffer(double* oldval, double* differ, vector<double>& wlist,vector<intvector>&verneighbr, HashMap& ver2edgehash);
	void inline computeDiffer(double* oldval, double* differ, vector<double>& wlist, vector<int>&vermark,	vector<intvector>&verneighbr, HashMap& ver2edgehash);
	void inline computeRefreshWeightList(vector<intvector>&verneighbr, vector<double>& w2list);
	void inline computeNorm(double* oldval, double* norms, double* mags);

	void gatherEdgeInfo(vector<int>& edgelist,vector<intvector>& edge2facelist,
		HashMap& ver2edgehash,vector<int>& ctredgelist,vector<int>&	nmedgelist,
		vector<int>& normedgelist,vector<double>& edgelenlist);
	//old one, no vermark
	void gatherVerInfo(vector<double>& verlist,vector<double>&verattrlist, vector<int>&edgelist, vector<double>&edgelenlist,
		vector<int>&nmedgelist, vector<int>&normedgelist, vector<int>&ctredgelist );
	//new one, with vermark
	void gatherVerInfo(vector<double>& verlist,vector<double>&verattrlist, vector<int>&vermark, vector<int>&edgelist, vector<double>&edgelenlist,
		vector<int>&nmedgelist, vector<int>&normedgelist, vector<int>&ctredgelist );
	void gatherFaceInfo(vector<int>& facelist);
	void splitOneNMEdge(int segnum, int i, vector<double>& verlist,vector<double>& verattrlist,vector<int>& edgelist,vector<intvector>& edge2facelist,
		vector<int>&nmedgelist, vector<int>& normedgelist,vector<int>&facelist, HashMap& ver2edgehash);
	void splitNMEdge(vector<double>& verlist,vector<double>& verattrlist,vector<int>& edgelist,vector<intvector>& edge2facelist,
		vector<int>&nmedgelist, vector<int>& normedgelist,vector<int>&facelist,		HashMap& ver2edgehash);
	void swapEdge(vector<double>&verlist,vector<int>& edgelist,vector<intvector>& edge2facelist,vector<int>&normedgelist, vector<int>&facelist,HashMap& ver2edgehash);
	bool splitTriangle(vector<double>&verlist, vector<double>&verattrilist, HashMap& ver2edgehash2,
		vector<int>&edgelist, vector<int>&normedgelist, vector<intvector>&edge2facelist,vector<int>&facelist, double alpha);

public:
	void InputData( const std::vector<double> &mver, const intvector &mface, const intvector &ctrmedge, const bool &doswap);
	void OutputData( std::vector<double> &mver, intvector &mface);
	void LiepaRefine(double alpha);
	void JUFair(double ratio, int times);
};