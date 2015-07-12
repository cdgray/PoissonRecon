/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior writften permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef MULTI_GRID_OCTREE_DATA_INCLUDED
#define MULTI_GRID_OCTREE_DATA_INCLUDED

#ifdef WIN32
#include <hash_map>
using stdext::hash_map;
#else
#include <ext/hash_map>
using namespace __gnu_cxx;

namespace __gnu_cxx
{
  template<> struct hash<long long> {
    size_t operator()(long long __x) const { return __x; }
  };
  template<> struct hash<const long long> {
    size_t operator()(const long long __x) const { return __x; }
  };
  
  
  template<> struct hash<unsigned long long> {
    size_t operator()(unsigned long long __x) const { return __x; }
  };
  template<> struct hash<const unsigned long long> {
    size_t operator()(const unsigned long long __x) const { return __x; }
  };
}

#endif

typedef float Real;
typedef OctNode<class TreeNodeData,Real> TreeOctNode;

class VertexData{
public:
	static long long EdgeIndex(const TreeOctNode* node,const int& eIndex,const int& maxDepth,int index[DIMENSION]);
	static long long EdgeIndex(const TreeOctNode* node,const int& eIndex,const int& maxDepth);
	static long long FaceIndex(const TreeOctNode* node,const int& fIndex,const int& maxDepth,int index[DIMENSION]);
	static long long FaceIndex(const TreeOctNode* node,const int& fIndex,const int& maxDepth);
	static long long CornerIndex(const int& depth,const int offSet[DIMENSION],const int& cIndex,const int& maxDepth,int index[DIMENSION]);
	static long long CornerIndex(const TreeOctNode* node,const int& cIndex,const int& maxDepth,int index[DIMENSION]);
	static long long CornerIndex(const TreeOctNode* node,const int& cIndex,const int& maxDepth);
	static long long CenterIndex(const int& depth,const int offSet[DIMENSION],const int& maxDepth,int index[DIMENSION]);
	static long long CenterIndex(const TreeOctNode* node,const int& maxDepth,int index[DIMENSION]);
	static long long CenterIndex(const TreeOctNode* node,const int& maxDepth);
};
class SortedTreeNodes{
public:
	TreeOctNode** treeNodes;
	int *nodeCount;
	int maxDepth;
	SortedTreeNodes(void);
	~SortedTreeNodes(void);
	void set(TreeOctNode& root,const int& setIndex);
};
class SortedTreeLeaves{
public:
	TreeOctNode** treeLeaves;
	int leafCount;
	SortedTreeLeaves(void);
	~SortedTreeLeaves(void);
	void set(TreeOctNode& root);
	void set(TreeOctNode& root,const int& maxDepth);
};
class IsoNodeData{
	static int UseAlloc;
public:
	static Allocator<IsoNodeData> Allocator;
	static int UseAllocator(void);
	static void SetAllocator(const int& blockSize);

	IsoNodeData(void);
	~IsoNodeData(void);

	union{
		Real cornerValues[Cube::CORNERS];
		struct{
			int mcIndex;
			int eSegmentCount;
			long long eSegments[2];
		};
	};
	inline int edgeCount(const int& faceIndex) const;
	inline int edgeIndex(const int& faceIndex,const int& e1,const int& e2) const;
	int addEdgeSegment(const int& edgeIndex1,const int& edgeIndex2);
};
class TreeNodeData{
public:
	static int UseIndex;
	union{
		IsoNodeData* isoNode;
		struct{
			int nodeIndex;
			Real centerWeightContribution;
		};
	};
	Real value;

	TreeNodeData(void);
	~TreeNodeData(void);
};

template<int Degree>
class Octree{
	TreeOctNode::NeighborKey neighborKey;	

	Real radius;
	void setNodeIndices(TreeOctNode& tree,int& idx);
	Real GetDotProduct(const int index[DIMENSION]) const;
	Real GetLaplacian(const int index[DIMENSION]) const;
	Real GetDivergence(const int index[DIMENSION],const Point3D<Real>& normal) const;
	void SetDivergence(const int index[DIMENSION],Point3D<Real>& div) const;

	class DivergenceFunction{
	public:
		Point3D<Real> normal;
		Octree<Degree>* ot;
		int index[DIMENSION],scratch[DIMENSION];
		void Function(TreeOctNode* node1,const TreeOctNode* node2);
	};

	class LaplacianProjectionFunction{
	public:
		double value;
		Octree<Degree>* ot;
		int index[DIMENSION],scratch[DIMENSION];
		void Function(TreeOctNode* node1,const TreeOctNode* node2);
	};
	class LaplacianMatrixFunction{
	public:
		int x2,y2,z2,d2;
		Octree<Degree>* ot;
		int index[DIMENSION],scratch[DIMENSION];
		int elementCount,offset;
		MatrixEntry<float>* rowElements;
		int Function(const TreeOctNode* node1,const TreeOctNode* node2);
	};
	class RestrictedLaplacianMatrixFunction{
	public:
		int depth,offset[3];
		Octree<Degree>* ot;
		Real radius;
		int index[DIMENSION],scratch[DIMENSION];
		int elementCount;
		MatrixEntry<float>* rowElements;
		int Function(const TreeOctNode* node1,const TreeOctNode* node2);
	};

	///////////////////////////
	// Evaluation Functions  //
	///////////////////////////
	class PointIndexValueFunction{
	public:
		int res2;
		double* valueTables;
		int index[DIMENSION];
		Real value;
		void Function(const TreeOctNode* node);
	};
	class PointIndexNormalFunction{
	public:
		int res2;
		double* valueTables;
		double* dValueTables;
		Point3D<Real> normal;
		int index[DIMENSION];
		void Function(const TreeOctNode* node);
	};
	class NodeMinMaxValueFunction{
	public:
		Real min,max;
		int maxDepth,res2;
		double* valueTable;
		int cIndex[DIMENSION][2];
		void Function(const TreeOctNode* node1,const TreeOctNode* node2);
	};

	class AdjacencyCountFunction{
	public:
		int adjacencyCount;
		void Function(const TreeOctNode* node1,const TreeOctNode* node2);
	};
	class AdjacencySetFunction{
	public:
		int *adjacencies,adjacencyCount;
		void Function(const TreeOctNode* node1,const TreeOctNode* node2);
	};

	class RefineFunction{
	public:
		int depth;
		void Function(TreeOctNode* node1,const TreeOctNode* node2);
	};
	class FaceEdgesFunction{
	public:
		int fIndex,maxDepth;
		std::vector<CoredEdgeIndex>* edges;
		hash_map<long long,int> *boundaryRoots,*interiorRoots;
		void Function(const TreeOctNode* node1,const TreeOctNode* node2);
	};

	int SolveFixedDepthMatrix(const int& depth,const SortedTreeNodes& sNodes);
	int SolveFixedDepthMatrix(const int& depth,const int& startingDepth,const SortedTreeNodes& sNodes);

	int GetFixedDepthLaplacian(SparseSymmetricMatrix<float>& matrix,const int& depth,const SortedTreeNodes& sNodes);
	int GetRestrictedFixedDepthLaplacian(SparseSymmetricMatrix<float>& matrix,const int* entries,const int& entryCount,const TreeOctNode* rNode,const Real& radius,const SortedTreeNodes& sNodes);

	void SetIsoSurfaceCorners(const Real& isoValue,const int& subdivisionDepth,const int& fullDepthIso);
	static int IsBoundaryFace(const TreeOctNode* node,const int& faceIndex,const int& subdivideDepth);
	static int IsBoundaryEdge(const TreeOctNode* node,const int& edgeIndex,const int& subdivideDepth);
	static int IsBoundaryEdge(const TreeOctNode* node,const int& dir,const int& x,const int& y,const int& subidivideDepth);
	void PreValidate(const Real& isoValue,const int& maxDepth,const int& subdivideDepth);
	void PreValidate(TreeOctNode* node,const Real& isoValue,const int& maxDepth,const int& subdivideDepth);
	void Validate(TreeOctNode* node,const Real& isoValue,const int& maxDepth,const int& fullDepthIso,const int& subdivideDepth);
	void Validate(TreeOctNode* node,const Real& isoValue,const int& maxDepth,const int& fullDepthIso);
	void Subdivide(TreeOctNode* node,const Real& isoValue,const int& maxDepth);

	int SetMCRootPositions(const Real& isoValue,hash_map<long long,int>& roots,hash_map<long long,Point3D<Real> >& normalHash,
		std::vector<Point3D<Real> >& positions,std::vector<Point3D<Real> >* normals,const int& nonLinearFit);
	// Gets the positions of all the iso-vertices
	int SetBoundaryMCRootPositions(const int& sDepth,const Real& isoValue,
		hash_map<long long,int>& boundaryRoots,hash_map<long long,Point3D<Real> >& boundaryNormalHash,CoredMeshData* mesh,const int& nonLinearFit);
	int SetMCRootPositions(TreeOctNode* node,const int& sDepth,const Real& isoValue,
		hash_map<long long,int>& boundaryRoots,hash_map<long long,int>* interiorRoots,
		hash_map<long long,Point3D<Real> >& boundaryNormalHash,hash_map<long long,Point3D<Real> >* interiorNormalHash,
		std::vector<Point3D<float> >* interiorPositions,
		CoredMeshData* mesh,const int& nonLinearFit);
	int SetMCRootPositions(TreeOctNode* node,const int& sDepth,const Real& isoValue,const Real* cornerValues,
		hash_map<long long,int>& boundaryRoots,hash_map<long long,int>* interiorRoots,
		hash_map<long long,Point3D<Real> >& boundaryNormalHash,hash_map<long long,Point3D<Real> >* interiorNormalHash,
		std::vector<Point3D<float> >* interiorPositions,
		CoredMeshData* mesh,const int& nonLinearFit);
	// Adds the triangles associates to the bounaries of the specified leaf node
	int GetMCIsoTriangles(TreeOctNode* node,CoredMeshData* mesh,hash_map<long long,int>& boundaryRoots,
		hash_map<long long,int>* interiorRoots,std::vector<Point3D<float> >* interiorPositions,const int& offSet,const int& sDepth);
	static int AddTriangles(CoredMeshData* mesh,std::vector<int>& edges);
	static int AddTriangles(CoredMeshData* mesh,std::vector<CoredPointIndex> edges[3],std::vector<Point3D<float> >* interiorPositions,const int& offSet);
	static int AddTriangles(CoredMeshData* mesh,std::vector<CoredPointIndex>& edges,std::vector<Point3D<float> >* interiorPositions,const int& offSet);
	static int SetMCFaceCurve(TreeOctNode* node,const int& edgeIndex1,const int& edgeIndex2,const int& maxDepth,std::vector<CoredPointIndex>& curve,
		std::vector<Point3D<float> >& boundaryPositions,hash_map<long long,int>& boundaryRoots,
		std::vector<Point3D<float> >* interiorPositions,const int& offSet,hash_map<long long,int>* interiorRoots,const int& sDepth);

	static int InteriorFaceRootCount(const TreeOctNode* node,const int &faceIndex,const Real& isoValue,const int& maxDepth);
	static int EdgeRootCount(const TreeOctNode* node,const int& edgeIndex,const Real& isoValue,const int& maxDepth);
	int GetRoot(const TreeOctNode* node,const int& edgeIndex,const Real& isoValue,const int& maxDepth,Point3D<Real> & position,hash_map<long long,Point3D<Real> >& normalHash,
		Point3D<Real>* normal,const int& nonLinearFit);
	int GetRoot(const TreeOctNode* node,const int& edgeIndex,const Real* cValues,const Real& isoValue,Point3D<Real> & position,
		hash_map<long long,Point3D<Real> >& normalHash,const int& nonLinearFit);
	static int GetRootIndex(const TreeOctNode* node,const int& edgeIndex,const Real& isoValue,const int& maxDepth,int eIndex[2],int& offset);
	static int GetRootIndex(const TreeOctNode* node,const int& edgeIndex,const int& maxDepth,int eIndex[2],int& offset);
	static int GetRootIndex(const TreeOctNode* node,const int& edgeIndex,const int& maxDepth,const int& sDepth,int eIndex[2],int& offset);
	static long long GetRootIndex(const TreeOctNode* node,const int& edgeIndex,const int& maxDepth);
	static int GetRootIndices(const TreeOctNode* node,const int& edgeIndex,const int& maxDepth,std::vector<long long>& indices);
	static int GetRootIndex(const TreeOctNode* node,const int& localEdgeIndex,const int& maxDepth,
		hash_map<long long,int>& boundaryRoots,hash_map<long long,int>* interiorRoots,CoredPointIndex& index);

	int NonLinearUpdateWeightContribution(TreeOctNode* node,const Point3D<Real>& position);
	Real NonLinearGetSampleWeight(TreeOctNode* node,const Point3D<Real>& position);
	void NonLinearGetSampleDepthAndWeight(TreeOctNode* node,const Point3D<Real>& position,const Real& samplesPerNode,Real& depth,Real& weight);
	int NonLinearSplatOrientedPoint(TreeOctNode* node,const Point3D<Real>& point,const Point3D<Real>& normal);
	void NonLinearSplatOrientedPoint(const Point3D<Real>& point,const Point3D<Real>& normal,const int& kernelDepth,const Real& samplesPerNode,const int& minDepth,const int& maxDepth);

	int HasNormals(TreeOctNode* node,const Real& epsilon);


public:
	static double maxMemoryUsage;
	static double MemoryUsage(void);
	std::vector< Point3D<Real> >* normals;
	Real postNormalSmooth;
	TreeOctNode tree;
	FunctionData<Degree,double> fData;
	Octree(void);

	void setFunctionData(const PPolynomial<Degree>& ReconstructionFunction,const int& maxDepth,const int& normalize,const Real& normalSmooth=-1);
	void finalize1(const int& refineNeighbors=-1);
	void finalize2(const int& refineNeighbors=-1);
	int setTree(char* fileName,const int& maxDepth,const int& binary,const int& kernelDepth,const Real& samplesPerNode,
		const Real& scaleFactor,Point3D<Real>& center,Real& scale,const int& resetSampleDepths=1);

	void SetLaplacianWeights(void);
	void ClipTree(void);
	int LaplacianMatrixIteration(const int& subdivideDepth);

	Real GetIsoValue(void);
	void GetMCIsoTriangles(const Real& isoValue,CoredMeshData* mesh,const int& fullDepthIso=0);
	void GetMCIsoTriangles(const Real& isoValue,const int& subdivideDepth,CoredMeshData* mesh,const int& fullDepthIso=0);
};

#include "MultiGridOctreeData.inl"
#endif // MULTI_GRID_OCTREE_DATA_INCLUDED
