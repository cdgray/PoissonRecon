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
prior written permission. 

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

#ifndef OCT_NODE_INCLUDED
#define OCT_NODE_INCLUDED

#include "Allocator.h"
#include "BinaryNode.h"
#include "MarchingCubes.h"

#define DIMENSION 3

template<class NodeData,class Real=float>
class OctNode
{
private:
	static int UseAlloc;

	class AdjacencyCountFunction{
	public:
		int count;
		void Function(const OctNode<NodeData,Real>* node1,const OctNode<NodeData,Real>* node2);
	};
	template<class NodeAdjacencyFunction>
	void __processNodeFaces(OctNode* node,NodeAdjacencyFunction* F,const int& cIndex1,const int& cIndex2,const int& cIndex3,const int& cIndex4);
	template<class NodeAdjacencyFunction>
	void __processNodeEdges(OctNode* node,NodeAdjacencyFunction* F,const int& cIndex1,const int& cIndex2);
	template<class NodeAdjacencyFunction>
	void __processNodeCorners(OctNode* node,NodeAdjacencyFunction* F,const int& cIndex);
	template<class NodeAdjacencyFunction>
	void __processNodeNodes(OctNode* node,NodeAdjacencyFunction* F);

	template<class TerminatingNodeAdjacencyFunction>
	static void __ProcessNodeAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const Real& cWidth2,TerminatingNodeAdjacencyFunction* F);
	template<class NodeAdjacencyFunction>
	static void __ProcessTerminatingNodeAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const Real& cWidth2,NodeAdjacencyFunction* F);
	template<class PointAdjacencyFunction>
	static void __ProcessPointAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node2,const Real& radius2,const Real& cWidth2,PointAdjacencyFunction* F);
	template<class PointAdjacencyFunction>
	static void __ProcessPointAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,const Real& radius1,OctNode* node2,const Real& radius2,const Real& cWidth2,PointAdjacencyFunction* F);
	template<class NodeAdjacencyFunction>
	static void __ProcessFixedDepthNodeAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const Real& cWidth2,const int& depth,NodeAdjacencyFunction* F);
	template<class NodeAdjacencyFunction>
	static void __ProcessMaxDepthNodeAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const Real& cWidth2,const int& depth,NodeAdjacencyFunction* F);

	// This is made private because the division by two has been pulled out.
	static inline int Overlap(const Point3D<Real>& center1,const Point3D<Real>& center2,const Real& dWidth);
	static inline int Overlap(const Real& c1,const Real& c2,const Real& c3,const Real& dWidth);
	inline int ChildOverlap(const Point3D<Real>& center1,const Real& radius1,const Real& radius2) const;
	inline int ChildOverlap(const Point3D<Real>& center1,const Real& radius2) const;
	inline int ChildOverlap(const Real& dx,const Real& dy,const Real& dz,const Real& radius) const;
	inline static int ChildOverlap(const Real& dx,const Real& dy,const Real& dz,const Real& d,const Real& cRadius2);
	inline int ChildOverlap2(const Point3D<Real>& center1,const Real& radius1,const Real& radius2,const Real& width2) const;

	const OctNode* __faceNeighbor(const int& dir,const int& off) const;
	const OctNode* __edgeNeighbor(const int& o,const int i[2],const int idx[2]) const;
	OctNode* __faceNeighbor(const int& dir,const int& off,const int& forceChildren);
	OctNode* __edgeNeighbor(const int& o,const int i[2],const int idx[2],const int& forceChildren);
public:
	static const int DepthShift,OffsetShift,OffsetShift1,OffsetShift2,OffsetShift3;
	static const int DepthMask,OffsetMask;

	static Allocator<OctNode> Allocator;
	static int UseAllocator(void);
	static void SetAllocator(int blockSize);

	OctNode* parent;
	OctNode* children;
	short d,off[3];
	NodeData nodeData;


	OctNode(void);
	~OctNode(void);
	int initChildren(void);

	void depthAndOffset(int& depth,int offset[3]) const; 
        int depth(void) const;
	static inline void DepthAndOffset(const long long& index,int& depth,int offset[3]);
	static inline void CenterAndWidth(const long long& index,Point3D<Real>& center,Real& width);
	static inline int Depth(const long long& index);
	static inline void Index(const int& depth,const int offset[3],short& d,short off[3]);
	void centerAndWidth(Point3D<Real>& center,Real& width) const;

	int leaves(void) const;
	int maxDepthLeaves(const int& maxDepth) const;
	int nodes(void) const;
	int maxDepth(void) const;

	const OctNode* root(void) const;

	const OctNode* nextLeaf(const OctNode* currentLeaf=NULL) const;
	OctNode* nextLeaf(OctNode* currentLeaf=NULL);
	const OctNode* nextNode(const OctNode* currentNode=NULL) const;
	OctNode* nextNode(OctNode* currentNode=NULL);
	const OctNode* nextBranch(const OctNode* current) const;
	OctNode* nextBranch(OctNode* current);

	void setFullDepth(const int& maxDepth);

	void printLeaves(void) const;
	void printRange(void) const;

	template<class NodeAdjacencyFunction>
	void processNodeFaces(OctNode* node,NodeAdjacencyFunction* F,const int& fIndex,const int& processCurrent=1);
	template<class NodeAdjacencyFunction>
	void processNodeEdges(OctNode* node,NodeAdjacencyFunction* F,const int& eIndex,const int& processCurrent=1);
	template<class NodeAdjacencyFunction>
	void processNodeCorners(OctNode* node,NodeAdjacencyFunction* F,const int& cIndex,const int& processCurrent=1);
	template<class NodeAdjacencyFunction>
	void processNodeNodes(OctNode* node,NodeAdjacencyFunction* F,const int& processCurrent=1);
	
	template<class NodeAdjacencyFunction>
	static void ProcessNodeAdjacentNodes(OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,NodeAdjacencyFunction* F,const int& processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessNodeAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const Real& width2,NodeAdjacencyFunction* F,const int& processCurrent=1);
	template<class TerminatingNodeAdjacencyFunction>
	static void ProcessTerminatingNodeAdjacentNodes(OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,TerminatingNodeAdjacencyFunction* F,const int& processCurrent=1);
	template<class TerminatingNodeAdjacencyFunction>
	static void ProcessTerminatingNodeAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const Real& width2,TerminatingNodeAdjacencyFunction* F,const int& processCurrent=1);
	template<class PointAdjacencyFunction>
	static void ProcessPointAdjacentNodes(const Point3D<Real>& center1,OctNode* node2,const Real& radius2,PointAdjacencyFunction* F,const int& processCurrent=1);
	template<class PointAdjacencyFunction>
	static void ProcessPointAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node2,const Real& radius2,const Real& width2,PointAdjacencyFunction* F,const int& processCurrent=1);
	template<class PointAdjacencyFunction>
	static void ProcessPointAdjacentNodes(const Point3D<Real>& center1,const Real& radius1,OctNode* node2,const Real& radius2,PointAdjacencyFunction* F,const int& processCurrent=1);
	template<class PointAdjacencyFunction>
	static void ProcessPointAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,const Real& radius1,OctNode* node2,const Real& radius2,const Real& width2,PointAdjacencyFunction* F,const int& processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessFixedDepthNodeAdjacentNodes(OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const int& depth,NodeAdjacencyFunction* F,const int& processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessFixedDepthNodeAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const Real& width2,const int& depth,NodeAdjacencyFunction* F,const int& processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessMaxDepthNodeAdjacentNodes(OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const int& depth,NodeAdjacencyFunction* F,const int& processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessMaxDepthNodeAdjacentNodes(const Real& dx,const Real& dy,const Real& dz,OctNode* node1,const Real& radius1,OctNode* node2,const Real& radius2,const Real& width2,const int& depth,NodeAdjacencyFunction* F,const int& processCurrent=1);

	static int CornerIndex(const Point3D<Real>& center,const Point3D<Real> &p);

	OctNode* faceNeighbor(const int& faceIndex,const int& forceChildren=0);
	const OctNode* faceNeighbor(const int& faceIndex) const;
	OctNode* edgeNeighbor(const int& edgeIndex,const int& forceChildren=0);
	const OctNode* edgeNeighbor(const int& edgeIndex) const;
	OctNode* cornerNeighbor(const int& cornerIndex,const int& forceChildren=0);
	const OctNode* cornerNeighbor(const int& cornerIndex) const;

	OctNode* getNearestLeaf(const Point3D<Real>& p);
	const OctNode* getNearestLeaf(const Point3D<Real>& p) const;

	static int CommonEdge(const OctNode* node1,const int& eIndex1,const OctNode* node2,const int& eIndex2);
	static int CompareForwardDepths(const void* v1,const void* v2);
	static int CompareForwardPointerDepths(const void* v1,const void* v2);
	static int CompareBackwardDepths(const void* v1,const void* v2);
	static int CompareBackwardPointerDepths(const void* v1,const void* v2);


	template<class NodeData2>
	OctNode& operator = (const OctNode<NodeData2,Real>& node);

	static inline int Overlap (const int &depth1,const int offSet1[DIMENSION],const Real& multiplier1,const int &depth2,const int offSet2[DIMENSION],const Real& multiplier2);
	static inline int Overlap2(const int &depth1,const int offSet1[DIMENSION],const Real& multiplier1,const int &depth2,const int offSet2[DIMENSION],const Real& multiplier2);

	class Neighbors{
	public:
		OctNode* neighbors[3][3][3];
		Neighbors(void);
		void clear(void);
	};
	class NeighborKey{
		Neighbors* neighbors;
	public:
		NeighborKey(void);
		~NeighborKey(void);

		void set(const int& depth);
		Neighbors& setNeighbors(OctNode* node);
		Neighbors& getNeighbors(OctNode* node);
	};

	int write(const char* fileName) const;
	int write(FILE* fp) const;
	int read(const char* fileName);
	int read(FILE* fp);
};


#include "Octree.inl"

#endif // OCT_NODE_INCLUDED
