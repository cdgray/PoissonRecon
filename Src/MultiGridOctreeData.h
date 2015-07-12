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

#ifndef MULTI_GRID_OCTREE_DATA_INCLUDED
#define MULTI_GRID_OCTREE_DATA_INCLUDED

#define BROODED_SLICES 1			// Work two slices at a time, so that children remain consecutive in memory.

#define GRADIENT_DOMAIN_SOLUTION 1	// Given the constraint vector-field V(p), there are two ways to solve for the coefficients, x, of the indicator function
									// with respect to the B-spline basis {B_i(p)}
									// 1] Find x minimizing:
									//			|| V(p) - \sum_i \nabla x_i B_i(p) ||^2
									//		which is solved by the system A_1x = b_1 where:
									//			A_1[i,j] = < \nabla B_i(p) , \nabla B_j(p) >
									//			b_1[i]   = < \nabla B_i(p) , V(p) >
									// 2] Formulate this as a Poisson equation:
									//			\sum_i x_i \Delta B_i(p) = \nabla \cdot V(p)
									//		which is solved by the system A_2x = b_2 where:
									//			A_2[i,j] = - < \Delta B_i(p) , B_j(p) >
									//			b_2[i]   = - < B_i(p) , \nabla \cdot V(p) >
									// Although the two system matrices should be the same (assuming that the B_i satisfy dirichlet/neumann boundary conditions)
									// the constraint vectors can differ when V does not satisfy the Neumann boundary conditions:
									//		A_1[i,j] = \int_R < \nabla B_i(p) , \nabla B_j(p) >
									//               = \int_R [ \nabla \cdot ( B_i(p) \nabla B_j(p) ) - B_i(p) \Delta B_j(p) ]
									//               = \int_dR < N(p) , B_i(p) \nabla B_j(p) > + A_2[i,j]
									// and the first integral is zero if either f_i is zero on the boundary dR or the derivative of B_i across the boundary is zero.
									// However, for the constraints we have:
									//		b_1(i)   = \int_R < \nabla B_i(p) , V(p) >
									//               = \int_R [ \nabla \cdot ( B_i(p) V(p) ) - B_i(p) \nabla \cdot V(p) ]
									//               = \int_dR < N(p) ,  B_i(p) V(p) > + b_2[i]
									// In particular, this implies that if the B_i satisfy the Neumann boundary conditions (rather than Dirichlet),
									// and V is not zero across the boundary, then the two constraints are different.
									// Forcing the < V(p) , N(p) > = 0 on the boundary, by killing off the component of the vector-field in the normal direction
									// (FORCE_NEUMANN_FIELD), makes the two systems equal, and the value of this flag should be immaterial.
									// Note that under interpretation 1, we have:
									//		\sum_i b_1(i) = < \nabla \sum_ i B_i(p) , V(p) > = 0
									// because the B_i's sum to one. However, in general, we could have
									//		\sum_i b_2(i) \neq 0.
									// This could cause trouble because the constant functions are in the kernel of the matrix A, so CG will misbehave if the constraint
									// has a non-zero DC term. (Again, forcing < V(p) , N(p) > = 0 along the boundary resolves this problem.)

#define FORCE_NEUMANN_FIELD 1		// This flag forces the normal component across the boundary of the integration domain to be zero.
									// This should be enabled if GRADIENT_DOMAIN_SOLUTION is not, so that CG doesn't run into trouble.

#define ROBERTO_TOLDO_FIX 1

#if !FORCE_NEUMANN_FIELD
#pragma message( "[WARNING] Not zeroing out normal component on boundary" )
#endif // !FORCE_NEUMANN_FIELD

#include "Hash.h"
#include "BSplineData.h"

template< class Real >
class TreeNodeData
{
public:
	static int NodeCount;
	int nodeIndex;

	TreeNodeData( void );
	~TreeNodeData( void );
};

template< class Real >
class RootInfo
{
	typedef OctNode< TreeNodeData< Real > , Real > TreeOctNode;
public:
	const TreeOctNode* node;
	int edgeIndex;
	long long key;
};

template< class Real >
class VertexData
{
	typedef OctNode< TreeNodeData< Real > , Real > TreeOctNode;
public:
	static const int VERTEX_COORDINATE_SHIFT = ( sizeof( long long ) * 8 ) / 3;
	static long long   EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth , int index[DIMENSION] );
	static long long   EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth );
	static long long   FaceIndex( const TreeOctNode* node , int fIndex , int maxDepth,int index[DIMENSION] );
	static long long   FaceIndex( const TreeOctNode* node , int fIndex , int maxDepth );
	static long long CornerIndex( const TreeOctNode* node , int cIndex , int maxDepth , int index[DIMENSION] );
	static long long CornerIndex( const TreeOctNode* node , int cIndex , int maxDepth );
	static long long CenterIndex( const TreeOctNode* node , int maxDepth , int index[DIMENSION] );
	static long long CenterIndex( const TreeOctNode* node , int maxDepth );
	static long long CornerIndex( int depth , const int offSet[DIMENSION] , int cIndex , int maxDepth , int index[DIMENSION] );
	static long long CenterIndex( int depth , const int offSet[DIMENSION] , int maxDepth , int index[DIMENSION] );
	static long long CornerIndexKey( const int index[DIMENSION] );
};
template< class Real >
class SortedTreeNodes
{
	typedef OctNode< TreeNodeData< Real > , Real > TreeOctNode;
protected:
	void _sortByZCoordinate( void );
public:
	Pointer( TreeOctNode* ) treeNodes;
	int *nodeCount;
	int maxDepth;
	SortedTreeNodes( void );
	~SortedTreeNodes( void );
	void set( TreeOctNode& root , std::vector< int >* map );
	Pointer( Pointer( int ) ) sliceOffsets;
	static int Slices( int depth );
	std::pair< int , int > sliceSpan( int depth , int off , int d ) const;
	struct CubeCornerIndices
	{
		int idx[Cube::CORNERS];
		CubeCornerIndices( void ) { memset( idx , -1 , sizeof( int ) * Cube::CORNERS ); }
		int& operator[] ( int i ) { return idx[i]; }
		const int& operator[] ( int i ) const { return idx[i]; }
	};
	struct SlabCornerTableData
	{
		SlabCornerTableData( void ){ cCount=0; }
		~SlabCornerTableData( void ){ clear(); }
		void clear( void ) { cTable.clear() ; cCount = 0; }
		CubeCornerIndices& operator[] ( const TreeOctNode* node );
		const CubeCornerIndices& operator[] ( const TreeOctNode* node ) const;
		CubeCornerIndices& cornerIndices( const TreeOctNode* node );
		const CubeCornerIndices& cornerIndices( const TreeOctNode* node ) const;
		int cCount;
		std::vector< CubeCornerIndices > cTable;
		std::vector< int > offsets;
	};
	void setCornerTable( SlabCornerTableData& cData , int depth , int offset , int maxDepth , int threads ) const;
	void setCornerTable( SlabCornerTableData& cData , int depth , int offset ,                int threads ) const { setCornerTable( cData , depth , offset , maxDepth-1 , threads ); }
	void setCornerTable( SlabCornerTableData& cData ,                                         int threads ) const { setCornerTable( cData , 0     , 0      , maxDepth-1 , threads ); }
	int getMaxSlabCornerCount( int depth , int maxDepth , int threads ) const;
	struct EdgeIndices
	{
		int idx[Cube::EDGES];
		EdgeIndices( void ) { memset( idx , -1 , sizeof( int ) * Cube::EDGES ); }
		int& operator[] ( int i ) { return idx[i]; }
		const int& operator[] ( int i ) const { return idx[i]; }
	};
	struct SlabEdgeTableData
	{
		SlabEdgeTableData( void ) { eCount=0; }
		~SlabEdgeTableData( void ) { clear(); }
		void clear( void ) { eTable.clear() , eCount=0; }
		EdgeIndices& operator[] ( const TreeOctNode* node );
		const EdgeIndices& operator[] ( const TreeOctNode* node ) const;
		EdgeIndices& edgeIndices( const TreeOctNode* node );
		const EdgeIndices& edgeIndices( const TreeOctNode* node ) const;
		int eCount;
		std::vector< EdgeIndices > eTable;
		std::vector< int > offsets;
	};
	void setEdgeTable( SlabEdgeTableData& eData , int depth , int offset , int maxDepth , int threads );
	void setEdgeTable( SlabEdgeTableData& eData , int depth , int offset ,                int threads ) { setEdgeTable( eData , depth , offset , maxDepth-1 , threads ); }
	void setEdgeTable( SlabEdgeTableData& eData ,                                         int threads ) { setEdgeTable( eData , 0     , 0      , maxDepth-1 , threads ); }
	int getMaxSlabEdgeCount( const TreeOctNode* rootNode , int depth , int threads ) const ;
};


template< class Real , int Degree >
class Octree
{
	typedef OctNode< TreeNodeData< Real > , Real > TreeOctNode;
	struct _PointData
	{
		Point3D< Real > position;
		Real weightedCoarserValue;
		Real weight;
		_PointData( Point3D< Real > p=Point3D< Real >() , Real w=0 ) { position = p , weight = w , weightedCoarserValue = Real(0); }
	};
public:
	struct NormalInfo
	{
		std::vector< int > normalIndices;
		std::vector< Point3D< Real > > normals;
		int normalIndex( const TreeOctNode* node ) const { return node->nodeData.nodeIndex>=normalIndices.size() ? -1 : normalIndices[ node->nodeData.nodeIndex ]; }
	};
	struct PointInfo
	{
		std::vector< int > pointIndices;
		std::vector< _PointData > points;
		int pointIndex( const TreeOctNode* node ) const { return node->nodeData.nodeIndex>=pointIndices.size() ? -1 : pointIndices[ node->nodeData.nodeIndex ]; }
	};
protected:
	SortedTreeNodes< Real > _sNodes;
	Real _samplesPerNode;
	int _splatDepth;
	int _minDepth;
	int _fullDepth;
	bool _constrainValues;
	int _boundaryType;
	Real _scale;
	Point3D< Real > _center;
	std::vector< int > _pointCount;
	Real _postDerivativeSmooth;
	BSplineData< Degree > _fData;

	bool _InBounds( Point3D< Real > ) const;

	double GetLaplacian  ( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[3] , const int off2[3] , bool childParent ) const;
	double GetDivergence1( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[3] , const int off2[3] , bool childParent , const Point3D< Real >& normal1 ) const;
	double GetDivergence2( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[3] , const int off2[3] , bool childParent , const Point3D< Real >& normal2 ) const;
	Point3D< double > GetDivergence1( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[3] , const int off2[3] , bool childParent ) const;
	Point3D< double > GetDivergence2( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[3] , const int off2[3] , bool childParent ) const;

	class FaceEdgesFunction
	{
	public:
		int fIndex , maxDepth;
		std::vector< std::pair< RootInfo< Real > , RootInfo< Real > > >* edges;
		hash_map< long long , std::pair< RootInfo< Real > , int > >* vertexCount;
		void Function( const TreeOctNode* node1 , const TreeOctNode* node2 );
		typename TreeOctNode::ConstNeighborKey3* neighborKey3;
		ConstPointer( unsigned char ) mcIndices;
	};
	template< class C , int N > struct Stencil{ C values[N][N][N]; };
	struct CenterValueStencil
	{
		Stencil< double , 3 > stencil;
		Stencil< double , 3 > stencils[8];
	};
	struct CornerValueStencil
	{
		Stencil< double , 3 > stencil[8];
		Stencil< double , 3 > stencils[8][8];
	};
	struct CornerNormalStencil
	{
		Stencil< Point3D< double > , 5 > stencil[8];
		Stencil< Point3D< double > , 5 > stencils[8][8];
	};

	void _setMultiColorIndices( int start , int end , std::vector< std::vector< int > >& indices ) const;
	int _SolveSystemGS( PointInfo& pointInfo , int depth , const typename BSplineData< Degree >::Integrator& integrator , const SortedTreeNodes< Real >& sNodes , Pointer( Real ) solution , Pointer( Real ) constraints , Pointer( Real ) metSolutionConstraints , int iters , bool coarseToFine , bool showResidual=false , double* bNorm2=NULL , double* inRNorm2=NULL , double* outRNorm2=NULL , bool forceSilent=false );
	int _SolveSystemCG( PointInfo& pointInfo , int depth , const typename BSplineData< Degree >::Integrator& integrator , const SortedTreeNodes< Real >& sNodes , Pointer( Real ) solution , Pointer( Real ) constraints , Pointer( Real ) metSolutionConstraints , int iters , bool coarseToFine , bool showResidual=false , double* bNorm2=NULL , double* inRNorm2=NULL , double* outRNorm2=NULL , double accuracy=0 );

	int GetMatrixRowSize( const typename TreeOctNode::Neighbors5& neighbors5 , bool symmetric ) const;
	int SetMatrixRow( const PointInfo& pointInfo , const typename TreeOctNode::Neighbors5& neighbors5 , Pointer( MatrixEntry< Real > ) row , int offset , const typename BSplineData< Degree >::Integrator& integrator , const Stencil< double , 5 >& stencil , bool symmetric ) const;

	void SetDivergenceStencil ( int depth , const typename BSplineData< Degree >::Integrator& integrator , Stencil< Point3D< double > , 5 >& stencil , bool scatter ) const;
	void SetDivergenceStencils( int depth , const typename BSplineData< Degree >::Integrator& integrator , Stencil< Point3D< double > , 5 > stencil[2][2][2] , bool scatter ) const;
	void SetLaplacianStencil  ( int depth , const typename BSplineData< Degree >::Integrator& integrator , Stencil< double , 5 >& stencil ) const;
	void SetLaplacianStencils ( int depth , const typename BSplineData< Degree >::Integrator& integrator , Stencil< double , 5 > stencil[2][2][2] ) const;
	void SetCenterEvaluationStencil ( const typename BSplineData< Degree >::template CenterEvaluator< 1 >& evaluator , int depth , Stencil< double , 3 >& stencil ) const;
	void SetCenterEvaluationStencils( const typename BSplineData< Degree >::template CenterEvaluator< 1 >& evaluator , int depth , Stencil< double , 3 > stencil[8] ) const;
	void SetCornerEvaluationStencil ( const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< double , 3 > stencil [8]    ) const;
	void SetCornerEvaluationStencils( const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< double , 3 > stencils[8][8] ) const;
	void SetCornerNormalEvaluationStencil ( const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 5 > stencil [8]    ) const;
	void SetCornerNormalEvaluationStencils( const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 5 > stencils[8][8] ) const;

	static void UpdateCoarserSupportBounds( const TreeOctNode* node , int& startX , int& endX , int& startY , int& endY , int& startZ , int& endZ );

	void UpdateConstraintsFromCoarser( const PointInfo& pointInfo , const typename TreeOctNode::Neighbors5& neighbors5 , const typename TreeOctNode::Neighbors5& pNeighbors5 , TreeOctNode* node , Pointer( Real ) constraints , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::Integrator& integrator , const Stencil< double , 5 >& stencil ) const;
	// Updates the constraints @(depth-1) based on the solution coefficients @(depth)
	void UpdateConstraintsFromFiner( const typename BSplineData< Degree >::Integrator& integrator , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real ) fineSolution , Pointer( Real ) coarseConstraints ) const;
	// Evaluate the points @(depth) using coefficients @(depth-1)
	void SetPointValuesFromCoarser( PointInfo& pointInfo , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real ) coarseCoefficients );
	// Evalutes the solution @(depth) at the points @(depth-1) and updates the met constraints @(depth-1)
	void SetPointConstraintsFromFiner( const PointInfo& pointInfo , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real )  finerCoefficients , Pointer( Real ) metConstraints) const;
	Real _WeightedCoarserFunctionValue( const _PointData& pointData , const typename TreeOctNode::NeighborKey3& neighborKey3 , const TreeOctNode* node , ConstPointer( Real ) coarseCoefficients ) const;
	Real _WeightedFinerFunctionValue  ( const _PointData& pointData , const typename TreeOctNode::NeighborKey3& neighborKey3 , const TreeOctNode* node , ConstPointer( Real )  finerCoefficients ) const;
	// Down samples constraints @(depth) to constraints @(depth-1)
	template< class C > void DownSample( int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( C ) fineConstraints    , Pointer( C ) coarseConstraints ) const;
	// Up samples solution @(depth-1) to solution @(depth)
	template< class C > void UpSample  ( int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( C ) coarseCoefficients , Pointer( C )  fineCoefficients ) const;
	int GetSliceMatrixAndUpdateConstraints( const PointInfo& pointInfo , SparseMatrix< Real >& matrix , Pointer( Real ) constraints , const typename BSplineData< Degree >::Integrator& integrator , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real ) metSolution , bool coarseToFine , int nStart , int nEnd );
	int GetMatrixAndUpdateConstraints( const PointInfo& pointInfo , SparseSymmetricMatrix< Real >& matrix , Pointer( Real ) constraints , const typename BSplineData< Degree >::Integrator& integrator , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real ) metSolution , bool coarseToFine );

	void SetIsoCorners( ConstPointer( Real ) solution , Real isoValue , Pointer( unsigned char ) mdIndices , TreeOctNode* leaf , typename SortedTreeNodes< Real >::SlabCornerTableData& cData , Pointer( char ) valuesSet , Pointer( Real ) values , typename TreeOctNode::ConstNeighborKey3& nKey , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< double , 3 > stencil[8] , const Stencil< double , 3 > stencils[8][8] );
	static bool IsBoundaryFace( const TreeOctNode* node , int faceIndex , int subdivideDepth );
	static bool IsBoundaryEdge( const TreeOctNode* node , int edgeIndex , int subdivideDepth );
	static bool IsBoundaryEdge( const TreeOctNode* node , int dir , int x , int y , int subidivideDepth );
	static bool IsBoundaryCorner( const TreeOctNode* node , int cornerIndex , int subdivideDepth );

	// For computing the iso-surface there is a lot of re-computation of information across shared geometry.
	// For function values we don't care so much.
	// For edges we need to be careful so that the mesh remains water-tight
	struct RootData : public SortedTreeNodes< Real >::SlabCornerTableData , public SortedTreeNodes< Real >::SlabEdgeTableData
	{
		// Edge to iso-vertex map
		hash_map< long long , int > boundaryRoots;
		// Vertex to ( value , normal ) map
		hash_map< long long , std::pair< Real , Point3D< Real > > > *boundaryValues;
		Pointer( int ) interiorRoots;
		Pointer( Real ) cornerValues;
		Pointer( Point3D< Real > ) cornerNormals;
		Pointer( char ) cornerValuesSet;
		Pointer( char ) cornerNormalsSet;
		Pointer( char ) edgesSet;
	};

	template< class Vertex >
	int SetMCRootPositions( const std::vector< Real >* kernelDensityWeights , ConstPointer( unsigned char ) mcIndices , TreeOctNode* node , int sDepth , Real isoValue , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , RootData& rootData ,
		std::vector< Vertex >* interiorVertices , CoredMeshData< Vertex >* mesh , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< Point3D< double > , 5 > stencil[8] , const Stencil< Point3D< double > , 5 > stencils[8][8] , int nonLinearFit );
	template< class Vertex >
	int GetMCIsoTriangles( ConstPointer( unsigned char ) mcIndices , TreeOctNode* node , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , CoredMeshData< Vertex >* mesh , RootData& rootData ,
		std::vector< Vertex >* interiorVertices , int offSet , int sDepth , bool polygonMesh , std::vector< Vertex >* barycenters );
	template< class Vertex >
	static int AddTriangles( CoredMeshData< Vertex >* mesh , std::vector< CoredPointIndex >& edges , std::vector< Vertex >* interiorVertices , int offSet , bool polygonMesh , std::vector< Vertex >* barycenters );

	void GetMCIsoEdges( ConstPointer( unsigned char ) mcIndices , TreeOctNode* node , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , int sDepth , std::vector< std::pair< RootInfo< Real > , RootInfo< Real > > >& edges );
	static int GetEdgeLoops( std::vector< std::pair< RootInfo< Real > , RootInfo< Real > > >& edges , std::vector< std::vector< std::pair< RootInfo< Real > , RootInfo< Real > > > >& loops);
	static void GetRootSpan( const RootInfo< Real >& ri , Point3D< Real >& start , Point3D< Real >& end );
	template< class Vertex >
	int GetRoot( const std::vector< Real >* kernelDensityWeights , const RootInfo< Real >& ri , Real isoValue , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , Vertex& vertex , RootData& rootData , int sDepth , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< Point3D< double > , 5 > nStencil[8] , const Stencil< Point3D< double > , 5 > nStencils[8][8] , int nonLinearFit );
	static int GetRootIndex( ConstPointer( unsigned char ) mcIndices , const TreeOctNode* node , int edgeIndex , int maxDepth , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , RootInfo< Real >& ri );
	static int GetRootIndex( const RootInfo< Real >& ri , RootData& rootData , CoredPointIndex& index );
	static int GetRootPair( ConstPointer( unsigned char ) mcIndices , const RootInfo< Real >& root , int maxDepth , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , RootInfo< Real >& pair );

	int UpdateWeightContribution( std::vector< Real >& kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , Real weight=Real(1.0) );
	Real GetSampleWeight( const std::vector< Real >& kernelDensityWeight , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , int splatDepth );
	Real GetSampleWeight( const std::vector< Real >& kernelDensityWeight , const TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::ConstNeighborKey3& neighborKey );
	void GetSampleDepthAndWeight( const std::vector< Real >& kernelDensityWeight , const TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::ConstNeighborKey3& neighborKey , Real samplesPerNode , Real& depth , Real& weight );
	Real GetSampleWeight( const std::vector< Real >& kernelDensityWeight , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey );
	void GetSampleDepthAndWeight( const std::vector< Real >& kernelDensityWeight , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , Real samplesPerNode , Real& depth , Real& weight );
	int SplatOrientedPoint( const std::vector< Real >& kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& point , const Point3D< Real >& normal , NormalInfo& normalInfo , typename TreeOctNode::NeighborKey3& neighborKey );
	Real SplatOrientedPoint( const std::vector< Real >& kernelDensityWeights , const Point3D<Real>& point , const Point3D<Real>& normal , NormalInfo& normalInfo , typename TreeOctNode::NeighborKey3& neighborKey , int kernelDepth , Real samplesPerNode , int minDepth , int maxDepth );

	int HasNormals( TreeOctNode* node , const NormalInfo& normalInfo );
	Real getCornerValue( const typename TreeOctNode::ConstNeighborKey3& neighborKey3 , const TreeOctNode* node , int corner , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< double , 3 >& stencil , const Stencil< double , 3 > stencils[8] , bool isInterior ) const;
	Point3D< Real > getCornerNormal( const typename TreeOctNode::ConstNeighbors5& neighbors5 , const typename TreeOctNode::ConstNeighbors5& pNeighbors5 , const TreeOctNode* node , int corner , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< Point3D< double > , 5 >& nStencil , const Stencil< Point3D< double > , 5 > nStencils[8] , bool isInterior ) const;
	Real getCenterValue( const typename TreeOctNode::ConstNeighborKey3& neighborKey3 , const TreeOctNode* node , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CenterEvaluator< 1 >& evaluator , const Stencil< double , 3 >& stencil , const Stencil< double , 3 >& pStencil , bool isInterior ) const;
	static bool _IsInset( const TreeOctNode* node );
	static bool _IsInsetSupported( const TreeOctNode* node );

	int refineBoundary( int subdivisionDepth , std::vector< int >* map );
public:
	int threads;
	static double maxMemoryUsage;
	TreeOctNode tree;

	static double MemoryUsage( void );
	Octree( void );

	void MakeComplete( std::vector< int >* map=NULL );
	void Finalize( int subdivisionDepth , std::vector< int >* map=NULL );
	void ClipTree( const NormalInfo& normalInfo );
	Real GetSolutionValue( Point3D< Real > p , const BSplineData< Degree >* fData=NULL ) const;
	Pointer( Real ) GetSolutionGrid( ConstPointer( Real ) solution , int& res , Real isoValue=0.f , int depth=-1 );
	template< class PointReal >
	int SetTree( char* fileName , int minDepth , int maxDepth , int fullDepth , int splatDepth , Real samplesPerNode ,
		Real scaleFactor , bool useConfidence , bool useNormalWeight , Real constraintWeight , int adaptiveExponent ,
		PointInfo& pointInfo ,
		NormalInfo& normalInfo ,
		std::vector< Real >& kernelDensityWeights , std::vector< Real >& centerWeights ,
		int boundaryType=BSplineElements< Degree >::NONE , XForm4x4< Real > xForm=XForm4x4< Real >::Identity );
	Pointer( Real ) SetLaplacianConstraints( const NormalInfo& normalInfo );
	Pointer( Real ) SolveSystem( PointInfo& pointInfo , Pointer( Real ) constraints , bool showResidual , int iters , int maxSolveDepth , int cgDepth=0 , double cgAccuracy=0 );

	Real GetIsoValue( ConstPointer( Real ) solution , const std::vector< Real >& centerWeights );
	template< class Vertex >
	void GetMCIsoTriangles( const std::vector< Real >* kernelDensityWeights , ConstPointer( Real ) solution , Real isoValue , int subdivideDepth , CoredMeshData< Vertex >* mesh , int fullDepthIso=0 , int nonLinearFit=1 , bool addBarycenter=false , bool polygonMesh=false );
};

#include "MultiGridOctreeData.inl"
#include "MultiGridOctreeData.SortedTreeNodes.inl"
#endif // MULTI_GRID_OCTREE_DATA_INCLUDED
