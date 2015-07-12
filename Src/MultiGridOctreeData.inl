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

#include "Octree.h"
#include "time.h"
#include "MemoryUsage.h"
#include "PointStream.h"
#include "MAT.h"

#define ITERATION_POWER 1.0/3
#define MEMORY_ALLOCATOR_BLOCK_SIZE 1<<12
#define SPLAT_ORDER 2

const double MATRIX_ENTRY_EPSILON = 0;
const double EPSILON              = 1e-6;
const double ROUND_EPS            = 1e-5;



//////////////////
// TreeNodeData //
//////////////////
template< class Real >
int TreeNodeData< Real >::NodeCount = 0;
template< class Real >
TreeNodeData< Real >::TreeNodeData( void ){ nodeIndex = NodeCount++; }
template< class Real >
TreeNodeData< Real >::~TreeNodeData( void ) { }


////////////
// Octree //
////////////
template< class Real , int Degree > double Octree< Real , Degree >::maxMemoryUsage=0;

template< class Real , int Degree >
double Octree< Real , Degree >::MemoryUsage(void)
{
	double mem = double( MemoryInfo::Usage() ) / (1<<20);
	if( mem>maxMemoryUsage ) maxMemoryUsage=mem;
	return mem;
}

template< class Real , int Degree >
Octree< Real , Degree >::Octree( void )
{
	threads = 1;
	_postDerivativeSmooth = 0;
	_constrainValues = false;
}

template< class Real , int Degree >
bool Octree< Real , Degree >::_IsInset( const TreeOctNode* node )
{
	int d , off[3];
	node->depthAndOffset( d , off );
	int res = 1<<d , o = 1<<(d-2);
	return ( off[0]>=o && off[0]<res-o && off[1]>=o && off[1]<res-o && off[2]>=o && off[2]<res-o );
}
template< class Real , int Degree >
bool Octree< Real , Degree >::_IsInsetSupported( const TreeOctNode* node )
{
	int d , off[3];
	node->depthAndOffset( d , off );
	int res = 1<<d , o = (1<<(d-2))-1;
	return ( off[0]>=o && off[0]<res-o && off[1]>=o && off[1]<res-o && off[2]>=o && off[2]<res-o );
}
template< class Real , int Degree >
int Octree< Real , Degree >::SplatOrientedPoint( const std::vector< Real >& kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , const Point3D<Real>& normal , NormalInfo& normalInfo , typename TreeOctNode::NeighborKey3& neighborKey )
{
	double x , dxdy , dxdydz , dx[DIMENSION][SPLAT_ORDER+1];
	double width;
	int off[3];
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node );
	Point3D<Real> center;
	Real w;
	node->centerAndWidth( center , w );
	width=w;
	for( int i=0 ; i<3 ; i++ )
	{
#if SPLAT_ORDER==2
		off[i] = 0;
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125+1.500*x+0.500*x*x;
		x = ( center[i] - position[i] ) / width;
		dx[i][1] = 0.750        -      x*x;

		dx[i][2] = 1. - dx[i][1] - dx[i][0];
#elif SPLAT_ORDER==1
		x = ( position[i] - center[i] ) / width;
		if( x<0 )
		{
			off[i] = 0;
			dx[i][0] = -x;
		}
		else
		{
			off[i] = 1;
			dx[i][0] = 1. - x;
		}
		dx[i][1] = 1. - dx[i][0];
#elif SPLAT_ORDER==0
		off[i] = 1;
		dx[i][0] = 1.;
#else
#     error Splat order not supported
#endif // SPLAT_ORDER
	}
	for( int i=off[0] ; i<=off[0]+SPLAT_ORDER ; i++ ) for( int j=off[1] ; j<=off[1]+SPLAT_ORDER ; j++ )
	{
		dxdy = dx[0][i] * dx[1][j];
		for( int k=off[2] ; k<=off[2]+SPLAT_ORDER ; k++ )
			if( neighbors.neighbors[i][j][k] )
			{
				dxdydz = dxdy * dx[2][k];
				TreeOctNode* _node = neighbors.neighbors[i][j][k];
				if( normalInfo.normalIndices.size()<TreeNodeData< Real >::NodeCount ) normalInfo.normalIndices.resize( TreeNodeData< Real >::NodeCount , -1 );
				int idx = normalInfo.normalIndex( _node );
				if( idx<0 )
				{
					Point3D<Real> n;
					n[0] = n[1] = n[2] = 0;
					idx = normalInfo.normalIndices[ _node->nodeData.nodeIndex ] = (int)normalInfo.normals.size();
					normalInfo.normals.push_back( n );
				}
				normalInfo.normals[idx] += normal * Real( dxdydz );
			}
	}
	return 0;
}
template< class Real , int Degree >
Real Octree< Real , Degree >::SplatOrientedPoint( const std::vector< Real >& kernelDensityWeights , const Point3D<Real>& position , const Point3D<Real>& normal , NormalInfo& normalInfo , typename TreeOctNode::NeighborKey3& neighborKey , int splatDepth , Real samplesPerNode , int minDepth , int maxDepth )
{
	double dx;
	Point3D<Real> n;
	TreeOctNode* temp;
	int cnt=0;
	double width;
	Point3D< Real > myCenter;
	Real myWidth;
	myCenter[0] = myCenter[1] = myCenter[2] = Real(0.5);
	myWidth = Real(1.0);

	temp = &tree;
	while( temp->depth()<splatDepth )
	{
		if( !temp->children )
		{
			fprintf( stderr , "Octree<Degree>::SplatOrientedPoint error\n" );
			return -1;
		}
		int cIndex=TreeOctNode::CornerIndex(myCenter,position);
		temp=&temp->children[cIndex];
		myWidth/=2;
		if(cIndex&1) myCenter[0] += myWidth/2;
		else		 myCenter[0] -= myWidth/2;
		if(cIndex&2) myCenter[1] += myWidth/2;
		else		 myCenter[1] -= myWidth/2;
		if(cIndex&4) myCenter[2] += myWidth/2;
		else		 myCenter[2] -= myWidth/2;
	}
	Real weight , depth;
	GetSampleDepthAndWeight( kernelDensityWeights , temp , position , neighborKey , samplesPerNode , depth , weight );

	if( depth<minDepth ) depth=Real(minDepth);
	if( depth>maxDepth ) depth=Real(maxDepth);
	int topDepth=int(ceil(depth));

	dx = 1.0-(topDepth-depth);
	if( topDepth<=minDepth )
	{
		topDepth=minDepth;
		dx=1;
	}
	else if( topDepth>maxDepth )
	{
		topDepth=maxDepth;
		dx=1;
	}
	while( temp->depth()>topDepth ) temp=temp->parent;
	while( temp->depth()<topDepth )
	{
		if(!temp->children) temp->initChildren();
		int cIndex=TreeOctNode::CornerIndex(myCenter,position);
		temp=&temp->children[cIndex];
		myWidth/=2;
		if(cIndex&1) myCenter[0] += myWidth/2;
		else		 myCenter[0] -= myWidth/2;
		if(cIndex&2) myCenter[1] += myWidth/2;
		else		 myCenter[1] -= myWidth/2;
		if(cIndex&4) myCenter[2] += myWidth/2;
		else		 myCenter[2] -= myWidth/2;
	}
	width = 1.0 / ( 1<<temp->depth() );
	n = normal * weight / Real( pow( width , 3 ) ) * Real( dx );
	SplatOrientedPoint( kernelDensityWeights , temp , position , n , normalInfo , neighborKey );
	if( fabs(1.0-dx) > EPSILON )
	{
		dx = Real(1.0-dx);
		temp = temp->parent;
		width = 1.0 / ( 1<<temp->depth() );

		n = normal * weight / Real( pow( width , 3 ) ) * Real( dx );
		SplatOrientedPoint( kernelDensityWeights , temp , position , n , normalInfo , neighborKey );
	}
	return weight;
}

template< class Real , int Degree >
void Octree< Real , Degree >::GetSampleDepthAndWeight( const std::vector< Real >& kernelDensityWeights , const TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::ConstNeighborKey3& neighborKey , Real samplesPerNode , Real& depth , Real& weight )
{
	const TreeOctNode* temp=node;
	weight = Real(1.0)/GetSampleWeight( kernelDensityWeights , temp , position , neighborKey );
	if( weight>=samplesPerNode ) depth = Real( temp->depth() + log( weight / samplesPerNode ) / log(double(1<<(DIMENSION-1))) );
	else
	{
		Real oldWeight , newWeight;
		oldWeight = newWeight = weight;
		while( newWeight<samplesPerNode && temp->parent )
		{
			temp=temp->parent;
			oldWeight = newWeight;
			newWeight = Real(1.0)/GetSampleWeight( kernelDensityWeights , temp , position, neighborKey );
		}
		depth = Real( temp->depth() + log( newWeight / samplesPerNode ) / log( newWeight / oldWeight ) );
	}
	weight = Real( pow( double(1<<(DIMENSION-1)) , -double(depth) ) );
}
template< class Real , int Degree >
void Octree< Real , Degree >::GetSampleDepthAndWeight( const std::vector< Real >& kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , Real samplesPerNode , Real& depth , Real& weight )
{
	TreeOctNode* temp=node;
	weight = Real(1.0)/GetSampleWeight( kernelDensityWeights , temp , position , neighborKey );
	if( weight>=samplesPerNode ) depth = Real( temp->depth() + log( weight / samplesPerNode ) / log(double(1<<(DIMENSION-1))) );
	else
	{
		Real oldWeight , newWeight;
		oldWeight = newWeight = weight;
		while( newWeight<samplesPerNode && temp->parent )
		{
			temp=temp->parent;
			oldWeight = newWeight;
			newWeight = Real(1.0)/GetSampleWeight( kernelDensityWeights , temp , position, neighborKey );
		}
		depth = Real( temp->depth() + log( newWeight / samplesPerNode ) / log( newWeight / oldWeight ) );
	}
	weight = Real( pow( double(1<<(DIMENSION-1)) , -double(depth) ) );
}
template< class Real , int Degree >
Real Octree< Real , Degree >::GetSampleWeight( const std::vector< Real >& kernelDensityWeights , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , int splatDepth )
{
	Point3D< Real > myCenter;
	Real myWidth;
	myCenter[0] = myCenter[1] = myCenter[2] = Real(0.5);
	myWidth = Real(1.0);

	TreeOctNode* temp = &tree;
	int d = 0;
	while( d<splatDepth )
	{
		if( !temp->children )
		{
			fprintf( stderr , "Octree<Degree>::SplatOrientedPoint error\n" );
			return -1;
		}
		int cIndex = TreeOctNode::CornerIndex( myCenter , position );
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if( cIndex&1 ) myCenter[0] += myWidth/2;
		else 		   myCenter[0] -= myWidth/2;
		if( cIndex&2 ) myCenter[1] += myWidth/2;
		else 		   myCenter[1] -= myWidth/2;
		if( cIndex&4 ) myCenter[2] += myWidth/2;
		else 		   myCenter[2] -= myWidth/2;
		d++;
	}
	return GetSampleWeight( kernelDensityWeights , temp , position , neighborKey );
}
template< class Real , int Degree >
Real Octree< Real , Degree >::GetSampleWeight( const std::vector< Real >& kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey )
{
	Real weight=0;
	double x , dxdy , dx[DIMENSION][3];
	double width;
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node );
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center,w);
	width=w;

	for( int i=0 ; i<DIMENSION ; i++ )
	{
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125 + 1.500*x + 0.500*x*x;
		x = ( center[i] - position[i] ) / width;
		dx[i][1] = 0.750           -       x*x;

		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
	}

	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		dxdy = dx[0][i] * dx[1][j];
		for( int k=0 ; k<3 ; k++ ) if( neighbors.neighbors[i][j][k] )
			weight += Real( dxdy * dx[2][k] * kernelDensityWeights[ neighbors.neighbors[i][j][k]->nodeData.nodeIndex ] );
	}
	return Real( 1.0 / weight );
}
template< class Real , int Degree >
Real Octree< Real , Degree >::GetSampleWeight( const std::vector< Real >& kernelDensityWeights , const TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::ConstNeighborKey3& neighborKey )
{
	Real weight=0;
	double x,dxdy,dx[DIMENSION][3];
	double width;
	typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors( node );
	Point3D<Real> center;
	Real w;
	node->centerAndWidth( center , w );
	width=w;

	for( int i=0 ; i<DIMENSION ; i++ )
	{
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125 + 1.500*x + 0.500*x*x;
		x = ( center[i] - position[i] ) / width;
		dx[i][1] = 0.750           -       x*x;

		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
	}

	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		dxdy = dx[0][i] * dx[1][j];
		for( int k=0 ; k<3 ; k++ ) if( neighbors.neighbors[i][j][k] )
			weight += Real( dxdy * dx[2][k] * kernelDensityWeights[ neighbors.neighbors[i][j][k]->nodeData.nodeIndex ] );
	}
	return Real( 1.0 / weight );
}
template< class Real , int Degree >
int Octree< Real , Degree >::UpdateWeightContribution( std::vector< Real >& kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , Real weight )
{
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node );
	if( kernelDensityWeights.size()<TreeNodeData< Real >::NodeCount ) kernelDensityWeights.resize( TreeNodeData< Real >::NodeCount , 0 );
	double x , dxdy , dx[DIMENSION][3] , width;
	Point3D< Real > center;
	Real w;
	node->centerAndWidth( center , w );
	width=w;
	const double SAMPLE_SCALE = 1. / ( 0.125 * 0.125 + 0.75 * 0.75 + 0.125 * 0.125 );

	for( int i=0 ; i<DIMENSION ; i++ )
	{
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125 + 1.500*x + 0.500*x*x;
		dx[i][1] = -0.25 - 2.*x - x*x;
		dx[i][2] = 1. - dx[i][1] - dx[i][0];
		// Note that we are splatting along a co-dimension one manifold, so uniform point samples
		// do not generate a unit sample weight.
		dx[i][0] *= SAMPLE_SCALE;
	}
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		dxdy = dx[0][i] * dx[1][j] * weight;
		TreeOctNode** _neighbors = neighbors.neighbors[i][j];
		for( int k=0 ; k<3 ; k++ ) if( _neighbors[k] ) kernelDensityWeights[ _neighbors[k]->nodeData.nodeIndex ] += Real( dxdy * dx[2][k] );
	}
	return 0;
}
template< class Real , int Degree >
bool Octree< Real , Degree >::_InBounds( Point3D< Real > p ) const
{
	if( _boundaryType==0 ){ if( p[0]<Real(0.25) || p[0]>Real(0.75) || p[1]<Real(0.25) || p[1]>Real(0.75) || p[2]<Real(0.25) || p[2]>Real(0.75) ) return false; }
	else                  { if( p[0]<Real(0.00) || p[0]>Real(1.00) || p[1]<Real(0.00) || p[1]>Real(1.00) || p[2]<Real(0.00) || p[2]>Real(1.00) ) return false; }
	return true;
}
template< class Real , int Degree >
template< class PointReal >
int Octree< Real , Degree >::SetTree( char* fileName , int minDepth , int maxDepth , int fullDepth , 
							int splatDepth , Real samplesPerNode , Real scaleFactor ,
							bool useConfidence , bool useNormalWeights , Real constraintWeight , int adaptiveExponent ,
							PointInfo& pointInfo ,
							NormalInfo& normalInfo ,
							std::vector< Real >& kernelDensityWeights , std::vector< Real >& centerWeights ,
							int boundaryType , XForm4x4< Real > xForm )
{
	if( splatDepth<0 ) splatDepth = 0;

	_boundaryType = boundaryType;
	if     ( _boundaryType<0 ) _boundaryType = -1;
	else if( _boundaryType>0 ) _boundaryType =  1;
	_samplesPerNode = samplesPerNode;
	_splatDepth = splatDepth;
	_constrainValues = (constraintWeight>0);

	XForm3x3< Real > xFormN;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) xFormN(i,j) = xForm(i,j);
	xFormN = xFormN.transpose().inverse();
	minDepth = std::min< int >( minDepth , maxDepth );	// minDepth <= maxDepth
	fullDepth = std::max< int >( minDepth , std::min< int >( fullDepth , maxDepth ) );	// minDepth <= fullDepth <= maxDepth
	// If _boundaryType==0, points are scaled to be in the [0.25,0.75]^3 cube so all depths have to be offset by
	// and the minDepth has to be 2.
	if( _boundaryType==0 )
	{
		minDepth++ , maxDepth++ , fullDepth++;
		if( splatDepth ) splatDepth++;
		minDepth = std::max< int >( minDepth , 2 );
	}
	// Otherwise the points are in the [0,1]^3 cube.
	// However, for Neumann constraints, the function at depth 0 is constant so the system matrix is zero if there
	// is no screening.
	else if( _boundaryType==1 && !_constrainValues ) minDepth = std::max< int >( minDepth , 1 );

	_postDerivativeSmooth = Real(1.0)/(1<<maxDepth);
	_fData.set( maxDepth , _boundaryType );

	_minDepth = minDepth;
	_fullDepth = fullDepth;
	double pointWeightSum = 0;
	Point3D< Real > min , max , myCenter;
	Real myWidth;
	int i , cnt=0;
	TreeOctNode* temp;

	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set( maxDepth );
	PointStream< PointReal >* pointStream;
	char* ext = GetFileExtension( fileName );
	if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryPointStream< PointReal >( fileName );
	else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYPointStream< PointReal >( fileName );
	else                                    pointStream = new  ASCIIPointStream< PointReal >( fileName );
	delete[] ext;

	tree.setFullDepth( _fullDepth );

	// Read through once to get the center and scale
	{
		double t = Time();
		Point3D< Real > p;
		Point3D< PointReal > _p , _n;
		while( pointStream->nextPoint( _p , _n ) )
		{
			p = xForm * Point3D< Real >(_p);
			for( i=0 ; i<DIMENSION ; i++ )
			{
				if( !cnt || p[i]<min[i] ) min[i] = p[i];
				if( !cnt || p[i]>max[i] ) max[i] = p[i];
			}
			cnt++;
		}

		if( _boundaryType==0 ) _scale = std::max< Real >( max[0]-min[0] , std::max< Real >( max[1]-min[1] , max[2]-min[2] ) ) * 2;
		else         _scale = std::max< Real >( max[0]-min[0] , std::max< Real >( max[1]-min[1] , max[2]-min[2] ) );
		_center = ( max+min ) /2;
	}

	_scale *= scaleFactor;
	for( i=0 ; i<DIMENSION ; i++ ) _center[i] -= _scale/2;
	if( splatDepth>0 )
	{
		double t = Time();
		cnt = 0;
		pointStream->reset();
		Point3D< Real > p , n;
		Point3D< PointReal > _p , _n;
		while( pointStream->nextPoint( _p , _n ) )
		{
			p = xForm * Point3D< Real >(_p) , n = xFormN * Point3D< Real >(_n);
			p = ( p - _center ) / _scale;
			if( !_InBounds(p) ) continue;
			myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
			myWidth = Real(1.0);
			Real weight=Real( 1. );
			if( useConfidence ) weight = Real( Length(n) );
			temp = &tree;
			int d=0;
			while( d<splatDepth )
			{
				UpdateWeightContribution( kernelDensityWeights , temp , p , neighborKey , weight );
				if( !temp->children ) temp->initChildren();
				int cIndex=TreeOctNode::CornerIndex( myCenter , p );
				temp = temp->children + cIndex;
				myWidth/=2;
				if( cIndex&1 ) myCenter[0] += myWidth/2;
				else           myCenter[0] -= myWidth/2;
				if( cIndex&2 ) myCenter[1] += myWidth/2;
				else           myCenter[1] -= myWidth/2;
				if( cIndex&4 ) myCenter[2] += myWidth/2;
				else           myCenter[2] -= myWidth/2;
				d++;
			}
			UpdateWeightContribution( kernelDensityWeights , temp , p , neighborKey , weight );
			cnt++;
		}
	}
	kernelDensityWeights.resize( TreeNodeData< Real >::NodeCount , 0 );

	std::vector< _PointData >& points = pointInfo.points;

	cnt = 0;
	pointStream->reset();
	Point3D< Real > p , n;
	Point3D< PointReal > _p , _n;
	while( pointStream->nextPoint( _p , _n ) )
	{
		p = xForm * Point3D< Real >(_p) , n = xFormN * Point3D< Real >(_n);
		n *= Real(-1.);
		p = ( p - _center ) / _scale;
		if( !_InBounds(p) ) continue;
		myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
		myWidth = Real(1.0);
		Real normalLength = Real( Length( n ) );
		if( normalLength!=normalLength || normalLength<=EPSILON ) continue;
		if( !useConfidence ) n /= normalLength;

		Real pointWeight = Real(1.f);
		if( samplesPerNode>0 && splatDepth ) pointWeight = SplatOrientedPoint( kernelDensityWeights , p , n , normalInfo , neighborKey , splatDepth , samplesPerNode , _minDepth , maxDepth );
		else
		{
			temp = &tree;
			int d=0;
			if( splatDepth )
			{
				while( d<splatDepth )
				{
					int cIndex=TreeOctNode::CornerIndex(myCenter,p);
					temp = &temp->children[cIndex];
					myWidth /= 2;
					if(cIndex&1) myCenter[0] += myWidth/2;
					else		 myCenter[0] -= myWidth/2;
					if(cIndex&2) myCenter[1] += myWidth/2;
					else		 myCenter[1] -= myWidth/2;
					if(cIndex&4) myCenter[2] += myWidth/2;
					else		 myCenter[2] -= myWidth/2;
					d++;
				}
				pointWeight = GetSampleWeight( kernelDensityWeights , temp , p , neighborKey );
			}
			for( i=0 ; i<DIMENSION ; i++ ) n[i] *= pointWeight;
			while( d<maxDepth )
			{
				if( !temp->children ) temp->initChildren();
				int cIndex=TreeOctNode::CornerIndex(myCenter,p);
				temp=&temp->children[cIndex];
				myWidth/=2;
				if(cIndex&1) myCenter[0] += myWidth/2;
				else		 myCenter[0] -= myWidth/2;
				if(cIndex&2) myCenter[1] += myWidth/2;
				else		 myCenter[1] -= myWidth/2;
				if(cIndex&4) myCenter[2] += myWidth/2;
				else		 myCenter[2] -= myWidth/2;
				d++;
			}
			SplatOrientedPoint( kernelDensityWeights , temp , p , n , normalInfo , neighborKey );
		}
		pointWeightSum += pointWeight;
		if( _constrainValues )
		{
			Real pointScreeningWeight = useNormalWeights ? Real( normalLength ) : Real(1.f);
			int d = 0;
			TreeOctNode* temp = &tree;
			myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
			myWidth = Real(1.0);
			while( 1 )
			{
				if( pointInfo.pointIndices.size()<TreeNodeData< Real >::NodeCount ) pointInfo.pointIndices.resize( TreeNodeData< Real >::NodeCount , -1 );
				int idx = pointInfo.pointIndex( temp );

				if( idx==-1 )
				{
					idx = (int)points.size();
					points.push_back( _PointData( p*pointScreeningWeight , pointScreeningWeight ) );
					pointInfo.pointIndices[ temp->nodeData.nodeIndex ] = idx;
				}
				else
				{
					points[idx].weight += pointScreeningWeight;
					points[idx].position += p*pointScreeningWeight;
				}

				int cIndex = TreeOctNode::CornerIndex( myCenter , p );
				if( !temp->children ) break;
				temp = &temp->children[cIndex];
				myWidth /= 2;
				if( cIndex&1 ) myCenter[0] += myWidth/2;
				else		   myCenter[0] -= myWidth/2;
				if( cIndex&2 ) myCenter[1] += myWidth/2;
				else		   myCenter[1] -= myWidth/2;
				if( cIndex&4 ) myCenter[2] += myWidth/2;
				else		   myCenter[2] -= myWidth/2;
				d++;
			}
		}
		cnt++;
	}

	if( _boundaryType==0 ) pointWeightSum *= Real(4.);
	constraintWeight *= Real( pointWeightSum );
	constraintWeight /= cnt;

	MemoryUsage( );
	delete pointStream;
	if( _constrainValues )
		// Set the average position and scale the weights
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode(node) )
			if( pointInfo.pointIndex( node )!=-1 )
			{
				int idx = pointInfo.pointIndex( node );
				points[idx].position /= points[idx].weight;
				int e = ( _boundaryType==0 ? node->depth()-1 : node->depth() ) * adaptiveExponent - ( _boundaryType==0 ? maxDepth-1 : maxDepth ) * (adaptiveExponent-1);
				if( e<0 ) points[idx].weight /= Real( 1<<(-e) );
				else      points[idx].weight *= Real( 1<<  e  );
				points[idx].weight *= Real( constraintWeight );
			}
#if FORCE_NEUMANN_FIELD
	if( _boundaryType==1 )
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode( node ) )
		{
			int d , off[3] , res;
			node->depthAndOffset( d , off );
			res = 1<<d;
			int idx = normalInfo.normalIndex( node );
			if( idx<0 ) continue;
			Point3D< Real >& normal = normalInfo.normals[ idx ];
			for( int d=0 ; d<3 ; d++ ) if( off[d]==0 || off[d]==res-1 ) normal[d] = 0;
		}
#endif // FORCE_NEUMANN_FIELD
	centerWeights.resize( tree.nodes() , 0 );
	kernelDensityWeights.resize( tree.nodes() , 0 );
	// Set the point weights for evaluating the iso-value
	for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode(node) )
	{
		int idx = normalInfo.normalIndex( node );
		if( idx<0 ) centerWeights[ node->nodeData.nodeIndex ] = 0;
		else        centerWeights[ node->nodeData.nodeIndex ] = Real( Length( normalInfo.normals[ idx ] ) );
	}
	MemoryUsage();
	return cnt;
}
template< class Real , int Degree >
void Octree< Real , Degree >::MakeComplete( std::vector< int >* map )
{
	tree.setFullDepth( tree.maxDepth() );
	refineBoundary( 0 , map );
	MemoryUsage();
}
template< class Real , int Degree >
void Octree< Real , Degree >::ClipTree( const NormalInfo& normalInfo )
{
	int maxDepth = tree.maxDepth();
	for( TreeOctNode* temp=tree.nextNode() ; temp ; temp=tree.nextNode(temp) )
		if( temp->children && temp->depth()>=_fullDepth )
		{
			int hasNormals=0;
			for( int i=0 ; i<Cube::CORNERS && !hasNormals ; i++ ) hasNormals = HasNormals( &temp->children[i] , normalInfo );
			if( !hasNormals ) temp->children=NULL;
		}
	MemoryUsage();
}

template< class Real , int Degree >
void Octree< Real , Degree >::Finalize( int subdivideDepth , std::vector< int >* map )
{
	int maxDepth = tree.maxDepth( );
	typename TreeOctNode::NeighborKey3 nKey;
	nKey.set( maxDepth );
	for( int d=maxDepth ; d>1 ; d-- )
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode( node ) ) if( node->depth()==d )
		{
			typename TreeOctNode::Neighbors3& neighbors = nKey.setNeighbors( node->parent->parent );
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( neighbors.neighbors[i][j][k] && !neighbors.neighbors[i][j][k]->children )
					neighbors.neighbors[i][j][k]->initChildren();
		}
	refineBoundary( subdivideDepth , map );
}
template< class Real , int Degree >
double Octree< Real , Degree >::GetLaplacian( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent ) const
{
	double vv[] =
	{
		integrator.dot( d , off1[0] , off2[0] , false , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , false , childParent )
	};
	double dd[] =
	{
		integrator.dot( d , off1[0] , off2[0] , true , true , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , true , true , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , true , true , childParent )
	};
	return dd[0]*vv[1]*vv[2] + vv[0]*dd[1]*vv[2] + vv[0]*vv[1]*dd[2];
}
template< class Real , int Degree >
double Octree< Real , Degree >::GetDivergence1( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent , const Point3D< Real >& normal1 ) const
{
	return Point3D< double >::Dot( GetDivergence1( integrator , d , off1 , off2 , childParent ) , normal1 );
}
template< class Real , int Degree > 
double Octree< Real , Degree >::GetDivergence2( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent , const Point3D< Real >& normal2 ) const
{
	return Point3D< double >::Dot( GetDivergence2( integrator , d , off1 , off2 , childParent ) , normal2 );
}
template< class Real , int Degree >
Point3D< double > Octree< Real , Degree >::GetDivergence1( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent ) const
{
	double vv[] =
	{
		integrator.dot( d , off1[0] , off2[0] , false , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , false , childParent )
	};
#if GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the vector-field with the gradient of the basis function
	double vd[] = 
	{
		integrator.dot( d , off1[0] , off2[0] , false , true , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , true , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , true , childParent )
	};
	return  Point3D< double >( vd[0]*vv[1]*vv[2] , vv[0]*vd[1]*vv[2] , vv[0]*vv[1]*vd[2] );
#else // !GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the divergence of the vector-field with the basis function
	double dv[] = 
	{
		integrator.dot( d , off1[0] , off2[0] , true , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , true , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , true , false , childParent )
	};
	return  -Point3D< double >( dv[0]*vv[1]*vv[2] , vv[0]*dv[1]*vv[2] , vv[0]*vv[1]*dv[2] );
#endif // GRADIENT_DOMAIN_SOLUTION
}
template< class Real , int Degree > 
Point3D< double > Octree< Real , Degree >::GetDivergence2( const typename BSplineData< Degree >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent ) const
{
	double vv[] =
	{
		integrator.dot( d , off1[0] , off2[0] , false , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , false , childParent )
	};
#if GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the vector-field with the gradient of the basis function
	double dv[] = 
	{
		integrator.dot( d , off1[0] , off2[0] , true , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , true , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , true , false , childParent )
	};
	return  Point3D< double >( dv[0]*vv[1]*vv[2] , vv[0]*dv[1]*vv[2] , vv[0]*vv[1]*dv[2] );
#else // !GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the divergence of the vector-field with the basis function
	double vd[] = 
	{
		integrator.dot( d , off1[0] , off2[0] , false , true , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , true , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , true , childParent )
	};
	return -Point3D< double >( vd[0]*vv[1]*vv[2] , vv[0]*vd[1]*vv[2] , vv[0]*vv[1]*vd[2] );
#endif // GRADIENT_DOMAIN_SOLUTION
}

template< class Real , int Degree >
int Octree< Real , Degree >::GetMatrixRowSize( const typename TreeOctNode::Neighbors5& neighbors5 , bool symmetric ) const
{
	int count = 0;
	int nodeIndex = neighbors5.neighbors[2][2][2]->nodeData.nodeIndex;
	const TreeOctNode* const * _nodes = &neighbors5.neighbors[0][0][0];
	if( symmetric )
	{
		for( int i=0 ; i<125 ; i++ ) if( _nodes[i] && _nodes[i]->nodeData.nodeIndex>=nodeIndex ) count++;
	}
	else
	{
		for( int i=0 ; i<125 ; i++ ) if( _nodes[i] ) count++;
	}
	return count;
}

template< class Real , int Degree >
int Octree< Real , Degree >::SetMatrixRow( const PointInfo& pointInfo , const typename TreeOctNode::Neighbors5& neighbors5 , Pointer( MatrixEntry< Real > ) row , int offset , const typename BSplineData< Degree >::Integrator& integrator , const Stencil< double , 5 >& stencil , bool symmetric ) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	bool hasYZPoints[3] , hasZPoints[3][3];
	Real diagonal = 0;
	Real splineValues[3*3*3*3*3];
	memset( splineValues , 0 , sizeof( Real ) * 3 * 3 * 3 * 3 * 3 );

	int count = 0;
	const TreeOctNode* node = neighbors5.neighbors[2][2][2];

	bool isInterior;
	int d , off[3];
	node->depthAndOffset( d , off );

	int o = _boundaryType==0 ? ( 1<<(d-2) ) : 0;
	int mn = 2+o , mx = (1<<d)-2-o;
	isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );


	if( _constrainValues )
	{
		int idx[3] ; node->centerIndex( idx );
		for( int j=0 ; j<3 ; j++ )
		{
			hasYZPoints[j] = false;
			for( int k=0 ; k<3 ; k++ )
			{
				hasZPoints[j][k] = false;
				for( int l=0 ; l<3 ; l++ )
				{
					const TreeOctNode* _node = neighbors5.neighbors[j+1][k+1][l+1];
					if( _node && pointInfo.pointIndex( _node )!=-1 )
					{
						const _PointData& pData = points[ pointInfo.pointIndex( _node ) ];
						Real* _splineValues = splineValues + 3*3*(3*(3*j+k)+l);
						Real weight = pData.weight;
						Point3D< Real > p = pData.position;
						for( int s=0 ; s<3 ; s++ )
						{
#if ROBERTO_TOLDO_FIX
							if( idx[0]+j-s>=0 && idx[0]+j-s<((2<<d)-1) ) _splineValues[3*0+s] = Real( _fData.baseBSplines[ idx[0]+j-s][s]( p[0] ) );
							if( idx[1]+k-s>=0 && idx[1]+k-s<((2<<d)-1) ) _splineValues[3*1+s] = Real( _fData.baseBSplines[ idx[1]+k-s][s]( p[1] ) );
							if( idx[2]+l-s>=0 && idx[2]+l-s<((2<<d)-1) ) _splineValues[3*2+s] = Real( _fData.baseBSplines[ idx[2]+l-s][s]( p[2] ) );
#else // !ROBERTO_TOLDO_FIX
							_splineValues[3*0+s] = Real( _fData.baseBSplines[ idx[0]+j-s][s]( p[0] ) );
							_splineValues[3*1+s] = Real( _fData.baseBSplines[ idx[1]+k-s][s]( p[1] ) );
							_splineValues[3*2+s] = Real( _fData.baseBSplines[ idx[2]+l-s][s]( p[2] ) );
#endif // ROBERTO_TOLDO_FIX
						}
						Real value = _splineValues[3*0+j] * _splineValues[3*1+k] * _splineValues[3*2+l];
						Real weightedValue = value * weight;
						for( int s=0 ; s<3 ; s++ ) _splineValues[3*0+s] *= weightedValue;
						diagonal += value * value * weight;
						hasYZPoints[j] = hasZPoints[j][k] = true;
					}
				}
			}
		}
	}

	Real pointValues[5][5][5];
	if( _constrainValues )
	{
		memset( pointValues , 0 , sizeof(Real)*5*5*5 );
		for( int i=0 ; i<3 ; i++ ) if( hasYZPoints[i] )
			for( int j=0 ; j<3 ; j++ ) if( hasZPoints[i][j] )
				for( int k=0 ; k<3 ; k++ )
				{
					const TreeOctNode* _node = neighbors5.neighbors[i+1][j+1][k+1];
					if( _node && pointInfo.pointIndex( _node )!=-1 )
						{
							const Real* _splineValuesX = splineValues + 3*(3*(3*(3*i+j)+k)+0)+2;
							const Real* _splineValuesY = splineValues + 3*(3*(3*(3*i+j)+k)+1)+2;
							const Real* _splineValuesZ = splineValues + 3*(3*(3*(3*i+j)+k)+2)+2;
							for( int ii=0 ; ii<=2 ; ii++ )
							{
								Real splineValue = _splineValuesX[-ii];
								for( int jj=0 ; jj<=2 ; jj++ )
								{
									Real* _pointValues = pointValues[i+ii][j+jj]+k;
									Real _splineValue = splineValue * _splineValuesY[-jj];
									for( int kk=0 ; kk<=2 ; kk++ ) _pointValues[kk] += _splineValue * _splineValuesZ[-kk];
								}
							}
						}
				}
	}

	pointValues[2][2][2] = diagonal;
	int nodeIndex = neighbors5.neighbors[2][2][2]->nodeData.nodeIndex;
	if( isInterior ) // General case, so try to make fast
	{
		const TreeOctNode* const * _nodes = &neighbors5.neighbors[0][0][0];
		const double* _stencil = &stencil.values[0][0][0];
		Real* _values = &pointValues[0][0][0];
		if( _constrainValues ) for( int i=0 ; i<125 ; i++ ) _values[i] = Real( _stencil[i] + _values[i] );
		else                   for( int i=0 ; i<125 ; i++ ) _values[i] = Real( _stencil[i] );
		if( symmetric ) pointValues[2][2][2] /= 2;
		row[count++] = MatrixEntry< Real >( nodeIndex-offset , _values[5*5*2+5*2+2] );
		if( symmetric )
		{
			for( int i=0 ; i<125 ; i++ ) if( i!=(5*5*2+5*2+2) && _nodes[i] && _nodes[i]->nodeData.nodeIndex>=nodeIndex )
				row[count++] = MatrixEntry< Real >( _nodes[i]->nodeData.nodeIndex-offset , _values[i] );
		}
		else
		{
			for( int i=0 ; i<125 ; i++ ) if( i!=(5*5*2+5*2+2) && _nodes[i] )
				row[count++] = MatrixEntry< Real >( _nodes[i]->nodeData.nodeIndex-offset , _values[i] );
		}
	}
	else
	{
		int d , off[3];
		node->depthAndOffset( d , off );
		Real temp = Real( GetLaplacian( integrator , d , off , off , false ) );
		if( _constrainValues ) temp += pointValues[2][2][2];
		if( symmetric ) temp /= 2;
		row[count++] = MatrixEntry< Real >( nodeIndex-offset , temp );
		for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
			if( (x!=2 || y!=2 || z!=2) && neighbors5.neighbors[x][y][z] && neighbors5.neighbors[x][y][z]->nodeData.nodeIndex>=0 && ( !symmetric || neighbors5.neighbors[x][y][z]->nodeData.nodeIndex>=nodeIndex ) )
			{
				const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
				int _d , _off[3];
				_node->depthAndOffset( _d , _off );
				Real temp = Real( GetLaplacian( integrator , d , off , _off , false ) );
				if( _constrainValues ) temp += pointValues[x][y][z];
				if( symmetric && x==2 && y==2 && z==2 ) temp /= 2;
				row[count++] = MatrixEntry< Real >( _node->nodeData.nodeIndex-offset , temp );
			}
	}
	return count;
}
// if( scatter ) normals come from the center ndoe
// else          normals come from the neighbors
template< class Real , int Degree >
void Octree< Real , Degree >::SetDivergenceStencil( int depth , const typename BSplineData< Degree >::Integrator& integrator , Stencil< Point3D< double > , 5 >& stencil , bool scatter ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	int offset[] = { center , center , center };
	for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
	{
		int _offset[] = { x+center-2 , y+center-2 , z+center-2 };
		if( scatter ) stencil.values[x][y][z] = GetDivergence1( integrator , depth , offset , _offset , false );
		else          stencil.values[x][y][z] = GetDivergence2( integrator , depth , offset , _offset , false );
	}
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetDivergenceStencils( int depth , const typename BSplineData< Degree >::Integrator& integrator , Stencil< Point3D< double > ,  5 > stencils[2][2][2] , bool scatter ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ )
	{
		int offset[] = { center+i , center+j , center+k };
		for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
		{
			int _offset[] = { x-2+center/2 , y-2+center/2 , z-2+center/2 };
			if( scatter ) stencils[i][j][k].values[x][y][z] = GetDivergence1( integrator , depth , offset , _offset , true );
			else          stencils[i][j][k].values[x][y][z] = GetDivergence2( integrator , depth , offset , _offset , true );
		}
	}
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetLaplacianStencil( int depth , const typename BSplineData< Degree >::Integrator& integrator , Stencil< double , 5 >& stencil ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	int offset[] = { center , center , center };
	for( int x=-2 ; x<=2 ; x++ ) for( int y=-2 ; y<=2 ; y++ ) for( int z=-2 ; z<=2 ; z++ )
	{
		int _offset[] = { x+center , y+center , z+center };
		stencil.values[x+2][y+2][z+2] = GetLaplacian( integrator , depth , offset , _offset , false );
	}
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetLaplacianStencils( int depth , const typename BSplineData< Degree >::Integrator& integrator , Stencil< double , 5 > stencils[2][2][2] ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ )
	{
		int offset[] = { center+i , center+j , center+k };
		for( int x=-2 ; x<=2 ; x++ ) for( int y=-2 ; y<=2 ; y++ ) for( int z=-2 ; z<=2 ; z++ )
		{
			int _offset[] = { x+center/2 , y+center/2 , z+center/2 };
			stencils[i][j][k].values[x+2][y+2][z+2] = GetLaplacian( integrator , depth , offset , _offset , true );
		}
	}
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetCenterEvaluationStencil( const typename BSplineData< Degree >::template CenterEvaluator< 1 >& evaluator , int depth , Stencil< double , 3 >& stencil ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
	{
		int off[] = { center+x-1 , center+y-1 , center+z-1 };
		stencil.values[x][y][z] = Real( evaluator.value( depth , center , off[0] , false , false ) * evaluator.value( depth , center , off[1] , false , false ) * evaluator.value( depth , center , off[2] , false , false ) );
	}
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetCenterEvaluationStencils( const typename BSplineData< Degree >::template CenterEvaluator< 1 >& evaluator , int depth , Stencil< double , 3 > stencils[8] ) const
{
	if( depth<3 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int idx[] = { center+cx , center+cy , center+cz };
		for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
		{
			int off[] = { center/2+x-1 , center/2+y-1 , center/2+z-1 };
			stencils[Cube::CornerIndex( cx , cy , cz ) ].values[x][y][z] = Real( evaluator.value( depth , idx[0] , off[0] , false , true ) * evaluator.value( depth , idx[1] , off[1] , false , true ) * evaluator.value( depth , idx[2] , off[2] , false , true ) );
		}
	}
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetCornerEvaluationStencil( const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< double , 3 > stencil[8] ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );
		for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
		{
			int off[] = { center+x-1 , center+y-1 , center+z-1 };
			stencil[c].values[x][y][z] = evaluator.value( depth , center , cx , off[0] , false , false ) * evaluator.value( depth , center , cy , off[1] , false , false ) * evaluator.value( depth , center , cz , off[2] , false , false );
		}
	}
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetCornerEvaluationStencils( const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< double , 3 > stencils[8][8] ) const
{
	if( depth<3 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );
		for( int _cx=0 ; _cx<2 ; _cx++ ) for( int _cy=0 ; _cy<2 ; _cy++ ) for( int _cz=0 ; _cz<2 ; _cz++ )
		{
			int _c = Cube::CornerIndex( _cx , _cy , _cz );
			int idx[] = { center+_cx , center+_cy , center+_cz };
			for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
			{
				int off[] = { center/2+x-1 , center/2+y-1 , center/2+z-1 };
				stencils[c][_c].values[x][y][z] = evaluator.value( depth , idx[0] , cx , off[0] , false , true ) * evaluator.value( depth , idx[1] , cy , off[1] , false , true ) * evaluator.value( depth , idx[2] , cz , off[2] , false , true );
			}
		}
	}
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetCornerNormalEvaluationStencil( const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 5 > stencil[8] ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );
		for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
		{
			int off[] = { center+x-2 , center+y-2 , center+z-2 };
			double v [] = { evaluator.value( depth , center , cx , off[0] , false , false ) , evaluator.value( depth , center , cy , off[1] , false , false ) , evaluator.value( depth , center , cz , off[2] , false , false ) };
			double dv[] = { evaluator.value( depth , center , cx , off[0] , true  , false ) , evaluator.value( depth , center , cy , off[1] , true  , false ) , evaluator.value( depth , center , cz , off[2] , true  , false ) };
			stencil[c].values[x][y][z] = Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] );
		}
	}
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetCornerNormalEvaluationStencils( const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 5 > stencils[8][8] ) const
{
	if( depth<3 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );	// Which corner of the finer cube
		for( int _cx=0 ; _cx<2 ; _cx++ ) for( int _cy=0 ; _cy<2 ; _cy++ ) for( int _cz=0 ; _cz<2 ; _cz++ )
		{
			int _c = Cube::CornerIndex( _cx , _cy , _cz );	// Which child node
			int idx[] = { center+_cx , center+_cy , center+_cz };
			for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
			{
				int off[] = { center/2+x-2 , center/2+y-2 , center/2+z-2 };
				double v [] = { evaluator.value( depth , idx[0] , cx , off[0] , false , true ) , evaluator.value( depth , idx[1] , cy , off[1] , false , true ) , evaluator.value( depth , idx[2] , cz , off[2] , false , true ) };
				double dv[] = { evaluator.value( depth , idx[0] , cx , off[0] , true  , true ) , evaluator.value( depth , idx[1] , cy , off[1] , true  , true ) , evaluator.value( depth , idx[2] , cz , off[2] , true  , true ) };
				stencils[c][_c].values[x][y][z] = Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] );
			}
		}
	}
}

template< class Real , int Degree >
void Octree< Real , Degree >::UpdateCoarserSupportBounds( const TreeOctNode* node , int& startX , int& endX , int& startY , int& endY , int& startZ , int& endZ )
{
	if( node->parent )
	{
		int x , y , z , c = int( node - node->parent->children );
		Cube::FactorCornerIndex( c , x , y , z );
		if( x==0 ) endX = 4;
		else     startX = 1;
		if( y==0 ) endY = 4;
		else     startY = 1;
		if( z==0 ) endZ = 4;
		else     startZ = 1;
	}
}
// Given the solution @( depth ) add to the met constraints @( depth-1 )
template< class Real , int Degree >
void Octree< Real , Degree >::UpdateConstraintsFromFiner( const typename BSplineData< Degree >::Integrator& integrator , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real ) fineSolution , Pointer( Real ) coarseConstraints ) const
{
	if( depth<=_minDepth ) return;
	Stencil< double , 5 > stencils[2][2][2];
	// Get the stencil describing the Laplacian relating coefficients @(depth) with coefficients @(depth-1)
	SetLaplacianStencils( depth , integrator , stencils );
	size_t start = sNodes.nodeCount[depth] , end = sNodes.nodeCount[depth+1] , range = end-start;
	int lStart = sNodes.nodeCount[depth-1];
	typename TreeOctNode::NeighborKey3 neighborKey3;
	neighborKey3.set( depth-1 );
	memset( coarseConstraints , 0 , sizeof(Real)*(sNodes.nodeCount[depth]-sNodes.nodeCount[depth-1]) );

	// Iterate over the nodes @( depth )
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey3 )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		TreeOctNode* node = sNodes.treeNodes[i];

		bool insetSupported = _boundaryType!=0 || _IsInsetSupported( node );

		// Offset the coarser constraints using the solution from the current resolutions.
		int x , y , z , c;
		c = int( node - node->parent->children );
		Cube::FactorCornerIndex( c , x , y , z );
		if( insetSupported )
		{
			typename TreeOctNode::Neighbors5 pNeighbors5;
			neighborKey3.getNeighbors( node->parent , pNeighbors5 );
			const Stencil< double , 5 >& lapStencil = stencils[x][y][z];

			Pointer( Real ) __coarseConstraints = coarseConstraints-lStart;
			bool isInterior;
			int d , off[3];
			{
				node->depthAndOffset( d , off );
				int o = _boundaryType==0 ? (1<<(d-2) ) : 0;
				int mn = 4+o , mx = (1<<d)-4-o;
				isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
			}
			// Offset the constraints using the solution from finer resolutions.
			int startX = 0 , endX = 5 , startY = 0 , endY = 5 , startZ = 0 , endZ = 5;
			UpdateCoarserSupportBounds( node , startX , endX , startY  , endY , startZ , endZ );

			Real solution = fineSolution[ node->nodeData.nodeIndex-start ];
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
				if( pNeighbors5.neighbors[x][y][z] && pNeighbors5.neighbors[x][y][z]->nodeData.nodeIndex>=0 )
				{
					const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
					if( isInterior )
#pragma omp atomic
						__coarseConstraints[ _node->nodeData.nodeIndex ] += Real( lapStencil.values[x][y][z] * solution );
					else
					{
						int _d , _off[3];
						_node->depthAndOffset( _d , _off );
#pragma omp atomic
						__coarseConstraints[ _node->nodeData.nodeIndex ] += Real( GetLaplacian( integrator , d , off , _off , true ) * solution );
					}
				}
		}
	}
}

template< class Real , int Degree >
void Octree< Real , Degree >::UpdateConstraintsFromCoarser( const PointInfo& pointInfo , const typename TreeOctNode::Neighbors5& neighbors5 , const typename TreeOctNode::Neighbors5& pNeighbors5 , TreeOctNode* node , Pointer( Real ) constraints , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::Integrator& integrator , const Stencil< double , 5 >& lapStencil ) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	if( node->depth()<=_minDepth ) return;
	bool isInterior;
	int d , off[3];
	{
		node->depthAndOffset( d , off );
		int o = _boundaryType==0 ? (1<<(d-2) ) : 0;
		int mn = 4+o , mx = (1<<d)-4-o;
		isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
	}
	Real constraint = Real( 0 );
	// Offset the constraints using the solution from lower resolutions.
	int startX = 0 , endX = 5 , startY = 0 , endY = 5 , startZ = 0 , endZ = 5;
	UpdateCoarserSupportBounds( node , startX , endX , startY  , endY , startZ , endZ );

	for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
		if( pNeighbors5.neighbors[x][y][z] && pNeighbors5.neighbors[x][y][z]->nodeData.nodeIndex>=0 )
		{
			const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
			Real _solution = metSolution[ _node->nodeData.nodeIndex ];
			{
				if( isInterior ) constraints[ node->nodeData.nodeIndex ] -= Real( lapStencil.values[x][y][z] * _solution );
				else
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					constraints[ node->nodeData.nodeIndex ] -= Real( GetLaplacian( integrator , d , off , _off , true ) * _solution );
				}
			}
		}
	if( _constrainValues )
	{
		double constraint = 0;
		int idx[3] ;
		node->centerIndex( idx );
		// Evaluate the current node's basis function at adjacent points
		for( int x=1 ; x<4 ; x++ ) for( int y=1 ; y<4 ; y++ ) for( int z=1 ; z<4 ; z++ )
			if( neighbors5.neighbors[x][y][z] && pointInfo.pointIndex( neighbors5.neighbors[x][y][z] )!=-1 )
			{
				const _PointData& pData = points[ pointInfo.pointIndex( neighbors5.neighbors[x][y][z] ) ];
				Real weightedPointValue = pData.weightedCoarserValue;
				Point3D< Real > p = pData.position;
				constraint += 
					_fData.baseBSplines[idx[0]][x-1]( p[0] ) *
					_fData.baseBSplines[idx[1]][y-1]( p[1] ) *
					_fData.baseBSplines[idx[2]][z-1]( p[2] ) * 
					weightedPointValue;
			}
		constraints[ node->nodeData.nodeIndex ] -= Real( constraint );
	}
}
struct UpSampleData
{
	int start;
	double v[2];
	UpSampleData( void ) { start = 0 , v[0] = v[1] = 0.; }
	UpSampleData( int s , double v1 , double v2 ) { start = s , v[0] = v1 , v[1] = v2; }
};
template< class Real , int Degree >
template< class C >
void Octree< Real , Degree >::DownSample( int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( C ) fineConstraints , Pointer( C ) coarseConstraints ) const
{
	if( depth==0 ) return;
	double cornerValue;
	if     ( _boundaryType==-1 ) cornerValue = 0.50;
	else if( _boundaryType== 1 ) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set( depth );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		int d , off[3];
		UpSampleData usData[3];
		sNodes.treeNodes[i]->depthAndOffset( d , off );
		for( int dd=0 ; dd<3 ; dd++ )
		{
			if     ( off[dd]  ==0          ) usData[dd] = UpSampleData( 1 , cornerValue , 0.00 );
			else if( off[dd]+1==(1<<depth) ) usData[dd] = UpSampleData( 0 , 0.00 , cornerValue );
			else if( off[dd]%2             ) usData[dd] = UpSampleData( 1 , 0.75 , 0.25 );
			else                             usData[dd] = UpSampleData( 0 , 0.25 , 0.75 );
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors( sNodes.treeNodes[i]->parent );
		C c = fineConstraints[ i-sNodes.nodeCount[depth] ];
		for( int ii=0 ; ii<2 ; ii++ )
		{
			int _ii = ii + usData[0].start;
			C cx = C( c*usData[0].v[ii] );
			for( int jj=0 ; jj<2 ; jj++ )
			{
				int _jj = jj + usData[1].start;
				C cxy = C( cx*usData[1].v[jj] );
				for( int kk=0 ; kk<2 ; kk++ )
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* pNode = neighbors.neighbors[_ii][_jj][_kk];
					if( pNode && pNode->nodeData.nodeIndex!=-1 )
#pragma omp atomic
						coarseConstraints[ pNode->nodeData.nodeIndex-sNodes.nodeCount[depth-1] ] += C( cxy*usData[2].v[kk] );
				}
			}
		}
	}
}
template< class Real , int Degree >
template< class C >
void Octree< Real , Degree >::UpSample( int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( C ) coarseCoefficients , Pointer( C ) fineCoefficients ) const
{
	double cornerValue;
	if     ( _boundaryType==-1 ) cornerValue = 0.50;
	else if( _boundaryType== 1 ) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	if( depth<=_minDepth ) return;

	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set( depth-1 );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		bool isInterior = true;
		TreeOctNode* node = sNodes.treeNodes[i];
		int d , off[3];
		UpSampleData usData[3];
		node->depthAndOffset( d , off );
		for( int d=0 ; d<3 ; d++ )
		{
			if     ( off[d]  ==0          ) usData[d] = UpSampleData( 1 , cornerValue , 0.00 ) , isInterior = false;
			else if( off[d]+1==(1<<depth) ) usData[d] = UpSampleData( 0 , 0.00 , cornerValue ) , isInterior = false;
			else if( off[d]%2             ) usData[d] = UpSampleData( 1 , 0.75 , 0.25 );
			else                            usData[d] = UpSampleData( 0 , 0.25 , 0.75 );
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors( node->parent );
		for( int ii=0 ; ii<2 ; ii++ )
		{
			int _ii = ii + usData[0].start;
			double dx = usData[0].v[ii];
			for( int jj=0 ; jj<2 ; jj++ )
			{
				int _jj = jj + usData[1].start;
				double dxy = dx * usData[1].v[jj];
				for( int kk=0 ; kk<2 ; kk++ )
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* node = neighbors.neighbors[_ii][_jj][_kk];
					if( node && node->nodeData.nodeIndex!=-1 )
					{
						double dxyz = dxy * usData[2].v[kk];
						int _i = node->nodeData.nodeIndex;
						fineCoefficients[ i-sNodes.nodeCount[depth] ] += coarseCoefficients[ _i-sNodes.nodeCount[depth-1] ] * Real( dxyz );
					}
				}
			}
		}
	}
}
// At each point @( depth ), evaluate the met solution @( depth-1 )
template< class Real , int Degree >
void Octree< Real , Degree >::SetPointValuesFromCoarser( PointInfo& pointInfo , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real ) coarseCoefficients )
{
	std::vector< _PointData >& points = pointInfo.points;
	// For every node at the current depth
	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set( depth );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		int pIdx = pointInfo.pointIndex( sNodes.treeNodes[i] );
		if( pIdx!=-1 )
		{
			neighborKey.getNeighbors( sNodes.treeNodes[i] );
			points[ pIdx ].weightedCoarserValue = _WeightedCoarserFunctionValue( points[pIdx] , neighborKey , sNodes.treeNodes[i] , coarseCoefficients-_sNodes.nodeCount[depth-1] );
		}
	}
}
template< class Real , int Degree >
Real Octree< Real , Degree >::_WeightedCoarserFunctionValue( const _PointData& pointData , const typename TreeOctNode::NeighborKey3& neighborKey , const TreeOctNode* pointNode , ConstPointer( Real ) coarseCoefficients ) const
{
	double pointValue = 0;
	int depth = pointNode->depth();
	if( _boundaryType==-1 && depth==0 ) return Real(-0.5) * pointData.weight;

	if( depth<=_minDepth ) return Real(0.);

	Real weight       = pointData.weight;
	Point3D< Real > p = pointData.position;

	// Iterate over all basis functions that overlap the point at the coarser resolutions
	{
		int d , _idx[3];
		const typename TreeOctNode::Neighbors3& neighbors = neighborKey.neighbors[depth-1];
		neighbors.neighbors[1][1][1]->depthAndOffset( d , _idx );
		_idx[0] = BinaryNode< double >::CenterIndex( d , _idx[0]-1 );
		_idx[1] = BinaryNode< double >::CenterIndex( d , _idx[1]-1 );
		_idx[2] = BinaryNode< double >::CenterIndex( d , _idx[2]-1 );

		for( int j=0 ; j<3 ; j++ )
		{
#if ROBERTO_TOLDO_FIX
			double xValue = 0;
			if( _idx[0]+j>=0 && _idx[0]+j<((1<<depth)-1) ) xValue = _fData.baseBSplines[ _idx[0]+j ][2-j]( p[0] );
			else continue;
#else // !ROBERTO_TOLDO_FIX
			double xValue = _fData.baseBSplines[ _idx[0]+j ][2-j]( p[0] );
#endif // ROBERTO_TOLDO_FIX
			for( int k=0 ; k<3 ; k++ )
			{
#if ROBERTO_TOLDO_FIX
				double xyValue = 0;
				if( _idx[1]+k>=0 && _idx[1]+k<((1<<depth)-1) ) xyValue = xValue * _fData.baseBSplines[ _idx[1]+k ][2-k]( p[1] );
				else continue;
#else // !ROBERTO_TOLDO_FIX
				double xyValue = xValue * _fData.baseBSplines[ _idx[1]+k ][2-k]( p[1] );
#endif // ROBERTO_TOLDO_FIX
				double _pointValue = 0;
				for( int l=0 ; l<3 ; l++ )
				{
					const TreeOctNode* basisNode = neighbors.neighbors[j][k][l];
#if ROBERTO_TOLDO_FIX
					if( basisNode && basisNode->nodeData.nodeIndex>=0 && _idx[2]+l>=0 && _idx[2]+l<((1<<depth)-1) )
						_pointValue += _fData.baseBSplines[ _idx[2]+l ][2-l]( p[2] ) * double( coarseCoefficients[basisNode->nodeData.nodeIndex] );
#else // !ROBERTO_TOLDO_FIX
					if( basisNode && basisNode->nodeData.nodeIndex>=0 )
						_pointValue += _fData.baseBSplines[ _idx[2]+l ][2-l]( p[2] ) * double( coarseCoefficients[basisNode->nodeData.nodeIndex] );
#endif // ROBERTO_TOLDO_FIX
				}
				pointValue += _pointValue * xyValue;
			}
		}
	}
	if( _boundaryType==-1 ) pointValue -= 0.5;
	return Real( pointValue * weight );
}
template< class Real , int Degree >
void Octree< Real , Degree >::SetPointConstraintsFromFiner( const PointInfo& pointInfo , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real ) finerCoefficients , Pointer( Real ) coarserConstraints ) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	// Note: We can't iterate over the finer point nodes as the point weights might be
	// scaled incorrectly, due to the adaptive exponent. So instead, we will iterate
	// over the coarser nodes and evaluate the finer solution at the associated points.
	if( !depth ) return;
	size_t start = sNodes.nodeCount[depth-1] , end = sNodes.nodeCount[depth] , range = end-start;
	memset( coarserConstraints , 0 , sizeof( Real ) * ( sNodes.nodeCount[depth]-sNodes.nodeCount[depth-1] ) );
	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set( depth-1 );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
	for( int i=sNodes.nodeCount[depth-1] ; i<sNodes.nodeCount[depth] ; i++ )
	{
		int pIdx = pointInfo.pointIndex( sNodes.treeNodes[i] );
		if( pIdx!=-1 )
		{
			typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors( sNodes.treeNodes[i] );
			// Evaluate the solution @( depth ) at the current point @( depth-1 )
			{
				Real finerPointValue = _WeightedFinerFunctionValue( points[pIdx] , neighborKey , sNodes.treeNodes[i] , finerCoefficients-sNodes.nodeCount[depth] );
				Point3D< Real > p = points[ pIdx ].position;
				// Update constraints for all nodes @( depth-1 ) that overlap the point
				int d , idx[3];
				neighbors.neighbors[1][1][1]->depthAndOffset( d, idx );
				// Set the (offset) index to the top-left-front corner of the 3x3x3 block of b-splines
				// overlapping the point.
				idx[0] = BinaryNode< double >::CenterIndex( d , idx[0]-1 );
				idx[1] = BinaryNode< double >::CenterIndex( d , idx[1]-1 );
				idx[2] = BinaryNode< double >::CenterIndex( d , idx[2]-1 );
				for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
					if( neighbors.neighbors[x][y][z] && neighbors.neighbors[x][y][z]->nodeData.nodeIndex!=-1 )
					{
#pragma omp atomic
						coarserConstraints[ neighbors.neighbors[x][y][z]->nodeData.nodeIndex - sNodes.nodeCount[depth-1] ] +=
							Real(
							_fData.baseBSplines[idx[0]+x][2-x]( p[0] ) *
							_fData.baseBSplines[idx[1]+y][2-y]( p[1] ) *
							_fData.baseBSplines[idx[2]+z][2-z]( p[2] ) * 
							finerPointValue
							);
					}
			}
		}
	}
}
template< class Real , int Degree >
Real Octree< Real , Degree >::_WeightedFinerFunctionValue( const _PointData& pointData , const typename TreeOctNode::NeighborKey3& neighborKey , const TreeOctNode* pointNode , ConstPointer( Real ) finerCoefficients ) const
{
	typename TreeOctNode::Neighbors3 childNeighbors;
	double pointValue = 0;
	int depth = pointNode->depth();
	Real weight       = pointData.weight;
	Point3D< Real > p = pointData.position;
	neighborKey.getChildNeighbors( p , depth , childNeighbors );
	// Iterate over all finer basis functions that overlap the point at the coarser resolutions
	int d , idx[3];
	{
		Point3D< Real > c;
		Real w;
		neighborKey.neighbors[depth].neighbors[1][1][1]->depthAndOffset( d , idx );
		neighborKey.neighbors[depth].neighbors[1][1][1]->centerAndWidth( c , w );
		d++;
		idx[0] *= 2 , idx[1] *= 2 , idx[2] *= 2;
		int cIndex=TreeOctNode::CornerIndex( c , p );
		if( cIndex&1 ) idx[0]++;
		if( cIndex&2 ) idx[1]++;
		if( cIndex&4 ) idx[2]++;
	}
	// Center the indexing at the top-left-front corner
	idx[0] = BinaryNode< double >::CenterIndex( d , idx[0]-1 );
	idx[1] = BinaryNode< double >::CenterIndex( d , idx[1]-1 );
	idx[2] = BinaryNode< double >::CenterIndex( d , idx[2]-1 );

	for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) for( int l=0 ; l<3 ; l++ )
	{
		const TreeOctNode* basisNode = childNeighbors.neighbors[j][k][l];
		if( basisNode && basisNode->nodeData.nodeIndex>=0 )
			pointValue += 
			_fData.baseBSplines[ idx[0]+j ][2-j]( p[0] ) *
			_fData.baseBSplines[ idx[1]+k ][2-k]( p[1] ) *
			_fData.baseBSplines[ idx[2]+l ][2-l]( p[2] ) *
			double( finerCoefficients[ basisNode->nodeData.nodeIndex ] );
	}
	if( _boundaryType==-1 ) pointValue -= Real(0.5);
	return Real( pointValue * weight );
}
template< class Real , int Degree >
int Octree< Real , Degree >::GetSliceMatrixAndUpdateConstraints( const PointInfo& pointInfo , SparseMatrix< Real >& matrix , Pointer( Real ) constraints , const typename BSplineData< Degree >::Integrator& integrator , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real ) metSolution , bool coarseToFine , int nStart , int nEnd )
{
	size_t range = nEnd-nStart;
	Stencil< double , 5 > stencil , stencils[2][2][2];
	SetLaplacianStencil ( depth , integrator , stencil );
	SetLaplacianStencils( depth , integrator , stencils );
	matrix.Resize( (int)range );
	typename TreeOctNode::NeighborKey3 neighborKey3;
	neighborKey3.set( depth );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey3 )
	for( int i=0 ; i<range ; i++ )
	{
		TreeOctNode* node = sNodes.treeNodes[i+nStart];
		// Get the matrix row size
		bool insetSupported = _boundaryType!=0 || _IsInsetSupported( node );
		typename TreeOctNode::Neighbors5 neighbors5;
		if( insetSupported ) neighborKey3.getNeighbors( node , neighbors5 );
		int count = insetSupported ? GetMatrixRowSize( neighbors5 , false ) : 1;

		// Allocate memory for the row
#pragma omp critical (matrix_set_row_size)
		{
			matrix.SetRowSize( i , count );
		}

		// Set the row entries
		if( insetSupported ) matrix.rowSizes[i] = SetMatrixRow( pointInfo , neighbors5 , matrix[i] , sNodes.nodeCount[depth] , integrator , stencil , false );
		else
		{
			matrix[i][0] = MatrixEntry< Real >( i , Real(1) );
			matrix.rowSizes[i] = 1;
		}

		if( depth>_minDepth )
		{
			// Offset the constraints using the solution from lower resolutions.
			int x , y , z , c;
			if( node->parent )
			{
				c = int( node - node->parent->children );
				Cube::FactorCornerIndex( c , x , y , z );
			}
			else x = y = z = 0;
			if( insetSupported && coarseToFine )
			{
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey3.getNeighbors( node->parent , pNeighbors5 );
				UpdateConstraintsFromCoarser( pointInfo , neighbors5 , pNeighbors5 , node , constraints , metSolution , integrator , stencils[x][y][z] );
			}
		}
	}
	return 1;
}
template< class Real , int Degree >
int Octree< Real , Degree >::GetMatrixAndUpdateConstraints( const PointInfo& pointInfo , SparseSymmetricMatrix< Real >& matrix , Pointer( Real ) constraints , const typename BSplineData< Degree >::Integrator& integrator , int depth , const SortedTreeNodes< Real >& sNodes , ConstPointer( Real ) metSolution , bool coarseToFine )
{
	size_t start = sNodes.nodeCount[depth] , end = sNodes.nodeCount[depth+1] , range = end-start;
	Stencil< double , 5 > stencil , stencils[2][2][2];
	SetLaplacianStencil ( depth , integrator , stencil );
	SetLaplacianStencils( depth , integrator , stencils );
	matrix.Resize( (int)range );
	typename TreeOctNode::NeighborKey3 neighborKey3;
	neighborKey3.set( depth );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey3 )
	for( int i=0 ; i<range ; i++ )
	{
		TreeOctNode* node = sNodes.treeNodes[i+start];
		// Get the matrix row size
		bool insetSupported = _boundaryType!=0 || _IsInsetSupported( node );
		typename TreeOctNode::Neighbors5 neighbors5;
		if( insetSupported ) neighborKey3.getNeighbors( node , neighbors5 );
		int count = insetSupported ? GetMatrixRowSize( neighbors5 , true ) : 1;

		// Allocate memory for the row
#pragma omp critical (matrix_set_row_size)
		matrix.SetRowSize( i , count );

		// Set the row entries
		if( insetSupported ) matrix.rowSizes[i] = SetMatrixRow( pointInfo , neighbors5 , matrix[i] , (int)start , integrator , stencil , true );
		else
		{
			matrix[i][0] = MatrixEntry< Real >( i , Real(1) );
			matrix.rowSizes[i] = 1;
		}
		if( depth>_minDepth )
		{
			// Offset the constraints using the solution from lower resolutions.
			int x , y , z , c;
			if( node->parent )
			{
				c = int( node - node->parent->children );
				Cube::FactorCornerIndex( c , x , y , z );
			}
			else x = y = z = 0;
			if( insetSupported && coarseToFine )
			{
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey3.getNeighbors( node->parent , pNeighbors5 );
				UpdateConstraintsFromCoarser( pointInfo , neighbors5 , pNeighbors5 , node , constraints , metSolution , integrator , stencils[x][y][z] );
			}
		}
	}
	return 1;
}

template< class Real , int Degree >
Pointer( Real ) Octree< Real , Degree >::SolveSystem( PointInfo& pointInfo , Pointer( Real ) constraints , bool showResidual , int iters , int maxSolveDepth , int cgDepth , double accuracy )
{
	int iter=0;
	typename BSplineData< Degree >::Integrator integrator;
	_fData.setIntegrator( integrator , _boundaryType==0 );
	iters = std::max< int >( 0 , iters );
	if( _boundaryType==0 ) maxSolveDepth++ , cgDepth++;

	Pointer( Real ) solution = AllocPointer< Real >( _sNodes.nodeCount[_sNodes.maxDepth] );
	memset( solution , 0 , sizeof(Real)*_sNodes.nodeCount[_sNodes.maxDepth] );

	solution[0] = 0;

	// The prolongation phase
	{
		std::vector< Real > metSolution( _sNodes.nodeCount[ _sNodes.maxDepth-1 ] , 0 );
		for( int d=_minDepth ; d<_sNodes.maxDepth ; d++ )
		{
			DumpOutput( "Depth[%d/%d]: %d\n" , _boundaryType==0 ? d-1 : d , _boundaryType==0 ? _sNodes.maxDepth-2 : _sNodes.maxDepth-1 , _sNodes.nodeCount[d+1]-_sNodes.nodeCount[d] );
			if( d==_minDepth )
				_SolveSystemCG( pointInfo , d , integrator , _sNodes , solution , constraints , GetPointer( metSolution ) , _sNodes.nodeCount[_minDepth+1]-_sNodes.nodeCount[_minDepth] , true , showResidual, NULL , NULL , NULL );
			else
			{
				if( d>cgDepth ) iter += _SolveSystemGS( pointInfo , d , integrator , _sNodes , solution , constraints , GetPointer( metSolution ) , d>maxSolveDepth ? 0 : iters , true , showResidual , NULL , NULL , NULL );
				else            iter += _SolveSystemCG( pointInfo , d , integrator , _sNodes , solution , constraints , GetPointer( metSolution ) , d>maxSolveDepth ? 0 : iters , true , showResidual , NULL , NULL , NULL , accuracy );
			}
		}
	}

	return solution;
}
template< class Real , int Degree >
void Octree< Real , Degree >::_setMultiColorIndices( int start , int end , std::vector< std::vector< int > >& indices ) const
{
	const int modulus = 3;
	indices.resize( modulus*modulus*modulus );
	int count[modulus*modulus*modulus];
	memset( count , 0 , sizeof(int)*modulus*modulus*modulus );
#pragma omp parallel for num_threads( threads )
	for( int i=start ; i<end ; i++ )
	{
		int d , off[3];
		_sNodes.treeNodes[i]->depthAndOffset( d , off );
		int idx = (modulus*modulus) * ( off[2]%modulus ) + modulus * ( off[1]%modulus ) + ( off[0]%modulus );
#pragma omp atomic
		count[idx]++;
	}

	for( int i=0 ; i<modulus*modulus*modulus ; i++ ) indices[i].reserve( count[i] ) , count[i]=0;

	for( int i=start ; i<end ; i++ )
	{
		int d , off[3];
		_sNodes.treeNodes[i]->depthAndOffset( d , off );
		int idx = (modulus*modulus) * ( off[2]%modulus ) + modulus * ( off[1]%modulus ) + ( off[0]%modulus );
		indices[idx].push_back( _sNodes.treeNodes[i]->nodeData.nodeIndex - start );
	}
}
template< class Real , int Degree >
int Octree< Real , Degree >::_SolveSystemGS( PointInfo& pointInfo , int depth , const typename BSplineData< Degree >::Integrator& integrator , const SortedTreeNodes< Real >& sNodes , Pointer( Real ) solution , Pointer( Real ) constraints , Pointer( Real ) metSolutionConstraints , int iters , bool coarseToFine , bool showResidual , double* bNorm2 , double* inRNorm2 , double* outRNorm2 , bool forceSilent )
{
	Pointer( Real ) metSolution = NullPointer< Real >();
	Pointer( Real ) metConstraints = NullPointer< Real >();
	if( coarseToFine ) metSolution    = metSolutionConstraints;	// This stores the up-sampled solution up to depth-2
	else               metConstraints = metSolutionConstraints; // This stores the down-sampled constraints up to depth

	double _maxMemoryUsage = maxMemoryUsage;
	maxMemoryUsage = 0;
	Vector< Real > X , B;
	int slices = 1<<depth;
	double systemTime=0. , solveTime=0. , updateTime=0. ,  evaluateTime = 0.;
#if BROODED_SLICES
	slices >>= 1;
#endif // BROODED_SLICES
	std::vector< int > offsets( slices+1 , 0 );
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		int d , off[3];
		sNodes.treeNodes[i]->depthAndOffset( d , off );
#if BROODED_SLICES
		offsets[ off[2]>>1 ]++;
#else // !BROODED_SLICES
		offsets[ off[2] ]++;
#endif // BROODED_SLICES
	}
	for( int i=1 ; i<slices ; i++ )  offsets[i] += offsets[i-1];
	for( int i=slices ; i>=1 ; i-- ) offsets[i]  = offsets[i-1];
	offsets[0] = 0;

	X.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );
	B.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );
	if( coarseToFine )
	{
		if( depth>_minDepth )
		{
			// Up-sample the cumulative change in solution @(depth-2) into the cumulative change in solution @(depth-1)
			if( depth-2>=_minDepth ) UpSample( depth-1 , sNodes , ( ConstPointer( Real ) )metSolution+_sNodes.nodeCount[depth-2] , metSolution+_sNodes.nodeCount[depth-1] );
			// Add in the change in solution @(depth-1)
#pragma omp parallel for num_threads( threads )
			for( int i=_sNodes.nodeCount[depth-1] ; i<_sNodes.nodeCount[depth] ; i++ ) metSolution[i] += solution[i];
			// Evaluate the points @(depth) using the cumulative change in solution @(depth-1)
			if( _constrainValues )
			{
				evaluateTime = Time();
				SetPointValuesFromCoarser( pointInfo , depth , sNodes , metSolution+_sNodes.nodeCount[depth-1] );
				evaluateTime = Time() - evaluateTime;
			}
		}
	}
	else if( depth<_sNodes.maxDepth-1 )
		for( int i=_sNodes.nodeCount[depth] ; i<_sNodes.nodeCount[depth+1] ; i++ ) constraints[i] -= metConstraints[i];
	// Initialize with the previously computed solution
#pragma omp parallel for num_threads( threads )
	for( int i=_sNodes.nodeCount[depth] ; i<_sNodes.nodeCount[depth+1] ; i++ ) X[ i-_sNodes.nodeCount[depth] ] = solution[i];
	double bNorm=0 , inRNorm=0 , outRNorm=0;
	if( depth>=_minDepth )
	{
#if BROODED_SLICES
		int frontOffset = ( showResidual || inRNorm2 ) ? 1 : 0;
		int backOffset = ( showResidual || outRNorm2 ) ? 1 : 0;
		int solveSlices = std::min< int >( iters , slices ) , matrixSlices = std::max< int >( 1 , std::min< int >( solveSlices+frontOffset+backOffset , slices ) );
#else // !BROODED_SLICES
		int frontOffset = ( showResidual || inRNorm2 ) ? 2 : 0;
		int backOffset = ( showResidual || outRNorm2 ) ? 2 : 0;
		int solveSlices = std::min< int >( 2*iters-1 , slices ) , matrixSlices = std::max< int >( 1 , std::min< int >( solveSlices+frontOffset+backOffset , slices ) );
#endif // BROODED_SLICES
		std::vector< SparseMatrix< Real > > _M( matrixSlices );
		std::vector< std::vector< std::vector< int > > > __mcIndices( std::max< int >( 0 , solveSlices ) );

		int dir = coarseToFine ? -1 : 1 , start = coarseToFine ? slices-1 : 0 , end = coarseToFine ? -1 : slices;
#if BROODED_SLICES
		for( int frontSlice=start-frontOffset*dir , backSlice = frontSlice-(iters-1)*dir ; backSlice!=end+backOffset*dir ; frontSlice+=dir , backSlice+=dir )
#else // !BROODED_SLICES
		for( int frontSlice=start-frontOffset*dir , backSlice = frontSlice-2*(iters-1)*dir ; backSlice!=end+backOffset*dir ; frontSlice+=dir , backSlice+=dir )
#endif // BROODED_SLICES
		{
			double t;
			if( frontSlice+frontOffset*dir>=0 && frontSlice+frontOffset*dir<slices )
			{
				int s = frontSlice+frontOffset*dir , _s = s % matrixSlices;
				t = Time();
				GetSliceMatrixAndUpdateConstraints( pointInfo , _M[_s] , constraints , integrator , depth , sNodes , metSolution , coarseToFine , offsets[s]+sNodes.nodeCount[depth] , offsets[s+1]+sNodes.nodeCount[depth] );
				systemTime += Time()-t;
				Pointer( TreeOctNode* ) const nodes = sNodes.treeNodes + sNodes.nodeCount[depth];
				for( int i=offsets[s] ; i<offsets[s+1] ; i++ )
				{
					if( _boundaryType!=0 || _IsInsetSupported( nodes[i] ) ) B[i] = constraints[ nodes[i]->nodeData.nodeIndex ];
					else                                                    B[i] = Real(0);
				}
				if( showResidual || inRNorm2 )
#pragma omp parallel for num_threads( threads ) reduction( + : bNorm , inRNorm )
					for( int j=0 ; j<_M[_s].rows ; j++ )
					{
						Real temp = Real(0);
						ConstPointer( MatrixEntry< Real > ) start = _M[_s][j];
						ConstPointer( MatrixEntry< Real > ) end = start + _M[_s].rowSizes[j];
						ConstPointer( MatrixEntry< Real > ) e;
						for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
						Real b = B[ j + offsets[s] ] ;
						bNorm += b*b;
						inRNorm += (temp-b) * (temp-b);
					}
				else if( bNorm2 )
#pragma omp parallel for num_threads( threads ) reduction( + : bNorm )
					for( int j=0 ; j<_M[_s].rows ; j++ )
					{
						Real b = B[ j + offsets[s] ] ;
						bNorm += b*b;
					}
			}
			t = Time();
			if( iters && frontSlice>=0 && frontSlice<slices )
			{
				int s = frontSlice , _s = s % matrixSlices , __s = s % solveSlices;
				for( int i=0 ; i<int( __mcIndices[__s].size() ) ; i++ ) __mcIndices[__s][i].clear();
				_setMultiColorIndices( sNodes.nodeCount[depth]+offsets[s] , sNodes.nodeCount[depth]+offsets[s+1] , __mcIndices[__s] );
			}
#if BROODED_SLICES
			for( int slice=frontSlice ; slice*dir>=backSlice*dir ; slice-=dir )
#else // !BROODED_SLICES
			for( int slice=frontSlice ; slice*dir>=backSlice*dir ; slice-=2*dir )
#endif // BROODED_SLICES
				if( slice>=0 && slice<slices )
				{
					int s = slice , _s = s % matrixSlices , __s = s % solveSlices;
					SparseMatrix< Real >::SolveGS( __mcIndices[__s] , _M[_s] , B , X , !coarseToFine , threads , offsets[s] );
				}
			solveTime += Time() - t;
			if( (showResidual || outRNorm2) && backSlice-backOffset*dir>=0 && backSlice-backOffset*dir<slices )
			{
				int s = backSlice-backOffset*dir , _s = s % matrixSlices;
#pragma omp parallel for num_threads( threads ) reduction( + : outRNorm )
				for( int j=0 ; j<_M[_s].rows ; j++ )
				{
					Real temp = Real(0);
					ConstPointer( MatrixEntry< Real > ) start = _M[_s][j];
					ConstPointer( MatrixEntry< Real > ) end = start + _M[_s].rowSizes[j];
					ConstPointer( MatrixEntry< Real > ) e;
					for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
					Real b = B[ j + offsets[s] ];
					outRNorm += (temp-b) * (temp-b);
				}
			}
		}
	}

	if( bNorm2 ) bNorm2[depth] = bNorm;
	if( inRNorm2 ) inRNorm2[depth] = inRNorm;
	if( outRNorm2 ) outRNorm2[depth] = outRNorm;
	if( showResidual && iters )
	{
		for( int i=0 ; i<depth ; i++ ) printf( "  " );
		printf( "GS: %.4e -> %.4e -> %.4e (%.2e) [%d]\n" , sqrt( bNorm ) , sqrt( inRNorm ) , sqrt( outRNorm ) , sqrt( outRNorm/bNorm ) , iters );
	}

	// Copy the old solution into the buffer, write in the new solution, compute the change, and update the met constraints
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ ) solution[i] = X[ i-sNodes.nodeCount[depth] ];
	if( !coarseToFine && depth>_minDepth )
	{
		// Explicitly compute the restriction of the met solution onto the coarser nodes
		// and down-sample the previous accumulation
		{
			UpdateConstraintsFromFiner( integrator , depth , sNodes , GetPointer( X ) , metConstraints+sNodes.nodeCount[depth-1] );
			if( _constrainValues ) SetPointConstraintsFromFiner( pointInfo , depth , sNodes , GetPointer( X ) , metConstraints+sNodes.nodeCount[depth-1] );
			if( depth<sNodes.maxDepth-1 ) DownSample( depth , sNodes , ( ConstPointer( Real ) )metConstraints+sNodes.nodeCount[depth] , metConstraints+sNodes.nodeCount[depth-1] );
		}
	}

	MemoryUsage();
	if( !forceSilent ) DumpOutput( "\tEvaluated / Got / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\n" , evaluateTime , systemTime , solveTime , float( maxMemoryUsage ) );
	maxMemoryUsage = std::max< double >( maxMemoryUsage , _maxMemoryUsage );

	return iters;
}
template< class Real , int Degree >
int Octree< Real , Degree >::_SolveSystemCG( PointInfo& pointInfo , int depth , const typename BSplineData< Degree >::Integrator& integrator , const SortedTreeNodes< Real >& sNodes , Pointer( Real ) solution , Pointer( Real ) constraints , Pointer( Real ) metSolutionConstraints , int iters , bool coarseToFine , bool showResidual , double* bNorm2 , double* inRNorm2 , double* outRNorm2 , double accuracy )
{
	Pointer( Real ) metSolution = NullPointer< Real >();
	Pointer( Real ) metConstraints = NullPointer< Real >();
	if( coarseToFine ) metSolution    = metSolutionConstraints;	// This stores the up-sampled solution up to depth-2
	else               metConstraints = metSolutionConstraints; // This stores the down-sampled constraints up to depth
	double _maxMemoryUsage = maxMemoryUsage;
	maxMemoryUsage = 0;
	int iter = 0;
	Vector< Real > X , B;
	SparseSymmetricMatrix< Real > M;
	double systemTime=0. , solveTime=0. , updateTime=0. ,  evaluateTime = 0.;
	X.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );
	if( coarseToFine )
	{
		if( depth>_minDepth )
		{
			// Up-sample the cumulative change in solution @(depth-2) into the cumulative change in solution @(depth-1)
			if( depth-2>=_minDepth ) UpSample( depth-1 , sNodes , ( ConstPointer( Real ) )metSolution+_sNodes.nodeCount[depth-2] , metSolution+_sNodes.nodeCount[depth-1] );
			// Add in the change in solution @(depth-1)
#pragma omp parallel for num_threads( threads )
			for( int i=_sNodes.nodeCount[depth-1] ; i<_sNodes.nodeCount[depth] ; i++ ) metSolution[i] += solution[i];
			// Evaluate the points @(depth) using the cumulative change in solution @(depth-1)
			if( _constrainValues )
			{
				evaluateTime = Time();
				SetPointValuesFromCoarser( pointInfo , depth , sNodes , metSolution+_sNodes.nodeCount[depth-1] );
				evaluateTime = Time() - evaluateTime;
			}
		}
	}
	else if( depth<_sNodes.maxDepth-1 )
		for( int i=_sNodes.nodeCount[depth] ; i<_sNodes.nodeCount[depth+1] ; i++ ) constraints[i] -= metConstraints[i];
	// Initialize with the previously computed solution
#pragma omp parallel for num_threads( threads )
	for( int i=_sNodes.nodeCount[depth] ; i<_sNodes.nodeCount[depth+1] ; i++ ) X[ i-_sNodes.nodeCount[depth] ] = solution[i];
	systemTime = Time();
	{
		// Get the system matrix (and adjust the right-hand-side based on the coarser solution if prolonging)
		if( coarseToFine ) GetMatrixAndUpdateConstraints( pointInfo , M , constraints , integrator , depth , sNodes , metSolution           , true  );
		else               GetMatrixAndUpdateConstraints( pointInfo , M , constraints , integrator , depth , sNodes , NullPointer< Real >() , false );
		// Set the constraint vector
		B.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );
		for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
			if( _boundaryType!=0 || _IsInsetSupported( sNodes.treeNodes[i] ) ) B[i-sNodes.nodeCount[depth]] = constraints[i];
			else                                                               B[i-sNodes.nodeCount[depth]] = Real(0);
	}
	systemTime = Time()-systemTime;

	solveTime = Time();
	// Solve the linear system
	accuracy = Real( accuracy / 100000 ) * M.rows;
	int res = 1<<depth;

	MapReduceVector< Real > mrVector;
	mrVector.resize( threads , M.rows );
	bool addDCTerm = (M.rows==res*res*res && !_constrainValues && _boundaryType!=-1);
	double bNorm , inRNorm , outRNorm;
	if( showResidual || bNorm2 ) bNorm = B.Norm( 2 );
	if( showResidual || inRNorm2 ) inRNorm = ( addDCTerm ? ( B - M * X - X.Average() ) : ( B - M * X ) ).Norm( 2 );

	if( _boundaryType==0 && depth>3 ) res -= 1<<(depth-2);
	if( iters ) iter += SparseSymmetricMatrix< Real >::SolveCG( M , B , iters , X , mrVector , Real( accuracy ) , 0 , addDCTerm );
	solveTime = Time()-solveTime;
	if( showResidual || outRNorm2 ) outRNorm = ( addDCTerm ? ( B - M * X - X.Average() ) : ( B - M * X ) ).Norm( 2 );
	if( bNorm2 ) bNorm2[depth] = bNorm * bNorm;
	if( inRNorm2 ) inRNorm2[depth] = inRNorm * inRNorm;
	if( outRNorm2 ) outRNorm2[depth] = outRNorm * outRNorm;
	if( showResidual && iters )
	{
		for( int i=0 ; i<depth ; i++ ) printf( "  " );
		printf( "CG: %.4e -> %.4e -> %.4e (%.2e) [%d]\n" , bNorm , inRNorm , outRNorm , outRNorm/bNorm , iter );
	}

	// Copy the old solution into the buffer, write in the new solution, compute the change, and update the met solution
	{
#pragma omp parallel for num_threads( threads )
		for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ ) solution[i] = X[ i-sNodes.nodeCount[depth] ];
		if( !coarseToFine && depth>_minDepth )
		{
			// Explicitly compute the restriction of the met solution onto the coarser nodes
			// and down-sample the previous accumulation
			{
				UpdateConstraintsFromFiner( integrator , depth , sNodes , GetPointer( X ) , metConstraints + sNodes.nodeCount[depth-1] );
				if( _constrainValues ) SetPointConstraintsFromFiner( pointInfo , depth , sNodes , GetPointer( X ) , metConstraints+sNodes.nodeCount[depth-1] );
				if( depth<sNodes.maxDepth-1 ) DownSample( depth , sNodes , ( ConstPointer( Real ) )metConstraints+sNodes.nodeCount[depth] , metConstraints+sNodes.nodeCount[depth-1] );
			}
		}
	}

	MemoryUsage();
	DumpOutput( "\tEvaluated / Got / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\n" , evaluateTime , systemTime , solveTime , float( maxMemoryUsage ) );
	maxMemoryUsage = std::max< double >( maxMemoryUsage , _maxMemoryUsage );
	return iter;
}
template< class Real , int Degree >
int Octree< Real , Degree >::HasNormals( TreeOctNode* node , const NormalInfo& normalInfo )
{
	int idx = normalInfo.normalIndex( node );
	if( idx>=0 )
	{
		const Point3D< Real >& normal = normalInfo.normals[ idx ];
		if( normal[0]!=0 || normal[1]!=0 || normal[2]!=0 ) return 1;
	}
	if( node->children ) for( int i=0 ; i<Cube::CORNERS ; i++ ) if( HasNormals( &node->children[i] , normalInfo ) ) return 1;
	return 0;
}
template< class Real , int Degree >
Pointer( Real ) Octree< Real , Degree >::SetLaplacianConstraints( const NormalInfo& normalInfo )
{
	// To set the Laplacian constraints, we iterate over the
	// splatted normals and compute the dot-product of the
	// divergence of the normal field with all the basis functions.
	// Within the same depth: set directly as a gather
	// Coarser depths 
	typename BSplineData< Degree >::Integrator integrator;
	_fData.setIntegrator( integrator , _boundaryType==0 );
	int maxDepth = _sNodes.maxDepth-1;
	Point3D< Real > zeroPoint;
	zeroPoint[0] = zeroPoint[1] = zeroPoint[2] = 0;
	Pointer( Real ) constraints = AllocPointer< Real >( _sNodes.nodeCount[_sNodes.maxDepth] );
	if( !constraints ) fprintf( stderr , "[ERROR] Failed to allocate constraints: %d * %d\n" , _sNodes.nodeCount[maxDepth+1] , sizeof( Real ) ) , exit( 0 );
	memset( constraints , 0 , sizeof(Real)*_sNodes.nodeCount[_sNodes.maxDepth] );
	Pointer( Real ) _constraints = AllocPointer< Real >( _sNodes.nodeCount[maxDepth] );
	if( !_constraints ) fprintf( stderr , "[ERROR] Failed to allocate _constraints: %d * %d\n" , _sNodes.nodeCount[maxDepth] , sizeof( Real ) ) , exit( 0 );
	memset( _constraints , 0 , sizeof(Real)*_sNodes.nodeCount[maxDepth] );

	for( int d=maxDepth ; d>=(_boundaryType==0?2:0) ; d-- )
	{
		int offset = d>0 ? _sNodes.treeNodes[ _sNodes.nodeCount[d-1] ]->nodeData.nodeIndex : 0;
		Stencil< Point3D< double > , 5 > stencil , stencils[2][2][2];
		SetDivergenceStencil ( d , integrator , stencil , false );
		SetDivergenceStencils( d , integrator , stencils , true );

		typename TreeOctNode::NeighborKey3 neighborKey3;
		neighborKey3.set( _fData.depth );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey3 )
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			TreeOctNode* node = _sNodes.treeNodes[i];
			int startX=0 , endX=5 , startY=0 , endY=5 , startZ=0 , endZ=5;
			int depth = node->depth();
			typename TreeOctNode::Neighbors5 neighbors5;
			neighborKey3.getNeighbors( node , neighbors5 );

			bool isInterior , isInterior2;
			{
				int d , off[3];
				node->depthAndOffset( d , off );
				int o = _boundaryType==0 ? (1<<(d-2)) : 0;
				int mn = 2+o , mx = (1<<d)-2-o;
				isInterior  = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
				mn += 2 , mx -= 2;
				isInterior2 = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
			}
			int cx , cy , cz;
			if( d )
			{
				int c = int( node - node->parent->children );
				Cube::FactorCornerIndex( c , cx , cy , cz );
			}
			else cx = cy = cz = 0;
			Stencil< Point3D< double > , 5 >& _stencil = stencils[cx][cy][cz];
			int d , off[3];
			node->depthAndOffset( d , off );
			// Set constraints from current depth
			// Gather the constraints from the vector-field at _node into the constraint stored with node
			{

				if( isInterior )
					for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
					{
						const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						if( _node )
						{
							int _idx = normalInfo.normalIndex( _node );
							if( _idx>=0 ) constraints[ node->nodeData.nodeIndex ] += Point3D< Real >::Dot( stencil.values[x][y][z] , normalInfo.normals[ _idx ] );
						}
					}
				else
					for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
					{
						const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						if( _node )
						{
							int _idx = normalInfo.normalIndex( _node );
							if( _idx>=0 )
							{
								int _d , _off[3];
								_node->depthAndOffset( _d , _off );
								constraints[ node->nodeData.nodeIndex ] += Real( GetDivergence2( integrator , d , off , _off , false , normalInfo.normals[ _idx ] ) );
							}
						}
					}
					UpdateCoarserSupportBounds( neighbors5.neighbors[2][2][2] , startX , endX , startY  , endY , startZ , endZ );
			}
			int idx = normalInfo.normalIndex( node );
			if( idx<0 ) continue;
			const Point3D< Real >& normal = normalInfo.normals[ idx ];
			if( normal[0]==0 && normal[1]==0 && normal[2]==0 ) continue;

			// Set the constraints for the parents
			if( depth>_minDepth )
			{
				neighborKey3.getNeighbors( node->parent , neighbors5 );

				for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
					if( neighbors5.neighbors[x][y][z] && neighbors5.neighbors[x][y][z]->nodeData.nodeIndex!=-1 )
					{
						TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						Real c;
						if( isInterior2 )
						{
							Point3D< double >& div = _stencil.values[x][y][z];
							c = Real( div[0] * normal[0] + div[1] * normal[1] + div[2] * normal[2] );
						}
						else
						{
							int _d , _off[3];
							_node->depthAndOffset( _d , _off );
							c = Real( GetDivergence1( integrator , d , off , _off , true , normal ) );
						}
#pragma omp atomic
						_constraints[ _node->nodeData.nodeIndex ] += c;
					}
			}
		}
		MemoryUsage();
	}
	std::vector< Point3D< Real > > coefficients( _sNodes.nodeCount[maxDepth] , zeroPoint );
	for( int d=maxDepth-1 ; d>=0 ; d-- )
	{
#pragma omp parallel for num_threads( threads )
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			TreeOctNode* node = _sNodes.treeNodes[i];
			int idx = normalInfo.normalIndex( node );
			if( idx<0 ) continue;
			const Point3D< Real >& normal = normalInfo.normals[ idx ];
			coefficients[i] += normal;
		}
	}
	// Fine-to-coarse down-sampling of constraints
	for( int d=maxDepth-1 ; d>=(_boundaryType==0?2:0) ; d-- ) DownSample( d , _sNodes , ( ConstPointer( Real ) )_constraints + _sNodes.nodeCount[d] , _constraints+_sNodes.nodeCount[d-1] );

	// Coarse-to-fine up-sampling of coefficients
	for( int d=(_boundaryType==0?2:0) ; d<maxDepth ; d++ ) UpSample( d , _sNodes , ( ConstPointer( Point3D< Real > ) ) GetPointer( coefficients ) + _sNodes.nodeCount[d-1] , GetPointer( coefficients ) + _sNodes.nodeCount[d] );

	// Add the accumulated constraints from all finer depths
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<_sNodes.nodeCount[maxDepth] ; i++ ) constraints[i] += _constraints[i];
	FreePointer( _constraints );

	// Compute the contribution from all coarser depths
	for( int d=0 ; d<=maxDepth ; d++ )
	{
		size_t start = _sNodes.nodeCount[d] , end = _sNodes.nodeCount[d+1] , range = end - start;
		Stencil< Point3D< double > , 5 > stencils[2][2][2];
		SetDivergenceStencils( d , integrator , stencils , false );
		typename TreeOctNode::NeighborKey3 neighborKey3;
		neighborKey3.set( maxDepth );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey3 )
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			TreeOctNode* node = _sNodes.treeNodes[i];
			int depth = node->depth();
			if( !depth ) continue;
			int startX=0 , endX=5 , startY=0 , endY=5 , startZ=0 , endZ=5;
			UpdateCoarserSupportBounds( node , startX , endX , startY  , endY , startZ , endZ );
			typename TreeOctNode::Neighbors5 neighbors5;
			neighborKey3.getNeighbors( node->parent , neighbors5 );

			bool isInterior;
			{
				int d , off[3];
				node->depthAndOffset( d , off );
				int o = _boundaryType==0 ? (1<<(d-2)) : 0;
				int mn = 4+o , mx = (1<<d)-4-o;
				isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
			}
			int cx , cy , cz;
			if( d )
			{
				int c = int( node - node->parent->children );
				Cube::FactorCornerIndex( c , cx , cy , cz );
			}
			else cx = cy = cz = 0;
			Stencil< Point3D< double > , 5 >& _stencil = stencils[cx][cy][cz];

			Real constraint = Real(0);
			int d , off[3];
			node->depthAndOffset( d , off );
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
				if( neighbors5.neighbors[x][y][z] && neighbors5.neighbors[x][y][z]->nodeData.nodeIndex!=-1 )
				{
					TreeOctNode* _node = neighbors5.neighbors[x][y][z];
					int _i = _node->nodeData.nodeIndex;
					if( isInterior )
					{
						Point3D< double >& div = _stencil.values[x][y][z];
						Point3D< Real >& normal = coefficients[_i];
						constraint += Real( div[0] * normal[0] + div[1] * normal[1] + div[2] * normal[2] );
					}
					else
					{
						int _d , _off[3];
						_node->depthAndOffset( _d , _off );
						constraint += Real( GetDivergence2( integrator , d , off , _off , true , coefficients[_i] ) );
					}
				}
				constraints[ node->nodeData.nodeIndex ] += constraint;
		}
	}

	MemoryUsage();
	return constraints;
}
template< class Real , int Degree >
void Octree< Real , Degree >::FaceEdgesFunction::Function( const TreeOctNode* node1 , const TreeOctNode* node2 )
{
	int mcIndex1 = mcIndices[ node1->nodeData.nodeIndex ];
	if( !node1->children && MarchingCubes::HasRoots( mcIndex1 ) )
	{
		RootInfo< Real > ri1 , ri2;
		typename hash_map< long long , std::pair< RootInfo< Real > , int > >::iterator iter;
		int isoTri[DIMENSION*MarchingCubes::MAX_TRIANGLES];
		int count=MarchingCubes::AddTriangleIndices( mcIndex1 , isoTri );

		for( int j=0 ; j<count ; j++ )
			for( int k=0 ; k<3 ; k++ )
				if( fIndex==Cube::FaceAdjacentToEdges( isoTri[j*3+k] , isoTri[j*3+((k+1)%3)] ) )
					if( GetRootIndex( mcIndices , node1 , isoTri[j*3+k] , maxDepth , *neighborKey3 , ri1 ) && GetRootIndex( mcIndices , node1 , isoTri[j*3+((k+1)%3)] , maxDepth , *neighborKey3 , ri2 ) )
					{
						long long key1=ri1.key , key2=ri2.key;
						edges->push_back( std::pair< RootInfo< Real > , RootInfo< Real > >( ri2 , ri1 ) );
						iter = vertexCount->find( key1 );
						if( iter==vertexCount->end() )
						{
							(*vertexCount)[key1].first = ri1;
							(*vertexCount)[key1].second=0;
						}
						iter=vertexCount->find(key2);
						if( iter==vertexCount->end() )
						{
							(*vertexCount)[key2].first = ri2;
							(*vertexCount)[key2].second=0;
						}
						(*vertexCount)[key1].second--;
						(*vertexCount)[key2].second++;
					}
					else fprintf( stderr , "Bad Edge 1: %d %d\n" , ri1.key , ri2.key );
	}
}
template< class Real , int Degree >
int Octree< Real , Degree >::refineBoundary( int subdivideDepth , std::vector< int >* map )
{
	// This implementation is somewhat tricky.
	// We would like to ensure that leaf-nodes across a subdivision boundary have the same depth.
	// We do this by calling the setNeighbors function.
	// The key is to implement this in a single pass through the leaves, ensuring that refinements don't propogate.
	// To this end, we do the minimal refinement that ensures that a cross boundary neighbor, and any of its cross-boundary
	// neighbors are all refined simultaneously.
	// For this reason, the implementation can only support nodes deeper than sDepth.
	bool flags[3][3][3];
	int maxDepth = tree.maxDepth();

	subdivideDepth = std::max< int >( subdivideDepth , 0 );
	if( _boundaryType==0 ) subdivideDepth += 2;
	subdivideDepth = std::min< int >( subdivideDepth , maxDepth );
	int sDepth = maxDepth - subdivideDepth;
	if( _boundaryType==0 ) sDepth = std::max< int >( 2 , sDepth );
	if( sDepth==0 )
	{
		_sNodes.set( tree , map );
		return sDepth;
	}

	// Ensure that face adjacent neighbors across the subdivision boundary exist to allow for
	// a consistent definition of the iso-surface
	typename TreeOctNode::NeighborKey3 nKey;
	nKey.set( maxDepth );
	for( TreeOctNode* leaf=tree.nextLeaf() ; leaf ; leaf=tree.nextLeaf( leaf ) )
		if( leaf->depth()>sDepth )
		{
			int d , off[3] , _off[3];
			leaf->depthAndOffset( d , off );
			int res = (1<<d)-1 , _res = ( 1<<(d-sDepth) )-1;
			_off[0] = off[0]&_res , _off[1] = off[1]&_res , _off[2] = off[2]&_res;
			bool boundary[2] = { (off[2]!=0 && _off[2]==0) , (off[2]!=res && _off[2]==_res) };

			if( boundary[0] || boundary[1] )
			{
				typename TreeOctNode::Neighbors3& neighbors = nKey.getNeighbors( leaf );
				for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) flags[i][j][k] = false;
				int z=0;
				if     ( boundary[0] && !neighbors.neighbors[1][1][0] ) z = -1;
				else if( boundary[1] && !neighbors.neighbors[1][1][2] ) z =  1;

				if( z )
				{
					flags[1][1][1+z] = true;
					nKey.setNeighbors( leaf , flags );
				}
			}
		}
	_sNodes.set( tree , map );
	MemoryUsage();
	return sDepth;
}
template< class Real , int Degree >
template< class Vertex >
void Octree< Real , Degree >::GetMCIsoTriangles( const std::vector< Real >* kernelDensityWeights , ConstPointer( Real ) solution , Real isoValue , int subdivideDepth , CoredMeshData< Vertex >* mesh , int fullDepthIso , int nonLinearFit , bool addBarycenter , bool polygonMesh )
{
	typename BSplineData< Degree >::template CornerEvaluator< 2 > evaluator;
	_fData.setCornerEvaluator( evaluator , 0 , _postDerivativeSmooth , _boundaryType==0 );

	int maxDepth = tree.maxDepth();
	int sDepth;
	{
		sDepth = std::max< int >( subdivideDepth , 0 );
		if( _boundaryType==0 ) sDepth += 2;
		sDepth = std::min< int >( sDepth , maxDepth );
		sDepth = maxDepth - sDepth;
		if( _boundaryType==0 ) sDepth = std::max< int >( 2 , sDepth );
	}

	RootData rootData , coarseRootData;
	std::vector< Vertex >* interiorVertices;
	std::vector< Real > metSolution( _sNodes.nodeCount[maxDepth] , 0 );
#pragma omp parallel for num_threads( threads )
	for( int i=_sNodes.nodeCount[_minDepth] ; i<_sNodes.nodeCount[maxDepth] ; i++ ) metSolution[i] = solution[i];
	for( int d=_minDepth ; d<maxDepth ; d++ ) UpSample( d , _sNodes , ( ConstPointer( Real ) )GetPointer( metSolution ) + _sNodes.nodeCount[d-1] , GetPointer( metSolution ) + _sNodes.nodeCount[d] );

	Pointer( unsigned char ) mcIndices = AllocPointer< unsigned char >( _sNodes.nodeCount[_sNodes.maxDepth] );
	if( !mcIndices ) fprintf( stderr , "[ERROR] Failed to allocate mcIndices\n" ) , exit( 0 );
	memset( mcIndices , 0 , sizeof(unsigned char) * _sNodes.nodeCount[_sNodes.maxDepth] );

	rootData.boundaryValues = new hash_map< long long , std::pair< Real , Point3D< Real > > >();
	int offSet = 0;
	int maxCCount = _sNodes.getMaxSlabCornerCount( sDepth , maxDepth , threads );
	int maxECount = _sNodes.getMaxSlabEdgeCount  ( &tree , sDepth , threads );
	rootData.cornerValues     = NewPointer< Real            >( maxCCount );
	rootData.cornerNormals    = NewPointer< Point3D< Real > >( maxCCount );
	rootData.interiorRoots    = NewPointer< int             >( maxECount );
	rootData.cornerValuesSet  = NewPointer< char            >( maxCCount );
	rootData.cornerNormalsSet = NewPointer< char            >( maxCCount );
	rootData.edgesSet         = NewPointer< char            >( maxECount );
	_sNodes.setCornerTable( coarseRootData , 0 , 0 , sDepth , threads );
	coarseRootData.cornerValues     = NewPointer< Real            >( coarseRootData.cCount );
	coarseRootData.cornerNormals    = NewPointer< Point3D< Real > >( coarseRootData.cCount );
	coarseRootData.cornerValuesSet  = NewPointer< char            >( coarseRootData.cCount );
	coarseRootData.cornerNormalsSet = NewPointer< char            >( coarseRootData.cCount );
	memset( coarseRootData.cornerValuesSet  , 0 , sizeof( char ) * coarseRootData.cCount );
	memset( coarseRootData.cornerNormalsSet , 0 , sizeof( char ) * coarseRootData.cCount );
	MemoryUsage();

	typename TreeOctNode::ConstNeighborKey3 nKey;
	nKey.set( maxDepth );
	// First process all leaf nodes at depths strictly finer than sDepth, one subtree at a time.
	std::vector< CornerValueStencil > vStencils( maxDepth+1 );
	std::vector< CornerNormalStencil > nStencils( maxDepth+1 );
	for( int d=_minDepth ; d<=maxDepth ; d++ )
	{
		SetCornerEvaluationStencil ( evaluator , d , vStencils[d].stencil  );
		SetCornerEvaluationStencils( evaluator , d , vStencils[d].stencils );
		SetCornerNormalEvaluationStencil ( evaluator , d , nStencils[d].stencil  );
		SetCornerNormalEvaluationStencils( evaluator , d , nStencils[d].stencils );
	}
	for( int i=0 ; i<SortedTreeNodes< Real >::Slices( sDepth ) ; i++ )
	{
		_sNodes.setCornerTable( rootData , sDepth , i , threads );
		_sNodes.setEdgeTable  ( rootData , sDepth , i , threads );
		memset( rootData.cornerValuesSet  , 0 , sizeof( char ) * rootData.cCount );
		memset( rootData.cornerNormalsSet , 0 , sizeof( char ) * rootData.cCount );
		memset( rootData.edgesSet         , 0 , sizeof( char ) * rootData.eCount );
		interiorVertices = new std::vector< Vertex >();
		for( int d=maxDepth ; d>sDepth && d>=_minDepth ; d-- )
		{
			std::pair< int , int > span = _sNodes.sliceSpan( sDepth , i , d );
			size_t leafNodeCount = 0;
			for( int i=span.first ; i<span.second ; i++ ) if( _sNodes.treeNodes[i]->nodeData.nodeIndex!=-1 && !_sNodes.treeNodes[i]->children ) leafNodeCount++;
			if( !leafNodeCount ) continue;
			std::vector< TreeOctNode* > leafNodes;
			leafNodes.reserve( leafNodeCount );
			for( int i=span.first ; i<span.second ; i++ ) if( _sNodes.treeNodes[i]->nodeData.nodeIndex!=-1 && !_sNodes.treeNodes[i]->children ) leafNodes.push_back( _sNodes.treeNodes[i] );

			// First set the corner values and associated marching-cube indices
#pragma omp parallel for num_threads( threads ) firstprivate( nKey )
			for( int i=0 ; i<leafNodeCount ; i++ )
			{
				TreeOctNode* leaf = leafNodes[i];
				SetIsoCorners( solution , isoValue , mcIndices , leaf , rootData , rootData.cornerValuesSet , rootData.cornerValues , nKey , GetPointer( metSolution ) , evaluator , vStencils[d].stencil , vStencils[d].stencils );

				// If this node shares a vertex with a node at the subdivision level, set the vertex value
				int d , off[3];
				leaf->depthAndOffset( d , off );
				int res = 1<<(d-sDepth);
				off[0] %= res , off[1] %=res , off[2] %= res;
				res--;
				if( !(off[0]%res) && !(off[1]%res) && !(off[2]%res) )
				{
					const TreeOctNode* temp = leaf;
					while( temp->depth()!=sDepth ) temp = temp->parent;
					int x = off[0]==0 ? 0 : 1 , y = off[1]==0 ? 0 : 1 , z = off[2]==0 ? 0 : 1;
					int c = Cube::CornerIndex( x , y , z );
					int idx = coarseRootData.cornerIndices( temp )[ c ];
					coarseRootData.cornerValues[ idx ] = rootData.cornerValues[ rootData.cornerIndices( leaf )[c] ];
					coarseRootData.cornerValuesSet[ idx ] = true;
				}

				// Compute the iso-vertices
				if( _boundaryType!=0 || _IsInset( leaf ) ) SetMCRootPositions( kernelDensityWeights , mcIndices , leaf , sDepth , isoValue , nKey , rootData , interiorVertices , mesh , solution , GetPointer( metSolution ) , evaluator , nStencils[d].stencil , nStencils[d].stencils , nonLinearFit );
			}
			// Note that this should be broken off for multi-threading as
			// the SetMCRootPositions writes to interiorPoints (with locking)
			// while GetMCIsoTriangles reads from interiorPoints (without locking)
			std::vector< Vertex > barycenters;
			std::vector< Vertex >* barycenterPtr = addBarycenter ? & barycenters : NULL;
#pragma omp parallel for num_threads( threads ) firstprivate( nKey )
			for( int i=0 ; i<leafNodeCount ; i++ )
			{
				TreeOctNode* leaf = leafNodes[i];
				if( _boundaryType!=0 || _IsInset( leaf ) ) GetMCIsoTriangles( mcIndices , leaf , nKey , mesh , rootData , interiorVertices , offSet , sDepth , polygonMesh , barycenterPtr );
			}
			for( int i=0 ; i<barycenters.size() ; i++ ) interiorVertices->push_back( barycenters[i] );
		}
		offSet = mesh->outOfCorePointCount();
		delete interiorVertices;
	}

	MemoryUsage();
	DeletePointer( rootData.cornerValues ) ; DeletePointer( rootData.cornerNormals );
	DeletePointer( rootData.cornerValuesSet ) ; DeletePointer( rootData.cornerNormalsSet );
	DeletePointer( rootData.interiorRoots );
	DeletePointer( rootData.edgesSet );
	coarseRootData.interiorRoots = NullPointer< int >();
	coarseRootData.boundaryValues = rootData.boundaryValues;
	for( hash_map< long long , int >::iterator iter=rootData.boundaryRoots.begin() ; iter!=rootData.boundaryRoots.end() ; iter++ ) 
		coarseRootData.boundaryRoots[iter->first] = iter->second;

	for( int d=sDepth ; d>=_minDepth ; d-- )
	{
		std::vector< Vertex > barycenters;
		std::vector< Vertex >* barycenterPtr = addBarycenter ? &barycenters : NULL;
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			TreeOctNode* leaf = _sNodes.treeNodes[i];
			if( leaf->children ) continue;
			// First set the corner values and associated marching-cube indices
			SetIsoCorners( solution , isoValue , mcIndices , leaf , coarseRootData , coarseRootData.cornerValuesSet , coarseRootData.cornerValues , nKey , GetPointer( metSolution ) , evaluator , vStencils[d].stencil , vStencils[d].stencils );

			// Now compute the iso-vertices
			{
				if( _boundaryType!=0 || _IsInset( leaf ) )
				{
					SetMCRootPositions< Vertex >( kernelDensityWeights , mcIndices , leaf , 0 , isoValue , nKey , coarseRootData , NULL , mesh , solution , GetPointer( metSolution ) , evaluator , nStencils[d].stencil , nStencils[d].stencils , nonLinearFit );
					GetMCIsoTriangles< Vertex >( mcIndices , leaf , nKey , mesh , coarseRootData , NULL , 0 , 0 , polygonMesh , barycenterPtr );
				}
			}
		}
	}
	MemoryUsage();
	DeletePointer( mcIndices );
	DeletePointer( coarseRootData.cornerValues ) ; DeletePointer( coarseRootData.cornerNormals );
	DeletePointer( coarseRootData.cornerValuesSet ) ; DeletePointer( coarseRootData.cornerNormalsSet );

	delete rootData.boundaryValues;
}
template< class Real , int Degree >
Real Octree< Real , Degree >::getCenterValue( const typename TreeOctNode::ConstNeighborKey3& neighborKey , const TreeOctNode* node , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CenterEvaluator< 1 >& evaluator , const Stencil< double , 3 >& stencil , const Stencil< double , 3 >& pStencil , bool isInterior ) const
{
	if( node->children ) fprintf( stderr , "[WARNING] getCenterValue assumes leaf node\n" );
	Real value=0;

	int d , off[3];
	node->depthAndOffset( d , off );

	if( isInterior )
	{
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
		{
			const TreeOctNode* n = neighborKey.neighbors[d].neighbors[i][j][k];
			if( n ) value += solution[ n->nodeData.nodeIndex ] * Real( stencil.values[i][j][k] );
		}
		if( d>_minDepth )
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
			{
				const TreeOctNode* n = neighborKey.neighbors[d-1].neighbors[i][j][k];
				if( n ) value += metSolution[n->nodeData.nodeIndex] * Real( pStencil.values[i][j][k] );
			}
	}
	else
	{
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
		{
			const TreeOctNode* n = neighborKey.neighbors[d].neighbors[i][j][k];
			if( n )
			{
				int _d , _off[3];
				n->depthAndOffset( _d , _off );
				value +=
					solution[ n->nodeData.nodeIndex ] * Real(
					evaluator.value( d , off[0] , _off[0] , false , false ) * evaluator.value( d , off[1] , _off[1] , false , false ) * evaluator.value( d , off[1] , _off[1] , false , false ) );
			}
		}
		if( d>_minDepth )
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
			{
				const TreeOctNode* n = neighborKey.neighbors[d-1].neighbors[i][j][k];
				if( n )
				{
					int _d , _off[3];
					n->depthAndOffset( _d , _off );
					value +=
						solution[ n->nodeData.nodeIndex ] * Real(
						evaluator.value( d , off[0] , _off[0] , false , false ) * evaluator.value( d , off[1] , _off[1] , false , false ) * evaluator.value( d , off[1] , _off[1] , false , false ) );
				}
			}
	}
	return value;
}
template< class Real , int Degree >
Real Octree< Real , Degree >::getCornerValue( const typename TreeOctNode::ConstNeighborKey3& neighborKey3 , const TreeOctNode* node , int corner , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< double , 3 >& stencil , const Stencil< double , 3 > stencils[8] , bool isInterior ) const
{
	double value = 0;
	if( _boundaryType==-1 ) value = -0.5;
	int d , off[3];
	node->depthAndOffset( d , off );

	int cx , cy , cz;
	int startX = 0 , endX = 3 , startY = 0 , endY = 3 , startZ = 0 , endZ = 3;
	Cube::FactorCornerIndex( corner , cx , cy , cz );
	{
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey3.neighbors[d];
		if( cx==0 ) endX = 2;
		else      startX = 1;
		if( cy==0 ) endY = 2;
		else      startY = 1;
		if( cz==0 ) endZ = 2;
		else      startZ = 1;
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( _node ) value += solution[ _node->nodeData.nodeIndex ] * stencil.values[x][y][z];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					value += solution[ _node->nodeData.nodeIndex ] * evaluator.value( d , off[0] , cx , _off[0] , false , false ) * evaluator.value( d , off[1] , cy , _off[1] , false , false ) * evaluator.value( d , off[2] , cz , _off[2] , false , false );
				}
			}
	}
	if( d>_minDepth )
	{
		int _corner = int( node - node->parent->children );
		int _cx , _cy , _cz;
		Cube::FactorCornerIndex( _corner , _cx , _cy , _cz );
		if( cx!=_cx ) startX = 0 , endX = 3;
		if( cy!=_cy ) startY = 0 , endY = 3;
		if( cz!=_cz ) startZ = 0 , endZ = 3;
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey3.neighbors[d-1];
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( _node ) value += metSolution[ _node->nodeData.nodeIndex ] * stencils[_corner].values[x][y][z];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					value += metSolution[ _node->nodeData.nodeIndex ] * evaluator.value( d , off[0] , cx , _off[0] , false , true ) * evaluator.value( d , off[1] , cy , _off[1] , false , true ) * evaluator.value( d , off[2] , cz , _off[2] , false , true );
				}
			}
	}
	return Real( value );
}
template< class Real , int Degree >
Point3D< Real > Octree< Real , Degree >::getCornerNormal( const typename TreeOctNode::ConstNeighbors5& neighbors5 , const typename TreeOctNode::ConstNeighbors5& pNeighbors5 , const TreeOctNode* node , int corner , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< Point3D< double > , 5 >& nStencil , const Stencil< Point3D< double > , 5 > nStencils[8] , bool isInterior ) const
{
	Point3D< double > normal;
	normal[0] = normal[1] = normal[2] = 0.;

	int d , off[3];
	node->depthAndOffset( d , off );

	int cx , cy , cz;
	int startX = 0 , endX = 5 , startY = 0 , endY = 5 , startZ = 0 , endZ = 5;
	Cube::FactorCornerIndex( corner , cx , cy , cz );
	{
		if( cx==0 ) endX = 4;
		else      startX = 1;
		if( cy==0 ) endY = 4;
		else      startY = 1;
		if( cz==0 ) endZ = 4;
		else      startZ = 1;
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors5.neighbors[x][y][z];
				if( _node ) normal += nStencil.values[x][y][z] * solution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					double v [] = { evaluator.value( d , off[0] , cx , _off[0] , false , false ) , evaluator.value( d , off[1] , cy , _off[1] , false , false ) , evaluator.value( d , off[2] , cz , _off[2] , false , false ) };
					double dv[] = { evaluator.value( d , off[0] , cx , _off[0] , true  , false ) , evaluator.value( d , off[1] , cy , _off[1] , true  , false ) , evaluator.value( d , off[2] , cz , _off[2] , true  , false ) };
					normal += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * solution[ _node->nodeData.nodeIndex ];
				}
			}
	}
	if( d>_minDepth )
	{
		int _cx , _cy , _cz , _corner = int( node - node->parent->children );
		Cube::FactorCornerIndex( _corner , _cx , _cy , _cz );
		if( cx!=_cx ) startX = 0 , endX = 5;
		if( cy!=_cy ) startY = 0 , endY = 5;
		if( cz!=_cz ) startZ = 0 , endZ = 5;
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=pNeighbors5.neighbors[x][y][z];
				if( _node ) normal += nStencils[_corner].values[x][y][z] * metSolution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					double v [] = { evaluator.value( d , off[0] , cx , _off[0] , false , true ) , evaluator.value( d , off[1] , cy , _off[1] , false , true ) , evaluator.value( d , off[2] , cz , _off[2] , false , true ) };
					double dv[] = { evaluator.value( d , off[0] , cx , _off[0] , true  , true ) , evaluator.value( d , off[1] , cy , _off[1] , true  , true ) , evaluator.value( d , off[2] , cz , _off[2] , true  , true ) };
					normal += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * metSolution[ _node->nodeData.nodeIndex ];
				}
			}
	}
	return Point3D< Real >( Real(normal[0]) , Real(normal[1]) , Real(normal[2]) );
}
template< class Real , int Degree >
Real Octree< Real , Degree >::GetIsoValue( ConstPointer( Real ) solution , const std::vector< Real >& centerWeights )
{
	Real isoValue=0 , weightSum=0;
	int maxDepth = tree.maxDepth();

	typename BSplineData< Degree >::template CenterEvaluator< 1 > evaluator;
	_fData.setCenterEvaluator( evaluator , 0 , 0 , _boundaryType==0 );
	std::vector< CenterValueStencil > vStencils( maxDepth+1 );
	for( int d=_minDepth ; d<=maxDepth ; d++ )
	{
		SetCenterEvaluationStencil ( evaluator , d , vStencils[d].stencil  );
		SetCenterEvaluationStencils( evaluator , d , vStencils[d].stencils );
	}
	std::vector< Real > metSolution( _sNodes.nodeCount[maxDepth] , 0 );
	std::vector< Real > centerValues( _sNodes.nodeCount[maxDepth+1] );
#pragma omp parallel for num_threads( threads )
	for( int i=_sNodes.nodeCount[_minDepth] ; i<_sNodes.nodeCount[maxDepth] ; i++ ) metSolution[i] = solution[i];
	for( int d=_minDepth ; d<maxDepth ; d++ ) UpSample( d , _sNodes , ( ConstPointer( Real ) )GetPointer( metSolution ) + _sNodes.nodeCount[d-1] , GetPointer( metSolution ) + _sNodes.nodeCount[d] );
	for( int d=maxDepth ; d>=_minDepth ; d-- )
	{
		typename TreeOctNode::ConstNeighborKey3 nKey;
		nKey.set( d );
#pragma omp parallel for num_threads( threads ) reduction( + : isoValue , weightSum ) firstprivate( nKey )
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			TreeOctNode* node = _sNodes.treeNodes[i];
			Real value = Real(0);
			if( node->children )
			{
				for( int c=0 ; c<Cube::CORNERS ; c++ ) value += centerValues[ node->children[c].nodeData.nodeIndex ];
				value /= Cube::CORNERS;
			}
			else
			{
				nKey.getNeighbors( node );
				int c=0 , x , y , z;
				if( node->parent ) c = int( node - node->parent->children );
				Cube::FactorCornerIndex( c , x , y , z );

				int d , off[3];
				node->depthAndOffset( d , off );
				int o = _boundaryType==0 ? (1<<(d-2)) : 0;
				int mn = 2+o , mx = (1<<d)-2-o;
				bool isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );

				value = getCenterValue( nKey , node , solution , GetPointer( metSolution ) , evaluator , vStencils[d].stencil , vStencils[d].stencils[c] , isInterior );
			}
			centerValues[i] = value;
			Real w = centerWeights[ node->nodeData.nodeIndex ];
			if( w!=0 ) isoValue += value * w , weightSum += w;
		}
	}
	if( _boundaryType==-1 ) return isoValue/weightSum - Real(0.5);
	else                    return isoValue/weightSum;
}

template< class Real , int Degree >
void Octree< Real , Degree >::SetIsoCorners( ConstPointer( Real ) solution , Real isoValue , Pointer( unsigned char ) mcIndices , TreeOctNode* leaf , typename SortedTreeNodes< Real >::SlabCornerTableData& cData , Pointer( char ) valuesSet , Pointer( Real ) values , typename TreeOctNode::ConstNeighborKey3& nKey , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< double , 3 > stencil[8] , const Stencil< double , 3 > stencils[8][8] )
{
	Real cornerValues[ Cube::CORNERS ];
	const typename SortedTreeNodes< Real >::CubeCornerIndices& cIndices = cData[ leaf ];

	bool isInterior;
	int d , off[3];
	leaf->depthAndOffset( d , off );
	int o = _boundaryType==0 ? (1<<(d-2)) : 0;
	int mn = 2+o , mx = (1<<d)-2-o;
	isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
	nKey.getNeighbors( leaf );
	for( int c=0 ; c<Cube::CORNERS ; c++ )
	{
		int vIndex = cIndices[c];
		if( valuesSet[vIndex] ) cornerValues[c] = values[vIndex];
		else
		{
			cornerValues[c] = getCornerValue( nKey , leaf , c , solution , metSolution , evaluator , stencil[c] , stencils[c] , isInterior );
			values[vIndex] = cornerValues[c];
			valuesSet[vIndex] = 1;
		}
	}
	mcIndices[ leaf->nodeData.nodeIndex ] = MarchingCubes::GetIndex( cornerValues , isoValue );

}
template< class Real , int Degree >
bool Octree< Real , Degree >::IsBoundaryCorner( const TreeOctNode* node , int cornerIndex , int subdivideDepth )
{
	if( subdivideDepth<0 ) return false;
	int d , off[3] , x , y , z;
	node->depthAndOffset( d , off );
	if( d<=subdivideDepth ) return true;
	Cube::FactorCornerIndex( cornerIndex , x , y , z );
	off[0] += x , off[1] += y , off[2] += z;
	int mask = 1<<( d - subdivideDepth );
	return !(off[0]%mask) && !(off[1]%mask) && !(off[2]%mask);
}
template< class Real , int Degree >
bool Octree< Real , Degree >::IsBoundaryFace( const TreeOctNode* node , int faceIndex , int subdivideDepth )
{
	int dir,offset,d,o[3],idx;

	if( subdivideDepth<0 ) return false;
	if( node->depth()<=subdivideDepth ) return true;
	Cube::FactorFaceIndex( faceIndex , dir , offset );
	if( dir!=2 ) return 0;
	node->depthAndOffset(d,o);

	idx=(int(o[dir])<<1) + (offset<<1);
	return !( idx%(2<<(d-subdivideDepth)) );
}
template< class Real , int Degree >
bool Octree< Real , Degree >::IsBoundaryEdge( const TreeOctNode* node , int edgeIndex , int subdivideDepth )
{
	int dir , x , y;
	Cube::FactorEdgeIndex( edgeIndex , dir , x , y );
	return IsBoundaryEdge( node , dir , x , y , subdivideDepth );
}
template< class Real , int Degree >
bool Octree< Real , Degree >::IsBoundaryEdge( const TreeOctNode* node , int dir , int x , int y , int subdivideDepth )
{

	if( subdivideDepth<0 ) return 0;
	if( node->depth()<=subdivideDepth ) return 1;
	int d , o[3] , idx1 , idx2 , mask;
	node->depthAndOffset( d , o );

	switch( dir )
	{
		case 0:
			idx1 = o[1] + x;
			idx2 = o[2] + y;
			break;
		case 1:
			idx1 = o[0] + x;
			idx2 = o[2] + y;
			break;
		case 2:
			idx1 = o[0] + x;
			idx2 = o[1] + y;
			break;
	}
	mask = 1<<( d - subdivideDepth );
	return !(idx1%(mask)) || !(idx2%(mask));
}
template< class Real , int Degree >
void Octree< Real , Degree >::GetRootSpan( const RootInfo< Real >& ri , Point3D< Real >& start , Point3D< Real >& end )
{
	int o , i1 , i2;
	Real width;
	Point3D< Real > c;

	Cube::FactorEdgeIndex( ri.edgeIndex , o , i1 , i2 );
	ri.node->centerAndWidth( c , width );
	switch(o)
	{
	case 0:
		start[0]          = c[0] - width/2;
		end  [0]          = c[0] + width/2;
		start[1] = end[1] = c[1] - width/2 + width*i1;
		start[2] = end[2] = c[2] - width/2 + width*i2;
		break;
	case 1:
		start[0] = end[0] = c[0] - width/2 + width*i1;
		start[1]          = c[1] - width/2;
		end  [1]          = c[1] + width/2;
		start[2] = end[2] = c[2] - width/2 + width*i2;
		break;
	case 2:
		start[0] = end[0] = c[0] - width/2 + width*i1;
		start[1] = end[1] = c[1] - width/2 + width*i2;
		start[2]          = c[2] - width/2;
		end  [2]          = c[2] + width/2;
		break;
	}
}
//////////////////////////////////////////////////////////////////////////////////////
// The assumption made when calling this code is that the edge has at most one root //
//////////////////////////////////////////////////////////////////////////////////////
template< class Real > void SetVertexValue(      PlyVertex< float >& vertex , Real value ){ ; }
template< class Real > void SetVertexValue( PlyValueVertex< float >& vertex , Real value ){ vertex.value = float(value); }
template< class Real , int Degree >
template< class Vertex >
int Octree< Real , Degree >::GetRoot( const std::vector< Real >* kernelDensityWeights , const RootInfo< Real >& ri , Real isoValue , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , Vertex& vertex , RootData& rootData , int sDepth , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< Point3D< double > , 5 > nStencil[8] , const Stencil< Point3D< double > , 5 > nStencils[8][8] , int nonLinearFit )
{
	Point3D< Real > position;
	int c1 , c2;
	Cube::EdgeCorners( ri.edgeIndex , c1 , c2 );

	long long key1 , key2;
	Point3D< Real > n[2];

	int i , o , i1 , i2 , rCount=0;
	Polynomial<2> P;
	double roots[2];
	double x0 , x1;
	Real center , width;
	Real averageRoot=0;
	Cube::FactorEdgeIndex( ri.edgeIndex , o , i1 , i2 );
	int idx1[3] , idx2[3];
	key1 = VertexData< Real >::CornerIndex( ri.node , c1 , _fData.depth , idx1 );
	key2 = VertexData< Real >::CornerIndex( ri.node , c2 , _fData.depth , idx2 );

	bool isBoundary1 = IsBoundaryCorner( ri.node , c1 , sDepth );
	bool isBoundary2 = IsBoundaryCorner( ri.node , c2 , sDepth );
	bool haveKey1 , haveKey2;
	std::pair< Real , Point3D< Real > > keyValue1 , keyValue2;
	int iter1 , iter2;
	{
		iter1 = rootData.cornerIndices( ri.node )[ c1 ];
		iter2 = rootData.cornerIndices( ri.node )[ c2 ];
		keyValue1.first = rootData.cornerValues[iter1];
		keyValue2.first = rootData.cornerValues[iter2];
		if( isBoundary1 )
		{
#pragma omp critical (normal_hash_access)
			{
				haveKey1 = ( rootData.boundaryValues->find( key1 )!=rootData.boundaryValues->end() );
				if( haveKey1 ) keyValue1 = (*rootData.boundaryValues)[key1];
			}
		}
		else
		{
			haveKey1 = ( rootData.cornerNormalsSet[iter1] != 0 );
			keyValue1.first = rootData.cornerValues[iter1];
			if( haveKey1 ) keyValue1.second = rootData.cornerNormals[iter1];
		}
		if( isBoundary2 )
		{
#pragma omp critical (normal_hash_access)
			{
				haveKey2 = ( rootData.boundaryValues->find( key2 )!=rootData.boundaryValues->end() );
				if( haveKey2 ) keyValue2 = (*rootData.boundaryValues)[key2];
			}
		}
		else
		{
			haveKey2 = ( rootData.cornerNormalsSet[iter2] != 0 );
			keyValue2.first = rootData.cornerValues[iter2];
			if( haveKey2 ) keyValue2.second = rootData.cornerNormals[iter2];
		}
	}
	typename TreeOctNode::ConstNeighbors5 neighbors5 , pNeighbors5;
	if( !haveKey1 || !haveKey2 )
	{
		neighborKey3.getNeighbors( ri.node , neighbors5 );
		if( ri.node->parent ) neighborKey3.getNeighbors( ri.node->parent , pNeighbors5 );
	}
	bool isInterior;
	if( !haveKey1 || !haveKey2 )
	{
		int d , off[3];
		ri.node->depthAndOffset( d , off );
		int o = _boundaryType==0 ? (1<<(d-2)) : 0;
		int mn = 2+o , mx = (1<<d)-2-o;
		isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
	}
	if( !haveKey1 ) keyValue1.second = getCornerNormal( neighbors5 , pNeighbors5 , ri.node , c1 , solution , metSolution , evaluator , nStencil[c1] , nStencils[c1] , isInterior );
	x0 = keyValue1.first;
	n[0] = keyValue1.second;

	if( !haveKey2 ) keyValue2.second = getCornerNormal( neighbors5 , pNeighbors5 , ri.node , c2 , solution , metSolution , evaluator , nStencil[c2] , nStencils[c2] , isInterior );
	x1 = keyValue2.first;
	n[1] = keyValue2.second;

	if( !haveKey1 )
		if( isBoundary1 )
		{
#pragma omp critical (normal_hash_access)
			{
				if( !haveKey1 ) (*rootData.boundaryValues)[key1] = keyValue1;
			}
		}
		else if( !haveKey1 ) rootData.cornerNormals[iter1] = keyValue1.second , rootData.cornerNormalsSet[iter1] = 1;
	if( !haveKey2 )
		if( isBoundary2 )
		{
#pragma omp critical (normal_hash_access)
			{
				if( !haveKey2 ) (*rootData.boundaryValues)[key2] = keyValue2;
			}
		}
		else if( !haveKey2 ) rootData.cornerNormals[iter2] = keyValue2.second , rootData.cornerNormalsSet[iter2] = 1;

	Point3D< Real > c;
	ri.node->centerAndWidth(c,width);
	center=c[o];
	for( i=0 ; i<DIMENSION ; i++ ) n[0][i] *= width , n[1][i] *= width;

	switch(o)
	{
	case 0:
		position[1] = c[1]-width/2+width*i1;
		position[2] = c[2]-width/2+width*i2;
		break;
	case 1:
		position[0] = c[0]-width/2+width*i1;
		position[2] = c[2]-width/2+width*i2;
		break;
	case 2:
		position[0] = c[0]-width/2+width*i1;
		position[1] = c[1]-width/2+width*i2;
		break;
	}
	double dx0,dx1;
	dx0 = n[0][o];
	dx1 = n[1][o];

	// The scaling will turn the Hermite Spline into a quadratic
	double scl=(x1-x0)/((dx1+dx0)/2);
	dx0 *= scl;
	dx1 *= scl;

	// Hermite Spline
	P.coefficients[0] = x0;
	P.coefficients[1] = dx0;
	P.coefficients[2] = 3*(x1-x0)-dx1-2*dx0;

	int rootCount = P.getSolutions( isoValue , roots , EPSILON );
	for( i=0 ; i<rootCount ; i++ )
		if( roots[i]>=0 && roots[i]<=1 )
		{
			averageRoot += Real( roots[i] );
			rCount++;
		}
	if( rCount && nonLinearFit ) averageRoot /= rCount;
	else					     averageRoot  = Real((x0-isoValue)/(x0-x1));
	if( averageRoot<0 || averageRoot>1 )
	{
		fprintf( stderr , "[WARNING] Bad average root: %f\n" , averageRoot );
		fprintf( stderr , "\t(%f %f) , (%f %f) (%f)\n" , x0 , x1 , dx0 , dx1 , isoValue );
		if( averageRoot<0 ) averageRoot = 0;
		if( averageRoot>1 ) averageRoot = 1;
	}
	position[o] = Real(center-width/2+width*averageRoot);
	vertex.point = position;
	if( kernelDensityWeights )
	{
		Real depth , weight;
		const TreeOctNode* temp = ri.node;
		while( temp->depth()>_splatDepth ) temp=temp->parent;
		GetSampleDepthAndWeight( *kernelDensityWeights , temp , position , neighborKey3 , _samplesPerNode , depth , weight );
		SetVertexValue( vertex , depth );
	}
	return 1;
}
template< class Real , int Degree >
int Octree< Real , Degree >::GetRootIndex( ConstPointer( unsigned char ) mcIndices , const TreeOctNode* node , int edgeIndex , int maxDepth , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , RootInfo< Real >& ri )
{
	int c1 , c2 , f1 , f2;
	const TreeOctNode *temp , *finest;
	int finestIndex;

	// The assumption is that the super-edge has a root along it.
	int mcIndex = mcIndices[ node->nodeData.nodeIndex ];
	if( !(MarchingCubes::edgeMask[mcIndex] & (1<<edgeIndex) ) ) return 0;

	Cube::FacesAdjacentToEdge( edgeIndex , f1 , f2 );

	finest = node;
	finestIndex = edgeIndex;
	if( node->depth()<maxDepth && !node->children )
	{
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey3.getNeighbors( node );
		int x , y , z;
		Cube::FactorFaceIndex( f1 , x , y , z );
		temp = neighbors.neighbors[x+1][y+1][z+1];
		if( temp && temp->nodeData.nodeIndex!=-1 && temp->children ) finest = temp , finestIndex = Cube::FaceReflectEdgeIndex( edgeIndex , f1 );
		else
		{
			int x , y , z;
			Cube::FactorFaceIndex( f2 , x , y , z );
			temp = neighbors.neighbors[x+1][y+1][z+1];
			if( temp && temp->nodeData.nodeIndex!=-1 && temp->children ) finest = temp , finestIndex = Cube::FaceReflectEdgeIndex( edgeIndex , f2 );
			else
			{
				int orientation , d1 , d2;
				Cube::FactorEdgeIndex( edgeIndex , orientation , d1 , d2 );
				switch( orientation )
				{
				case 0: temp = neighbors.neighbors[1][d1<<1][d2<<1] ; break;
				case 1: temp = neighbors.neighbors[d1<<1][1][d2<<1] ; break;
				case 2: temp = neighbors.neighbors[d1<<1][d2<<1][1] ; break;
				}
				if( temp && temp->nodeData.nodeIndex!=-1 && temp->children ) finest = temp , finestIndex = Cube::EdgeReflectEdgeIndex( edgeIndex );
			}
		}
	}

	Cube::EdgeCorners( finestIndex , c1 , c2 );
	if( finest->children )
	{
		if     ( GetRootIndex( mcIndices , finest->children + c1 , finestIndex , maxDepth , neighborKey3 , ri ) ) return 1;
		else if( GetRootIndex( mcIndices , finest->children + c2 , finestIndex , maxDepth , neighborKey3 , ri ) ) return 1;
		else
		{
			int d1 , off1[3] , d2 , off2[3];
			node->depthAndOffset( d1 , off1 );
			finest->depthAndOffset( d2 , off2 );
			int c1 , c2;
			Cube::EdgeCorners( edgeIndex , c1 , c2 );
			int x1 , x2 , y1 , y2 , z1 , z2;
			Cube::FactorCornerIndex( c1 , x1 , y1 , z1 );
			Cube::FactorCornerIndex( c2 , x2 , y2 , z2 );
			fprintf( stderr , "[WARNING] Couldn't find root index with either child [%d] (%d %d %d) { %d %d %d - %d %d %d } -> [%d] (%d %d %d) (%d %d)\n" , d1 , off1[0] , off1[1] , off1[2] , x1 , y1 , z1 , x2 , y2 , z2 , d2 , off2[0] , off2[1] , off2[2] , node->children!=NULL , finest->children!=NULL );
			printf( "\t" );
			for( int i=0 ; i<8 ; i++ ) if( mcIndices[ node->nodeData.nodeIndex ] & (1<<i) ) printf( "1" ); else printf( "0" );
			printf( "\t" );
			for( int i=0 ; i<8 ; i++ ) if( mcIndices[ finest->nodeData.nodeIndex ] & (1<<i) ) printf( "1" ); else printf( "0" );
			printf( "\n" );
			return 0;
		}
	}
	else
	{
		int o,i1,i2;
		Cube::FactorEdgeIndex(finestIndex,o,i1,i2);
		int d,off[3];
		finest->depthAndOffset(d,off);
		ri.node=finest;
		ri.edgeIndex=finestIndex;
		int offset,eIndex[2];
		offset = BinaryNode< Real >::CenterIndex( d , off[o] );
		switch(o)
		{
		case 0:
			eIndex[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[1],i1);
			eIndex[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[2],i2);
			break;
		case 1:
			eIndex[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[0],i1);
			eIndex[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[2],i2);
			break;
		case 2:
			eIndex[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[0],i1);
			eIndex[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[1],i2);
			break;
		}
		ri.key= (long long)(o) | (long long)(eIndex[0])<<5 | (long long)(eIndex[1])<<25 | (long long)(offset)<<45;
		return 1;
	}
}
template< class Real , int Degree >
int Octree< Real , Degree >::GetRootPair( ConstPointer( unsigned char ) mcIndices , const RootInfo< Real >& ri , int maxDepth , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , RootInfo< Real >& pair )
{
	const TreeOctNode* node = ri.node;
	int c1 , c2 , c;
	Cube::EdgeCorners( ri.edgeIndex , c1 , c2 );
	while( node->parent )
	{
		c = int(node-node->parent->children);
		if( c!=c1 && c!=c2 ) return 0;
		if( !MarchingCubes::HasEdgeRoots( mcIndices[ node->parent->nodeData.nodeIndex ] , ri.edgeIndex ) )
		{
			if( c==c1 ) return GetRootIndex( mcIndices , &node->parent->children[c2] , ri.edgeIndex , maxDepth , neighborKey3 , pair );
			else        return GetRootIndex( mcIndices , &node->parent->children[c1] , ri.edgeIndex , maxDepth , neighborKey3 , pair );
		}
		node = node->parent;
	}
	return 0;

}
template< class Real , int Degree >
int Octree< Real , Degree >::GetRootIndex( const RootInfo< Real >& ri , RootData& rootData , CoredPointIndex& index )
{
	long long key = ri.key;
	hash_map< long long , int >::iterator rootIter;
	rootIter = rootData.boundaryRoots.find( key );
	if( rootIter!=rootData.boundaryRoots.end() )
	{
		index.inCore = 1;
		index.index = rootIter->second;
		return 1;
	}
	else if( rootData.interiorRoots )
	{
		int eIndex = rootData.edgeIndices( ri.node )[ ri.edgeIndex ];
		if( rootData.edgesSet[ eIndex ] )
		{
			index.inCore = 0;
			index.index = rootData.interiorRoots[ eIndex ];
			return 1;
		}
	}
	return 0;
}
template< class Real , int Degree >
template< class Vertex >
int Octree< Real , Degree >::SetMCRootPositions( const std::vector< Real >* kernelDensityWeights , ConstPointer( unsigned char ) mcIndices , TreeOctNode* node , int sDepth , Real isoValue , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , RootData& rootData , 
	std::vector< Vertex >* interiorVertices , CoredMeshData< Vertex >* mesh , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< Degree >::template CornerEvaluator< 2 >& evaluator , const Stencil< Point3D< double > , 5 > nStencil[8] , const Stencil< Point3D< double > , 5 > nStencils[8][8] , int nonLinearFit )
{
	Vertex vertex;
	int eIndex;
	RootInfo< Real > ri;
	int count=0;
	if( !MarchingCubes::HasRoots( mcIndices[ node->nodeData.nodeIndex ] ) ) return 0;
	for( int i=0 ; i<DIMENSION ; i++ ) for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ )
	{
		long long key;
		eIndex = Cube::EdgeIndex( i , j , k );
		if( GetRootIndex( mcIndices , node , eIndex , _fData.depth , neighborKey3 , ri ) )
		{
			key = ri.key;
			if( !rootData.interiorRoots || IsBoundaryEdge( node , i , j , k , sDepth ) )
			{
				hash_map< long long , int >::iterator iter , end;
				// Check if the root has already been set
#pragma omp critical (boundary_roots_hash_access)
				{
					iter = rootData.boundaryRoots.find( key );
					end  = rootData.boundaryRoots.end();
				}
				if( iter==end )
				{
					// Get the root information
					GetRoot( kernelDensityWeights , ri , isoValue , neighborKey3 , vertex , rootData , sDepth , solution , metSolution , evaluator , nStencil , nStencils , nonLinearFit );
					vertex.point = vertex.point * _scale + _center;
					// Add the root if it hasn't been added already
#pragma omp critical (boundary_roots_hash_access)
					{
						iter = rootData.boundaryRoots.find( key );
						end  = rootData.boundaryRoots.end();
						if( iter==end )
						{
							mesh->inCorePoints.push_back( vertex );
							rootData.boundaryRoots[key] = int( mesh->inCorePoints.size() ) - 1;
						}
					}
					if( iter==end ) count++;
				}
			}
			else
			{
				int nodeEdgeIndex = rootData.edgeIndices( ri.node )[ ri.edgeIndex ];
				if( !rootData.edgesSet[ nodeEdgeIndex ] )
				{
					// Get the root information
					GetRoot( kernelDensityWeights , ri , isoValue , neighborKey3 , vertex , rootData , sDepth , solution , metSolution , evaluator , nStencil , nStencils , nonLinearFit );
					vertex.point = vertex.point * _scale + _center;
					// Add the root if it hasn't been added already
#pragma omp critical (add_point_access)
					{
						if( !rootData.edgesSet[ nodeEdgeIndex ] )
						{
							rootData.interiorRoots[ nodeEdgeIndex ] = mesh->addOutOfCorePoint_s( vertex );
							interiorVertices->push_back( vertex );
							rootData.edgesSet[ nodeEdgeIndex ] = 1;
							count++;
						}
					}
				}
			}
		}
	}
	return count;
}
template< class Real , int Degree >
void Octree< Real , Degree >::GetMCIsoEdges( ConstPointer( unsigned char ) mcIndices , TreeOctNode* node , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , int sDepth , std::vector< std::pair< RootInfo< Real > , RootInfo< Real > > >& edges )
{
	const TreeOctNode* temp;
	int count=0 , tris=0;
	int isoTri[ DIMENSION * MarchingCubes::MAX_TRIANGLES ];
	FaceEdgesFunction fef;
	int ref , fIndex;
	typename hash_map< long long , std::pair< RootInfo< Real > , int > >::iterator iter;
	hash_map< long long , std::pair< RootInfo< Real > , int > > vertexCount;

	fef.edges = &edges;
	fef.maxDepth = _fData.depth;
	fef.vertexCount = &vertexCount;
	fef.neighborKey3 = &neighborKey3;
	fef.mcIndices = mcIndices;
	const TreeOctNode* _temp[Cube::NEIGHBORS];
	{
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey3.getNeighbors( node );
		for( int f=0 ; f<Cube::NEIGHBORS ; f++ )
		{
			int x , y , z;
			Cube::FactorFaceIndex( f , x , y , z );
			_temp[f] = neighbors.neighbors[x+1][y+1][z+1];
		}
	}
	int mcIndex = mcIndices[ node->nodeData.nodeIndex ];
	count = MarchingCubes::AddTriangleIndices( mcIndex , isoTri );
	for( fIndex=0 ; fIndex<Cube::NEIGHBORS ; fIndex++ )
	{
		ref = Cube::FaceReflectFaceIndex( fIndex , fIndex );
		fef.fIndex = ref;
		temp = _temp[fIndex];
		// If the face neighbor exists and has higher resolution than the current node,
		// get the iso-curve from the neighbor
		if( temp && temp->nodeData.nodeIndex!=-1 && temp->children && !IsBoundaryFace( node , fIndex , sDepth ) ) temp->processNodeFaces( temp , &fef , ref );
		// Otherwise, get it from the node
		else
		{
			RootInfo< Real > ri1 , ri2;
			for( int j=0 ; j<count ; j++ )
				for( int k=0 ; k<3 ; k++ )
					if( fIndex==Cube::FaceAdjacentToEdges( isoTri[j*3+k] , isoTri[j*3+((k+1)%3)] ) )
						if( GetRootIndex( mcIndices , node , isoTri[j*3+k] , _fData.depth , neighborKey3 , ri1 ) && GetRootIndex( mcIndices , node , isoTri[j*3+((k+1)%3)] , _fData.depth , neighborKey3 , ri2 ) )
						{
							long long key1 = ri1.key , key2 = ri2.key;
							edges.push_back( std::pair< RootInfo< Real > , RootInfo< Real > >( ri1 , ri2 ) );
							iter=vertexCount.find( key1 );
							if( iter==vertexCount.end() )
							{
								vertexCount[key1].first = ri1;
								vertexCount[key1].second = 0;
							}
							iter=vertexCount.find( key2 );
							if( iter==vertexCount.end() )
							{
								vertexCount[key2].first = ri2;
								vertexCount[key2].second = 0;
							}
							vertexCount[key1].second++;
							vertexCount[key2].second--;
						}
						else
						{
							int r1 = MarchingCubes::HasEdgeRoots( mcIndex , isoTri[j*3+k] );
							int r2 = MarchingCubes::HasEdgeRoots( mcIndex , isoTri[j*3+((k+1)%3)] );
							fprintf( stderr , "Bad Edge 2: %d %d\t%d %d\n" , ri1.key , ri2.key , r1 , r2 );
						}
		}
	}
	for( int i=0 ; i<int(edges.size()) ; i++ )
	{
		iter = vertexCount.find( edges[i].first.key );
		if( iter==vertexCount.end() ) printf( "Could not find vertex: %lld\n" , edges[i].first );
		else if( vertexCount[ edges[i].first.key ].second )
		{
			RootInfo< Real > ri;
			GetRootPair( mcIndices , vertexCount[edges[i].first.key].first , _fData.depth , neighborKey3 , ri );
			long long key = ri.key;
			iter = vertexCount.find( key );
			if( iter==vertexCount.end() )
			{
				int d , off[3];
				node->depthAndOffset( d , off );
				printf( "Vertex pair not in list 1 (%lld) %d\t[%d] (%d %d %d)\n" , key , IsBoundaryEdge( ri.node , ri.edgeIndex , sDepth ) , d , off[0] , off[1] , off[2] );
			}
			else
			{
				edges.push_back( std::pair< RootInfo< Real > , RootInfo< Real > >( ri , edges[i].first ) );
				vertexCount[ key ].second++;
				vertexCount[ edges[i].first.key ].second--;
			}
		}

		iter = vertexCount.find( edges[i].second.key );
		if( iter==vertexCount.end() ) printf( "Could not find vertex: %lld\n" , edges[i].second );
		else if( vertexCount[edges[i].second.key].second )
		{
			RootInfo< Real > ri;
			GetRootPair( mcIndices , vertexCount[edges[i].second.key].first , _fData.depth , neighborKey3 , ri );
			long long key = ri.key;
			iter=vertexCount.find( key );
			if( iter==vertexCount.end() )
			{
				int d , off[3];
				node->depthAndOffset( d , off );
				printf( "Vertex pair not in list 2\t[%d] (%d %d %d)\n" , d , off[0] , off[1] , off[2] );
			}
			else
			{
				edges.push_back( std::pair< RootInfo< Real > , RootInfo< Real > >( edges[i].second , ri ) );
				vertexCount[key].second--;
				vertexCount[ edges[i].second.key ].second++;
			}
		}
	}
}
template< class Real , int Degree >
template< class Vertex >
int Octree< Real , Degree >::GetMCIsoTriangles( ConstPointer( unsigned char ) mcIndices , TreeOctNode* node , typename TreeOctNode::ConstNeighborKey3& neighborKey3 , CoredMeshData< Vertex >* mesh , RootData& rootData , std::vector< Vertex >* interiorVertices , int offSet , int sDepth , bool polygonMesh , std::vector< Vertex >* barycenters )
{
	int tris=0;
	std::vector< std::pair< RootInfo< Real > , RootInfo< Real > > > edges;
	std::vector< std::vector< std::pair< RootInfo< Real > , RootInfo< Real > > > > edgeLoops;
	GetMCIsoEdges( mcIndices , node , neighborKey3 , sDepth , edges );

	GetEdgeLoops( edges , edgeLoops );
	for( int i=0 ; i<int(edgeLoops.size()) ; i++ )
	{
		CoredPointIndex p;
		std::vector<CoredPointIndex> edgeIndices;
		for( int j=int(edgeLoops[i].size())-1 ; j>=0 ; j-- )
		{
			if( !GetRootIndex( edgeLoops[i][j].first , rootData , p ) ) printf( "Bad Point Index\n" );
			else edgeIndices.push_back( p );
		}
		tris += AddTriangles( mesh , edgeIndices , interiorVertices , offSet , polygonMesh , barycenters );
	}
	return tris;
}

template< class Real , int Degree >
int Octree< Real , Degree >::GetEdgeLoops( std::vector< std::pair< RootInfo< Real > , RootInfo< Real > > >& edges , std::vector< std::vector< std::pair< RootInfo< Real > , RootInfo< Real > > > >& loops )
{
	int loopSize=0;
	long long frontIdx , backIdx;
	std::pair< RootInfo< Real > , RootInfo< Real > > e , temp;
	loops.clear();

	while( edges.size() )
	{
		std::vector< std::pair< RootInfo< Real > ,  RootInfo< Real > > > front , back;
		e = edges[0];
		loops.resize( loopSize+1 );
		edges[0] = edges.back();
		edges.pop_back();
		frontIdx = e.second.key;
		backIdx = e.first.key;
		for( int j=int(edges.size())-1 ; j>=0 ; j-- )
		{
			if( edges[j].first.key==frontIdx || edges[j].second.key==frontIdx )
			{
				if( edges[j].first.key==frontIdx ) temp = edges[j];
				else temp.first = edges[j].second , temp.second = edges[j].first;
				frontIdx = temp.second.key;
				front.push_back(temp);
				edges[j] = edges.back();
				edges.pop_back();
				j = int(edges.size());
			}
			else if( edges[j].first.key==backIdx || edges[j].second.key==backIdx )
			{
				if( edges[j].second.key==backIdx ) temp = edges[j];
				else temp.first = edges[j].second , temp.second = edges[j].first;
				backIdx = temp.first.key;
				back.push_back(temp);
				edges[j] = edges.back();
				edges.pop_back();
				j = int(edges.size());
			}
		}
		for( int j=int(back.size())-1 ; j>=0 ; j-- ) loops[loopSize].push_back( back[j] );
		loops[loopSize].push_back(e);
		for( int j=0 ; j<int(front.size()) ; j++ ) loops[loopSize].push_back( front[j] );
		loopSize++;
	}
	return int(loops.size());
}
template< class Real , int Degree >
template< class Vertex >
int Octree< Real , Degree >::AddTriangles( CoredMeshData< Vertex >* mesh , std::vector< CoredPointIndex >& edges , std::vector< Vertex >* interiorVertices , int offSet , bool polygonMesh , std::vector< Vertex >* barycenters )
{
	MinimalAreaTriangulation< Real > MAT;
	std::vector< Point3D< Real > > vertices;
	std::vector< TriangleIndex > triangles;
	if( polygonMesh )
	{
		std::vector< CoredVertexIndex > vertices( edges.size() );
		for( int i=0 ; i<int(edges.size()) ; i++ )
		{
			vertices[i].idx    =  edges[i].index;
			vertices[i].inCore = (edges[i].inCore!=0);
		}
		mesh->addPolygon_s( vertices );
		return 1;
	}
	if( edges.size()>3 )
	{
		bool isCoplanar = false;

		if( barycenters )
			for( int i=0 ; i<int(edges.size()) ; i++ )
				for( int j=0 ; j<i ; j++ )
					if( (i+1)%edges.size()!=j && (j+1)%edges.size()!=i )
					{
						Vertex v1 , v2;
						if( edges[i].inCore ) v1 =  mesh->inCorePoints[ edges[i].index        ];
						else                  v1 = (*interiorVertices)[ edges[i].index-offSet ];
						if( edges[j].inCore ) v2 =  mesh->inCorePoints[ edges[j].index        ];
						else                  v2 = (*interiorVertices)[ edges[j].index-offSet ];
						for( int k=0 ; k<3 ; k++ ) if( v1.point[k]==v2.point[k] ) isCoplanar = true;
					}
		if( isCoplanar )
		{
			Vertex c;
			c *= 0;
			for( int i=0 ; i<int(edges.size()) ; i++ )
			{
				Vertex p;
				if(edges[i].inCore)	p =  mesh->inCorePoints[edges[i].index       ];
				else				p = (*interiorVertices)[edges[i].index-offSet];
				c += p;
			}
			c /= Real( edges.size() );
			int cIdx;
			cIdx = mesh->addOutOfCorePoint_s( c );
#pragma omp critical (add_barycenter_access)
			barycenters->push_back( c );
			for( int i=0 ; i<int(edges.size()) ; i++ )
			{
				std::vector< CoredVertexIndex > vertices( 3 );
				vertices[0].idx = edges[i                 ].index;
				vertices[1].idx = edges[(i+1)%edges.size()].index;
				vertices[2].idx = cIdx;
				vertices[0].inCore = (edges[i                 ].inCore!=0);
				vertices[1].inCore = (edges[(i+1)%edges.size()].inCore!=0);
				vertices[2].inCore = 0;
				mesh->addPolygon_s( vertices );
			}
			return int( edges.size() );
		}
		else
		{
			vertices.resize( edges.size() );
			// Add the points
			for( int i=0 ; i<int(edges.size()) ; i++ )
			{
				Vertex p;
				if( edges[i].inCore ) p =  mesh->inCorePoints[edges[i].index       ];
				else                  p = (*interiorVertices)[edges[i].index-offSet];
				vertices[i] = p.point;
			}
			MAT.GetTriangulation( vertices , triangles );
			for( int i=0 ; i<int(triangles.size()) ; i++ )
			{
				std::vector< CoredVertexIndex > _vertices( 3 );
				for( int j=0 ; j<3 ; j++ )
				{
					_vertices[j].idx    =  edges[ triangles[i].idx[j] ].index;
					_vertices[j].inCore = (edges[ triangles[i].idx[j] ].inCore!=0);
				}
				mesh->addPolygon_s( _vertices );
			}
		}
	}
	else if( edges.size()==3 )
	{
		std::vector< CoredVertexIndex > vertices( 3 );
		for( int i=0 ; i<3 ; i++ )
		{
			vertices[i].idx    =  edges[i].index;
			vertices[i].inCore = (edges[i].inCore!=0);
		}
		mesh->addPolygon_s( vertices );
	}
	return int(edges.size())-2;
}
template< class Real , int Degree >
Real Octree< Real , Degree >::GetSolutionValue( Point3D< Real > p , const BSplineData< Degree >* fData ) const
{
	Real value = Real(0);
	BSplineData< Degree , Real > _fData;
	if( !fData ) _fData.set( tree.maxDepth() , true , _boundaryType ) , fData = &_fData;
	const TreeOctNode* n = tree.nextNode();
	while( n )
	{
		Point3D< Real > c;
		Real w;
		n->centerAndWidth( c , w );
		c -= p , w *= Real(1.5);
		if( fabs(c[0])>w || fabs(c[1])>w || fabs(c[2])>w )
		{
			n = tree.nextBranch( n );
			continue;
		}
		int d , off[3];
		n->depthAndOffset( d , off );
		value += (Real)
			(
			n->nodeData.solution *
			fData->baseFunctions[ BinaryNode< Real >::CenterIndex( d , off[0] ) ]( p[0] ) *
			fData->baseFunctions[ BinaryNode< Real >::CenterIndex( d , off[1] ) ]( p[1] ) *
			fData->baseFunctions[ BinaryNode< Real >::CenterIndex( d , off[2] ) ]( p[2] )
			);
		n = tree.nextNode( n );
	}
	if( _boundaryType==-1 ) value -= Real(0.5);
	return value;
}
template< class Real , int Degree >
Pointer( Real ) Octree< Real , Degree >::GetSolutionGrid( ConstPointer( Real ) solution , int& res , Real isoValue , int depth )
{
	int maxDepth = _boundaryType==0 ? tree.maxDepth()-1 : tree.maxDepth();
	if( depth<=0 || depth>maxDepth ) depth = maxDepth;
	res = 1<<depth;
	typename BSplineData< Degree >::ValueTables< Real > vTables = _fData.getValueTables< Real >( _fData.VALUE_FLAG );
	Pointer( Real ) values = NewPointer< Real >( res * res * res );
	memset( values , 0 , sizeof( Real ) * res  * res * res );

	for( TreeOctNode* n=tree.nextNode() ; n ; n=tree.nextNode( n ) )
	{
		if( n->depth()>(_boundaryType==0?depth+1:depth) ) continue;
		if( n->depth()<_minDepth ) continue;
		int d , idx[3] , start[3] , end[3];
		n->depthAndOffset( d , idx );
		bool skip=false;
		for( int i=0 ; i<3 ; i++ )
		{
			// Get the index of the functions
			idx[i] = BinaryNode< double >::CenterIndex( d , idx[i] );
			// Figure out which samples fall into the range
			vTables.setSampleSpan( idx[i] , start[i] , end[i] );
			// We only care about the odd indices
			if( !(start[i]&1) ) start[i]++;
			if( !(  end[i]&1) )   end[i]--;
			if( _boundaryType==0 )
			{
				// (start[i]-1)>>1 >=   res/2 
				// (  end[i]-1)<<1 <  3*res/2
				start[i] = std::max< int >( start[i] ,   res+1 );
				end  [i] = std::min< int >( end  [i] , 3*res-1 );
			}
		}
		if( skip ) continue;
		Real coefficient = solution[ n->nodeData.nodeIndex ];
		for( int x=start[0] ; x<=end[0] ; x+=2 )
			for( int y=start[1] ; y<=end[1] ; y+=2 )
				for( int z=start[2] ; z<=end[2] ; z+=2 )
				{
					int xx = (x-1)>>1 , yy=(y-1)>>1 , zz = (z-1)>>1;
					if( _boundaryType==0 ) xx -= res/2 , yy -= res/2 , zz -= res/2;
					values[ zz*res*res + yy*res + xx ] +=
						coefficient *
						vTables.valueTable[ idx[0] + x*vTables.functionCount ] *
						vTables.valueTable[ idx[1] + y*vTables.functionCount ] *
						vTables.valueTable[ idx[2] + z*vTables.functionCount ];
				}
	}
	if( _boundaryType==-1 ) for( int i=0 ; i<res*res*res ; i++ ) values[i] -= Real(0.5);
	for( int i=0 ; i<res*res*res ; i++ ) values[i] -= isoValue;

	return values;
}

////////////////
// VertexData //
////////////////
template< class Real >
long long VertexData< Real >::CenterIndex(const TreeOctNode* node,int maxDepth)
{
	int idx[DIMENSION];
	return CenterIndex(node,maxDepth,idx);
}
template< class Real >
long long VertexData< Real >::CenterIndex(const TreeOctNode* node,int maxDepth,int idx[DIMENSION])
{
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++) idx[i]=BinaryNode<Real>::CornerIndex( maxDepth+1 , d+1 , o[i]<<1 , 1 );
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
template< class Real >
long long VertexData< Real >::CenterIndex( int depth , const int offSet[DIMENSION] , int maxDepth , int idx[DIMENSION] )
{
	for(int i=0;i<DIMENSION;i++) idx[i]=BinaryNode<Real>::CornerIndex( maxDepth+1 , depth+1 , offSet[i]<<1 , 1 );
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
template< class Real >
long long VertexData< Real >::CornerIndex(const TreeOctNode* node,int cIndex,int maxDepth)
{
	int idx[DIMENSION];
	return CornerIndex(node,cIndex,maxDepth,idx);
}
template< class Real >
long long VertexData< Real >::CornerIndex( const TreeOctNode* node , int cIndex , int maxDepth , int idx[DIMENSION] )
{
	int x[DIMENSION];
	Cube::FactorCornerIndex( cIndex , x[0] , x[1] , x[2] );
	int d , o[3];
	node->depthAndOffset( d , o );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode<Real>::CornerIndex( maxDepth+1 , d , o[i] , x[i] );
	return CornerIndexKey( idx );
}
template< class Real >
long long VertexData< Real >::CornerIndex( int depth , const int offSet[DIMENSION] , int cIndex , int maxDepth , int idx[DIMENSION] )
{
	int x[DIMENSION];
	Cube::FactorCornerIndex( cIndex , x[0] , x[1] , x[2] );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode<Real>::CornerIndex( maxDepth+1 , depth , offSet[i] , x[i] );
	return CornerIndexKey( idx );
}
template< class Real >
long long VertexData< Real >::CornerIndexKey( const int idx[DIMENSION] )
{
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
template< class Real >
long long VertexData< Real >::FaceIndex(const TreeOctNode* node,int fIndex,int maxDepth){
	int idx[DIMENSION];
	return FaceIndex(node,fIndex,maxDepth,idx);
}
template< class Real >
long long VertexData< Real >::FaceIndex(const TreeOctNode* node,int fIndex,int maxDepth,int idx[DIMENSION])
{
	int dir,offset;
	Cube::FactorFaceIndex(fIndex,dir,offset);
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);}
	idx[dir]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,o[dir],offset);
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
template< class Real >
long long VertexData< Real >::EdgeIndex(const TreeOctNode* node,int eIndex,int maxDepth)
{
	int idx[DIMENSION];
	return EdgeIndex(node,eIndex,maxDepth,idx);
}
template< class Real >
long long VertexData< Real >::EdgeIndex(const TreeOctNode* node,int eIndex,int maxDepth,int idx[DIMENSION])
{
	int o,i1,i2;
	int d,off[3];
	node->depthAndOffset(d,off);
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,d+1,off[i]<<1,1);}
	Cube::FactorEdgeIndex(eIndex,o,i1,i2);
	switch(o){
		case 0:
			idx[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[1],i1);
			idx[2]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[2],i2);
			break;
		case 1:
			idx[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[0],i1);
			idx[2]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[2],i2);
			break;
		case 2:
			idx[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[0],i1);
			idx[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,off[1],i2);
			break;
	};
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
