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

#define LOW_MEMORY 1

/////////////////////
// SortedTreeNodes //
/////////////////////
template< class Real >
SortedTreeNodes< Real >::SortedTreeNodes( void )
{
	nodeCount = NULL;
	treeNodes = NullPointer< TreeOctNode* >();
	maxDepth = 0;
	sliceOffsets = NullPointer< Pointer( int ) >();
}
template< class Real >
SortedTreeNodes< Real >::~SortedTreeNodes( void )
{
	if( nodeCount ) delete[] nodeCount , nodeCount = NULL;
	if( treeNodes ) DeletePointer( treeNodes );
	if( sliceOffsets )
	{
		for( int d=0 ; d<maxDepth ; d++ ) FreePointer( sliceOffsets[d] );
		FreePointer( sliceOffsets );
	}
}
template< class Real >
void SortedTreeNodes< Real >::set( TreeOctNode& root , std::vector< int >* map )
{
	if( nodeCount ) delete[] nodeCount;
	if( treeNodes ) DeletePointer( treeNodes );
	if( sliceOffsets )
	{
		for( int d=0 ; d<maxDepth ; d++ ) FreePointer( sliceOffsets[d] );
		FreePointer( sliceOffsets );
	}
	maxDepth = root.maxDepth()+1;
	nodeCount = new int[ maxDepth+1 ];
	treeNodes = NewPointer< TreeOctNode* >( root.nodes() );

	int startDepth = 0;
	nodeCount[0] = 0 , nodeCount[1] = 1;
	treeNodes[0] = &root;
	for( int d=startDepth+1 ; d<maxDepth ; d++ )
	{
		nodeCount[d+1] = nodeCount[d];
		for( int i=nodeCount[d-1] ; i<nodeCount[d] ; i++ )
		{
			TreeOctNode* temp = treeNodes[i];
			if( temp->children ) for( int c=0 ; c<8 ; c++ ) treeNodes[ nodeCount[d+1]++ ] = temp->children + c;
		}
	}
	_sortByZCoordinate();
	if( map )
	{
		map->resize( nodeCount[maxDepth] );
		for( int i=0 ; i<nodeCount[maxDepth] ; i++ ) (*map)[i] = treeNodes[i]->nodeData.nodeIndex;
	}
	for( int i=0 ; i<nodeCount[maxDepth] ; i++ ) treeNodes[i]->nodeData.nodeIndex = i;
}
#if BROODED_SLICES
template< class Real >
int SortedTreeNodes< Real >::Slices( int depth ){ return depth<=0 ? 1 :1<<(depth-1); }
template< class Real >
std::pair< int , int > SortedTreeNodes< Real >::sliceSpan( int depth , int off , int d ) const
{
	int dd;
	if( !d || depth ) dd = d-depth;
	else dd = d-depth-1;
	return std::pair< int , int >( nodeCount[d] + sliceOffsets[d][off<<dd] , nodeCount[d] + sliceOffsets[d][(off+1)<<dd] );
}
#else // !BROODED_SLICES
template< class Real >
int SortedTreeNodes< Real >::Slices( int depth ){ return 1<<depth; }
template< class Real >
std::pair< int , int > SortedTreeNodes< Real >::sliceSpan( int depth , int off , int d ) const
{
	int dd = d-depth;
	return std::pair< int , int >( nodeCount[d] + sliceOffsets[d][off<<dd] , nodeCount[d] + sliceOffsets[d][(off+1)<<dd] );
}
#endif // BROODED_SLICES
template< class Real >
void SortedTreeNodes< Real >::_sortByZCoordinate( void )
{
	sliceOffsets = AllocPointer< Pointer( int ) >( maxDepth );
	for( int d=0 ; d<maxDepth ;  d++ )
	{
		int slices = Slices( d );
		sliceOffsets[d] = AllocPointer< int >( slices+1 );
		memset( sliceOffsets[d] , 0 , sizeof(int)*(slices+1) );
		for( int i=nodeCount[d] ; i<nodeCount[d+1] ; i++ )
		{
			int _d , _off[3];
			treeNodes[i]->depthAndOffset( _d , _off );
#if BROODED_SLICES
			sliceOffsets[d][ _off[2]>>1 ]++;
#else // !BROODED_SLICES
			sliceOffsets[d][ _off[2] ]++;
#endif // BROODED_SLICES
		}
		for( int i=1 ; i<slices ; i++ ) sliceOffsets[d][i] += sliceOffsets[d][i-1];
		for( int i=slices ; i>=1 ; i-- ) sliceOffsets[d][i] = sliceOffsets[d][i-1];
		sliceOffsets[d][0] = 0;
	}
	for( TreeOctNode* node=treeNodes[0]->nextNode() ; node ; node=treeNodes[0]->nextNode( node ) )
	{
		int d , off[3];
		node->depthAndOffset( d , off );
#if BROODED_SLICES
		treeNodes[ nodeCount[d] + sliceOffsets[d][ off[2]>>1 ] ] = node;
		sliceOffsets[d][ off[2]>>1 ]++;
#else // !BROODED_SLICES
		treeNodes[ nodeCount[d] + sliceOffsets[d][ off[2] ] ] = node;
		sliceOffsets[d][ off[2] ]++;
#endif // BROODED_SLICES
	}
	for( int d=0 ; d<maxDepth ; d++ )
	{
		for( int i=Slices(d) ; i>=1 ; i-- ) sliceOffsets[d][i] = sliceOffsets[d][i-1];
		sliceOffsets[d][0] = 0;
	}
}
template< class Real >
typename SortedTreeNodes< Real >::CubeCornerIndices& SortedTreeNodes< Real >::SlabCornerTableData::operator[] ( const TreeOctNode* node ) { return cTable[ node->nodeData.nodeIndex + offsets[node->depth()] ]; }
template< class Real >
const typename SortedTreeNodes< Real >::CubeCornerIndices& SortedTreeNodes< Real >::SlabCornerTableData::operator[] ( const TreeOctNode* node ) const { return cTable[ node->nodeData.nodeIndex + offsets[node->depth()] ]; }
template< class Real >
typename SortedTreeNodes< Real >::CubeCornerIndices& SortedTreeNodes< Real >::SlabCornerTableData::cornerIndices( const TreeOctNode* node ) { return cTable[ node->nodeData.nodeIndex + offsets[node->depth()] ]; }
template< class Real >
const typename SortedTreeNodes< Real >::CubeCornerIndices& SortedTreeNodes< Real >::SlabCornerTableData::cornerIndices( const TreeOctNode* node ) const { return cTable[ node->nodeData.nodeIndex + offsets[node->depth()] ]; }
template< class Real >
void SortedTreeNodes< Real >::setCornerTable( SlabCornerTableData& cData , int minDepth , int off , int maxDepth , int threads ) const
{
	if( threads<=0 ) threads = 1;
	// The vector of per-depth node spans
	std::vector< std::pair< int , int > > spans( this->maxDepth , std::pair< int , int >( -1 , -1 ) );
	cData.offsets.resize( this->maxDepth , -1 );
	int nodeCount = 0;
	for( int d=minDepth ; d<=maxDepth ; d++ )
	{
		spans[d] = sliceSpan( minDepth , off , d );
		cData.offsets[d] = nodeCount - spans[d].first;
		nodeCount += spans[d].second - spans[d].first;
	}
	cData.cTable.resize( nodeCount );
	int count = 0;
	typename TreeOctNode::ConstNeighborKey3 neighborKey;
	neighborKey.set( maxDepth );
#if LOW_MEMORY
	for( int d=minDepth ; d<=maxDepth ; d++ )
	{
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
		for( int i=spans[d].first ; i<spans[d].second ; i++ )
		{
			TreeOctNode* node = treeNodes[i];
			if( d<maxDepth && node->children ) continue;
			const typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors( node );
			int _d , _off[3];
			node->depthAndOffset( _d , _off );
			int dd = (d-minDepth);
#if BROODED_SLICES
			dd++;
#endif // BROODED_SLICES
			for( int c=0 ; c<Cube::CORNERS ; c++ )	// Iterate over the cell's corners
			{
				bool cornerOwner = true;
				int x , y , z;
				int ac = Cube::AntipodalCornerIndex( c ); // The index of the node relative to the corner
				Cube::FactorCornerIndex( c , x , y , z );
				for( int cc=0 ; cc<Cube::CORNERS ; cc++ ) // Iterate over the corner's cells
				{
					int xx , yy , zz;
					Cube::FactorCornerIndex( cc , xx , yy , zz );
					xx += x , yy += y , zz += z;
					if( neighbors.neighbors[xx][yy][zz] && ( cc<ac || ( d<maxDepth && neighbors.neighbors[xx][yy][zz]->children ) ) )
						if( ( _off[2]>>dd )==( (_off[2]+zz-1)>>dd ) ){ cornerOwner = false ; break; }
				}
				if( cornerOwner )
				{
					int myCount = ( i + cData.offsets[d] ) * Cube::CORNERS + c;

					const TreeOctNode* n = node;
					int d = n->depth();
					do
					{
						const typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d];
						// Set all the corner indices at the current depth
						for( int cc=0 ; cc<Cube::CORNERS ; cc++ )
						{
							int xx , yy , zz;
							Cube::FactorCornerIndex( cc , xx , yy , zz );
							xx += x , yy += y , zz += z;
							if( neighborKey.neighbors[d].neighbors[xx][yy][zz] )
								if( ( _off[2]>>dd )==( (_off[2]+zz-1)>>dd ) ) cData[ neighbors.neighbors[xx][yy][zz] ][ Cube::AntipodalCornerIndex(cc) ] = myCount;
						}
						// If we are not at the root and the parent also has the corner
						if( d==minDepth || n!=(n->parent->children+c) ) break;
						n = n->parent;
						d--;
					}
					while( 1 );
				}
			}
		}
	}
	std::vector< int > map( nodeCount );
	for( int c=0 ; c<Cube::CORNERS ; c++ )
	{
		memset( &map[0] , 0 , sizeof(int)*nodeCount );
		int start = c*nodeCount , end = (c+1)*nodeCount;
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<nodeCount ; i++ ) for( int cc=0 ; cc<Cube::CORNERS ; cc++ )
		{
			int idx = cData.cTable[i][cc]; 
			if( idx>=start && idx<end ) map[idx%nodeCount] = 1;
		}
		for( int i=0 ; i<nodeCount ; i++ ) if( map[i] ) map[i] = count++;
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<nodeCount ; i++ ) for( int cc=0 ; cc<Cube::CORNERS ; cc++ )
		{
			int idx = cData.cTable[i][cc]; 
			if( idx>=start && idx<end ) cData.cTable[i][cc] = map[idx%nodeCount];
		}
	}
#else // !LOW_MEMORY
	std::vector< int > cIndices( nodeCount*Cube::CORNERS , 0 );
	for( int d=minDepth ; d<=maxDepth ; d++ )
	{
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
		for( int i=spans[d].first ; i<spans[d].second ; i++ )
		{
			TreeOctNode* node = treeNodes[i];
			if( d<maxDepth && node->children ) continue;
			const typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors( node );
			int _d , _off[3];
			node->depthAndOffset( _d , _off );
			int dd = (d-minDepth);
#if BROODED_SLICES
			dd++;
#endif // BROODED_SLICES
			for( int c=0 ; c<Cube::CORNERS ; c++ )	// Iterate over the cell's corners
			{
				bool cornerOwner = true;
				int x , y , z;
				int ac = Cube::AntipodalCornerIndex( c ); // The index of the node relative to the corner
				Cube::FactorCornerIndex( c , x , y , z );
				for( int cc=0 ; cc<Cube::CORNERS ; cc++ ) // Iterate over the corner's cells
				{
					int xx , yy , zz;
					Cube::FactorCornerIndex( cc , xx , yy , zz );
					xx += x , yy += y , zz += z;
					if( neighbors.neighbors[xx][yy][zz] && ( cc<ac || ( d<maxDepth && neighbors.neighbors[xx][yy][zz]->children ) ) )
						if( ( _off[2]>>dd )==( (_off[2]+zz-1)>>dd ) ){ cornerOwner = false ; break; }
				}
				if( cornerOwner )
				{
					int myCount = ( treeNodes[i]->nodeData.nodeIndex+cData.offsets[d] )*Cube::CORNERS + c;
					cIndices[myCount] = 1;

					const TreeOctNode* n = node;
					int d = n->depth();
					do
					{
						const typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d];
						// Set all the corner indices at the current depth
						for( int cc=0 ; cc<Cube::CORNERS ; cc++ )
						{
							int xx , yy , zz;
							Cube::FactorCornerIndex( cc , xx , yy , zz );
							xx += x , yy += y , zz += z;
							if( neighborKey.neighbors[d].neighbors[xx][yy][zz] )
								if( ( _off[2]>>dd )==( (_off[2]+zz-1)>>dd ) ) cData[ neighbors.neighbors[xx][yy][zz] ][ Cube::AntipodalCornerIndex(cc) ] = myCount;
						}
						// If we are not at the root and the parent also has the corner
						if( d==minDepth || n!=(n->parent->children+c) ) break;
						n = n->parent;
						d--;
					}
					while( 1 );
				}
			}
		}
	}
	for( int i=0 ; i<cIndices.size() ; i++ ) if( cIndices[i] ) cIndices[i] = count++;
	for( int d=minDepth ; d<=maxDepth ; d++ )
#pragma omp parallel for num_threads( threads )
		for( int i=spans[d].first ; i<spans[d].second ; i++ ) for( int j=0 ; j<Cube::CORNERS ; j++ ) cData[ treeNodes[i] ][j] = cIndices[ cData[ treeNodes[i] ][j] ];
#endif // LOW_MEMORY
	cData.cCount = count;
}
template< class Real >
int SortedTreeNodes< Real >::getMaxSlabCornerCount( int depth , int maxDepth , int threads ) const
{
	if( threads<=0 ) threads = 1;
	int slabs = 1<< depth;
	std::vector< int > cornerCount( slabs , 0 );
	typename TreeOctNode::ConstNeighborKey3 neighborKey;
	neighborKey.set( maxDepth );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
	for( int i=nodeCount[depth] ; i<nodeCount[maxDepth+1] ; i++ )
	{
		TreeOctNode* node = treeNodes[i];
		int d , off[3];
		node->depthAndOffset( d , off );
		if( d<maxDepth && node->children ) continue;

		const typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors( node );
		int dd = (d-depth);
#if BROODED_SLICES
		dd++;
#endif // BROODED_SLICES
		for( int c=0 ; c<Cube::CORNERS ; c++ )	// Iterate over the cell's corners
		{
			bool cornerOwner = true;
			int x , y , z;
			int ac = Cube::AntipodalCornerIndex( c ); // The index of the node relative to the corner
			Cube::FactorCornerIndex( c , x , y , z );
			for( int cc=0 ; cc<Cube::CORNERS ; cc++ ) // Iterate over the corner's cells
			{
				int xx , yy , zz;
				Cube::FactorCornerIndex( cc , xx , yy , zz );
				xx += x , yy += y , zz += z;
				if( neighbors.neighbors[xx][yy][zz] && ( cc<ac || ( d<maxDepth && neighbors.neighbors[xx][yy][zz]->children ) ) )
					if( ( off[2]>>dd )==( (off[2]+zz-1)>>dd ) ){ cornerOwner = false ; break; }
			}
			if( cornerOwner )
#pragma omp atomic
				cornerCount[ ( off[2]>>(d-depth) ) ]++;
		}
	}
	int maxCount = 0;
#if BROODED_SLICES
	if( depth==0 ) maxCount = std::max< int >( maxCount , cornerCount[0] );
	for( int i=0 ; i<slabs/2 ; i++ ) maxCount = std::max< int >( maxCount , cornerCount[2*i] + cornerCount[2*i+1] );
#else // !BROODED_SLICES
	for( int i=0 ; i<slabs ; i++ ) maxCount = std::max< int >( maxCount , cornerCount[i] );
#endif // BROODED_SLICES
	return maxCount;
}
template< class Real >
typename SortedTreeNodes< Real >::EdgeIndices& SortedTreeNodes< Real >::SlabEdgeTableData::operator[] ( const TreeOctNode* node ) { return eTable[ node->nodeData.nodeIndex + offsets[node->depth()] ]; }
template< class Real >
const typename SortedTreeNodes< Real >::EdgeIndices& SortedTreeNodes< Real >::SlabEdgeTableData::operator[] ( const TreeOctNode* node ) const { return eTable[ node->nodeData.nodeIndex + offsets[node->depth()] ]; }
template< class Real >
typename SortedTreeNodes< Real >::EdgeIndices& SortedTreeNodes< Real >::SlabEdgeTableData::edgeIndices( const TreeOctNode* node ) { return eTable[ node->nodeData.nodeIndex + offsets[node->depth()] ]; }
template< class Real >
const typename SortedTreeNodes< Real >::EdgeIndices& SortedTreeNodes< Real >::SlabEdgeTableData::edgeIndices( const TreeOctNode* node ) const { return eTable[ node->nodeData.nodeIndex + offsets[node->depth()] ]; }
template< class Real >
void SortedTreeNodes< Real >::setEdgeTable( SlabEdgeTableData& eData , int minDepth , int off , int maxDepth , int threads )
{
	if( threads<=0 ) threads = 1;
	std::vector< std::pair< int , int > > spans( this->maxDepth , std::pair< int , int >( -1 , -1 ) );

	eData.offsets.resize( this->maxDepth , -1 );
	int nodeCount = 0;
	for( int d=minDepth ; d<=maxDepth ; d++ )
	{
		spans[d] = sliceSpan( minDepth , off , d );
		eData.offsets[d] = nodeCount - spans[d].first;
		nodeCount += spans[d].second - spans[d].first;
	}
	eData.eTable.resize( nodeCount );
	int count = 0;
	typename TreeOctNode::ConstNeighborKey3 neighborKey;
	neighborKey.set( maxDepth );
#if LOW_MEMORY
	for( int d=minDepth ; d<=maxDepth ; d++ )
	{
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
		for( int i=spans[d].first ; i<spans[d].second ; i++ )
		{
			TreeOctNode* node = treeNodes[i];
			const typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors( node );
			int _d , _off[3];
			node->depthAndOffset( _d , _off );
			int dd = (d-minDepth);
#if BROODED_SLICES
			dd++;
#endif // BROODED_SLICES

			for( int e=0 ; e<Cube::EDGES ; e++ )
			{
				bool edgeOwner = true;
				int o , _i , _j;
				Cube::FactorEdgeIndex( e , o , _i , _j );
				int ac = Square::AntipodalCornerIndex( Square::CornerIndex( _i , _j ) );
				for( int cc=0 ; cc<Square::CORNERS ; cc++ )
				{
					int ii , jj , x , y , z;
					Square::FactorCornerIndex( cc , ii , jj );
					ii += _i , jj += _j;
					switch( o )
					{
					case 0: y = ii , z = jj , x = 1 ; break;
					case 1: x = ii , z = jj , y = 1 ; break;
					case 2: x = ii , y = jj , z = 1 ; break;
					}
					if( neighbors.neighbors[x][y][z] && cc<ac )
						if( ( _off[2]>>dd )==( (_off[2]+z-1)>>dd ) ){ edgeOwner = false ; break; }
				}
				if( edgeOwner )
				{
					int myCount = ( i + eData.offsets[d] ) * Cube::EDGES + e;
					// Set all edge indices
					for( int cc=0 ; cc<Square::CORNERS ; cc++ )
					{
						int ii , jj , aii , ajj , x , y , z;
						Square::FactorCornerIndex( cc , ii , jj );
						Square::FactorCornerIndex( Square::AntipodalCornerIndex( cc ) , aii , ajj );
						ii += _i , jj += _j;
						switch( o )
						{
						case 0: y = ii , z = jj , x = 1 ; break;
						case 1: x = ii , z = jj , y = 1 ; break;
						case 2: x = ii , y = jj , z = 1 ; break;
						}
						if( neighbors.neighbors[x][y][z] )
							if( ( _off[2]>>dd )==( (_off[2]+z-1)>>dd ) ) eData[ neighbors.neighbors[x][y][z] ][ Cube::EdgeIndex( o , aii , ajj ) ] = myCount;
					}
				}
			}
		}
	}
	std::vector< int > map( nodeCount );
	for( int e=0 ; e<Cube::EDGES ; e++ )
	{
		memset( &map[0] , 0 , sizeof(int)*nodeCount );
		int start = e*nodeCount , end = (e+1)*nodeCount;
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<nodeCount ; i++ ) for( int ee=0 ; ee<Cube::EDGES ; ee++ )
		{
			int idx = eData.eTable[i][ee]; 
			if( idx>=start && idx<end ) map[idx%nodeCount] = 1;
		}
		for( int i=0 ; i<nodeCount ; i++ ) if( map[i] ) map[i] = count++;
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<nodeCount ; i++ ) for( int ee=0 ; ee<Cube::EDGES ; ee++ )
		{
			int idx = eData.eTable[i][ee]; 
			if( idx>=start && idx<end ) eData.eTable[i][ee] = map[idx%nodeCount];
		}
	}
#else // !LOW_MEMORY
	std::vector< int > eIndices( nodeCount*Cube::EDGES , 0 );
	for( int d=minDepth ; d<=maxDepth ; d++ )
	{
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
		for( int i=spans[d].first ; i<spans[d].second ; i++ )
		{
			TreeOctNode* node = treeNodes[i];
			const typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors( node );
			int _d , _off[3];
			node->depthAndOffset( _d , _off );
			int dd = (d-minDepth);
#if BROODED_SLICES
			dd++;
#endif // BROODED_SLICES

			for( int e=0 ; e<Cube::EDGES ; e++ )
			{
				bool edgeOwner = true;
				int o , _i , _j;
				Cube::FactorEdgeIndex( e , o , _i , _j );
				int ac = Square::AntipodalCornerIndex( Square::CornerIndex( _i , _j ) );
				for( int cc=0 ; cc<Square::CORNERS ; cc++ )
				{
					int ii , jj , x , y , z;
					Square::FactorCornerIndex( cc , ii , jj );
					ii += _i , jj += _j;
					switch( o )
					{
					case 0: y = ii , z = jj , x = 1 ; break;
					case 1: x = ii , z = jj , y = 1 ; break;
					case 2: x = ii , y = jj , z = 1 ; break;
					}
					if( neighbors.neighbors[x][y][z] && cc<ac )
						if( ( _off[2]>>dd )==( (_off[2]+z-1)>>dd ) ){ edgeOwner = false ; break; }
				}
				if( edgeOwner )
				{
					int myCount = (treeNodes[i]->nodeData.nodeIndex + eData.offsets[d])*Cube::EDGES + e;
					eIndices[myCount] = 1;
					// Set all edge indices
					for( int cc=0 ; cc<Square::CORNERS ; cc++ )
					{
						int ii , jj , aii , ajj , x , y , z;
						Square::FactorCornerIndex( cc , ii , jj );
						Square::FactorCornerIndex( Square::AntipodalCornerIndex( cc ) , aii , ajj );
						ii += _i , jj += _j;
						switch( o )
						{
						case 0: y = ii , z = jj , x = 1 ; break;
						case 1: x = ii , z = jj , y = 1 ; break;
						case 2: x = ii , y = jj , z = 1 ; break;
						}
						if( neighbors.neighbors[x][y][z] )
							if( ( _off[2]>>dd )==( (_off[2]+z-1)>>dd ) ) eData[ neighbors.neighbors[x][y][z] ][ Cube::EdgeIndex( o , aii , ajj ) ] = myCount;
					}
				}
			}
		}
	}
	for( int i=0 ; i<eIndices.size() ; i++ ) if( eIndices[i] ) eIndices[i] = count++;
	for( int d=minDepth ; d<=maxDepth ; d++ )
#pragma omp parallel for num_threads( threads )
		for( int i=spans[d].first ; i<spans[d].second ; i++ ) for( int j=0 ; j<Cube::EDGES ; j++ ) eData[ treeNodes[i] ][j] = eIndices[ eData[ treeNodes[i] ][j] ];
#endif // LOW_MEMORY
	eData.eCount = count;
}
template< class Real >
int SortedTreeNodes< Real >::getMaxSlabEdgeCount( const TreeOctNode* rootNode , int depth , int threads ) const
{
	if( threads<=0 ) threads = 1;
	int slices = 1<<depth;
	std::vector< int > edgeCount( slices , 0 );
	typename TreeOctNode::ConstNeighborKey3 neighborKey;
	neighborKey.set( maxDepth-1 );
#pragma omp parallel for num_threads( threads ) firstprivate( neighborKey )
	for( int i=nodeCount[depth] ; i<nodeCount[maxDepth] ; i++ )
	{
		TreeOctNode* node = treeNodes[i];
		const typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors( node );
		int d , off[3];
		node->depthAndOffset( d , off );
		int dd = (d-depth);
#if BROODED_SLICES
		dd++;
#endif // BROODED_SLICES

		for( int e=0 ; e<Cube::EDGES ; e++ )
		{
			bool edgeOwner = true;
			int o , i , j;
			Cube::FactorEdgeIndex( e , o , i , j );
			int ac = Square::AntipodalCornerIndex( Square::CornerIndex( i , j ) );
			for( int cc=0 ; cc<Square::CORNERS ; cc++ )
			{
				int ii , jj , x , y , z;
				Square::FactorCornerIndex( cc , ii , jj );
				ii += i , jj += j;
				switch( o )
				{
				case 0: y = ii , z = jj , x = 1 ; break;
				case 1: x = ii , z = jj , y = 1 ; break;
				case 2: x = ii , y = jj , z = 1 ; break;
				}
				if( neighbors.neighbors[x][y][z] && cc<ac )
					if( ( off[2]>>dd )==( (off[2]+z-1)>>dd ) ){ edgeOwner = false ; break; } 
			}
			if( edgeOwner )
#pragma omp atomic
				edgeCount[ off[2]>>(d-depth) ]++;
		}
	}
	int maxCount = 0;
#if BROODED_SLICES
	if( depth==0 ) maxCount = std::max< int >( maxCount , edgeCount[0] );
	else for( int i=0 ; i<slices/2 ; i++ ) maxCount = std::max< int >( maxCount , edgeCount[2*i]+edgeCount[2*i+1] );
#else // !BRODED_SLICES
	for( int i=0 ; i<slices ; i++ ) maxCount = std::max< int >( maxCount , edgeCount[i] );
#endif // BROODED_SLICES
	return maxCount;
}
