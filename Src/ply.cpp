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

#include "ply.h"

//
// PLY data structures
//
char *elem_names[] = { "vertex", "face" };

typedef struct PlyVertex
{
	float x, y, z;
} PlyVertex;
typedef struct PlyVertex2
{
	float x1 , x2 , y1 , y2 , z1 , z2 , s;
} PlyVertex2;

typedef struct PlyOrientedVertex
{
	float x, y, z , nx, ny, nz;
} PlyOrientedVertex;

typedef struct PlyFace
{
	unsigned char nr_vertices;
	int *vertices;
	int segment;
} PlyFace;

static PlyProperty vert_props[] = {
	{"x", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,x), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,y), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,z), 0, 0, 0, 0}
};
static PlyProperty vert_props_2[] =
{
	{"x1", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex2,x1), 0, 0, 0, 0},
	{"x2", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex2,x2), 0, 0, 0, 0},
	{"y1", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex2,y1), 0, 0, 0, 0},
	{"y2", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex2,y2), 0, 0, 0, 0},
	{"z1", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex2,z1), 0, 0, 0, 0},
	{"z2", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex2,z2), 0, 0, 0, 0},
	{"s" , PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex2,s ), 0, 0, 0, 0},
};
static PlyProperty oriented_vert_props[] = {
	{"x",  PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,x ), 0, 0, 0, 0},
	{"y",  PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,y ), 0, 0, 0, 0},
	{"z",  PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,z ), 0, 0, 0, 0},
	{"nx", PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,nx), 0, 0, 0, 0},
	{"ny", PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,ny), 0, 0, 0, 0},
	{"nz", PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,nz), 0, 0, 0, 0}
};

// List of property information for a vertex
static PlyProperty face_props[] = {
	{ "vertex_indices" , PLY_INT, PLY_INT, offsetof(PlyFace,vertices),
		1, PLY_UCHAR, PLY_UCHAR, offsetof(PlyFace,nr_vertices)},
};

int PlyWritePolygons( char* fileName , CoredMeshData* mesh , int file_type , const Point3D<float>& translate , float scale , char** comments , int commentNum , XForm4x4< float > xForm )
{
	XForm3x3< float > xFormN;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) xFormN(i,j) = xForm(i,j);
	xFormN = xFormN.transpose().inverse();
	int i;
	int nr_vertices=int(mesh->outOfCorePointCount()+mesh->inCorePoints.size());
	int nr_faces=mesh->polygonCount();
	float version;
	PlyFile *ply = ply_open_for_writing(fileName, 2, elem_names, file_type, &version);
	if( !ply ) return 0;

	mesh->resetIterator();
	
	//
	// describe vertex and face properties
	//
	ply_element_count( ply , "vertex" , nr_vertices );
	ply_describe_property( ply , "vertex" , &vert_props[0] );
	ply_describe_property( ply , "vertex" , &vert_props[1] );
	ply_describe_property( ply , "vertex" , &vert_props[2] );
	
	ply_element_count( ply , "face" , nr_faces );
	ply_describe_property( ply , "face" , &face_props[0] );
	
	// Write in the comments
	for( i=0 ; i<commentNum ; i++ ) ply_put_comment( ply , comments[i] );

	ply_header_complete( ply );
	
	// write vertices
	ply_put_element_setup( ply , "vertex" );
	Point3D< float > p;
	for( i=0 ; i<int( mesh->inCorePoints.size() ) ; i++ )
	{
		PlyVertex ply_vertex;
		p = mesh->inCorePoints[i] * scale + translate;
		p = xForm * p;
		ply_vertex.x = p[0];
		ply_vertex.y = p[1];
		ply_vertex.z = p[2];
		ply_put_element(ply, (void *) &ply_vertex);
	}
	for( i=0; i<mesh->outOfCorePointCount() ; i++ )
	{
		PlyVertex ply_vertex;
		mesh->nextOutOfCorePoint(p);
		p = p * scale +translate;
		p = xForm * p;
		ply_vertex.x = p[0];
		ply_vertex.y = p[1];
		ply_vertex.z = p[2];
		ply_put_element(ply, (void *) &ply_vertex);		
	}  // for, write vertices
	
	// write faces
	std::vector< CoredVertexIndex > polygon;
	ply_put_element_setup( ply , "face" );
	for( i=0 ; i<nr_faces ; i++ )
	{
		//
		// create and fill a struct that the ply code can handle
		//
		PlyFace ply_face;
		mesh->nextPolygon( polygon );
		ply_face.nr_vertices = int( polygon.size() );
		ply_face.vertices = new int[ polygon.size() ];
		for( int i=0 ; i<int(polygon.size()) ; i++ )
			if( polygon[i].inCore ) ply_face.vertices[i] = polygon[i].idx;
			else                    ply_face.vertices[i] = polygon[i].idx + int( mesh->inCorePoints.size() );
		ply_put_element( ply, (void *) &ply_face );
		delete[] ply_face.vertices;
	}  // for, write faces
	
	ply_close( ply );
	return 1;
}
int PlyWritePolygons( char* fileName , CoredMeshData* mesh , int file_type , char** comments , int commentNum , XForm4x4< float > xForm )
{
	XForm3x3< float > xFormN;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) xFormN(i,j) = xForm(i,j);
	xFormN = xFormN.transpose().inverse();
	int i;
	int nr_vertices=int(mesh->outOfCorePointCount()+mesh->inCorePoints.size());
	int nr_faces=mesh->polygonCount();
	float version;
	PlyFile *ply = ply_open_for_writing(fileName, 2, elem_names, file_type, &version);
	if( !ply ) return 0;

	mesh->resetIterator();
	
	//
	// describe vertex and face properties
	//
	ply_element_count( ply , "vertex" , nr_vertices );
	ply_describe_property( ply , "vertex" , &vert_props[0] );
	ply_describe_property( ply , "vertex" , &vert_props[1] );
	ply_describe_property( ply , "vertex" , &vert_props[2] );
	
	ply_element_count( ply , "face" , nr_faces );
	ply_describe_property( ply , "face" , &face_props[0] );
	
	// Write in the comments
	for( i=0 ; i<commentNum ; i++ ) ply_put_comment( ply , comments[i] );

	ply_header_complete( ply );
	
	// write vertices
	ply_put_element_setup( ply , "vertex" );
	Point3D< float > p;
	for( i=0 ; i<int( mesh->inCorePoints.size() ) ; i++ )
	{
		PlyVertex ply_vertex;
		p = xForm * mesh->inCorePoints[i];
		ply_vertex.x = p[0];
		ply_vertex.y = p[1];
		ply_vertex.z = p[2];
		ply_put_element(ply, (void *) &ply_vertex);
	}
	for( i=0; i<mesh->outOfCorePointCount() ; i++ )
	{
		PlyVertex ply_vertex;
		mesh->nextOutOfCorePoint(p);
		p = xForm * p;
		ply_vertex.x = p[0];
		ply_vertex.y = p[1];
		ply_vertex.z = p[2];
		ply_put_element(ply, (void *) &ply_vertex);		
	}  // for, write vertices
	
	// write faces
	std::vector< CoredVertexIndex > polygon;
	ply_put_element_setup( ply , "face" );
	for( i=0 ; i<nr_faces ; i++ )
	{
		//
		// create and fill a struct that the ply code can handle
		//
		PlyFace ply_face;
		mesh->nextPolygon( polygon );
		ply_face.nr_vertices = int( polygon.size() );
		ply_face.vertices = new int[ polygon.size() ];
		for( int i=0 ; i<int(polygon.size()) ; i++ )
			if( polygon[i].inCore ) ply_face.vertices[i] = polygon[i].idx;
			else                    ply_face.vertices[i] = polygon[i].idx + int( mesh->inCorePoints.size() );
		ply_put_element( ply, (void *) &ply_face );
		delete[] ply_face.vertices;
	}  // for, write faces
	
	ply_close( ply );
	return 1;
}


int PlyWritePolygons( char* fileName , CoredMeshData2* mesh , int file_type , const Point3D<float>& translate , float scale , char** comments , int commentNum )
{
	int i;
	int nr_vertices=int(mesh->outOfCorePointCount()+mesh->inCorePoints.size());
	int nr_faces=mesh->polygonCount();
	float version;
	PlyFile *ply = ply_open_for_writing(fileName, 2, elem_names, file_type, &version);
	if( !ply ) return 0;

	mesh->resetIterator();
	
	//
	// describe vertex and face properties
	//
	ply_element_count( ply , "vertex" , nr_vertices );
	ply_describe_property( ply , "vertex" , &vert_props_2[0] );
	ply_describe_property( ply , "vertex" , &vert_props_2[1] );
	ply_describe_property( ply , "vertex" , &vert_props_2[2] );
	ply_describe_property( ply , "vertex" , &vert_props_2[3] );
	ply_describe_property( ply , "vertex" , &vert_props_2[4] );
	ply_describe_property( ply , "vertex" , &vert_props_2[5] );
	ply_describe_property( ply , "vertex" , &vert_props_2[6] );
	
	ply_element_count( ply , "face" , nr_faces );
	ply_describe_property( ply , "face" , &face_props[0] );
	
	// Write in the comments
	for( i=0 ; i<commentNum ; i++ ) ply_put_comment( ply , comments[i] );

	ply_header_complete( ply );
	
	// write vertices
	ply_put_element_setup( ply , "vertex" );
	for( i=0 ; i<int( mesh->inCorePoints.size() ) ; i++ )
	{
		CoredMeshData2::Vertex p = mesh->inCorePoints[i];
		PlyVertex2 ply_vertex;
		ply_vertex.x1 = p.start[0]*scale+translate[0];
		ply_vertex.y1 = p.start[1]*scale+translate[1];
		ply_vertex.z1 = p.start[2]*scale+translate[2];
		ply_vertex.x2 = p.end  [0]*scale+translate[0];
		ply_vertex.y2 = p.end  [1]*scale+translate[1];
		ply_vertex.z2 = p.end  [2]*scale+translate[2];
		ply_vertex.s  = p.value;
		ply_put_element(ply, (void *) &ply_vertex);
	}
	for( i=0; i<mesh->outOfCorePointCount() ; i++ )
	{
		CoredMeshData2::Vertex p;
		PlyVertex2 ply_vertex;
		mesh->nextOutOfCorePoint( p );
		ply_vertex.x1 = p.start[0]*scale+translate[0];
		ply_vertex.y1 = p.start[1]*scale+translate[1];
		ply_vertex.z1 = p.start[2]*scale+translate[2];
		ply_vertex.x2 = p.end  [0]*scale+translate[0];
		ply_vertex.y2 = p.end  [1]*scale+translate[1];
		ply_vertex.z2 = p.end  [2]*scale+translate[2];
		ply_vertex.s  = p.value;
		ply_put_element(ply, (void *) &ply_vertex);		
	}  // for, write vertices
	
	// write faces
	std::vector< CoredVertexIndex > polygon;
	ply_put_element_setup( ply , "face" );
	for( i=0 ; i<nr_faces ; i++ )
	{
		//
		// create and fill a struct that the ply code can handle
		//
		PlyFace ply_face;
		mesh->nextPolygon( polygon );
		ply_face.nr_vertices = int( polygon.size() );
		ply_face.vertices = new int[ polygon.size() ];
		for( int i=0 ; i<int(polygon.size()) ; i++ )
			if( polygon[i].inCore ) ply_face.vertices[i] = polygon[i].idx;
			else                    ply_face.vertices[i] = polygon[i].idx + int( mesh->inCorePoints.size() );
		ply_put_element( ply, (void *) &ply_face );
		delete[] ply_face.vertices;
	}  // for, write faces
	
	ply_close( ply );
	return 1;
}
