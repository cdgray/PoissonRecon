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

#define ITERATION_POWER 1.0/3
#define MEMORY_ALLOCATOR_BLOCK_SIZE 1<<12

#define FORCE_UNIT_NORMALS 1
#define USE_DOT_RATIOS 1

#define PAD_SIZE (Real(1.0))

const Real EPSILON=Real(1e-6);
const Real ROUND_EPS=Real(1e-5);
/////////////////////
// SortedTreeNodes //
/////////////////////

SortedTreeNodes::SortedTreeNodes(void){
	nodeCount=NULL;
	treeNodes=NULL;
	maxDepth=0;
}
SortedTreeNodes::~SortedTreeNodes(void){
	if(nodeCount){delete[] nodeCount;}
	if(treeNodes){delete[] treeNodes;}
}
void SortedTreeNodes::set(TreeOctNode& root,const int& setIndex){
	if(nodeCount){delete[] nodeCount;}
	if(treeNodes){delete[] treeNodes;}
	maxDepth=root.maxDepth()+1;
	nodeCount=new int[maxDepth+1];
	treeNodes=new TreeOctNode*[root.nodes()];

	TreeOctNode* temp=root.nextNode();
	int i,cnt=0;
	while(temp){
		treeNodes[cnt++]=temp;
		temp=root.nextNode(temp);
	}
	qsort(treeNodes,cnt,sizeof(const TreeOctNode*),TreeOctNode::CompareForwardPointerDepths);
	for(i=0;i<=maxDepth;i++){nodeCount[i]=0;}
	for(i=0;i<cnt;i++){
		if(setIndex){treeNodes[i]->nodeData.nodeIndex=i;}
		nodeCount[treeNodes[i]->depth()+1]++;
	}
	for(i=1;i<=maxDepth;i++){nodeCount[i]+=nodeCount[i-1];}
}


//////////////////////
// SortedTreeLeaves //
//////////////////////
SortedTreeLeaves::SortedTreeLeaves(void){treeLeaves=NULL;}
SortedTreeLeaves::~SortedTreeLeaves(void){if(treeLeaves){delete[] treeLeaves;}}
void SortedTreeLeaves::set(TreeOctNode& root){
	if(treeLeaves){delete[] treeLeaves;}
	leafCount=root.leaves();
	treeLeaves=new TreeOctNode*[root.leaves()];

	TreeOctNode* temp=root.nextLeaf();
	int cnt=0;
	while(temp){
		treeLeaves[cnt++]=temp;
		temp=root.nextLeaf(temp);
	}
	qsort(treeLeaves,cnt,sizeof(const TreeOctNode*),TreeOctNode::CompareBackwardPointerDepths);
}
void SortedTreeLeaves::set(TreeOctNode& root,const int& maxDepth){
	if(treeLeaves){delete[] treeLeaves;}
	leafCount=root.maxDepthLeaves(maxDepth);
	treeLeaves=new TreeOctNode*[leafCount];

	TreeOctNode* temp=root.nextLeaf();
	int cnt=0;
	while(temp){
		if(temp->depth()<=maxDepth){treeLeaves[cnt++]=temp;}
		temp=root.nextLeaf(temp);
	}
	qsort(treeLeaves,cnt,sizeof(const TreeOctNode*),TreeOctNode::CompareBackwardPointerDepths);
}

//////////////////
// TreeNodeData //
//////////////////
int TreeNodeData::UseIndex=1;
TreeNodeData::TreeNodeData(void){
	if(UseIndex){
		nodeIndex=-1;
		centerWeightContribution=0;
	}
	else{isoNode=NULL;}
	value=0;
}
TreeNodeData::~TreeNodeData(void){
	if(!IsoNodeData::UseAllocator()){if(isoNode){delete isoNode;}}
	isoNode=NULL;
}

/////////////////
// IsoNodeData //
/////////////////
int IsoNodeData::UseAlloc=0;
Allocator<IsoNodeData> IsoNodeData::Allocator;
void IsoNodeData::SetAllocator(const int& blockSize){
	if(blockSize>0){
		UseAlloc=1;
		Allocator.set(blockSize);
	}
	else{UseAlloc=0;}
}
int IsoNodeData::UseAllocator(void){return UseAlloc;}


IsoNodeData::IsoNodeData(void){
	eSegmentCount=0;
	eSegments[0]=eSegments[1]=0;
}
IsoNodeData::~IsoNodeData(void){;}
int IsoNodeData::addEdgeSegment(const int& edgeIndex1,const int& edgeIndex2){
	int faceIndex=Cube::FaceAdjacentToEdges(edgeIndex1,edgeIndex2);
	if(faceIndex<0){return -1;}
	int eCount=edgeCount(faceIndex);
	eSegments[0]|=(long long)(edgeIndex1)<<(4*(faceIndex*2+eCount));
	eSegments[1]|=(long long)(edgeIndex2)<<(4*(faceIndex*2+eCount));
	eCount++;
	eSegmentCount &= ~(3<<(faceIndex*2));
	eSegmentCount |= eCount<<(faceIndex*2);
	return faceIndex;
}
inline int IsoNodeData::edgeCount(const int& faceIndex) const{return (eSegmentCount>>(faceIndex*2))&3;}
inline int IsoNodeData::edgeIndex(const int& faceIndex,const int& e1,const int& e2) const{return int((eSegments[e2]>>(4*(faceIndex*2+e1)))&15);}
////////////
// Octree //
////////////
template<int Degree>
double Octree<Degree>::maxMemoryUsage=0;

template<int Degree>
double Octree<Degree>::MemoryUsage(void){
	double mem=MemoryInfo::Usage()/(1<<20);
	if(mem>maxMemoryUsage){maxMemoryUsage=mem;}
	return mem;
}

template<int Degree>
Octree<Degree>::Octree(void){
	radius=0;
	postNormalSmooth=0;
}

template<int Degree>
void Octree<Degree>::setNodeIndices(TreeOctNode& node,int& idx){
	node.nodeData.nodeIndex=idx;
	idx++;
	if(node.children){for(int i=0;i<Cube::CORNERS;i++){setNodeIndices(node.children[i],idx);}}
}
template<int Degree>
int Octree<Degree>::NonLinearSplatOrientedPoint(TreeOctNode* node,const Point3D<Real>& position,const Point3D<Real>& normal){
	double x,dxdy,dxdydz,dx[DIMENSION][3];
	int i,j,k;
	TreeOctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);
	double width;
	Point3D<Real> center;
	Real w;

	node->centerAndWidth(center,w);
	width=w;
	for(int i=0;i<3;i++){
		x=(center.coords[i]-position.coords[i]-width)/width;
		dx[i][0]=1.125+1.500*x+0.500*x*x;
		x=(center.coords[i]-position.coords[i])/width;
		dx[i][1]=0.750        -      x*x;
		dx[i][2]=1.0-dx[i][1]-dx[i][0];
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			dxdy=dx[0][i]*dx[1][j];
			for(k=0;k<3;k++){
				if(neighbors.neighbors[i][j][k]){
					dxdydz=dxdy*dx[2][k];
					int idx=neighbors.neighbors[i][j][k]->nodeData.nodeIndex;
					if(idx<0){
						Point3D<Real> n;
						n.coords[0]=n.coords[1]=n.coords[2]=0;
						idx=neighbors.neighbors[i][j][k]->nodeData.nodeIndex=int(normals->size());
						normals->push_back(n);
					}
					(*normals)[idx].coords[0]+=Real(normal.coords[0]*dxdydz);
					(*normals)[idx].coords[1]+=Real(normal.coords[1]*dxdydz);
					(*normals)[idx].coords[2]+=Real(normal.coords[2]*dxdydz);
				}
			}
		}
	}
	return 0;
}
template<int Degree>
void Octree<Degree>::NonLinearSplatOrientedPoint(const Point3D<Real>& position,const Point3D<Real>& normal,const int& splatDepth,const Real& samplesPerNode,
												 const int& minDepth,const int& maxDepth){
	double dx;
	Point3D<Real> n;
	TreeOctNode* temp;
	int i;
	double width;
	Point3D<Real> myCenter;
	Real myWidth;
	myCenter.coords[0]=myCenter.coords[1]=myCenter.coords[2]=Real(0.5);
	myWidth=Real(1.0);

	temp=&tree;
	while(temp->depth()<splatDepth){
		if(!temp->children){
			printf("Octree<Degree>::NonLinearSplatOrientedPoint error\n");
			return;
		}
		int cIndex=TreeOctNode::CornerIndex(myCenter,position);
		temp=&temp->children[cIndex];
		myWidth/=2;
		if(cIndex&1){myCenter.coords[0]+=myWidth/2;}
		else		{myCenter.coords[0]-=myWidth/2;}
		if(cIndex&2){myCenter.coords[1]+=myWidth/2;}
		else		{myCenter.coords[1]-=myWidth/2;}
		if(cIndex&4){myCenter.coords[2]+=myWidth/2;}
		else		{myCenter.coords[2]-=myWidth/2;}
	}
	Real alpha,newDepth;
	NonLinearGetSampleDepthAndWeight(temp,position,samplesPerNode,newDepth,alpha);

	if(newDepth<minDepth){newDepth=Real(minDepth);}
	if(newDepth>maxDepth){newDepth=Real(maxDepth);}
	int topDepth=int(ceil(newDepth));

	dx=1.0-(topDepth-newDepth);
	if(topDepth<=minDepth){
		topDepth=minDepth;
		dx=1;
	}
	else if(topDepth>maxDepth){
		topDepth=maxDepth;
		dx=1;
	}
	while(temp->depth()>topDepth){temp=temp->parent;}
	while(temp->depth()<topDepth){
		if(!temp->children){temp->initChildren();}
		int cIndex=TreeOctNode::CornerIndex(myCenter,position);
		temp=&temp->children[cIndex];
		myWidth/=2;
		if(cIndex&1){myCenter.coords[0]+=myWidth/2;}
		else		{myCenter.coords[0]-=myWidth/2;}
		if(cIndex&2){myCenter.coords[1]+=myWidth/2;}
		else		{myCenter.coords[1]-=myWidth/2;}
		if(cIndex&4){myCenter.coords[2]+=myWidth/2;}
		else		{myCenter.coords[2]-=myWidth/2;}
	}
	width=1.0/(1<<temp->depth());
	for(i=0;i<DIMENSION;i++){n.coords[i]=normal.coords[i]*alpha/Real(pow(width,3))*Real(dx);}
	NonLinearSplatOrientedPoint(temp,position,n);
	if(fabs(1.0-dx)>EPSILON){
		dx=Real(1.0-dx);
		temp=temp->parent;
		width=1.0/(1<<temp->depth());

		for(i=0;i<DIMENSION;i++){n.coords[i]=normal.coords[i]*alpha/Real(pow(width,3))*Real(dx);}
		NonLinearSplatOrientedPoint(temp,position,n);
	}
}
template<int Degree>
void Octree<Degree>::NonLinearGetSampleDepthAndWeight(TreeOctNode* node,const Point3D<Real>& position,const Real& samplesPerNode,Real& depth,Real& weight){
	TreeOctNode* temp=node;
	weight=Real(1.0)/NonLinearGetSampleWeight(temp,position);
	if(weight>=samplesPerNode+1){depth=Real(temp->depth()+log(weight/(samplesPerNode+1))/log(double(1<<(DIMENSION-1))));}
	else{
		Real oldAlpha,newAlpha;
		oldAlpha=newAlpha=weight;
		while(newAlpha<(samplesPerNode+1) && temp->parent){
			temp=temp->parent;
			oldAlpha=newAlpha;
			newAlpha=Real(1.0)/NonLinearGetSampleWeight(temp,position);
		}
		depth=Real(temp->depth()+log(newAlpha/(samplesPerNode+1))/log(newAlpha/oldAlpha));
	}
	weight=Real(pow(double(1<<(DIMENSION-1)),-double(depth)));
}

template<int Degree>
Real Octree<Degree>::NonLinearGetSampleWeight(TreeOctNode* node,const Point3D<Real>& position){
	Real weight=0;
	double x,dxdy,dx[DIMENSION][3];
	int i,j,k;
	TreeOctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);
	double width;
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center,w);
	width=w;

	for(i=0;i<DIMENSION;i++){
		x=(center.coords[i]-position.coords[i]-width)/width;
		dx[i][0]=1.125+1.500*x+0.500*x*x;
		x=(center.coords[i]-position.coords[i])/width;
		dx[i][1]=0.750        -      x*x;
		dx[i][2]=1.0-dx[i][1]-dx[i][0];
	}

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			dxdy=dx[0][i]*dx[1][j];
			for(k=0;k<3;k++){
				if(neighbors.neighbors[i][j][k]){weight+=Real(dxdy*dx[2][k]*neighbors.neighbors[i][j][k]->nodeData.centerWeightContribution);}
			}
		}
	}
	return Real(1.0/weight);
}
template<int Degree>
int Octree<Degree>::NonLinearUpdateWeightContribution(TreeOctNode* node,const Point3D<Real>& position){
	int i,j,k;
	TreeOctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);
	double x,dxdy,dx[DIMENSION][3];
	double width;
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center,w);
	width=w;

	for(i=0;i<DIMENSION;i++){
		x=(center.coords[i]-position.coords[i]-width)/width;
		dx[i][0]=1.125+1.500*x+0.500*x*x;
		x=(center.coords[i]-position.coords[i])/width;
		dx[i][1]=0.750        -      x*x;
		dx[i][2]=1.0-dx[i][1]-dx[i][0];
	}

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			dxdy=dx[0][i]*dx[1][j];
			for(k=0;k<3;k++){
				if(neighbors.neighbors[i][j][k]){neighbors.neighbors[i][j][k]->nodeData.centerWeightContribution+=Real(dxdy*dx[2][k]);}
			}
		}
	}
	return 0;
}

template<int Degree>
int Octree<Degree>::setTree(char* fileName,const int& maxDepth,const int& binary,
							const int& kernelDepth,
							const Real& samplesPerNode,const Real& scaleFactor,Point3D<Real>& center,Real& scale,const int& resetSamples){

	FILE* fp;
	Point3D<Real> min,max,position,normal,myCenter;
	Real myWidth;
	int i,cnt=0;
	float c[2*DIMENSION];
	TreeOctNode* temp;
	int splatDepth=0;

	TreeNodeData::UseIndex=1;
	neighborKey.set(maxDepth);
	splatDepth=kernelDepth;
	if(splatDepth<0){splatDepth=0;}
	if(binary){fp=fopen(fileName,"rb");}
	else{fp=fopen(fileName,"r");}
	if(!fp){return 0;}

	DumpOutput("Setting bounding box\n");
	// Read through once to get the center and scale
	while(1){
		if(binary){if(fread(c,sizeof(float),2*DIMENSION,fp)!=6){break;}}
		else{if(fscanf(fp," %f %f %f %f %f %f ",&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])!=2*DIMENSION){break;}}
		for(i=0;i<DIMENSION;i++){
			if(!cnt || c[i]<min.coords[i]){min.coords[i]=c[i];}
			if(!cnt || c[i]>max.coords[i]){max.coords[i]=c[i];}
		}
		cnt++;
	}
	for(i=0;i<DIMENSION;i++){
		if(!i || scale<max.coords[i]-min.coords[i]){scale=Real(max.coords[i]-min.coords[i]);}
		center.coords[i]=Real(max.coords[i]+min.coords[i])/2;
	}
	DumpOutput("Samples: %d\n",cnt);
	scale*=scaleFactor;
	for(i=0;i<DIMENSION;i++){center.coords[i]-=scale/2;}
	if(splatDepth>0){
		DumpOutput("Setting sample weights\n");
		fseek(fp,SEEK_SET,0);
		cnt=0;
		while(1){
			if(binary){if(fread(c,sizeof(float),2*DIMENSION,fp)!=2*DIMENSION){break;}}
			else{if(fscanf(fp," %f %f %f %f %f %f ",&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])!=2*DIMENSION){break;}}
			for(i=0;i<DIMENSION;i++){position.coords[i]=(c[i]-center.coords[i])/scale;}
			myCenter.coords[0]=myCenter.coords[1]=myCenter.coords[2]=Real(0.5);
			myWidth=Real(1.0);
			for(i=0;i<DIMENSION;i++){if(position.coords[i]<myCenter.coords[i]-myWidth/2 || position.coords[i]>myCenter.coords[i]+myWidth/2){break;}}
			if(i!=DIMENSION){continue;}
			temp=&tree;
			int d=0;
			while(d<splatDepth){
				NonLinearUpdateWeightContribution(temp,position);
				if(!temp->children){temp->initChildren();}
				int cIndex=TreeOctNode::CornerIndex(myCenter,position);
				temp=&temp->children[cIndex];
				myWidth/=2;
				if(cIndex&1){myCenter.coords[0]+=myWidth/2;}
				else		{myCenter.coords[0]-=myWidth/2;}
				if(cIndex&2){myCenter.coords[1]+=myWidth/2;}
				else		{myCenter.coords[1]-=myWidth/2;}
				if(cIndex&4){myCenter.coords[2]+=myWidth/2;}
				else		{myCenter.coords[2]-=myWidth/2;}
				d++;
			}
			NonLinearUpdateWeightContribution(temp,position);
			cnt++;
		}
	}

	DumpOutput("Adding Points and Normals\n");
	normals=new std::vector<Point3D<Real> >();
	fseek(fp,SEEK_SET,0);
	cnt=0;
	while(1){
		if(binary){if(fread(c,sizeof(float),2*DIMENSION,fp)!=2*DIMENSION){break;}}
		else{if(fscanf(fp," %f %f %f %f %f %f ",&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])!=2*DIMENSION){break;}}
		for(i=0;i<DIMENSION;i++){
			position.coords[i]=(c[i]-center.coords[i])/scale;
			normal.coords[i]=c[DIMENSION+i];
		}
		myCenter.coords[0]=myCenter.coords[1]=myCenter.coords[2]=Real(0.5);
		myWidth=Real(1.0);
		for(i=0;i<DIMENSION;i++){if(position.coords[i]<myCenter.coords[i]-myWidth/2 || position.coords[i]>myCenter.coords[i]+myWidth/2){break;}}
		if(i!=DIMENSION){continue;}
#if FORCE_UNIT_NORMALS
		Real l=Real(Length(normal));
		if(l>EPSILON){l=1.f/l;}
		else{l=0;}
		l*=(2<<maxDepth);
		normal.coords[0]*=l;
		normal.coords[1]*=l;
		normal.coords[2]*=l;
#endif // FORCE_UNIT_NORMALS
		if(resetSamples && samplesPerNode>0 && splatDepth){
			NonLinearSplatOrientedPoint(position,normal,splatDepth,samplesPerNode,1,maxDepth);
		}
		else{
			Real alpha=1;
			temp=&tree;
			if(splatDepth){
				int d=0;
				while(d<splatDepth){
					int cIndex=TreeOctNode::CornerIndex(myCenter,position);
					temp=&temp->children[cIndex];
					myWidth/=2;
					if(cIndex&1){myCenter.coords[0]+=myWidth/2;}
					else		{myCenter.coords[0]-=myWidth/2;}
					if(cIndex&2){myCenter.coords[1]+=myWidth/2;}
					else		{myCenter.coords[1]-=myWidth/2;}
					if(cIndex&4){myCenter.coords[2]+=myWidth/2;}
					else		{myCenter.coords[2]-=myWidth/2;}
					d++;
				}
				alpha=NonLinearGetSampleWeight(temp,position);

			}
			for(i=0;i<DIMENSION;i++){normal.coords[i]*=alpha;}
			int d=0;
			while(d<maxDepth){
				if(!temp->children){temp->initChildren();}
				int cIndex=TreeOctNode::CornerIndex(myCenter,position);
				temp=&temp->children[cIndex];
				myWidth/=2;
				if(cIndex&1){myCenter.coords[0]+=myWidth/2;}
				else		{myCenter.coords[0]-=myWidth/2;}
				if(cIndex&2){myCenter.coords[1]+=myWidth/2;}
				else		{myCenter.coords[1]-=myWidth/2;}
				if(cIndex&4){myCenter.coords[2]+=myWidth/2;}
				else		{myCenter.coords[2]-=myWidth/2;}
				d++;
			}
			NonLinearSplatOrientedPoint(temp,position,normal);
		}
	}

	DumpOutput("Memory Usage: %.3f MB\n",float(MemoryUsage()));
	fclose(fp);
	return cnt;
}

template<int Degree>
void Octree<Degree>::setFunctionData(const PPolynomial<Degree>& ReconstructionFunction,	const int& maxDepth,const int& normalize,const Real& normalSmooth){

	radius=Real(fabs(ReconstructionFunction.polys[0].start));
	if(normalSmooth>0){postNormalSmooth=normalSmooth;}
	fData.set(maxDepth,ReconstructionFunction,normalize,USE_DOT_RATIOS);
}

template<int Degree>
void Octree<Degree>::finalize1(const int& refineNeighbors)
{
	TreeOctNode* temp;

	if(refineNeighbors>=0){
		RefineFunction rf;
		temp=tree.nextNode();
		while(temp){
			if(temp->nodeData.nodeIndex>=0 && Length((*normals)[temp->nodeData.nodeIndex])>EPSILON){
				rf.depth=temp->depth()-refineNeighbors;
				TreeOctNode::ProcessMaxDepthNodeAdjacentNodes(temp,2*radius,&tree,Real(0.5),temp->depth()-refineNeighbors,&rf);
			}
			temp=tree.nextNode(temp);
		}
	}
	else if(refineNeighbors==-1234){
		temp=tree.nextLeaf();
		while(temp){
			if(!temp->children && temp->depth()<fData.depth){temp->initChildren();}
			temp=tree.nextLeaf(temp);
		}
	}
}
template<int Degree>
void Octree<Degree>::finalize2(const int& refineNeighbors)
{
	TreeOctNode* temp;

	if(refineNeighbors>=0){
		RefineFunction rf;
		temp=tree.nextNode();
		while(temp){
			if(fabs(temp->nodeData.value)>EPSILON){
				rf.depth=temp->depth()-refineNeighbors;
				TreeOctNode::ProcessMaxDepthNodeAdjacentNodes(temp,2*radius,&tree,Real(0.5),temp->depth()-refineNeighbors,&rf);
			}
			temp=tree.nextNode(temp);
		}
	}
}
template <int Degree>
void Octree<Degree>::SetDivergence(const int idx[DIMENSION],Point3D<Real>& div) const
{
#if USE_DOT_RATIOS
	double dot=fData.dotTable[idx[0]]*fData.dotTable[idx[1]]*fData.dotTable[idx[2]];
	div.coords[0]=fData.dDotTable[idx[0]]*dot;
	div.coords[1]=fData.dDotTable[idx[1]]*dot;
	div.coords[2]=fData.dDotTable[idx[2]]*dot;
#else // !USE_DOT_RATIOS
	double dot[DIMENSION];
	dot[0]=fData.dotTable[idx[0]];
	dot[1]=fData.dotTable[idx[1]];
	dot[2]=fData.dotTable[idx[2]];
	return Real(
		fData.dDotTable[idx[0]]*normal.coords[0]*dot[1]*dot[2]+
		fData.dDotTable[idx[1]]*normal.coords[1]*dot[0]*dot[2]+
		fData.dDotTable[idx[2]]*normal.coords[2]*dot[0]*dot[1]
	);
#endif // !USE_DOT_RATIOS
}

template <int Degree>
Real Octree<Degree>::GetDivergence(const int idx[DIMENSION],const Point3D<Real>& normal) const
{
#if USE_DOT_RATIOS
	double dot=fData.dotTable[idx[0]]*fData.dotTable[idx[1]]*fData.dotTable[idx[2]];
	return Real(dot*(fData.dDotTable[idx[0]]*normal.coords[0]+fData.dDotTable[idx[1]]*normal.coords[1]+fData.dDotTable[idx[2]]*normal.coords[2]));
#else // !USE_DOT_RATIOS
	double dot[DIMENSION];
	dot[0]=fData.dotTable[idx[0]];
	dot[1]=fData.dotTable[idx[1]];
	dot[2]=fData.dotTable[idx[2]];
	return Real(
		fData.dDotTable[idx[0]]*normal.coords[0]*dot[1]*dot[2]+
		fData.dDotTable[idx[1]]*normal.coords[1]*dot[0]*dot[2]+
		fData.dDotTable[idx[2]]*normal.coords[2]*dot[0]*dot[1]
	);
#endif // !USE_DOT_RATIOS
}
template<int Degree>
Real Octree<Degree>::GetLaplacian(const int idx[DIMENSION]) const
{
#if USE_DOT_RATIOS
	return Real(fData.dotTable[idx[0]]*fData.dotTable[idx[1]]*fData.dotTable[idx[2]]*(fData.d2DotTable[idx[0]]+fData.d2DotTable[idx[1]]+fData.d2DotTable[idx[2]]));
#else // !USE_DOT_RATIOS
	double dot[3];
	dot[0]=fData.dotTable[idx[0]];
	dot[1]=fData.dotTable[idx[1]];
	dot[2]=fData.dotTable[idx[2]];
	return Real(
		fData.d2DotTable[idx[0]]*dot[1]*dot[2]+
		fData.d2DotTable[idx[1]]*dot[0]*dot[2]+
		fData.d2DotTable[idx[2]]*dot[0]*dot[1]
	);
#endif // USE_DOT_RATIOS
}
template<int Degree>
Real Octree<Degree>::GetDotProduct(const int idx[DIMENSION]) const
{
	return Real(fData.dotTable[idx[0]]*fData.dotTable[idx[1]]*fData.dotTable[idx[2]]);
}
template<int Degree>
int Octree<Degree>::GetFixedDepthLaplacian(SparseSymmetricMatrix<float>& matrix,const int& depth,const SortedTreeNodes& sNodes)
{
	LaplacianMatrixFunction mf;
	mf.ot=this;
	mf.offset=sNodes.nodeCount[depth];
	Real myRadius=int(2*radius-ROUND_EPS)+ROUND_EPS;
	matrix.Resize(sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth]);
	mf.rowElements=(MatrixEntry<float>*)malloc(sizeof(MatrixEntry<float>)*matrix.rows);
	for(int i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];i++){
		mf.elementCount=0;
		mf.d2=int(sNodes.treeNodes[i]->d);
		mf.x2=int(sNodes.treeNodes[i]->off[0]);
		mf.y2=int(sNodes.treeNodes[i]->off[1]);
		mf.z2=int(sNodes.treeNodes[i]->off[2]);
		mf.index[0]=mf.x2*fData.res;
		mf.index[1]=mf.y2*fData.res;
		mf.index[2]=mf.z2*fData.res;
		TreeOctNode::ProcessTerminatingNodeAdjacentNodes(sNodes.treeNodes[i],myRadius-Real(0.5),&tree,Real(0.5),&mf);
		matrix.SetRowSize(i-sNodes.nodeCount[depth],mf.elementCount);
		memcpy(matrix.m_ppElements[i-sNodes.nodeCount[depth]],mf.rowElements,sizeof(MatrixEntry<float>)*mf.elementCount);
	}
	free(mf.rowElements);
	return 1;
}
template<int Degree>
int Octree<Degree>::GetRestrictedFixedDepthLaplacian(SparseSymmetricMatrix<float>& matrix,const int* entries,const int& entryCount,
													 const TreeOctNode* rNode,const Real& radius,
													 const SortedTreeNodes& sNodes){
	int i;
	RestrictedLaplacianMatrixFunction mf;
	Real myRadius=int(2*radius-ROUND_EPS)+ROUND_EPS;
	mf.ot=this;
	mf.radius=radius;
	rNode->depthAndOffset(mf.depth,mf.offset);
	matrix.Resize(entryCount);
	mf.rowElements=(MatrixEntry<float>*)malloc(sizeof(MatrixEntry<float>)*matrix.rows);
	for(i=0;i<entryCount;i++){sNodes.treeNodes[entries[i]]->nodeData.nodeIndex=i;}
	for(i=0;i<entryCount;i++){
		mf.elementCount=0;
		mf.index[0]=int(sNodes.treeNodes[entries[i]]->off[0])*fData.res;
		mf.index[1]=int(sNodes.treeNodes[entries[i]]->off[1])*fData.res;
		mf.index[2]=int(sNodes.treeNodes[entries[i]]->off[2])*fData.res;
		TreeOctNode::ProcessTerminatingNodeAdjacentNodes(sNodes.treeNodes[entries[i]],myRadius-Real(0.5),&tree,Real(0.5),&mf);
		matrix.SetRowSize(i,mf.elementCount);
		memcpy(matrix.m_ppElements[i],mf.rowElements,sizeof(MatrixEntry<float>)*mf.elementCount);
	}
	for(i=0;i<entryCount;i++){sNodes.treeNodes[entries[i]]->nodeData.nodeIndex=entries[i];}
	free(mf.rowElements);
	return 1;
}

template<int Degree>
int Octree<Degree>::LaplacianMatrixIteration(const int& subdivideDepth){
	int i,iter=0;
	SortedTreeNodes sNodes;
	double t;
	fData.setDotTables(fData.D2_DOT_FLAG);
	sNodes.set(tree,1);

	SparseMatrix<float>::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);

	sNodes.treeNodes[0]->nodeData.value=0;
	for(i=1;i<sNodes.maxDepth;i++){
		DumpOutput("Depth: %d/%d\n",i,sNodes.maxDepth-1);
		t=Time();
		if(subdivideDepth>0){iter+=SolveFixedDepthMatrix(i,subdivideDepth,sNodes);}
		else{iter+=SolveFixedDepthMatrix(i,sNodes);}
	}
	SparseMatrix<float>::Allocator.reset();
	fData.clearDotTables(fData.DOT_FLAG | fData.D_DOT_FLAG | fData.D2_DOT_FLAG);
	return iter;
}

template<int Degree>
int Octree<Degree>::SolveFixedDepthMatrix(const int& depth,const SortedTreeNodes& sNodes){
	int i,iter=0;
	Vector<double> V,Solution;
	SparseSymmetricMatrix<Real> matrix;
	Real myRadius,myRadius1,myRadius2;
	double gTime,sTime,uTime;
	Real dx,dy,dz;
	int x1,x2,y1,y2,z1,z2;


	gTime=Time();
	V.Resize(sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth]);
	for(i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];i++){V[i-sNodes.nodeCount[depth]]=sNodes.treeNodes[i]->nodeData.value;}
	SparseSymmetricMatrix<float>::Allocator.rollBack();
	GetFixedDepthLaplacian(matrix,depth,sNodes);
	gTime=Time()-gTime;
	DumpOutput("\tMatrix entries: %d / %d^2 = %.4f%%\n",matrix.Entries(),matrix.rows,100.0*(matrix.Entries()/double(matrix.rows))/matrix.rows);
	DumpOutput("\tMemory Usage: %.3f MB\n",float(MemoryUsage()));
	sTime=Time();
	iter+=SparseSymmetricMatrix<Real>::Solve(matrix,V,int(pow(matrix.rows,ITERATION_POWER)),Solution,double(EPSILON),1);
	sTime=Time()-sTime;
	uTime=Time();
	for(i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];i++){sNodes.treeNodes[i]->nodeData.value    =Real(Solution[i-sNodes.nodeCount[depth]]);}

	myRadius=Real(radius+ROUND_EPS-0.5);
	myRadius1=myRadius2=radius;
	myRadius1=Real(int(2*myRadius1+1-ROUND_EPS))/2-ROUND_EPS;
	myRadius2=Real(int(2*myRadius2  -ROUND_EPS))/2+ROUND_EPS;

	myRadius /=(1<<depth);
	myRadius1/=(1<<depth);
	myRadius2/=(1<<depth);

	if(depth<sNodes.maxDepth-1){
		LaplacianProjectionFunction pf;
		TreeOctNode *node1,*node2;
		pf.ot=this;
		int idx1,idx2,off=sNodes.nodeCount[depth];
		// First pass: idx2 is the solution coefficient propogated
		for(i=0;i<matrix.rows;i++){
			idx1=i;
			node1=sNodes.treeNodes[idx1+off];
			if(!node1->children){continue;}
			x1=int(node1->off[0]);
			y1=int(node1->off[1]);
			z1=int(node1->off[2]);
			for(int j=0;j<matrix.rowSizes[i];j++){
				idx2=matrix.m_ppElements[i][j].N;
				node2=sNodes.treeNodes[idx2+off];
				x2=int(node2->off[0]);
				y2=int(node2->off[1]);
				z2=int(node2->off[2]);
				pf.value=Solution[idx2];
				pf.index[0]=fData.res*x2;
				pf.index[1]=fData.res*y2;
				pf.index[2]=fData.res*z2;
				dx=Real(x2-x1)/(1<<depth);
				dy=Real(y2-y1)/(1<<depth);
				dz=Real(z2-z1)/(1<<depth);
				if(fabs(dx)<myRadius && fabs(dy)<myRadius && fabs(dz)<myRadius){node1->processNodeNodes(node2,&pf,0);}
				else{TreeOctNode::ProcessNodeAdjacentNodes(dx,dy,dz,node2,myRadius1,node1,myRadius2,Real(1.0/(1<<depth)),&pf,0);}
			}
		}
		// Second pass: idx1 is the solution coefficient propogated
		for(i=0;i<matrix.rows;i++){
			idx1=i;
			node1=sNodes.treeNodes[idx1+off];
			x1=int(node1->off[0]);
			y1=int(node1->off[1]);
			z1=int(node1->off[2]);
			pf.value=Solution[idx1];
			pf.index[0]=fData.res*x1;
			pf.index[1]=fData.res*y1;
			pf.index[2]=fData.res*z1;
			for(int j=0;j<matrix.rowSizes[i];j++){
				idx2=matrix.m_ppElements[i][j].N;
				node2=sNodes.treeNodes[idx2+off];
				if(idx1!=idx2 && node2->children){
					x2=int(node2->off[0]);
					y2=int(node2->off[1]);
					z2=int(node2->off[2]);
					dx=Real(x1-x2)/(1<<depth);
					dy=Real(y1-y2)/(1<<depth);
					dz=Real(z1-z2)/(1<<depth);
					if(fabs(dx)<myRadius && fabs(dy)<myRadius && fabs(dz)<myRadius){node2->processNodeNodes(node1,&pf,0);}
					else{TreeOctNode::ProcessNodeAdjacentNodes(dx,dy,dz,node1,myRadius1,node2,myRadius2,Real(1.0/(1<<depth)),&pf,0);}
				}
			}
		}
	}
	uTime=Time()-uTime;
	DumpOutput("\tGot / Solved / Updated in: %6.3f / %6.3f / %6.3f\n",gTime,sTime,uTime);
	return iter;
}
template<int Degree>
int Octree<Degree>::SolveFixedDepthMatrix(const int& depth,const int& startingDepth,const SortedTreeNodes& sNodes){
	int i,j,d,iter=0;
	SparseSymmetricMatrix<Real> matrix;
	AdjacencySetFunction asf;
	AdjacencyCountFunction acf;
	Vector<Real> Values;
	Vector<double> SubValues,SubSolution;
	double gTime,sTime,uTime;
	Real myRadius,myRadius2;
	Real dx,dy,dz;
	Real sRadius=radius/(1<<depth);

	if(startingDepth>=depth){return SolveFixedDepthMatrix(depth,sNodes);}

	Values.Resize(sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth]);

	for(i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];i++){
		Values[i-sNodes.nodeCount[depth]]=sNodes.treeNodes[i]->nodeData.value;
		sNodes.treeNodes[i]->nodeData.value=0;
	}

	myRadius=2*radius-Real(0.5);
	myRadius=int(myRadius-ROUND_EPS)+ROUND_EPS;
	myRadius2=Real(radius+ROUND_EPS-0.5);
	d=depth-startingDepth;
	for(i=sNodes.nodeCount[d];i<sNodes.nodeCount[d+1];i++){
		gTime=Time();
		TreeOctNode* temp;
		// Get all of the entries associated to the subspace
		acf.adjacencyCount=0;
		temp=sNodes.treeNodes[i]->nextNode();
		while(temp){
			if(temp->depth()==depth){
				acf.Function(temp,temp);
				temp=sNodes.treeNodes[i]->nextBranch(temp);
			}
			else{temp=sNodes.treeNodes[i]->nextNode(temp);}
		}
		for(j=sNodes.nodeCount[d];j<sNodes.nodeCount[d+1];j++){
			if(i==j){continue;}
			TreeOctNode::ProcessFixedDepthNodeAdjacentNodes(sNodes.treeNodes[i],Real(0.5),sNodes.treeNodes[j],myRadius,depth,&acf);
		}
		if(!acf.adjacencyCount){continue;}
		asf.adjacencies=new int[acf.adjacencyCount];
		asf.adjacencyCount=0;
		temp=sNodes.treeNodes[i]->nextNode();
		while(temp){
			if(temp->depth()==depth){
				asf.Function(temp,temp);
				temp=sNodes.treeNodes[i]->nextBranch(temp);
			}
			else{temp=sNodes.treeNodes[i]->nextNode(temp);}
		}
		for(j=sNodes.nodeCount[d];j<sNodes.nodeCount[d+1];j++){
			if(i==j){continue;}
			TreeOctNode::ProcessFixedDepthNodeAdjacentNodes(sNodes.treeNodes[i],Real(0.5),sNodes.treeNodes[j],myRadius,depth,&asf);
		}

		DumpOutput("\tNodes[%d/%d]: %d\n",i-sNodes.nodeCount[d]+1,sNodes.nodeCount[d+1]-sNodes.nodeCount[d],asf.adjacencyCount);
		// Get the associated vector
		SubValues.Resize(asf.adjacencyCount);
		for(j=0;j<asf.adjacencyCount;j++){SubValues[j]=Values[asf.adjacencies[j]-sNodes.nodeCount[depth]];}
		SubSolution.Resize(asf.adjacencyCount);
		for(j=0;j<asf.adjacencyCount;j++){SubSolution[j]=sNodes.treeNodes[asf.adjacencies[j]]->nodeData.value;}
		// Get the associated matrix
		SparseSymmetricMatrix<float>::Allocator.rollBack();
		GetRestrictedFixedDepthLaplacian(matrix,asf.adjacencies,asf.adjacencyCount,sNodes.treeNodes[i],myRadius,sNodes);
		gTime=Time()-gTime;
		DumpOutput("\t\tMatrix entries: %d / %d^2 = %.4f%%\n",matrix.Entries(),matrix.rows,100.0*(matrix.Entries()/double(matrix.rows))/matrix.rows);
		DumpOutput("\t\tMemory Usage: %.3f MB\n",float(MemoryUsage()));

		// Solve the matrix
		sTime=Time();
		iter+=SparseSymmetricMatrix<Real>::Solve(matrix,SubValues,int(pow(matrix.rows,ITERATION_POWER)),SubSolution,double(EPSILON),0);
		sTime=Time()-sTime;

		uTime=Time();
		LaplacianProjectionFunction lpf;
		lpf.ot=this;

		// Update the solution for all nodes in the sub-tree
		for(j=0;j<asf.adjacencyCount;j++){
			temp=sNodes.treeNodes[asf.adjacencies[j]];
			while(temp->depth()>sNodes.treeNodes[i]->depth()){temp=temp->parent;}
			if(temp->nodeData.nodeIndex>=sNodes.treeNodes[i]->nodeData.nodeIndex){sNodes.treeNodes[asf.adjacencies[j]]->nodeData.value=Real(SubSolution[j]);}
		}

		// Update the values in the next depth
		int x1,x2,y1,y2,z1,z2;
		if(depth<sNodes.maxDepth-1){
			int idx1,idx2;
			TreeOctNode *node1,*node2;
			// First pass: idx2 is the solution coefficient propogated
			for(j=0;j<matrix.rows;j++){
				idx1=asf.adjacencies[j];
				node1=sNodes.treeNodes[idx1];
				if(!node1->children){continue;}
				x1=int(node1->off[0]);
				y1=int(node1->off[1]);
				z1=int(node1->off[2]);

				for(int k=0;k<matrix.rowSizes[j];k++){
					idx2=asf.adjacencies[matrix.m_ppElements[j][k].N];
					node2=sNodes.treeNodes[idx2];
					temp=node2;
					while(temp->depth()>d){temp=temp->parent;}
					if(temp!=sNodes.treeNodes[i]){continue;}
					lpf.value=Real(SubSolution[matrix.m_ppElements[j][k].N]);
					x2=int(node2->off[0]);
					y2=int(node2->off[1]);
					z2=int(node2->off[2]);
					lpf.index[0]=fData.res*x2;
					lpf.index[1]=fData.res*y2;
					lpf.index[2]=fData.res*z2;
					dx=Real(x2-x1)/(1<<depth);
					dy=Real(y2-y1)/(1<<depth);
					dz=Real(z2-z1)/(1<<depth);
					if(fabs(dx)<myRadius2 && fabs(dy)<myRadius2 && fabs(dz)<myRadius2){node1->processNodeNodes(node2,&lpf,0);}
					else{TreeOctNode::ProcessNodeAdjacentNodes(dx,dy,dz,node2,sRadius,node1,sRadius,Real(1.0/(1<<depth)),&lpf,0);}
				}
			}
			// Second pass: idx1 is the solution coefficient propogated
			for(j=0;j<matrix.rows;j++){
				idx1=asf.adjacencies[j];
				node1=sNodes.treeNodes[idx1];
				temp=node1;
				while(temp->depth()>d){temp=temp->parent;}
				if(temp!=sNodes.treeNodes[i]){continue;}
				x1=int(node1->off[0]);
				y1=int(node1->off[1]);
				z1=int(node1->off[2]);

				lpf.value=Real(SubSolution[j]);
				lpf.index[0]=fData.res*x1;
				lpf.index[1]=fData.res*y1;
				lpf.index[2]=fData.res*z1;
				for(int k=0;k<matrix.rowSizes[j];k++){
					idx2=asf.adjacencies[matrix.m_ppElements[j][k].N];
					node2=sNodes.treeNodes[idx2];
					if(!node2->children){continue;}

					if(idx1!=idx2){
						x2=int(node2->off[0]);
						y2=int(node2->off[1]);
						z2=int(node2->off[2]);
						dx=Real(x1-x2)/(1<<depth);
						dy=Real(y1-y2)/(1<<depth);
						dz=Real(z1-z2)/(1<<depth);
						if(fabs(dx)<myRadius2 && fabs(dy)<myRadius2 && fabs(dz)<myRadius2){node2->processNodeNodes(node1,&lpf,0);}
						else{TreeOctNode::ProcessNodeAdjacentNodes(dx,dy,dz,node1,sRadius,node2,sRadius,Real(1.0/(1<<depth)),&lpf,0);}
					}
				}
			}
		}
		uTime=Time()-uTime;
		DumpOutput("\t\tGot / Solved / Updated in: %6.3f / %6.3f / %6.3f\n",gTime,sTime,uTime);
		delete[] asf.adjacencies;
	}

	return iter;
}
template<int Degree>
int Octree<Degree>::HasNormals(TreeOctNode* node,const Real& epsilon){
	int hasNormals=0;
	if(node->nodeData.nodeIndex>=0 && Length((*normals)[node->nodeData.nodeIndex])>epsilon){hasNormals=1;}
	if(node->children){for(int i=0;i<Cube::CORNERS && !hasNormals;i++){hasNormals|=HasNormals(&node->children[i],epsilon);}}

	return hasNormals;
}
template<int Degree>
void Octree<Degree>::ClipTree(void){
	TreeOctNode* temp;
	temp=tree.nextNode();
	while(temp){
		if(temp->children){
			int hasNormals=0;
			for(int i=0;i<Cube::CORNERS && !hasNormals;i++){hasNormals=HasNormals(&temp->children[i],EPSILON);}
			if(!hasNormals){temp->children=NULL;}
		}
		temp=tree.nextNode(temp);
	}
}

template<int Degree>
void Octree<Degree>::SetLaplacianWeights(void){
	TreeOctNode* temp;

	fData.setDotTables(fData.DOT_FLAG | fData.D_DOT_FLAG);
	DivergenceFunction df;
	df.ot=this;
	temp=tree.nextNode();
	while(temp){
		if(temp->nodeData.nodeIndex<0 || Length((*normals)[temp->nodeData.nodeIndex])<=EPSILON){
			temp=tree.nextNode(temp);
			continue;
		}
		df.normal=(*normals)[temp->nodeData.nodeIndex];
		df.index[0]=int(temp->off[0])*fData.res;
		df.index[1]=int(temp->off[1])*fData.res;
		df.index[2]=int(temp->off[2])*fData.res;
		TreeOctNode::ProcessNodeAdjacentNodes(temp,radius,&tree,radius,&df);
		temp=tree.nextNode(temp);
	}
	fData.clearDotTables(fData.D_DOT_FLAG);
	temp=tree.nextNode();
	while(temp){
		if(temp->nodeData.nodeIndex<0){temp->nodeData.centerWeightContribution=0;}
		else{temp->nodeData.centerWeightContribution=Real(Length((*normals)[temp->nodeData.nodeIndex]));}
		temp=tree.nextNode(temp);
	}
	MemoryUsage();

	delete normals;
	normals=NULL;
}
template<int Degree>
void Octree<Degree>::DivergenceFunction::Function(TreeOctNode* node1,const TreeOctNode*){
	scratch[0]=index[0]+int(node1->off[0]);
	scratch[1]=index[1]+int(node1->off[1]);
	scratch[2]=index[2]+int(node1->off[2]);
	node1->nodeData.value+=ot->GetDivergence(scratch,normal);
}
template<int Degree>
void Octree<Degree>::LaplacianProjectionFunction::Function(TreeOctNode* node1,const TreeOctNode*){
	scratch[0]=index[0]+int(node1->off[0]);
	scratch[1]=index[1]+int(node1->off[1]);
	scratch[2]=index[2]+int(node1->off[2]);
	node1->nodeData.value-=Real(ot->GetLaplacian(scratch)*value);
}
template<int Degree>
void Octree<Degree>::AdjacencyCountFunction::Function(const TreeOctNode*,const TreeOctNode*){adjacencyCount++;}
template<int Degree>
void Octree<Degree>::AdjacencySetFunction::Function(const TreeOctNode* node1,const TreeOctNode*){adjacencies[adjacencyCount++]=node1->nodeData.nodeIndex;}
template<int Degree>
void Octree<Degree>::RefineFunction::Function(TreeOctNode* node1,const TreeOctNode*){
	if(!node1->children && node1->depth()<depth){node1->initChildren();}
}
template<int Degree>
void Octree<Degree>::FaceEdgesFunction::Function(const TreeOctNode* node1,const TreeOctNode* node2){
	if(node1->nodeData.isoNode){
		CoredEdgeIndex e;
		hash_map<long long,int>::iterator rootIter;
		int cnt=node1->nodeData.isoNode->edgeCount(fIndex);
		for(int i=0;i<cnt;i++){
			if(	GetRootIndex(node1,node1->nodeData.isoNode->edgeIndex(fIndex,i,0),maxDepth,*boundaryRoots,interiorRoots,e.idx[1]) &&
				GetRootIndex(node1,node1->nodeData.isoNode->edgeIndex(fIndex,i,1),maxDepth,*boundaryRoots,interiorRoots,e.idx[0])){
					edges->push_back(e);
			}
//			else{fprintf(stderr,"Bad Edge: %d %d\n",e.idx[0],e.idx[1]);}
		}
	}
}
template<int Degree>
void Octree<Degree>::PointIndexValueFunction::Function(const TreeOctNode* node){
	int idx[DIMENSION];
	idx[0]=index[0]+int(node->off[0])*res2;
	idx[1]=index[1]+int(node->off[1])*res2;
	idx[2]=index[2]+int(node->off[2])*res2;
	value+=node->nodeData.value*   Real( valueTables[idx[0]]* valueTables[idx[1]]* valueTables[idx[2]]);
}
template<int Degree>
void Octree<Degree>::PointIndexNormalFunction::Function(const TreeOctNode* node){
	int idx[DIMENSION];
	idx[0]=index[0]+int(node->off[0])*res2;
	idx[1]=index[1]+int(node->off[1])*res2;
	idx[2]=index[2]+int(node->off[2])*res2;
	normal.coords[0]+=	node->nodeData.value*   Real(dValueTables[idx[0]]* valueTables[idx[1]]* valueTables[idx[2]]);
	normal.coords[1]+=	node->nodeData.value*   Real( valueTables[idx[0]]*dValueTables[idx[1]]* valueTables[idx[2]]);
	normal.coords[2]+=	node->nodeData.value*   Real( valueTables[idx[0]]* valueTables[idx[1]]*dValueTables[idx[2]]);
}
template<int Degree>
void Octree<Degree>::NodeMinMaxValueFunction::Function(const TreeOctNode* node1,const TreeOctNode* node2){
	double mn[DIMENSION],mx[DIMENSION];
	int idx[DIMENSION];
	VertexData::CenterIndex(node1,maxDepth,idx);
	for(int i=0;i<DIMENSION;i++){
		int off=int(node1->off[i])*res2;
		int idx1=off+cIndex[i][0];
		int idx2=off+cIndex[i][1];
		if(idx[i]<=cIndex[i][0]){
			mn[i]=valueTable[idx2];
			mx[i]=valueTable[idx1];
		}
		else if(idx[i]>=cIndex[i][1]){
			mn[i]=valueTable[idx1];
			mx[i]=valueTable[idx2];
		}
		else{
			mx[i]=valueTable[off+idx[i]];
			if(valueTable[idx1]<valueTable[idx2])	{mn[i]=valueTable[idx1];}
			else									{mn[i]=valueTable[idx2];}
		}
	}
	if(node1->nodeData.value>0){
		min+=Real(node1->nodeData.value*mn[0]*mn[1]*mn[2]);
		max+=Real(node1->nodeData.value*mx[0]*mx[1]*mx[2]);
	}
	else{
		max+=Real(node1->nodeData.value*mn[0]*mn[1]*mn[2]);
		min+=Real(node1->nodeData.value*mx[0]*mx[1]*mx[2]);
	}
}

template<int Degree>
int Octree<Degree>::LaplacianMatrixFunction::Function(const TreeOctNode* node1,const TreeOctNode* node2){
	Real temp;
	int d1=int(node1->d);
	int x1,y1,z1;
	x1=int(node1->off[0]);
	y1=int(node1->off[1]);
	z1=int(node1->off[2]);
	int dDepth=d2-d1;
	int d;
	d=(x2>>dDepth)-x1;
	if(d<0){return 0;}
	if(!dDepth){
		if(!d){
			d=y2-y1;
			if(d<0){return 0;}
			else if(!d){
				d=z2-z1;
				if(d<0){return 0;}
			}
		}
		scratch[0]=index[0]+x1;
		scratch[1]=index[1]+y1;
		scratch[2]=index[2]+z1;
		temp=ot->GetLaplacian(scratch);
		if(node1==node2){temp/=2;}
		if(fabs(temp)>EPSILON){
			rowElements[elementCount].Value=temp;
			rowElements[elementCount].N=node1->nodeData.nodeIndex-offset;
			elementCount++;
		}
		return 0;
	}
	return 1;
}
template<int Degree>
int Octree<Degree>::RestrictedLaplacianMatrixFunction::Function(const TreeOctNode* node1,const TreeOctNode* node2){
	int d1,d2,off1[3],off2[3];
	node1->depthAndOffset(d1,off1);
	node2->depthAndOffset(d2,off2);
	int dDepth=d2-d1;
	int d;
	d=(off2[0]>>dDepth)-off1[0];
	if(d<0){return 0;}

	if(!dDepth){
		if(!d){
			d=off2[1]-off1[1];
			if(d<0){return 0;}
			else if(!d){
				d=off2[2]-off1[2];
				if(d<0){return 0;}
			}
		}
		// Since we are getting the restricted matrix, we don't want to propogate out to terms that don't contribute...
		if(!TreeOctNode::Overlap2(depth,offset,0.5,d1,off1,radius)){return 0;}
		scratch[0]=index[0]+BinaryNode<Real>::Index(d1,off1[0]);
		scratch[1]=index[1]+BinaryNode<Real>::Index(d1,off1[1]);
		scratch[2]=index[2]+BinaryNode<Real>::Index(d1,off1[2]);
		Real temp=ot->GetLaplacian(scratch);
		if(node1==node2){temp/=2;}
		if(fabs(temp)>EPSILON){
			rowElements[elementCount].Value=temp;
			rowElements[elementCount].N=node1->nodeData.nodeIndex;
			elementCount++;
		}
		return 0;
	}
	return 1;
}
template<int Degree>
void Octree<Degree>::GetMCIsoTriangles(const Real& isoValue,CoredMeshData* mesh,const int& fullDepthIso){
	double t;
	TreeOctNode* temp;
	SortedTreeLeaves sLeaves;
	hash_map<long long,int> roots;
	hash_map<long long,Point3D<Real> > *normalHash=new hash_map<long long,Point3D<Real> >();

	SetIsoSurfaceCorners(isoValue,0,fullDepthIso);
	// At the point all of the corner values have been set and all nodes are valid. Now it's just a matter
	// of running marching cubes.

	t=Time();
	fData.setValueTables(fData.VALUE_FLAG | fData.D_VALUE_FLAG,postNormalSmooth);
	temp=tree.nextLeaf();
	while(temp){
		SetMCRootPositions(temp,0,isoValue,roots,NULL,*normalHash,NULL,NULL,mesh,1);
		temp=tree.nextLeaf(temp);
	}
	MemoryUsage();

	DumpOutput("Normal Size: %.2f MB\n",double(sizeof(Point3D<Real>)*normalHash->size())/1000000);
	DumpOutput("Set %d root positions in: %f\n",mesh->inCorePoints.size(),Time()-t);
	DumpOutput("Memory Usage: %.3f MB\n",float(MemoryUsage()));

	fData.clearValueTables();
	delete normalHash;

	DumpOutput("Post deletion size: %.3f MB\n",float(MemoryUsage()));

	temp=tree.nextNode();
	while(temp){
		if(temp->nodeData.isoNode){
			temp->nodeData.isoNode->mcIndex=MarchingCubes::GetIndex(temp->nodeData.isoNode->cornerValues,isoValue);
			temp->nodeData.isoNode->eSegmentCount=0;
			temp->nodeData.isoNode->eSegments[0]=temp->nodeData.isoNode->eSegments[1]=0;
		}
		temp=tree.nextNode(temp);
	}

	t=Time();

	// Now get the iso-surfaces, running from finest nodes to coarsest in order to allow for edge propogation from
	// finer faces to coarser ones.
	sLeaves.set(tree);

	for(int i=0;i<sLeaves.leafCount;i++){GetMCIsoTriangles(sLeaves.treeLeaves[i],mesh,roots,NULL,NULL,0,0);}
	DumpOutput("Added triangles in: %f\n",Time()-t);
	DumpOutput("Memory Usage: %.3f MB\n",float(MemoryUsage()));
}
template<int Degree>
void Octree<Degree>::GetMCIsoTriangles(const Real& isoValue,const int& subdivideDepth,CoredMeshData* mesh,const int& fullDepthIso){
	TreeOctNode* temp;
	hash_map<long long,int> boundaryRoots,*interiorRoots;
	hash_map<long long,Point3D<Real> > *boundaryNormalHash,*interiorNormalHash;
	std::vector<Point3D<float> >* interiorPoints;

	int sDepth;
	if(subdivideDepth<=0){sDepth=0;}
	else{sDepth=fData.depth-subdivideDepth;}
	if(sDepth<0){sDepth=0;}

	SetIsoSurfaceCorners(isoValue,sDepth,fullDepthIso);
	// At this point all of the corner values have been set and all nodes are valid. Now it's just a matter
	// of running marching cubes.

	boundaryNormalHash=new hash_map<long long,Point3D<Real> >();
	int offSet=0;
	SortedTreeNodes sNodes;
	sNodes.set(tree,0);
	fData.setValueTables(fData.VALUE_FLAG | fData.D_VALUE_FLAG,postNormalSmooth);

	// Set the root positions for all leaf nodes below the subdivide threshold
	SetBoundaryMCRootPositions(sDepth,isoValue,boundaryRoots,*boundaryNormalHash,mesh,1);

	for(int i=sNodes.nodeCount[sDepth];i<sNodes.nodeCount[sDepth+1];i++){
		SortedTreeLeaves sLeaves;
		interiorRoots=new hash_map<long long,int>();
		interiorNormalHash=new hash_map<long long,Point3D<Real> >();
		interiorPoints=new std::vector<Point3D<float> >();

		temp=sNodes.treeNodes[i]->nextLeaf();
		while(temp){
			if(temp->nodeData.isoNode){
				SetMCRootPositions(temp,sDepth,isoValue,temp->nodeData.isoNode->cornerValues,
					boundaryRoots,interiorRoots,*boundaryNormalHash,interiorNormalHash,interiorPoints,mesh,1);
				temp->nodeData.isoNode->mcIndex=MarchingCubes::GetIndex(temp->nodeData.isoNode->cornerValues,isoValue);
				temp->nodeData.isoNode->eSegmentCount=0;
				temp->nodeData.isoNode->eSegments[0]=temp->nodeData.isoNode->eSegments[1]=0;
			}
			temp=sNodes.treeNodes[i]->nextLeaf(temp);
		}
		delete interiorNormalHash;

		sLeaves.set(*sNodes.treeNodes[i]);
		for(int j=0;j<sLeaves.leafCount;j++){GetMCIsoTriangles(sLeaves.treeLeaves[j],mesh,boundaryRoots,interiorRoots,interiorPoints,offSet,sDepth);}
		delete interiorRoots;
		delete interiorPoints;
		offSet=mesh->outOfCorePointCount();
	}
	delete boundaryNormalHash;

	SortedTreeLeaves sLeaves;
	sLeaves.set(tree,sDepth-1);
	for(int j=0;j<sLeaves.leafCount;j++){GetMCIsoTriangles(sLeaves.treeLeaves[j],mesh,boundaryRoots,NULL,NULL,0,0);}
}
template<int Degree>
Real Octree<Degree>::GetIsoValue(void){
	const TreeOctNode* temp;
	Real isoValue,weightSum,w;
	Real myRadius;
	PointIndexValueFunction cf;
	Point3D<Real> center;
	Real width;

	fData.setValueTables(fData.VALUE_FLAG,0);
	cf.valueTables=fData.valueTables;
	cf.res2=fData.res2;
	myRadius=radius;
	isoValue=weightSum=0;
	temp=tree.nextNode();
	while(temp){
		w=temp->nodeData.centerWeightContribution;
		if(w>EPSILON){
			cf.value=0;
			VertexData::CenterIndex(temp,fData.depth,cf.index);
			temp->centerAndWidth(center,width);
			TreeOctNode::ProcessPointAdjacentNodes(center,&tree,myRadius,&cf);
			isoValue+=cf.value*w;
			weightSum+=w;
		}
		temp=tree.nextNode(temp);
	}
	return isoValue/weightSum;
}
template<int Degree>
void Octree<Degree>::SetIsoSurfaceCorners(const Real& isoValue,const int& subdivideDepth,const int& fullDepthIso){
	double t=Time();
	int i,j;
	hash_map<long long,Real> values;
	Real cornerValues[Cube::CORNERS];
	PointIndexValueFunction cf;
	Point3D<Real> position;
	TreeOctNode* temp;
	long long key;
	SortedTreeNodes *sNodes=new SortedTreeNodes();
	sNodes->set(tree,0);

	temp=tree.nextNode();
	while(temp){
		temp->nodeData.isoNode=NULL;
		temp=tree.nextNode(temp);
	}
	TreeNodeData::UseIndex=0;	
	// Start by setting the corner values of all the nodes
	cf.valueTables=fData.valueTables;
	cf.res2=fData.res2;
	for(i=0;i<sNodes->nodeCount[subdivideDepth];i++){
		temp=sNodes->treeNodes[i];
		if(!temp->children){
			for(j=0;j<Cube::CORNERS;j++){
				cf.value=0;
				VertexData::CornerIndex(temp,j,fData.depth,cf.index);
				position.coords[0]=BinaryNode<Real>::CornerIndexPosition(cf.index[0],fData.depth+1);
				position.coords[1]=BinaryNode<Real>::CornerIndexPosition(cf.index[1],fData.depth+1);
				position.coords[2]=BinaryNode<Real>::CornerIndexPosition(cf.index[2],fData.depth+1);
				TreeOctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);
				cornerValues[j]=cf.value;
			}
			if(temp->depth()<fData.depth || MarchingCubes::HasRoots(cornerValues,isoValue)){
				if(IsoNodeData::UseAllocator()){temp->nodeData.isoNode=IsoNodeData::Allocator.newElements(1);}
				else{temp->nodeData.isoNode=new IsoNodeData();}
				memcpy(temp->nodeData.isoNode->cornerValues,cornerValues,sizeof(Real)*Cube::CORNERS);
			}
		}
	}

	MemoryUsage();

	for(i=sNodes->nodeCount[subdivideDepth];i<sNodes->nodeCount[subdivideDepth+1];i++){
		temp=sNodes->treeNodes[i]->nextLeaf();
		while(temp){
			for(j=0;j<Cube::CORNERS;j++){
				key=VertexData::CornerIndex(temp,j,fData.depth,cf.index);
				if(values.find(key)!=values.end()){cornerValues[j]=values[key];}
				else{
					cf.value=0;
					position.coords[0]=BinaryNode<Real>::CornerIndexPosition(cf.index[0],fData.depth+1);
					position.coords[1]=BinaryNode<Real>::CornerIndexPosition(cf.index[1],fData.depth+1);
					position.coords[2]=BinaryNode<Real>::CornerIndexPosition(cf.index[2],fData.depth+1);
					TreeOctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);
					values[key]=cf.value;
					cornerValues[j]=cf.value;
				}
			}
			if(temp->depth()<fData.depth || MarchingCubes::HasRoots(cornerValues,isoValue)){
				if(IsoNodeData::UseAllocator()){temp->nodeData.isoNode=IsoNodeData::Allocator.newElements(1);}
				else{temp->nodeData.isoNode=new IsoNodeData();}
				memcpy(temp->nodeData.isoNode->cornerValues,cornerValues,sizeof(Real)*Cube::CORNERS);
			}
			temp=sNodes->treeNodes[i]->nextLeaf(temp);
		}
		MemoryUsage();
		values.clear();
	}
	delete sNodes;
	DumpOutput("Set corner values in: %f\n",Time()-t);
	DumpOutput("Memory Usage: %.3f MB\n",float(MemoryUsage()));

	if(subdivideDepth){
		printf("Pre-Validating\n");
		PreValidate(isoValue,fData.depth,subdivideDepth);
		printf("Done\n");
	}
	// Now validate all of the node edges
	t=Time();
	temp=tree.nextLeaf();
	while(temp){
		if(subdivideDepth)	{Validate(temp,isoValue,fData.depth,fullDepthIso,subdivideDepth);}
		else				{Validate(temp,isoValue,fData.depth,fullDepthIso);}
		temp=tree.nextLeaf(temp);
	}
	DumpOutput("Validated in: %f\n",Time()-t);
	DumpOutput("Memory Usage: %.3f MB\n",float(MemoryUsage()));
}
template<int Degree>
void Octree<Degree>::Subdivide(TreeOctNode* node,const Real& isoValue,const int& maxDepth){
	int i,j,c[4];
	Real value,cornerValues[Cube::CORNERS];
	Real cornerValues2[Cube::CORNERS][Cube::CORNERS];
	PointIndexValueFunction cf;
	Point3D<Real> position;
	cf.valueTables=fData.valueTables;
	cf.res2=fData.res2;

	for(i=0;i<Cube::CORNERS;i++){cornerValues[i]=node->nodeData.isoNode->cornerValues[i];}
	node->initChildren();
	// Since we are allocating blocks, it is possible that some of the memory was pre-allocated with
	// the wrong initialization
	for(i=0;i<Cube::CORNERS;i++){node->children[i].nodeData.isoNode=NULL;}

	// Now set the corner values for the new children
	// Copy old corner values
	for(i=0;i<Cube::CORNERS;i++){cornerValues2[i][i]=cornerValues[i];}
	// 8 of 27 corners set

	// Set center corner
	cf.value=0;
	VertexData::CenterIndex(node,maxDepth,cf.index);
	position.coords[0]=BinaryNode<Real>::CornerIndexPosition(cf.index[0],maxDepth+1);
	position.coords[1]=BinaryNode<Real>::CornerIndexPosition(cf.index[1],maxDepth+1);
	position.coords[2]=BinaryNode<Real>::CornerIndexPosition(cf.index[2],maxDepth+1);
	TreeOctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);
	value=cf.value;
	for(i=0;i<Cube::CORNERS;i++){cornerValues2[i][Cube::AntipodalCornerIndex(i)]=value;}
	// 9 of 27 set

	// Set face corners
	for(i=0;i<Cube::NEIGHBORS;i++){
		int dir,offset,e;
		Cube::FactorFaceIndex(i,dir,offset);
		cf.value=0;
		VertexData::FaceIndex(node,i,maxDepth,cf.index);
		position.coords[0]=BinaryNode<Real>::CornerIndexPosition(cf.index[0],maxDepth+1);
		position.coords[1]=BinaryNode<Real>::CornerIndexPosition(cf.index[1],maxDepth+1);
		position.coords[2]=BinaryNode<Real>::CornerIndexPosition(cf.index[2],maxDepth+1);
		TreeOctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);
		value=cf.value;
		Cube::FaceCorners(i,c[0],c[1],c[2],c[3]);
		e=Cube::EdgeIndex(dir,0,0);
		for(j=0;j<4;j++){cornerValues2[c[j]][Cube::EdgeReflectCornerIndex(c[j],e)]=value;}
	}
	// 15 of 27 set

	// Set edge corners
	for(i=0;i<Cube::EDGES;i++){
		int o,i1,i2,f;
		Cube::FactorEdgeIndex(i,o,i1,i2);
		cf.value=0;
		VertexData::EdgeIndex(node,i,maxDepth,cf.index);
		position.coords[0]=BinaryNode<Real>::CornerIndexPosition(cf.index[0],maxDepth+1);
		position.coords[1]=BinaryNode<Real>::CornerIndexPosition(cf.index[1],maxDepth+1);
		position.coords[2]=BinaryNode<Real>::CornerIndexPosition(cf.index[2],maxDepth+1);
		TreeOctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);
		value=cf.value;
		Cube::EdgeCorners(i,c[0],c[1]);
		f=Cube::FaceIndex(o,0);
		for(j=0;j<2;j++){cornerValues2[c[j]][Cube::FaceReflectCornerIndex(c[j],f)]=value;}
	}
	// 27 of 27 set

	for(i=0;i<Cube::CORNERS;i++){
		if(node->children[i].depth()<maxDepth || MarchingCubes::HasRoots(cornerValues2[i],isoValue)){
			if(IsoNodeData::UseAllocator()){node->children[i].nodeData.isoNode=IsoNodeData::Allocator.newElements(1);}
			else{node->children[i].nodeData.isoNode=new IsoNodeData();}
			memcpy(node->children[i].nodeData.isoNode->cornerValues,cornerValues2[i],sizeof(Real)*Cube::CORNERS);
		}
	}
	if(node->nodeData.isoNode){
		if(!IsoNodeData::UseAllocator()){delete node->nodeData.isoNode;}
		node->nodeData.isoNode=NULL;
	}
}

template<int Degree>
int Octree<Degree>::InteriorFaceRootCount(const TreeOctNode* node,const int &faceIndex,const Real& isoValue,const int& maxDepth){
	int c1,c2,e1,e2,dir,off,cnt=0;
	int corners[Cube::CORNERS/2];
	if(node->children){
		Cube::FaceCorners(faceIndex,corners[0],corners[1],corners[2],corners[3]);
		Cube::FactorFaceIndex(faceIndex,dir,off);
		c1=corners[0];
		c2=corners[3];
		switch(dir){
			case 0:
				e1=Cube::EdgeIndex(1,off,1);
				e2=Cube::EdgeIndex(2,off,1);
				break;
			case 1:
				e1=Cube::EdgeIndex(0,off,1);
				e2=Cube::EdgeIndex(2,1,off);
				break;
			case 2:
				e1=Cube::EdgeIndex(0,1,off);
				e2=Cube::EdgeIndex(1,1,off);
				break;
		};
		cnt+=EdgeRootCount(&node->children[c1],e1,isoValue,maxDepth)+EdgeRootCount(&node->children[c1],e2,isoValue,maxDepth);
		switch(dir){
			case 0:
				e1=Cube::EdgeIndex(1,off,0);
				e2=Cube::EdgeIndex(2,off,0);
				break;
			case 1:
				e1=Cube::EdgeIndex(0,off,0);
				e2=Cube::EdgeIndex(2,0,off);
				break;
			case 2:
				e1=Cube::EdgeIndex(0,0,off);
				e2=Cube::EdgeIndex(1,0,off);
				break;
		};
		cnt+=EdgeRootCount(&node->children[c2],e1,isoValue,maxDepth)+EdgeRootCount(&node->children[c2],e2,isoValue,maxDepth);
		for(int i=0;i<Cube::CORNERS/2;i++){if(node->children[corners[i]].children){cnt+=InteriorFaceRootCount(&node->children[corners[i]],faceIndex,isoValue,maxDepth);}}
	}
	return cnt;
}

template<int Degree>
int Octree<Degree>::EdgeRootCount(const TreeOctNode* node,const int& edgeIndex,const Real& isoValue,const int& maxDepth){
	int f1,f2,c1,c2;
	const TreeOctNode* temp;
	Cube::FacesAdjacentToEdge(edgeIndex,f1,f2);

	int eIndex;
	const TreeOctNode* finest=node;
	eIndex=edgeIndex;
	if(node->depth()<maxDepth){
		temp=node->faceNeighbor(f1);
		if(temp && temp->children){
			finest=temp;
			eIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f1);
		}
		else{
			temp=node->faceNeighbor(f2);
			if(temp && temp->children){
				finest=temp;
				eIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f2);
			}
			else{
				temp=node->edgeNeighbor(edgeIndex);
				if(temp && temp->children){
					finest=temp;
					eIndex=Cube::EdgeReflectEdgeIndex(edgeIndex);
				}
			}
		}
	}

	Cube::EdgeCorners(eIndex,c1,c2);
	if(finest->children){return EdgeRootCount(&finest->children[c1],eIndex,isoValue,maxDepth)+EdgeRootCount(&finest->children[c2],eIndex,isoValue,maxDepth);}
	else{
		if(!finest->nodeData.isoNode){return 0;}
		Real* cValues=finest->nodeData.isoNode->cornerValues;
		if((cValues[c1]<=isoValue && cValues[c2]<=isoValue) || (cValues[c1]>=isoValue && cValues[c2]>=isoValue)){return 0;}
		else{return 1;}
	}
}
template<int Degree>
int Octree<Degree>::IsBoundaryFace(const TreeOctNode* node,const int& faceIndex,const int& subdivideDepth){
	int dir,offset,d,o[3],idx;

	if(subdivideDepth<0){return 0;}
	if(node->d<=subdivideDepth){return 1;}
	Cube::FactorFaceIndex(faceIndex,dir,offset);
	node->depthAndOffset(d,o);

	idx=(int(o[dir])<<1) + (offset<<1);
	return !(idx%(2<<(int(node->d)-subdivideDepth)));
}
template<int Degree>
int Octree<Degree>::IsBoundaryEdge(const TreeOctNode* node,const int& edgeIndex,const int& subdivideDepth){
	int dir,x,y;
	Cube::FactorEdgeIndex(edgeIndex,dir,x,y);
	return IsBoundaryEdge(node,dir,x,y,subdivideDepth);
}
template<int Degree>
int Octree<Degree>::IsBoundaryEdge(const TreeOctNode* node,const int& dir,const int& x,const int& y,const int& subdivideDepth){
	int d,o[3],idx1,idx2,mask;

	if(subdivideDepth<0){return 0;}
	if(node->d<=subdivideDepth){return 1;}
	node->depthAndOffset(d,o);

	switch(dir){
		case 0:
			idx1=(int(o[1])<<1) + (x<<1);
			idx2=(int(o[2])<<1) + (y<<1);
			break;
		case 1:
			idx1=(int(o[0])<<1) + (x<<1);
			idx2=(int(o[2])<<1) + (y<<1);
			break;
		case 2:
			idx1=(int(o[0])<<1) + (x<<1);
			idx2=(int(o[1])<<1) + (y<<1);
			break;
		default:
			fprintf(stderr,"Bad Code\n");
			exit(0);
	}
	mask=2<<(int(node->d)-subdivideDepth);
	return !(idx1%(mask)) || !(idx2%(mask));
}

template<int Degree>
void Octree<Degree>::PreValidate(TreeOctNode* node,const Real& isoValue,const int& maxDepth,const int& subdivideDepth){
	int sub=0;
//	if(int(node->d)<subdivideDepth){sub=1;}
	for(int i=0;i<Cube::NEIGHBORS && !sub;i++){
		TreeOctNode* neighbor=node->faceNeighbor(i);
		if(neighbor && neighbor->children){
			if(IsBoundaryFace(node,i,subdivideDepth)){sub=1;}
		}
	}
	if(sub){
		Subdivide(node,isoValue,maxDepth);
		for(int i=0;i<Cube::NEIGHBORS;i++){
			if(IsBoundaryFace(node,i,subdivideDepth)){
				TreeOctNode* neighbor=node->faceNeighbor(i);
				while(neighbor && !neighbor->children){
					PreValidate(neighbor,isoValue,maxDepth,subdivideDepth);
					neighbor=node->faceNeighbor(i);
				}
			}
		}
	}
}

template<int Degree>
void Octree<Degree>::PreValidate(const Real& isoValue,const int& maxDepth,const int& subdivideDepth){
	TreeOctNode* temp;

	temp=tree.nextLeaf();
	while(temp){
		PreValidate(temp,isoValue,maxDepth,subdivideDepth);
		temp=tree.nextLeaf(temp);
	}

}
template<int Degree>
void Octree<Degree>::Validate(TreeOctNode* node,const Real& isoValue,const int& maxDepth,const int& fullDepthIso){
	int i,sub=0;
	TreeOctNode* treeNode=node;
	TreeOctNode* neighbor;
	if(node->depth()>=maxDepth){return;}

	// Check if full-depth extraction is enabled and we have an iso-node that is not at maximum depth
	if(!sub && fullDepthIso && MarchingCubes::HasRoots(node->nodeData.isoNode->cornerValues,isoValue)){sub=1;}

	// Check if the node has faces that are ambiguous and are adjacent to finer neighbors
	for(i=0;i<Cube::NEIGHBORS && !sub;i++){
		neighbor=treeNode->faceNeighbor(i);
		if(neighbor && neighbor->children){if(MarchingCubes::IsAmbiguous(node->nodeData.isoNode->cornerValues,isoValue,i)){sub=1;}}
	}

	// Check if the node has edges with more than one root
	for(i=0;i<Cube::EDGES && !sub;i++){if(EdgeRootCount(node,i,isoValue,maxDepth)>1){sub=1;}}

	for(i=0;i<Cube::NEIGHBORS && !sub;i++){
		neighbor=node->faceNeighbor(i);
		if(	neighbor && neighbor->children &&
			!MarchingCubes::HasRoots(node->nodeData.isoNode->cornerValues,isoValue,i) &&
			InteriorFaceRootCount(neighbor,Cube::FaceReflectFaceIndex(i,i),isoValue,maxDepth)){sub=1;}
	}
	if(sub){
		Subdivide(node,isoValue,maxDepth);
		for(i=0;i<Cube::NEIGHBORS;i++){
			neighbor=treeNode->faceNeighbor(i);
			if(neighbor && !neighbor->children){Validate(neighbor,isoValue,maxDepth,fullDepthIso);}
		}
		for(i=0;i<Cube::EDGES;i++){
			neighbor=treeNode->edgeNeighbor(i);
			if(neighbor && !neighbor->children){Validate(neighbor,isoValue,maxDepth,fullDepthIso);}
		}
		for(i=0;i<Cube::CORNERS;i++){if(!node->children[i].children){Validate(&node->children[i],isoValue,maxDepth,fullDepthIso);}}
	}
}
template<int Degree>
void Octree<Degree>::Validate(TreeOctNode* node,const Real& isoValue,const int& maxDepth,const int& fullDepthIso,const int& subdivideDepth){
	int i,sub=0;
	TreeOctNode* treeNode=node;
	TreeOctNode* neighbor;
	if(node->depth()>=maxDepth){return;}

	// Check if full-depth extraction is enabled and we have an iso-node that is not at maximum depth
	if(!sub && fullDepthIso && MarchingCubes::HasRoots(node->nodeData.isoNode->cornerValues,isoValue)){sub=1;}

	// Check if the node has faces that are ambiguous and are adjacent to finer neighbors
	for(i=0;i<Cube::NEIGHBORS && !sub;i++){
		neighbor=treeNode->faceNeighbor(i);
		if(neighbor && neighbor->children){if(MarchingCubes::IsAmbiguous(node->nodeData.isoNode->cornerValues,isoValue,i) || IsBoundaryFace(node,i,subdivideDepth)){sub=1;}}
	}

	// Check if the node has edges with more than one root
	for(i=0;i<Cube::EDGES && !sub;i++){if(EdgeRootCount(node,i,isoValue,maxDepth)>1){sub=1;}}

	for(i=0;i<Cube::NEIGHBORS && !sub;i++){
		neighbor=node->faceNeighbor(i);
		if(	neighbor && neighbor->children &&
			!MarchingCubes::HasRoots(node->nodeData.isoNode->cornerValues,isoValue,i) &&
			InteriorFaceRootCount(neighbor,Cube::FaceReflectFaceIndex(i,i),isoValue,maxDepth)){sub=1;}
	}
	if(sub){
		Subdivide(node,isoValue,maxDepth);
		for(i=0;i<Cube::NEIGHBORS;i++){
			neighbor=treeNode->faceNeighbor(i);
			if(neighbor && !neighbor->children){Validate(neighbor,isoValue,maxDepth,fullDepthIso,subdivideDepth);}
		}
		for(i=0;i<Cube::EDGES;i++){
			neighbor=treeNode->edgeNeighbor(i);
			if(neighbor && !neighbor->children){Validate(neighbor,isoValue,maxDepth,fullDepthIso,subdivideDepth);}
		}
		for(i=0;i<Cube::CORNERS;i++){if(!node->children[i].children){Validate(&node->children[i],isoValue,maxDepth,fullDepthIso,subdivideDepth);}}
	}
}
//////////////////////////////////////////////////////////////////////////////////////
// The assumption made when calling this code is that the edge has at most one root //
//////////////////////////////////////////////////////////////////////////////////////
template<int Degree>
int Octree<Degree>::GetRoot(const TreeOctNode* node,const int& edgeIndex,const Real* cValues,const Real& isoValue,Point3D<Real> & position,
							hash_map<long long,Point3D<Real> >& normalHash,const int& nonLinearFit){
	int c1,c2;
	Cube::EdgeCorners(edgeIndex,c1,c2);

	if((cValues[c1]<=isoValue && cValues[c2]<=isoValue) || (cValues[c1]>=isoValue && cValues[c2]>=isoValue)){return 0;}

	long long key;
	Point3D<Real> n[2];
	PointIndexNormalFunction cnf;
	Point3D<Real> p;
	cnf.valueTables=fData.valueTables;
	cnf.dValueTables=fData.dValueTables;
	cnf.res2=fData.res2;

	int i,o,i1,i2,rCount=0;
	Polynomial<2> P;
	std::vector<double> roots;
	double x0,x1;
	Real center,width;
	Real averageRoot=0;
	Cube::FactorEdgeIndex(edgeIndex,o,i1,i2);
	key=VertexData::CornerIndex(node,c1,fData.depth,cnf.index);

	if(normalHash.find(key)==normalHash.end()){
		cnf.normal.coords[0]=cnf.normal.coords[1]=cnf.normal.coords[2]=0;
		p.coords[0]=BinaryNode<Real>::CornerIndexPosition(cnf.index[0],fData.depth+1);
		p.coords[1]=BinaryNode<Real>::CornerIndexPosition(cnf.index[1],fData.depth+1);
		p.coords[2]=BinaryNode<Real>::CornerIndexPosition(cnf.index[2],fData.depth+1);
		TreeOctNode::ProcessPointAdjacentNodes(p,postNormalSmooth,&tree,radius,&cnf);
		normalHash[key]=cnf.normal;
	}
	n[0]=normalHash[key];

	key=VertexData::CornerIndex(node,c2,fData.depth,cnf.index);
	if(normalHash.find(key)==normalHash.end()){
		cnf.normal.coords[0]=cnf.normal.coords[1]=cnf.normal.coords[2]=0;
		p.coords[0]=BinaryNode<Real>::CornerIndexPosition(cnf.index[0],fData.depth+1);
		p.coords[1]=BinaryNode<Real>::CornerIndexPosition(cnf.index[1],fData.depth+1);
		p.coords[2]=BinaryNode<Real>::CornerIndexPosition(cnf.index[2],fData.depth+1);
		TreeOctNode::ProcessPointAdjacentNodes(p,postNormalSmooth,&tree,radius,&cnf);
		normalHash[key]=cnf.normal;
	}
	n[1]=normalHash[key];

	Point3D<Real> c;
	node->centerAndWidth(c,width);
	center=c.coords[o];
	for(i=0;i<DIMENSION;i++){
		n[0].coords[i]*=width;
		n[1].coords[i]*=width;
	}

	switch(o){
				case 0:
					position.coords[1]=c.coords[1]-width/2+width*i1;
					position.coords[2]=c.coords[2]-width/2+width*i2;
					break;
				case 1:
					position.coords[0]=c.coords[0]-width/2+width*i1;
					position.coords[2]=c.coords[2]-width/2+width*i2;
					break;
				case 2:
					position.coords[0]=c.coords[0]-width/2+width*i1;
					position.coords[1]=c.coords[1]-width/2+width*i2;
					break;
	}
	double dx0,dx1;
	x0=cValues[c1];
	x1=cValues[c2];
	dx0=n[0].coords[o];
	dx1=n[1].coords[o];

	// The scaling will turn the Hermite Spline into a quadratic
	double scl=(x1-x0)/((dx1+dx0)/2);
	dx0*=scl;
	dx1*=scl;

	// Hermite Spline
	P.coefficients[0]=x0;
	P.coefficients[1]=dx0;
	P.coefficients[2]=3*(x1-x0)-dx1-2*dx0;

	P.getSolutions(isoValue,roots,EPSILON);
	for(i=0;i<int(roots.size());i++){
		if(roots[i]>=0 && roots[i]<=1){
			averageRoot+=Real(roots[i]);
			rCount++;
		}
	}
	if(rCount && nonLinearFit)	{averageRoot/=rCount;}
	else						{averageRoot=Real((x0-isoValue)/(x0-x1));}

	position.coords[o]=Real(center-width/2+width*averageRoot);
	return 1;
}

template<int Degree>
int Octree<Degree>::GetRoot(const TreeOctNode* node,const int& edgeIndex,const Real& isoValue,const int& maxDepth,Point3D<Real>& position,hash_map<long long,Point3D<Real> >& normals,
							Point3D<Real>* normal,const int& nonLinearFit){
	int c1,c2,f1,f2;
	const TreeOctNode *temp,*finest;
	int eIndex;

	Cube::FacesAdjacentToEdge(edgeIndex,f1,f2);

	finest=node;
	eIndex=edgeIndex;
	if(node->depth()<fData.depth){
		temp=node->faceNeighbor(f1);
		if(temp && temp->children){
			finest=temp;
			eIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f1);
		}
		else{
			temp=node->faceNeighbor(f2);
			if(temp && temp->children){
				finest=temp;
				eIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f2);
			}
			else{
				temp=node->edgeNeighbor(edgeIndex);
				if(temp && temp->children){
					finest=temp;
					eIndex=Cube::EdgeReflectEdgeIndex(edgeIndex);
				}
			}
		}
	}

	Cube::EdgeCorners(eIndex,c1,c2);
	if(finest->children){
		if		(GetRoot(&finest->children[c1],eIndex,isoValue,maxDepth,position,normals,normal,nonLinearFit))	{return 1;}
		else if	(GetRoot(&finest->children[c2],eIndex,isoValue,maxDepth,position,normals,normal,nonLinearFit))	{return 1;}
		else																									{return 0;}
	}
	else{
		if(!finest->nodeData.isoNode){return 0;}
		return GetRoot(finest,eIndex,finest->nodeData.isoNode->cornerValues,isoValue,position,normals,nonLinearFit);
	}
}
template<int Degree>
int Octree<Degree>::GetRootIndex(const TreeOctNode* node,const int& edgeIndex,const Real& isoValue,const int& maxDepth,int eIndex[2],int& offset){
	int c1,c2,f1,f2;
	const TreeOctNode *temp,*finest;
	int finestIndex;

	Cube::FacesAdjacentToEdge(edgeIndex,f1,f2);

	finest=node;
	finestIndex=edgeIndex;
	if(node->depth()<maxDepth){
		temp=node->faceNeighbor(f1);
		if(temp && temp->children){
			finest=temp;
			finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f1);
		}
		else{
			temp=node->faceNeighbor(f2);
			if(temp && temp->children){
				finest=temp;
				finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f2);
			}
			else{
				temp=node->edgeNeighbor(edgeIndex);
				if(temp && temp->children){
					finest=temp;
					finestIndex=Cube::EdgeReflectEdgeIndex(edgeIndex);
				}
			}
		}
	}

	Cube::EdgeCorners(finestIndex,c1,c2);
	if(finest->children){
		if		(GetRootIndex(&finest->children[c1],finestIndex,isoValue,maxDepth,eIndex,offset))	{return 1;}
		else if	(GetRootIndex(&finest->children[c2],finestIndex,isoValue,maxDepth,eIndex,offset))	{return 1;}
		else																							{return 0;}
	}
	else{
		if(!finest->nodeData.isoNode){return 0;}
		Real* cValues=finest->nodeData.isoNode->cornerValues;
		if((cValues[c1]<=isoValue && cValues[c2]<=isoValue) || (cValues[c1]>=isoValue && cValues[c2]>=isoValue)){return 0;}

		int o,i1,i2;
		Cube::FactorEdgeIndex(finestIndex,o,i1,i2);
		int d,off[3];
		finest->depthAndOffset(d,off);
		offset=BinaryNode<Real>::Index(d,off[o]);
		switch(o){
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
		return 1;
	}
}
template<int Degree>
int Octree<Degree>::GetRootIndex(const TreeOctNode* node,const int& edgeIndex,const int& maxDepth,const int& sDepth,int eIndex[2],int& offset){
	int c1,c2,f1,f2;
	const TreeOctNode *temp,*finest;
	int finestIndex;

	Cube::FacesAdjacentToEdge(edgeIndex,f1,f2);

	finest=node;
	finestIndex=edgeIndex;
	if(node->depth()<maxDepth){
		if(IsBoundaryFace(node,f1,sDepth)){temp=NULL;}
		else{temp=node->faceNeighbor(f1);}
		if(temp && temp->children){
			finest=temp;
			finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f1);
		}
		else{
			if(IsBoundaryFace(node,f2,sDepth)){temp=NULL;}
			else{temp=node->faceNeighbor(f2);}
			if(temp && temp->children){
				finest=temp;
				finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f2);
			}
			else{
				if(IsBoundaryEdge(node,edgeIndex,sDepth)){temp=NULL;}
				else{temp=node->edgeNeighbor(edgeIndex);}
				if(temp && temp->children){
					finest=temp;
					finestIndex=Cube::EdgeReflectEdgeIndex(edgeIndex);
				}
			}
		}
	}

	Cube::EdgeCorners(finestIndex,c1,c2);
	if(finest->children){
		if		(GetRootIndex(&finest->children[c1],finestIndex,maxDepth,sDepth,eIndex,offset))	{return 1;}
		else if	(GetRootIndex(&finest->children[c2],finestIndex,maxDepth,sDepth,eIndex,offset))	{return 1;}
		else																							{return 0;}
	}
	else{
		if(!finest->nodeData.isoNode){return 0;}
		if(!(MarchingCubes::edgeMask[finest->nodeData.isoNode->mcIndex] & (1<<finestIndex))){return 0;}

		int o,i1,i2;
		Cube::FactorEdgeIndex(finestIndex,o,i1,i2);
		int d,off[3];
		finest->depthAndOffset(d,off);
		offset=BinaryNode<Real>::Index(d,off[o]);
		switch(o){
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
		return 1;
	}
}

template<int Degree>
int Octree<Degree>::GetRootIndex(const TreeOctNode* node,const int& edgeIndex,const int& maxDepth,int eIndex[2],int& offset){
	int c1,c2,f1,f2;
	const TreeOctNode *temp,*finest;
	int finestIndex;


	// The assumption is that the super-edge has a root along it. 
	if(!node->nodeData.isoNode || !(MarchingCubes::edgeMask[node->nodeData.isoNode->mcIndex] & (1<<edgeIndex))){return 0;}

	Cube::FacesAdjacentToEdge(edgeIndex,f1,f2);

	finest=node;
	finestIndex=edgeIndex;
	if(node->depth()<maxDepth){
		temp=node->faceNeighbor(f1);
		if(temp && temp->children){
			finest=temp;
			finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f1);
		}
		else{
			temp=node->faceNeighbor(f2);
			if(temp && temp->children){
				finest=temp;
				finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f2);
			}
			else{
				temp=node->edgeNeighbor(edgeIndex);
				if(temp && temp->children){
					finest=temp;
					finestIndex=Cube::EdgeReflectEdgeIndex(edgeIndex);
				}
			}
		}
	}

	Cube::EdgeCorners(finestIndex,c1,c2);
	if(finest->children){
		if		(GetRootIndex(&finest->children[c1],finestIndex,maxDepth,eIndex,offset))	{return 1;}
		else if	(GetRootIndex(&finest->children[c2],finestIndex,maxDepth,eIndex,offset))	{return 1;}
		else																							{return 0;}
	}
	else{
//		if(!finest->nodeData.isoNode || (MarchingCubes::edgeMask[finest->nodeData.isoNode->mcIndex] & (1<<finestIndex))){return 0;}

		int o,i1,i2;
		Cube::FactorEdgeIndex(finestIndex,o,i1,i2);
		int d,off[3];
		finest->depthAndOffset(d,off);
		offset=BinaryNode<Real>::Index(d,off[o]);
		switch(o){
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
		return 1;
	}
}
template<int Degree>
long long Octree<Degree>::GetRootIndex(const TreeOctNode* node,const int& edgeIndex,const int& maxDepth){
	int o,i1,i2;
	Cube::FactorEdgeIndex(edgeIndex,o,i1,i2);
	int d,off[3],eIndex[2],offset;
	node->depthAndOffset(d,off);
	offset=BinaryNode<Real>::Index(d,off[o]);
	switch(o){
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
	return (long long)(o) | (long long)(eIndex[0])<<5 | (long long)(eIndex[1])<<25 | (long long)(offset)<<45;
}
template<int Degree>
int Octree<Degree>::GetRootIndices(const TreeOctNode* node,const int& edgeIndex,const int& maxDepth,std::vector<long long>& indices){
	int c1,c2,f1,f2;
	const TreeOctNode *temp,*finest;
	int finestIndex;

	Cube::FacesAdjacentToEdge(edgeIndex,f1,f2);

	finest=node;
	finestIndex=edgeIndex;
	if(node->depth()<maxDepth){
		temp=node->faceNeighbor(f1);
		if(temp && temp->children){
			finest=temp;
			finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f1);
		}
		else{
			temp=node->faceNeighbor(f2);
			if(temp && temp->children){
				finest=temp;
				finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f2);
			}
			else{
				temp=node->edgeNeighbor(edgeIndex);
				if(temp && temp->children){
					finest=temp;
					finestIndex=Cube::EdgeReflectEdgeIndex(edgeIndex);
				}
			}
		}
	}

	Cube::EdgeCorners(finestIndex,c1,c2);
	if(finest->children){return GetRootIndices(&finest->children[c1],finestIndex,maxDepth,indices)+GetRootIndices(&finest->children[c2],finestIndex,maxDepth,indices);}
	else{
		if(!finest->nodeData.isoNode || !(MarchingCubes::edgeMask[finest->nodeData.isoNode->mcIndex] & (1<<finestIndex))){return 0;}
		else{
			indices.push_back(GetRootIndex(finest,finestIndex,maxDepth));
			return 1;
		}
	}
}

template<int Degree>
int Octree<Degree>::GetRootIndex(const TreeOctNode* node,const int& localEdgeIndex,const int& maxDepth,
								 hash_map<long long,int>& boundaryRoots,hash_map<long long,int>* interiorRoots,CoredPointIndex& index){
	int edgeIndex[2],offset,o,i1,i2;
	if(GetRootIndex(node,localEdgeIndex,maxDepth,edgeIndex,offset)){
		Cube::FactorEdgeIndex(localEdgeIndex,o,i1,i2);
		long long key= (long long)(o) | (long long)(edgeIndex[0])<<5 | (long long)(edgeIndex[1])<<25 | (long long)(offset)<<45;
		hash_map<long long,int>::iterator rootIter=boundaryRoots.find(key);
		if(rootIter!=boundaryRoots.end()){
			index.inCore=1;
			index.index=rootIter->second;
			return 1;
		}
		else if(interiorRoots){
			rootIter=interiorRoots->find(key);
			if(rootIter!=interiorRoots->end()){
				index.inCore=0;
				index.index=rootIter->second;
				return 1;
			}
		}
	}
	return 0;
}
template<int Degree>
int Octree<Degree>::SetMCRootPositions(const Real& isoValue,hash_map<long long,int>& roots,hash_map<long long,Point3D<Real> >& normalHash,
									   std::vector<Point3D<Real> >& positions,std::vector<Point3D<Real> >* normals,const int& nonLinearFit){
	TreeOctNode* temp;
	Point3D<Real> position,normal;
	int i,j,k,eIndex;
	int edgeIndex[2],offset;

	temp=tree.nextLeaf();
	while(temp){
		if(!temp->nodeData.isoNode){
			temp=tree.nextLeaf(temp);
			continue;
		}
		for(i=0;i<DIMENSION;i++){
			for(j=0;j<2;j++){
				for(k=0;k<2;k++){
					long long key;
					int rData;
					eIndex=Cube::EdgeIndex(i,j,k);
					if(GetRootIndex(temp,eIndex,isoValue,fData.depth,edgeIndex,offset)){
						key= (long long)(i) | (long long)(edgeIndex[0])<<5 | (long long)(edgeIndex[1])<<25 | (long long)(offset)<<45;
						if(roots.find(key)==roots.end()){
							if(normals)	{GetRoot(temp,eIndex,isoValue,fData.depth,position,normalHash,&normal,nonLinearFit);}
							else		{GetRoot(temp,eIndex,isoValue,fData.depth,position,normalHash,NULL,nonLinearFit);}
							rData=int(positions.size());
							positions.push_back(position);
							if(normals){normals->push_back(normal);}
							roots[key]=rData;
						}
					}
				}
			}
		}
		temp=tree.nextLeaf(temp);
	}
	return int(positions.size());
}
template<int Degree>
int Octree<Degree>::SetMCRootPositions(TreeOctNode* node,const int& sDepth,const Real& isoValue,const Real* cornerValues,
									   hash_map<long long,int>& boundaryRoots,hash_map<long long,int>* interiorRoots,
									   hash_map<long long,Point3D<Real> >& boundaryNormalHash,hash_map<long long,Point3D<Real> >* interiorNormalHash,
									   std::vector<Point3D<float> >* interiorPositions,
									   CoredMeshData* mesh,const int& nonLinearFit){
	Point3D<Real> position;
	int i,j,k,eIndex,f1,f2;
	int count=0;
	TreeOctNode* temp;

	for(i=0;i<DIMENSION;i++){
		for(j=0;j<2;j++){
			for(k=0;k<2;k++){
				long long key;
				eIndex=Cube::EdgeIndex(i,j,k);
				Cube::FacesAdjacentToEdge(eIndex,f1,f2);

				temp=node->faceNeighbor(f1);
				if(temp && temp->children){continue;}
				temp=node->faceNeighbor(f2);
				if(temp && temp->children){continue;}
				temp=node->edgeNeighbor(eIndex);
				if(temp && temp->children){continue;}


				key=GetRootIndex(node,eIndex,fData.depth);
				if(!interiorRoots || IsBoundaryEdge(node,i,j,k,sDepth)){
					if(boundaryRoots.find(key)==boundaryRoots.end()){
						if(GetRoot(node,eIndex,cornerValues,isoValue,position,boundaryNormalHash,nonLinearFit)){
							mesh->inCorePoints.push_back(position);
							boundaryRoots[key]=int(mesh->inCorePoints.size())-1;
							count++;
						}
					}
				}
				else{
					if(interiorRoots->find(key)==interiorRoots->end()){
						if(GetRoot(node,eIndex,cornerValues,isoValue,position,*interiorNormalHash,nonLinearFit)){
							(*interiorRoots)[key]=mesh->addOutOfCorePoint(position);
							interiorPositions->push_back(position);
							count++;
						}
					}
				}
			}
		}
	}
	return count;
}

template<int Degree>
int Octree<Degree>::SetMCRootPositions(TreeOctNode* node,const int& sDepth,const Real& isoValue,
									   hash_map<long long,int>& boundaryRoots,hash_map<long long,int>* interiorRoots,
									   hash_map<long long,Point3D<Real> >& boundaryNormalHash,hash_map<long long,Point3D<Real> >* interiorNormalHash,
									   std::vector<Point3D<float> >* interiorPositions,
									   CoredMeshData* mesh,const int& nonLinearFit){
	Point3D<Real> position;
	int i,j,k,eIndex;
	int edgeIndex[2],offset;
	int count=0;

	if(!node->nodeData.isoNode){return 0;}
	for(i=0;i<DIMENSION;i++){
		for(j=0;j<2;j++){
			for(k=0;k<2;k++){
				long long key;
				eIndex=Cube::EdgeIndex(i,j,k);
				if(GetRootIndex(node,eIndex,isoValue,fData.depth,edgeIndex,offset)){
					key= (long long)(i) | (long long)(edgeIndex[0])<<5 | (long long)(edgeIndex[1])<<25 | (long long)(offset)<<45;
					if(!interiorRoots || IsBoundaryEdge(node,i,j,k,sDepth)){
						if(boundaryRoots.find(key)==boundaryRoots.end()){
							GetRoot(node,eIndex,isoValue,fData.depth,position,boundaryNormalHash,NULL,nonLinearFit);
							mesh->inCorePoints.push_back(position);
							boundaryRoots[key]=int(mesh->inCorePoints.size())-1;
							count++;
						}
					}
					else{
						if(interiorRoots->find(key)==interiorRoots->end()){
							GetRoot(node,eIndex,isoValue,fData.depth,position,*interiorNormalHash,NULL,nonLinearFit);
							(*interiorRoots)[key]=mesh->addOutOfCorePoint(position);
							interiorPositions->push_back(position);
							count++;
						}
					}
				}
			}
		}
	}
	return count;
}
template<int Degree>
int Octree<Degree>::SetBoundaryMCRootPositions(const int& sDepth,const Real& isoValue,
											   hash_map<long long,int>& boundaryRoots,hash_map<long long,Point3D<Real> >& boundaryNormalHash,
											   CoredMeshData* mesh,const int& nonLinearFit){
	Point3D<Real> position;
	int i,j,k,eIndex,hits;
	int edgeIndex[2],offset;
	int count=0;
	TreeOctNode* node;

	node=tree.nextLeaf();
	while(node){
		hits=0;
		if(node->nodeData.isoNode){
			for(i=0;i<DIMENSION;i++){
				for(j=0;j<2;j++){
					for(k=0;k<2;k++){
						if(IsBoundaryEdge(node,i,j,k,sDepth)){
							hits++;
							long long key;
							eIndex=Cube::EdgeIndex(i,j,k);
							if(GetRootIndex(node,eIndex,isoValue,fData.depth,edgeIndex,offset)){
								key= (long long)(i) | (long long)(edgeIndex[0])<<5 | (long long)(edgeIndex[1])<<25 | (long long)(offset)<<45;
								if(boundaryRoots.find(key)==boundaryRoots.end()){
									GetRoot(node,eIndex,isoValue,fData.depth,position,boundaryNormalHash,NULL,nonLinearFit);
									mesh->inCorePoints.push_back(position);
									boundaryRoots[key]=int(mesh->inCorePoints.size())-1;
									count++;
								}
							}
						}
					}
				}
			}
		}
		if(hits){node=tree.nextLeaf(node);}
		else{node=tree.nextBranch(node);}
	}
	return count;
}
template<int Degree>
int Octree<Degree>::GetMCIsoTriangles(TreeOctNode* node,CoredMeshData* mesh,hash_map<long long,int>& boundaryRoots,
									  hash_map<long long,int>* interiorRoots,std::vector<Point3D<float> >* interiorPositions,const int& offSet,const int& sDepth)
{
	int isoTri[DIMENSION*MarchingCubes::MAX_TRIANGLES];
	size_t i,j,count,cnt=0;
	std::vector<CoredPointIndex> curves[3];
	if(node->nodeData.isoNode){count=MarchingCubes::AddTriangleIndices(node->nodeData.isoNode->mcIndex,isoTri);}
	else{count=0;}
	if(!count){return 0;}
	for(i=0;i<count;i++){
		for(j=0;j<3;j++){SetMCFaceCurve(node,isoTri[3*i+j],isoTri[3*i+(j+1)%3],fData.depth,curves[j],mesh->inCorePoints,boundaryRoots,interiorPositions,offSet,interiorRoots,sDepth);}
		cnt+=AddTriangles(mesh,curves,interiorPositions,offSet);
	}
	return int(cnt);
}

template<int Degree>
int Octree<Degree>::SetMCFaceCurve(TreeOctNode* node,const int& edgeIndex1,const int& edgeIndex2,const int& maxDepth,std::vector<CoredPointIndex>& curve,
								   std::vector<Point3D<float> >& boundaryPositions,hash_map<long long,int>& boundaryRoots,
								   std::vector<Point3D<float> >* interiorPositions,const int& offSet,hash_map<long long,int>* interiorRoots,const int& sDepth){
	TreeOctNode *neighbor,*subNeighbor;
	curve.clear();
	int faceIndex=node->nodeData.isoNode->addEdgeSegment(edgeIndex1,edgeIndex2);
	CoredPointIndex pIndex;
	long long key;
	int eCount;
	int o,i1,i2;
	int edgeIndex[2],offset;
	hash_map<long long,int>::iterator rootIter;
	if(!GetRootIndex(node,edgeIndex1,maxDepth,sDepth,edgeIndex,offset)){
//		fprintf(stderr,"Bad Edge Node 1!!!\n");
		return 0;
	}
	Cube::FactorEdgeIndex(edgeIndex1,o,i1,i2);
	key= (long long)(o) | (long long)(edgeIndex[0])<<5 | (long long)(edgeIndex[1])<<25 | (long long)(offset)<<45;
	rootIter=boundaryRoots.find(key);
	if(rootIter==boundaryRoots.end()){
		if(!interiorRoots){
//			fprintf(stderr,"Bad Edge Node 2a\n");
			return 0;
		}
		else{
			rootIter=interiorRoots->find(key);
			if(rootIter==interiorRoots->end()){
//				fprintf(stderr,"Bad Edge Node 2b\n");
				return 0;
			}
			else{pIndex.inCore=0;}
		}
	}
	else{pIndex.inCore=1;}
	pIndex.index=rootIter->second;
	neighbor=node->faceNeighbor(faceIndex);
	if(faceIndex<0 || !neighbor || !neighbor->children || IsBoundaryFace(node,faceIndex,sDepth)){curve.push_back(pIndex);}
	else{
		int eIndex=Cube::FaceReflectEdgeIndex(edgeIndex1,faceIndex);
		int fIndex=Cube::FaceReflectFaceIndex(faceIndex,faceIndex);
		TreeOctNode* subNode;
		if(pIndex.inCore)	{subNode=neighbor->getNearestLeaf(boundaryPositions[pIndex.index]);}
		else				{subNode=neighbor->getNearestLeaf((*interiorPositions)[pIndex.index-offSet]);}
		while(1){
			if(!subNode){
//				fprintf(stderr,"Bad SubNode\n");
				return 0;
			}
			else if(!subNode->nodeData.isoNode){
//				fprintf(stderr,"Bad SubNode->NodeData.IsoNode\n");
				return 0;
			}
			curve.push_back(pIndex);
			eCount=subNode->nodeData.isoNode->edgeCount(fIndex);
			for(int i=0;i<eCount;i++){
				for(int j=0;j<2;j++){if(subNode->nodeData.isoNode->edgeIndex(fIndex,i,j)==eIndex){eIndex=subNode->nodeData.isoNode->edgeIndex(fIndex,i,(j+1)%2);}}
			}
			if(TreeOctNode::CommonEdge(node,edgeIndex2,subNode,eIndex)){break;}
			if(!GetRootIndex(subNode,eIndex,maxDepth,edgeIndex,offset)){
//				fprintf(stderr,"Bad Edge Node 3\n");
				return 0;
			}
			Cube::FactorEdgeIndex(eIndex,o,i1,i2);
			key= (long long)(o) | (long long)(edgeIndex[0])<<5 | (long long)(edgeIndex[1])<<25 | (long long)(offset)<<45;
			rootIter=boundaryRoots.find(key);
			if(rootIter==boundaryRoots.end()){
				if(!interiorRoots){
//					fprintf(stderr,"Bad Edge Node 3\n");
					return 0;
				}
				else{
					rootIter=interiorRoots->find(key);
					if(rootIter==interiorRoots->end()){
//						fprintf(stderr,"Bad Edge Node 3\n");
						return 0;
					}
					else{pIndex.inCore=0;}
				}
			}
			else{pIndex.inCore=1;}
			pIndex.index=rootIter->second;

			int f,f1,f2;
			Cube::FacesAdjacentToEdge(eIndex,f1,f2);
			if(f1==fIndex)	{f=f2;}
			else			{f=f1;}
			subNeighbor=subNode->faceNeighbor(f);
			if(pIndex.inCore)	{subNode=subNeighbor->getNearestLeaf(boundaryPositions[pIndex.index]);}
			else				{subNode=subNeighbor->getNearestLeaf((*interiorPositions)[pIndex.index-offSet]);}
			eIndex=Cube::FaceReflectEdgeIndex(eIndex,f);
		}
	}
	return 1;
}
template<int Degree>
int Octree<Degree>::AddTriangles(CoredMeshData* mesh,std::vector<CoredPointIndex> edges[3],std::vector<Point3D<float> >* interiorPositions,const int& offSet){
	std::vector<CoredPointIndex> e;
	for(int i=0;i<3;i++){for(size_t j=0;j<edges[i].size();j++){e.push_back(edges[i][j]);}}
	return AddTriangles(mesh,e,interiorPositions,offSet);
}
template<int Degree>
int Octree<Degree>::AddTriangles(CoredMeshData* mesh,std::vector<int>& edges){
	int i,j,k,startIndex,endIndex;
	double ar;
	Triangle t;
	int pCount=int(edges.size());

	if(pCount>3){
		std::vector<int> start,end;
		for(i=0;i<pCount;i++){
			for(j=0;j<3;j++){for(k=0;k<3;k++){t.p[j][k]=mesh->inCorePoints[edges[(i+j)%pCount]].coords[k];}}
			double temp=t.AspectRatio();
			if(!i || temp<ar){
				ar=temp;
				startIndex=i;
			}
		}
		startIndex=(startIndex+1)%pCount;
		endIndex=startIndex+pCount/2;
		for(j=startIndex;j<=endIndex;j++){start.push_back(edges[j%pCount]);}
		for(j=endIndex;j<=startIndex+pCount;j++){end.push_back(edges[j%pCount]);}
		AddTriangles(mesh,start);
		AddTriangles(mesh,end);
	}
	if(pCount==3){
		TriangleIndex tri;
		for(j=0;j<3;j++){tri.idx[j]=edges[j];}
		mesh->addTriangle(tri);
	}
	return int(edges.size())-2;
}
template<int Degree>
int Octree<Degree>::AddTriangles(CoredMeshData* mesh,std::vector<CoredPointIndex>& edges,std::vector<Point3D<float> >* interiorPositions,const int& offSet){
	int i,j,k,startIndex,endIndex;
	double ar;
	Triangle t;
	int pCount=int(edges.size());

	if(pCount>3){
		std::vector<CoredPointIndex> start,end;
		startIndex=0;
		ar=0;
		for(i=0;i<pCount;i++){
			for(j=0;j<3;j++){
				if(edges[(i+j)%pCount].inCore)	{for(k=0;k<3;k++){t.p[j][k]=mesh->inCorePoints[edges[(i+j)%pCount].index].coords[k];}}
				else							{for(k=0;k<3;k++){t.p[j][k]=(*interiorPositions)[edges[(i+j)%pCount].index-offSet].coords[k];}}
			}
			double temp=t.AspectRatio();
			if(!i || temp<ar){
				ar=temp;
				startIndex=i;
			}
		}
		startIndex=(startIndex+1)%pCount;
		endIndex=startIndex+pCount/2;
		for(j=startIndex;j<=endIndex;j++){start.push_back(edges[j%pCount]);}
		for(j=endIndex;j<=startIndex+pCount;j++){end.push_back(edges[j%pCount]);}
		AddTriangles(mesh,start,interiorPositions,offSet);
		AddTriangles(mesh,end,interiorPositions,offSet);
	}
	if(pCount==3){
		TriangleIndex tri;
		int inCoreFlag=0;
		for(j=0;j<3;j++){
			tri.idx[j]=edges[j].index;
			if(edges[j].inCore){inCoreFlag|=CoredMeshData::IN_CORE_FLAG[j];}
		}
		mesh->addTriangle(tri,inCoreFlag);
	}
	return int(edges.size())-2;
}
////////////////
// VertexData //
////////////////
long long VertexData::CenterIndex(const TreeOctNode* node,const int& maxDepth){
	int idx[DIMENSION];
	return CenterIndex(node,maxDepth,idx);
}
long long VertexData::CenterIndex(const TreeOctNode* node,const int& maxDepth,int idx[DIMENSION]){
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);}
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
long long VertexData::CenterIndex(const int& depth,const int offSet[DIMENSION],const int& maxDepth,int idx[DIMENSION]){
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,depth+1,offSet[i]<<1,1);}
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
long long VertexData::CornerIndex(const TreeOctNode* node,const int& cIndex,const int& maxDepth){
	int idx[DIMENSION];
	return CornerIndex(node,cIndex,maxDepth,idx);
}
long long VertexData::CornerIndex(const TreeOctNode* node,const int& cIndex,const int& maxDepth,int idx[DIMENSION]){
	int x[DIMENSION];
	Cube::FactorCornerIndex(cIndex,x[0],x[1],x[2]);
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,o[i],x[i]);}
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
long long VertexData::CornerIndex(const int& depth,const int offSet[DIMENSION],const int& cIndex,const int& maxDepth,int idx[DIMENSION]){
	int x[DIMENSION];
	Cube::FactorCornerIndex(cIndex,x[0],x[1],x[2]);
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,depth,offSet[i],x[i]);}
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
long long VertexData::FaceIndex(const TreeOctNode* node,const int& fIndex,const int& maxDepth){
	int idx[DIMENSION];
	return FaceIndex(node,fIndex,maxDepth,idx);
}
long long VertexData::FaceIndex(const TreeOctNode* node,const int& fIndex,const int& maxDepth,int idx[DIMENSION]){
	int dir,offset;
	Cube::FactorFaceIndex(fIndex,dir,offset);
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);}
	idx[dir]=BinaryNode<Real>::CornerIndex(maxDepth+1,d,o[dir],offset);
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
long long VertexData::EdgeIndex(const TreeOctNode* node,const int& eIndex,const int& maxDepth){
	int idx[DIMENSION];
	return FaceIndex(node,eIndex,maxDepth,idx);
}
long long VertexData::EdgeIndex(const TreeOctNode* node,const int& eIndex,const int& maxDepth,int idx[DIMENSION]){
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
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
