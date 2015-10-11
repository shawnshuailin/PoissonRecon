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

#include <stdlib.h>
#include <math.h>
#include <algorithm>

/////////////
// OctNode //
/////////////
template< class NodeData > const int OctNode< NodeData >::DepthShift=5;
template< class NodeData > const int OctNode< NodeData >::OffsetShift = ( sizeof(long long)*8 - DepthShift ) / 3;
template< class NodeData > const int OctNode< NodeData >::DepthMask=(1<<DepthShift)-1;
template< class NodeData > const int OctNode< NodeData >::OffsetMask=(1<<OffsetShift)-1;
template< class NodeData > const int OctNode< NodeData >::OffsetShift1=DepthShift;//低位是depth，从depthshift到depthshift+offsetshift是第一个offset
template< class NodeData > const int OctNode< NodeData >::OffsetShift2=OffsetShift1+OffsetShift;//这是第二个offset
template< class NodeData > const int OctNode< NodeData >::OffsetShift3=OffsetShift2+OffsetShift;//这是第三个offset

template< class NodeData > int OctNode< NodeData >::UseAlloc=0;
template< class NodeData > Allocator<OctNode< NodeData > > OctNode< NodeData >::NodeAllocator;

template< class NodeData >
void OctNode< NodeData >::SetAllocator(int blockSize)
{
	if(blockSize>0)
	{
		UseAlloc=1;
		NodeAllocator.set(blockSize);
	}
	else{UseAlloc=0;}
}
template< class NodeData >
int OctNode< NodeData >::UseAllocator(void){return UseAlloc;}

template< class NodeData >
OctNode< NodeData >::OctNode(void){
	parent=children=NULL;
	_depthAndOffset = 0;
}

template< class NodeData >
OctNode< NodeData >::~OctNode(void){
	if(!UseAlloc){if(children){delete[] children;}}
	parent=children=NULL;
}
template< class NodeData >
void OctNode< NodeData >::setFullDepth( int maxDepth )
{
	//maxDepth代表了完全八叉树的深度，因此这里要对所有maxDepth之上的child node进行初始化，迭代进行
	if( maxDepth )
	{
		if( !children ) initChildren();
		for( int i=0 ; i<8 ; i++ ) children[i].setFullDepth( maxDepth-1 );//迭代设定child node
	}
}

template< class NodeData >
int OctNode< NodeData >::initChildren( void )
{
	if( UseAlloc ) children=NodeAllocator.newElements(8);
	else
	{
		if( children ) delete[] children;
		children = NULL;
		children = new OctNode[Cube::CORNERS];
	}
	if( !children )
	{
		fprintf(stderr,"Failed to initialize children in OctNode::initChildren\n");
		exit(0);
		return 0;
	}
	int d , off[3];
	depthAndOffset( d , off );//当前node的深度和offset
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ )
	{
		int idx=Cube::CornerIndex(i,j,k);
		children[idx].parent = this;
		children[idx].children = NULL;
		int off2[3];
		off2[0] = (off[0]<<1)+i;//感觉这个off所代表的偏移记录了node本身以及其所有parent节点在每一次octree划分时所处的位置
		off2[1] = (off[1]<<1)+j;//根据这个off可以迅速定位在每一层划分时，这个节点当时的父节点所处的划分位置
		off2[2] = (off[2]<<1)+k;//这套offset的index系统因为在不同depth下会产生重复，因此必须与depth相配合，才能找出正确的node
		//如果size(long long)等于16，那么offset中每一个index的长度是41bit，足够作为每一层的index
		children[idx]._depthAndOffset = Index( d+1 , off2 );
	}
	return 1;
}
template< class NodeData >
inline void OctNode< NodeData >::Index(int depth,const int offset[3],short& d,short off[3]){
	d=short(depth);
	off[0]=short((1<<depth)+offset[0]-1);//这个是绝对索引, 1<<depth是上面所有层的node总数，本层node的index初始值就是1<<depth-1
	off[1]=short((1<<depth)+offset[1]-1);
	off[2]=short((1<<depth)+offset[2]-1);
}

template< class NodeData >
inline void OctNode< NodeData >::depthAndOffset( int& depth , int offset[DIMENSION] ) const
{
	//这个变量在节点声明时就存储了depth和offset，现在只是进行逆操作把对应值取出来
	depth = int( _depthAndOffset & DepthMask );
	offset[0] = int( (_depthAndOffset>>OffsetShift1) & OffsetMask );
	offset[1] = int( (_depthAndOffset>>OffsetShift2) & OffsetMask );
	offset[2] = int( (_depthAndOffset>>OffsetShift3) & OffsetMask );
}
template< class NodeData >
inline void OctNode< NodeData >::centerIndex( int index[DIMENSION] ) const
{
	//根据深度和offset计算出当前节点在binaryNode情况下的center index
	int d , off[DIMENSION];
	depthAndOffset( d , off );
	for( int i=0 ; i<DIMENSION ; i++ ) index[i] = BinaryNode::CenterIndex( d , off[i] );
}
template< class NodeData >
inline unsigned long long OctNode< NodeData >::Index( int depth , const int offset[3] )
{
	//这个函数是把node的depth和offset经偏移操作后放在变量_depthAndOffset中
	//这些偏移的位数是提前约定好的
	//返回值应该与_depthAndOffset相等，offset属于相对偏移
	unsigned long long idx=0;
	idx |= ( ( (unsigned long long)(depth    ) ) & DepthMask  );
	idx |= ( ( (unsigned long long)(offset[0]) ) & OffsetMask ) << OffsetShift1;
	idx |= ( ( (unsigned long long)(offset[1]) ) & OffsetMask ) << OffsetShift2;
	idx |= ( ( (unsigned long long)(offset[2]) ) & OffsetMask ) << OffsetShift3;
	return idx;
}
template< class NodeData >
inline int OctNode< NodeData >::depth( void ) const {return int( _depthAndOffset & DepthMask );}
template< class NodeData >
inline void OctNode< NodeData >::DepthAndOffset(const long long& index,int& depth,int offset[3]){
	depth=int(index&DepthMask);//先求出当前节点的depth
	offset[0]=(int((index>>OffsetShift1)&OffsetMask)+1)&(~(1<<depth));//后一个位与操作相当于减去1<<depth
	offset[1]=(int((index>>OffsetShift2)&OffsetMask)+1)&(~(1<<depth));//那么这里的index应该记录的是当前node在octree中的绝对索引，现在要转化成相对的offset
	offset[2]=(int((index>>OffsetShift3)&OffsetMask)+1)&(~(1<<depth));
}
template< class NodeData >
inline int OctNode< NodeData >::Depth(const long long& index){return int(index&DepthMask);}
template< class NodeData >
template< class Real >
void OctNode< NodeData >::centerAndWidth( Point3D<Real>& center , Real& width ) const
{
	int depth , offset[3];
	depthAndOffset( depth , offset );//找出当前节点的深度和offset
	width = Real( 1.0 / (1<<depth) );//因为进行了归一化操作，因此节点的width与深度关系可以直接计算
	for( int dim=0 ; dim<DIMENSION ; dim++ ) center.coords[dim] = Real( 0.5+offset[dim] ) * width;//各个坐标轴的中心点也与offset相关
}
template< class NodeData >
template< class Real >
bool OctNode< NodeData >::isInside( Point3D< Real > p ) const//根据中心点位置和width判断给定数据点是否在OctNode内
{
	Point3D< Real > c;
	Real w;
	centerAndWidth( c , w );//中心点位置和width
	w /= 2;
	return (c[0]-w)<p[0] && p[0]<=(c[0]+w) && (c[1]-w)<p[1] && p[1]<=(c[1]+w) && (c[2]-w)<p[2] && p[2]<=(c[2]+w);
}
template< class NodeData >
template< class Real >
inline void OctNode< NodeData >::CenterAndWidth(const long long& index,Point3D<Real>& center,Real& width){
	int depth,offset[3];
	depth=index&DepthMask;//利用绝对索引获得center和width
	offset[0]=(int((index>>OffsetShift1)&OffsetMask)+1)&(~(1<<depth));
	offset[1]=(int((index>>OffsetShift2)&OffsetMask)+1)&(~(1<<depth));
	offset[2]=(int((index>>OffsetShift3)&OffsetMask)+1)&(~(1<<depth));
	width=Real(1.0/(1<<depth));
	for(int dim=0;dim<DIMENSION;dim++){center.coords[dim]=Real(0.5+offset[dim])*width;}
}

template< class NodeData >
int OctNode< NodeData >::maxDepth(void) const{
	//返回当前节点的所有子节点的最大深度
	if(!children){return 0;}
	else{
		int c,d;
		for(int i=0;i<Cube::CORNERS;i++){
			d=children[i].maxDepth();
			if(!i || d>c){c=d;}
		}
		return c+1;
	}
}
template< class NodeData >
size_t OctNode< NodeData >::nodes( void ) const
{
	//返回当前节点的所有子节点数目，包括自己
	if( !children ) return 1;
	else
	{
		size_t c=0;
		for( int i=0 ; i<Cube::CORNERS ; i++ ) c += children[i].nodes();
		return c+1;
	}
}
template< class NodeData >
size_t OctNode< NodeData >::leaves( void ) const
{
	if( !children ) return 1;
	else
	{
		size_t c=0;
		for( int i=0 ; i<Cube::CORNERS ; i++ ) c += children[i].leaves();
		return c;
	}
}
template< class NodeData >
size_t OctNode< NodeData >::maxDepthLeaves( int maxDepth ) const
{
	//给出在深度小于maxDepth时当前节点包含的叶节点数目
	if( depth()>maxDepth ) return 0;
	if( !children ) return 1;
	else
	{
		size_t c=0;
		for( int i=0 ; i<Cube::CORNERS ; i++ ) c += children[i].maxDepthLeaves(maxDepth);
		return c;
	}
}
template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::root(void) const{//返回根节点
	const OctNode* temp=this;
	while(temp->parent){temp=temp->parent;}
	return temp;
}


template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::nextBranch( const OctNode* current ) const
{
	//如果current不存在父节点或者current就是当前节点，就返回NULL
	//不太明白为什么current等于当前节点就返回NULL
	if( !current->parent || current==this ) return NULL;
	if(current-current->parent->children==Cube::CORNERS-1) return nextBranch( current->parent );//有可能会返回其父节点的下一个sibling
	else return current+1;//如果current不是其父节点的最后一个子节点，就返回其下一个sibling
}
template< class NodeData >
OctNode< NodeData >* OctNode< NodeData >::nextBranch(OctNode* current){
	if(!current->parent || current==this){return NULL;}
	if(current-current->parent->children==Cube::CORNERS-1){return nextBranch(current->parent);}
	else{return current+1;}
}
template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::prevBranch( const OctNode* current ) const
{
	//如果current不存在父节点或者current就是当前节点，就返回NULL
	if( !current->parent || current==this ) return NULL;
	if( current-current->parent->children==0 ) return prevBranch( current->parent );//有可能会返回其父节点的上一个sibling
	else return current-1;//如果current不是其父节点的第一个子节点，就返回其上一个sibling
}
template< class NodeData >
OctNode< NodeData >* OctNode< NodeData >::prevBranch( OctNode* current )
{
	if( !current->parent || current==this ) return NULL;
	if( current-current->parent->children==0 ) return prevBranch( current->parent );
	else return current-1;
}
template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::nextLeaf(const OctNode* current) const{
	if(!current){//如果没指明current，则返回当前节点的第一个子节点，有可能是自己
		const OctNode< NodeData >* temp=this;
		while(temp->children){temp=&temp->children[0];}
		return temp;
	}
	if(current->children){return current->nextLeaf();}//如果当期节点存在子节点，则返回子节点中的第一个
	const OctNode* temp=nextBranch(current);//如果current为根节点或者current节点就是当前节点，返回值为NULL
	if(!temp){return NULL;}
	else{return temp->nextLeaf();}//否则找出下一个sibling的叶节点并返回
}
template< class NodeData >
OctNode< NodeData >* OctNode< NodeData >::nextLeaf(OctNode* current){//注释同上
	if(!current){
		OctNode< NodeData >* temp=this;
		while(temp->children){temp=&temp->children[0];}
		return temp;
	}
	if(current->children){return current->nextLeaf();}
	OctNode* temp=nextBranch(current);
	if(!temp){return NULL;}
	else{return temp->nextLeaf();}
}

template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::nextNode( const OctNode* current ) const
{
	//找出下一个会被访问到的节点
	if( !current ) return this;
	else if( current->children ) return &current->children[0];
	else return nextBranch(current);
}
template< class NodeData >
OctNode< NodeData >* OctNode< NodeData >::nextNode( OctNode* current )
{
	//找出下一个会被访问到的节点
	if( !current ) return this;
	else if( current->children ) return &current->children[0];
	else return nextBranch( current );
}

template< class NodeData >
void OctNode< NodeData >::printRange(void) const
{
	Point3D< float > center;
	float width;
	centerAndWidth(center,width);//当前节点的center和width
	for(int dim=0;dim<DIMENSION;dim++){
		printf("%[%f,%f]",center.coords[dim]-width/2,center.coords[dim]+width/2);//当前节点在归一化octree中的距离范围
		if(dim<DIMENSION-1){printf("x");}
		else printf("\n");
	}
}

template< class NodeData >
void OctNode< NodeData >::AdjacencyCountFunction::Function(const OctNode< NodeData >* node1,const OctNode< NodeData >* node2){count++;}

template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::processNodeNodes(OctNode* node,NodeAdjacencyFunction* F,int processCurrent){
	if(processCurrent){F->Function(this,node);}
	if(children){__processNodeNodes(node,F);}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::processNodeFaces(OctNode* node,NodeAdjacencyFunction* F,int fIndex,int processCurrent){
	if(processCurrent){F->Function(this,node);}
	if(children){
		int c1,c2,c3,c4;
		Cube::FaceCorners(fIndex,c1,c2,c3,c4);
		__processNodeFaces(node,F,c1,c2,c3,c4);
	}
}
template< class NodeData >
template< class NodeAdjacencyFunction >
void OctNode< NodeData >::processNodeFaces( const OctNode* node , NodeAdjacencyFunction* F , int fIndex , int processCurrent ) const
{
	if( processCurrent ) F->Function( this , node );
	if(children)
	{
		int c1 , c2 , c3 , c4;
		Cube::FaceCorners( fIndex , c1 , c2 , c3 , c4 );
		__processNodeFaces( node , F , c1 , c2 , c3 , c4 );
	}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::processNodeEdges(OctNode* node,NodeAdjacencyFunction* F,int eIndex,int processCurrent){
	if(processCurrent){F->Function(this,node);}
	if(children){
		int c1,c2;
		Cube::EdgeCorners(eIndex,c1,c2);
		__processNodeEdges(node,F,c1,c2);
	}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::processNodeCorners(OctNode* node,NodeAdjacencyFunction* F,int cIndex,int processCurrent){
	if(processCurrent){F->Function(this,node);}
	OctNode< NodeData >* temp=this;
	while(temp->children){
		temp=&temp->children[cIndex];
		F->Function(temp,node);
	}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::__processNodeNodes(OctNode* node,NodeAdjacencyFunction* F){
	F->Function(&children[0],node);
	F->Function(&children[1],node);
	F->Function(&children[2],node);
	F->Function(&children[3],node);
	F->Function(&children[4],node);
	F->Function(&children[5],node);
	F->Function(&children[6],node);
	F->Function(&children[7],node);
	if(children[0].children){children[0].__processNodeNodes(node,F);}
	if(children[1].children){children[1].__processNodeNodes(node,F);}
	if(children[2].children){children[2].__processNodeNodes(node,F);}
	if(children[3].children){children[3].__processNodeNodes(node,F);}
	if(children[4].children){children[4].__processNodeNodes(node,F);}
	if(children[5].children){children[5].__processNodeNodes(node,F);}
	if(children[6].children){children[6].__processNodeNodes(node,F);}
	if(children[7].children){children[7].__processNodeNodes(node,F);}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::__processNodeEdges(OctNode* node,NodeAdjacencyFunction* F,int cIndex1,int cIndex2){
	F->Function(&children[cIndex1],node);
	F->Function(&children[cIndex2],node);
	if(children[cIndex1].children){children[cIndex1].__processNodeEdges(node,F,cIndex1,cIndex2);}
	if(children[cIndex2].children){children[cIndex2].__processNodeEdges(node,F,cIndex1,cIndex2);}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::__processNodeFaces(OctNode* node,NodeAdjacencyFunction* F,int cIndex1,int cIndex2,int cIndex3,int cIndex4){
	F->Function(&children[cIndex1],node);
	F->Function(&children[cIndex2],node);
	F->Function(&children[cIndex3],node);
	F->Function(&children[cIndex4],node);
	if(children[cIndex1].children){children[cIndex1].__processNodeFaces(node,F,cIndex1,cIndex2,cIndex3,cIndex4);}
	if(children[cIndex2].children){children[cIndex2].__processNodeFaces(node,F,cIndex1,cIndex2,cIndex3,cIndex4);}
	if(children[cIndex3].children){children[cIndex3].__processNodeFaces(node,F,cIndex1,cIndex2,cIndex3,cIndex4);}
	if(children[cIndex4].children){children[cIndex4].__processNodeFaces(node,F,cIndex1,cIndex2,cIndex3,cIndex4);}
}
template< class NodeData >
template< class NodeAdjacencyFunction >
void OctNode< NodeData >::__processNodeFaces( const OctNode* node , NodeAdjacencyFunction* F , int cIndex1 , int cIndex2 , int cIndex3 , int cIndex4 ) const
{
	F->Function( &children[cIndex1] , node );
	F->Function( &children[cIndex2] , node );
	F->Function( &children[cIndex3] , node );
	F->Function( &children[cIndex4] , node );
	if( children[cIndex1].children ) children[cIndex1].__processNodeFaces( node , F , cIndex1 , cIndex2 , cIndex3 , cIndex4 );
	if( children[cIndex2].children ) children[cIndex2].__processNodeFaces( node , F , cIndex1 , cIndex2 , cIndex3 , cIndex4 );
	if( children[cIndex3].children ) children[cIndex3].__processNodeFaces( node , F , cIndex1 , cIndex2 , cIndex3 , cIndex4 );
	if( children[cIndex4].children ) children[cIndex4].__processNodeFaces( node , F , cIndex1 , cIndex2 , cIndex3 , cIndex4 );
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::ProcessNodeAdjacentNodes(int maxDepth,OctNode* node1,int width1,OctNode* node2,int width2,NodeAdjacencyFunction* F,int processCurrent){
	int c1[3],c2[3],w1,w2;
	node1->centerIndex(maxDepth+1,c1);
	node2->centerIndex(maxDepth+1,c2);
	w1=node1->width(maxDepth+1);
	w2=node2->width(maxDepth+1);

	ProcessNodeAdjacentNodes(c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2],node1,(width1*w1)>>1,node2,(width2*w2)>>1,w2,F,processCurrent);
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::ProcessNodeAdjacentNodes(int dx,int dy,int dz,
													  OctNode< NodeData >* node1,int radius1,
													  OctNode< NodeData >* node2,int radius2,int width2,
													  NodeAdjacencyFunction* F,int processCurrent){
	if(!Overlap(dx,dy,dz,radius1+radius2)){return;}
	if(processCurrent){F->Function(node2,node1);}
	if(!node2->children){return;}
	__ProcessNodeAdjacentNodes(-dx,-dy,-dz,node1,radius1,node2,radius2,width2/2,F);
}
template< class NodeData >
template<class TerminatingNodeAdjacencyFunction>
void OctNode< NodeData >::ProcessTerminatingNodeAdjacentNodes(int maxDepth,OctNode* node1,int width1,OctNode* node2,int width2,TerminatingNodeAdjacencyFunction* F,int processCurrent){
	int c1[3],c2[3],w1,w2;
	node1->centerIndex(maxDepth+1,c1);
	node2->centerIndex(maxDepth+1,c2);
	w1=node1->width(maxDepth+1);
	w2=node2->width(maxDepth+1);

	ProcessTerminatingNodeAdjacentNodes(c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2],node1,(width1*w1)>>1,node2,(width2*w2)>>1,w2,F,processCurrent);
}
template< class NodeData >
template<class TerminatingNodeAdjacencyFunction>
void OctNode< NodeData >::ProcessTerminatingNodeAdjacentNodes(int dx,int dy,int dz,
																 OctNode< NodeData >* node1,int radius1,
																 OctNode< NodeData >* node2,int radius2,int width2,
																 TerminatingNodeAdjacencyFunction* F,int processCurrent)
{
	if(!Overlap(dx,dy,dz,radius1+radius2)){return;}
	if(processCurrent){F->Function(node2,node1);}
	if(!node2->children){return;}
	__ProcessTerminatingNodeAdjacentNodes(-dx,-dy,-dz,node1,radius1,node2,radius2,width2/2,F);
}
template< class NodeData >
template<class PointAdjacencyFunction>
void OctNode< NodeData >::ProcessPointAdjacentNodes( int maxDepth , const int c1[3] , OctNode* node2 , int width2 , PointAdjacencyFunction* F , int processCurrent )
{
	int c2[3] , w2;
	node2->centerIndex( maxDepth+1 , c2 );
	w2 = node2->width( maxDepth+1 );
	ProcessPointAdjacentNodes( c1[0]-c2[0] , c1[1]-c2[1] , c1[2]-c2[2] , node2 , (width2*w2)>>1 , w2 , F , processCurrent );
}
template< class NodeData >
template<class PointAdjacencyFunction>
void OctNode< NodeData >::ProcessPointAdjacentNodes(int dx,int dy,int dz,
													   OctNode< NodeData >* node2,int radius2,int width2,
													   PointAdjacencyFunction* F,int processCurrent)
{
	if( !Overlap(dx,dy,dz,radius2) ) return;
	if( processCurrent ) F->Function(node2);
	if( !node2->children ) return;
	__ProcessPointAdjacentNodes( -dx , -dy , -dz , node2 , radius2 , width2>>1 , F );
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::ProcessFixedDepthNodeAdjacentNodes(int maxDepth,
																OctNode< NodeData >* node1,int width1,
																OctNode< NodeData >* node2,int width2,
																int depth,NodeAdjacencyFunction* F,int processCurrent)
{
	int c1[3],c2[3],w1,w2;
	node1->centerIndex(maxDepth+1,c1);
	node2->centerIndex(maxDepth+1,c2);
	w1=node1->width(maxDepth+1);
	w2=node2->width(maxDepth+1);

	ProcessFixedDepthNodeAdjacentNodes(c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2],node1,(width1*w1)>>1,node2,(width2*w2)>>1,w2,depth,F,processCurrent);
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::ProcessFixedDepthNodeAdjacentNodes(int dx,int dy,int dz,
																OctNode< NodeData >* node1,int radius1,
																OctNode< NodeData >* node2,int radius2,int width2,
																int depth,NodeAdjacencyFunction* F,int processCurrent)
{
	int d=node2->depth();
	if(d>depth){return;}
	if(!Overlap(dx,dy,dz,radius1+radius2)){return;}
	if(d==depth){if(processCurrent){F->Function(node2,node1);}}
	else{
		if(!node2->children){return;}
		__ProcessFixedDepthNodeAdjacentNodes(-dx,-dy,-dz,node1,radius1,node2,radius2,width2/2,depth-1,F);
	}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::ProcessMaxDepthNodeAdjacentNodes(int maxDepth,
															  OctNode< NodeData >* node1,int width1,
															  OctNode< NodeData >* node2,int width2,
															  int depth,NodeAdjacencyFunction* F,int processCurrent)
{
	int c1[3],c2[3],w1,w2;
	node1->centerIndex(maxDepth+1,c1);
	node2->centerIndex(maxDepth+1,c2);
	w1=node1->width(maxDepth+1);
	w2=node2->width(maxDepth+1);
	ProcessMaxDepthNodeAdjacentNodes(c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2],node1,(width1*w1)>>1,node2,(width2*w2)>>1,w2,depth,F,processCurrent);
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::ProcessMaxDepthNodeAdjacentNodes(int dx,int dy,int dz,
															  OctNode< NodeData >* node1,int radius1,
															  OctNode< NodeData >* node2,int radius2,int width2,
															  int depth,NodeAdjacencyFunction* F,int processCurrent)
{
	int d=node2->depth();
	if(d>depth){return;}
	if(!Overlap(dx,dy,dz,radius1+radius2)){return;}
	if(processCurrent){F->Function(node2,node1);}
	if(d<depth && node2->children){__ProcessMaxDepthNodeAdjacentNodes(-dx,-dy,-dz,node1,radius1,node2,radius2,width2>>1,depth-1,F);}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::__ProcessNodeAdjacentNodes(int dx,int dy,int dz,
														OctNode* node1,int radius1,
														OctNode* node2,int radius2,int cWidth2,
														NodeAdjacencyFunction* F)
{
	int cWidth=cWidth2>>1;
	int radius=radius2>>1;
	int o=ChildOverlap(dx,dy,dz,radius1+radius,cWidth);
	if(o){
		int dx1=dx-cWidth;
		int dx2=dx+cWidth;
		int dy1=dy-cWidth;
		int dy2=dy+cWidth;
		int dz1=dz-cWidth;
		int dz2=dz+cWidth;
		if(o&  1){F->Function(&node2->children[0],node1);if(node2->children[0].children){__ProcessNodeAdjacentNodes(dx1,dy1,dz1,node1,radius1,&node2->children[0],radius,cWidth,F);}}
		if(o&  2){F->Function(&node2->children[1],node1);if(node2->children[1].children){__ProcessNodeAdjacentNodes(dx2,dy1,dz1,node1,radius1,&node2->children[1],radius,cWidth,F);}}
		if(o&  4){F->Function(&node2->children[2],node1);if(node2->children[2].children){__ProcessNodeAdjacentNodes(dx1,dy2,dz1,node1,radius1,&node2->children[2],radius,cWidth,F);}}
		if(o&  8){F->Function(&node2->children[3],node1);if(node2->children[3].children){__ProcessNodeAdjacentNodes(dx2,dy2,dz1,node1,radius1,&node2->children[3],radius,cWidth,F);}}
		if(o& 16){F->Function(&node2->children[4],node1);if(node2->children[4].children){__ProcessNodeAdjacentNodes(dx1,dy1,dz2,node1,radius1,&node2->children[4],radius,cWidth,F);}}
		if(o& 32){F->Function(&node2->children[5],node1);if(node2->children[5].children){__ProcessNodeAdjacentNodes(dx2,dy1,dz2,node1,radius1,&node2->children[5],radius,cWidth,F);}}
		if(o& 64){F->Function(&node2->children[6],node1);if(node2->children[6].children){__ProcessNodeAdjacentNodes(dx1,dy2,dz2,node1,radius1,&node2->children[6],radius,cWidth,F);}}
		if(o&128){F->Function(&node2->children[7],node1);if(node2->children[7].children){__ProcessNodeAdjacentNodes(dx2,dy2,dz2,node1,radius1,&node2->children[7],radius,cWidth,F);}}
	}
}
template< class NodeData >
template<class TerminatingNodeAdjacencyFunction>
void OctNode< NodeData >::__ProcessTerminatingNodeAdjacentNodes(int dx,int dy,int dz,
																   OctNode* node1,int radius1,
																   OctNode* node2,int radius2,int cWidth2,
																   TerminatingNodeAdjacencyFunction* F)
{
	int cWidth=cWidth2>>1;
	int radius=radius2>>1;
	int o=ChildOverlap(dx,dy,dz,radius1+radius,cWidth);
	if(o){
		int dx1=dx-cWidth;
		int dx2=dx+cWidth;
		int dy1=dy-cWidth;
		int dy2=dy+cWidth;
		int dz1=dz-cWidth;
		int dz2=dz+cWidth;
		if(o&  1){if(F->Function(&node2->children[0],node1) && node2->children[0].children){__ProcessTerminatingNodeAdjacentNodes(dx1,dy1,dz1,node1,radius1,&node2->children[0],radius,cWidth,F);}}
		if(o&  2){if(F->Function(&node2->children[1],node1) && node2->children[1].children){__ProcessTerminatingNodeAdjacentNodes(dx2,dy1,dz1,node1,radius1,&node2->children[1],radius,cWidth,F);}}
		if(o&  4){if(F->Function(&node2->children[2],node1) && node2->children[2].children){__ProcessTerminatingNodeAdjacentNodes(dx1,dy2,dz1,node1,radius1,&node2->children[2],radius,cWidth,F);}}
		if(o&  8){if(F->Function(&node2->children[3],node1) && node2->children[3].children){__ProcessTerminatingNodeAdjacentNodes(dx2,dy2,dz1,node1,radius1,&node2->children[3],radius,cWidth,F);}}
		if(o& 16){if(F->Function(&node2->children[4],node1) && node2->children[4].children){__ProcessTerminatingNodeAdjacentNodes(dx1,dy1,dz2,node1,radius1,&node2->children[4],radius,cWidth,F);}}
		if(o& 32){if(F->Function(&node2->children[5],node1) && node2->children[5].children){__ProcessTerminatingNodeAdjacentNodes(dx2,dy1,dz2,node1,radius1,&node2->children[5],radius,cWidth,F);}}
		if(o& 64){if(F->Function(&node2->children[6],node1) && node2->children[6].children){__ProcessTerminatingNodeAdjacentNodes(dx1,dy2,dz2,node1,radius1,&node2->children[6],radius,cWidth,F);}}
		if(o&128){if(F->Function(&node2->children[7],node1) && node2->children[7].children){__ProcessTerminatingNodeAdjacentNodes(dx2,dy2,dz2,node1,radius1,&node2->children[7],radius,cWidth,F);}}
	}
}
template< class NodeData >
template<class PointAdjacencyFunction>
void OctNode< NodeData >::__ProcessPointAdjacentNodes(int dx,int dy,int dz,
														 OctNode* node2,int radius2,int cWidth2,
														 PointAdjacencyFunction* F)
{
	int cWidth=cWidth2>>1;
	int radius=radius2>>1;
	int o=ChildOverlap(dx,dy,dz,radius,cWidth);
	if( o )
	{
		int dx1=dx-cWidth;
		int dx2=dx+cWidth;
		int dy1=dy-cWidth;
		int dy2=dy+cWidth;
		int dz1=dz-cWidth;
		int dz2=dz+cWidth;
		if(o&  1){F->Function(&node2->children[0]);if(node2->children[0].children){__ProcessPointAdjacentNodes(dx1,dy1,dz1,&node2->children[0],radius,cWidth,F);}}
		if(o&  2){F->Function(&node2->children[1]);if(node2->children[1].children){__ProcessPointAdjacentNodes(dx2,dy1,dz1,&node2->children[1],radius,cWidth,F);}}
		if(o&  4){F->Function(&node2->children[2]);if(node2->children[2].children){__ProcessPointAdjacentNodes(dx1,dy2,dz1,&node2->children[2],radius,cWidth,F);}}
		if(o&  8){F->Function(&node2->children[3]);if(node2->children[3].children){__ProcessPointAdjacentNodes(dx2,dy2,dz1,&node2->children[3],radius,cWidth,F);}}
		if(o& 16){F->Function(&node2->children[4]);if(node2->children[4].children){__ProcessPointAdjacentNodes(dx1,dy1,dz2,&node2->children[4],radius,cWidth,F);}}
		if(o& 32){F->Function(&node2->children[5]);if(node2->children[5].children){__ProcessPointAdjacentNodes(dx2,dy1,dz2,&node2->children[5],radius,cWidth,F);}}
		if(o& 64){F->Function(&node2->children[6]);if(node2->children[6].children){__ProcessPointAdjacentNodes(dx1,dy2,dz2,&node2->children[6],radius,cWidth,F);}}
		if(o&128){F->Function(&node2->children[7]);if(node2->children[7].children){__ProcessPointAdjacentNodes(dx2,dy2,dz2,&node2->children[7],radius,cWidth,F);}}
	}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::__ProcessFixedDepthNodeAdjacentNodes(int dx,int dy,int dz,
																  OctNode* node1,int radius1,
																  OctNode* node2,int radius2,int cWidth2,
																  int depth,NodeAdjacencyFunction* F)
{
	int cWidth=cWidth2>>1;
	int radius=radius2>>1;
	int o=ChildOverlap(dx,dy,dz,radius1+radius,cWidth);
	if(o){
		int dx1=dx-cWidth;
		int dx2=dx+cWidth;
		int dy1=dy-cWidth;
		int dy2=dy+cWidth;
		int dz1=dz-cWidth;
		int dz2=dz+cWidth;
		if(node2->depth()==depth){
			if(o&  1){F->Function(&node2->children[0],node1);}
			if(o&  2){F->Function(&node2->children[1],node1);}
			if(o&  4){F->Function(&node2->children[2],node1);}
			if(o&  8){F->Function(&node2->children[3],node1);}
			if(o& 16){F->Function(&node2->children[4],node1);}
			if(o& 32){F->Function(&node2->children[5],node1);}
			if(o& 64){F->Function(&node2->children[6],node1);}
			if(o&128){F->Function(&node2->children[7],node1);}
		}
		else{
			if(o&  1){if(node2->children[0].children){__ProcessFixedDepthNodeAdjacentNodes(dx1,dy1,dz1,node1,radius1,&node2->children[0],radius,cWidth,depth,F);}}
			if(o&  2){if(node2->children[1].children){__ProcessFixedDepthNodeAdjacentNodes(dx2,dy1,dz1,node1,radius1,&node2->children[1],radius,cWidth,depth,F);}}
			if(o&  4){if(node2->children[2].children){__ProcessFixedDepthNodeAdjacentNodes(dx1,dy2,dz1,node1,radius1,&node2->children[2],radius,cWidth,depth,F);}}
			if(o&  8){if(node2->children[3].children){__ProcessFixedDepthNodeAdjacentNodes(dx2,dy2,dz1,node1,radius1,&node2->children[3],radius,cWidth,depth,F);}}
			if(o& 16){if(node2->children[4].children){__ProcessFixedDepthNodeAdjacentNodes(dx1,dy1,dz2,node1,radius1,&node2->children[4],radius,cWidth,depth,F);}}
			if(o& 32){if(node2->children[5].children){__ProcessFixedDepthNodeAdjacentNodes(dx2,dy1,dz2,node1,radius1,&node2->children[5],radius,cWidth,depth,F);}}
			if(o& 64){if(node2->children[6].children){__ProcessFixedDepthNodeAdjacentNodes(dx1,dy2,dz2,node1,radius1,&node2->children[6],radius,cWidth,depth,F);}}
			if(o&128){if(node2->children[7].children){__ProcessFixedDepthNodeAdjacentNodes(dx2,dy2,dz2,node1,radius1,&node2->children[7],radius,cWidth,depth,F);}}
		}
	}
}
template< class NodeData >
template<class NodeAdjacencyFunction>
void OctNode< NodeData >::__ProcessMaxDepthNodeAdjacentNodes(int dx,int dy,int dz,
																OctNode* node1,int radius1,
																OctNode* node2,int radius2,int cWidth2,
																int depth,NodeAdjacencyFunction* F)
{
	int cWidth=cWidth2>>1;
	int radius=radius2>>1;
	int o=ChildOverlap(dx,dy,dz,radius1+radius,cWidth);
	if(o){
		int dx1=dx-cWidth;
		int dx2=dx+cWidth;
		int dy1=dy-cWidth;
		int dy2=dy+cWidth;
		int dz1=dz-cWidth;
		int dz2=dz+cWidth;
		if(node2->depth()<=depth){
			if(o&  1){F->Function(&node2->children[0],node1);}
			if(o&  2){F->Function(&node2->children[1],node1);}
			if(o&  4){F->Function(&node2->children[2],node1);}
			if(o&  8){F->Function(&node2->children[3],node1);}
			if(o& 16){F->Function(&node2->children[4],node1);}
			if(o& 32){F->Function(&node2->children[5],node1);}
			if(o& 64){F->Function(&node2->children[6],node1);}
			if(o&128){F->Function(&node2->children[7],node1);}
		}
		if(node2->depth()<depth){
			if(o&  1){if(node2->children[0].children){__ProcessMaxDepthNodeAdjacentNodes(dx1,dy1,dz1,node1,radius1,&node2->children[0],radius,cWidth,depth,F);}}
			if(o&  2){if(node2->children[1].children){__ProcessMaxDepthNodeAdjacentNodes(dx2,dy1,dz1,node1,radius1,&node2->children[1],radius,cWidth,depth,F);}}
			if(o&  4){if(node2->children[2].children){__ProcessMaxDepthNodeAdjacentNodes(dx1,dy2,dz1,node1,radius1,&node2->children[2],radius,cWidth,depth,F);}}
			if(o&  8){if(node2->children[3].children){__ProcessMaxDepthNodeAdjacentNodes(dx2,dy2,dz1,node1,radius1,&node2->children[3],radius,cWidth,depth,F);}}
			if(o& 16){if(node2->children[4].children){__ProcessMaxDepthNodeAdjacentNodes(dx1,dy1,dz2,node1,radius1,&node2->children[4],radius,cWidth,depth,F);}}
			if(o& 32){if(node2->children[5].children){__ProcessMaxDepthNodeAdjacentNodes(dx2,dy1,dz2,node1,radius1,&node2->children[5],radius,cWidth,depth,F);}}
			if(o& 64){if(node2->children[6].children){__ProcessMaxDepthNodeAdjacentNodes(dx1,dy2,dz2,node1,radius1,&node2->children[6],radius,cWidth,depth,F);}}
			if(o&128){if(node2->children[7].children){__ProcessMaxDepthNodeAdjacentNodes(dx2,dy2,dz2,node1,radius1,&node2->children[7],radius,cWidth,depth,F);}}
		}
	}
}
template< class NodeData >
inline int OctNode< NodeData >::ChildOverlap(int dx,int dy,int dz,int d,int cRadius2)
{
	int w1=d-cRadius2;
	int w2=d+cRadius2;
	int overlap=0;

	int test=0,test1=0;
	if(dx<w2 && dx>-w1){test =1;}
	if(dx<w1 && dx>-w2){test|=2;}

	if(!test){return 0;}
	if(dz<w2 && dz>-w1){test1 =test;}
	if(dz<w1 && dz>-w2){test1|=test<<4;}

	if(!test1){return 0;}
	if(dy<w2 && dy>-w1){overlap =test1;}
	if(dy<w1 && dy>-w2){overlap|=test1<<2;}
	return overlap;
}

template< class NodeData >
template< class Real >
OctNode< NodeData >* OctNode< NodeData >::getNearestLeaf(const Point3D<Real>& p){
	//找出离p最近的当前节点下的子节点，有可能是自己
	Point3D<Real> center;
	Real width;
	OctNode< NodeData >* temp;
	int cIndex;
	if(!children){return this;}
	centerAndWidth(center,width);
	temp=this;
	while(temp->children){
		cIndex=CornerIndex(center,p);//返回离p最近的corner index
		temp=&temp->children[cIndex];
		width/=2;
		if(cIndex&1){center.coords[0]+=width/2;}
		else		{center.coords[0]-=width/2;}
		if(cIndex&2){center.coords[1]+=width/2;}
		else		{center.coords[1]-=width/2;}
		if(cIndex&4){center.coords[2]+=width/2;}
		else		{center.coords[2]-=width/2;}
	}
	return temp;
}
template< class NodeData >
template< class Real >
const OctNode< NodeData >* OctNode< NodeData >::getNearestLeaf(const Point3D<Real>& p) const{
	int nearest;
	Real temp,dist2;
	if(!children){return this;}
	for(int i=0;i<Cube::CORNERS;i++){
		temp=SquareDistance(children[i].center,p);//这个通过欧式距离来查找，也是迭代查找
		if(!i || temp<dist2){
			dist2=temp;
			nearest=i;
		}
	}
	return children[nearest].getNearestLeaf(p);
}

template< class NodeData >
int OctNode< NodeData >::CommonEdge(const OctNode< NodeData >* node1,int eIndex1,const OctNode< NodeData >* node2,int eIndex2){
	//经验证是正确的，但是思路有点复杂，目的是判断给出的两条边是否在共同的父节点上处于同一条公共边上，而不是判断两条边是否是重合的边
	int o1,o2,i1,i2,j1,j2;

	Cube::FactorEdgeIndex(eIndex1,o1,i1,j1);//找出node1上eIndex1指定的边的orientation等信息
	Cube::FactorEdgeIndex(eIndex2,o2,i2,j2);//同上
	if(o1!=o2){return 0;}//如果两个边在不同的坐标轴上，肯定不是公共边，返回0

	int dir[2];
	int idx1[2];
	int idx2[2];
	switch(o1){
		case 0:	dir[0]=1;	dir[1]=2;	break;
		case 1:	dir[0]=0;	dir[1]=2;	break;
		case 2:	dir[0]=0;	dir[1]=1;	break;
	};
	int d1,d2,off1[3],off2[3];
	node1->depthAndOffset(d1,off1);//返回node1的depth和offset
	node2->depthAndOffset(d2,off2);//同上
	idx1[0]=off1[dir[0]]+(1<<d1)+i1;
	idx1[1]=off1[dir[1]]+(1<<d1)+j1;
	idx2[0]=off2[dir[0]]+(1<<d2)+i2;
	idx2[1]=off2[dir[1]]+(1<<d2)+j2;
	if(d1>d2){//证明两个node还有可能不在同一个depth上
		idx2[0]<<=(d1-d2);
		idx2[1]<<=(d1-d2);
	}
	else{
		idx1[0]<<=(d2-d1);
		idx1[1]<<=(d2-d1);
	}
	if(idx1[0]==idx2[0] && idx1[1]==idx2[1]){return 1;}
	else									{return 0;}
}
template< class NodeData >
template< class Real >
int OctNode< NodeData >::CornerIndex(const Point3D<Real>& center,const Point3D<Real>& p){
	int cIndex=0;//根据点的位置与中心点位置关系来确定corner index，找出当前节点中离p最近的corner
	if(p.coords[0]>center.coords[0]){cIndex|=1;}
	if(p.coords[1]>center.coords[1]){cIndex|=2;}
	if(p.coords[2]>center.coords[2]){cIndex|=4;}
	return cIndex;
}
template< class NodeData >
template< class NodeData2 >
OctNode< NodeData >& OctNode< NodeData >::operator = ( const OctNode< NodeData2 >& node )
{
	int i;
	if(children){delete[] children;}
	children=NULL;

	this->depth = node.depth;
	for(i=0;i<DIMENSION;i++){this->offset[i] = node.offset[i];}
	if(node.children){
		initChildren();
		for(i=0;i<Cube::CORNERS;i++){children[i] = node.children[i];}
	}
	return *this;
}
template< class NodeData >
int OctNode< NodeData >::CompareForwardDepths(const void* v1,const void* v2){
	//判断depth是否相同，难道对于形参是void的情况在调用时可以不用带括号
	return ((const OctNode< NodeData >*)v1)->depth-((const OctNode< NodeData >*)v2)->depth;
}
template< class NodeData >
int OctNode< NodeData >::CompareByDepthAndXYZ( const void* v1 , const void* v2 )
{
	//根据节点的深度以及offset的值判断两个node是否相同，不同还需要返回差异值
	const OctNode< NodeData > *n1 = (*(const OctNode< NodeData >**)v1);
	const OctNode< NodeData > *n2 = (*(const OctNode< NodeData >**)v2);
	if( n1->d!=n2->d ) return int(n1->d)-int(n2->d);
	else if( n1->off[0]!=n2->off[0] ) return int(n1->off[0]) - int(n2->off[0]);
	else if( n1->off[1]!=n2->off[1] ) return int(n1->off[1]) - int(n2->off[1]);
	else if( n1->off[2]!=n2->off[2] ) return int(n1->off[2]) - int(n2->off[2]);
	return 0;
}

long long _InterleaveBits( int p[3] )
{
	long long key = 0;
	for( int i=0 ; i<32 ; i++ ) key |= ( ( p[0] & (1<<i) )<<(2*i) ) | ( ( p[1] & (1<<i) )<<(2*i+1) ) | ( ( p[2] & (1<<i) )<<(2*i+2) );//这样p[0]下一位的取值移位后会与p[2]当前位数字重合
	return key;
}
template< class NodeData >
int OctNode< NodeData >::CompareByDepthAndZIndex( const void* v1 , const void* v2 )
{
	//暂时不明白这个函数的用意何在???
	const OctNode< NodeData >* n1 = (*(const OctNode< NodeData >**)v1);
	const OctNode< NodeData >* n2 = (*(const OctNode< NodeData >**)v2);
	int d1 , off1[3] , d2 , off2[3];
	n1->depthAndOffset( d1 , off1 ) , n2->depthAndOffset( d2 , off2 );
	if     ( d1>d2 ) return  1;
	else if( d1<d2 ) return -1;
	long long k1 = _InterleaveBits( off1 ) , k2 = _InterleaveBits( off2 );
	if     ( k1>k2 ) return  1;
	else if( k1<k2 ) return -1;
	else             return  0;
}

template< class NodeData >
int OctNode< NodeData >::CompareForwardPointerDepths( const void* v1 , const void* v2 )
{
	const OctNode< NodeData >* n1 = (*(const OctNode< NodeData >**)v1);
	const OctNode< NodeData >* n2 = (*(const OctNode< NodeData >**)v2);
	if(n1->d!=n2->d){return int(n1->d)-int(n2->d);}//先判断是否在同一深度上
	while( n1->parent!=n2->parent )//回溯父节点，直到遇到公共父节点
	{
		n1=n1->parent;
		n2=n2->parent;
	}
	//找出二者在父节点相同时的offset差异
	if(n1->off[0]!=n2->off[0]){return int(n1->off[0])-int(n2->off[0]);}
	if(n1->off[1]!=n2->off[1]){return int(n1->off[1])-int(n2->off[1]);}
	return int(n1->off[2])-int(n2->off[2]);
	return 0;
}
template< class NodeData >
int OctNode< NodeData >::CompareBackwardDepths(const void* v1,const void* v2){
	return ((const OctNode< NodeData >*)v2)->depth-((const OctNode< NodeData >*)v1)->depth;
}
template< class NodeData >
int OctNode< NodeData >::CompareBackwardPointerDepths(const void* v1,const void* v2){
	return (*(const OctNode< NodeData >**)v2)->depth()-(*(const OctNode< NodeData >**)v1)->depth();
}
template< class NodeData >
template< class Real >
inline int OctNode< NodeData >::Overlap2(const int &depth1,const int offSet1[DIMENSION],const Real& multiplier1,const int &depth2,const int offSet2[DIMENSION],const Real& multiplier2){
	int d=depth2-depth1;//depth2一定比depth1要大???
	Real w=multiplier2+multiplier1*(1<<d);
	Real w2=Real((1<<(d-1))-0.5);
	if(
		fabs(Real(offSet2[0]-(offSet1[0]<<d))-w2)>=w ||
		fabs(Real(offSet2[1]-(offSet1[1]<<d))-w2)>=w ||
		fabs(Real(offSet2[2]-(offSet1[2]<<d))-w2)>=w
		){return 0;}
	return 1;
}
template< class NodeData >
inline int OctNode< NodeData >::Overlap(int c1,int c2,int c3,int dWidth){
	if(c1>=dWidth || c1<=-dWidth || c2>=dWidth || c2<=-dWidth || c3>=dWidth || c3<=-dWidth){return 0;}
	else{return 1;}
}
template< class NodeData >
OctNode< NodeData >* OctNode< NodeData >::faceNeighbor(int faceIndex,int forceChildren){return __faceNeighbor(faceIndex>>1,faceIndex&1,forceChildren);}
template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::faceNeighbor(int faceIndex) const {return __faceNeighbor(faceIndex>>1,faceIndex&1);}
template< class NodeData >
OctNode< NodeData >* OctNode< NodeData >::__faceNeighbor(int dir,int off,int forceChildren){
	if(!parent){return NULL;}
	int pIndex=int(this-parent->children);
	pIndex^=(1<<dir);
	if((pIndex & (1<<dir))==(off<<dir)){return &parent->children[pIndex];}
	else{
		OctNode* temp=parent->__faceNeighbor(dir,off,forceChildren);
		if(!temp){return NULL;}
		if(!temp->children){
			if(forceChildren){temp->initChildren();}
			else{return temp;}
		}
		return &temp->children[pIndex];
	}
}
template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::__faceNeighbor(int dir,int off) const {
	if(!parent){return NULL;}
	int pIndex=int(this-parent->children);
	pIndex^=(1<<dir);
	if((pIndex & (1<<dir))==(off<<dir)){return &parent->children[pIndex];}
	else{
		const OctNode* temp=parent->__faceNeighbor(dir,off);
		if(!temp || !temp->children){return temp;}
		else{return &temp->children[pIndex];}
	}
}

template< class NodeData >
OctNode< NodeData >* OctNode< NodeData >::edgeNeighbor(int edgeIndex,int forceChildren){
	int idx[2],o,i[2];
	Cube::FactorEdgeIndex(edgeIndex,o,i[0],i[1]);//现根据edgeindex找出它的三个index坐标
	switch(o){
		case 0:	idx[0]=1;	idx[1]=2;	break;//假设分别代表xyz
		case 1:	idx[0]=0;	idx[1]=2;	break;
		case 2:	idx[0]=0;	idx[1]=1;	break;
	};
	return __edgeNeighbor(o,i,idx,forceChildren);
}
template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::edgeNeighbor(int edgeIndex) const {
	int idx[2],o,i[2];
	Cube::FactorEdgeIndex(edgeIndex,o,i[0],i[1]);
	switch(o){
		case 0:	idx[0]=1;	idx[1]=2;	break;
		case 1:	idx[0]=0;	idx[1]=2;	break;
		case 2:	idx[0]=0;	idx[1]=1;	break;
	};
	return __edgeNeighbor(o,i,idx);
}
template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::__edgeNeighbor(int o,const int i[2],const int idx[2]) const{
	if(!parent){return NULL;}
	int pIndex=int(this-parent->children);
	int aIndex,x[DIMENSION];

	Cube::FactorCornerIndex(pIndex,x[0],x[1],x[2]);
	aIndex=(~((i[0] ^ x[idx[0]]) | ((i[1] ^ x[idx[1]])<<1))) & 3;
	pIndex^=(7 ^ (1<<o));
	if(aIndex==1)	{	// I can get the neighbor from the parent's face adjacent neighbor
		const OctNode* temp=parent->__faceNeighbor(idx[0],i[0]);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==2)	{	// I can get the neighbor from the parent's face adjacent neighbor
		const OctNode* temp=parent->__faceNeighbor(idx[1],i[1]);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==0)	{	// I can get the neighbor from the parent
		return &parent->children[pIndex];
	}
	else if(aIndex==3)	{	// I can get the neighbor from the parent's edge adjacent neighbor
		const OctNode* temp=parent->__edgeNeighbor(o,i,idx);
		if(!temp || !temp->children){return temp;}
		else{return &temp->children[pIndex];}
	}
	else{return NULL;}
}
template< class NodeData >
OctNode< NodeData >* OctNode< NodeData >::__edgeNeighbor(int o,const int i[2],const int idx[2],int forceChildren){
	if(!parent){return NULL;}//没有父节点，没有Neighbor
	int pIndex=int(this-parent->children);//这是个什么操作，没看懂
	int aIndex,x[DIMENSION];

	Cube::FactorCornerIndex(pIndex,x[0],x[1],x[2]);//得到角点位置
	aIndex=(~((i[0] ^ x[idx[0]]) | ((i[1] ^ x[idx[1]])<<1))) & 3;//^代表异或
	pIndex^=(7 ^ (1<<o));
	if(aIndex==1)	{	// I can get the neighbor from the parent's face adjacent neighbor
		OctNode* temp=parent->__faceNeighbor(idx[0],i[0],0);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==2)	{	// I can get the neighbor from the parent's face adjacent neighbor
		OctNode* temp=parent->__faceNeighbor(idx[1],i[1],0);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==0)	{	// I can get the neighbor from the parent
		return &parent->children[pIndex];
	}
	else if(aIndex==3)	{	// I can get the neighbor from the parent's edge adjacent neighbor
		OctNode* temp=parent->__edgeNeighbor(o,i,idx,forceChildren);
		if(!temp){return NULL;}
		if(!temp->children){
			if(forceChildren){temp->initChildren();}
			else{return temp;}
		}
		return &temp->children[pIndex];
	}
	else{return NULL;}
}

template< class NodeData >
const OctNode< NodeData >* OctNode< NodeData >::cornerNeighbor(int cornerIndex) const {
	int pIndex,aIndex=0;
	if(!parent){return NULL;}

	pIndex=int(this-parent->children);
	aIndex=(cornerIndex ^ pIndex);	// The disagreement bits
	pIndex=(~pIndex)&7;				// The antipodal point
	if(aIndex==7){					// Agree on no bits
		return &parent->children[pIndex];
	}
	else if(aIndex==0){				// Agree on all bits
		const OctNode* temp=((const OctNode*)parent)->cornerNeighbor(cornerIndex);
		if(!temp || !temp->children){return temp;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==6){				// Agree on face 0
		const OctNode* temp=((const OctNode*)parent)->__faceNeighbor(0,cornerIndex & 1);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==5){				// Agree on face 1
		const OctNode* temp=((const OctNode*)parent)->__faceNeighbor(1,(cornerIndex & 2)>>1);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==3){				// Agree on face 2
		const OctNode* temp=((const OctNode*)parent)->__faceNeighbor(2,(cornerIndex & 4)>>2);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==4){				// Agree on edge 2
		const OctNode* temp=((const OctNode*)parent)->edgeNeighbor(8 | (cornerIndex & 1) | (cornerIndex & 2) );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==2){				// Agree on edge 1
		const OctNode* temp=((const OctNode*)parent)->edgeNeighbor(4 | (cornerIndex & 1) | ((cornerIndex & 4)>>1) );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==1){				// Agree on edge 0
		const OctNode* temp=((const OctNode*)parent)->edgeNeighbor(((cornerIndex & 2) | (cornerIndex & 4))>>1 );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else{return NULL;}
}
template< class NodeData >
OctNode< NodeData >* OctNode< NodeData >::cornerNeighbor(int cornerIndex,int forceChildren){
	int pIndex,aIndex=0;
	if(!parent){return NULL;}

	pIndex=int(this-parent->children);
	aIndex=(cornerIndex ^ pIndex);	// The disagreement bits
	pIndex=(~pIndex)&7;				// The antipodal point
	if(aIndex==7){					// Agree on no bits
		return &parent->children[pIndex];
	}
	else if(aIndex==0){				// Agree on all bits
		OctNode* temp=((OctNode*)parent)->cornerNeighbor(cornerIndex,forceChildren);
		if(!temp){return NULL;}
		if(!temp->children){
			if(forceChildren){temp->initChildren();}
			else{return temp;}
		}
		return &temp->children[pIndex];
	}
	else if(aIndex==6){				// Agree on face 0
		OctNode* temp=((OctNode*)parent)->__faceNeighbor(0,cornerIndex & 1,0);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==5){				// Agree on face 1
		OctNode* temp=((OctNode*)parent)->__faceNeighbor(1,(cornerIndex & 2)>>1,0);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==3){				// Agree on face 2
		OctNode* temp=((OctNode*)parent)->__faceNeighbor(2,(cornerIndex & 4)>>2,0);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==4){				// Agree on edge 2
		OctNode* temp=((OctNode*)parent)->edgeNeighbor(8 | (cornerIndex & 1) | (cornerIndex & 2) );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==2){				// Agree on edge 1
		OctNode* temp=((OctNode*)parent)->edgeNeighbor(4 | (cornerIndex & 1) | ((cornerIndex & 4)>>1) );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==1){				// Agree on edge 0
		OctNode* temp=((OctNode*)parent)->edgeNeighbor(((cornerIndex & 2) | (cornerIndex & 4))>>1 );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else{return NULL;}
}
////////////////////////
// OctNodeNeighborKey //
////////////////////////
template< class NodeData >
OctNode< NodeData >::Neighbors3::Neighbors3(void){clear();}
template< class NodeData >
void OctNode< NodeData >::Neighbors3::clear(void){
	for(int i=0;i<3;i++){for(int j=0;j<3;j++){for(int k=0;k<3;k++){neighbors[i][j][k]=NULL;}}}
}
template< class NodeData >
OctNode< NodeData >::NeighborKey3::NeighborKey3( void ){ _depth=-1 , neighbors=NULL; }
template< class NodeData >
OctNode< NodeData >::NeighborKey3::NeighborKey3( const NeighborKey3& nKey3 )
{
	_depth = 0 , neighbors = NULL;
	set( nKey3._depth );
	for( int d=0 ; d<=_depth ; d++ ) memcpy( &neighbors[d] , &nKey3.neighbors[d] , sizeof(Neighbors3) );
}
template< class NodeData >
OctNode< NodeData >::NeighborKey3::~NeighborKey3(void)
{
	if( neighbors ) delete[] neighbors;
	neighbors = NULL;
}

template< class NodeData >
void OctNode< NodeData >::NeighborKey3::set( int d )
{
	if( neighbors ) delete[] neighbors;//猜对了，neighbors果然是一个数组
	neighbors = NULL;
	_depth = d;
	if( d<0 ) return;
	neighbors = new Neighbors3[d+1];//声明时要多声明一组，因为初始深度为0。应该是这样，上面对d<0情况进行了判断，没有d==0
}
template< class NodeData >
template< class Real >
bool OctNode< NodeData >::NeighborKey3::setChildNeighbors( Point3D< Real > p , int d , typename OctNode< NodeData >::Neighbors3& childNeighbors ) const
{
	if( !neighbors[d].neighbors[1][1][1] ) return false;
	int i , j , k , x1 , y1 , z1 , x2 , y2 , z2;
	Point3D< Real > c;
	Real w;
	neighbors[d].neighbors[1][1][1]->centerAndWidth( c , w );
	int idx = CornerIndex( c , p );
	Cube::FactorCornerIndex(   idx    , x1 , y1 , z1 );
	Cube::FactorCornerIndex( (~idx)&7 , x2 , y2 , z2 );

	if( !neighbors[d].neighbors[1][1][1]->children ) neighbors[d].neighbors[1][1][1]->initChildren();
	for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
		childNeighbors.neighbors[x2+i][y2+j][z2+k] = &neighbors[d].neighbors[1][1][1]->children[Cube::CornerIndex(i,j,k)];

	// Set the neighbors from across the faces
	i=x1<<1;
	if( neighbors[d].neighbors[i][1][1] )
	{
		if( !neighbors[d].neighbors[i][1][1]->children ) neighbors[d].neighbors[i][1][1]->initChildren();
		for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ ) childNeighbors.neighbors[i][y2+j][z2+k] = &neighbors[d].neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];
	}
	j=y1<<1;
	if( neighbors[d].neighbors[1][j][1] )
	{
		if( !neighbors[d].neighbors[1][j][1]->children ) neighbors[d].neighbors[1][j][1]->initChildren();
		for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ ) childNeighbors.neighbors[x2+i][j][z2+k] = &neighbors[d].neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];
	}
	k=z1<<1;
	if( neighbors[d].neighbors[1][1][k] )
	{
		if( !neighbors[d].neighbors[1][1][k]->children ) neighbors[d].neighbors[1][1][k]->initChildren();
		for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) childNeighbors.neighbors[x2+i][y2+j][k] = &neighbors[d].neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];
	}

	// Set the neighbors from across the edges
	i=x1<<1 , j=y1<<1;
	if( neighbors[d].neighbors[i][j][1] )
	{
		if( !neighbors[d].neighbors[i][j][1]->children ) neighbors[d].neighbors[i][j][1]->initChildren();
		for( k=0 ; k<2 ; k++ ) childNeighbors.neighbors[i][j][z2+k] = &neighbors[d].neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];
	}
	i=x1<<1 , k=z1<<1;
	if( neighbors[d].neighbors[i][1][k] )
	{
		if( !neighbors[d].neighbors[i][1][k]->children ) neighbors[d].neighbors[i][1][k]->initChildren();
		for( j=0 ; j<2 ; j++ ) childNeighbors.neighbors[i][y2+j][k] = &neighbors[d].neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];
	}
	j=y1<<1 , k=z1<<1;
	if( neighbors[d].neighbors[1][j][k] )
	{
		if( !neighbors[d].neighbors[1][j][k]->children ) neighbors[d].neighbors[1][j][k]->initChildren();
		for( i=0 ; i<2 ; i++ ) childNeighbors.neighbors[x2+i][j][k] = &neighbors[d].neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];
	}

	// Set the neighbor from across the corner
	i=x1<<1 , j=y1<<1 , k=z1<<1;
	if( neighbors[d].neighbors[i][j][k] )
	{
		if( !neighbors[d].neighbors[i][j][k]->children ) neighbors[d].neighbors[i][j][k]->initChildren();
		childNeighbors.neighbors[i][j][k] = &neighbors[d].neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
	}
	return true;
}
template< class NodeData >
template< class Real >
bool OctNode< NodeData >::NeighborKey3::getChildNeighbors( Point3D< Real > p , int d , typename OctNode< NodeData >::Neighbors3& childNeighbors ) const
{
	if( !neighbors[d].neighbors[1][1][1] ) return false;
	int i , j , k , x1 , y1 , z1 , x2 , y2 , z2;
	Point3D< Real > c;
	Real w;
	neighbors[d].neighbors[1][1][1]->centerAndWidth( c , w );
	int idx = CornerIndex( c , p );
	Cube::FactorCornerIndex(   idx    , x1 , y1 , z1 );
	Cube::FactorCornerIndex( (~idx)&7 , x2 , y2 , z2 );

	// Set the neighbors of the center cell
	if( neighbors[d].neighbors[1][1][1] && neighbors[d].neighbors[1][1][1]->children )
		for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
			childNeighbors.neighbors[x2+i][y2+j][z2+k] = &neighbors[d].neighbors[1][1][1]->children[Cube::CornerIndex(i,j,k)];

	// Set the neighbors from across the faces
	i=x1<<1;
	if( neighbors[d].neighbors[i][1][1] && neighbors[d].neighbors[i][1][1]->children )
		for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ ) childNeighbors.neighbors[i][y2+j][z2+k] = &neighbors[d].neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];
	else 
		for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ ) childNeighbors.neighbors[i][y2+j][z2+k] = NULL;
	j=y1<<1;
	if( neighbors[d].neighbors[1][j][1] && neighbors[d].neighbors[1][j][1]->children )
		for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ ) childNeighbors.neighbors[x2+i][j][z2+k] = &neighbors[d].neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];
	else
		for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ ) childNeighbors.neighbors[x2+i][j][z2+k] = NULL;
	k=z1<<1;
	if( neighbors[d].neighbors[1][1][k] && neighbors[d].neighbors[1][1][k]->children )
		for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) childNeighbors.neighbors[x2+i][y2+j][k] = &neighbors[d].neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];
	else
		for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) childNeighbors.neighbors[x2+i][y2+j][k] = NULL;

	// Set the neighbors from across the edges
	i=x1<<1 , j=y1<<1;
	if( neighbors[d].neighbors[i][j][1] && neighbors[d].neighbors[i][j][1]->children )
		for( k=0 ; k<2 ; k++ ) childNeighbors.neighbors[i][j][z2+k] = &neighbors[d].neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];
	else
		for( k=0 ; k<2 ; k++ ) childNeighbors.neighbors[i][j][z2+k] = NULL;
	i=x1<<1 , k=z1<<1;
	if( neighbors[d].neighbors[i][1][k] && neighbors[d].neighbors[i][1][k]->children )
		for( j=0 ; j<2 ; j++ ) childNeighbors.neighbors[i][y2+j][k] = &neighbors[d].neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];
	else
		for( j=0 ; j<2 ; j++ ) childNeighbors.neighbors[i][y2+j][k] = NULL;
	j=y1<<1 , k=z1<<1;
	if( neighbors[d].neighbors[1][j][k] && neighbors[d].neighbors[1][j][k]->children )
		for( i=0 ; i<2 ; i++ ) childNeighbors.neighbors[x2+i][j][k] = &neighbors[d].neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];
	else
		for( i=0 ; i<2 ; i++ ) childNeighbors.neighbors[x2+i][j][k] = NULL;

	// Set the neighbor from across the corner
	i=x1<<1 , j=y1<<1 , k=z1<<1;
	if( neighbors[d].neighbors[i][j][k] && neighbors[d].neighbors[i][j][k]->children )
		childNeighbors.neighbors[i][j][k] = &neighbors[d].neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
	else
		childNeighbors.neighbors[i][j][k] = NULL;

	return true;
}
template< class NodeData >
template< class Real >
typename OctNode< NodeData >::Neighbors3& OctNode< NodeData >::NeighborKey3::setNeighbors( OctNode< NodeData >* root , Point3D< Real > p , int d )
{
	//这里设置Neighbor的操作与Neighbor5非常相似，但是不明白为什么要加上一个数据点p作为输入参数
	//每次找的角点都是离p最近的一个，其对应的child node才会被更新Neighbor，那么其他的没有遍历到是不是就不管了
	if( !neighbors[d].neighbors[1][1][1] || !neighbors[d].neighbors[1][1][1]->isInside( p ) )//neighbor[1][1][1]实际上就是node本身吧
	{
		neighbors[d].clear();

		if( !d ) neighbors[d].neighbors[1][1][1] = root;//深度为0时，就设置为root
		else
		{
			Neighbors3& temp = setNeighbors( root , p , d-1 );//还是递归设置

			int i , j , k , x1 , y1 , z1 , x2 , y2 , z2;
			Point3D< Real > c;
			Real w;
			temp.neighbors[1][1][1]->centerAndWidth( c , w );//也就是中心节点的center和width
			int idx = CornerIndex( c , p );//得到距离p最近的角点的index
			Cube::FactorCornerIndex(   idx    , x1 , y1 , z1 );
			Cube::FactorCornerIndex( (~idx)&7 , x2 , y2 , z2 );

			//为这个角点代表的child node设置Neighbor
			if( !temp.neighbors[1][1][1]->children ) temp.neighbors[1][1][1]->initChildren();//细分到child node
			for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
				neighbors[d].neighbors[x2+i][y2+j][z2+k] = &temp.neighbors[1][1][1]->children[Cube::CornerIndex(i,j,k)];

			//下面分别针对面、边、角点重合的node先进行了child node划分，然后找出其中相邻node
			// Set the neighbors from across the faces
			i=x1<<1;//下面这六个node与当前node是面重合
			if( temp.neighbors[i][1][1] )
			{
				if( !temp.neighbors[i][1][1]->children ) temp.neighbors[i][1][1]->initChildren();
				for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][y2+j][z2+k] = &temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];
			}
			j=y1<<1;
			if( temp.neighbors[1][j][1] )
			{
				if( !temp.neighbors[1][j][1]->children ) temp.neighbors[1][j][1]->initChildren();
				for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[x2+i][j][z2+k] = &temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];
			}
			k=z1<<1;
			if( temp.neighbors[1][1][k] )
			{
				if( !temp.neighbors[1][1][k]->children ) temp.neighbors[1][1][k]->initChildren();
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) neighbors[d].neighbors[x2+i][y2+j][k] = &temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];
			}

			// Set the neighbors from across the edges
			i=x1<<1 , j=y1<<1;//下面这些node与中心node是边重合
			if( temp.neighbors[i][j][1] )
			{
				if( !temp.neighbors[i][j][1]->children ) temp.neighbors[i][j][1]->initChildren();
				for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][j][z2+k] = &temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];
			}
			i=x1<<1 , k=z1<<1;
			if( temp.neighbors[i][1][k] )
			{
				if( !temp.neighbors[i][1][k]->children ) temp.neighbors[i][1][k]->initChildren();
				for( j=0 ; j<2 ; j++ ) neighbors[d].neighbors[i][y2+j][k] = &temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];
			}
			j=y1<<1 , k=z1<<1;
			if( temp.neighbors[1][j][k] )
			{
				if( !temp.neighbors[1][j][k]->children ) temp.neighbors[1][j][k]->initChildren();
				for( i=0 ; i<2 ; i++ ) neighbors[d].neighbors[x2+i][j][k] = &temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];
			}

			// Set the neighbor from across the corner
			i=x1<<1 , j=y1<<1 , k=z1<<1;//下面这些node与中心node是角点重合
			if( temp.neighbors[i][j][k] )
			{
				if( !temp.neighbors[i][j][k]->children ) temp.neighbors[i][j][k]->initChildren();
				neighbors[d].neighbors[i][j][k] = &temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
			}
		}
	}
	return neighbors[d];
}
template< class NodeData >
template< class Real >
typename OctNode< NodeData >::Neighbors3& OctNode< NodeData >::NeighborKey3::getNeighbors( OctNode< NodeData >* root , Point3D< Real > p , int d )
{
	//这个get函数除了没有initChild操作，以及多加了一些对child node存在性的判断之外，其他的与set完全相同，不明白其用意???
	if( !neighbors[d].neighbors[1][1][1] || !neighbors[d].neighbors[1][1][1]->isInside( p ) )
	{
		neighbors[d].clear();

		if( !d ) neighbors[d].neighbors[1][1][1] = root;
		else
		{
			Neighbors3& temp = getNeighbors( root , p , d-1 );

			int i , j , k , x1 , y1 , z1 , x2 , y2 , z2;
			Point3D< Real > c;
			Real w;
			temp.neighbors[1][1][1]->centerAndWidth( c , w );
			int idx = CornerIndex( c , p );
			Cube::FactorCornerIndex(   idx    , x1 , y1 , z1 );
			Cube::FactorCornerIndex( (~idx)&7 , x2 , y2 , z2 );

			if( !temp.neighbors[1][1][1] || !temp.neighbors[1][1][1]->children )
			{
				fprintf( stderr , "[ERROR] Couldn't find node at appropriate depth\n" );
				exit( 0 );
			}
			for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
				neighbors[d].neighbors[x2+i][y2+j][z2+k] = &temp.neighbors[1][1][1]->children[Cube::CornerIndex(i,j,k)];


			// Set the neighbors from across the faces
			i=x1<<1;
			if( temp.neighbors[i][1][1] && temp.neighbors[i][1][1]->children )
				for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][y2+j][z2+k] = &temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];
			j=y1<<1;
			if( temp.neighbors[1][j][1] && temp.neighbors[1][j][1]->children )
				for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[x2+i][j][z2+k] = &temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];
			k=z1<<1;
			if( temp.neighbors[1][1][k] && temp.neighbors[1][1][k]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) neighbors[d].neighbors[x2+i][y2+j][k] = &temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];

			// Set the neighbors from across the edges
			i=x1<<1 , j=y1<<1;
			if( temp.neighbors[i][j][1] && temp.neighbors[i][j][1]->children )
				for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][j][z2+k] = &temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];
			i=x1<<1 , k=z1<<1;
			if( temp.neighbors[i][1][k] && temp.neighbors[i][1][k]->children )
				for( j=0 ; j<2 ; j++ ) neighbors[d].neighbors[i][y2+j][k] = &temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];
			j=y1<<1 , k=z1<<1;
			if( temp.neighbors[1][j][k] && temp.neighbors[1][j][k]->children )
				for( i=0 ; i<2 ; i++ ) neighbors[d].neighbors[x2+i][j][k] = &temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];

			// Set the neighbor from across the corner
			i=x1<<1 , j=y1<<1 , k=z1<<1;
			if( temp.neighbors[i][j][k] && temp.neighbors[i][j][k]->children )
				neighbors[d].neighbors[i][j][k] = &temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
		}
	}
	return neighbors[d];
}

template< class NodeData >
typename OctNode< NodeData >::Neighbors3& OctNode< NodeData >::NeighborKey3::setNeighbors( OctNode< NodeData >* node )
{
	//内部也进行迭代，对node的深度之上的所有Neighbor进行了设置，但设定的Neighbor只与当前node相关吧
	int d = node->depth();
	if( node==neighbors[d].neighbors[1][1][1] )
	{
		bool reset = false;
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) if( !neighbors[d].neighbors[i][j][k] ) reset = true;
		if( reset ) neighbors[d].neighbors[1][1][1] = NULL;
	}
	if( node!=neighbors[d].neighbors[1][1][1] )
	{
		neighbors[d].clear();

		if( !node->parent ) neighbors[d].neighbors[1][1][1] = node;
		else
		{
			int i,j,k,x1,y1,z1,x2,y2,z2;
			int idx=int(node-node->parent->children);
			Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
			Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);
			for(i=0;i<2;i++){
				for(j=0;j<2;j++){
					for(k=0;k<2;k++){
						neighbors[d].neighbors[x2+i][y2+j][z2+k]=&node->parent->children[Cube::CornerIndex(i,j,k)];
					}
				}
			}
			Neighbors3& temp=setNeighbors(node->parent);

			// Set the neighbors from across the faces
			i=x1<<1;
			if(temp.neighbors[i][1][1]){
				if(!temp.neighbors[i][1][1]->children) temp.neighbors[i][1][1]->initChildren();
				for(j=0;j<2;j++) for(k=0;k<2;k++) neighbors[d].neighbors[i][y2+j][z2+k]=&temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];
			}
			j=y1<<1;
			if(temp.neighbors[1][j][1]){
				if(!temp.neighbors[1][j][1]->children){temp.neighbors[1][j][1]->initChildren();}
				for(i=0;i<2;i++){for(k=0;k<2;k++){neighbors[d].neighbors[x2+i][j][z2+k]=&temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];}}
			}
			k=z1<<1;
			if(temp.neighbors[1][1][k]){
				if(!temp.neighbors[1][1][k]->children){temp.neighbors[1][1][k]->initChildren();}
				for(i=0;i<2;i++){for(j=0;j<2;j++){neighbors[d].neighbors[x2+i][y2+j][k]=&temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];}}
			}

			// Set the neighbors from across the edges
			i=x1<<1;	j=y1<<1;
			if(temp.neighbors[i][j][1]){
				if(!temp.neighbors[i][j][1]->children){temp.neighbors[i][j][1]->initChildren();}
				for(k=0;k<2;k++){neighbors[d].neighbors[i][j][z2+k]=&temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];}
			}
			i=x1<<1;	k=z1<<1;
			if(temp.neighbors[i][1][k]){
				if(!temp.neighbors[i][1][k]->children){temp.neighbors[i][1][k]->initChildren();}
				for(j=0;j<2;j++){neighbors[d].neighbors[i][y2+j][k]=&temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];}
			}
			j=y1<<1;	k=z1<<1;
			if(temp.neighbors[1][j][k]){
				if(!temp.neighbors[1][j][k]->children){temp.neighbors[1][j][k]->initChildren();}
				for(i=0;i<2;i++){neighbors[d].neighbors[x2+i][j][k]=&temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];}
			}

			// Set the neighbor from across the corner
			i=x1<<1;	j=y1<<1;	k=z1<<1;
			if(temp.neighbors[i][j][k]){
				if(!temp.neighbors[i][j][k]->children){temp.neighbors[i][j][k]->initChildren();}
				neighbors[d].neighbors[i][j][k]=&temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
			}
		}
	}
	return neighbors[d];
}
// Note the assumption is that if you enable an edge, you also enable adjacent faces.
// And, if you enable a corner, you enable adjacent edges and faces.
template< class NodeData >
typename OctNode< NodeData >::Neighbors3& OctNode< NodeData >::NeighborKey3::setNeighbors( OctNode< NodeData >* node , bool flags[3][3][3] )
{
	//与上面的函数类似，只是假如了flags作为标志位
	int d = node->depth();
	if( node==neighbors[d].neighbors[1][1][1] )
	{
		bool reset = false;
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) if( flags[i][j][k] && !neighbors[d].neighbors[i][j][k] ) reset = true;
		if( reset ) neighbors[d].neighbors[1][1][1] = NULL;
	}
	if( node!=neighbors[d].neighbors[1][1][1] )
	{
		neighbors[d].clear();

		if( !node->parent ) neighbors[d].neighbors[1][1][1] = node;
		else
		{
			int x1,y1,z1,x2,y2,z2;
			int idx=int(node-node->parent->children);
			Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
			Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);
			for( int i=0 ; i<2 ; i++ )
				for( int j=0 ; j<2 ; j++ )
					for( int k=0 ; k<2 ; k++ )
						neighbors[d].neighbors[x2+i][y2+j][z2+k]=&node->parent->children[Cube::CornerIndex(i,j,k)];

			Neighbors3& temp=setNeighbors( node->parent , flags );

			// Set the neighbors from across the faces
			{
				int i=x1<<1;
				if( temp.neighbors[i][1][1] )
				{
					if( flags[i][1][1] && !temp.neighbors[i][1][1]->children ) temp.neighbors[i][1][1]->initChildren();
					if( temp.neighbors[i][1][1]->children ) for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][y2+j][z2+k] = &temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];
				}
			}
			{
				int j = y1<<1;
				if( temp.neighbors[1][j][1] )
				{
					if( flags[1][j][1] && !temp.neighbors[1][j][1]->children ) temp.neighbors[1][j][1]->initChildren();
					if( temp.neighbors[1][j][1]->children ) for( int i=0 ; i<2 ; i++ ) for( int k=0 ; k<2 ; k++ ) neighbors[d].neighbors[x2+i][j][z2+k] = &temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];
				}
			}
			{
				int k = z1<<1;
				if( temp.neighbors[1][1][k] )
				{
					if( flags[1][1][k] && !temp.neighbors[1][1][k]->children ) temp.neighbors[1][1][k]->initChildren();
					if( temp.neighbors[1][1][k]->children ) for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) neighbors[d].neighbors[x2+i][y2+j][k] = &temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];
				}
			}

			// Set the neighbors from across the edges
			{
				int i=x1<<1 , j=y1<<1;
				if( temp.neighbors[i][j][1] )
				{
					if( flags[i][j][1] && !temp.neighbors[i][j][1]->children ) temp.neighbors[i][j][1]->initChildren();
					if( temp.neighbors[i][j][1]->children ) for( int k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][j][z2+k] = &temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];
				}
			}
			{
				int i=x1<<1 , k=z1<<1;
				if( temp.neighbors[i][1][k] )
				{
					if( flags[i][1][k] && !temp.neighbors[i][1][k]->children ) temp.neighbors[i][1][k]->initChildren();
					if( temp.neighbors[i][1][k]->children ) for( int j=0 ; j<2 ; j++ ) neighbors[d].neighbors[i][y2+j][k] = &temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];
				}
			}
			{
				int j=y1<<1 , k=z1<<1;
				if( temp.neighbors[1][j][k] )
				{
					if( flags[1][j][k] && !temp.neighbors[1][j][k]->children ) temp.neighbors[1][j][k]->initChildren();
					if( temp.neighbors[1][j][k]->children ) for( int i=0 ; i<2 ; i++ ) neighbors[d].neighbors[x2+i][j][k] = &temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];
				}
			}

			// Set the neighbor from across the corner
			{
				int i=x1<<1 , j=y1<<1 , k=z1<<1;
				if( temp.neighbors[i][j][k] )
				{
					if( flags[i][j][k] && !temp.neighbors[i][j][k]->children ) temp.neighbors[i][j][k]->initChildren();
					if( temp.neighbors[i][j][k]->children ) neighbors[d].neighbors[i][j][k] = &temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
				}
			}
		}
	}
	return neighbors[d];
}

template< class NodeData >
typename OctNode< NodeData >::Neighbors3& OctNode< NodeData >::NeighborKey3::getNeighbors( OctNode< NodeData >* node )
{
	//相较于setNeighbor函数，也是多加入了对childnode的判断
	int d=node->depth();
	if(node!=neighbors[d].neighbors[1][1][1])
	{
		neighbors[d].clear();

		if( !node->parent ) neighbors[d].neighbors[1][1][1] = node;
		else
		{
			int i,j,k,x1,y1,z1,x2,y2,z2;
			int idx=int(node-node->parent->children);
			Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
			Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);
			for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
				neighbors[d].neighbors[x2+i][y2+j][z2+k] = node->parent->children + Cube::CornerIndex(i,j,k);

			Neighbors3& temp=getNeighbors(node->parent);

			// Set the neighbors from across the faces
			i=x1<<1;
			if( temp.neighbors[i][1][1] && temp.neighbors[i][1][1]->children )
				for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][y2+j][z2+k] = temp.neighbors[i][1][1]->children + Cube::CornerIndex(x2,j,k);
			j=y1<<1;
			if( temp.neighbors[1][j][1] && temp.neighbors[1][j][1]->children )
				for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[x2+i][j][z2+k] = temp.neighbors[1][j][1]->children + Cube::CornerIndex(i,y2,k);
			k=z1<<1;
			if( temp.neighbors[1][1][k] && temp.neighbors[1][1][k]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++) neighbors[d].neighbors[x2+i][y2+j][k] = temp.neighbors[1][1][k]->children + Cube::CornerIndex(i,j,z2);

			// Set the neighbors from across the edges
			i=x1<<1;	j=y1<<1;
			if( temp.neighbors[i][j][1] && temp.neighbors[i][j][1]->children )
				for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][j][z2+k] = temp.neighbors[i][j][1]->children + Cube::CornerIndex(x2,y2,k);
			i=x1<<1;	k=z1<<1;
			if( temp.neighbors[i][1][k] && temp.neighbors[i][1][k]->children )
				for( j=0 ; j<2 ; j++ ) neighbors[d].neighbors[i][y2+j][k] = temp.neighbors[i][1][k]->children + Cube::CornerIndex(x2,j,z2);
			j=y1<<1;	k=z1<<1;
			if( temp.neighbors[1][j][k] && temp.neighbors[1][j][k]->children )
				for( i=0 ; i<2 ; i++ ) neighbors[d].neighbors[x2+i][j][k] = temp.neighbors[1][j][k]->children + Cube::CornerIndex(i,y2,z2);

			// Set the neighbor from across the corner
			i=x1<<1;	j=y1<<1;	k=z1<<1;
			if( temp.neighbors[i][j][k] && temp.neighbors[i][j][k]->children )
				neighbors[d].neighbors[i][j][k] = temp.neighbors[i][j][k]->children + Cube::CornerIndex(x2,y2,z2);
		}
	}
	return neighbors[node->depth()];
}
template< class NodeData >
void OctNode< NodeData >::NeighborKey3::setNeighbors( OctNode< NodeData >* node , typename OctNode< NodeData >::Neighbors5& neighbors )
{
	//利用一环邻域Neighbor来设置二环邻域Neighbor???
	neighbors.clear();
	if( !node ) return;
	if( !node->parent ) neighbors.neighbors[2][2][2] = node;
	else
	{
		int c = int( node - node->parent->children );
		const OctNode< NodeData >::Neighbors3& _neighbors = setNeighbors( node->parent );
		switch( c )
		{
		case 0:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] )
				{
					if( !_neighbors.neighbors[i][j][k]->children ) _neighbors.neighbors[i][j][k]->initChildren();
					for( int ii=0 ; ii<( (i==2) ? 1 : 2 ) ; ii++ )
						for( int jj=0 ; jj<( (j==2) ? 1 : 2 ) ; jj++ )
							for( int kk=0 ; kk<( (k==2) ? 1 : 2 ) ; kk++ )
								neighbors.neighbors[2*i+ii][2*j+jj][2*k+kk] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
				}
			break;
		case 1:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] )
				{
					if( !_neighbors.neighbors[i][j][k]->children ) _neighbors.neighbors[i][j][k]->initChildren();
					for( int ii=( (i==0) ? 1 : 0 ) ; ii<2 ; ii++ )
						for( int jj=0 ; jj<( (j==2) ? 1 : 2 ) ; jj++ )
							for( int kk=0 ; kk<( (k==2) ? 1 : 2 ) ; kk++ )
								neighbors.neighbors[2*i+ii-1][2*j+jj][2*k+kk] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
				}
			break;
		case 2:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] )
				{
					if( !_neighbors.neighbors[i][j][k]->children ) _neighbors.neighbors[i][j][k]->initChildren();
					for( int ii=0 ; ii<( (i==2) ? 1 : 2 ) ; ii++ )
						for( int jj=( (j==0) ? 1 : 0 ) ; jj<2 ; jj++ )
							for( int kk=0 ; kk<( (k==2) ? 1 : 2 ) ; kk++ )
								neighbors.neighbors[2*i+ii][2*j+jj-1][2*k+kk] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
				}
			break;
		case 3:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] )
				{
					if( !_neighbors.neighbors[i][j][k]->children ) _neighbors.neighbors[i][j][k]->initChildren();
					for( int ii=( (i==0) ? 1 : 0 ) ; ii<2 ; ii++ )
						for( int jj=( (j==0) ? 1 : 0 ) ; jj<2 ; jj++ )
							for( int kk=0 ; kk<( (k==2) ? 1 : 2 ) ; kk++ )
								neighbors.neighbors[2*i+ii-1][2*j+jj-1][2*k+kk] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
				}
			break;
		case 4:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] )
				{
					if( !_neighbors.neighbors[i][j][k]->children ) _neighbors.neighbors[i][j][k]->initChildren();
					for( int ii=0 ; ii<( (i==2) ? 1 : 2 ) ; ii++ )
						for( int jj=0 ; jj<( (j==2) ? 1 : 2 ) ; jj++ )
							for( int kk=( (k==0) ? 1 : 0 ) ; kk<2 ; kk++ )
								neighbors.neighbors[2*i+ii][2*j+jj][2*k+kk-1] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
				}
			break;
		case 5:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] )
				{
					if( !_neighbors.neighbors[i][j][k]->children ) _neighbors.neighbors[i][j][k]->initChildren();
					for( int ii=( (i==0) ? 1 : 0 ) ; ii<2 ; ii++ )
						for( int jj=0 ; jj<( (j==2) ? 1 : 2 ) ; jj++ )
							for( int kk=( (k==0) ? 1 : 0 ) ; kk<2 ; kk++ )
								neighbors.neighbors[2*i+ii-1][2*j+jj][2*k+kk-1] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
				}
			break;
		case 6:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] )
				{
					if( !_neighbors.neighbors[i][j][k]->children ) _neighbors.neighbors[i][j][k]->initChildren();
					for( int ii=0 ; ii<( (i==2) ? 1 : 2 ) ; ii++ )
						for( int jj=( (j==0) ? 1 : 0 ) ; jj<2 ; jj++ )
							for( int kk=( (k==0) ? 1 : 0 ) ; kk<2 ; kk++ )
								neighbors.neighbors[2*i+ii][2*j+jj-1][2*k+kk-1] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
				}
			break;
		case 7:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] )
				{
					if( !_neighbors.neighbors[i][j][k]->children ) _neighbors.neighbors[i][j][k]->initChildren();
					for( int ii=( (i==0) ? 1 : 0 ) ; ii<2 ; ii++ )
						for( int jj=( (j==0) ? 1 : 0 ) ; jj<2 ; jj++ )
							for( int kk=( (k==0) ? 1 : 0 ) ; kk<2 ; kk++ )
								neighbors.neighbors[2*i+ii-1][2*j+jj-1][2*k+kk-1] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
				}
			break;
		}
	}
}
template< class NodeData >
void OctNode< NodeData >::NeighborKey3::getNeighbors( OctNode< NodeData >* node , typename OctNode< NodeData >::Neighbors5& neighbors )
{
	neighbors.clear();
	if( !node ) return;
	if( !node->parent ) neighbors.neighbors[2][2][2] = node;
	else
	{
		int c = int( node - node->parent->children );
		const OctNode< NodeData >::Neighbors3& _neighbors = getNeighbors( node->parent );
		OctNode< NodeData >* const * _nodes = &_neighbors.neighbors[0][0][0];
		const OctNode< NodeData >* const * _node;
		const OctNode< NodeData >* __node;
		int iS , iE , jS , jE , kS , kE;
#define _S( i ) ( (i==0) ? 1 : 0 )
#define _E( i ) ( (i==2) ? 1 : 2 )
		switch( c )
		{
		case 0:
			_node = _nodes;
			for( int i=0 ; i<3 ; i++ ){ iE=_E(i) ; for( int j=0 ; j<3 ; j++ ){ jE=_E(j) ; for( int k=0 ; k<3 ; k++ , _node++ ){ kE=_E(k) , __node = *_node;
			if( __node && __node->children ) for( int ii=0 , iii=2*i ; ii<iE ; ii++ , iii++ ) for( int jj=0 , jjj=2*j ; jj<jE ; jj++ , jjj++ ) for( int kk=0 , kkk=2*k ; kk<kE ; kk++ , kkk++ )
				neighbors.neighbors[iii][jjj][kkk] = __node->children + Cube::CornerIndex( ii , jj , kk );
			} } }
			break;
		case 1:
			_node = _nodes;
			for( int i=0 ; i<3 ; i++ ){ iS=_S(i) ; for( int j=0 ; j<3 ; j++ ){ jE=_E(j) ; for( int k=0 ; k<3 ; k++ , _node++ ){ kE=_E(k) , __node = *_node;
			if( __node && __node->children ) for( int ii=iS , iii=2*i+iS-1 ; ii<2 ; ii++ , iii++ ) for( int jj=0 , jjj=2*j ; jj<jE ; jj++ , jjj++ ) for( int kk=0 , kkk=2*k ; kk<kE ; kk++ , kkk++ )
				neighbors.neighbors[iii][jjj][kkk] = __node->children + Cube::CornerIndex( ii , jj , kk );
			} } }
			break;
		case 2:
			_node = _nodes;
			for( int i=0 ; i<3 ; i++ ){ iE=_E(i) ; for( int j=0 ; j<3 ; j++ ){ jS=_S(j) ; for( int k=0 ; k<3 ; k++ , _node++ ){ kE=_E(k) , __node = *_node;
			if( __node && __node->children ) for( int ii=0 , iii=2*i ; ii<iE ; ii++ , iii++ ) for( int jj=jS , jjj=2*j+jS-1 ; jj<2 ; jj++ , jjj++ ) for( int kk=0 , kkk=2*k ; kk<kE; kk++ , kkk++ )
				neighbors.neighbors[iii][jjj][kkk] = __node->children + Cube::CornerIndex( ii , jj , kk );
			} } }
			break;
		case 3:
			_node = _nodes;
			for( int i=0  ; i<3 ; i++ ){ iS=_S(i) ; for( int j=0 ; j<3 ; j++ ){ jS=_S(j) ; for( int k=0 ; k<3 ; k++ , _node++ ){ kE=_E(k) , __node = *_node;
			if( __node && __node->children ) for( int ii=iS , iii=2*i+iS-1 ; ii<2 ; ii++ , iii++ ) for( int jj=jS , jjj=2*j+jS-1 ; jj<2 ; jj++ , jjj++ ) for( int kk=0 , kkk=2*k ; kk<kE ; kk++ , kkk++ )
				neighbors.neighbors[iii][jjj][kkk] = __node->children + Cube::CornerIndex( ii , jj , kk );
			} } }
			break;
		case 4:
			_node = _nodes;
			for( int i=0 ; i<3 ; i++ ){ iE=_E(i) ; for( int j=0 ; j<3 ; j++ ){ jE=_E(j) ; for( int k=0 ; k<3 ; k++ , _node++ ){ kS=_S(k) , __node = *_node;
			if( __node && __node->children ) for( int ii=0 , iii=2*i ; ii<iE ; ii++ , iii++ ) for( int jj=0 , jjj=2*j ; jj<jE ; jj++ , jjj++ ) for( int kk=kS , kkk=2*k+kS-1 ; kk<2 ; kk++ , kkk++ )
				neighbors.neighbors[iii][jjj][kkk] = __node->children + Cube::CornerIndex( ii , jj , kk );
			} } }
			break;
		case 5:
			_node = _nodes;
			for( int i=0 ; i<3 ; i++ ){ iS=_S(i) ; for( int j=0 ; j<3 ; j++ ){ jE=_E(j) ; for( int k=0 ; k<3 ; k++ , _node++ ){ kS=_S(k) , __node = *_node;
			if( __node && __node->children ) for( int ii=iS , iii=2*i+iS-1 ; ii<2 ; ii++ , iii++ ) for( int jj=0 , jjj=2*j ; jj<jE ; jj++ , jjj++ ) for( int kk=kS , kkk=2*k+kS-1 ; kk<2 ; kk++ , kkk++ )
				neighbors.neighbors[iii][jjj][kkk] = __node->children + Cube::CornerIndex( ii , jj , kk );
			} } }
			break;
		case 6:
			_node = _nodes;
			for( int i=0 ; i<3 ; i++ ){ iE=_E(i) ; for( int j=0 ; j<3 ; j++ ){ jS=_S(j) ; for( int k=0 ; k<3 ; k++ , _node++ ){ kS=_S(k) , __node = *_node;
			if( __node && __node->children ) for( int ii=0 , iii=2*i ; ii<iE ; ii++ , iii++ ) for( int jj=jS , jjj=2*j+jS-1 ; jj<2 ; jj++ , jjj++ ) for( int kk=kS , kkk=2*k+kS-1 ; kk<2 ; kk++ , kkk++ )
				neighbors.neighbors[iii][jjj][kkk] = __node->children + Cube::CornerIndex( ii , jj , kk );
			} } }
			break;
		case 7:
			_node = _nodes;
			for( int i=0 ; i<3 ; i++ ){ iS=_S(i) ; for( int j=0 ; j<3 ; j++ ){ jS=_S(j) ; for( int k=0 ; k<3 ; k++ , _node++ ){ kS=_S(k) , __node = *_node;
			if( __node && __node->children ) for( int ii=iS , iii=2*i+iS-1 ; ii<2 ; ii++ , iii++ ) for( int jj=jS , jjj=2*j+jS-1 ; jj<2 ; jj++ , jjj++ ) for( int kk=kS , kkk=2*k+kS-1 ; kk<2 ; kk++ , kkk++ )
				neighbors.neighbors[iii][jjj][kkk] = __node->children + Cube::CornerIndex( ii , jj , kk );
			} } }
			break;
		}
#undef _S
#undef _E
	}
}
template< class NodeData >
void OctNode< NodeData >::ConstNeighborKey3::getNeighbors( const OctNode< NodeData >* node , typename OctNode< NodeData >::ConstNeighbors5& neighbors )
{
	neighbors.clear();
	if( !node ) return;
	if( !node->parent ) neighbors.neighbors[2][2][2] = node;
	else
	{
		int c = int( node - node->parent->children );
		const OctNode< NodeData >::ConstNeighbors3& _neighbors = getNeighbors( node->parent );
		switch( c )
		{
		case 0:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] && _neighbors.neighbors[i][j][k]->children )
					for( int ii=0 ; ii<( (i==2) ? 1 : 2 ) ; ii++ )
						for( int jj=0 ; jj<( (j==2) ? 1 : 2 ) ; jj++ )
							for( int kk=0 ; kk<( (k==2) ? 1 : 2 ) ; kk++ )
								neighbors.neighbors[2*i+ii][2*j+jj][2*k+kk] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
			break;
		case 1:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] && _neighbors.neighbors[i][j][k]->children )
					for( int ii=( (i==0) ? 1 : 0 ) ; ii<2 ; ii++ )
						for( int jj=0 ; jj<( (j==2) ? 1 : 2 ) ; jj++ )
							for( int kk=0 ; kk<( (k==2) ? 1 : 2 ) ; kk++ )
								neighbors.neighbors[2*i+ii-1][2*j+jj][2*k+kk] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
			break;
		case 2:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] && _neighbors.neighbors[i][j][k]->children )
					for( int ii=0 ; ii<( (i==2) ? 1 : 2 ) ; ii++ )
						for( int jj=( (j==0) ? 1 : 0 ) ; jj<2 ; jj++ )
							for( int kk=0 ; kk<( (k==2) ? 1 : 2 ) ; kk++ )
								neighbors.neighbors[2*i+ii][2*j+jj-1][2*k+kk] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
			break;
		case 3:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] && _neighbors.neighbors[i][j][k]->children )
					for( int ii=( (i==0) ? 1 : 0 ) ; ii<2 ; ii++ )
						for( int jj=( (j==0) ? 1 : 0 ) ; jj<2 ; jj++ )
							for( int kk=0 ; kk<( (k==2) ? 1 : 2 ) ; kk++ )
								neighbors.neighbors[2*i+ii-1][2*j+jj-1][2*k+kk] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
			break;
		case 4:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] && _neighbors.neighbors[i][j][k]->children )
					for( int ii=0 ; ii<( (i==2) ? 1 : 2 ) ; ii++ )
						for( int jj=0 ; jj<( (j==2) ? 1 : 2 ) ; jj++ )
							for( int kk=( (k==0) ? 1 : 0 ) ; kk<2 ; kk++ )
								neighbors.neighbors[2*i+ii][2*j+jj][2*k+kk-1] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
			break;
		case 5:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] && _neighbors.neighbors[i][j][k]->children )
					for( int ii=( (i==0) ? 1 : 0 ) ; ii<2 ; ii++ )
						for( int jj=0 ; jj<( (j==2) ? 1 : 2 ) ; jj++ )
							for( int kk=( (k==0) ? 1 : 0 ) ; kk<2 ; kk++ )
								neighbors.neighbors[2*i+ii-1][2*j+jj][2*k+kk-1] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
			break;
		case 6:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] && _neighbors.neighbors[i][j][k]->children )
					for( int ii=0 ; ii<( (i==2) ? 1 : 2 ) ; ii++ )
						for( int jj=( (j==0) ? 1 : 0 ) ; jj<2 ; jj++ )
							for( int kk=( (k==0) ? 1 : 0 ) ; kk<2 ; kk++ )
								neighbors.neighbors[2*i+ii][2*j+jj-1][2*k+kk-1] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
			break;
		case 7:
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( _neighbors.neighbors[i][j][k] && _neighbors.neighbors[i][j][k]->children )
					for( int ii=( (i==0) ? 1 : 0 ) ; ii<2 ; ii++ )
						for( int jj=( (j==0) ? 1 : 0 ) ; jj<2 ; jj++ )
							for( int kk=( (k==0) ? 1 : 0 ) ; kk<2 ; kk++ )
								neighbors.neighbors[2*i+ii-1][2*j+jj-1][2*k+kk-1] = _neighbors.neighbors[i][j][k]->children + Cube::CornerIndex( ii , jj , kk );
			break;
		}
	}
}

///////////////////////
// ConstNeighborKey3 //
///////////////////////
template< class NodeData >
OctNode< NodeData >::ConstNeighbors3::ConstNeighbors3(void){clear();}
template< class NodeData >
void OctNode< NodeData >::ConstNeighbors3::clear(void){
	for(int i=0;i<3;i++){for(int j=0;j<3;j++){for(int k=0;k<3;k++){neighbors[i][j][k]=NULL;}}}
}
template< class NodeData >
OctNode< NodeData >::ConstNeighborKey3::ConstNeighborKey3( void ){ _depth=-1 , neighbors=NULL; }
template< class NodeData >
OctNode< NodeData >::ConstNeighborKey3::ConstNeighborKey3( const ConstNeighborKey3& key3 )
{
	_depth = 0 , neighbors = NULL;
	set( key3._depth );
	for( int d=0 ; d<=_depth ; d++ ) memcpy( &neighbors[d] , &key3.neighbors[d] , sizeof(ConstNeighbors3) );
}
template< class NodeData >
OctNode< NodeData >::ConstNeighborKey3::~ConstNeighborKey3(void){
	if( neighbors ) delete[] neighbors;
	neighbors=NULL;
}

template< class NodeData >
void OctNode< NodeData >::ConstNeighborKey3::set( int d )
{
	if( neighbors ) delete[] neighbors;
	neighbors = NULL;
	_depth = d;
	if( d<0 ) return;
	neighbors = new ConstNeighbors3[d+1];
}
template< class NodeData >
typename OctNode< NodeData >::ConstNeighbors3& OctNode< NodeData >::ConstNeighborKey3::getNeighbors(const OctNode< NodeData >* node)
{
	int d=node->depth();
	if( node!=neighbors[d].neighbors[1][1][1] )
	{
		neighbors[d].clear();

		if(!node->parent) neighbors[d].neighbors[1][1][1]=node;
		else
		{
			int i,j,k,x1,y1,z1,x2,y2,z2;
			int idx=int(node-node->parent->children);
			Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
			Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);
			for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
				neighbors[d].neighbors[x2+i][y2+j][z2+k] = node->parent->children + Cube::CornerIndex(i,j,k);

			ConstNeighbors3& temp=getNeighbors(node->parent);

			// Set the neighbors from across the faces
			i=x1<<1;
			if( temp.neighbors[i][1][1] && temp.neighbors[i][1][1]->children )
				for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][y2+j][z2+k] = temp.neighbors[i][1][1]->children + Cube::CornerIndex(x2,j,k);
			j=y1<<1;
			if( temp.neighbors[1][j][1] && temp.neighbors[1][j][1]->children )
				for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[x2+i][j][z2+k] = temp.neighbors[1][j][1]->children + Cube::CornerIndex(i,y2,k);
			k=z1<<1;
			if( temp.neighbors[1][1][k] && temp.neighbors[1][1][k]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) neighbors[d].neighbors[x2+i][y2+j][k] = temp.neighbors[1][1][k]->children + Cube::CornerIndex(i,j,z2);

			// Set the neighbors from across the edges
			i=x1<<1;	j=y1<<1;
			if( temp.neighbors[i][j][1] && temp.neighbors[i][j][1]->children )
				for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][j][z2+k] = temp.neighbors[i][j][1]->children + Cube::CornerIndex(x2,y2,k);
			i=x1<<1;	k=z1<<1;
			if( temp.neighbors[i][1][k] && temp.neighbors[i][1][k]->children )
				for( j=0 ; j<2 ; j++ ) neighbors[d].neighbors[i][y2+j][k] = temp.neighbors[i][1][k]->children + Cube::CornerIndex(x2,j,z2);
			j=y1<<1;	k=z1<<1;
			if( temp.neighbors[1][j][k] && temp.neighbors[1][j][k]->children )
				for( i=0 ; i<2 ; i++ ) neighbors[d].neighbors[x2+i][j][k] = temp.neighbors[1][j][k]->children + Cube::CornerIndex(i,y2,z2);

			// Set the neighbor from across the corner
			i=x1<<1;	j=y1<<1;	k=z1<<1;
			if( temp.neighbors[i][j][k] && temp.neighbors[i][j][k]->children )
				neighbors[d].neighbors[i][j][k] = temp.neighbors[i][j][k]->children + Cube::CornerIndex(x2,y2,z2);
		}
	}
	return neighbors[node->depth()];
}
template< class NodeData >
typename OctNode< NodeData >::ConstNeighbors3& OctNode< NodeData >::ConstNeighborKey3::getNeighbors( const OctNode< NodeData >* node , int minDepth )
{
	int d=node->depth();
	if( d<minDepth ) fprintf( stderr , "[ERROR] Node depth lower than min-depth: %d < %d\n" , d , minDepth ) , exit( 0 );
	if( node!=neighbors[d].neighbors[1][1][1] )
	{
		neighbors[d].clear();

		if( d==minDepth ) neighbors[d].neighbors[1][1][1]=node;
		else
		{
			int i,j,k,x1,y1,z1,x2,y2,z2;
			int idx = int(node-node->parent->children);
			Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
			Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);

			ConstNeighbors3& temp=getNeighbors( node->parent , minDepth );

			//  Set the syblings
			for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
				neighbors[d].neighbors[x2+i][y2+j][z2+k] = node->parent->children + Cube::CornerIndex(i,j,k);

			// Set the neighbors from across the faces
			i=x1<<1;
			if( temp.neighbors[i][1][1] && temp.neighbors[i][1][1]->children )
				for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][y2+j][z2+k] = temp.neighbors[i][1][1]->children + Cube::CornerIndex(x2,j,k);

			j=y1<<1;
			if( temp.neighbors[1][j][1] && temp.neighbors[1][j][1]->children )
				for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[x2+i][j][z2+k] = temp.neighbors[1][j][1]->children + Cube::CornerIndex(i,y2,k);

			k=z1<<1;
			if( temp.neighbors[1][1][k] && temp.neighbors[1][1][k]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) neighbors[d].neighbors[x2+i][y2+j][k] = temp.neighbors[1][1][k]->children + Cube::CornerIndex(i,j,z2);

			// Set the neighbors from across the edges
			i=x1<<1 , j=y1<<1;
			if( temp.neighbors[i][j][1] && temp.neighbors[i][j][1]->children )
				for( k=0 ; k<2 ; k++ ) neighbors[d].neighbors[i][j][z2+k] = temp.neighbors[i][j][1]->children + Cube::CornerIndex(x2,y2,k);

			i=x1<<1 , k=z1<<1;
			if( temp.neighbors[i][1][k] && temp.neighbors[i][1][k]->children )
				for( j=0 ; j<2 ; j++ ) neighbors[d].neighbors[i][y2+j][k] = temp.neighbors[i][1][k]->children + Cube::CornerIndex(x2,j,z2);

			j=y1<<1 , k=z1<<1;
			if( temp.neighbors[1][j][k] && temp.neighbors[1][j][k]->children )
				for( i=0 ; i<2 ; i++ ) neighbors[d].neighbors[x2+i][j][k] = temp.neighbors[1][j][k]->children + Cube::CornerIndex(i,y2,z2);

			// Set the neighbor from across the corner
			i=x1<<1 , j=y1<<1 , k=z1<<1;
			if( temp.neighbors[i][j][k] && temp.neighbors[i][j][k]->children )
				neighbors[d].neighbors[i][j][k] = temp.neighbors[i][j][k]->children + Cube::CornerIndex(x2,y2,z2);
		}
	}
	return neighbors[node->depth()];
}

template< class NodeData > OctNode< NodeData >::Neighbors5::Neighbors5( void ){ clear(); }
template< class NodeData > OctNode< NodeData >::ConstNeighbors5::ConstNeighbors5( void ){ clear(); }
template< class NodeData >
void OctNode< NodeData >::Neighbors5::clear( void )
{
	for( int i=0 ; i<5 ; i++ ) for( int j=0 ; j<5 ; j++ ) for( int k=0 ; k<5 ; k++ ) neighbors[i][j][k] = NULL;
}
template< class NodeData >
void OctNode< NodeData >::ConstNeighbors5::clear( void )
{
	for( int i=0 ; i<5 ; i++ ) for( int j=0 ; j<5 ; j++ ) for( int k=0 ; k<5 ; k++ ) neighbors[i][j][k] = NULL;
}
template< class NodeData >
OctNode< NodeData >::NeighborKey5::NeighborKey5( void )
{
	_depth = -1;
	neighbors = NULL;
}
template< class NodeData >
OctNode< NodeData >::ConstNeighborKey5::ConstNeighborKey5( void )
{
	_depth = -1;
	neighbors = NULL;
}
template< class NodeData >
OctNode< NodeData >::NeighborKey5::~NeighborKey5( void )
{
	if( neighbors ) delete[] neighbors;
	neighbors = NULL;
}
template< class NodeData >
OctNode< NodeData >::ConstNeighborKey5::~ConstNeighborKey5( void )
{
	if( neighbors ) delete[] neighbors;
	neighbors = NULL;
}
template< class NodeData >
void OctNode< NodeData >::NeighborKey5::set( int d )
{
	if( neighbors ) delete[] neighbors;
	neighbors = NULL;
	if(d<0) return;
	_depth = d;
	neighbors=new Neighbors5[d+1];
}
template< class NodeData >
void OctNode< NodeData >::ConstNeighborKey5::set( int d )
{
	if( neighbors ) delete[] neighbors;
	neighbors = NULL;
	if(d<0) return;
	_depth = d;
	neighbors=new ConstNeighbors5[d+1];
}
template< class NodeData >
typename OctNode< NodeData >::Neighbors5& OctNode< NodeData >::NeighborKey5::getNeighbors( OctNode* node )
{
	int d=node->depth();
	if( node!=neighbors[d].neighbors[2][2][2] )
	{
		//这算是getNeighbors?如果要重新设置neighbors的值，应该算是setNeighbors
		neighbors[d].clear();

		if( !node->parent ) neighbors[d].neighbors[2][2][2]=node;
		else
		{
			getNeighbors( node->parent );
			Neighbors5& temp = neighbors[d-1];//父节点的Neighbor
			int x1 , y1 , z1 , x2 , y2 , z2;
			int idx = int( node - node->parent->children );
			Cube::FactorCornerIndex( idx , x1 , y1 , z1 );//找出当前node在parent node的child中的具体排位

			Neighbors5& n = neighbors[d];//当前节点的Neighbor
			Cube::FactorCornerIndex( (~idx)&7 , x2 , y2 , z2 );//(~idx)&7效果应该等同于(7-idx)
			int i , j , k;
			int fx0 = x2+1 , fy0 = y2+1 , fz0 = z2+1;	// Indices of the bottom left corner of the parent within the 5x5x5
			int cx1 = x1*2+1 , cy1 = y1*2+1 , cz1 = z1*2+1;
			int cx2 = x2*2+1 , cy2 = y2*2+1 , cz2 = z2*2+1;
			int fx1 = x1*3 , fy1 = y1*3 , fz1 = z1*3;
			int fx2 = x2*4 , fy2 = y2*4 , fz2 = z2*4;

			//  Set the syblings
			for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
				n.neighbors[fx0+i][fy0+j][fz0+k] = node->parent->children + Cube::CornerIndex( i , j , k );

			//因为当前只设置二环邻域的Neighbor，因此只需要在父节点一环邻域的子节点中查找就可以满足要求
			// Set the neighbors from across the faces
			//across face意味着在parent层面，两个parent有一个相邻面，六个
			//简单验证了一下，应该是对的
			if( temp.neighbors[cx1][2][2] && temp.neighbors[cx1][2][2]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx1+i][fy0+j][fz0+k] = temp.neighbors[cx1][2][2]->children + Cube::CornerIndex( i , j , k );
			if( temp.neighbors[2][cy1][2] && temp.neighbors[2][cy1][2]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx0+i][fy1+j][fz0+k] = temp.neighbors[2][cy1][2]->children + Cube::CornerIndex( i , j , k );
			if( temp.neighbors[2][2][cz1] && temp.neighbors[2][2][cz1]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx0+i][fy0+j][fz1+k] = temp.neighbors[2][2][cz1]->children + Cube::CornerIndex( i , j , k );
			if( temp.neighbors[cx2][2][2] && temp.neighbors[cx2][2][2]->children )
				for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx2  ][fy0+j][fz0+k] = temp.neighbors[cx2][2][2]->children + Cube::CornerIndex( x1 , j , k );
			if( temp.neighbors[2][cy2][2] && temp.neighbors[2][cy2][2]->children )
				for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx0+i][fy2  ][fz0+k] = temp.neighbors[2][cy2][2]->children + Cube::CornerIndex( i , y1 , k );
			if( temp.neighbors[2][2][cz2] && temp.neighbors[2][2][cz2]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ )
					n.neighbors[fx0+i][fy0+j][fz2  ] = temp.neighbors[2][2][cz2]->children + Cube::CornerIndex( i , j , z1 );

			// Set the neighbors from across the edges
			//parent之间只有相邻边，12个???
			if( temp.neighbors[cx1][cy1][2] && temp.neighbors[cx1][cy1][2]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx1+i][fy1+j][fz0+k] = temp.neighbors[cx1][cy1][2]->children + Cube::CornerIndex( i , j , k );
			if( temp.neighbors[cx1][2][cz1] && temp.neighbors[cx1][2][cz1]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx1+i][fy0+j][fz1+k] = temp.neighbors[cx1][2][cz1]->children + Cube::CornerIndex( i , j , k );
			if( temp.neighbors[2][cy1][cz1] && temp.neighbors[2][cy1][cz1]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx0+i][fy1+j][fz1+k] = temp.neighbors[2][cy1][cz1]->children + Cube::CornerIndex( i , j , k );
			if( temp.neighbors[cx1][cy2][2] && temp.neighbors[cx1][cy2][2]->children )
				for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx1+i][fy2  ][fz0+k] = temp.neighbors[cx1][cy2][2]->children + Cube::CornerIndex( i , y1 , k );
			if( temp.neighbors[cx1][2][cz2] && temp.neighbors[cx1][2][cz2]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ )
					n.neighbors[fx1+i][fy0+j][fz2  ] = temp.neighbors[cx1][2][cz2]->children + Cube::CornerIndex( i , j , z1 );
			if( temp.neighbors[cx2][cy1][2] && temp.neighbors[cx2][cy1][2]->children )
				for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx2  ][fy1+j][fz0+k] = temp.neighbors[cx2][cy1][2]->children + Cube::CornerIndex( x1 , j , k );
			if( temp.neighbors[2][cy1][cz2] && temp.neighbors[2][cy1][cz2]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ )
					n.neighbors[fx0+i][fy1+j][fz2  ] = temp.neighbors[2][cy1][cz2]->children + Cube::CornerIndex( i , j , z1 );
			if( temp.neighbors[cx2][2][cz1] && temp.neighbors[cx2][2][cz1]->children )
				for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx2  ][fy0+j][fz1+k] = temp.neighbors[cx2][2][cz1]->children + Cube::CornerIndex( x1 , j , k );
			if( temp.neighbors[2][cy2][cz1] && temp.neighbors[2][cy2][cz1]->children )
				for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx0+i][fy2  ][fz1+k] = temp.neighbors[2][cy2][cz1]->children + Cube::CornerIndex( i , y1 , k );
			if( temp.neighbors[cx2][cy2][2] && temp.neighbors[cx2][cy2][2]->children )
				for( k=0 ; k<2 ; k++ )
					n.neighbors[fx2  ][fy2  ][fz0+k] = temp.neighbors[cx2][cy2][2]->children + Cube::CornerIndex( x1 , y1 , k );
			if( temp.neighbors[cx2][2][cz2] && temp.neighbors[cx2][2][cz2]->children )
				for( j=0 ; j<2 ; j++ )
					n.neighbors[fx2  ][fy0+j][fz2  ] = temp.neighbors[cx2][2][cz2]->children + Cube::CornerIndex( x1 , j , z1 );
			if( temp.neighbors[2][cy2][cz2] && temp.neighbors[2][cy2][cz2]->children )
				for( i=0 ; i<2 ; i++ )
					n.neighbors[fx0+i][fy2  ][fz2  ] = temp.neighbors[2][cy2][cz2]->children + Cube::CornerIndex( i , y1 , z1 );

			// Set the neighbor from across the corners
			//parent之间只有相邻点，8个
			if( temp.neighbors[cx1][cy1][cz1] && temp.neighbors[cx1][cy1][cz1]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx1+i][fy1+j][fz1+k] = temp.neighbors[cx1][cy1][cz1]->children + Cube::CornerIndex( i , j , k );
			if( temp.neighbors[cx1][cy1][cz2] && temp.neighbors[cx1][cy1][cz2]->children )
				for( i=0 ; i<2 ; i++ ) for( j=0 ; j<2 ; j++ )
					n.neighbors[fx1+i][fy1+j][fz2  ] = temp.neighbors[cx1][cy1][cz2]->children + Cube::CornerIndex( i , j , z1 );
			if( temp.neighbors[cx1][cy2][cz1] && temp.neighbors[cx1][cy2][cz1]->children )
				for( i=0 ; i<2 ; i++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx1+i][fy2  ][fz1+k] = temp.neighbors[cx1][cy2][cz1]->children + Cube::CornerIndex( i , y1 , k );
			if( temp.neighbors[cx2][cy1][cz1] && temp.neighbors[cx2][cy1][cz1]->children )
				for( j=0 ; j<2 ; j++ ) for( k=0 ; k<2 ; k++ )
					n.neighbors[fx2  ][fy1+j][fz1+k] = temp.neighbors[cx2][cy1][cz1]->children + Cube::CornerIndex( x1 , j , k );
			if( temp.neighbors[cx1][cy2][cz2] && temp.neighbors[cx1][cy2][cz2]->children )
				for( i=0 ; i<2 ; i++ )
					n.neighbors[fx1+i][fy2  ][fz2  ] = temp.neighbors[cx1][cy2][cz2]->children + Cube::CornerIndex( i , y1 , z1 );
			if( temp.neighbors[cx2][cy1][cz2] && temp.neighbors[cx2][cy1][cz2]->children )
				for( j=0 ; j<2 ; j++ )
					n.neighbors[fx2  ][fy1+j][fz2  ] = temp.neighbors[cx2][cy1][cz2]->children + Cube::CornerIndex( x1 , j , z1 );
			if( temp.neighbors[cx2][cy2][cz1] && temp.neighbors[cx2][cy2][cz1]->children )
				for( k=0 ; k<2 ; k++ )
					n.neighbors[fx2  ][fy2  ][fz1+k] = temp.neighbors[cx2][cy2][cz1]->children + Cube::CornerIndex( x1 , y1 , k );
			if( temp.neighbors[cx2][cy2][cz2] && temp.neighbors[cx2][cy2][cz2]->children )
				n.neighbors[fx2  ][fy2  ][fz2  ] = temp.neighbors[cx2][cy2][cz2]->children + Cube::CornerIndex( x1 , y1 , z1 );
		}
	}
	return neighbors[d];
}
template< class NodeData >
typename OctNode< NodeData >::Neighbors5& OctNode< NodeData >::NeighborKey5::setNeighbors( OctNode* node , int xStart , int xEnd , int yStart , int yEnd , int zStart , int zEnd )
{
	int d=node->depth();//先确定当前node在octree的哪一层
	if( node!=neighbors[d].neighbors[2][2][2] )//为什么感觉Neighbor5与当前参数node没有明显的从属关系，否则为什么要先在此进行判断
	{
		neighbors[d].clear();

		if( !node->parent ) neighbors[d].neighbors[2][2][2] = node;//如果是根结点，那么直接把位置(2,2,2)指定为根结点
		else
		{
			setNeighbors( node->parent , xStart , xEnd , yStart , yEnd , zStart , zEnd );
			Neighbors5& temp = neighbors[d-1];
			int x1 , y1 , z1 , x2 , y2 , z2 , ii , jj , kk;
			int idx = int( node-node->parent->children );//得到当前node在其parent的所有child中的index偏移
			Cube::FactorCornerIndex( idx , x1 , y1 , z1 );//利用index计算出当前node在parent划分中属于八块中的哪一块

			for( int i=xStart ; i<xEnd ; i++ )//简单推算了一下，xStart(yStart, zStart)的取值可能为0，xEnd(yEnd, zEnd)的取值可能为5
			{
				x2 = i+x1;
				ii = x2&1;//保留x2最低位
				x2 = 1+(x2>>1);//经推算
				for( int j=yStart ; j<yEnd ; j++ )
				{
					y2 = j+y1;//y2同理，同x2一样，为什么要强行更改最低位
					jj = y2&1;
					y2 = 1+(y2>>1);
					for( int k=zStart ; k<zEnd ; k++ )
					{
						z2 = k+z1;
						kk = z2&1;
						z2 = 1+(z2>>1);
						if(temp.neighbors[x2][y2][z2] )
						{
							if( !temp.neighbors[x2][y2][z2]->children ) temp.neighbors[x2][y2][z2]->initChildren();
							neighbors[d].neighbors[i][j][k] = temp.neighbors[x2][y2][z2]->children + Cube::CornerIndex(ii,jj,kk);
							//这里不但涉及了当前node的子节点，还包括了当前node的neighbor的子节点
							//这里还涉及一个问题，就是所有的node都是有序编号的吗?那么不同depth之间的序号是怎么实现的，深度优先还是宽度优先???
							//Neighbor node是属于同一个depth之间的关系
						}
					}
				}
			}
		}
	}
	return neighbors[d];
}
template< class NodeData >
typename OctNode< NodeData >::ConstNeighbors5& OctNode< NodeData >::ConstNeighborKey5::getNeighbors( const OctNode* node )
{
	int d=node->depth();
	if( node!=neighbors[d].neighbors[2][2][2] )
	{
		neighbors[d].clear();

		if(!node->parent) neighbors[d].neighbors[2][2][2]=node;
		else
		{
			getNeighbors( node->parent );
			ConstNeighbors5& temp = neighbors[d-1];
			int x1,y1,z1,x2,y2,z2,ii,jj,kk;
			int idx=int(node-node->parent->children);
			Cube::FactorCornerIndex(idx,x1,y1,z1);

			for(int i=0;i<5;i++)
			{
				x2=i+x1;
				ii=x2&1;
				x2=1+(x2>>1);
				for(int j=0;j<5;j++)
				{
					y2=j+y1;
					jj=y2&1;
					y2=1+(y2>>1);
					for(int k=0;k<5;k++)
					{
						z2=k+z1;
						kk=z2&1;
						z2=1+(z2>>1);
						if(temp.neighbors[x2][y2][z2] && temp.neighbors[x2][y2][z2]->children)
							neighbors[d].neighbors[i][j][k] = temp.neighbors[x2][y2][z2]->children + Cube::CornerIndex(ii,jj,kk);
					}
				}
			}
		}
	}
	return neighbors[d];
}


template< class NodeData >
int OctNode< NodeData >::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int ret=write(fp);
	fclose(fp);
	return ret;
}
template< class NodeData >
int OctNode< NodeData >::write(FILE* fp) const{
	fwrite(this,sizeof(OctNode< NodeData >),1,fp);//这写进去的是指针???
	if(children){for(int i=0;i<Cube::CORNERS;i++){children[i].write(fp);}}
	return 1;
}
template< class NodeData >
int OctNode< NodeData >::read(const char* fileName){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int ret=read(fp);
	fclose(fp);
	return ret;
}
template< class NodeData >
int OctNode< NodeData >::read(FILE* fp){
	fread(this,sizeof(OctNode< NodeData >),1,fp);//不太明白读出来的是什么信息
	parent=NULL;
	if(children){
		children=NULL;
		initChildren();
		for(int i=0;i<Cube::CORNERS;i++){
			children[i].read(fp);
			children[i].parent=this;
		}
	}
	return 1;
}
template< class NodeData >
int OctNode< NodeData >::width(int maxDepth) const {
	int d=depth();//返回当前节点的深度
	return 1<<(maxDepth-d); //但是这个width所代表的意义是???代表当前节点包含的最大深度时的子节点个数???
}
template< class NodeData >
void OctNode< NodeData >::centerIndex(int maxDepth,int index[DIMENSION]) const
{
	//具体意思还不清楚，后面再看
	int d,o[3];
	depthAndOffset(d,o);//当前的深度和offset
	for(int i=0;i<DIMENSION;i++) index[i]=BinaryNode::CornerIndex( maxDepth , d+1 , o[i]<<1 , 1 );
}
