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


template< class NodeData >
class OctNode
{
private:
	static int UseAlloc;//只有一个副本
	unsigned long long _depthAndOffset;//在OctNode中，把当前节点的depth和offset存在一个变量中，用特定的偏移位来区分
	//虽然不明白为什么要这么做，难道是为了节省变量

	class AdjacencyCountFunction
	{
	public:
		int count;
		void Function( const OctNode< NodeData >* node1 , const OctNode< NodeData >* node2 );//function的作用就是为了count++???
	};
	template<class NodeAdjacencyFunction>
	void __processNodeFaces(OctNode* node,NodeAdjacencyFunction* F,int cIndex1,int cIndex2,int cIndex3,int cIndex4);
	template< class NodeAdjacencyFunction >
	void __processNodeFaces( const OctNode* node , NodeAdjacencyFunction* F , int cIndex1 , int cIndex2 , int cIndex3 , int cIndex4 ) const;
	template<class NodeAdjacencyFunction>
	void __processNodeEdges(OctNode* node,NodeAdjacencyFunction* F,int cIndex1,int cIndex2);
	template<class NodeAdjacencyFunction>
	void __processNodeNodes(OctNode* node,NodeAdjacencyFunction* F);
	template<class NodeAdjacencyFunction>
	static void __ProcessNodeAdjacentNodes(int dx,int dy,int dz,OctNode* node1,int radius1,OctNode* node2,int radius2,int cWidth2,NodeAdjacencyFunction* F);
	template<class TerminatingNodeAdjacencyFunction>
	static void __ProcessTerminatingNodeAdjacentNodes(int dx,int dy,int dz,OctNode* node1,int radius1,OctNode* node2,int radius2,int cWidth2,TerminatingNodeAdjacencyFunction* F);
	template<class PointAdjacencyFunction>
	static void __ProcessPointAdjacentNodes(int dx,int dy,int dz,OctNode* node2,int radius2,int cWidth2,PointAdjacencyFunction* F);
	template<class NodeAdjacencyFunction>
	static void __ProcessFixedDepthNodeAdjacentNodes(int dx,int dy,int dz,OctNode* node1,int radius1,OctNode* node2,int radius2,int cWidth2,int depth,NodeAdjacencyFunction* F);
	template<class NodeAdjacencyFunction>
	static void __ProcessMaxDepthNodeAdjacentNodes(int dx,int dy,int dz,OctNode* node1,int radius1,OctNode* node2,int radius2,int cWidth2,int depth,NodeAdjacencyFunction* F);

	// This is made private because the division by two has been pulled out.
	//暂时不太清楚下面两个函数的作用???
	static inline int Overlap(int c1,int c2,int c3,int dWidth);//static可以通过类名来调用
	inline static int ChildOverlap(int dx,int dy,int dz,int d,int cRadius2);

	const OctNode* __faceNeighbor(int dir,int off) const;
	const OctNode* __edgeNeighbor(int o,const int i[2],const int idx[2]) const;
	OctNode* __faceNeighbor(int dir,int off,int forceChildren);
	OctNode* __edgeNeighbor(int o,const int i[2],const int idx[2],int forceChildren);
public:
	static const int DepthShift,OffsetShift,OffsetShift1,OffsetShift2,OffsetShift3;
	static const int DepthMask,OffsetMask;

	static Allocator< OctNode > NodeAllocator;
	static int UseAllocator( void );
	static void SetAllocator( int blockSize );

	OctNode* parent;
	OctNode* children;
	NodeData nodeData;

	OctNode(void);
	~OctNode(void);
	int initChildren( void );//初始化，为child节点分配空间，设置子节点的depth和offset

	void depthAndOffset( int& depth , int offset[DIMENSION] ) const; //这个变量在节点声明时就存储了depth和offset，现在只是进行逆操作把对应值取出来
	void centerIndex( int index[DIMENSION] ) const;//根据深度和offset计算出当前节点在binaryNode情况下的center index
	int depth( void ) const;//当前节点的深度
	
	static inline void DepthAndOffset( const long long& index , int& depth , int offset[DIMENSION] );//index记录的是当前node在octree中的绝对索引，现在要转化成相对的offset

	template< class Real > static inline void CenterAndWidth( const long long& index , Point3D< Real >& center , Real& width );//利用绝对索引获得center和width
	static inline int Depth( const long long& index );
	static inline void Index( int depth , const int offset[3] , short& d , short off[DIMENSION] );//这个函数给出当前node在octree中的绝对索引
	static inline unsigned long long Index( int depth , const int offset[3] );//返回值应该与_depthAndOffset相等，offset属于相对偏移
	template< class Real > void centerAndWidth( Point3D<Real>& center , Real& width ) const;//返回当前节点的center和width
	template< class Real > bool isInside( Point3D< Real > p ) const;//判断p是否在当前节点内部

	size_t leaves( void ) const;//当前节点的叶节点数目
	size_t maxDepthLeaves( int maxDepth ) const;//给出在深度小于maxDepth时当前节点包含的叶节点数目
	size_t nodes( void ) const;//返回当前节点的所有子节点数目，包括自己
	int maxDepth( void ) const;//返回当前节点的所有子节点的最大深度

	const OctNode* root( void ) const;//返回根节点

	//下面的八个函数比较奇怪，在调用时并没有指明给定参数与当前节点的关系，在允许赋值为NULL的情况下会出现多种不同的处理情况

	//如果没指明current，则返回当前节点的第一个子节点，有可能是自己
	//如果当期节点存在子节点，则返回子节点中的第一个，如果current为根节点或者current节点就是当前节点，返回值为NULL
	//否则找出下一个sibling的叶节点并返回
	const OctNode* nextLeaf(const OctNode* currentLeaf=NULL) const;
	OctNode* nextLeaf(OctNode* currentLeaf=NULL);
	//返回下一个可能访问到的节点，看代码应该是宽度优先遍历
	const OctNode* nextNode(const OctNode* currentNode=NULL) const;
	OctNode* nextNode(OctNode* currentNode=NULL);
	//下面四个函数存在疑问，不明白当current等于当前节点时，为什么直接返回NULL
	//找出给定current节点的下一个宽度优先遍历节点
	//如果current不存在父节点或者current就是当前节点，就返回NULL
	//有可能会返回其父节点的下一个sibling，如果current不是其父节点的最后一个子节点，就返回其下一个sibling
	const OctNode* nextBranch(const OctNode* current) const;
	OctNode* nextBranch(OctNode* current);
	//找出current节点上一个宽度优先遍历节点
	//如果current不存在父节点或者current就是当前节点，就返回NULL
	//有可能会返回其父节点的上一个sibling，如果current不是其父节点的第一个子节点，就返回其上一个sibling
	const OctNode* prevBranch(const OctNode* current) const;
	OctNode* prevBranch(OctNode* current);

	//maxDepth代表了完全八叉树的深度，因此这里要对所有maxDepth之上的child node进行初始化，迭代进行
	void setFullDepth(int maxDepth);

	void printLeaves(void) const;
	void printRange(void) const;//当前节点在归一化octree中的距离范围

	template<class NodeAdjacencyFunction>
	void processNodeFaces(OctNode* node,NodeAdjacencyFunction* F,int fIndex,int processCurrent=1);
	template< class NodeAdjacencyFunction >
	void processNodeFaces( const OctNode* node , NodeAdjacencyFunction* F , int fIndex , int processCurrent=1 ) const;
	template<class NodeAdjacencyFunction>
	void processNodeEdges(OctNode* node,NodeAdjacencyFunction* F,int eIndex,int processCurrent=1);
	template<class NodeAdjacencyFunction>
	void processNodeCorners(OctNode* node,NodeAdjacencyFunction* F,int cIndex,int processCurrent=1);
	template<class NodeAdjacencyFunction>
	void processNodeNodes(OctNode* node,NodeAdjacencyFunction* F,int processCurrent=1);
	
	template<class NodeAdjacencyFunction>
	static void ProcessNodeAdjacentNodes(int maxDepth,OctNode* node1,int width1,OctNode* node2,int width2,NodeAdjacencyFunction* F,int processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessNodeAdjacentNodes(int dx,int dy,int dz,OctNode* node1,int radius1,OctNode* node2,int radius2,int width2,NodeAdjacencyFunction* F,int processCurrent=1);
	template<class TerminatingNodeAdjacencyFunction>
	static void ProcessTerminatingNodeAdjacentNodes(int maxDepth,OctNode* node1,int width1,OctNode* node2,int width2,TerminatingNodeAdjacencyFunction* F,int processCurrent=1);
	template<class TerminatingNodeAdjacencyFunction>
	static void ProcessTerminatingNodeAdjacentNodes(int dx,int dy,int dz,OctNode* node1,int radius1,OctNode* node2,int radius2,int width2,TerminatingNodeAdjacencyFunction* F,int processCurrent=1);
	template<class PointAdjacencyFunction>
	static void ProcessPointAdjacentNodes(int maxDepth,const int center1[3],OctNode* node2,int width2,PointAdjacencyFunction* F,int processCurrent=1);
	template<class PointAdjacencyFunction>
	static void ProcessPointAdjacentNodes(int dx,int dy,int dz,OctNode* node2,int radius2,int width2,PointAdjacencyFunction* F,int processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessFixedDepthNodeAdjacentNodes(int maxDepth,OctNode* node1,int width1,OctNode* node2,int width2,int depth,NodeAdjacencyFunction* F,int processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessFixedDepthNodeAdjacentNodes(int dx,int dy,int dz,OctNode* node1,int radius1,OctNode* node2,int radius2,int width2,int depth,NodeAdjacencyFunction* F,int processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessMaxDepthNodeAdjacentNodes(int maxDepth,OctNode* node1,int width1,OctNode* node2,int width2,int depth,NodeAdjacencyFunction* F,int processCurrent=1);
	template<class NodeAdjacencyFunction>
	static void ProcessMaxDepthNodeAdjacentNodes(int dx,int dy,int dz,OctNode* node1,int radius1,OctNode* node2,int radius2,int width2,int depth,NodeAdjacencyFunction* F,int processCurrent=1);

	//根据点的位置与中心点位置关系来确定corner index，找出当前节点中离p最近的corner
	template< class Real > static int CornerIndex( const Point3D<Real>& center , const Point3D<Real> &p );

	//返回指定(face/edge/corner)index下的Neighbor，face是6个，edge是12个，corner是8个，加起来就是一环邻域内的26个Neighbor node
	OctNode* faceNeighbor(int faceIndex,int forceChildren=0);
	const OctNode* faceNeighbor(int faceIndex) const;
	OctNode* edgeNeighbor(int edgeIndex,int forceChildren=0);
	const OctNode* edgeNeighbor(int edgeIndex) const;
	OctNode* cornerNeighbor(int cornerIndex,int forceChildren=0);
	const OctNode* cornerNeighbor(int cornerIndex) const;
	//找出离p最近的当前节点下的子节点，有可能是自己
	template< class Real > OctNode* getNearestLeaf(const Point3D<Real>& p);
	template< class Real > const OctNode* getNearestLeaf(const Point3D<Real>& p) const;

	//目的是判断给出的两条边是否在共同的父节点上处于同一条公共边上，而不是判断两条边是否是重合的边
	static int CommonEdge(const OctNode* node1,int eIndex1,const OctNode* node2,int eIndex2);
	//判断depth是否相同
	static int CompareForwardDepths(const void* v1,const void* v2);
	//根据节点的深度以及offset的值判断两个node是否相同，不同还需要返回差异值
	static int CompareByDepthAndXYZ( const void* v1 , const void* v2 );
	//暂时不明白这个函数的用意何在???
	static int CompareByDepthAndZIndex( const void* v1 , const void* v2 );
	//先判断是否在同一深度上，再回溯父节点，直到遇到公共父节点，找出二者在父节点相同时的offset差异
	static int CompareForwardPointerDepths(const void* v1,const void* v2);
	//没看明白下面两个函数的用处
	static int CompareBackwardDepths(const void* v1,const void* v2);
	static int CompareBackwardPointerDepths(const void* v1,const void* v2);


	template<class NodeData2>
	OctNode& operator = ( const OctNode< NodeData2 >& node );

	//没看懂multiplier1/multiplier2的意思???
	template< class Real >
	static inline int Overlap2(const int &depth1,const int offSet1[DIMENSION],const Real& multiplier1,const int &depth2,const int offSet2[DIMENSION],const Real& multiplier2);


	int write(const char* fileName) const;
	int write(FILE* fp) const;
	int read(const char* fileName);
	int read(FILE* fp);

	class Neighbors5//比一环邻域多了一层，考虑了二环邻域125个neighbor node
	{
	public:
		OctNode* neighbors[5][5][5];
		Neighbors5( void );
		void clear( void );
	};
	class ConstNeighbors5
	{
	public:
		const OctNode* neighbors[5][5][5];
		ConstNeighbors5( void );
		void clear( void );
	};

	class NeighborKey5
	{
		int _depth;
	public:
		Neighbors5* neighbors;//存储在[0,_depth]深度范围内每一层的neighbors，而且这个neighbors对象为Neighbor5，表示二环邻域节点

		NeighborKey5( void );
		~NeighborKey5( void );

		void set( int depth );
		Neighbors5& getNeighbors( OctNode* node );//下面两个函数都应该算作是设置Neighbor
		Neighbors5& setNeighbors( OctNode* node ,  int xStart=0 , int xEnd=5 , int yStart=0 , int yEnd=5 , int zStart=0 , int zEnd=5 );
	};
	class ConstNeighborKey5
	{
		int _depth;
	public:
		ConstNeighbors5* neighbors;

		ConstNeighborKey5( void );
		~ConstNeighborKey5( void );

		void set( int depth );
		ConstNeighbors5& getNeighbors( const OctNode* node );
	};

	class Neighbors3//27个neighbors
	{
	public:
		OctNode* neighbors[3][3][3];
		Neighbors3( void );
		void clear( void );
	};
	class ConstNeighbors3
	{
	public:
		const OctNode* neighbors[3][3][3];
		ConstNeighbors3( void );
		void clear( void );
	};
	class NeighborKey3
	{
		int _depth;
	public:
		Neighbors3* neighbors;//octree node的27个neighbor，实际应该是26个，而且根据set的内容，这应该是一个数组
		//数组中每一项代表不同的depth下的27个neighbor

		NeighborKey3( void );
		NeighborKey3( const NeighborKey3& key3 );
		~NeighborKey3( void );
		int depth( void ) const { return _depth; }

		void set( int depth );
		//每次找的角点都是离p最近的一个，其对应的child node才会被更新Neighbor，那么其他的没有遍历到是不是就不管了
		template< class Real > Neighbors3& setNeighbors( OctNode* root , Point3D< Real > p , int d );
		//这个get函数除了没有initChild操作，以及多加了一些对child node存在性的判断之外，其他的与set完全相同，不明白其用意???
		template< class Real > Neighbors3& getNeighbors( OctNode* root , Point3D< Real > p , int d );		
		Neighbors3& setNeighbors( OctNode* node , bool flags[3][3][3] );
		Neighbors3& setNeighbors( OctNode* node );
		Neighbors3& getNeighbors( OctNode* node );
		void setNeighbors( OctNode* node , typename OctNode< NodeData >::Neighbors5& neighbors );
		void getNeighbors( OctNode* node , typename OctNode< NodeData >::Neighbors5& neighbors );

		template< class Real > bool setChildNeighbors( Point3D< Real > p , int d , Neighbors3& childNeighbors ) const;
		template< class Real > bool getChildNeighbors( Point3D< Real > p , int d , Neighbors3& childNeighbors ) const;
	};
	class ConstNeighborKey3
	{
		int _depth;
	public:
		ConstNeighbors3* neighbors;

		ConstNeighborKey3( void );
		ConstNeighborKey3( const ConstNeighborKey3& key3 );
		~ConstNeighborKey3( void );
		int depth( void ) const { return _depth; }

		void set( int depth );
		ConstNeighbors3& getNeighbors( const OctNode* node );
		ConstNeighbors3& getNeighbors( const OctNode* node , int minDepth );
		void getNeighbors( const OctNode* node , typename OctNode< NodeData >::ConstNeighbors5& neighbors );
	};

	void centerIndex(int maxDepth,int index[DIMENSION]) const;
	int width(int maxDepth) const;
};


#include "Octree.inl"

#endif // OCT_NODE
