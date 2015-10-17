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
	static int UseAlloc;//ֻ��һ������
	unsigned long long _depthAndOffset;//��OctNode�У��ѵ�ǰ�ڵ��depth��offset����һ�������У����ض���ƫ��λ������
	//��Ȼ������ΪʲôҪ��ô�����ѵ���Ϊ�˽�ʡ����

	class AdjacencyCountFunction
	{
	public:
		int count;
		void Function( const OctNode< NodeData >* node1 , const OctNode< NodeData >* node2 );//function�����þ���Ϊ��count++???
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
	//��ʱ��̫���������������������???
	static inline int Overlap(int c1,int c2,int c3,int dWidth);//static����ͨ������������
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
	int initChildren( void );//��ʼ����Ϊchild�ڵ����ռ䣬�����ӽڵ��depth��offset

	void depthAndOffset( int& depth , int offset[DIMENSION] ) const; //��������ڽڵ�����ʱ�ʹ洢��depth��offset������ֻ�ǽ���������Ѷ�Ӧֵȡ����
	void centerIndex( int index[DIMENSION] ) const;//������Ⱥ�offset�������ǰ�ڵ���binaryNode����µ�center index
	int depth( void ) const;//��ǰ�ڵ�����
	
	static inline void DepthAndOffset( const long long& index , int& depth , int offset[DIMENSION] );//index��¼���ǵ�ǰnode��octree�еľ�������������Ҫת������Ե�offset

	template< class Real > static inline void CenterAndWidth( const long long& index , Point3D< Real >& center , Real& width );//���þ����������center��width
	static inline int Depth( const long long& index );
	static inline void Index( int depth , const int offset[3] , short& d , short off[DIMENSION] );//�������������ǰnode��octree�еľ�������
	static inline unsigned long long Index( int depth , const int offset[3] );//����ֵӦ����_depthAndOffset��ȣ�offset�������ƫ��
	template< class Real > void centerAndWidth( Point3D<Real>& center , Real& width ) const;//���ص�ǰ�ڵ��center��width
	template< class Real > bool isInside( Point3D< Real > p ) const;//�ж�p�Ƿ��ڵ�ǰ�ڵ��ڲ�

	size_t leaves( void ) const;//��ǰ�ڵ��Ҷ�ڵ���Ŀ
	size_t maxDepthLeaves( int maxDepth ) const;//���������С��maxDepthʱ��ǰ�ڵ������Ҷ�ڵ���Ŀ
	size_t nodes( void ) const;//���ص�ǰ�ڵ�������ӽڵ���Ŀ�������Լ�
	int maxDepth( void ) const;//���ص�ǰ�ڵ�������ӽڵ��������

	const OctNode* root( void ) const;//���ظ��ڵ�

	//����İ˸������Ƚ���֣��ڵ���ʱ��û��ָ�����������뵱ǰ�ڵ�Ĺ�ϵ��������ֵΪNULL������»���ֶ��ֲ�ͬ�Ĵ������

	//���ûָ��current���򷵻ص�ǰ�ڵ�ĵ�һ���ӽڵ㣬�п������Լ�
	//������ڽڵ�����ӽڵ㣬�򷵻��ӽڵ��еĵ�һ�������currentΪ���ڵ����current�ڵ���ǵ�ǰ�ڵ㣬����ֵΪNULL
	//�����ҳ���һ��sibling��Ҷ�ڵ㲢����
	const OctNode* nextLeaf(const OctNode* currentLeaf=NULL) const;
	OctNode* nextLeaf(OctNode* currentLeaf=NULL);
	//������һ�����ܷ��ʵ��Ľڵ㣬������Ӧ���ǿ�����ȱ���
	const OctNode* nextNode(const OctNode* currentNode=NULL) const;
	OctNode* nextNode(OctNode* currentNode=NULL);
	//�����ĸ������������ʣ������׵�current���ڵ�ǰ�ڵ�ʱ��Ϊʲôֱ�ӷ���NULL
	//�ҳ�����current�ڵ����һ��������ȱ����ڵ�
	//���current�����ڸ��ڵ����current���ǵ�ǰ�ڵ㣬�ͷ���NULL
	//�п��ܻ᷵���丸�ڵ����һ��sibling�����current�����丸�ڵ�����һ���ӽڵ㣬�ͷ�������һ��sibling
	const OctNode* nextBranch(const OctNode* current) const;
	OctNode* nextBranch(OctNode* current);
	//�ҳ�current�ڵ���һ��������ȱ����ڵ�
	//���current�����ڸ��ڵ����current���ǵ�ǰ�ڵ㣬�ͷ���NULL
	//�п��ܻ᷵���丸�ڵ����һ��sibling�����current�����丸�ڵ�ĵ�һ���ӽڵ㣬�ͷ�������һ��sibling
	const OctNode* prevBranch(const OctNode* current) const;
	OctNode* prevBranch(OctNode* current);

	//maxDepth��������ȫ�˲�������ȣ��������Ҫ������maxDepth֮�ϵ�child node���г�ʼ������������
	void setFullDepth(int maxDepth);

	void printLeaves(void) const;
	void printRange(void) const;//��ǰ�ڵ��ڹ�һ��octree�еľ��뷶Χ

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

	//���ݵ��λ�������ĵ�λ�ù�ϵ��ȷ��corner index���ҳ���ǰ�ڵ�����p�����corner
	template< class Real > static int CornerIndex( const Point3D<Real>& center , const Point3D<Real> &p );

	//����ָ��(face/edge/corner)index�µ�Neighbor��face��6����edge��12����corner��8��������������һ�������ڵ�26��Neighbor node
	OctNode* faceNeighbor(int faceIndex,int forceChildren=0);
	const OctNode* faceNeighbor(int faceIndex) const;
	OctNode* edgeNeighbor(int edgeIndex,int forceChildren=0);
	const OctNode* edgeNeighbor(int edgeIndex) const;
	OctNode* cornerNeighbor(int cornerIndex,int forceChildren=0);
	const OctNode* cornerNeighbor(int cornerIndex) const;
	//�ҳ���p����ĵ�ǰ�ڵ��µ��ӽڵ㣬�п������Լ�
	template< class Real > OctNode* getNearestLeaf(const Point3D<Real>& p);
	template< class Real > const OctNode* getNearestLeaf(const Point3D<Real>& p) const;

	//Ŀ�����жϸ������������Ƿ��ڹ�ͬ�ĸ��ڵ��ϴ���ͬһ���������ϣ��������ж��������Ƿ����غϵı�
	static int CommonEdge(const OctNode* node1,int eIndex1,const OctNode* node2,int eIndex2);
	//�ж�depth�Ƿ���ͬ
	static int CompareForwardDepths(const void* v1,const void* v2);
	//���ݽڵ������Լ�offset��ֵ�ж�����node�Ƿ���ͬ����ͬ����Ҫ���ز���ֵ
	static int CompareByDepthAndXYZ( const void* v1 , const void* v2 );
	//��ʱ����������������������???
	static int CompareByDepthAndZIndex( const void* v1 , const void* v2 );
	//���ж��Ƿ���ͬһ����ϣ��ٻ��ݸ��ڵ㣬ֱ�������������ڵ㣬�ҳ������ڸ��ڵ���ͬʱ��offset����
	static int CompareForwardPointerDepths(const void* v1,const void* v2);
	//û���������������������ô�
	static int CompareBackwardDepths(const void* v1,const void* v2);
	static int CompareBackwardPointerDepths(const void* v1,const void* v2);


	template<class NodeData2>
	OctNode& operator = ( const OctNode< NodeData2 >& node );

	//û����multiplier1/multiplier2����˼???
	template< class Real >
	static inline int Overlap2(const int &depth1,const int offSet1[DIMENSION],const Real& multiplier1,const int &depth2,const int offSet2[DIMENSION],const Real& multiplier2);


	int write(const char* fileName) const;
	int write(FILE* fp) const;
	int read(const char* fileName);
	int read(FILE* fp);

	class Neighbors5//��һ���������һ�㣬�����˶�������125��neighbor node
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
		Neighbors5* neighbors;//�洢��[0,_depth]��ȷ�Χ��ÿһ���neighbors���������neighbors����ΪNeighbor5����ʾ��������ڵ�

		NeighborKey5( void );
		~NeighborKey5( void );

		void set( int depth );
		Neighbors5& getNeighbors( OctNode* node );//��������������Ӧ������������Neighbor
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

	class Neighbors3//27��neighbors
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
		Neighbors3* neighbors;//octree node��27��neighbor��ʵ��Ӧ����26�������Ҹ���set�����ݣ���Ӧ����һ������
		//������ÿһ�����ͬ��depth�µ�27��neighbor

		NeighborKey3( void );
		NeighborKey3( const NeighborKey3& key3 );
		~NeighborKey3( void );
		int depth( void ) const { return _depth; }

		void set( int depth );
		//ÿ���ҵĽǵ㶼����p�����һ�������Ӧ��child node�Żᱻ����Neighbor����ô������û�б������ǲ��ǾͲ�����
		template< class Real > Neighbors3& setNeighbors( OctNode* root , Point3D< Real > p , int d );
		//���get��������û��initChild�������Լ������һЩ��child node�����Ե��ж�֮�⣬��������set��ȫ��ͬ��������������???
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
