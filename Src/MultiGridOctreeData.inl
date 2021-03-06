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
#include "MyTime.h"
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
int TreeNodeData::NodeCount = 0;
TreeNodeData::TreeNodeData( void ){ nodeIndex = NodeCount++; }
TreeNodeData::~TreeNodeData( void ) { }


////////////
// Octree //
////////////
template< class Real > double Octree< Real >::maxMemoryUsage=0;

template< class Real >
double Octree< Real >::MemoryUsage(void)
{
	double mem = double( MemoryInfo::Usage() ) / (1<<20);
	if( mem>maxMemoryUsage ) maxMemoryUsage=mem;
	return mem;
}

template< class Real >
Octree< Real >::Octree( void )
{
	threads = 1;
	_constrainValues = false;
}

template< class Real >
bool Octree< Real >::_IsInset( const TreeOctNode* node )
{
	int d , off[3];
	node->depthAndOffset( d , off );
	int res = 1<<d , o = 1<<(d-2);
	return ( off[0]>=o && off[0]<res-o && off[1]>=o && off[1]<res-o && off[2]>=o && off[2]<res-o );
}
template< class Real >
bool Octree< Real >::_IsInsetSupported( const TreeOctNode* node )
{
	int d , off[3];
	node->depthAndOffset( d , off );
	int res = 1<<d , o = (1<<(d-2))-1;
	return ( off[0]>=o && off[0]<res-o && off[1]>=o && off[1]<res-o && off[2]>=o && off[2]<res-o );
}
template< class Real >
template< class V >
int Octree< Real >::SplatPointData( TreeOctNode* node , const Point3D< Real >& position , const V& v , SparseNodeData< V >& dataInfo , typename TreeOctNode::ConstNeighborKey3& neighborKey )
{
	double x , dxdy , dxdydz , dx[DIMENSION][SPLAT_ORDER+1];
	double width;
	int off[3];
	typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors( node );
	Point3D<Real> center;
	Real w;
	node->centerAndWidth( center , w );
	width=w;
	for( int i=0 ; i<3 ; i++ )
	{
#if SPLAT_ORDER==2//SPLAT算法计算权重是根据二次幂计算，参考函数UpdateWeightContribution
		off[i] = 0;
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125+1.500*x+0.500*x*x;
		x = ( center[i] - position[i] ) / width;
		dx[i][1] = 0.750        -      x*x;

		dx[i][2] = 1. - dx[i][1] - dx[i][0];
#elif SPLAT_ORDER==1//一次幂
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
				const TreeOctNode* _node = neighbors.neighbors[i][j][k];
				if( _node )
				{
					if( dataInfo.indices.size()<TreeNodeData::NodeCount ) dataInfo.indices.resize( TreeNodeData::NodeCount , -1 );
					int idx = dataInfo.index( _node );
					if( idx<0 )
					{
						dataInfo.indices[ _node->nodeData.nodeIndex ] = (int)dataInfo.data.size();
						dataInfo.data.push_back( v * Real(dxdydz) );
					}
					else dataInfo.data[idx] += v * Real( dxdydz );
				}
			}
	}
	return 0;
}
template< class Real >
template< class V >
int Octree< Real >::SplatPointData( TreeOctNode* node , const Point3D< Real >& position , const V& v , SparseNodeData< V >& dataInfo , typename TreeOctNode::NeighborKey3& neighborKey )
{
	//参考函数UpdateWeightContribution
	double x , dxdy , dxdydz , dx[DIMENSION][SPLAT_ORDER+1];
	double width;
	int off[3];
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node );//找出node的一环邻域octree node
	Point3D<Real> center;
	Real w;
	node->centerAndWidth( center , w );
	width=w;
	for( int i=0 ; i<3 ; i++ )
	{
#if SPLAT_ORDER==2//扩展到周围二环邻域，二环邻域是二次权重
		off[i] = 0;
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125+1.500*x+0.500*x*x;
		x = ( center[i] - position[i] ) / width;
		dx[i][1] = 0.750        -      x*x;

		dx[i][2] = 1. - dx[i][1] - dx[i][0];
#elif SPLAT_ORDER==1//扩展到周围一环邻域，一环邻域是线性权重
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
				if( (int)dataInfo.indices.size()<TreeNodeData::NodeCount ) dataInfo.indices.resize( TreeNodeData::NodeCount , -1 );
				int idx = dataInfo.index( _node );
				if( idx<0 )
				{
					dataInfo.indices[ _node->nodeData.nodeIndex ] = (int)dataInfo.data.size();
					dataInfo.data.push_back( v * Real(dxdydz) );
				}
				else dataInfo.data[idx] += v * Real( dxdydz );
			}
	}
	return 0;
}
template< class Real >
template< class V >
Real Octree< Real >::SplatPointData( ConstPointer( Real ) kernelDensityWeights , const Point3D< Real >& position , const V& v , SparseNodeData< V >& dataInfo , typename TreeOctNode::NeighborKey3& neighborKey , int minDepth , int maxDepth , int dim )
{
	//v存储了被翻转的normal信息，dataInfo是用来存储splatted normal的
	double dx;
	V _v;
	TreeOctNode* temp;
	int cnt=0;
	double width;
	Point3D< Real > myCenter( (Real)0.5 , (Real)0.5 , (Real)0.5 );
	Real myWidth = (Real)1.;

	temp = &tree;
	while( temp->depth()<_splatDepth )//在node Depth小于splatDepth之前进行此操作
	{
		if( !temp->children ) fprintf( stderr , "[ERROR] Octree::SplatPointData\n" ) , exit( 0 );
		int cIndex=TreeOctNode::CornerIndex( myCenter , position );//根据position和center的相对位置判断最近的cube corner index
		temp=&temp->children[cIndex];//移动到这个corner相关的child node
		myWidth/=2;//width减半
		if( cIndex&1 ) myCenter[0] += myWidth/2;//更新center
		else		   myCenter[0] -= myWidth/2;
		if( cIndex&2 ) myCenter[1] += myWidth/2;
		else 	  	   myCenter[1] -= myWidth/2;
		if( cIndex&4 ) myCenter[2] += myWidth/2;
		else 		   myCenter[2] -= myWidth/2;
	}//width和center更新直到node Depth == splatDepth为止
	Real weight , depth;
	//temp当前是在splatDepth下包含position的节点，这里计算的是W_D(s.p)和Depth(s.p)，根据point density调整到合适的深度
	GetSampleDepthAndWeight( kernelDensityWeights , temp , position , neighborKey , depth , weight );

	if( depth<minDepth ) depth = Real(minDepth);
	if( depth>maxDepth ) depth = Real(maxDepth);
	int topDepth=int(ceil(depth));

	dx = 1.0-(topDepth-depth);
	if( topDepth<=minDepth )
	{
		topDepth = minDepth;
		dx=1;
	}
	else if( topDepth>maxDepth )
	{
		topDepth = maxDepth;
		dx=1;
	}
	while( temp->depth()>topDepth ) temp=temp->parent;
	while( temp->depth()<topDepth )
	{
		if(!temp->children) temp->initChildren();
		int cIndex=TreeOctNode::CornerIndex( myCenter , position );
		temp=&temp->children[cIndex];
		myWidth/=2;
		if( cIndex&1 ) myCenter[0] += myWidth/2;
		else		   myCenter[0] -= myWidth/2;
		if( cIndex&2 ) myCenter[1] += myWidth/2;
		else		   myCenter[1] -= myWidth/2;
		if( cIndex&4 ) myCenter[2] += myWidth/2;
		else		   myCenter[2] -= myWidth/2;
	}
	width = 1.0 / ( 1<<temp->depth() );
	_v = v * weight / Real( pow( width , dim ) ) * Real( dx );//v是当前position的normal，
	//上面的式子等于version 1中的for(i=0;i<DIMENSION;i++){n.coords[i]=normal.coords[i]*alpha/Real(pow(width,3))*Real(dx);}
	SplatPointData( temp , position , _v , dataInfo , neighborKey );//dataInfo原来是用来存储splatted normal的
	if( fabs(1.0-dx) > EPSILON )//这里好像是考虑到GetSampleDepthAndWeight中求出的Depth卡在两层之间时，两个层都要进行normal splatting，权重又dx决定
	{
		dx = Real(1.0-dx);
		temp = temp->parent;
		width = 1.0 / ( 1<<temp->depth() );

		_v = v * weight / Real( pow( width , dim ) ) * Real( dx );
		SplatPointData( temp , position , _v , dataInfo , neighborKey );
	}
	return weight;
}
template< class Real >
template< class V >
void Octree< Real >::MultiSplatPointData( ConstPointer( Real ) kernelDensityWeights , const Point3D< Real >& position , const V& v , SparseNodeData< V >& dataInfo , typename TreeOctNode::NeighborKey3& neighborKey , int maxDepth , int dim )
{
	Real _depth , weight;
	if( kernelDensityWeights ) GetSampleDepthAndWeight( kernelDensityWeights , position , neighborKey , _depth , weight );
	else weight = (Real)1. , _depth = (Real)maxDepth;
	int depth = std::min< int >( maxDepth , (int)ceil( _depth ) );
	V _v = v * weight;

	Point3D< Real > myCenter( (Real)0.5 , (Real)0.5 , (Real)0.5 );
	Real myWidth = (Real)1.;

	TreeOctNode* temp = &tree;
	while( temp->depth()<=depth )
	{
		SplatPointData( temp , position , _v * Real( pow( (float)(1<<temp->depth()) , dim ) ) , dataInfo , neighborKey );
		if( temp->depth()<depth )
		{
			if( !temp->children ) temp->initChildren();
			int cIndex = TreeOctNode::CornerIndex( myCenter , position );
			temp = &temp->children[cIndex];
			myWidth /= 2;
			if( cIndex&1 ) myCenter[0] += myWidth/2;
			else		   myCenter[0] -= myWidth/2;
			if( cIndex&2 ) myCenter[1] += myWidth/2;
			else 	  	   myCenter[1] -= myWidth/2;
			if( cIndex&4 ) myCenter[2] += myWidth/2;
			else 		   myCenter[2] -= myWidth/2;
		}
		else break;
	}
}
template< class Real >
template< class V >
void Octree< Real >::MultiSplatPointData( ConstPointer( Real ) kernelDensityWeights , const Point3D< Real >& position , const V& v , SparseNodeData< V >& dataInfo , typename TreeOctNode::ConstNeighborKey3& neighborKey , int dim )
{
	Real depth , weight;
	if( kernelDensityWeights ) GetSampleDepthAndWeight( kernelDensityWeights , position , neighborKey , depth , weight );
	else weight = (Real)1.;
	V _v = v * weight;

	Point3D< Real > myCenter( (Real)0.5 , (Real)0.5 , (Real)0.5 );
	Real myWidth = (Real)1.;

	TreeOctNode* temp = &tree;
	while( true )
	{
		SplatPointData( temp , position , _v * Real( pow( 1<<temp->depth() , dim ) ) , dataInfo , neighborKey );
		if( !temp->children ) break;
		int cIndex = TreeOctNode::CornerIndex( myCenter , position );
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if( cIndex&1 ) myCenter[0] += myWidth/2;
		else		   myCenter[0] -= myWidth/2;
		if( cIndex&2 ) myCenter[1] += myWidth/2;
		else 	  	   myCenter[1] -= myWidth/2;
		if( cIndex&4 ) myCenter[2] += myWidth/2;
		else 		   myCenter[2] -= myWidth/2;
	}
}

template< class Real >
void Octree< Real >::GetSampleDepthAndWeight( ConstPointer( Real ) kernelDensityWeights , const Point3D< Real >& position , typename TreeOctNode::ConstNeighborKey3& neighborKey , Real& depth , Real& weight )
{
	TreeOctNode* temp;
	Point3D< Real > myCenter( (Real)0.5 , (Real)0.5 , (Real)0.5 );
	Real myWidth = Real( 1. );

	// Get the finest node with depth less than or equal to the splat depth that contains the point
	temp = &tree;
	while( temp->depth()<_splatDepth )
	{
		if( !temp->children ) break;// fprintf( stderr , "[ERROR] Octree::GetSampleDepthAndWeight\n" ) , exit( 0 );
		int cIndex = TreeOctNode::CornerIndex( myCenter , position );
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if( cIndex&1 ) myCenter[0] += myWidth/2;
		else		   myCenter[0] -= myWidth/2;
		if( cIndex&2 ) myCenter[1] += myWidth/2;
		else		   myCenter[1] -= myWidth/2;
		if( cIndex&4 ) myCenter[2] += myWidth/2;
		else		   myCenter[2] -= myWidth/2;
	}
	return GetSampleDepthAndWeight( kernelDensityWeights , temp , position , neighborKey , depth , weight );
}
template< class Real >
void Octree< Real >::GetSampleDepthAndWeight( ConstPointer( Real ) kernelDensityWeights , const Point3D< Real >& position , typename TreeOctNode::NeighborKey3& neighborKey , Real& depth , Real& weight )
{
	TreeOctNode* temp;
	Point3D< Real > myCenter( (Real)0.5 , (Real)0.5 , (Real)0.5 );
	Real myWidth = Real( 1. );

	// Get the finest node with depth less than or equal to the splat depth that contains the point
	temp = &tree;
	while( temp->depth()<_splatDepth )
	{
		if( !temp->children ) break;// fprintf( stderr , "[ERROR] Octree::GetSampleDepthAndWeight\n" ) , exit( 0 );
		int cIndex = TreeOctNode::CornerIndex( myCenter , position );
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if( cIndex&1 ) myCenter[0] += myWidth/2;
		else		   myCenter[0] -= myWidth/2;
		if( cIndex&2 ) myCenter[1] += myWidth/2;
		else		   myCenter[1] -= myWidth/2;
		if( cIndex&4 ) myCenter[2] += myWidth/2;
		else		   myCenter[2] -= myWidth/2;
	}
	return GetSampleDepthAndWeight( kernelDensityWeights , temp , position , neighborKey , depth , weight );
}
template< class Real >
void Octree< Real >::GetSampleDepthAndWeight( ConstPointer( Real ) kernelDensityWeights , const TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::ConstNeighborKey3& neighborKey , Real& depth , Real& weight )
{
	const TreeOctNode* temp=node;
	weight = GetSamplesPerNode( kernelDensityWeights , temp , position , neighborKey );
	if( weight>=(Real)1. ) depth = Real( temp->depth() + log( weight ) / log(double(1<<(DIMENSION-1))) );
	else
	{
		Real oldWeight , newWeight;
		oldWeight = newWeight = weight;
		while( newWeight<(Real)1. && temp->parent )
		{
			temp=temp->parent;
			oldWeight = newWeight;
			newWeight = GetSamplesPerNode( kernelDensityWeights , temp , position, neighborKey );
		}
		depth = Real( temp->depth() + log( newWeight ) / log( newWeight / oldWeight ) );
	}
	weight = Real( pow( double(1<<(DIMENSION-1)) , -double(depth) ) );
}
template< class Real >
void Octree< Real >::GetSampleDepthAndWeight( ConstPointer( Real ) kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , Real& depth , Real& weight )
{
	TreeOctNode* temp=node;
	//返回了当前point position在node的加权weight，相当于评估了每个neighbor node的point分布密度后，根据不同node的远近关系，
	//对当前position的影响又依靠再一次的splat weight决定，所以就好像把splat weight乘了两次，
	//实际上表示的意思是position周围point越多，靠近的距离越近，该返回值越大，position对于当前的贡献应该会越小
	weight = GetSamplesPerNode( kernelDensityWeights , temp , position , neighborKey );
	//SGP07文中的Depth_s.p，根据position的density调整其计算的Depth
	//公式是log(4)(W'/W)，但由于在计算node weight时已经用smaplesPerNode等进行了缩放，这里只用了log(weight)，而不是version1中的log(weight/(samplesPerNode+1)
	//所以density较大的情况下要向更大深度拓展，在较小的情况下要向更小深度收缩，如下所示
	if( weight>=(Real)1. ) depth = Real( temp->depth() + log( weight ) / log(double(1<<(DIMENSION-1))) );
	else
	{
		Real oldWeight , newWeight;
		oldWeight = newWeight = weight;
		while( newWeight<(Real)1. && temp->parent )
		{
			temp=temp->parent;
			oldWeight = newWeight;//不断向上刷新point在整个branch上的加权weight
			newWeight = GetSamplesPerNode( kernelDensityWeights , temp , position, neighborKey );//GetSamplesPerNode越往上走得到的weight会越大
		}
		depth = Real( temp->depth() + log( newWeight ) / log( newWeight / oldWeight ) );
	}
	weight = Real( pow( double(1<<(DIMENSION-1)) , -double(depth) ) );
}
template< class Real >
Real Octree< Real >::GetSamplesPerNode( ConstPointer( Real ) kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey )
{
	Real weight=0;
	double x , dxdy , dx[DIMENSION][3];
	double width;
	//node是当前position在splatDepth层存在的node，根据node找出position在splatDepth层的27个Neighbor
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node );
	Point3D<Real> center;
	Real w;
	node->centerAndWidth( center , w );
	width = w;

	for( int i=0 ; i<DIMENSION ; i++ )//0.125/0.75/0.125权重
	{
		x = ( center[i] - position[i] - width ) / width;//参考函数UpdateWeightContribution
		dx[i][0] = 1.125 + 1.500*x + 0.500*x*x;
		x = ( center[i] - position[i] ) / width;
		dx[i][1] = 0.750           -       x*x;

		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
	}
	//把node得到的权重累加通过再一次权重累加的形式返回给position代表的point，推测position周围的density
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		dxdy = dx[0][i] * dx[1][j];
		for( int k=0 ; k<3 ; k++ ) if( neighbors.neighbors[i][j][k] )
			weight += Real( dxdy * dx[2][k] * kernelDensityWeights[ neighbors.neighbors[i][j][k]->nodeData.nodeIndex ] );
	}
	return weight;//返回了当前node的加权weight
}
template< class Real >
Real Octree< Real >::GetSamplesPerNode( ConstPointer( Real ) kernelDensityWeights , const TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::ConstNeighborKey3& neighborKey )
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
		x = ( center[i] - position[i] - width ) / width;//参考函数UpdateWeightContribution
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
	return weight;
}
template< class Real >
int Octree< Real >::UpdateWeightContribution( std::vector< Real >& kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , Real weight )
{
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node );//设置了单层的一环邻域neighbor
	if( (int)kernelDensityWeights.size()<TreeNodeData::NodeCount ) kernelDensityWeights.resize( TreeNodeData::NodeCount , 0 );//原来的weight保持不变，新分配的空间被初始化为0
	double x , dxdy , dx[DIMENSION][3] , width;
	Point3D< Real > center;
	Real w;
	node->centerAndWidth( center , w );//当前node的center position和width
	width = w;

	for( int i=0 ; i<DIMENSION ; i++ )
	{
		x = ( center[i] - position[i] - width ) / width;//经验证，就是注释中所谓的0.125/0.75/0.125的划分，至于权重为什么这么设定，论文中似乎没有提到
		//而且，position点是在当前center node中的，下面的第一个x的值是p的某一维坐标到最左侧Neighbor中心点的距离，第二个x是
		//到自己所在node的中心点的距离，根据下面公式计算，当p的位置在自己所在node的center时，可以得到0.125/0.75/0.125的默认权重，
		//偏向左侧Neighbor时，权重逐渐变成0.5/0.5/0，偏向右侧Neighbor时，权重逐渐变成0/0.5/0.5，两侧权重都呈现出单调增减趋势
		//中间node的权重会根据position到center的距离出现抛物线变化
		dx[i][0] = 1.125 + 1.500*x + 0.500*x*x;
		dx[i][1] = -0.25 - 2.*x - x*x;
		dx[i][2] = 1. - dx[i][1] - dx[i][0];
	}

	// Suppose that the samples are uniformly placed along the middle of the three slices.
	// Then splatting the points we get coefficients:
	//		0.125 / 0.75 / 0.125 across the three slices.
	// Sampling at the center slice we get:
	//		0.125^2 + 0.75^2 + 0.125^2 = 19/32
	//这块对weight的缩放实在是不太理解，所有weight加起来和为1，为什么还要放大这些计算好的weight
	//splat算法就是投影到二维上，所以这里考虑二维九宫格的计算，中间格为这几个weight平方和
	const double SAMPLE_SCALE = 1. / ( 0.125 * 0.125 + 0.75 * 0.75 + 0.125 * 0.125 );
	weight *= (Real)SAMPLE_SCALE;

	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		//早期版本的poisson中并没有weight的进一步缩放，这里的weight考虑了normal大小，考虑了samplesPerNode的值，还有上面的缩放
		dxdy = dx[0][i] * dx[1][j] * weight;
		TreeOctNode** _neighbors = neighbors.neighbors[i][j];
		//要注意splat算法中其含义是将三维投射到二维中，所以在最后一维的weight好像是多个point sampleweight的叠加，具体可Google splatting(雪球算法)
		for( int k=0 ; k<3 ; k++ ) if( _neighbors[k] ) kernelDensityWeights[ _neighbors[k]->nodeData.nodeIndex ] += Real( dxdy * dx[2][k] );//第三个轴的叠加作用在这里区分
	}
	return 0;
}
template< class Real >
bool Octree< Real >::_InBounds( Point3D< Real > p ) const
{
	if( _boundaryType==0 ){ if( p[0]<Real(0.25) || p[0]>Real(0.75) || p[1]<Real(0.25) || p[1]>Real(0.75) || p[2]<Real(0.25) || p[2]>Real(0.75) ) return false; }
	else                  { if( p[0]<Real(0.00) || p[0]>Real(1.00) || p[1]<Real(0.00) || p[1]>Real(1.00) || p[2]<Real(0.00) || p[2]>Real(1.00) ) return false; }
	return true;
}
template< class Real >
template< class PointReal >
int Octree< Real >::SetTree( OrientedPointStream< PointReal >* pointStream , int minDepth , int maxDepth , int fullDepth , 
							int splatDepth , Real samplesPerNode , Real scaleFactor ,
							bool useConfidence , bool useNormalWeights , Real constraintWeight , int adaptiveExponent ,
							std::vector< Real >& kernelDensityWeights , SparseNodeData< PointData >& pointInfo , SparseNodeData< Point3D< Real > >& normalInfo , std::vector< Real >& centerWeights ,
							XForm4x4< Real >& xForm , int boundaryType , bool makeComplete )
{
	//kernel Depth在函数中被替换为splatDepth，看来kernel 就是用来计算point density的
	if( splatDepth<0 ) splatDepth = 0;//按照文中所写，splatDepth是比Depth要小，这是为了估算point density设置的
	//看Kazhdan的注释，在早期版本中minDepth就是fullDepth，不知道这里为什么还要区分开
	_boundaryType = boundaryType;
	if     ( _boundaryType<0 ) _boundaryType = -1;//看样子boundaryType有三种不同状态
	else if( _boundaryType>0 ) _boundaryType =  1;
	_splatDepth = splatDepth;//后面splat的时候会用到
	_constrainValues = (constraintWeight>0);//表征是否需要用原点云中数据点作为固定约束点

	XForm3x3< Real > xFormN;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) xFormN(i,j) = xForm(i,j);//应该只是旋转部分
	xFormN = xFormN.transpose().inverse();//这两个操作在xFormN是规范旋转矩阵的情况下互相抵消，这里在干啥???
	minDepth = std::min< int >( minDepth , maxDepth );	// minDepth <= maxDepth
	fullDepth = std::max< int >( minDepth , std::min< int >( fullDepth , maxDepth ) );	// minDepth <= fullDepth <= maxDepth
	//实在不明白为什么会有这么多的depth，fullDepth和maxDepth的区别以及fullDepth的用处???
	//按照Kazhdan的注释，fullDepth是指过了这个深度，octree就可以adaptive，maxDepth是指octree的最大深度，
	//早期版本使用minDepth来完成fullDepth的功能，现在可能用法已经完全不一样了
	//Depth---Note that since the reconstructor adapts the octree to the sampling density, the specified reconstruction depth is only an upper bound.
	// If _boundaryType==0, points are scaled to be in the [0.25,0.75]^3 cube so all depths have to be offset by
	// and the minDepth has to be 2.<---------------WHY???
	if( _boundaryType==0 )//为什么boundaryType等于0这么特殊
	{
		minDepth++ , maxDepth++ , fullDepth++;//深度增加是什么意思???
		if( splatDepth ) splatDepth++;
		minDepth = std::max< int >( minDepth , 2 );
	}
	// Otherwise the points are in the [0,1]^3 cube.
	// However, for Neumann constraints, the function at depth 0 is constant so the system matrix is zero if there
	// is no screening.
#if 0 //if 0是屏蔽掉一些当前用不到但后面可能会用的代码，用的时候换成if 1就可以了
	else if( _boundaryType==1 && !_constrainValues ) minDepth = std::max< int >( minDepth , 1 );
#endif

	//二次样条的某一曲线段只与相应的两段折线段，三个控制多边形顶点有关，改变其中一个顶点，将影响三段样条曲线段。
	//同样的，对三次样条，某一曲线段由相应的三段折线段，四个控制点决定。
	//spline is a numeric function that is piecewise-defined by polynomial functions, 
	//and which possesses a high degree of smoothness at the places where the polynomial pieces connect
	//In the computer-aided design and computer graphics, spline functions are constructed as linear combinations of B-splines with a set of control points.
	_fData.set( maxDepth , _boundaryType );//B样条曲线设置参数，边界条件类型一共有三种，查查NEUMANN boundary type，这种boundary约束只要求
	//边界的梯度保持为0就可以了
	_minDepth = minDepth;
	_fullDepth = fullDepth;
	double pointWeightSum = 0;
	Point3D< Real > min , max , myCenter;
	Real myWidth;
	int i , cnt=0;
	TreeOctNode* temp;

	typename TreeOctNode::NeighborKey3 neighborKey;//同上
	neighborKey.set( maxDepth );//这里只是初始化空间，没有实际的neighbor node赋值操作

	tree.setFullDepth( _fullDepth );//初始化设置child node，fullDepth下的所有child都有

	// Read through once to get the center and scale
	{
		double t = Time();
		Point3D< Real > p;
		OrientedPoint3D< PointReal > _p;
		while( pointStream->nextPoint( _p ) )//取点
		{
			p = xForm * Point3D< Real >(_p.p);//变换
			for( i=0 ; i<DIMENSION ; i++ )
			{
				if( !cnt || p[i]<min[i] ) 
					min[i] = p[i];//找每个维度的最大最小值
				if( !cnt || p[i]>max[i] ) 
					max[i] = p[i];
			}
			cnt++;//point count
		}

		if( _boundaryType==0 ) _scale = std::max< Real >( max[0]-min[0] , std::max< Real >( max[1]-min[1] , max[2]-min[2] ) ) * 2;//最大一个维度的2倍，
		//在boundaryType等于0时，points are scaled to be in the [0.25,0.75]^3 cube，不明白为什么
		else         _scale = std::max< Real >( max[0]-min[0] , std::max< Real >( max[1]-min[1] , max[2]-min[2] ) );//最大的维度，这为什么要有区别
		_center = ( max+min ) /2;//得到中心点坐标
	}

	_scale *= scaleFactor;//scaleFactor是单位cube与包围当前点云cube的比值
	for( i=0 ; i<DIMENSION ; i++ ) _center[i] -= _scale/2;//利用这个操作，再加上后面的缩放操作，可以把点云中心点调整到(0.5)

	// Update the transformation
	{
		XForm4x4< Real > trans = XForm4x4< Real >::Identity() , scale = XForm4x4< Real >::Identity();
		for( int i=0 ; i<3 ; i++ ) scale(i,i) = (Real)(1./_scale ) , trans(3,i) = -_center[i];
		xForm = scale * trans * xForm;//旋转、平移、缩放
	}

	if( splatDepth>0 )
	{
		double t = Time();
		cnt = 0;
		pointStream->reset();//reset搞啥
		Point3D< Real > p , n;
		OrientedPoint3D< PointReal > _p;
		while( pointStream->nextPoint( _p ) )
		{
			p = xForm * Point3D< Real >(_p.p) , n = xFormN * Point3D< Real >(_p.n);//对点的位置和法向进行变换
			if( !_InBounds(p) ) continue;//把boundary之外的点干掉
			myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
			myWidth = Real(1.0);//这个是总的cube宽度，boundaryType为0时根据上面的scale已经被缩放到[0.25,0.75]范围内
			Real weight=Real( 1. );
			if( useConfidence ) weight = Real( Length(n) );//如果有weight，weight与normal长度有关系，参考Kazhdan的解释
			if( samplesPerNode>0 ) weight /= (Real)samplesPerNode;//参照kazhdan的解释，在每个node的sample大于1时，相当于每个sample贡献值再次分散
			temp = &tree;
			int d=0;
			while( d<splatDepth )//在splatDepth之前，这个点会被反复用于计算point density
			{
				//这里的weight在传递进函数后并不改变，它只是在函数执行前记录normal以及samplePerNode是否对当前计算kernelDensityWeight有影响
				UpdateWeightContribution( kernelDensityWeights , temp , p , neighborKey , weight );//这里kernelDensityWeights貌似只是个指针，还没有内容
				if( !temp->children ) temp->initChildren();//初始化child node，树的子节点初始化到了fullDepth这一次层，当splatDepth大于fullDepth时需要继续生成子节点
				//但此时应该只要生成靠近P一侧的子节点的子树，而不是全部子节点都要生成
				int cIndex=TreeOctNode::CornerIndex( myCenter , p );//判断出当前位置点p在哪一个child node上
				temp = temp->children + cIndex;//移动指针到该child node
				myWidth/=2;//细分，width要更新
				if( cIndex&1 ) myCenter[0] += myWidth/2;//更新child node中心点的位置，重新计算p对于child node的weight和contribution
				else           myCenter[0] -= myWidth/2;
				if( cIndex&2 ) myCenter[1] += myWidth/2;
				else           myCenter[1] -= myWidth/2;
				if( cIndex&4 ) myCenter[2] += myWidth/2;
				else           myCenter[2] -= myWidth/2;
				d++;
			}
			UpdateWeightContribution( kernelDensityWeights , temp , p , neighborKey , weight );//在d==splatDepth时再计算一次
			cnt++;//
		}
	}
	kernelDensityWeights.resize( TreeNodeData::NodeCount , 0 );//这里还要resize是为了???

	std::vector< PointData >& points = pointInfo.data;

	cnt = 0;
	pointStream->reset();
	Point3D< Real > p , n;
	OrientedPoint3D< PointReal > _p;
	while( pointStream->nextPoint( _p ) )//splat normal，准备方程Ax=b的右侧部分
	{
		p = xForm * Point3D< Real >(_p.p) , n = xFormN * Point3D< Real >(_p.n);
		n *= Real(-1.);//normal被翻转了，为啥???
		if( !_InBounds(p) ) continue;
		myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
		myWidth = Real(1.0);
		Real normalLength = Real( Length( n ) );
		if( /*isnan( normalLength ) || !isfinite( normalLength ) || */normalLength<=EPSILON ) continue;
		if( !useConfidence ) n /= normalLength;//如果不用weight的话就把normal单位化

		Real pointWeight = Real(1.f);
		//normalInfo是用来存储splatted normal的，pointWeight是p在其周围node进行splat normal后在调整后的Depth得到的point density，是splat normal的主要权重
		if( samplesPerNode>0 && splatDepth ) pointWeight = SplatPointData( GetPointer( kernelDensityWeights ) , p , n , normalInfo , neighborKey , _minDepth , maxDepth );
		else
		{
			temp = &tree;
			int d=0;
			if( splatDepth )
			{
				while( d<splatDepth )//找到splatDepth处的octree node
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
				pointWeight = (Real)1.0/GetSamplesPerNode( GetPointer( kernelDensityWeights ) , temp , p , neighborKey );
			}
			for( i=0 ; i<DIMENSION ; i++ ) n[i] *= pointWeight;//把pointweight加载到normal上，后面重建会用到
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
			}//找到maxDepth处的octree node
			SplatPointData( temp , p , n , normalInfo , neighborKey );//论文中的splat point过程，要细看
		}
		pointWeightSum += pointWeight;
		if( _constrainValues )//点位置的screen约束，而且是从根节点向下迭代设置的
		{
			Real pointScreeningWeight = useNormalWeights ? Real( normalLength ) : Real(1.f);//设置screen weight
			int d = 0;
			TreeOctNode* temp = &tree;
			myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
			myWidth = Real(1.0);
			while( 1 )
			{
				if( (int)pointInfo.indices.size()<TreeNodeData::NodeCount ) pointInfo.indices.resize( TreeNodeData::NodeCount , -1 );
				int idx = pointInfo.index( temp );

				if( idx==-1 )
				{
					idx = (int)points.size();//点位置的权重约束
					points.push_back( PointData( p*pointScreeningWeight , pointScreeningWeight ) );
					pointInfo.indices[ temp->nodeData.nodeIndex ] = idx;
				}
				else
				{
					points[idx].weight += pointScreeningWeight;//还是按照node来绑定这种point position约束
					points[idx].position += p*pointScreeningWeight;
				}

				int cIndex = TreeOctNode::CornerIndex( myCenter , p );
				if( !temp->children ) break;//没有child node时跳出
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

	if( _boundaryType==0 ) pointWeightSum *= Real(4.);//这是为什么，不仅对pointWeightSum进行调整，还要与constraintWeight放在一起计算
	constraintWeight *= Real( pointWeightSum );
	constraintWeight /= cnt;

	MemoryUsage( );
	if( _constrainValues )//如果要用point position进行screen约束的话
	{
		// Set the average position and scale the weights
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode(node) )
		{
			if( pointInfo.index( node )!=-1 )//检查一个node是否包含待约束的point
			{
				int idx = pointInfo.index( node );
				points[idx].position /= points[idx].weight;//计算weight average
				int e = ( _boundaryType==0 ? node->depth()-1 : node->depth() ) * adaptiveExponent - ( _boundaryType==0 ? maxDepth-1 : maxDepth ) * (adaptiveExponent-1);
				if( e<0 ) points[idx].weight /= Real( 1<<(-e) );
				else      points[idx].weight *= Real( 1<<  e  );
				points[idx].weight *= Real( constraintWeight );//权重调整???
			}
		}
	}
#if FORCE_NEUMANN_FIELD
	if( _boundaryType==1 )
	{	
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode( node ) )
		{
			int d , off[3] , res;
			node->depthAndOffset( d , off );
			res = 1<<d;
			int idx = normalInfo.index( node );
			if( idx<0 ) continue;
			Point3D< Real >& normal = normalInfo.data[ idx ];
			for( int d=0 ; d<3 ; d++ ) if( off[d]==0 || off[d]==res-1 ) normal[d] = 0;//调整Neumann边界的normal约束
		}
	}
#endif // FORCE_NEUMANN_FIELD
	centerWeights.resize( tree.nodes() , 0 );
	kernelDensityWeights.resize( tree.nodes() , 0 );
	// Set the point weights for evaluating the iso-value
	for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode(node) )
	{
		int idx = normalInfo.index( node );
		if( idx<0 ) centerWeights[ node->nodeData.nodeIndex ] = 0;
		else        centerWeights[ node->nodeData.nodeIndex ] = Real( Length( normalInfo.data[ idx ] ) );//所谓的node centerWeight就是node 对应的normalInfo的值
	}
	MemoryUsage();
	{
		std::vector< int > indexMap;
		if( makeComplete ) MakeComplete( &indexMap );
		else ClipTree( normalInfo ) , Finalize( &indexMap );//测试node的子节点是否有normal，没有就置为NULL，然后生成indexMap

		{
			std::vector< int > temp = pointInfo.indices;
			pointInfo.indices.resize( indexMap.size() );
			for( size_t i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<(int)temp.size() ) pointInfo.indices[i] = temp[ indexMap[i] ];//基于上一步利用normal调整octree node的结果来删除某些pointInfo
				else                               pointInfo.indices[i] = -1;
		}
		{
			std::vector< int > temp = normalInfo.indices;
			normalInfo.indices.resize( indexMap.size() );
			for( size_t i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<(int)temp.size() ) normalInfo.indices[i] = temp[ indexMap[i] ];//基于上一步利用normal调整octree node的结果来删除某些normalInfo
				else                               normalInfo.indices[i] = -1;
		}
		{
			std::vector< Real > temp = centerWeights;
			centerWeights.resize( indexMap.size() );
			for( size_t i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<(int)temp.size() ) centerWeights[i] = temp[ indexMap[i] ];//基于上一步利用normal调整octree node的结果来删除某些centerWeights
				else                               centerWeights[i] = (Real)0;
		}
		{
			std::vector< Real > temp = kernelDensityWeights;
			kernelDensityWeights.resize( indexMap.size() );
			for( size_t i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<(int)temp.size() ) kernelDensityWeights[i] = temp[ indexMap[i] ];//基于上一步利用normal调整octree node的结果来删除某些kernelDensityWeights
				else                               kernelDensityWeights[i] = (Real)0;
		}
	}
	return cnt;
}

template< class Real >
template< class Vertex >
void Octree< Real >::GetAdaptiveOctreeGrid(std::vector<Vertex>& vec_octreeGridVertex, std::vector< int >& vec_octreeGridFace)
{
	int x, y, z;
	int iCurVertexSize = 0;
	int iCubeCornerTable [] = {0, 1, 2, 1, 3, 2, 1, 3, 5, 3, 7, 5, 5, 7, 6, 6, 4, 5, 4, 6, 0, 2, 0, 6, 3, 2, 7, 2, 6, 7, 0, 1, 4, 1, 5, 4};
	for( TreeOctNode* temp=tree.nextLeaf() ; temp ; temp=tree.nextLeaf(temp) )
	{
		Point3D<Real> center;
		Real w;
		temp->centerAndWidth( center , w );
		center -= Point3D<Real>(w/2, w/2, w/2);
		for (int i = 0; i < Cube::CORNERS; ++i)
		{
			Cube::FactorCornerIndex(i, x, y, z);
			Point3D<Real> corner = center + Point3D<Real>(w * x, w * y, w * z);
			Vertex newVertex;
			newVertex.point = corner * _scale + _center;
			vec_octreeGridVertex.push_back(newVertex);
		}
		for (int i = 0; i < 36; ++i)
		{
			vec_octreeGridFace.push_back(iCubeCornerTable[i] + iCurVertexSize);
		}
		iCurVertexSize += Cube::CORNERS;
	}
}

template< class Real >
template< class PointReal , class Data , class _Data >
int Octree< Real >::SetTree( OrientedPointStreamWithData< PointReal , Data >* pointStream , int minDepth , int maxDepth , int fullDepth , 
							int splatDepth , Real samplesPerNode , Real scaleFactor ,
							bool useConfidence , bool useNormalWeights , Real constraintWeight , int adaptiveExponent ,
							std::vector< Real >& kernelDensityWeights , SparseNodeData< PointData >& pointInfo , SparseNodeData< Point3D< Real > >& normalInfo , std::vector< Real >& centerWeights ,
							SparseNodeData< ProjectiveData< _Data > >& dataValues ,
							XForm4x4< Real >& xForm , int boundaryType , bool makeComplete )
{
	if( splatDepth<0 ) splatDepth = 0;

	_boundaryType = boundaryType;
	if     ( _boundaryType<0 ) _boundaryType = -1;
	else if( _boundaryType>0 ) _boundaryType =  1;

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
#if 0
	else if( _boundaryType==1 && !_constrainValues ) minDepth = std::max< int >( minDepth , 1 );
#endif

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

	tree.setFullDepth( _fullDepth );

	// Read through once to get the center and scale
	{
		double t = Time();
		Point3D< Real > p;
		OrientedPoint3D< PointReal > _p;
		while( pointStream->nextPoint( _p ) )
		{
			p = xForm * Point3D< Real >(_p.p);
			for( i=0 ; i<DIMENSION ; i++ )
			{
				if( !cnt || p[i]<min[i] ) min[i] = p[i];
				if( !cnt || p[i]>max[i] ) max[i] = p[i];
			}
			cnt++;
		}

		if( _boundaryType==0 ) _scale = std::max< Real >( max[0]-min[0] , std::max< Real >( max[1]-min[1] , max[2]-min[2] ) ) * 2;
		else                   _scale = std::max< Real >( max[0]-min[0] , std::max< Real >( max[1]-min[1] , max[2]-min[2] ) );
		_center = ( max+min ) /2;
	}

	_scale *= scaleFactor;
	for( i=0 ; i<DIMENSION ; i++ ) _center[i] -= _scale/2;

	// Update the transformation
	{
		XForm4x4< Real > trans = XForm4x4< Real >::Identity() , scale = XForm4x4< Real >::Identity();
		for( int i=0 ; i<3 ; i++ ) scale(i,i) = (Real)(1./_scale ) , trans(3,i) = -_center[i];
		xForm = scale * trans * xForm;
	}

	if( splatDepth>0 )
	{
		double t = Time();
		cnt = 0;
		pointStream->reset();
		Point3D< Real > p , n;
		OrientedPoint3D< PointReal > _p;
		while( pointStream->nextPoint( _p ) )
		{
			p = xForm * Point3D< Real >(_p.p) , n = xFormN * Point3D< Real >(_p.n);
			if( !_InBounds(p) ) continue;
			myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
			myWidth = Real(1.0);
			Real weight=Real( 1. );
			if( useConfidence ) weight = Real( Length(n) );
			if( samplesPerNode>0 ) weight /= (Real)samplesPerNode;
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
	kernelDensityWeights.resize( TreeNodeData::NodeCount , 0 );

	std::vector< PointData >& points = pointInfo.data;

	cnt = 0;
	pointStream->reset();
	Point3D< Real > p , n;
	OrientedPoint3D< PointReal > _p;
	Data _d;
	while( pointStream->nextPoint( _p , _d ) )
	{
		p = xForm * Point3D< Real >(_p.p) , n = xFormN * Point3D< Real >(_p.n);
		n *= Real(-1.);
		if( !_InBounds(p) ) continue;
		myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
		myWidth = Real(1.0);
		Real normalLength = Real( Length( n ) );
		if( normalLength!=normalLength || normalLength<=EPSILON ) continue;
		if( !useConfidence ) n /= normalLength;

		Real pointWeight = Real(1.f);
		if( samplesPerNode>0 && splatDepth )
		{
			MultiSplatPointData( GetPointer( kernelDensityWeights ) , p , ProjectiveData< _Data >( _Data( _d ) , (Real)1. ) , dataValues , neighborKey , maxDepth , 2 );
			pointWeight = SplatPointData( GetPointer( kernelDensityWeights ) , p , n , normalInfo , neighborKey , _minDepth , maxDepth , 3 );
		}
		else
		{
			MultiSplatPointData( NullPointer( Real ) , p , ProjectiveData< _Data >( _Data( _d ) , (Real)1. ) , dataValues , neighborKey , maxDepth , 2 );

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
				pointWeight = (Real)1.0/GetSamplesPerNode( GetPointer( kernelDensityWeights ) , temp , p , neighborKey );
			}
			for( i=0 ; i<DIMENSION ; i++ ) n[i] *= pointWeight;
			while( d<maxDepth )
			{
				if( !temp->children ) temp->initChildren();
				int cIndex=TreeOctNode::CornerIndex( myCenter , p );
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
			SplatPointData( temp , p , n , normalInfo , neighborKey );
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
				if( (int)pointInfo.indices.size()<TreeNodeData::NodeCount ) pointInfo.indices.resize( TreeNodeData::NodeCount , -1 );
				int idx = pointInfo.index( temp );

#if POINT_DATA_RES
				if( idx==-1 )
				{
					idx = (int)points.size();
					points.resize( idx+1 );
					pointInfo.indices[ temp->nodeData.nodeIndex ] = idx;
				}
				points[idx].addPoint( p*pointScreeningWeight , myCenter , myWidth , pointScreeningWeight );
#else // !POINT_DATA_RES
				if( idx==-1 )
				{
					idx = (int)points.size();
					points.push_back( PointData( p*pointScreeningWeight , pointScreeningWeight ) );
					pointInfo.indices[ temp->nodeData.nodeIndex ] = idx;
				}
				else
				{
					points[idx].weight += pointScreeningWeight;
					points[idx].position += p*pointScreeningWeight;
				}
#endif // POINT_DATA_RES

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
	if( _constrainValues )
		// Set the average position and scale the weights
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode(node) )
			if( pointInfo.index( node )!=-1 )
			{
				int idx = pointInfo.index( node );
#if POINT_DATA_RES
				for( int c=0 ; c<_PointData::SAMPLES ; c++ )
					if( points[idx].weights[c] )
					{
						points[idx].positions[c] /= points[idx].weights[c];
						int e = ( _boundaryType==0 ? node->depth()-1 : node->depth() ) * adaptiveExponent - ( _boundaryType==0 ? maxDepth-1 : maxDepth ) * (adaptiveExponent-1);
						if( e<0 ) points[idx].weights[c] /= Real( 1<<(-e) );
						else      points[idx].weights[c] *= Real( 1<<  e  );
						points[idx].weights[c] *= Real( constraintWeight );
					}
#else // !POINT_DATA_RES
				points[idx].position /= points[idx].weight;
				int e = ( _boundaryType==0 ? node->depth()-1 : node->depth() ) * adaptiveExponent - ( _boundaryType==0 ? maxDepth-1 : maxDepth ) * (adaptiveExponent-1);
				if( e<0 ) points[idx].weight /= Real( 1<<(-e) );
				else      points[idx].weight *= Real( 1<<  e  );
				points[idx].weight *= Real( constraintWeight );
#endif // POINT_DATA_RES
			}
#if FORCE_NEUMANN_FIELD
	if( _boundaryType==1 )
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode( node ) )
		{
			int d , off[3] , res;
			node->depthAndOffset( d , off );
			res = 1<<d;
			int idx = normalInfo.index( node );
			if( idx<0 ) continue;
			Point3D< Real >& normal = normalInfo.data[ idx ];
			for( int d=0 ; d<3 ; d++ ) if( off[d]==0 || off[d]==res-1 ) normal[d] = 0;
		}
#endif // FORCE_NEUMANN_FIELD
	centerWeights.resize( tree.nodes() , 0 );
	kernelDensityWeights.resize( tree.nodes() , 0 );
	// Set the point weights for evaluating the iso-value
	for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode(node) )
	{
		int idx = normalInfo.index( node );
		if( idx<0 ) centerWeights[ node->nodeData.nodeIndex ] = 0;
		else        centerWeights[ node->nodeData.nodeIndex ] = Real( Length( normalInfo.data[ idx ] ) );
	}
	MemoryUsage();
	{
		std::vector< int > indexMap;
		if( makeComplete ) MakeComplete( &indexMap );
		else ClipTree( normalInfo ) , Finalize( &indexMap );

		{
			std::vector< int > temp = pointInfo.indices;
			pointInfo.indices.resize( indexMap.size() );
			for( size_t i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<(int)temp.size() ) pointInfo.indices[i] = temp[ indexMap[i] ];
				else                               pointInfo.indices[i] = -1;
		}
		{
			std::vector< int > temp = normalInfo.indices;
			normalInfo.indices.resize( indexMap.size() );
			for( size_t i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<(int)temp.size() ) normalInfo.indices[i] = temp[ indexMap[i] ];
				else                               normalInfo.indices[i] = -1;
		}
		{
			std::vector< Real > temp = centerWeights;
			centerWeights.resize( indexMap.size() );
			for( size_t i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<(int)temp.size() ) centerWeights[i] = temp[ indexMap[i] ];
				else                               centerWeights[i] = (Real)0;
		}
		{
			std::vector< Real > temp = kernelDensityWeights;
			kernelDensityWeights.resize( indexMap.size() );
			for( size_t i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<(int)temp.size() ) kernelDensityWeights[i] = temp[ indexMap[i] ];
				else                               kernelDensityWeights[i] = (Real)0;
		}
		{
			std::vector< int > temp = dataValues.indices;
			dataValues.indices.resize( indexMap.size() );
			for( size_t i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<(int)temp.size() ) dataValues.indices[i] = temp[ indexMap[i] ];
				else                               dataValues.indices[i] = -1;
		}
	}
	return cnt;
}

template< class Real >
void Octree< Real >::MakeComplete( std::vector< int >* map )
{
	tree.setFullDepth( tree.maxDepth() );
	refineBoundary( map );
	MemoryUsage();
}
template< class Real >
void Octree< Real >::ClipTree( const SparseNodeData< Point3D< Real > >& normalInfo )//测试node的子节点是否有normal，没有就置为NULL
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

template< class Real >
void Octree< Real >::Finalize( std::vector< int >* map )
{
	int maxDepth = tree.maxDepth( );
	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set( maxDepth );
	for( int d=maxDepth ; d>1 ; d-- )
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode( node ) ) if( node->depth()==d )
		{
			typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node->parent->parent );
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( neighbors.neighbors[i][j][k] && !neighbors.neighbors[i][j][k]->children )
					neighbors.neighbors[i][j][k]->initChildren();
		}
	refineBoundary( map );
}
template< class Real >
double Octree< Real >::GetLaplacian( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent ) const
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
template< class Real >
double Octree< Real >::GetDivergence1( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent , const Point3D< Real >& normal1 ) const
{
	return Point3D< double >::Dot( GetDivergence1( integrator , d , off1 , off2 , childParent ) , normal1 );
}
template< class Real > 
double Octree< Real >::GetDivergence2( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent , const Point3D< Real >& normal2 ) const
{
	return Point3D< double >::Dot( GetDivergence2( integrator , d , off1 , off2 , childParent ) , normal2 );
}
template< class Real >
Point3D< double > Octree< Real >::GetDivergence1( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent ) const
{
	double vv[] =//根据下文提示，这里v应该代表vector field
	{
		//取出在off1与off2定义的BSpline的dot integral结果
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
	return  Point3D< double >( vd[0]*vv[1]*vv[2] , vv[0]*vd[1]*vv[2] , vv[0]*vv[1]*vd[2] );//divergence就是任意坐标轴的梯度乘以其他两个坐标轴的scalar value，然后相加
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
template< class Real > 
Point3D< double > Octree< Real >::GetDivergence2( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent ) const
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

template< class Real >
int Octree< Real >::GetMatrixRowSize( const typename TreeOctNode::Neighbors5& neighbors5 , bool symmetric ) const
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

template< class Real >
int Octree< Real >::SetMatrixRow( const SparseNodeData< PointData >& pointInfo , const typename TreeOctNode::Neighbors5& neighbors5 , Pointer( MatrixEntry< Real > ) row , int offset , const typename BSplineData< 2 >::Integrator& integrator , const Stencil< double , 5 >& stencil , bool symmetric ) const
{
	const std::vector< PointData >& points = pointInfo.data;
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
					if( _node && pointInfo.index( _node )!=-1 )
					{
						const PointData& pData = points[ pointInfo.index( _node ) ];
						Real* _splineValues = splineValues + 3*3*(3*(3*j+k)+l);
						Real weight = pData.weight;
						Point3D< Real > p = pData.position;
						for( int s=0 ; s<3 ; s++ )
						{
							if( idx[0]+j-s>=0 && idx[0]+j-s<((2<<d)-1) ) _splineValues[3*0+s] = Real( _fData.baseBSplines[ idx[0]+j-s][s]( p[0] ) );
							if( idx[1]+k-s>=0 && idx[1]+k-s<((2<<d)-1) ) _splineValues[3*1+s] = Real( _fData.baseBSplines[ idx[1]+k-s][s]( p[1] ) );
							if( idx[2]+l-s>=0 && idx[2]+l-s<((2<<d)-1) ) _splineValues[3*2+s] = Real( _fData.baseBSplines[ idx[2]+l-s][s]( p[2] ) );
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
					if( _node && pointInfo.index( _node )!=-1 )
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
// if( scatter ) normals come from the center node
// else          normals come from the neighbors
template< class Real >
void Octree< Real >::SetDivergenceStencil( int depth , const typename BSplineData< 2 >::Integrator& integrator , Stencil< Point3D< double > , 5 >& stencil , bool scatter ) const
{
	if( depth<2 ) return;//depth太小的情况下都没有足够的Neighbor用来计算
	int center = 1<<(depth-1);//不明白为什么要取center
	int offset[] = { center , center , center };
	for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
	{
		int _offset[] = { x+center-2 , y+center-2 , z+center-2 };//左右各两个，加上中间自己
		//scatter的时候以center node normal为主，因此offset用了vector field，而Neighbor用了gradient进行dot求解divergence，用了divergence1，具体区别见函数内部
		if( scatter ) stencil.values[x][y][z] = GetDivergence1( integrator , depth , offset , _offset , false );
		//反之，以Neighbor node normal为主，因此_offset用了vector field，center node用了gradient，用了divergence2
		else          stencil.values[x][y][z] = GetDivergence2( integrator , depth , offset , _offset , false );//GetDivergence1是vd，GetDivergence2返回的是dv
	}
}
template< class Real >
void Octree< Real >::SetDivergenceStencils( int depth , const typename BSplineData< 2 >::Integrator& integrator , Stencil< Point3D< double > ,  5 > stencils[2][2][2] , bool scatter ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ )
	{
		int offset[] = { center+i , center+j , center+k };
		for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
		{
			int _offset[] = { x-2+center/2 , y-2+center/2 , z-2+center/2 };//这个函数是用来计算child node和parent node之间的dot，所以_offset要在center/2周围变动
			if( scatter ) stencils[i][j][k].values[x][y][z] = GetDivergence1( integrator , depth , offset , _offset , true );
			else          stencils[i][j][k].values[x][y][z] = GetDivergence2( integrator , depth , offset , _offset , true );//childParent flag被置为true
		}
	}
}
template< class Real >
void Octree< Real >::SetLaplacianStencil( int depth , const typename BSplineData< 2 >::Integrator& integrator , Stencil< double , 5 >& stencil ) const
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
template< class Real >
void Octree< Real >::SetLaplacianStencils( int depth , const typename BSplineData< 2 >::Integrator& integrator , Stencil< double , 5 > stencils[2][2][2] ) const
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
template< class Real >
void Octree< Real >::SetCenterEvaluationStencil( const typename BSplineData< 2 >::template CenterEvaluator< 1 >& evaluator , int depth , Stencil< double , 3 >& stencil ) const
{
	//这种stencil属于正常情况下的模板，左右位置的neighbor都存在
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
	{
		int off[] = { center+x-1 , center+y-1 , center+z-1 };//左右各偏移一个位置
		//每个维度上都算了一遍，然后相乘，center与neighbor属于同一层，而且求的是value乘积不是gradient value乘积
		stencil.values[x][y][z] = Real( evaluator.value( depth , center , off[0] , false , false ) * evaluator.value( depth , center , off[1] , false , false ) * evaluator.value( depth , center , off[2] , false , false ) );
	}
}
template< class Real >
void Octree< Real >::SetCenterEvaluationStencils( const typename BSplineData< 2 >::template CenterEvaluator< 1 >& evaluator , int depth , Stencil< double , 3 > stencils[8] ) const
{
	//这个函数产生的child parent的value stencil，根据子节点与parent node位置的不同分为八种情况
	if( depth<3 ) return;
	int center = 1<<(depth-1);//depth与center是同一层
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )//center由上面的式子赋值，一定是偶数，位于其parent的左侧，这里加0或1就是要遍历同一个parent下的所有child node，将这不同的八种情况都作为stencil
	{
		int idx[] = { center+cx , center+cy , center+cz };
		for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
		{
			int off[] = { center/2+x-1 , center/2+y-1 , center/2+z-1 };//off是父节点那一层的neighbor，与center的parent node左右偏离各一个node
			stencils[Cube::CornerIndex( cx , cy , cz ) ].values[x][y][z] = Real( evaluator.value( depth , idx[0] , off[0] , false , true ) * evaluator.value( depth , idx[1] , off[1] , false , true ) * evaluator.value( depth , idx[2] , off[2] , false , true ) );
		}
	}
}
template< class Real >
void Octree< Real >::SetCornerEvaluationStencil( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< double , 3 > stencil[8] ) const
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
template< class Real >
void Octree< Real >::SetCornerEvaluationStencils( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< double , 3 > stencils[8][8] ) const
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
template< class Real >
void Octree< Real >::SetCornerNormalEvaluationStencil( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 3 > stencil[8] ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );
		for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
		{
			int off[] = { center+x-1 , center+y-1 , center+z-1 };
			double v [] = { evaluator.value( depth , center , cx , off[0] , false , false ) , evaluator.value( depth , center , cy , off[1] , false , false ) , evaluator.value( depth , center , cz , off[2] , false , false ) };
			double dv[] = { evaluator.value( depth , center , cx , off[0] , true  , false ) , evaluator.value( depth , center , cy , off[1] , true  , false ) , evaluator.value( depth , center , cz , off[2] , true  , false ) };
			stencil[c].values[x][y][z] = Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] );
		}
	}
}
template< class Real >
void Octree< Real >::SetCornerNormalEvaluationStencils( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 3 > stencils[8][8] ) const
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
			for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
			{
				int off[] = { center/2+x-1 , center/2+y-1 , center/2+z-1 };
				double v [] = { evaluator.value( depth , idx[0] , cx , off[0] , false , true ) , evaluator.value( depth , idx[1] , cy , off[1] , false , true ) , evaluator.value( depth , idx[2] , cz , off[2] , false , true ) };
				double dv[] = { evaluator.value( depth , idx[0] , cx , off[0] , true  , true ) , evaluator.value( depth , idx[1] , cy , off[1] , true  , true ) , evaluator.value( depth , idx[2] , cz , off[2] , true  , true ) };
				stencils[c][_c].values[x][y][z] = Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] );
			}
		}
	}
}
template< class Real >
void Octree< Real >::SetCornerNormalEvaluationStencil( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 5 > stencil[8] ) const
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
template< class Real >
void Octree< Real >::SetCornerNormalEvaluationStencils( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 5 > stencils[8][8] ) const
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

template< class Real >
void Octree< Real >::UpdateCoarserSupportBounds( const TreeOctNode* node , int& startX , int& endX , int& startY , int& endY , int& startZ , int& endZ )
{
	//这个函数的意思在于，当前node处于parent位置不同，那么在parent Neighbor层面，能与其BSpline产生重合的parent Neighbor也不同
	//从一维角度来看，如果node处于parent左侧，即x=0，那么能对其BSpline产生影响的(child node层面)只有左右各两个node，而这四个node分别从属于
	//parent neighbor从0到3四个node中，因此要修改endX=4(后面使用时取了小于号，而不是小于等于)；同理，如果x=1，则四个overlapping child node从属于
	//parent neighbor从1到4四个node中，因此要修改StartX=1
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
template< class Real >
void Octree< Real >::UpdateConstraintsFromFiner( const typename BSplineData< 2 >::Integrator& integrator , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) fineSolution , Pointer( Real ) coarseConstraints ) const
{
	if( depth<=_minDepth ) return;
	Stencil< double , 5 > stencils[2][2][2];
	// Get the stencil describing the Laplacian relating coefficients @(depth) with coefficients @(depth-1)
	SetLaplacianStencils( depth , integrator , stencils );
	size_t start = sNodes.nodeCount[depth] , end = sNodes.nodeCount[depth+1] , range = end-start;
	int lStart = sNodes.nodeCount[depth-1];
	memset( coarseConstraints , 0 , sizeof(Real)*(sNodes.nodeCount[depth]-sNodes.nodeCount[depth-1]) );

	// Iterate over the nodes @( depth )
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth-1 );
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		TreeOctNode* node = sNodes.treeNodes[i];

		bool insetSupported = _boundaryType!=0 || _IsInsetSupported( node );

		// Offset the coarser constraints using the solution from the current resolutions.
		int x , y , z , c;
		c = int( node - node->parent->children );
		Cube::FactorCornerIndex( c , x , y , z );
		if( insetSupported )
		{
			typename TreeOctNode::Neighbors5 pNeighbors5;
			neighborKey.getNeighbors( node->parent , pNeighbors5 );
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

template< class Real >
void Octree< Real >::UpdateConstraintsFromCoarser( const SparseNodeData< PointData >& pointInfo , const typename TreeOctNode::Neighbors5& neighbors5 , const typename TreeOctNode::Neighbors5& pNeighbors5 , TreeOctNode* node , Pointer( Real ) constraints , ConstPointer( Real ) metSolution , const typename BSplineData< 2 >::Integrator& integrator , const Stencil< double , 5 >& lapStencil ) const
{
	const std::vector< PointData >& points = pointInfo.data;
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
			if( neighbors5.neighbors[x][y][z] && pointInfo.index( neighbors5.neighbors[x][y][z] )!=-1 )
			{
				const PointData& pData = points[ pointInfo.index( neighbors5.neighbors[x][y][z] ) ];
				Point3D< Real > p = pData.position;
				constraint += 
					_fData.baseBSplines[idx[0]][x-1]( p[0] ) *
					_fData.baseBSplines[idx[1]][y-1]( p[1] ) *
					_fData.baseBSplines[idx[2]][z-1]( p[2] ) * 
					pData.weightedCoarserDValue;
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
template< class Real >
template< class C >
void Octree< Real >::DownSample( int depth , const SortedTreeNodes& sNodes , ConstPointer( C ) fineConstraints , Pointer( C ) coarseConstraints ) const
{
	if( depth==0 ) return;
	double cornerValue;
	if     ( _boundaryType==-1 ) cornerValue = 0.50;
	else if( _boundaryType== 1 ) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth );//NeighborKey3代表从0到Depth的Neighbor3
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )//第d层sNodes的起始和结束index
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		int d , off[3];
		UpSampleData usData[3];
		sNodes.treeNodes[i]->depthAndOffset( d , off );
		for( int dd=0 ; dd<3 ; dd++ )//x/y/z三个dimension，按照如下权重指定downSample时的Weight
		{
			if     ( off[dd]  ==0          ) usData[dd] = UpSampleData( 1 , cornerValue , 0.00 );//左边界
			else if( off[dd]+1==(1<<depth) ) usData[dd] = UpSampleData( 0 , 0.00 , cornerValue );//右边界
			else if( off[dd]%2             ) usData[dd] = UpSampleData( 1 , 0.75 , 0.25 );//奇数项，在父节点中偏右，为什么在不是边界的情况下也只考虑父节点的两个Neighbor
			else                             usData[dd] = UpSampleData( 0 , 0.25 , 0.75 );//偶数项，在父节点中偏左
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors( sNodes.treeNodes[i]->parent );//父节点的Neighbor
		C c = fineConstraints[ i-sNodes.nodeCount[depth] ];//fineConstraints在函数外已经做了指针偏移，现在调用时偏移量为i-sNodes.nodeCount[depth]
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
					TreeOctNode* pNode = neighbors.neighbors[_ii][_jj][_kk];//找出父节点的Neighbor
					if( pNode )
#pragma omp atomic
						//coarseConstraints在函数外已经做了指针偏移，现在调用时偏移量为pNode->nodeData.nodeIndex-sNodes.nodeCount[depth-1]
						coarseConstraints[ pNode->nodeData.nodeIndex-sNodes.nodeCount[depth-1] ] += C( cxy*usData[2].v[kk] );//给父节点及其Neighbor加上weighted constraint
				}
			}
		}
	}
}
template< class Real >
template< class C >
void Octree< Real >::UpSample( int depth , const SortedTreeNodes& sNodes , ConstPointer( C ) coarseCoefficients , Pointer( C ) fineCoefficients ) const
{
	//原理跟DownSample差不多，都是用child parent之间的位置关系进行加权
	double cornerValue;
	if     ( _boundaryType==-1 ) cornerValue = 0.50;
	else if( _boundaryType== 1 ) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	if( depth<=_minDepth ) return;

	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth-1 );
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		bool isInterior = true;
		TreeOctNode* node = sNodes.treeNodes[i];
		int d , off[3];
		UpSampleData usData[3];
		node->depthAndOffset( d , off );
		for( int d=0 ; d<3 ; d++ )
		{
			if     ( off[d]  ==0          ) usData[d] = UpSampleData( 1 , cornerValue , 0.00 ) , isInterior = false;//左边界
			else if( off[d]+1==(1<<depth) ) usData[d] = UpSampleData( 0 , 0.00 , cornerValue ) , isInterior = false;//右边界
			else if( off[d]%2             ) usData[d] = UpSampleData( 1 , 0.75 , 0.25 );//奇数
			else                            usData[d] = UpSampleData( 0 , 0.25 , 0.75 );//偶数
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
					if( node )
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
template< class Real >
void Octree< Real >::SetPointValuesFromCoarser( SparseNodeData< PointData >& pointInfo , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) coarseCoefficients )
{
	std::vector< PointData >& points = pointInfo.data;
	// For every node at the current depth
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth );
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		int pIdx = pointInfo.index( sNodes.treeNodes[i] );
		if( pIdx!=-1 )
		{
			neighborKey.getNeighbors( sNodes.treeNodes[i] );
			points[ pIdx ].weightedCoarserDValue = (Real)( _CoarserFunctionValue( points[pIdx] , neighborKey , sNodes.treeNodes[i] , coarseCoefficients-_sNodes.nodeCount[depth-1] ) - 0.5 ) * points[pIdx].weight;
		}
	}
}
template< class Real >
Real Octree< Real >::_CoarserFunctionValue( const PointData& pointData , const typename TreeOctNode::NeighborKey3& neighborKey , const TreeOctNode* pointNode , ConstPointer( Real ) coarseCoefficients ) const
{
	double pointValue = 0;
	int depth = pointNode->depth();

	if( depth<=_minDepth ) return Real(0.);

	Point3D< Real > p = pointData.position;

	// Iterate over all basis functions that overlap the point at the coarser resolutions
	{
		int d , _idx[3];
		const typename TreeOctNode::Neighbors3& neighbors = neighborKey.neighbors[depth-1];
		neighbors.neighbors[1][1][1]->depthAndOffset( d , _idx );
		_idx[0] = BinaryNode::CenterIndex( d , _idx[0]-1 );
		_idx[1] = BinaryNode::CenterIndex( d , _idx[1]-1 );
		_idx[2] = BinaryNode::CenterIndex( d , _idx[2]-1 );

		for( int j=0 ; j<3 ; j++ )
		{
			double xValue = 0;
			if( _idx[0]+j>=0 && _idx[0]+j<((1<<depth)-1) ) xValue = _fData.baseBSplines[ _idx[0]+j ][2-j]( p[0] );
			else continue;
			for( int k=0 ; k<3 ; k++ )
			{
				double xyValue = 0;
				if( _idx[1]+k>=0 && _idx[1]+k<((1<<depth)-1) ) xyValue = xValue * _fData.baseBSplines[ _idx[1]+k ][2-k]( p[1] );
				else continue;
				double _pointValue = 0;
				for( int l=0 ; l<3 ; l++ )
				{
					const TreeOctNode* basisNode = neighbors.neighbors[j][k][l];
					if( basisNode && basisNode->nodeData.nodeIndex>=0 && _idx[2]+l>=0 && _idx[2]+l<((1<<depth)-1) )
						_pointValue += _fData.baseBSplines[ _idx[2]+l ][2-l]( p[2] ) * double( coarseCoefficients[basisNode->nodeData.nodeIndex] );
				}
				pointValue += _pointValue * xyValue;
			}
		}
	}
	return Real( pointValue );
}
template< class Real >
void Octree< Real >::SetPointConstraintsFromFiner( const SparseNodeData< PointData >& pointInfo , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) finerCoefficients , Pointer( Real ) coarserConstraints ) const
{
	const std::vector< PointData >& points = pointInfo.data;
	// Note: We can't iterate over the finer point nodes as the point weights might be
	// scaled incorrectly, due to the adaptive exponent. So instead, we will iterate
	// over the coarser nodes and evaluate the finer solution at the associated points.
	if( !depth ) return;
	size_t start = sNodes.nodeCount[depth-1] , end = sNodes.nodeCount[depth] , range = end-start;
	memset( coarserConstraints , 0 , sizeof( Real ) * ( sNodes.nodeCount[depth]-sNodes.nodeCount[depth-1] ) );
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth-1 );
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth-1] ; i<sNodes.nodeCount[depth] ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		int pIdx = pointInfo.index( sNodes.treeNodes[i] );
		if( pIdx!=-1 )
		{
			typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors( sNodes.treeNodes[i] );
			// Evaluate the solution @( depth ) at the current point @( depth-1 )
			{
				Real finerPointDValue = (Real)( _FinerFunctionValue( points[pIdx] , neighborKey , sNodes.treeNodes[i] , finerCoefficients-sNodes.nodeCount[depth] ) - 0.5 ) * points[pIdx].weight;
				Point3D< Real > p = points[ pIdx ].position;
				// Update constraints for all nodes @( depth-1 ) that overlap the point
				int d , idx[3];
				neighbors.neighbors[1][1][1]->depthAndOffset( d, idx );
				// Set the (offset) index to the top-left-front corner of the 3x3x3 block of b-splines
				// overlapping the point.
				idx[0] = BinaryNode::CenterIndex( d , idx[0]-1 );
				idx[1] = BinaryNode::CenterIndex( d , idx[1]-1 );
				idx[2] = BinaryNode::CenterIndex( d , idx[2]-1 );
				for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
					if( neighbors.neighbors[x][y][z] )
					{
#pragma omp atomic
						coarserConstraints[ neighbors.neighbors[x][y][z]->nodeData.nodeIndex - sNodes.nodeCount[depth-1] ] +=
							Real(
							_fData.baseBSplines[idx[0]+x][2-x]( p[0] ) *
							_fData.baseBSplines[idx[1]+y][2-y]( p[1] ) *
							_fData.baseBSplines[idx[2]+z][2-z]( p[2] ) * 
							finerPointDValue
							);
					}
			}
		}
	}
}
template< class Real >
Real Octree< Real >::_FinerFunctionValue( const PointData& pointData , const typename TreeOctNode::NeighborKey3& neighborKey , const TreeOctNode* pointNode , ConstPointer( Real ) finerCoefficients ) const
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
	idx[0] = BinaryNode::CenterIndex( d , idx[0]-1 );
	idx[1] = BinaryNode::CenterIndex( d , idx[1]-1 );
	idx[2] = BinaryNode::CenterIndex( d , idx[2]-1 );

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
	return Real( pointValue );
}
template< class Real >
int Octree< Real >::GetSliceMatrixAndUpdateConstraints( const SparseNodeData< PointData >& pointInfo , SparseMatrix< Real >& matrix , Pointer( Real ) constraints , const typename BSplineData< 2 >::Integrator& integrator , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) metSolution , bool coarseToFine , int nStart , int nEnd )
{
	size_t range = nEnd-nStart;
	Stencil< double , 5 > stencil , stencils[2][2][2];
	SetLaplacianStencil ( depth , integrator , stencil );
	SetLaplacianStencils( depth , integrator , stencils );
	matrix.Resize( (int)range );
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<(int)range ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		TreeOctNode* node = sNodes.treeNodes[i+nStart];
		// Get the matrix row size
		bool insetSupported = _boundaryType!=0 || _IsInsetSupported( node );
		typename TreeOctNode::Neighbors5 neighbors5;
		if( insetSupported ) neighborKey.getNeighbors( node , neighbors5 );
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
				neighborKey.getNeighbors( node->parent , pNeighbors5 );
				UpdateConstraintsFromCoarser( pointInfo , neighbors5 , pNeighbors5 , node , constraints , metSolution , integrator , stencils[x][y][z] );
			}
		}
	}
	return 1;
}
template< class Real >
int Octree< Real >::GetMatrixAndUpdateConstraints( const SparseNodeData< PointData >& pointInfo , SparseSymmetricMatrix< Real >& matrix , Pointer( Real ) constraints , const typename BSplineData< 2 >::Integrator& integrator , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) metSolution , bool coarseToFine )
{
	size_t start = sNodes.nodeCount[depth] , end = sNodes.nodeCount[depth+1] , range = end-start;
	Stencil< double , 5 > stencil , stencils[2][2][2];
	SetLaplacianStencil ( depth , integrator , stencil );
	SetLaplacianStencils( depth , integrator , stencils );
	matrix.Resize( (int)range );
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<(int)range ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		TreeOctNode* node = sNodes.treeNodes[i+start];
		// Get the matrix row size
		bool insetSupported = _boundaryType!=0 || _IsInsetSupported( node );
		typename TreeOctNode::Neighbors5 neighbors5;
		if( insetSupported ) neighborKey.getNeighbors( node , neighbors5 );
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
				neighborKey.getNeighbors( node->parent , pNeighbors5 );
				UpdateConstraintsFromCoarser( pointInfo , neighbors5 , pNeighbors5 , node , constraints , metSolution , integrator , stencils[x][y][z] );
			}
		}
	}
	return 1;
}

template< class Real >
Pointer( Real ) Octree< Real >::SolveSystem( SparseNodeData< PointData >& pointInfo , Pointer( Real ) constraints , bool showResidual , int iters , int maxSolveDepth , int cgDepth , double accuracy )
{
	int iter=0;
	typename BSplineData< 2 >::Integrator integrator;
	_fData.setIntegrator( integrator , _boundaryType==0 );
	iters = std::max< int >( 0 , iters );
	if( _boundaryType==0 ) maxSolveDepth++ , cgDepth++;

	Pointer( Real ) solution = AllocPointer< Real >( _sNodes.nodeCount[_sNodes.maxDepth] );//Octree node总数
	memset( solution , 0 , sizeof(Real)*_sNodes.nodeCount[_sNodes.maxDepth] );

	solution[0] = 0;

	std::vector< Real > metSolution( _sNodes.nodeCount[ _sNodes.maxDepth-1 ] , 0 );
	for( int d=_minDepth ; d<_sNodes.maxDepth ; d++ )
	{
		DumpOutput( "Depth[%d/%d]: %d\n" , _boundaryType==0 ? d-1 : d , _boundaryType==0 ? _sNodes.maxDepth-2 : _sNodes.maxDepth-1 , _sNodes.nodeCount[d+1]-_sNodes.nodeCount[d] );
		if( d==_minDepth )//这一层求解用了conjugate-gradient
			_SolveSystemCG( pointInfo , d , integrator , _sNodes , solution , constraints , GetPointer( metSolution ) , _sNodes.nodeCount[_minDepth+1]-_sNodes.nodeCount[_minDepth] , true , showResidual, NULL , NULL , NULL );//minDepth层的node数目
		else
		{	//这些求解用到了block gaussian seidel iteration
			if( d>cgDepth ) iter += _SolveSystemGS( pointInfo , d , integrator , _sNodes , solution , constraints , GetPointer( metSolution ) , d>maxSolveDepth ? 0 : iters , true , showResidual , NULL , NULL , NULL );
			else            iter += _SolveSystemCG( pointInfo , d , integrator , _sNodes , solution , constraints , GetPointer( metSolution ) , d>maxSolveDepth ? 0 : iters , true , showResidual , NULL , NULL , NULL , accuracy );//如果没大于cgDepth，求解还是用conjugate-gradient
		}
	}

	return solution;
}
template< class Real >
void Octree< Real >::_setMultiColorIndices( int start , int end , std::vector< std::vector< int > >& indices ) const
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
template< class Real >
int Octree< Real >::_SolveSystemGS( SparseNodeData< PointData >& pointInfo , int depth , const typename BSplineData< 2 >::Integrator& integrator , const SortedTreeNodes& sNodes , Pointer( Real ) solution , Pointer( Real ) constraints , Pointer( Real ) metSolutionConstraints , int iters , bool coarseToFine , bool showResidual , double* bNorm2 , double* inRNorm2 , double* outRNorm2 , bool forceSilent )
{
	Pointer( Real ) metSolution = NullPointer( Real );
	Pointer( Real ) metConstraints = NullPointer( Real );
	if( coarseToFine ) metSolution    = metSolutionConstraints;	// This stores the up-sampled solution up to depth-2
	else               metConstraints = metSolutionConstraints; // This stores the down-sampled constraints up to depth

	double _maxMemoryUsage = maxMemoryUsage;
	maxMemoryUsage = 0;
	Vector< Real > X , B;
	int slices = 1<<depth;//binary node在depth下的划分数目
	double systemTime=0. , solveTime=0. , updateTime=0. ,  evaluateTime = 0.;
	std::vector< int > offsets( slices+1 , 0 );
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )//当前层求解的node index起始和截止
	{
		int d , off[3];
		sNodes.treeNodes[i]->depthAndOffset( d , off );
		offsets[ off[2] ]++;//在z轴统计
	}
	for( int i=1 ; i<slices ; i++ )  offsets[i] += offsets[i-1];//叠加
	for( int i=slices ; i>=1 ; i-- ) offsets[i]  = offsets[i-1];//offset[i+1] - offset[i]代表node的z轴坐标等于i时在depth层的node size
	offsets[0] = 0;

	X.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );//分配x和b还是按照depth层的节点数目
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
	for( int i=_sNodes.nodeCount[depth] ; i<_sNodes.nodeCount[depth+1] ; i++ ) X[ i-_sNodes.nodeCount[depth] ] = solution[i];//用前一步求解结果进行迭代
	double bNorm=0 , inRNorm=0 , outRNorm=0;
	if( depth>=_minDepth )
	{
		int frontOffset = ( showResidual || inRNorm2 ) ? 2 : 0;
		int backOffset = ( showResidual || outRNorm2 ) ? 2 : 0;
		int solveSlices = std::min< int >( 2*iters-1 , slices ) , matrixSlices = std::max< int >( 1 , std::min< int >( solveSlices+frontOffset+backOffset , slices ) );
		std::vector< SparseMatrix< Real > > _M( matrixSlices );
		std::vector< std::vector< std::vector< int > > > __mcIndices( std::max< int >( 0 , solveSlices ) );

		int dir = coarseToFine ? -1 : 1 , start = coarseToFine ? slices-1 : 0 , end = coarseToFine ? -1 : slices;
		for( int frontSlice=start-frontOffset*dir , backSlice = frontSlice-2*(iters-1)*dir ; backSlice!=end+backOffset*dir ; frontSlice+=dir , backSlice+=dir )
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
			for( int slice=frontSlice ; slice*dir>=backSlice*dir ; slice-=2*dir )
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
template< class Real >
int Octree< Real >::_SolveSystemCG( SparseNodeData< PointData >& pointInfo , int depth , const typename BSplineData< 2 >::Integrator& integrator , const SortedTreeNodes& sNodes , Pointer( Real ) solution , Pointer( Real ) constraints , Pointer( Real ) metSolutionConstraints , int iters , bool coarseToFine , bool showResidual , double* bNorm2 , double* inRNorm2 , double* outRNorm2 , double accuracy )
{
	Pointer( Real ) metSolution = NullPointer( Real );
	Pointer( Real ) metConstraints = NullPointer( Real );
	if( coarseToFine ) metSolution    = metSolutionConstraints;	// This stores the up-sampled solution up to depth-2
	else               metConstraints = metSolutionConstraints; // This stores the down-sampled constraints up to depth
	double _maxMemoryUsage = maxMemoryUsage;
	maxMemoryUsage = 0;
	int iter = 0;
	Vector< Real > X , B;
	SparseSymmetricMatrix< Real > M;
	double systemTime=0. , solveTime=0. , updateTime=0. ,  evaluateTime = 0.;
	X.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );//当前层的node数目
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
		if( coarseToFine ) GetMatrixAndUpdateConstraints( pointInfo , M , constraints , integrator , depth , sNodes , metSolution         , true  );
		else               GetMatrixAndUpdateConstraints( pointInfo , M , constraints , integrator , depth , sNodes , NullPointer( Real ) , false );
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
template< class Real >
int Octree< Real >::HasNormals( TreeOctNode* node , const SparseNodeData< Point3D< Real > >& normalInfo )
{
	int idx = normalInfo.index( node );
	if( idx>=0 )
	{
		const Point3D< Real >& normal = normalInfo.data[ idx ];
		if( normal[0]!=0 || normal[1]!=0 || normal[2]!=0 ) return 1;
	}
	if( node->children ) for( int i=0 ; i<Cube::CORNERS ; i++ ) if( HasNormals( &node->children[i] , normalInfo ) ) return 1;
	return 0;
}
template< class Real >
Pointer( Real ) Octree< Real >::SetLaplacianConstraints( const SparseNodeData< Point3D< Real > >& normalInfo )
{
	// To set the Laplacian constraints, we iterate over the
	// splatted normals and compute the dot-product of the
	// divergence of the normal field with all the basis functions.
	// Within the same depth: set directly as a gather
	// Coarser depths 
	//constraint作为最后返回值，一共经历了三次修改，第一次是同层之间的overlapping node dot，第二次是用finer node downsample coarser node，
	//第三次是用所有的coarser node计算 finer node；第一次用了center node的同层neighbor node normal和gradient，第二次计算时用了center node的normal和
	//gradient，第三次用了upsample的normal field
	typename BSplineData< 2 >::Integrator integrator;
	_fData.setIntegrator( integrator , _boundaryType==0 );
	int maxDepth = _sNodes.maxDepth-1;//_sNodes的maxDepth比Octree大一层，因此局部变量maxDepth代表了Octree的实际深度
	Point3D< Real > zeroPoint;
	zeroPoint[0] = zeroPoint[1] = zeroPoint[2] = 0;
	Pointer( Real ) constraints = AllocPointer< Real >( _sNodes.nodeCount[_sNodes.maxDepth] );//_sNodes.nodeCount[_sNodes.maxDepth]代表Octree中所有node的数目
	if( !constraints ) fprintf( stderr , "[ERROR] Failed to allocate constraints: %d * %zu\n" , _sNodes.nodeCount[_sNodes.maxDepth] , sizeof( Real ) ) , exit( 0 );
	memset( constraints , 0 , sizeof(Real)*_sNodes.nodeCount[_sNodes.maxDepth] );
	//比上面的constraint少了一层的node，这里maxDepth = _sNodes.maxDepth-1，局部变量
	Pointer( Real ) _constraints = AllocPointer< Real >( _sNodes.nodeCount[maxDepth] );//表示了除最后一层以外其余层node数目之和
	if( !_constraints ) fprintf( stderr , "[ERROR] Failed to allocate _constraints: %d * %zu\n" , _sNodes.nodeCount[maxDepth] , sizeof( Real ) ) , exit( 0 );
	memset( _constraints , 0 , sizeof(Real)*_sNodes.nodeCount[maxDepth] );
	MemoryUsage();

	for( int d=maxDepth ; d>=(_boundaryType==0?2:0) ; d-- )//maxDepth = _sNodes.maxDepth -1
	{
		//sorted node在某一depth下的node Count是由nodeCount[d]和nodeCount[d+1]决定的，nodeIndex应该是在splatting normal时候分配的
		//这个offset就是d-1层nodeIndex的开始
		int offset = d>0 ? _sNodes.treeNodes[ _sNodes.nodeCount[d-1] ]->nodeData.nodeIndex : 0;
		Stencil< Point3D< double > , 5 > stencil , stencils[2][2][2];//Stencil功能类似于version 1中的dotTable, dDotTable，Neighbor=5因为任意BSpline有五个BSpline与之重合
		SetDivergenceStencil ( d , integrator , stencil , false );//scatter为false，以d层center node的Neighbor为主，准备利用其normal field赋值vector field
		SetDivergenceStencils( d , integrator , stencils , true );//以d层center node为主，准备利用center node的normal field赋值vector field

		std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _fData.depth );//neighborKeys记录了从0到depth的neighbor3
#pragma omp parallel for num_threads( threads )
		//这是第d层nodeCount的开始和结束，当d=maxDepth=_sNodes.maxDepth-1时，是_sNodes最后一层的node，因此这个循环包含了Octree最深一层的node
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];//线程变量
			TreeOctNode* node = _sNodes.treeNodes[i];
			int startX=0 , endX=5 , startY=0 , endY=5 , startZ=0 , endZ=5;
			int depth = node->depth();
			typename TreeOctNode::Neighbors5 neighbors5;
			neighborKey.getNeighbors( node , neighbors5 );//通过回溯node的parent节点，找出其三环邻域内的所有节点

			bool isInterior , isInterior2;
			{
				int d , off[3];
				node->depthAndOffset( d , off );
				int o = _boundaryType==0 ? (1<<(d-2)) : 0;
				int mn = 2+o , mx = (1<<d)-2-o;
				//不在两边，既不在小于2的位置，也不在最后两个位置，其实mx=res-2-o
				isInterior  = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );//在这个范围内可以保证左右两侧的overlapping BSpline都存在
				mn += 2 , mx -= 2;//又向内缩进2个node，用在了计算与parent neighbor的重叠关系中
				isInterior2 = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
			}
			int cx , cy , cz;
			if( d )
			{
				int c = int( node - node->parent->children );//node在其所有兄弟节点中的偏移
				Cube::FactorCornerIndex( c , cx , cy , cz );//判断出corner index
			}
			else cx = cy = cz = 0;
			Stencil< Point3D< double > , 5 >& _stencil = stencils[cx][cy][cz];//根据在父节点下充当子节点位置的不同取stencils
			int d , off[3];
			node->depthAndOffset( d , off );
			// Set constraints from current depth
			// Gather the constraints from the vector-field at _node into the constraint stored with node
			{
				//在考虑同一层node间的重合关系时，直接对constraint变量进行赋值，而不是_constraint
				if( isInterior )
				{
					for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
					{
						const TreeOctNode* _node = neighbors5.neighbors[x][y][z];//取三环邻域node
						if( _node )
						{
							int _idx = normalInfo.index( _node );//取出normal index
							//因为stencil取值时利用Neighbor node的normal field为目标，这里的normalInfo自然来源于邻域node
							if( _idx>=0 ) constraints[ node->nodeData.nodeIndex ] += Point3D< Real >::Dot( stencil.values[x][y][z] , normalInfo.data[ _idx ] );
						}
					}
				}
				else
				{	
					for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
					{
						const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						if( _node )
						{
							int _idx = normalInfo.index( _node );
							if( _idx>=0 )
							{
								int _d , _off[3];
								_node->depthAndOffset( _d , _off );
								//GetDivergence2是off的gradient乘以_off的vector field，然后再点乘_off的normal，跟上面的计算过程一样，只是没用到标准的stencil
								constraints[ node->nodeData.nodeIndex ] += Real( GetDivergence2( integrator , d , off , _off , false , normalInfo.data[ _idx ] ) );
							}
						}
					}
				}
				//neighbors[2][2][2]就是current node，这是根据current node在parent node中的位置修改startX等的值，因为parent neighbor[0-4]并不都与current node重合
				UpdateCoarserSupportBounds( neighbors5.neighbors[2][2][2] , startX , endX , startY  , endY , startZ , endZ );
			}
			int idx = normalInfo.index( node );
			if( idx<0 ) continue;
			const Point3D< Real >& normal = normalInfo.data[ idx ];
			if( normal[0]==0 && normal[1]==0 && normal[2]==0 ) continue;

			// Set the constraints for the parents
			if( depth>_minDepth )
			{
				neighborKey.getNeighbors( node->parent , neighbors5 );//取出了父节点的neighbor

				for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
				{
					if( neighbors5.neighbors[x][y][z] )
					{
						TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						Real c;
						if( isInterior2 )
						{
							Point3D< double >& div = _stencil.values[x][y][z];//根据在父节点下位置的不同取出的stencil
							c = Real( div[0] * normal[0] + div[1] * normal[1] + div[2] * normal[2] );
						}
						else
						{
							int _d , _off[3];
							_node->depthAndOffset( _d , _off );
							//off的vector field和_off的gradient vector field，由于是child node和parent node的dot，childParent flag被置为true
							c = Real( GetDivergence1( integrator , d , off , _off , true , normal ) );
						}
#pragma omp atomic
						_constraints[ _node->nodeData.nodeIndex ] += c;//修改parent node的constraint，这里对_constraint赋值而不是constraint
					}
				}
			}
		}
		MemoryUsage();
	}

	// Fine-to-coarse down-sampling of constraints，maxDepth-1，比实际Octree的深度又少一层，因为是在parent node间进行_constraint下采样
	//这里采用了指针偏移的方式，那么_constraints + _sNodes.nodeCount[d]应该代表第d层_constraints的开始，下同理
	for( int d=maxDepth-1 ; d>=(_boundaryType==0?2:0) ; d-- ) DownSample( d , _sNodes , ( ConstPointer( Real ) )_constraints + _sNodes.nodeCount[d] , _constraints+_sNodes.nodeCount[d-1] );

	// Add the accumulated constraints from all finer depths
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<_sNodes.nodeCount[maxDepth] ; i++ ) constraints[i] += _constraints[i];
	//线程累加，maxDepth=_sNodes.maxDepth - 1，少了_sNodes最后一层的node，把刚才DownSample的结果更新到全局变量constraint中
	FreePointer( _constraints );


	std::vector< Point3D< Real > > coefficients( _sNodes.nodeCount[maxDepth] , zeroPoint );//_sNodes.nodeCount[maxDepth]代表所有父节点的node数目
	for( int d=maxDepth-1 ; d>=0 ; d-- )//maxDepth等于_sNodes.MaxDepth -1
	{
#pragma omp parallel for num_threads( threads )
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			int idx = normalInfo.index( _sNodes.treeNodes[i] );
			if( idx<0 ) continue;
			coefficients[i] = normalInfo.data[ idx ];//coefficients赋值为normalInfo，应该是加入Laplacian之后去赋值系数矩阵
		}
	}

	// Coarse-to-fine up-sampling of coefficients，用child和parent之间的位置关系来进行加权，更新finer resolution coefficients
	for( int d=(_boundaryType==0?2:0) ; d<maxDepth ; d++ ) UpSample( d , _sNodes , ( ConstPointer( Point3D< Real > ) ) GetPointer( coefficients ) + _sNodes.nodeCount[d-1] , GetPointer( coefficients ) + _sNodes.nodeCount[d] );//这里没有问题???d=0时d-1小于零；如果认为d代表finer depth，那么d应该从1或者minDepth+1开始吧

	// Compute the contribution from all coarser depths，为什么最后又重新更新了coarse depth对finer的constraint
	for( int d=0 ; d<=maxDepth ; d++ )//maxDepth = _sNodes.maxDepth - 1;
	{
		size_t start = _sNodes.nodeCount[d] , end = _sNodes.nodeCount[d+1] , range = end - start;
		Stencil< Point3D< double > , 5 > stencils[2][2][2];
		SetDivergenceStencils( d , integrator , stencils , false );//false表明以Neighbor node的normal field进行计算
		std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( maxDepth );//这个不能放在循环外吗???
#pragma omp parallel for num_threads( threads )
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
			TreeOctNode* node = _sNodes.treeNodes[i];
			int depth = node->depth();
			if( !depth ) continue;
			int startX=0 , endX=5 , startY=0 , endY=5 , startZ=0 , endZ=5;
			UpdateCoarserSupportBounds( node , startX , endX , startY  , endY , startZ , endZ );//修改StartX等
			typename TreeOctNode::Neighbors5 neighbors5;
			neighborKey.getNeighbors( node->parent , neighbors5 );

			bool isInterior;
			{
				int d , off[3];
				node->depthAndOffset( d , off );
				int o = _boundaryType==0 ? (1<<(d-2)) : 0;
				int mn = 4+o , mx = (1<<d)-4-o;
				//在这个范围内可以保证左右两侧的overlapping parent BSpline都存在
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
			{
				if( neighbors5.neighbors[x][y][z] )//parent node neighbor
				{
					TreeOctNode* _node = neighbors5.neighbors[x][y][z];
					int _i = _node->nodeData.nodeIndex;//parent node index
					if( isInterior )//isInterior为真，直接使用模板数据
					{
						Point3D< double >& div = _stencil.values[x][y][z];
						Point3D< Real >& normal = coefficients[_i];//这里用的是upsample的normal field
						constraint += Real( div[0] * normal[0] + div[1] * normal[1] + div[2] * normal[2] );
					}
					else
					{
						int _d , _off[3];
						_node->depthAndOffset( _d , _off );
						constraint += Real( GetDivergence2( integrator , d , off , _off , true , coefficients[_i] ) );
					}
				}
			}
			constraints[ node->nodeData.nodeIndex ] += constraint;
		}
	}
	MemoryUsage();
	return constraints;
}
template< class Real >
void Octree< Real >::refineBoundary( std::vector< int >* map ){ _sNodes.set( tree , map ); }


template< class Real >
template< class V >
V Octree< Real >::getCenterValue( const typename TreeOctNode::ConstNeighborKey3& neighborKey , const TreeOctNode* node , ConstPointer( V ) solution , ConstPointer( V ) metSolution , const typename BSplineData< 2 >::template CenterEvaluator< 1 >& evaluator , const Stencil< double , 3 >& stencil , const Stencil< double , 3 >& pStencil , bool isInterior ) const
{
	//stencil是相同Depth下的stencil，pStencil是child与parent的stencil，isInterior表示边界flag
	//metSolution记录了自上而下的upSample过的求解结果，solution记录原始求解结果
	if( node->children ) fprintf( stderr , "[WARNING] getCenterValue assumes leaf node\n" );//只允许leaf node，也就是最后一层的node
	V value(0);

	int d , off[3];
	node->depthAndOffset( d , off );

	if( isInterior )
	{
		//isInterior代表非边界node，处理的是正常情况，两侧的neighbor node都存在
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
		{
			const TreeOctNode* n = neighborKey.neighbors[d].neighbors[i][j][k];
			if( n ) value += solution[ n->nodeData.nodeIndex ] * Real( stencil.values[i][j][k] );//同一层之间node center value的weighted average
		}
		if( d>_minDepth )
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
			{
				const TreeOctNode* n = neighborKey.neighbors[d-1].neighbors[i][j][k];
				if( n ) value += metSolution[n->nodeData.nodeIndex] * Real( pStencil.values[i][j][k] );//从parent处upsample过的solution的weighted average
			}
	}
	else
	{
		//非标准型
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
		{
			const TreeOctNode* n = neighborKey.neighbors[d].neighbors[i][j][k];
			if( n )
			{
				int _d , _off[3];
				n->depthAndOffset( _d , _off );
				//无标准模式，临时单独求解
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
template< class Real >
template< class V >
V Octree< Real >::getCornerValue( const typename TreeOctNode::ConstNeighborKey3& neighborKey , const TreeOctNode* node , int corner , ConstPointer( V ) solution , ConstPointer( V ) metSolution , const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , const Stencil< double , 3 >& stencil , const Stencil< double , 3 > stencils[8] , bool isInterior ) const
{
	V value(0);
	int d , off[3];
	node->depthAndOffset( d , off );

	int cx , cy , cz;
	int startX = 0 , endX = 3 , startY = 0 , endY = 3 , startZ = 0 , endZ = 3;
	Cube::FactorCornerIndex( corner , cx , cy , cz );
	{
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d];
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
				if( _node ) value += solution[ _node->nodeData.nodeIndex ] * Real( stencil.values[x][y][z] );
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					value += solution[ _node->nodeData.nodeIndex ] * Real( evaluator.value( d , off[0] , cx , _off[0] , false , false ) * evaluator.value( d , off[1] , cy , _off[1] , false , false ) * evaluator.value( d , off[2] , cz , _off[2] , false , false ) );
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
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d-1];
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( _node ) value += metSolution[ _node->nodeData.nodeIndex ] * Real( stencils[_corner].values[x][y][z] );
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					value += metSolution[ _node->nodeData.nodeIndex ] * Real( evaluator.value( d , off[0] , cx , _off[0] , false , true ) * evaluator.value( d , off[1] , cy , _off[1] , false , true ) * evaluator.value( d , off[2] , cz , _off[2] , false , true ) );
				}
			}
	}
	return Real( value );
}
template< class Real >
std::pair< Real , Point3D< Real > > Octree< Real >::getCornerValueAndNormal( const typename TreeOctNode::ConstNeighborKey3& neighborKey , const TreeOctNode* node , int corner , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , const Stencil< double , 3 >& vStencil , const Stencil< double , 3 > vStencils[8] , const Stencil< Point3D< double > , 3 >& nStencil , const Stencil< Point3D< double > , 3 > nStencils[8] , bool isInterior ) const
{
	double value = 0;
	Point3D< double > normal;
	int d , off[3];
	node->depthAndOffset( d , off );

	int cx , cy , cz;
	int startX = 0 , endX = 3 , startY = 0 , endY = 3 , startZ = 0 , endZ = 3;
	Cube::FactorCornerIndex( corner , cx , cy , cz );
	{
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d];
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
				if( _node ) value += solution[ _node->nodeData.nodeIndex ] * vStencil.values[x][y][z] , normal += nStencil.values[x][y][z] * solution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					double v [] = { evaluator.value( d , off[0] , cx , _off[0] , false , false ) , evaluator.value( d , off[1] , cy , _off[1] , false , false ) , evaluator.value( d , off[2] , cz , _off[2] , false , false ) };
					double dv[] = { evaluator.value( d , off[0] , cx , _off[0] , true  , false ) , evaluator.value( d , off[1] , cy , _off[1] , true  , false ) , evaluator.value( d , off[2] , cz , _off[2] , true  , false ) };
					value += solution[ _node->nodeData.nodeIndex ] * evaluator.value( d , off[0] , cx , _off[0] , false , false ) * evaluator.value( d , off[1] , cy , _off[1] , false , false ) * evaluator.value( d , off[2] , cz , _off[2] , false , false );
					normal += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * solution[ _node->nodeData.nodeIndex ];
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
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d-1];
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( _node ) value += metSolution[ _node->nodeData.nodeIndex ] * vStencils[_corner].values[x][y][z] , normal += nStencils[_corner].values[x][y][z] * metSolution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					double v [] = { evaluator.value( d , off[0] , cx , _off[0] , false , true ) , evaluator.value( d , off[1] , cy , _off[1] , false , true ) , evaluator.value( d , off[2] , cz , _off[2] , false , true ) };
					double dv[] = { evaluator.value( d , off[0] , cx , _off[0] , true  , true ) , evaluator.value( d , off[1] , cy , _off[1] , true  , true ) , evaluator.value( d , off[2] , cz , _off[2] , true  , true ) };
					value += metSolution[ _node->nodeData.nodeIndex ] * evaluator.value( d , off[0] , cx , _off[0] , false , true ) * evaluator.value( d , off[1] , cy , _off[1] , false , true ) * evaluator.value( d , off[2] , cz , _off[2] , false , true );
					normal += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * metSolution[ _node->nodeData.nodeIndex ];
				}
			}
	}
	return std::pair< Real , Point3D< Real > >( Real( value ) , Point3D< Real >( normal ) );
}
template< class Real >
Point3D< Real > Octree< Real >::getCornerNormal( const typename TreeOctNode::ConstNeighbors5& neighbors5 , const typename TreeOctNode::ConstNeighbors5& pNeighbors5 , const TreeOctNode* node , int corner , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , const Stencil< Point3D< double > , 5 >& nStencil , const Stencil< Point3D< double > , 5 > nStencils[8] , bool isInterior ) const
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

template< class Real >
template< class V >
V Octree< Real >::_Evaluate( ConstPointer( V ) coefficients , Point3D< Real > p , typename TreeOctNode::ConstNeighborKey3& neighborKey3 ) const
{
	V value = V(0);

	for( int d=0 ; d<=neighborKey3.depth() ; d++ ) for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
	{
		const TreeOctNode* n = neighborKey3.neighbors[d].neighbors[i][j][k];
		if( n )
		{
			int dd , off[3];
			n->depthAndOffset( dd , off );
			value += 
				(
				coefficients[ n->nodeData.nodeIndex ] *
				(Real)
				(
					_fData.baseFunctions[ BinaryNode::CenterIndex( dd , off[0] ) ]( p[0] ) *
					_fData.baseFunctions[ BinaryNode::CenterIndex( dd , off[1] ) ]( p[1] ) *
					_fData.baseFunctions[ BinaryNode::CenterIndex( dd , off[2] ) ]( p[2] )
				)
			);
		}
	}
	return value;
}
template< class Real >
template< class V >
V Octree< Real >::_Evaluate( const SparseNodeData< V >& coefficients , Point3D< Real > p , typename TreeOctNode::ConstNeighborKey3& neighborKey3 ) const
{
	V value = V(0);

	double functionValues[3][3];
	for( int d=0 ; d<=neighborKey3.depth() ; d++ )
	{
		{
			const TreeOctNode* n = neighborKey3.neighbors[d].neighbors[1][1][1];
			if( !n ) fprintf( stderr , "[ERROR] _Evaluate: Center node needs to exist\n" ) , exit( 0 );
			int dd , off[3];
			n->depthAndOffset( dd , off );
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) if( off[i]-1+j>=0 && off[i]-1+j<(1<<dd) )
				functionValues[i][j] = _fData.baseFunctions[ BinaryNode::CenterIndex( dd , off[i]-1+j ) ]( p[i] );
		}

		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
		{
			double functionValue = functionValues[0][i] * functionValues[1][j];
			for( int k=0 ; k<3 ; k++ )
			{
				const TreeOctNode* n = neighborKey3.neighbors[d].neighbors[i][j][k];
				if( n )
				{
					int idx = coefficients.index( n );
					if( idx>=0 ) value += coefficients.data[idx] * (Real)( functionValue * functionValues[2][k] );
				}
			}
		}
	}
	return value;
}

template< class Real >
template< class V >
V Octree< Real >::Evaluate( ConstPointer( V ) coefficients , Point3D< Real > p , const BSplineData< 2 >* fData , int depth ) const
{
	int _depth;
	if( depth<0 ) _depth = -1;
	else _depth = _boundaryType==0 ? depth+1 : depth;

	V value = V(0);
	BSplineData< 2 > _fData;
	if( !fData ) _fData.set( tree.maxDepth() , _boundaryType ) , fData = &_fData;
	const TreeOctNode* n = tree.nextNode();
	while( n )
	{
		Point3D< Real > c;
		Real w;
		n->centerAndWidth( c , w );
		c -= p , w *= Real(1.5);
		if( fabs(c[0])>w || fabs(c[1])>w || fabs(c[2])>w || ( _depth>=0 && n->depth()>_depth ) )
		{
			n = tree.nextBranch( n );
			continue;
		}
		int d , off[3];
		n->depthAndOffset( d , off );
		value += 
			(
			coefficients[ n->nodeData.nodeIndex ] *
			(Real)
			(
				fData->baseFunctions[ BinaryNode::CenterIndex( d , off[0] ) ]( p[0] ) *
				fData->baseFunctions[ BinaryNode::CenterIndex( d , off[1] ) ]( p[1] ) *
				fData->baseFunctions[ BinaryNode::CenterIndex( d , off[2] ) ]( p[2] )
				)
			);
		n = tree.nextNode( n );
	}
	return value;
}
template< class Real >
template< class V >
V Octree< Real >::Evaluate( const SparseNodeData< V >& coefficients , Point3D< Real > p , const BSplineData< 2 >* fData , int depth ) const
{
	int _depth;
	if( depth<0 ) _depth = -1;
	else _depth = _boundaryType==0 ? depth+1 : depth;
	V value = V(0);
	BSplineData< 2 > _fData;
	if( !fData ) _fData.set( tree.maxDepth() , _boundaryType ) , fData = &_fData;
	const TreeOctNode* n = tree.nextNode();
	while( n )
	{
		Point3D< Real > c;
		Real w;
		n->centerAndWidth( c , w );
		c -= p , w *= Real(1.5);
		if( fabs(c[0])>w || fabs(c[1])>w || fabs(c[2])>w || (_depth>=0 && n->depth()>_depth ) ) n = tree.nextBranch( n );
		else
		{
			int d , off[3];
			n->depthAndOffset( d , off );

			int idx = coefficients.index( n );
			if( idx>=0 )
				value += 
				(
				coefficients.data[idx] *
				(Real)
				(
					fData->baseFunctions[ BinaryNode::CenterIndex( d , off[0] ) ]( p[0] ) *
					fData->baseFunctions[ BinaryNode::CenterIndex( d , off[1] ) ]( p[1] ) *
					fData->baseFunctions[ BinaryNode::CenterIndex( d , off[2] ) ]( p[2] )
					)
				);
			n = tree.nextNode( n );
		}
	}
	return value;
}
template< class Real >
template< class V >
Pointer( V ) Octree< Real >::Evaluate( ConstPointer( V ) coefficients , int& res , Real isoValue , int depth )
{
	int maxDepth = _boundaryType==0 ? tree.maxDepth()-1 : tree.maxDepth();
	if( depth<=0 || depth>maxDepth ) depth = maxDepth;
	res = 1<<depth;
	typename BSplineData< 2 >::template ValueTables< Real > vTables = _fData.template getValueTables< Real >( _fData.VALUE_FLAG );
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
			idx[i] = BinaryNode::CenterIndex( d , idx[i] );
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
		V coefficient = coefficients[ n->nodeData.nodeIndex ];
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
	for( int i=0 ; i<res*res*res ; i++ ) values[i] -= isoValue;

	return values;
}

////////////////
// VertexData //
////////////////
long long VertexData::CenterIndex(const TreeOctNode* node,int maxDepth)
{
	//指定node的中心点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	int idx[DIMENSION];
	return CenterIndex(node,maxDepth,idx);
}
long long VertexData::CenterIndex(const TreeOctNode* node,int maxDepth,int idx[DIMENSION])
{
	//指定node的中心点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++) idx[i]=BinaryNode::CornerIndex( maxDepth+1 , d+1 , o[i]<<1 , 1 );
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::CenterIndex( int depth , const int offSet[DIMENSION] , int maxDepth , int idx[DIMENSION] )
{
	//指定node(知道depth和offset)的中心点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	for(int i=0;i<DIMENSION;i++) idx[i]=BinaryNode::CornerIndex( maxDepth+1 , depth+1 , offSet[i]<<1 , 1 );
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::CornerIndex(const TreeOctNode* node,int cIndex,int maxDepth)
{
	//指定node的指定corner(cIndex)点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	int idx[DIMENSION];
	return CornerIndex(node,cIndex,maxDepth,idx);
}
long long VertexData::CornerIndex( const TreeOctNode* node , int cIndex , int maxDepth , int idx[DIMENSION] )
{
	//指定node的指定corner(cIndex)点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	int x[DIMENSION];
	Cube::FactorCornerIndex( cIndex , x[0] , x[1] , x[2] );
	int d , o[3];
	node->depthAndOffset( d , o );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , d , o[i] , x[i] );
	return CornerIndexKey( idx );
}
long long VertexData::CornerIndex( int depth , const int offSet[DIMENSION] , int cIndex , int maxDepth , int idx[DIMENSION] )
{
	//指定node(知道depth和offset)的指定corner(cIndex)点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	int x[DIMENSION];
	Cube::FactorCornerIndex( cIndex , x[0] , x[1] , x[2] );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , depth , offSet[i] , x[i] );
	return CornerIndexKey( idx );
}
long long VertexData::CornerIndexKey( const int idx[DIMENSION] )
{
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::FaceIndex(const TreeOctNode* node,int fIndex,int maxDepth){
	//指定node，指定face(fIndex)的中心点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	int idx[DIMENSION];
	return FaceIndex(node,fIndex,maxDepth,idx);
}
long long VertexData::FaceIndex(const TreeOctNode* node,int fIndex,int maxDepth,int idx[DIMENSION])
{
	//指定node，指定face(fIndex)的中心点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	int dir,offset;
	Cube::FactorFaceIndex(fIndex,dir,offset);
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);}
	idx[dir]=BinaryNode::CornerIndex(maxDepth+1,d,o[dir],offset);
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth )
{ 
	//指定node，指定edge(eIndex)的中心点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	int idx[DIMENSION] ; 
	return EdgeIndex( node , eIndex , maxDepth , idx ); 
}
long long VertexData::EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth , int idx[DIMENSION] )
{
	//指定node，指定edge(eIndex)的中心点，在深度达到maxDepth的全划分(center和corner)点序列中的index
	int o , i1 , i2;
	int d , off[3];
	node->depthAndOffset( d ,off );//得到depth和offset的具体值
	Cube::FactorEdgeIndex( eIndex , o , i1 , i2 );//根据edgeindex得到orientation等信息
	//难道是把每一维都看作是一个binary node
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , d+1 , off[i]<<1 , 1 );
	switch(o)
	{
		case 0:
			idx[1] = BinaryNode::CornerIndex( maxDepth+1 , d , off[1] , i1 );
			idx[2] = BinaryNode::CornerIndex( maxDepth+1 , d , off[2] , i2 );
			break;
		case 1:
			idx[0] = BinaryNode::CornerIndex( maxDepth+1 , d , off[0] , i1 );
			idx[2] = BinaryNode::CornerIndex( maxDepth+1 , d , off[2] , i2 );
			break;
		case 2:
			idx[0] = BinaryNode::CornerIndex( maxDepth+1 , d , off[0] , i1 );
			idx[1] = BinaryNode::CornerIndex( maxDepth+1 , d , off[1] , i2 );
			break;
	};
	//最后还是通过移位操作来存储edge index
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
