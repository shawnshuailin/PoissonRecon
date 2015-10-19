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

#ifndef BINARY_NODE_INCLUDED
#define BINARY_NODE_INCLUDED

#define MSVC_2010_FIX 1


class BinaryNode//binarynode貌似是用来表示在一个维度情况下进行划分时的center和corner情况，因为octree在每个维度上划分时都是binary tree
{//然后才会有相应的width，center position等信息
public:
	static inline int CenterCount( int depth ) { return  1<<depth; }//深度为depth时的center数量，depth应该是从0开始的
	static inline int CornerCount( int depth ) { return (1<<depth)+1; }//深度为depth时的corner数量
	static inline int CumulativeCenterCount( int maxDepth ) { return (1<<(maxDepth+1))-1; }//累积的center数量，查了定义，树的深度是从零开始算的，计算到根节点的距离，高度是到叶节点的距离
	static inline int CumulativeCornerCount( int maxDepth ) { return (1<<(maxDepth+1))+maxDepth; }//累积的corner数量
	//下面应该是center或者corner在cumulative情况下的index
	static inline int CenterIndex( int depth , int offSet ) { return (1<<depth)+offSet-1; }//当深度为depth，偏移为offset时的center index
	static inline int CornerIndex( int depth , int offSet ) { return (1<<depth)+offSet+depth; }//当深度为depth，偏移为offset时的corner index
	//这个函数是计算在给定maxDepth下的划分时，如果把binarynode上所有中心点和角点列出来，
	//那么指定depth下offset的点在这个完整序列中的index，具体用法可以见其VertexData类中的调用函数
	static inline int CornerIndex( int maxDepth , int depth , int offSet , int forwardCorner ){ return (offSet+forwardCorner)<<(maxDepth-depth); }
	//在maxDepth下，单位1长度的node被划分成1<<maxDepth份，因此每个corner的position就应该可以得到
	template< class Real > static inline Real CornerIndexPosition(int index,int maxDepth){ return Real(index)/(1<<maxDepth); }
	template< class Real > static inline Real Width(int depth){ return Real(1.0/(1<<depth)); }//标准二分下当前深度的单位node宽度
	//根据深度求出当前offset下node的中心点坐标和width
	template< class Real > static inline void CenterAndWidth( int depth , int offset , Real& center , Real& width )
	  {
	    width=Real (1.0/(1<<depth) );
	    center=Real((0.5+offset)*width);
	  }
	//现根据index计算出depth和offset，然后计算center和width
	template< class Real > static inline void CenterAndWidth( int idx , Real& center , Real& width )
	  {
	    int depth , offset;
	    DepthAndOffset( idx , depth , offset );
	    CenterAndWidth( depth , offset , center , width );
	  }
	//根据index计算depth和offset，这里的index指的是当前划分层次下的这一段(或者说这个node)的累加index，
	//也就是从depth=0时开始计数的index
	static inline void DepthAndOffset( int idx , int& depth , int& offset )
	  {
	    int i=idx+1;
#if MSVC_2010_FIX
		depth = 0;
#else // !MSVC_2010_FIX
	    depth = -1;
#endif // MSVC_2010_FIX
	    while( i )
		{
	      i >>= 1;
	      depth++;
	    }
#if MSVC_2010_FIX
		depth--;
#endif // MSVC_2010_FIX
	    offset = ( idx+1 ) - (1<<depth);
	  }
};
#endif // BINARY_NODE_INCLUDED
