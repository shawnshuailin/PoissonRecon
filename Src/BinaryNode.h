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


class BinaryNode//binarynodeò����������ʾ��һ��ά������½��л���ʱ��center��corner�������Ϊoctree��ÿ��ά���ϻ���ʱ����binary tree
{//Ȼ��Ż�����Ӧ��width��center position����Ϣ
public:
	static inline int CenterCount( int depth ) { return  1<<depth; }//���Ϊdepthʱ��center������depthӦ���Ǵ�0��ʼ��
	static inline int CornerCount( int depth ) { return (1<<depth)+1; }//���Ϊdepthʱ��corner����
	static inline int CumulativeCenterCount( int maxDepth ) { return (1<<(maxDepth+1))-1; }//�ۻ���center���������˶��壬��������Ǵ��㿪ʼ��ģ����㵽���ڵ�ľ��룬�߶��ǵ�Ҷ�ڵ�ľ���
	static inline int CumulativeCornerCount( int maxDepth ) { return (1<<(maxDepth+1))+maxDepth; }//�ۻ���corner����
	//����Ӧ����center����corner��cumulative����µ�index
	static inline int CenterIndex( int depth , int offSet ) { return (1<<depth)+offSet-1; }//�����Ϊdepth��ƫ��Ϊoffsetʱ��center index
	static inline int CornerIndex( int depth , int offSet ) { return (1<<depth)+offSet+depth; }//�����Ϊdepth��ƫ��Ϊoffsetʱ��corner index
	//��������Ǽ����ڸ���maxDepth�µĻ���ʱ�������binarynode���������ĵ�ͽǵ��г�����
	//��ôָ��depth��offset�ĵ���������������е�index�������÷����Լ���VertexData���еĵ��ú���
	static inline int CornerIndex( int maxDepth , int depth , int offSet , int forwardCorner ){ return (offSet+forwardCorner)<<(maxDepth-depth); }
	//��maxDepth�£���λ1���ȵ�node�����ֳ�1<<maxDepth�ݣ����ÿ��corner��position��Ӧ�ÿ��Եõ�
	template< class Real > static inline Real CornerIndexPosition(int index,int maxDepth){ return Real(index)/(1<<maxDepth); }
	template< class Real > static inline Real Width(int depth){ return Real(1.0/(1<<depth)); }//��׼�����µ�ǰ��ȵĵ�λnode���
	//������������ǰoffset��node�����ĵ������width
	template< class Real > static inline void CenterAndWidth( int depth , int offset , Real& center , Real& width )
	  {
	    width=Real (1.0/(1<<depth) );
	    center=Real((0.5+offset)*width);
	  }
	//�ָ���index�����depth��offset��Ȼ�����center��width
	template< class Real > static inline void CenterAndWidth( int idx , Real& center , Real& width )
	  {
	    int depth , offset;
	    DepthAndOffset( idx , depth , offset );
	    CenterAndWidth( depth , offset , center , width );
	  }
	//����index����depth��offset�������indexָ���ǵ�ǰ���ֲ���µ���һ��(����˵���node)���ۼ�index��
	//Ҳ���Ǵ�depth=0ʱ��ʼ������index
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
