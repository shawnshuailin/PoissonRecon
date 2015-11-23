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

#ifndef BSPLINE_DATA_INCLUDED
#define BSPLINE_DATA_INCLUDED

#include "PPolynomial.h"
#include "Array.h"

//The BSplineElementCoefficients represents the linear combination of B-splines for a given interval
template< int Degree >
struct BSplineElementCoefficients//表示系数的struct
{
	int coeffs[Degree+1];//degree+1个项，每一项的系数，但这些系数属于不同的BSpline，只是在同一个interval上
	BSplineElementCoefficients( void ){ memset( coeffs , 0 , sizeof( int ) * ( Degree+1 ) ); }
	int& operator[]( int idx ){ return coeffs[idx]; }
	const int& operator[]( int idx ) const { return coeffs[idx]; }
};

//BSplineElements is a representation of the linear combination of B-splines. Consider a function expressed as a linear combination of B-splines of degree d. 
//On any interval, the restriction of the function can be expressed as the linear combination d+1 elements (polynomials) deriving from the d+1 different 
//B-splines that overlap the interval. The BSplineElementCoefficients represents the linear combination for a given interval and 
//the BSplineElements represents these linear combinations over all intervals.
//Note that with the BSplineElements, you can represent all linear combination of B-Splines as well as other functions. 
//(E.g. The linear combination of B-splines has to be smooth of degree d-1, while BSplineElements can represent discontinuous functions.) 
//Thus, by defining operations like differentiation, integration, and prolongation over BSplineElements, 
//we automatically get these for functions that are linear combinations of B-splines.
//也就是说BSplineElements首先要有多个interval来定义这些函数，而每个interval上函数叠加的weight就存储在一个BSplineElementCoefficients中
template< int Degree >
struct BSplineElements : public std::vector< BSplineElementCoefficients< Degree > >
{
	static const int _off = (Degree+1)/2;//_off是幂次为Degree时，第i个interval上定义的BSpline的最左侧覆盖区间就是i-_off，最右侧是i-_off+Degree
	void _addLeft ( int offset , int boundary );//addLeft和addRight都是循环迭代设置的
	void _addRight( int offset , int boundary );
public:
	enum
	{
		NONE      =  0,
		DIRICHLET = -1,
		NEUMANN   =  1
	};
	// Coefficients are ordered as "/" "-" "\"
	int denominator;//分母，什么的分母???

	BSplineElements( void ) { denominator = 1; }
	BSplineElements( int res , int offset , int boundary=NONE , int inset=0 );

	void upSample( BSplineElements& high ) const;
	void differentiate( BSplineElements< Degree-1 >& d ) const;

	void print( FILE* fp=stdout ) const
	{
		for( int i=0 ; i<std::vector< BSplineElementCoefficients< Degree > >::size() ; i++ )//这里调用size的意义在哪???，为什么不是this->size()
		{
			printf( "%d]" , i );
			for( int j=0 ; j<=Degree ; j++ ) printf( " %d" , (*this)[i][j] );//既然能用this调用operator []，为什么不能调用size函数
			printf( " (%d)\n" , denominator );
		}
	}
};

//BSplineData是用来处理各种BSpline的dot，dDot以及d2Dot的吧
//BSplineData作用相当于version1中FunctionData，用法也近似
template< int Degree >
class BSplineData
{
	int _boundaryType;
	//第一个Degree+1代表第一个B样条曲线幂次为Degree，第二个代表第二个B样条曲线幂次为Degree，因此这就是两个B样条dot的结果
	double _vvIntegrals[Degree+1][Degree+1];
	double _vdIntegrals[Degree+1][Degree  ];//每一项存储的是多项式乘积并积分结果，因此下面的_vdIntegrals等才出现项数减少的情况
	double _dvIntegrals[Degree  ][Degree+1];
	double _ddIntegrals[Degree  ][Degree  ];

public:
	struct Integrator
	{
		struct IntegralTables
		{
			//根据vv_ccIntegrals等在后面被调用的情况，可以得出在offset符合区间[0, Degree)和[res-Degree, res-1]之间时，BSpline需要特殊处理
			//因此vector filed有2*Degree个特殊情况，剩下一个是正常的，总共2*Degree+1个。这里的offset指的是BSpline中心node所在的位置
			//每个BSpline可以覆盖Degree+1个interval，在这些interval上可能与之发生区间重叠的有2Degree+1个BSpline，这是后一个2*Degree+1的原因
			double vv_ccIntegrals[2*Degree+1][2*Degree+1] , vv_cpIntegrals[(2*Degree+1)*2][2*Degree+1];//p代表parent node
			double dv_ccIntegrals[2*Degree+1][2*Degree+1] , dv_cpIntegrals[(2*Degree+1)*2][2*Degree+1];
			double vd_ccIntegrals[2*Degree+1][2*Degree+1] , vd_cpIntegrals[(2*Degree+1)*2][2*Degree+1];
			double dd_ccIntegrals[2*Degree+1][2*Degree+1] , dd_cpIntegrals[(2*Degree+1)*2][2*Degree+1];//d代表derivative
		};
		std::vector< IntegralTables > iTables;
		double dot( int depth , int off1 , int off2 , bool d1 , bool d2 , bool childParent=false ) const;
	};
	double dot( int depth1 , int off1 , int depth2 , int off2 , bool d1 , bool d2 , bool inset=false ) const;
	void setIntegrator( Integrator& integrator , bool inset , bool useDotRatios=false ) const;
	template< int Radius >
	struct CenterEvaluator//centerEvaluator应该是计算中心点的scalar function的值
	{
		struct ValueTables
		{
			double vValues[2*Degree+1][ 3*(2*Radius+1) ];//如果称作是radius，区间为[-radius..., 0, ..., radius]，
			//而乘以3是因为index等于1是跟当前层neighbor的数据，剩下的两个是跟两个child neighbor的数据
			double dValues[2*Degree+1][ 3*(2*Radius+1) ];//Degree是2*Degree+1是因为在两侧的边缘进行了特殊处理吧
		};
		std::vector< ValueTables > vTables;
		double value( int depth , int off1 , int off2 , bool d , bool childParent=false ) const;
	};
	template< int Radius >
	void setCenterEvaluator( CenterEvaluator< Radius >& evaluator , double smoothingRadius , double dSmoothingRadius, bool inset ) const;
	double value( int depth , int off , double smoothingRadius , double s , bool d , bool inset=false ) const;
	template< int Radius >
	struct CornerEvaluator//corner evaluators是用来计算corner点的scalar function的值
	{
		struct ValueTables
		{
			double vValues[2*Degree+1][4*Radius+3];
			double dValues[2*Degree+1][4*Radius+3];
		};
		std::vector< ValueTables > vTables;
		double value( int depth , int off1 , int c1 , int off2 , bool d , bool childParent=false ) const;
	};
	template< int Radius >
	void setCornerEvaluator( CornerEvaluator< Radius >& evaluator , double smoothingRadius , double dSmoothingRadius, bool inset ) const;

	//BSplineComponents is a representation of an individual B-spline. In general, if the B-spline is of degree d, it will be made up of d+1 piecewise polynomial pieces.
	//That’s what stored in the “_polys” member of the component.
	struct BSplineComponents//当BSpline系数为Degree时，一共有Degree+1段折线，需要有Degree+1个多项式来表示这些折线，而且每个多项式的系数是Degree+1个
	{
		Polynomial< Degree > polys[Degree+1];
		Polynomial< Degree >& operator[] ( int idx ) { return polys[idx]; }
		const Polynomial< Degree >& operator[] ( int idx ) const { return polys[idx]; }
		void printnl( void ) const  { for( int d=0 ; d<=Degree ; d++ ) polys[d].printnl(); }
		BSplineComponents scale( double s ) const { BSplineComponents b ; for( int d=0 ; d<=Degree ; d++ ) b[d] = polys[d].scale(s) ; return b; }
		BSplineComponents shift( double s ) const { BSplineComponents b ; for( int d=0 ; d<=Degree ; d++ ) b[d] = polys[d].shift(s) ; return b; }
	};

	int depth;
	size_t functionCount , sampleCount;
	PPolynomial< Degree   >  baseFunction ,  leftBaseFunction ,  rightBaseFunction ,  leftRightBaseFunction;
	PPolynomial< Degree-1 > dBaseFunction , dLeftBaseFunction , dRightBaseFunction , dLeftRightBaseFunction;
	BSplineComponents baseBSpline , leftBSpline , rightBSpline , leftRightBSpline;
	Pointer( PPolynomial< Degree > ) baseFunctions;
	Pointer( BSplineComponents ) baseBSplines;

	BSplineData( void );

	const static int  VV_DOT_FLAG = 1;
	const static int  DV_DOT_FLAG = 2;//derivative * value在set divergence时会用到
	const static int  DD_DOT_FLAG = 4;//derivative * derivative在set laplacian时会用到
	const static int   VALUE_FLAG = 1;
	const static int D_VALUE_FLAG = 2;
	template< class Real >
	struct DotTables//dottable在设置Ax=b的矩阵中会用到，因为涉及set laplacian以及set divergence，但在这一版本中好像没用到，getDotTables函数也没有被调用
	{
		size_t functionCount;
		Pointer( Real ) vvDotTable;
		Pointer( Real ) dvDotTable;
		Pointer( Real ) ddDotTable;

		DotTables( void );
		~DotTables( void );

		inline size_t Index( int i1 , int i2 ) const;
		static inline size_t SymmetricIndex( int i1 , int i2 );
		static inline int SymmetricIndex( int i1 , int i2 , size_t& index );
	};
	template< class Real >
	struct ValueTables//valueTable是缓存，用来存储所计算出的function value及derivative function value，供后面计算corner value或者center value使用
	{
		size_t functionCount , sampleCount;
		Pointer( Real ) valueTable;
		Pointer( Real ) dValueTable;

		ValueTables( void );
		~ValueTables( void );

		inline size_t Index( int i1 , int i2 ) const;
		void setSampleSpan( int idx , int& start , int& end , double smooth=0 ) const;
	};
	void set( int maxDepth , int boundaryType=BSplineElements< Degree >::NONE );
	template< class Real >
	typename BSplineData< Degree >::template DotTables< Real > getDotTables( int flags , bool useDotRatios=true , bool inset=false ) const;
	template< class Real >
	typename BSplineData< Degree >::template ValueTables< Real > getValueTables( int flags , double valueSmooth=0 , double normalSmooth=0 ) const;
};

template< int Degree1 , int Degree2 > void SetBSplineElementIntegrals( double integrals[Degree1+1][Degree2+1] );

#include "BSplineData.inl"
#endif // BSPLINE_DATA_INCLUDED