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
struct BSplineElementCoefficients//��ʾϵ����struct
{
	int coeffs[Degree+1];//degree+1���ÿһ���ϵ��������Щϵ�����ڲ�ͬ��BSpline��ֻ����ͬһ��interval��
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
//Ҳ����˵BSplineElements����Ҫ�ж��interval��������Щ��������ÿ��interval�Ϻ������ӵ�weight�ʹ洢��һ��BSplineElementCoefficients��
template< int Degree >
struct BSplineElements : public std::vector< BSplineElementCoefficients< Degree > >
{
	static const int _off = (Degree+1)/2;//_off���ݴ�ΪDegreeʱ����i��interval�϶����BSpline������า���������i-_off�����Ҳ���i-_off+Degree
	void _addLeft ( int offset , int boundary );//addLeft��addRight����ѭ���������õ�
	void _addRight( int offset , int boundary );
public:
	enum
	{
		NONE      =  0,
		DIRICHLET = -1,
		NEUMANN   =  1
	};
	// Coefficients are ordered as "/" "-" "\"
	int denominator;//��ĸ��ʲô�ķ�ĸ???

	BSplineElements( void ) { denominator = 1; }
	BSplineElements( int res , int offset , int boundary=NONE , int inset=0 );

	void upSample( BSplineElements& high ) const;
	void differentiate( BSplineElements< Degree-1 >& d ) const;

	void print( FILE* fp=stdout ) const
	{
		for( int i=0 ; i<std::vector< BSplineElementCoefficients< Degree > >::size() ; i++ )//�������size����������???��Ϊʲô����this->size()
		{
			printf( "%d]" , i );
			for( int j=0 ; j<=Degree ; j++ ) printf( " %d" , (*this)[i][j] );//��Ȼ����this����operator []��Ϊʲô���ܵ���size����
			printf( " (%d)\n" , denominator );
		}
	}
};

//BSplineData�������������BSpline��dot��dDot�Լ�d2Dot�İ�
//BSplineData�����൱��version1��FunctionData���÷�Ҳ����
template< int Degree >
class BSplineData
{
	int _boundaryType;
	//��һ��Degree+1�����һ��B���������ݴ�ΪDegree���ڶ�������ڶ���B���������ݴ�ΪDegree��������������B����dot�Ľ��
	double _vvIntegrals[Degree+1][Degree+1];
	double _vdIntegrals[Degree+1][Degree  ];//ÿһ��洢���Ƕ���ʽ�˻������ֽ������������_vdIntegrals�Ȳų����������ٵ����
	double _dvIntegrals[Degree  ][Degree+1];
	double _ddIntegrals[Degree  ][Degree  ];

public:
	struct Integrator
	{
		struct IntegralTables
		{
			//����vv_ccIntegrals���ں��汻���õ���������Եó���offset��������[0, Degree)��[res-Degree, res-1]֮��ʱ��BSpline��Ҫ���⴦��
			//���vector filed��2*Degree�����������ʣ��һ���������ģ��ܹ�2*Degree+1���������offsetָ����BSpline����node���ڵ�λ��
			//ÿ��BSpline���Ը���Degree+1��interval������Щinterval�Ͽ�����֮���������ص�����2Degree+1��BSpline�����Ǻ�һ��2*Degree+1��ԭ��
			double vv_ccIntegrals[2*Degree+1][2*Degree+1] , vv_cpIntegrals[(2*Degree+1)*2][2*Degree+1];//p����parent node
			double dv_ccIntegrals[2*Degree+1][2*Degree+1] , dv_cpIntegrals[(2*Degree+1)*2][2*Degree+1];
			double vd_ccIntegrals[2*Degree+1][2*Degree+1] , vd_cpIntegrals[(2*Degree+1)*2][2*Degree+1];
			double dd_ccIntegrals[2*Degree+1][2*Degree+1] , dd_cpIntegrals[(2*Degree+1)*2][2*Degree+1];//d����derivative
		};
		std::vector< IntegralTables > iTables;
		double dot( int depth , int off1 , int off2 , bool d1 , bool d2 , bool childParent=false ) const;
	};
	double dot( int depth1 , int off1 , int depth2 , int off2 , bool d1 , bool d2 , bool inset=false ) const;
	void setIntegrator( Integrator& integrator , bool inset , bool useDotRatios=false ) const;
	template< int Radius >
	struct CenterEvaluator//centerEvaluatorӦ���Ǽ������ĵ��scalar function��ֵ
	{
		struct ValueTables
		{
			double vValues[2*Degree+1][ 3*(2*Radius+1) ];//���������radius������Ϊ[-radius..., 0, ..., radius]��
			//������3����Ϊindex����1�Ǹ���ǰ��neighbor�����ݣ�ʣ�µ������Ǹ�����child neighbor������
			double dValues[2*Degree+1][ 3*(2*Radius+1) ];//Degree��2*Degree+1����Ϊ������ı�Ե���������⴦���
		};
		std::vector< ValueTables > vTables;
		double value( int depth , int off1 , int off2 , bool d , bool childParent=false ) const;
	};
	template< int Radius >
	void setCenterEvaluator( CenterEvaluator< Radius >& evaluator , double smoothingRadius , double dSmoothingRadius, bool inset ) const;
	double value( int depth , int off , double smoothingRadius , double s , bool d , bool inset=false ) const;
	template< int Radius >
	struct CornerEvaluator//corner evaluators����������corner���scalar function��ֵ
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
	//That��s what stored in the ��_polys�� member of the component.
	struct BSplineComponents//��BSplineϵ��ΪDegreeʱ��һ����Degree+1�����ߣ���Ҫ��Degree+1������ʽ����ʾ��Щ���ߣ�����ÿ������ʽ��ϵ����Degree+1��
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
	const static int  DV_DOT_FLAG = 2;//derivative * value��set divergenceʱ���õ�
	const static int  DD_DOT_FLAG = 4;//derivative * derivative��set laplacianʱ���õ�
	const static int   VALUE_FLAG = 1;
	const static int D_VALUE_FLAG = 2;
	template< class Real >
	struct DotTables//dottable������Ax=b�ľ����л��õ�����Ϊ�漰set laplacian�Լ�set divergence��������һ�汾�к���û�õ���getDotTables����Ҳû�б�����
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
	struct ValueTables//valueTable�ǻ��棬�����洢���������function value��derivative function value�����������corner value����center valueʹ��
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