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

#ifndef P_POLYNOMIAL_INCLUDED
#define P_POLYNOMIAL_INCLUDED
#include <vector>
#include "Polynomial.h"

template<int Degree>
class StartingPolynomial{
public:
	Polynomial<Degree> p;//最高幂次为Degree的多项式
	double start;//start好像是用来表示多项式的未知数赋值空间的最小值，所以在scale和shift操作中都受到了影响
	//这些start只表示了每个多项式的最小取值，而各区间之间好像是连续定义的，但后一个start并不是上一个polynomial的结束，
	//所有polynomial的定义域都是[start, infinity]
	template<int Degree2>
	StartingPolynomial<Degree+Degree2>  operator * (const StartingPolynomial<Degree2>& p) const;//多项式相乘
	StartingPolynomial scale(double s) const;//参数缩放
	StartingPolynomial shift(double t) const;//参数平移
	int operator < (const StartingPolynomial& sp) const;//啥玩意，比较大小只比较start，没太看懂用法
	static int Compare(const void* v1,const void* v2);//小于符号和这里的compare函数都是在container中对StartingPolynomial对象进行排序时需要的弱排序比较
};

template<int Degree>
class PPolynomial//查查PPolynomial和StaringPolynomial
{
//所谓的PPolynomial是指在表示B样条曲线时，对于不同定义域，其多项式表示并不相同，因此要用一个多项式vector来存储多个多项式
public:
	size_t polyCount;//当前对象中有多少个StartingPolynomial
	//这些start只表示了每个多项式的最小取值
	//polys是一个array，存储了多个多项式，赋值空间连续，以各个多项式的start为分隔划分区间
	//这里确实证明了PPolynomial中不同多项式的start的前后关系，但是下一个start并不是上一个多项式阈值范围的上限，所有的多项式上限都是infinity
	StartingPolynomial<Degree>* polys;

	PPolynomial(void);
	PPolynomial(const PPolynomial<Degree>& p);
	~PPolynomial(void);

	PPolynomial& operator = (const PPolynomial& p);

	int size(void) const;//返回实际占用的内存空间

	void set( size_t size );
	// Note: this method will sort the elements in sps
	void set( StartingPolynomial<Degree>* sps , int count );
	void reset( size_t newSize );


	double operator()( double t ) const;
	double integral( double tMin , double tMax ) const;
	double Integral( void ) const;

	//两个表示了多个多项式的对象之间的加减乘运算以及赋值操作，要按照start大小排序，难道两个中都没有start相等的情况出现，进一步验证start的意义
	template<int Degree2>
	PPolynomial<Degree>& operator = (const PPolynomial<Degree2>& p);

	PPolynomial  operator + (const PPolynomial& p) const;
	PPolynomial  operator - (const PPolynomial& p) const;

	template<int Degree2>
	PPolynomial<Degree+Degree2> operator * (const Polynomial<Degree2>& p) const;

	template<int Degree2>
	PPolynomial<Degree+Degree2> operator * (const PPolynomial<Degree2>& p) const;


	PPolynomial& operator += ( double s );
	PPolynomial& operator -= ( double s );
	PPolynomial& operator *= ( double s );
	PPolynomial& operator /= ( double s );
	PPolynomial  operator +  ( double s ) const;
	PPolynomial  operator -  ( double s ) const;
	PPolynomial  operator *  ( double s ) const;
	PPolynomial  operator /  ( double s ) const;

	PPolynomial& addScaled(const PPolynomial& poly,double scale);

	PPolynomial scale( double s ) const;
	PPolynomial shift( double t ) const;

	PPolynomial< Degree-1 > derivative(void) const;
	PPolynomial< Degree+1 > integral(void) const;

	void getSolutions(double c,std::vector<double>& roots,double EPS,double min=-DBL_MAX,double max=DBL_MAX) const;

	void printnl( void ) const;

	PPolynomial< Degree+1 > MovingAverage( double radius ) const;
	static PPolynomial BSpline( double radius=0.5 );

	void write( FILE* fp , int samples , double min , double max ) const;
};
#include "PPolynomial.inl"
#endif // P_POLYNOMIAL_INCLUDED
