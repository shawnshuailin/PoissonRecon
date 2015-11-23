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

#ifndef POLYNOMIAL_INCLUDED
#define POLYNOMIAL_INCLUDED

#include <vector>

template< int Degree >
class Polynomial{
public:
	double coefficients[Degree+1];//Degree代表了多项式的幂次，因此一共有Degree+1个

	Polynomial(void);
	template<int Degree2>
	Polynomial(const Polynomial<Degree2>& P);
	double operator()( double t ) const;//给定未知数的值，算出多项式的值
	double integral( double tMin , double tMax ) const;//在指定上下界范围内对多项式进行积分

	int operator == (const Polynomial& p) const;//判断多项式相等或不等
	int operator != (const Polynomial& p) const;
	int isZero(void) const;//判零或置零
	void setZero(void);

	template<int Degree2>
	Polynomial& operator  = (const Polynomial<Degree2> &p);
	Polynomial& operator += (const Polynomial& p);
	Polynomial& operator -= (const Polynomial& p);
	Polynomial  operator -  (void) const;
	Polynomial  operator +  (const Polynomial& p) const;
	Polynomial  operator -  (const Polynomial& p) const;
	template<int Degree2>
	Polynomial<Degree+Degree2>  operator *  (const Polynomial<Degree2>& p) const;

	Polynomial& operator += ( double s );
	Polynomial& operator -= ( double s );
	Polynomial& operator *= ( double s );
	Polynomial& operator /= ( double s );
	Polynomial  operator +  ( double s ) const;
	Polynomial  operator -  ( double s ) const;
	Polynomial  operator *  ( double s ) const;
	Polynomial  operator /  ( double s ) const;

	Polynomial scale( double s ) const;//声明一个新多项式，多项式的参数除以scale
	Polynomial shift( double t ) const;//所谓的shift就是把原变量移动了-t，然后所有的系数会发生变化，推导一下就可以看出规律

	Polynomial<Degree-1> derivative(void) const;//多项式的导数
	Polynomial<Degree+1> integral(void) const;//多项式积分

	void printnl(void) const;

	Polynomial& addScaled(const Polynomial& p,double scale);

	static void Negate(const Polynomial& in,Polynomial& out);
	static void Subtract(const Polynomial& p1,const Polynomial& p2,Polynomial& q);
	static void Scale(const Polynomial& p,double w,Polynomial& q);
	static void AddScaled(const Polynomial& p1,double w1,const Polynomial& p2,double w2,Polynomial& q);
	static void AddScaled(const Polynomial& p1,const Polynomial& p2,double w2,Polynomial& q);
	static void AddScaled(const Polynomial& p1,double w1,const Polynomial& p2,Polynomial& q);

	void getSolutions(double c,std::vector<double>& roots,double EPS) const;//多项式求解，c是一个常量，roots用来保存结果，EPS表示极小值，用来判断系数是否接近零
	int getSolutions( double c , double* roots , double EPS ) const;//同上

	static Polynomial BSplineComponent( int i );//计算幂次为Degree的B样条曲线的每一个多项式
};

#include "Polynomial.inl"
#endif // POLYNOMIAL_INCLUDED
