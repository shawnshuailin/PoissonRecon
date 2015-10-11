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
	Polynomial<Degree> p;//����ݴ�ΪDegree�Ķ���ʽ
	double start;//start������������ʾ����ʽ��δ֪����ֵ�ռ����Сֵ��������scale��shift�����ж��ܵ���Ӱ��
	//��Щstartֻ��ʾ��ÿ������ʽ����Сȡֵ����������֮���������������ģ�����һ��start��������һ��polynomial�Ľ�����
	//����polynomial�Ķ�������[start, infinity]
	template<int Degree2>
	StartingPolynomial<Degree+Degree2>  operator * (const StartingPolynomial<Degree2>& p) const;//����ʽ���
	StartingPolynomial scale(double s) const;//��������
	StartingPolynomial shift(double t) const;//����ƽ��
	int operator < (const StartingPolynomial& sp) const;//ɶ���⣬�Ƚϴ�Сֻ�Ƚ�start��û̫�����÷�
	static int Compare(const void* v1,const void* v2);//С�ڷ��ź������compare����������container�ж�StartingPolynomial�����������ʱ��Ҫ��������Ƚ�
};

template<int Degree>
class PPolynomial//���PPolynomial��StaringPolynomial
{
//��ν��PPolynomial��ָ�ڱ�ʾB��������ʱ�����ڲ�ͬ�����������ʽ��ʾ������ͬ�����Ҫ��һ������ʽvector���洢�������ʽ
public:
	size_t polyCount;//��ǰ�������ж��ٸ�StartingPolynomial
	//��Щstartֻ��ʾ��ÿ������ʽ����Сȡֵ
	//polys��һ��array���洢�˶������ʽ����ֵ�ռ��������Ը�������ʽ��startΪ�ָ���������
	//����ȷʵ֤����PPolynomial�в�ͬ����ʽ��start��ǰ���ϵ��������һ��start��������һ������ʽ��ֵ��Χ�����ޣ����еĶ���ʽ���޶���infinity
	StartingPolynomial<Degree>* polys;

	PPolynomial(void);
	PPolynomial(const PPolynomial<Degree>& p);
	~PPolynomial(void);

	PPolynomial& operator = (const PPolynomial& p);

	int size(void) const;//����ʵ��ռ�õ��ڴ�ռ�

	void set( size_t size );
	// Note: this method will sort the elements in sps
	void set( StartingPolynomial<Degree>* sps , int count );
	void reset( size_t newSize );


	double operator()( double t ) const;
	double integral( double tMin , double tMax ) const;
	double Integral( void ) const;

	//������ʾ�˶������ʽ�Ķ���֮��ļӼ��������Լ���ֵ������Ҫ����start��С�����ѵ������ж�û��start��ȵ�������֣���һ����֤start������
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
