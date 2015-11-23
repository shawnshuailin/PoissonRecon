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

#include "Factor.h"

////////////////////////
// StartingPolynomial //
////////////////////////
template<int Degree>
template<int Degree2>
StartingPolynomial<Degree+Degree2> StartingPolynomial<Degree>::operator * (const StartingPolynomial<Degree2>& p) const{
	StartingPolynomial<Degree+Degree2> sp;
	if(start>p.start){sp.start=start;}
	else{sp.start=p.start;}
	sp.p=this->p*p.p;
	return sp;
}
template<int Degree>
StartingPolynomial<Degree> StartingPolynomial<Degree>::scale(double s) const{
	StartingPolynomial q;
	q.start=start*s;
	q.p=p.scale(s);
	return q;
}
template<int Degree>
StartingPolynomial<Degree> StartingPolynomial<Degree>::shift(double s) const{
	StartingPolynomial q;
	q.start=start+s;
	q.p=p.shift(s);
	return q;
}


template<int Degree>
int StartingPolynomial<Degree>::operator < (const StartingPolynomial<Degree>& sp) const{
	if(start<sp.start){return 1;}
	else{return 0;}
}
template<int Degree>
int StartingPolynomial<Degree>::Compare(const void* v1,const void* v2){
	double d=((StartingPolynomial*)(v1))->start-((StartingPolynomial*)(v2))->start;
	if		(d<0)	{return -1;}
	else if	(d>0)	{return  1;}
	else			{return  0;}
}

/////////////////
// PPolynomial //
/////////////////
template<int Degree>
PPolynomial<Degree>::PPolynomial(void){
	polyCount=0;
	polys=NULL;
}
template<int Degree>
PPolynomial<Degree>::PPolynomial(const PPolynomial<Degree>& p){
	polyCount=0;
	polys=NULL;
	set(p.polyCount);
	memcpy(polys,p.polys,sizeof(StartingPolynomial<Degree>)*p.polyCount);
}

template<int Degree>
PPolynomial<Degree>::~PPolynomial(void){
	if(polyCount){free(polys);}
	polyCount=0;
	polys=NULL;
}
template<int Degree>
void PPolynomial<Degree>::set(StartingPolynomial<Degree>* sps,int count){
	int i,c=0;
	set(count);
	qsort(sps,count,sizeof(StartingPolynomial<Degree>),StartingPolynomial<Degree>::Compare);//按照start大小进行排序
	for( i=0 ; i<count ; i++ )
	{
		if( !c || sps[i].start!=polys[c-1].start ) polys[c++] = sps[i];
		else{polys[c-1].p+=sps[i].p;}//如果不是第一个放进array的而且两个start相等，就直接把两个多项式相加
	}
	reset( c );
}
template <int Degree>
int PPolynomial<Degree>::size(void) const{return int(sizeof(StartingPolynomial<Degree>)*polyCount);}

template<int Degree>
void PPolynomial<Degree>::set( size_t size )
{
	if(polyCount){free(polys);}
	polyCount=0;
	polys=NULL;
	polyCount=size;
	if(size){
		polys=(StartingPolynomial<Degree>*)malloc(sizeof(StartingPolynomial<Degree>)*size);
		memset(polys,0,sizeof(StartingPolynomial<Degree>)*size);
	}
}
template<int Degree>
void PPolynomial<Degree>::reset( size_t newSize )
{
	polyCount=newSize;
	polys=(StartingPolynomial<Degree>*)realloc(polys,sizeof(StartingPolynomial<Degree>)*newSize);
}

template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator = (const PPolynomial<Degree>& p){
	set(p.polyCount);
	memcpy(polys,p.polys,sizeof(StartingPolynomial<Degree>)*p.polyCount);
	return *this;
}

template<int Degree>
template<int Degree2>
PPolynomial<Degree>& PPolynomial<Degree>::operator  = (const PPolynomial<Degree2>& p){
	set(p.polyCount);
	for(int i=0;i<int(polyCount);i++){
		polys[i].start=p.polys[i].start;//单纯的复制行为
		polys[i].p=p.polys[i].p;
	}
	return *this;
}

template<int Degree>
double PPolynomial<Degree>::operator ()( double t ) const
{
	double v=0;
	for( int i=0 ; i<int(polyCount) && t>polys[i].start ; i++ ) v+=polys[i].p(t);//累加多项式的值，start必须小于t才行
	return v;
}

template<int Degree>
double PPolynomial<Degree>::integral( double tMin , double tMax ) const
{
	int m=1;
	double start,end,s,v=0;
	start=tMin;
	end=tMax;
	if(tMin>tMax){
		m=-1;
		start=tMax;
		end=tMin;
	}
	for(int i=0;i<int(polyCount) && polys[i].start<end;i++){
		if(start<polys[i].start){s=polys[i].start;}//从这个代码来看，start应该是未知数阈值空间的最小值，即阈值范围的起点
		else{s=start;}
		v+=polys[i].p.integral(s,end);//为什么没有对end进行限制，为什么用了end而不是polys[i+1].start
	}
	return v*m;
}

//这些start只表示了每个多项式的最小取值，而各区间之间好像是连续定义的，如同BSpline的base function一样，后一个的start就是上一个结束
//如果后一个start代表了上一个函数区间定义的最大值，为什么积分的时候还需要积分到最后一个poly的start
template<int Degree>
double PPolynomial<Degree>::Integral(void) const{return integral(polys[0].start,polys[polyCount-1].start);}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator + (const PPolynomial<Degree>& p) const{
	PPolynomial q;
	int i,j;
	size_t idx=0;
	q.set(polyCount+p.polyCount);
	i=j=-1;

	while(idx<q.polyCount){
		if		(j>=int(p.polyCount)-1)				{q.polys[idx]=  polys[++i];}//先自加后结算，这里如果i大于原来的polyCount，说明原来的多项式已经全有序复制了
		//剩下的都在p表示的PPolynomial中，而且是有序的(start 序)，可以直接全复制过来，j大于p.polyCount时也一样，在此之前，要比较当前剩余多项式
		//中start与p剩余多项式中start的顺序
		else if	(i>=int(  polyCount)-1)				{q.polys[idx]=p.polys[++j];}
		else if(polys[i+1].start<p.polys[j+1].start){q.polys[idx]=  polys[++i];}//谁小谁先来
		else										{q.polys[idx]=p.polys[++j];}
		idx++;
	}
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator - (const PPolynomial<Degree>& p) const{
	PPolynomial q;//两个表示了多个多项式之间的相减，要按照start大小排序，难道两个中都没有start相等的情况出现
	int i,j;
	size_t idx=0;
	q.set(polyCount+p.polyCount);
	i=j=-1;

	while(idx<q.polyCount){
		if		(j>=int(p.polyCount)-1)				{q.polys[idx]=  polys[++i];}
		else if	(i>=int(  polyCount)-1)				{q.polys[idx].start=p.polys[++j].start;q.polys[idx].p=p.polys[j].p*(-1.0);}
		else if(polys[i+1].start<p.polys[j+1].start){q.polys[idx]=  polys[++i];}
		else										{q.polys[idx].start=p.polys[++j].start;q.polys[idx].p=p.polys[j].p*(-1.0);}
		idx++;
	}
	return q;
}
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::addScaled(const PPolynomial<Degree>& p,double scale){
	int i,j;//为什么是先比较，再scale
	StartingPolynomial<Degree>* oldPolys=polys;
	size_t idx=0,cnt=0,oldPolyCount=polyCount;
	polyCount=0;
	polys=NULL;
	set(oldPolyCount+p.polyCount);
	i=j=-1;
	while(cnt<polyCount){
		if		(j>=int( p.polyCount)-1)				{polys[idx]=oldPolys[++i];}
		else if	(i>=int(oldPolyCount)-1)				{polys[idx].start= p.polys[++j].start;polys[idx].p=p.polys[j].p*scale;}
		else if	(oldPolys[i+1].start<p.polys[j+1].start){polys[idx]=oldPolys[++i];}
		else											{polys[idx].start= p.polys[++j].start;polys[idx].p=p.polys[j].p*scale;}
		if(idx && polys[idx].start==polys[idx-1].start)	{polys[idx-1].p+=polys[idx].p;}//当start相等时还要进行合并，为什么在加减中没看到
		else{idx++;}
		cnt++;
	}
	free(oldPolys);
	reset(idx);
	return *this;
}
template<int Degree>
template<int Degree2>
PPolynomial<Degree+Degree2> PPolynomial<Degree>::operator * (const PPolynomial<Degree2>& p) const{
	PPolynomial<Degree+Degree2> q;
	StartingPolynomial<Degree+Degree2> *sp;
	int i,j,spCount=int(polyCount*p.polyCount);

	sp=(StartingPolynomial<Degree+Degree2>*)malloc(sizeof(StartingPolynomial<Degree+Degree2>)*spCount);
	for(i=0;i<int(polyCount);i++){
		for(j=0;j<int(p.polyCount);j++){
			sp[i*p.polyCount+j]=polys[i]*p.polys[j];//竟然不用比start大小，乘的时候，start取了交集，这就可以用来生成文章中所说的sparse coefficient Matrix
			//但感觉这种排序下来还是会存在start的乱序
		}
	}
	q.set(sp,spCount);//所以在这里对乱序进行了整理，用了qsort
	free(sp);
	return q;
}
template<int Degree>
template<int Degree2>
PPolynomial<Degree+Degree2> PPolynomial<Degree>::operator * (const Polynomial<Degree2>& p) const{//这块与一个单独的多项式相乘
	PPolynomial<Degree+Degree2> q;
	q.set(polyCount);
	for(int i=0;i<int(polyCount);i++){
		q.polys[i].start=polys[i].start;
		q.polys[i].p=polys[i].p*p;
	}
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::scale( double s ) const
{
	PPolynomial q;
	q.set(polyCount);
	for(size_t i=0;i<polyCount;i++){q.polys[i]=polys[i].scale(s);}//直接照搬startingPolynomial
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::shift( double s ) const
{
	PPolynomial q;
	q.set(polyCount);
	for(size_t i=0;i<polyCount;i++){q.polys[i]=polys[i].shift(s);}//直接照搬startingPolynomial
	return q;
}
template<int Degree>
PPolynomial<Degree-1> PPolynomial<Degree>::derivative(void) const{
	PPolynomial<Degree-1> q;
	q.set(polyCount);
	for(size_t i=0;i<polyCount;i++){
		q.polys[i].start=polys[i].start;
		q.polys[i].p=polys[i].p.derivative();//直接照搬Polynomial
	}
	return q;
}
template<int Degree>
PPolynomial<Degree+1> PPolynomial<Degree>::integral(void) const{
	int i;
	PPolynomial<Degree+1> q;
	q.set(polyCount);
	for(i=0;i<int(polyCount);i++){
		q.polys[i].start=polys[i].start;
		q.polys[i].p=polys[i].p.integral();//直接照搬Polynomial
		q.polys[i].p-=q.polys[i].p(q.polys[i].start);//还减去了这个多项式在start处的值
	}
	return q;
}
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator  += ( double s ) {polys[0].p+=s;}//为什么只在第一个多项式上操作???，只修改了第一个多项式上的常数项
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator  -= ( double s ) {polys[0].p-=s;}//为什么只在第一个多项式上操作???，同上
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator  *= ( double s )
{
	for(int i=0;i<int(polyCount);i++){polys[i].p*=s;}//直接照搬Polynomial
	return *this;
}
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator  /= ( double s )
{
	for(size_t i=0;i<polyCount;i++){polys[i].p/=s;}//直接照搬Polynomial
	return *this;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator + ( double s ) const
{
	PPolynomial q=*this;
	q+=s;
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator - ( double s ) const
{
	PPolynomial q=*this;
	q-=s;
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator * ( double s ) const
{
	PPolynomial q=*this;
	q*=s;
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator / ( double s ) const
{
	PPolynomial q=*this;
	q/=s;
	return q;
}

template<int Degree>
void PPolynomial<Degree>::printnl(void) const{
	Polynomial<Degree> p;

	if(!polyCount){
		Polynomial<Degree> p;
		printf("[-Infinity,Infinity]\n");
	}
	else{
		for(size_t i=0;i<polyCount;i++){
			printf("[");
			if		(polys[i  ].start== DBL_MAX){printf("Infinity,");}//足以证明两个相邻多项式之间在取值空间上的前后关系
			else if	(polys[i  ].start==-DBL_MAX){printf("-Infinity,");}
			else								{printf("%f,",polys[i].start);}
			if(i+1==polyCount)					{printf("Infinity]\t");}
			//最后一个多项式的上限是infinity，在这个区间内累加了所有的多项式，而在前面的区间中，只包含了PPolynomial中的部分多项式
			//这里确实证明了PPolynomial中不同多项式的start的前后关系，但是下一个start并不是上一个多项式阈值范围的上限，所有的上限都是infinity
			else if (polys[i+1].start== DBL_MAX){printf("Infinity]\t");}
			else if	(polys[i+1].start==-DBL_MAX){printf("-Infinity]\t");}
			else								{printf("%f]\t",polys[i+1].start);}
			p=p+polys[i].p;
			p.printnl();
		}
	}
	printf("\n");
}
template< >
PPolynomial< 0 > PPolynomial< 0 >::BSpline( double radius )//怎么保证在Degree等于0会优先调用这个函数
{
	//这个函数式Degree=0时的boxFilter，后续Degree>0的smoothing filter需要通过不断调用MovingAverage函数来获得，如下所示
	//这个函数感觉就是07文章中写的box Filter，在|r|<0.5时函数值为1，其他地方为0
	PPolynomial q;
	q.set(2);//使得polyCount等于2

	q.polys[0].start=-radius;
	q.polys[1].start= radius;

	q.polys[0].p.coefficients[0]= 1.0;
	q.polys[1].p.coefficients[0]=-1.0;
	return q;
}
template< int Degree >
PPolynomial< Degree > PPolynomial<Degree>::BSpline( double radius )//radius默认为0.5，radius设定了smoothing filter中的半径大小
{
	return PPolynomial< Degree-1 >::BSpline().MovingAverage( radius );//BSpline采用了默认radius
}
template<int Degree>
PPolynomial<Degree+1> PPolynomial<Degree>::MovingAverage( double radius ) const//就是泊松重建中BoxFilter的卷积过程，每调用一次就是卷积一次
{
	PPolynomial<Degree+1> A;
	Polynomial<Degree+1> p;
	StartingPolynomial<Degree+1>* sps;

	sps=(StartingPolynomial<Degree+1>*)malloc(sizeof(StartingPolynomial<Degree+1>)*polyCount*2);//这里的数目也没有问题

	for(int i=0;i<int(polyCount);i++){//通过与第一版泊松重建代码的比较，这个函数应该是用来BoxFilter的卷积过程，卷积结果以StartingPolynomial表示，所以并不直观
		sps[2*i  ].start=polys[i].start-radius;
		sps[2*i+1].start=polys[i].start+radius;
		p=polys[i].p.integral()-polys[i].p.integral()(polys[i].start);
		sps[2*i  ].p=p.shift(-radius);
		sps[2*i+1].p=p.shift( radius)*-1;
	}
	A.set(sps,int(polyCount*2));//重新整理顺序，整理后polyCount有可能会变小，不再是polyCount*2
	free(sps);
	return A*1.0/(2*radius);
}
template<int Degree>
void PPolynomial<Degree>::getSolutions(double c,std::vector<double>& roots,double EPS,double min,double max) const{
	Polynomial<Degree> p;//c应该代表常数项
	std::vector<double> tempRoots;

	p.setZero();
	for(size_t i=0;i<polyCount;i++){
		p+=polys[i].p;//为什么要累加求解，不懂???好像是因为阈值区间越大，覆盖的polynomial就越多，例如在[polys[0].start, polys[1].start]上只有polys[0]有定义
		//而在[polyCount-1].start, infinity]区间上，所有的polys都有定义，与定义在八叉树上不同尺度下的BSpline函数很相似
		if(polys[i].start>max){break;}//先检查当前的start是否已经超出了max的范围，如果超出，后面其他的poly肯定也超出，所以break
		if(i<polyCount-1 && polys[i+1].start<min){continue;}//如果当前的poly start还没到范围内，就先攒着，加上后面的一块进行求解
		p.getSolutions(c,tempRoots,EPS);
		for(size_t j=0;j<tempRoots.size();j++){
			if(tempRoots[j]>polys[i].start && (i+1==polyCount || tempRoots[j]<=polys[i+1].start)){
				if(tempRoots[j]>min && tempRoots[j]<max){roots.push_back(tempRoots[j]);}//记录有效解
			}
		}
	}
}

template<int Degree>
void PPolynomial<Degree>::write(FILE* fp,int samples,double min,double max) const{
	fwrite(&samples,sizeof(int),1,fp);
	for(int i=0;i<samples;i++){
		double x=min+i*(max-min)/(samples-1);
		float v=(*this)(x);//求值，写入文件
		fwrite(&v,sizeof(float),1,fp);
	}
}
