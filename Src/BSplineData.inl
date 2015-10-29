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

/////////////////
// BSplineData //
/////////////////
// Support[i]:
//		Odd:  i +/- 0.5 * ( 1 + Degree )
//			i - 0.5 * ( 1 + Degree ) < 0
// <=>		i < 0.5 * ( 1 + Degree )
//			i + 0.5 * ( 1 + Degree ) > 0
// <=>		i > - 0.5 * ( 1 + Degree )
//			i + 0.5 * ( 1 + Degree ) > r
// <=>      i > r - 0.5 * ( 1 + Degree )
//			i - 0.5 * ( 1 + Degree ) < r
// <=>      i < r + 0.5 * ( 1 + Degree )
//		Even: i + 0.5 +/- 0.5 * ( 1 + Degree )
//			i - 0.5 * Degree < 0
// <=>		i < 0.5 * Degree
//			i + 1 + 0.5 * Degree > 0
// <=>		i > -1 - 0.5 * Degree
//			i + 1 + 0.5 * Degree > r
// <=>		i > r - 1 - 0.5 * Degree
//			i - 0.5 * Degree < r
// <=>		i < r + 0.5 * Degree
template< int Degree > inline bool LeftOverlap( unsigned int depth , int offset )//没用到
{
	offset <<= 1;
	if( Degree & 1 ) return (offset < 1+Degree) && (offset > -1-Degree );//Degree&1根据上面的解释应该在判断Degree的奇偶
	else             return (offset <   Degree) && (offset > -2-Degree );
}
template< int Degree > inline bool RightOverlap( unsigned int depth , int offset )//没用到
{
	offset <<= 1;
	int r = 1<<(depth+1);
	if( Degree & 1 ) return (offset > 2-1-Degree) && (offset < 2+1+Degree );
	else             return (offset > 2-2-Degree) && (offset < 2+  Degree );
}
template< int Degree > inline int ReflectLeft( unsigned int depth , int offset )//没用到
{
	if( Degree&1 ) return   -offset;
	else           return -1-offset;
}
template< int Degree > inline int ReflectRight( unsigned int depth , int offset )//没用到
{
	int r = 1<<(depth+1);
	if( Degree&1 ) return r  -offset;
	else           return r-1-offset;
}

template< int Degree >
BSplineData< Degree >::BSplineData( void )
{
	functionCount = sampleCount = 0;
	SetBSplineElementIntegrals< Degree   , Degree   >( _vvIntegrals );
	SetBSplineElementIntegrals< Degree   , Degree-1 >( _vdIntegrals );//出现d的地方就下降一个等级
	SetBSplineElementIntegrals< Degree-1 , Degree   >( _dvIntegrals );
	SetBSplineElementIntegrals< Degree-1 , Degree-1 >( _ddIntegrals );
}

template< int Degree >
double BSplineData< Degree >::Integrator::dot( int depth , int off1 , int off2 , bool d1 , bool d2 , bool childParent ) const
{
	if( depth<0 || depth>=int( iTables.size() ) ) return 0.;
	const typename Integrator::IntegralTables& iTable = iTables[depth];
	if( childParent )
	{
		int c = off1&1;
		off1 >>= 1 , depth--;
		int ii , d = off2-off1 , res = (1<<depth);
		if( depth<0 || off1<0 || off2<0 || off1>=res || off2>=res || d<-Degree || d>Degree ) return 0;
		if     ( off1<     Degree ) ii = off1;
		else if( off1>=res-Degree ) ii = 2*Degree + off1 - (res-1);
		else                        ii = Degree;
		if     ( d1 && d2 ) return iTable.dd_cpIntegrals[2*ii+c][d+Degree];
		else if( d1       ) return iTable.dv_cpIntegrals[2*ii+c][d+Degree];
		else if(       d2 ) return iTable.vd_cpIntegrals[2*ii+c][d+Degree];
		else                return iTable.vv_cpIntegrals[2*ii+c][d+Degree];
	}
	else
	{
		int ii , d = off2-off1 , res = (1<<depth);
		if( off1<0 || off2<0 || off1>=res || off2>=res || d<-Degree || d>Degree ) return 0;//首先检查offset，res是当前深度下所有node的最大值吧，
		//还要保证这两个BSpline function相差在Degree之间，否则二者没有交集
		if     ( off1<     Degree ) ii = off1;//小于Degree使得ii=off1，一共有Degree个
		else if( off1>=res-Degree ) ii = 2*Degree + off1 - (res-1);//这里相当于在判断off1是否在两个边上，换算后相当于ii>=Degree+1，从res-Degree到res-1一共有Degree个
		else                        ii = Degree;//这有单独一个，加起来就是iTable的数目2*Degree+1个
		if     ( d1 && d2 ) return iTable.dd_ccIntegrals[ii][d+Degree];//两个都是derivative???
		else if( d1       ) return iTable.dv_ccIntegrals[ii][d+Degree];//一个是derivative
		else if(       d2 ) return iTable.vd_ccIntegrals[ii][d+Degree];
		else                return iTable.vv_ccIntegrals[ii][d+Degree];//两个都不是derivative
	}
}
template< int Degree >
template< int Radius >
double BSplineData< Degree >::CenterEvaluator< Radius >::value( int depth , int off1 , int off2 , bool d , bool childParent ) const
{
	//这里的childParent是一个Flag，标志两个node是否在两个不同的Depth下
	if( depth<0 || depth>=int( vTables.size() ) ) return 0.;
	if( childParent )
	{
		int c = off1&1;//保留off1的最后一位，这一位代表了在当前Depth划分时的offset，下面除以2求其parent index的时候就会把这个信息抹掉
		off1 >>= 1 , depth--;//为什么肯定off1就比off2多了一层，非要在深度上减一
		const typename CenterEvaluator::ValueTables& vTable = vTables[depth];
		int ii , dd = off1-off2 , res = (1<<depth);
		if( depth<0 || off1<0 || off2<0 || off1>=res || off2>=res || dd<-Radius || dd>Radius ) return 0;
		if     ( off2<     Degree ) ii = off2;
		else if( off2>=res-Degree ) ii = 2*Degree + off2 - (res-1);//估计中间符合条件的只占了ii=Degree那一个数组，剩下都分散了
		else                        ii = Degree;
		if( d ) return vTable.dValues[ii][(dd+Radius)*3+2*c];//明白了，之所以有*3的约束，index为1时放了在同一个Depth下的数据，见下面的code
		//而index为0或者2时放了两个child的数据
		else    return vTable.vValues[ii][(dd+Radius)*3+2*c];
	}
	else
	{
		const typename CenterEvaluator::ValueTables& vTable = vTables[depth];
		int ii , dd = off1-off2 , res = (1<<depth);
		if( off1<0 || off2<0 || off1>=res || off2>=res || dd<-Radius || dd>Radius ) return 0;
		if     ( off2<     Degree ) ii = off2;
		else if( off2>=res-Degree ) ii = 2*Degree + off2 - (res-1);
		else                        ii = Degree;
		if( d ) return vTable.dValues[ii][(dd+Radius)*3+1];//乘以3之后，只取第一个，剩下的两个是两个child的数据
		else    return vTable.vValues[ii][(dd+Radius)*3+1];
	}
}
template< int Degree >
template< int Radius >
double BSplineData< Degree >::CornerEvaluator< Radius >::value( int depth , int off1 , int c1 , int off2 , bool d , bool childParent ) const
{
	if( c1<0 || c1>=2 )
	{
		fprintf( stderr , "[WARNING] Clamping corner to {0,1}\n" );
		c1 = std::max< int >( 0 , std::min< int >( c1 , 1 ) );
	}
	if( depth<0 || depth>=int( vTables.size() ) ) return 0.;
	if( childParent )
	{
		int c = off1&1;
		off1 >>= 1 , depth--;
		const typename CornerEvaluator::ValueTables& vTable = vTables[depth];
		int ii , dd = off1-off2 , res = (1<<depth);
		if( depth<0 || off1<0 || off2<0 || off1>=res || off2>=res || dd<-Radius || dd>Radius ) return 0;
		if     ( off2<     Degree ) ii = off2;
		else if( off2>=res-Degree ) ii = 2*Degree + off2 - (res-1);
		else                        ii = Degree;
		if( d ) return vTable.dValues[ii][(dd+Radius)*2+c+c1];
		else    return vTable.vValues[ii][(dd+Radius)*2+c+c1];
	}
	else
	{
		const typename CornerEvaluator::ValueTables& vTable = vTables[depth];
		int ii , dd = off1-off2 , res = (1<<depth);
		if( off1<0 || off2<0 || off1>=res || off2>=res || dd<-Radius || dd>Radius ) return 0;
		if     ( off2<     Degree ) ii = off2;
		else if( off2>=res-Degree ) ii = 2*Degree + off2 - (res-1);
		else                        ii = Degree;
		if( d ) return vTable.dValues[ii][(dd+Radius)*2+2*c1];
		else    return vTable.vValues[ii][(dd+Radius)*2+2*c1];
	}
}
template< int Degree >
void BSplineData< Degree >::set( int maxDepth , int boundaryType )
{
	_boundaryType = boundaryType>0 ? 1 : ( boundaryType<0 ? -1 : 0 );

	depth = maxDepth;
	// [Warning] This assumes that the functions spacing is dual，不明白
	//难道在binaryNode的每个维度上都构建了functionCount个baseFunction
	functionCount = BinaryNode::CumulativeCenterCount( depth );//binary node是octree在每个单一维度上的简化情况，而function count是累积的???
	sampleCount   = BinaryNode::CenterCount( depth ) + BinaryNode::CornerCount( depth );//B样条曲线的可采样点就是这些center加上corner，这些点有vertex value吧
	baseFunctions = NewPointer< PPolynomial< Degree > >( functionCount );//base function，预分配空间，PPolynomial本身就表示多个多项式的叠加
	//这里的空间应该分配给了每一个center节点，每个节点都是下面baseFunction计算出的结果
	baseBSplines = NewPointer< BSplineComponents >( functionCount );//预分配空间，BSplineComponents与PPolynomial有啥区别
	//BSplineComponents没有指明自变量空间

	baseFunction = PPolynomial< Degree >::BSpline();//Degree等于2的BSpline，作为平滑方程，经过推导，应该就是box Filter进行三次卷积后形成的平滑方程
	for( int i=0 ; i<=Degree ; i++ ) baseBSpline[i] = Polynomial< Degree >::BSplineComponent( i ).shift( double(-(Degree+1)/2) + i - 0.5 );//为什么shift这么多???分别是-1.5，-0.5，0.5
	dBaseFunction = baseFunction.derivative();//basefunction的导数
	StartingPolynomial< Degree > sPolys[Degree+4];

	for( int i=0 ; i<Degree+3 ; i++ )
	{
		sPolys[i].start = double(-(Degree+1)/2) + i - 1.5;
		sPolys[i].p *= 0;
		if(         i<=Degree   )  sPolys[i].p += baseBSpline[i  ].shift( -1 ) * _boundaryType;
		if( i>=1 && i<=Degree+1 )  sPolys[i].p += baseBSpline[i-1];
		for( int j=0 ; j<i ; j++ ) sPolys[i].p -= sPolys[j].p;
	}
	leftBaseFunction.set( sPolys , Degree+3 );
	for( int i=0 ; i<Degree+3 ; i++ )
	{
		sPolys[i].start = double(-(Degree+1)/2) + i - 0.5;
		sPolys[i].p *= 0;
		if(         i<=Degree   )  sPolys[i].p += baseBSpline[i  ];
		if( i>=1 && i<=Degree+1 )  sPolys[i].p += baseBSpline[i-1].shift( 1 ) * _boundaryType;
		for( int j=0 ; j<i ; j++ ) sPolys[i].p -= sPolys[j].p;
	}
	rightBaseFunction.set( sPolys , Degree+3 );
	for( int i=0 ; i<Degree+4 ; i++ )
	{
		sPolys[i].start = double(-(Degree+1)/2) + i - 1.5;
		sPolys[i].p *= 0;
		if(         i<=Degree   )  sPolys[i].p += baseBSpline[i  ].shift( -1 ) * _boundaryType; // The left-shifted B-spline
		if( i>=1 && i<=Degree+1 )  sPolys[i].p += baseBSpline[i-1];             // The centered B-Spline
		if( i>=2 && i<=Degree+2 )  sPolys[i].p += baseBSpline[i-2].shift(  1 ) * _boundaryType; // The right-shifted B-spline
		for( int j=0 ; j<i ; j++ ) sPolys[i].p -= sPolys[j].p;
	}
	leftRightBaseFunction.set( sPolys , Degree+4 );

	dLeftBaseFunction  =  leftBaseFunction.derivative();
	dRightBaseFunction = rightBaseFunction.derivative();
	dLeftRightBaseFunction = leftRightBaseFunction.derivative();
	leftRightBSpline = leftBSpline = rightBSpline = baseBSpline;
	leftBSpline [1] +=  leftBSpline[2].shift( -1 ) ,  leftBSpline[0] *= 0;
	rightBSpline[1] += rightBSpline[0].shift(  1 ) , rightBSpline[2] *= 0;
	leftRightBSpline[1] += leftRightBSpline[2].shift( -1 ) + leftRightBSpline[0].shift( 1 ) , leftRightBSpline[0] *= 0 , leftRightBSpline[2] *= 0 ;

	double c , w;
	for( size_t i=0 ; i<functionCount ; i++ )
	{
		BinaryNode::CenterAndWidth( int(i) , c , w );//找出center和width
		//下面这两种属于正常情况下的baseFunction和baseSpline
		baseFunctions[i] = baseFunction.scale(w).shift(c);//base function，width是w，center是c，对基函数进行缩放和平移，这里相当于把平滑方程平移到每一个中心点位置，且平滑尺度大小不一
		baseBSplines[i] = baseBSpline.scale(w).shift(c);//base BSpline function也一样
		if( _boundaryType )//if(-1)是返回true
		{
			//这里应该是对边界情况单独进行处理，off=r-1，而且r=1<<d，那么r-1就是深度为d下的最后一个node
			int d , off , r;
			BinaryNode::DepthAndOffset( int(i) , d , off );//找出划分深度和在同一深度下的offset
			r = 1<<d;
			if     ( off==0 && off==r-1 ) baseFunctions[i] = leftRightBaseFunction.scale(w).shift(c);//那是不是意味着r=1，d=0，leftright是用在这种情况下的
			else if( off==0             ) baseFunctions[i] =      leftBaseFunction.scale(w).shift(c);//最左边，只有leftBase
			else if(           off==r-1 ) baseFunctions[i] =     rightBaseFunction.scale(w).shift(c);//最右边，只有rightBase
			if     ( off==0 && off==r-1 ) baseBSplines [i] = leftRightBSpline.scale(w).shift(c);//同上???
			else if( off==0             ) baseBSplines [i] =      leftBSpline.scale(w).shift(c);
			else if(           off==r-1 ) baseBSplines [i] =     rightBSpline.scale(w).shift(c);
		}
	}
}
template< int Degree >
double BSplineData< Degree >::dot( int depth1 ,  int off1 , int depth2 , int off2 , bool d1 , bool d2 , bool inset ) const//两个BSplineData的点乘
{
	const int _Degree1 = (d1 ? (Degree-1) : Degree) , _Degree2 = (d2 ? (Degree-1) : Degree);//d1或者d2为true代表是为derivative积分，所以degree要减一
	int sums[ Degree+1 ][ Degree+1 ];

	int depth = std::max< int >( depth1 , depth2 );
	//inset在这里代表边界条件
	BSplineElements< Degree > b1( 1<<depth1 , off1 , _boundaryType , inset ? ( 1<<(depth1-2) ) : 0 ) , b2( 1<<depth2 , off2 , _boundaryType , inset ? ( 1<<(depth2-2) ) : 0 );

	BSplineElements< Degree > b;//depth不够，要做自循环
	while( depth1<depth ) b=b1 , b.upSample( b1 ) , depth1++;
	while( depth2<depth ) b=b2 , b.upSample( b2 ) , depth2++;

	BSplineElements< Degree-1 > db1 , db2;
	b1.differentiate( db1 ) , b2.differentiate( db2 );

	int start1=-1 , end1=-1 , start2=-1 , end2=-1;
	for( int i=0 ; i<int( b1.size() ) ; i++ ) for( int j=0 ; j<=Degree ; j++ )
	{
		if( b1[i][j] && start1==-1 ) start1 = i;
		if( b1[i][j] ) end1 = i+1;
		if( b2[i][j] && start2==-1 ) start2 = i;
		if( b2[i][j] ) end2 = i+1;
	}
	if( start1==end1 || start2==end2 || start1>=end2 || start2>=end1 ) return 0.;
	int start = std::max< int >( start1 , start2 ) , end = std::min< int >( end1 , end2 );
	memset( sums , 0 , sizeof( sums ) );
	for( int i=start ; i<end ; i++ ) for( int j=0 ; j<=_Degree1 ; j++ ) for( int k=0 ; k<=_Degree2 ; k++ ) sums[j][k] += ( d1 ?  db1[i][j] : b1[i][j] ) * ( d2 ? db2[i][k] : b2[i][k] );
	double _dot = 0;
	if     ( d1 && d2 ) for( int j=0 ; j<=_Degree1 ; j++ ) for( int k=0 ; k<=_Degree2 ; k++ ) _dot += _ddIntegrals[j][k] * sums[j][k];
	else if( d1       ) for( int j=0 ; j<=_Degree1 ; j++ ) for( int k=0 ; k<=_Degree2 ; k++ ) _dot += _dvIntegrals[j][k] * sums[j][k];
	else if(       d2 ) for( int j=0 ; j<=_Degree1 ; j++ ) for( int k=0 ; k<=_Degree2 ; k++ ) _dot += _vdIntegrals[j][k] * sums[j][k];
	else                for( int j=0 ; j<=_Degree1 ; j++ ) for( int k=0 ; k<=_Degree2 ; k++ ) _dot += _vvIntegrals[j][k] * sums[j][k];
	_dot /= b1.denominator;
	_dot /= b2.denominator;
	if     ( d1 && d2 ) return _dot * (1<<depth);
	else if( d1 || d2 ) return _dot;
	else                return _dot / (1<<depth);
}
template< int Degree >
double BSplineData< Degree >::value( int depth ,  int off , double smoothingRadius ,  double s , bool d , bool inset ) const
{
	PPolynomial< Degree+1 >  function;
	PPolynomial< Degree   > dFunction;

	if( off<0 || off>=(1<<depth) ) return 0;
	int idx = BinaryNode::CenterIndex( depth , off );

	if( smoothingRadius>0 ) function = baseFunctions[idx].MovingAverage( smoothingRadius );
	else                    function = baseFunctions[idx];
	dFunction = function.derivative();

	if( d ) return dFunction(s);
	else    return  function(s);
}
template< int Degree >
void BSplineData< Degree >::setIntegrator( Integrator& integrator , bool inset , bool useDotRatios ) const
{
	integrator.iTables.resize( depth+1 );//估计要为每一层设置integrator table，i和j一个是从[0,2*Degree]，另一个是从[-Degree, Degree]
	for( int d=0 ; d<=depth ; d++ ) for( int i=0 ; i<=2*Degree ; i++ ) for( int j=-Degree ; j<=Degree ; j++ )
	{
		int res = 1<<d , ii = (i<=Degree ? i : i+res-1 - 2*Degree );
		integrator.iTables[d].vv_ccIntegrals[i][j+Degree] = dot( d , ii , d , ii+j , false , false , inset );
		integrator.iTables[d].dv_ccIntegrals[i][j+Degree] = dot( d , ii , d , ii+j , true  , false , inset );
		integrator.iTables[d].vd_ccIntegrals[i][j+Degree] = dot( d , ii , d , ii+j , false , true  , inset );
		integrator.iTables[d].dd_ccIntegrals[i][j+Degree] = dot( d , ii , d , ii+j , true  , true  , inset );
		if( useDotRatios )
		{
			integrator.iTables[d].dv_ccIntegrals[i][j+Degree] /= integrator.iTables[d].vv_ccIntegrals[i][j+Degree];
			integrator.iTables[d].vd_ccIntegrals[i][j+Degree] /= integrator.iTables[d].vv_ccIntegrals[i][j+Degree];
			integrator.iTables[d].dd_ccIntegrals[i][j+Degree] /= integrator.iTables[d].vv_ccIntegrals[i][j+Degree];
		}
	}
	for( int d=1 ; d<=depth ; d++ ) for( int i=0 ; i<=2*Degree ; i++ ) for( int j=-Degree ; j<=Degree ; j++ )
	{
		int res = 1<<d , ii = (i<=Degree ? i : i+(res/2)-1 - 2*Degree );
		for( int c=0 ; c<2 ; c++ )//当前child node加0或1
		{
			integrator.iTables[d].vv_cpIntegrals[2*i+c][j+Degree] = dot( d , 2*ii+c , d-1 , ii+j , false , false , inset );
			integrator.iTables[d].dv_cpIntegrals[2*i+c][j+Degree] = dot( d , 2*ii+c , d-1 , ii+j , true  , false , inset );
			integrator.iTables[d].vd_cpIntegrals[2*i+c][j+Degree] = dot( d , 2*ii+c , d-1 , ii+j , false , true  , inset );
			integrator.iTables[d].dd_cpIntegrals[2*i+c][j+Degree] = dot( d , 2*ii+c , d-1 , ii+j , true  , true  , inset );
			if( useDotRatios )
			{
				integrator.iTables[d].dv_cpIntegrals[2*i+c][j+Degree] /= integrator.iTables[d].vv_cpIntegrals[2*i+c][j+Degree];
				integrator.iTables[d].vd_cpIntegrals[2*i+c][j+Degree] /= integrator.iTables[d].vv_cpIntegrals[2*i+c][j+Degree];
				integrator.iTables[d].dd_cpIntegrals[2*i+c][j+Degree] /= integrator.iTables[d].vv_cpIntegrals[2*i+c][j+Degree];
			}
		}
	}
}
template< int Degree >
template< int Radius >
void BSplineData< Degree >::setCenterEvaluator( CenterEvaluator< Radius >& evaluator , double smoothingRadius , double dSmoothingRadius , bool inset ) const
{
	evaluator.vTables.resize( depth+1 );
	for( int d=0 ; d<=depth ; d++ ) for( int i=0 ; i<=2*Degree ; i++ ) for( int j=-Radius ; j<=Radius ; j++ ) for( int k=-1 ; k<=1 ; k++ )
	{
		int res = 1<<d , ii = (i<=Degree ? i : i+res-1 - 2*Degree );
		double s = 0.5+ii+j+0.25*k;//这是什么公式???
		evaluator.vTables[d].vValues[i][(j+Radius)*3+(k+1)] = value( d , ii ,  smoothingRadius , s/res , false , inset );
		evaluator.vTables[d].dValues[i][(j+Radius)*3+(k+1)] = value( d , ii , dSmoothingRadius , s/res , true  , inset );
	}
}
template< int Degree >
template< int Radius >
void BSplineData< Degree >::setCornerEvaluator( CornerEvaluator< Radius >& evaluator , double smoothingRadius , double dSmoothingRadius , bool inset ) const
{
	evaluator.vTables.resize( depth+1 );
	for( int d=0 ; d<=depth ; d++ ) for( int i=0 ; i<=2*Degree ; i++ ) for( int j=-Radius ; j<=Radius ; j++ ) for( int k=0 ; k<=2 ; k++ )
	{
		int res = 1<<d , ii = (i<=Degree ? i : i+res-1 - 2*Degree );
		double s = ii+j+0.5*k;
		evaluator.vTables[d].vValues[i][(j+Radius)*2+k] = value( d , ii ,  smoothingRadius , s/res , false , inset );
		evaluator.vTables[d].dValues[i][(j+Radius)*2+k] = value( d , ii , dSmoothingRadius , s/res , true  , inset );
	}
}


template< int Degree >
template< class Real >
BSplineData< Degree >::DotTables< Real >::DotTables( void )
{
	vvDotTable = NullPointer( Real );	
	dvDotTable = NullPointer( Real );	
	ddDotTable = NullPointer( Real );	
}
template< int Degree >
template< class Real >
BSplineData< Degree >::DotTables< Real >::~DotTables( void )
{
	DeletePointer( vvDotTable ); 
	DeletePointer( dvDotTable ); 
	DeletePointer( ddDotTable ); 
}
template< int Degree >
template< class Real >
inline size_t BSplineData< Degree >::DotTables< Real >::Index( int i1 , int i2 ) const { return size_t(i1)*functionCount + size_t(i2); }
template< int Degree >
template< class Real >
inline size_t BSplineData< Degree >::DotTables< Real >::SymmetricIndex( int i1 , int i2 )
{
	size_t _i1 = i1 , _i2 = i2;
	if( i1>i2 ) return ((_i1*_i1+i1)>>1)+_i2;
	else        return ((_i2*_i2+i2)>>1)+_i1;
}
template< int Degree >
template< class Real >
inline int BSplineData< Degree >::DotTables< Real >::SymmetricIndex( int i1 , int i2 , size_t& index )
{
	size_t _i1 = i1 , _i2 = i2;
	if( i1<i2 )
	{
		index = ((_i2*_i2+_i2)>>1)+_i1;
		return 1;
	}
	else
	{
		index = ((_i1*_i1+_i1)>>1)+_i2;
		return 0;
	}
}
template< int Degree >
template< class Real >
typename BSplineData< Degree >::template DotTables< Real > BSplineData< Degree >::getDotTables( int flags , bool useDotRatios , bool inset ) const
{
	typename BSplineData< Degree >::template DotTables< Real > dTables;
	dTables.functionCount = functionCount;

	size_t size = ( functionCount*functionCount + functionCount )>>1;
	size_t fullSize = functionCount*functionCount;
	if( flags & VV_DOT_FLAG )
	{
		dTables.vvDotTable = NewPointer< Real >( size );
		memset( dTables.vvDotTable , 0 , sizeof(Real)*size );
	}
	if( flags & DV_DOT_FLAG )
	{
		dTables.dvDotTable = NewPointer< Real >( fullSize );
		memset( dTables.dvDotTable , 0 , sizeof(Real)*fullSize );
	}
	if( flags & DD_DOT_FLAG )
	{
		dTables.ddDotTable = NewPointer< Real >( size );
		memset( dTables.ddDotTable , 0 , sizeof(Real)*size );
	}
	int vvSums[Degree+1][Degree+1];
	int vdSums[Degree+1][Degree  ];
	int dvSums[Degree  ][Degree+1];
	int ddSums[Degree  ][Degree  ];
	double vvIntegrals[Degree+1][Degree+1];
	double vdIntegrals[Degree+1][Degree  ];
	double dvIntegrals[Degree  ][Degree+1];
	double ddIntegrals[Degree  ][Degree  ];
	SetBSplineElementIntegrals< Degree   , Degree   >( vvIntegrals );
	SetBSplineElementIntegrals< Degree   , Degree-1 >( vdIntegrals );
	SetBSplineElementIntegrals< Degree-1 , Degree   >( dvIntegrals );
	SetBSplineElementIntegrals< Degree-1 , Degree-1 >( ddIntegrals );

	for( int d1=0 ; d1<=depth ; d1++ ) for( int off1=0 ; off1<(1<<d1) ; off1++ )
	{
		int ii = BinaryNode::CenterIndex( d1 , off1 );
		BSplineElements< Degree > b1( 1<<d1 , off1 , _boundaryType , inset ? ( 1<<(d1-2) ) : 0 );
		BSplineElements< Degree-1 > db1;
		b1.differentiate( db1 );
		int start1 , end1;

		start1 = -1 , end1 = -1;
		for( int i=0 ; i<int(b1.size()) ; i++ ) for( int j=0 ; j<=Degree ; j++ )
		{
			if( b1[i][j] && start1==-1 ) start1 = i;
			if( b1[i][j] ) end1 = i+1;
		}
		if( start1==end1 ) continue;
		for( int d2=d1 ; d2<=depth ; d2++ )
		{
			for( int off2=0 ; off2<(1<<d2) ; off2++ )
			{
				int start2 = off2-Degree;
				int end2   = off2+Degree+1;
				if( start2>=end1 || start1>=end2 ) continue;
				start2 = std::max< int >( start1 , start2 );
				end2   = std::min< int >(   end1 ,   end2 );
				if( d1==d2 && off2<off1 ) continue;
				int jj = BinaryNode::CenterIndex( d2 , off2 );
				BSplineElements< Degree > b2( 1<<d2 , off2 , _boundaryType , inset ? ( 1<<(d2-2) ) : 0 );
				BSplineElements< Degree-1 > db2;
				b2.differentiate( db2 );

				size_t idx = DotTables< Real >::SymmetricIndex( ii , jj );
				size_t idx1 = DotTables< Real >::Index( ii , jj ) , idx2 = DotTables< Real >::Index( jj , ii );

				memset( vvSums , 0 , sizeof( int ) * ( Degree+1 ) * ( Degree+1 ) );
				memset( vdSums , 0 , sizeof( int ) * ( Degree+1 ) * ( Degree   ) );
				memset( dvSums , 0 , sizeof( int ) * ( Degree   ) * ( Degree+1 ) );
				memset( ddSums , 0 , sizeof( int ) * ( Degree   ) * ( Degree   ) );
				for( int i=start2 ; i<end2 ; i++ )
				{
					for( int j=0 ; j<=Degree ; j++ ) for( int k=0 ; k<=Degree ; k++ ) vvSums[j][k] +=  b1[i][j] *  b2[i][k];
					for( int j=0 ; j<=Degree ; j++ ) for( int k=0 ; k< Degree ; k++ ) vdSums[j][k] +=  b1[i][j] * db2[i][k];
					for( int j=0 ; j< Degree ; j++ ) for( int k=0 ; k<=Degree ; k++ ) dvSums[j][k] += db1[i][j] *  b2[i][k];
					for( int j=0 ; j< Degree ; j++ ) for( int k=0 ; k< Degree ; k++ ) ddSums[j][k] += db1[i][j] * db2[i][k];
				}
				double vvDot = 0 , dvDot = 0 , vdDot = 0 , ddDot = 0;
				for( int j=0 ; j<=Degree ; j++ ) for( int k=0 ; k<=Degree ; k++ ) vvDot += vvIntegrals[j][k] * vvSums[j][k];
				for( int j=0 ; j<=Degree ; j++ ) for( int k=0 ; k< Degree ; k++ ) vdDot += vdIntegrals[j][k] * vdSums[j][k];
				for( int j=0 ; j< Degree ; j++ ) for( int k=0 ; k<=Degree ; k++ ) dvDot += dvIntegrals[j][k] * dvSums[j][k];
				for( int j=0 ; j< Degree ; j++ ) for( int k=0 ; k< Degree ; k++ ) ddDot += ddIntegrals[j][k] * ddSums[j][k];
				vvDot /= (1<<d2);
				ddDot *= (1<<d2);
				vvDot /= ( b1.denominator * b2.denominator );
				dvDot /= ( b1.denominator * b2.denominator );
				vdDot /= ( b1.denominator * b2.denominator );
				ddDot /= ( b1.denominator * b2.denominator );
				if( fabs(vvDot)<1e-15 ) continue;
				if( flags & VV_DOT_FLAG ) dTables.vvDotTable[idx] = Real( vvDot );
				if( flags & DV_DOT_FLAG ) dTables.dvDotTable[idx1] = Real( dvDot );
				if( flags & DV_DOT_FLAG ) dTables.dvDotTable[idx2] = Real( dvDot );
				if( flags & DD_DOT_FLAG ) dTables.ddDotTable[idx ] = Real( ddDot );
			}
			BSplineElements< Degree > b;
			b = b1;
			b.upSample( b1 );
			b1.differentiate( db1 );
			start1 = -1;
			for( int i=0 ; i<int(b1.size()) ; i++ ) for( int j=0 ; j<=Degree ; j++ )
			{
				if( b1[i][j] && start1==-1 ) start1 = i;
				if( b1[i][j] ) end1 = i+1;
			}
		}
	}
	return dTables;
}
template< int Degree >
template< class Real >
BSplineData< Degree >::ValueTables< Real >::ValueTables( void )
{
	valueTable = NullPointer( Real );
	dValueTable = NullPointer( Real );
}
template< int Degree >
template< class Real >
BSplineData< Degree >::ValueTables< Real >::~ValueTables( void )
{
	DeletePointer( valueTable ); 
	DeletePointer( dValueTable ); 
}
template< int Degree >
template< class Real >
inline size_t BSplineData< Degree >::ValueTables< Real >::Index( int i1 , int i2 ) const { return size_t(i1)*functionCount + size_t(i2); }
template< int Degree >
template< class Real >
typename BSplineData< Degree >::template ValueTables< Real > BSplineData< Degree >::getValueTables( int flags , double valueSmooth , double derivativeSmooth ) const
{
	typename BSplineData< Degree >::template ValueTables< Real > vTables;
	vTables.functionCount = functionCount;
	vTables.sampleCount = sampleCount;

	if( flags &   VALUE_FLAG ) vTables.valueTable = NewPointer< Real >( functionCount*sampleCount );
	if( flags & D_VALUE_FLAG ) vTables.dValueTable = NewPointer< Real >( functionCount*sampleCount );
	PPolynomial< Degree+1 >  function;
	PPolynomial< Degree   > dFunction;
	for( size_t i=0 ; i<functionCount ; i++ )
	{
		if( valueSmooth>0 )      function=baseFunctions[i].MovingAverage( valueSmooth );
		else                     function=baseFunctions[i];
		if( derivativeSmooth>0 ) dFunction=baseFunctions[i].derivative().MovingAverage( derivativeSmooth );
		else                     dFunction=baseFunctions[i].derivative();

		for( size_t j=0 ; j<sampleCount ; j++ )
		{
			double x=double(j)/(sampleCount-1);
			if( flags &   VALUE_FLAG ) vTables.valueTable[j*functionCount+i] = Real( function(x));
			if( flags & D_VALUE_FLAG ) vTables.dValueTable[j*functionCount+i] = Real(dFunction(x));
		}
	}
	return vTables;
}
template< int Degree >
template< class Real >
void BSplineData< Degree >::ValueTables< Real >::setSampleSpan( int idx , int& start , int& end , double smooth ) const
{
	int d , off , res;
	BinaryNode::DepthAndOffset( idx , d , off );
	res = 1<<d;
	double _start = ( off + 0.5 - 0.5*(Degree+1) ) / res - smooth;
	double _end   = ( off + 0.5 + 0.5*(Degree+1) ) / res + smooth;
	//   (start)/(sampleCount-1) >_start && (start-1)/(sampleCount-1)<=_start
	// => start > _start * (sampleCount-1 ) && start <= _start*(sampleCount-1) + 1
	// => _start * (sampleCount-1) + 1 >= start > _start * (sampleCount-1)
	start = int( floor( _start * (sampleCount-1) + 1 ) );
	if( start<0 ) start = 0;
	//   (end)/(sampleCount-1)<_end && (end+1)/(sampleCount-1)>=_end
	// => end < _end * (sampleCount-1 ) && end >= _end*(sampleCount-1) - 1
	// => _end * (sampleCount-1) > end >= _end * (sampleCount-1) - 1
	end = int( ceil( _end * (sampleCount-1) - 1 ) );
	if( end>=int(sampleCount) ) end = int(sampleCount)-1;
}


/////////////////////
// BSplineElements //
/////////////////////
template< int Degree >
BSplineElements< Degree >::BSplineElements( int res , int offset , int boundary , int inset )
{
	denominator = 1;
	std::vector< BSplineElementCoefficients< Degree > >::resize( res , BSplineElementCoefficients< Degree >() );//每一个独立的BSplineElementCoefficients都代表一个独立的多项式
	//是否也是BSpline中一个独立的Element???

	for( int i=0 ; i<=Degree ; i++ )
	{
		int idx = -_off + offset + i;//_off = (Degree+1)/2 = 1，则idx= offset + i - 1
		if( idx>=0 && idx<res ) (*this)[idx][i] = 1;//初始化???
	}
	if( boundary!=0 )
	{
		_addLeft( offset-2*res , boundary ) , _addRight( offset+2*res , boundary );
		if( Degree&1 ) _addLeft( offset-res , boundary ) , _addRight(  offset+res     , boundary );//Degree=2的话Degree&1等于0
		else           _addLeft( -offset-1  , boundary ) , _addRight( -offset-1+2*res , boundary );
	}
	if( inset ) for( int i=0 ; i<inset && i<res ; i++ ) for( int j=0 ; j<=Degree ; j++ ) (*this)[i][j] = (*this)[res-1-i][j] = 0;//好像把两侧边界的系数全置为0
}
template< int Degree >
void BSplineElements< Degree >::_addLeft( int offset , int boundary )
{
	int res = int( std::vector< BSplineElementCoefficients< Degree > >::size() );//为啥直接这样调用
	bool set = false;
	for( int i=0 ; i<=Degree ; i++ )
	{
		int idx = -_off + offset + i;
		if( idx>=0 && idx<res ) (*this)[idx][i] += boundary , set = true;//+=boundary是什么意思???
	}
	if( set ) _addLeft( offset-2*res , boundary );
}
template< int Degree >
void BSplineElements< Degree >::_addRight( int offset , int boundary )
{
	int res = int( std::vector< BSplineElementCoefficients< Degree > >::size() );
	bool set = false;
	for( int i=0 ; i<=Degree ; i++ )
	{
		int idx = -_off + offset + i;
		if( idx>=0 && idx<res ) (*this)[idx][i] += boundary , set = true;
	}
	if( set ) _addRight( offset+2*res , boundary );
}
template< int Degree >
void BSplineElements< Degree >::upSample( BSplineElements< Degree >& high ) const
{
	fprintf( stderr , "[ERROR] B-spline up-sampling not supported for degree %d\n" , Degree );
	exit( 0 );
}
template<>
void BSplineElements< 1 >::upSample( BSplineElements< 1 >& high ) const//这...，没看懂在上采样什么东西???
{
	high.resize( size()*2 );
	high.assign( high.size() , BSplineElementCoefficients<1>() );
	for( int i=0 ; i<int(size()) ; i++ )
	{
		high[2*i+0][0] += 1 * (*this)[i][0];
		high[2*i+0][1] += 0 * (*this)[i][0];
		high[2*i+1][0] += 2 * (*this)[i][0];
		high[2*i+1][1] += 1 * (*this)[i][0];

		high[2*i+0][0] += 1 * (*this)[i][1];
		high[2*i+0][1] += 2 * (*this)[i][1];
		high[2*i+1][0] += 0 * (*this)[i][1];
		high[2*i+1][1] += 1 * (*this)[i][1];
	}
	high.denominator = denominator * 2;
}
template<>
void BSplineElements< 2 >::upSample( BSplineElements< 2 >& high ) const
{
	//    /----\
	//   /      \
	//  /        \  = 1  /--\       +3    /--\     +3      /--\   +1        /--\
	// /          \     /    \           /    \           /    \           /    \
	// |----------|     |----------|   |----------|   |----------|   |----------|

	high.resize( size()*2 );
	high.assign( high.size() , BSplineElementCoefficients<2>() );
	for( int i=0 ; i<int(size()) ; i++ )
	{
		high[2*i+0][0] += 1 * (*this)[i][0];
		high[2*i+0][1] += 0 * (*this)[i][0];
		high[2*i+0][2] += 0 * (*this)[i][0];
		high[2*i+1][0] += 3 * (*this)[i][0];
		high[2*i+1][1] += 1 * (*this)[i][0];
		high[2*i+1][2] += 0 * (*this)[i][0];

		high[2*i+0][0] += 3 * (*this)[i][1];
		high[2*i+0][1] += 3 * (*this)[i][1];
		high[2*i+0][2] += 1 * (*this)[i][1];
		high[2*i+1][0] += 1 * (*this)[i][1];
		high[2*i+1][1] += 3 * (*this)[i][1];
		high[2*i+1][2] += 3 * (*this)[i][1];

		high[2*i+0][0] += 0 * (*this)[i][2];
		high[2*i+0][1] += 1 * (*this)[i][2];
		high[2*i+0][2] += 3 * (*this)[i][2];
		high[2*i+1][0] += 0 * (*this)[i][2];
		high[2*i+1][1] += 0 * (*this)[i][2];
		high[2*i+1][2] += 1 * (*this)[i][2];
	}
	high.denominator = denominator * 4;
}

template< int Degree >
void BSplineElements< Degree >::differentiate( BSplineElements< Degree-1 >& d ) const
{
	//这里利用非静态函数的静态调用与直接用this->size()有什么区别
	d.resize( std::vector< BSplineElementCoefficients< Degree > >::size() );
	d.assign( d.size()  , BSplineElementCoefficients< Degree-1 >() );
	for( int i=0 ; i<int(std::vector< BSplineElementCoefficients< Degree > >::size()) ; i++ ) for( int j=0 ; j<=Degree ; j++ )
	{
		if( j-1>=0 )   d[i][j-1] -= (*this)[i][j];//这里同样出现了this指针，为什么上面调用size函数时不使用this指针
		if( j<Degree ) d[i][j  ] += (*this)[i][j];
	}
	d.denominator = denominator;
}
// If we were really good, we would implement this integral table to store
// rational values to improve precision...
template< int Degree1 , int Degree2 >
void SetBSplineElementIntegrals( double integrals[Degree1+1][Degree2+1] )
{
	for( int i=0 ; i<=Degree1 ; i++ )
	{
		Polynomial< Degree1 > p1 = Polynomial< Degree1 >::BSplineComponent( i );//p1的幂次为Degree1，range index为i
		for( int j=0 ; j<=Degree2 ; j++ )
		{
			Polynomial< Degree2 > p2 = Polynomial< Degree2 >::BSplineComponent( j );//p2的幂次为Degree2，range index为j
			integrals[i][j] = ( p1 * p2 ).integral( 0 , 1 );//先进行多项式相乘，然后积分，区间为[0,1]，不明白为啥，按照二维数组进行存储，i+j代表了某一乘积项的degree
		}
	}
}
