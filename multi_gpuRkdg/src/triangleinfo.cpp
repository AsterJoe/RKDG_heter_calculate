#include "../inc/triangleinfo.h"

CTriangleInfo::CTriangleInfo():
area(NULL),
perimeter(NULL),
outer_normal_vector(NULL),
mass_coeff(NULL),
vol_bf_value(NULL),
vol_bdf_value(NULL),
edge_bf_value(NULL),
vol_gauss_weight(NULL),
edge_gauss_weight(NULL),
_cell_num(0)
{

}

// 以下函数是本文件使用的临时函数

static inline double Hanshu1(double x, double y)
{
	return x * x * x * x;
}
static inline double Hanshu2(double x, double y)
{
	return x * x * x * y;
}
static inline double Hanshu3(double x, double y)
{
	return x * x * y * y;
}
static inline double Hanshu4(double x, double y)
{
	return x * y * y * y;
}
static inline double Hanshu5(double x, double y)
{
	return y * y * y * y;
}
static inline double Hanshu6(double x,double y)
{
	return x * x * x;
}
static inline double Hanshu7(double x, double y)
{
	return x * x * y;
}
static inline double Hanshu8(double x, double y)
{
	return x * y * y;
}
static inline double Hanshu9(double x, double y)
{
	return y * y * y;
}
static inline double Hanshu10(double x, double y)
{
	return x * x;
}
static inline double Hanshu11(double x, double y)
{
	return x * y;
}
static inline double Hanshu12(double x, double y)
{
	return y * y;
}
static inline double Hanshu13(double x, double y)
{
	return x;
}
static inline double Hanshu14(double x, double y)
{
	return y;
}

double CTriangleInfo::basisFunction(int tindex, int findex, double x, double y)
{
	assert(tindex>=0 && tindex<_cell_num && findex>=0 && findex<BASIS_FUNCTIONS);

	// 该三角形上的全局索引
	int gindex = tindex*BASIS_FUNCTIONS*BASIS_FUNCTION_COEFFS+findex*BASIS_FUNCTION_COEFFS;  


	double bx = barycenter[tindex].getX();
	double by = barycenter[tindex].getY();

	return (
			basis_fun_coeff[gindex]
		 +  basis_fun_coeff[gindex+1]*(x-bx)
		 +  basis_fun_coeff[gindex+2]*(y-by)
		 +  basis_fun_coeff[gindex+3]*(x-bx)*(x-bx)
		 +  basis_fun_coeff[gindex+4]*(x-bx)*(y-by)
		 +  basis_fun_coeff[gindex+5]*(y-by)*(y-by)
		 );
}

double CTriangleInfo::derivativeFunction(int tindex, int findex, double x, double y)
{
	assert(tindex>=0 && tindex<_cell_num && findex>=0 && findex<2*BASIS_FUNCTIONS);

	int bi = findex/2;	// 第几个基函数
	int dxy = findex%2; // 对x还是y求偏导

	double bx = barycenter[tindex].getX();
	double by = barycenter[tindex].getY();

	int gindex = tindex*BASIS_FUNCTIONS*BASIS_FUNCTION_COEFFS+bi*BASIS_FUNCTION_COEFFS;


	if ( 0==dxy )
	{
		return	basis_fun_coeff[gindex+1]
			  + 2*basis_fun_coeff[gindex+3]*(x-bx)
			  + basis_fun_coeff[gindex+4]*(y-by);
	}

	return    basis_fun_coeff[gindex+2]
			+ basis_fun_coeff[gindex+4]*(x-bx)
			+ 2*basis_fun_coeff[gindex+5]*(y-by);
}

CVertice2D& CTriangleInfo::getBarycenter(int index)
{
	assert(index>=0 && index<_cell_num);
	return barycenter[index];
}


CVertice2D& CTriangleInfo::getEdgeMiddleVertice(int tindex, int eindex)
{
	assert(tindex>=0 && tindex<_cell_num && eindex>=0 && eindex<TRIANGLE_EDGES);

	return edge_middle_vertice[tindex*TRIANGLE_EDGES+eindex];
}


double CTriangleInfo::getRadiusOfInscribedCircle(int index)
{
	assert(index>=0 && index<_cell_num);
	return radius_of_inscribed_circle[index];
}



void CTriangleInfo::setRadiusOfInscribedCircle(int index, double radius)
{
	assert(index>=0 && index<_cell_num);
	radius_of_inscribed_circle[index] = radius;
}


void CTriangleInfo::initInformation( int index, double x[3], double y[3] )
{
	assert( index>=0 && index<_cell_num );
	double l[TRIANGLE_VERTICES], Sxishu[14];

	//Gauss积分节点权重常数
	double w1 = 0.797426985353087;
	double w2 = 0.101286507323456;

	double w3 = 0.059715871789770;
	double w4 = 0.470142064105115;

	double cp = 0.5 + sqrt(15.0)/10.0; //线Gauss点比例系数

	int i, j;


	// 三角形面积
	area[index] = 0.5*fabs((x[2]-x[0])*(y[1]-y[0])-(y[2]-y[0])*(x[1]-x[0]));

	// 三角形三条边长度和周长
	l[0] = sqrt((x[1]-x[2])*(x[1]-x[2])+(y[1]-y[2])*(y[1]-y[2]));
	l[1] = sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]));
	l[2] = sqrt((x[0]-x[1])*(x[0]-x[1])+(y[0]-y[1])*(y[0]-y[1]));

	perimeter[index] = l[0]+l[1]+l[2];

	barycenter[index].setVertice((x[0]+x[1]+x[2])/3.0, (y[0]+y[1]+y[2])/3.0);
	radius_of_inscribed_circle[index] = 2.0*area[index]/perimeter[index];

	// 外法向量
	outer_normal_vector[index] = (y[2]-y[1])/l[0];
	outer_normal_vector[index+1*_cell_num] = (x[1]-x[2])/l[0];
	outer_normal_vector[index+2*_cell_num] = (y[0]-y[2])/l[1];
	outer_normal_vector[index+3*_cell_num] = (x[2]-x[0])/l[1];
	outer_normal_vector[index+4*_cell_num] = (y[1]-y[0])/l[2];
	outer_normal_vector[index+5*_cell_num] = (x[0]-x[1])/l[2];

	int gindex = TRIANGLE_EDGES*index;
	// 边中点
	edge_middle_vertice[gindex].setVertice((x[1]+x[2])/2.0, (y[1]+y[2])/2.0);
	edge_middle_vertice[gindex+1].setVertice((x[0]+x[2])/2.0, (y[0]+y[2])/2.0);
	edge_middle_vertice[gindex+2].setVertice((x[1]+x[0])/2.0, (y[1]+y[0])/2.0);

	gindex = VOLUME_GPOINTS*index;
	// 体高斯积分节点
	vol_gauss_vertice[gindex].setVertice((x[0]+x[1]+x[2])/3.0, (y[0]+y[1]+y[2])/3.0);
	vol_gauss_vertice[gindex+1].setVertice(w1*x[0]+w2*(x[1]+x[2]), w1*y[0]+w2*(y[1]+y[2]));
	vol_gauss_vertice[gindex+2].setVertice(w1*x[1]+w2*(x[0]+x[2]), w1*y[1]+w2*(y[0]+y[2]));
	vol_gauss_vertice[gindex+3].setVertice(w1*x[2]+w2*(x[0]+x[1]), w1*y[2]+w2*(y[0]+y[1]));
	vol_gauss_vertice[gindex+4].setVertice(w3*x[0]+w4*(x[1]+x[2]), w3*y[0]+w4*(y[1]+y[2]));
	vol_gauss_vertice[gindex+5].setVertice(w3*x[1]+w4*(x[0]+x[2]), w3*y[1]+w4*(y[0]+y[2]));
	vol_gauss_vertice[gindex+6].setVertice(w3*x[2]+w4*(x[0]+x[1]), w3*y[2]+w4*(y[0]+y[1]));

	gindex = EDGE_GPOINTS*TRIANGLE_EDGES*index;
	// 边上高斯积分节点
	edge_gauss_vertice[gindex].setVertice(cp*x[1]+(1-cp)*x[2], cp*y[1]+(1-cp)*y[2]);
	edge_gauss_vertice[gindex+1].setVertice(0.5*x[1]+0.5*x[2], 0.5*y[1]+0.5*y[2]);			
	edge_gauss_vertice[gindex+2].setVertice((1-cp)*x[1]+cp*x[2], (1-cp)*y[1]+cp*y[2]);
	edge_gauss_vertice[gindex+3].setVertice(cp*x[2]+(1-cp)*x[0], cp*y[2]+(1-cp)*y[0]);
	edge_gauss_vertice[gindex+4].setVertice(0.5*x[2]+0.5*x[0], 0.5*y[2]+0.5*y[0]);
	edge_gauss_vertice[gindex+5].setVertice((1-cp)*x[2]+cp*x[0], (1-cp)*y[2]+cp*y[0]);
	edge_gauss_vertice[gindex+6].setVertice(cp*x[0]+(1-cp)*x[1], cp*y[0]+(1-cp)*y[1]);
	edge_gauss_vertice[gindex+7].setVertice(0.5*x[0]+0.5*x[1], 0.5*y[0]+0.5*y[1]);
	edge_gauss_vertice[gindex+8].setVertice((1-cp)*x[0]+cp*x[1], (1-cp)*y[0]+cp*y[1]);

	vol_gauss_weight[index] = area[index] * 0.225000000000000;
	vol_gauss_weight[index+1*_cell_num] = area[index] * 0.125939180544827;
	vol_gauss_weight[index+2*_cell_num] = area[index] * 0.125939180544827;
	vol_gauss_weight[index+3*_cell_num] = area[index] * 0.125939180544827;
	vol_gauss_weight[index+4*_cell_num] = area[index] * 0.132394152788506;
	vol_gauss_weight[index+5*_cell_num] = area[index] * 0.132394152788506;
	vol_gauss_weight[index+6*_cell_num] = area[index] * 0.132394152788506;

	for ( i=0; i<TRIANGLE_EDGES; ++i )
	{
		edge_gauss_weight[index+(EDGE_GPOINTS*i)*_cell_num] = 0.55555555555556 * 0.5 * l[i];
		edge_gauss_weight[index+(EDGE_GPOINTS*i+1)*_cell_num] = 0.88888888888889 * 0.5 * l[i];
		edge_gauss_weight[index+(EDGE_GPOINTS*i+2)*_cell_num] = 0.55555555555556 * 0.5 * l[i];
	}

	initXishuMatrix(index, x, y, Sxishu);

	solveCoefficient(index, Sxishu);

	initMassMatrix(index);

	double xv, yv;
	// 初始化基函数在边上高斯积分点的值
	for ( i=0; i<BASIS_FUNCTIONS; ++i )
	{
		for ( j=0; j<TRIANGLE_EDGES*EDGE_GPOINTS; ++j )
		{
			xv = edge_gauss_vertice[gindex+j].getX();
			yv = edge_gauss_vertice[gindex+j].getY();

			edge_bf_value[index+j*_cell_num+i*TRIANGLE_EDGES*EDGE_GPOINTS*_cell_num] 
				= basisFunction(index, i,xv, yv);
		}
	}

	// 初始化基函数导数在体高斯积分点的值
	for ( i=0; i<2*BASIS_FUNCTIONS; ++i )
	{
		for ( j=0; j<VOLUME_GPOINTS; ++j )
		{
			xv = vol_gauss_vertice[index*VOLUME_GPOINTS+j].getX();
			yv = vol_gauss_vertice[index*VOLUME_GPOINTS+j].getY();

			vol_bdf_value[index+j*_cell_num+i*VOLUME_GPOINTS*_cell_num] = derivativeFunction(index, i, xv, yv);
		}
	}

	// 初始化基函数在高斯积分节点上的值
	for ( i=0; i<BASIS_FUNCTIONS; ++i )
	{
		for ( j=0; j<VOLUME_GPOINTS; ++j )
		{
			xv = vol_gauss_vertice[index*VOLUME_GPOINTS+j].getX();
			yv = vol_gauss_vertice[index*VOLUME_GPOINTS+j].getY();

			vol_bf_value[index+j*_cell_num+i*VOLUME_GPOINTS*_cell_num] = basisFunction(index, i, xv, yv);
		}
	}

}


void CTriangleInfo::initXishuMatrix(int index, double x[], double y[], double Sxishu[])
{
	double A,B,C,D,E,F,G;
	double w1,w2,w3,w4;
	double ss,  m_GNode[13][2];

	A = 0.065130102902216;
	B = 0.869739794195568;
	C = 0.312865496004875;
	D = 0.638444188569809;
	E = 0.048690315425316;
	F = 0.260345966079038;
	G = 0.479308067841923;

	w1 = -0.149570044467670;
	w2 =  0.053347235608839;
	w3 =  0.175615257433204;
	w4 =  0.077113760890257;

	ss = 1.0 / (sqrt(area[index]));

	//计算13个面积分Gauss节点(减重心坐标后除以面积系数)

	m_GNode[0][0] = 0;
	m_GNode[0][1] = 0;

	m_GNode[1][0] = (B * x[0] + A * (x[1] + x[2]) -barycenter[index].getX()) * ss;
	m_GNode[1][1] = (B * y[0] + A * (y[1] + y[2]) -barycenter[index].getY()) * ss;

	m_GNode[2][0] = (B * x[1] + A * (x[0] + x[2]) -barycenter[index].getX()) * ss;
	m_GNode[2][1] = (B * y[1] + A * (y[0] + y[2]) -barycenter[index].getY()) * ss;

	m_GNode[3][0] = (B *x[2] + A * (x[0] + x[1]) - barycenter[index].getX()) * ss;
	m_GNode[3][1] = (B *y[2] + A * (y[0] + y[1]) - barycenter[index].getY()) * ss;

	m_GNode[4][0] = (G * x[0] + F * (x[1] + x[2]) - barycenter[index].getX()) * ss;
	m_GNode[4][1] = (G * y[0] + F * (y[1] + y[2]) - barycenter[index].getY()) * ss;

	m_GNode[5][0] = (G * x[1] + F * (x[0] + x[2]) - barycenter[index].getX()) * ss;
	m_GNode[5][1] = (G * y[1] + F * (y[0] + y[2]) - barycenter[index].getY()) * ss;

	m_GNode[6][0] = (G * x[2] + F * (x[0] + x[1]) - barycenter[index].getX()) * ss;
	m_GNode[6][1] = (G * y[2] + F * (y[0] + y[1]) - barycenter[index].getY()) * ss;

	m_GNode[7][0] = (C * x[0] + D * x[1] + E * x[2] - barycenter[index].getX()) * ss;
	m_GNode[7][1] = (C * y[0] + D * y[1] + E * y[2] - barycenter[index].getY()) * ss;

	m_GNode[8][0] = (D * x[0] + C * x[1] + E * x[2] - barycenter[index].getX()) * ss;
	m_GNode[8][1] = (D * y[0] + C * y[1] + E * y[2] - barycenter[index].getY()) * ss;

	m_GNode[9][0] = (D * x[0] + E * x[1] + C * x[2] - barycenter[index].getX()) * ss;
	m_GNode[9][1] = (D * y[0] + E * y[1] + C * y[2] - barycenter[index].getY()) * ss;

	m_GNode[10][0] = (E * x[0] + D * x[1] + C * x[2] - barycenter[index].getX()) * ss;
	m_GNode[10][1] = (E * y[0] + D * y[1] + C * y[2] - barycenter[index].getY()) * ss;

	m_GNode[11][0] = (E * x[0] + C * x[1] + D * x[2] - barycenter[index].getX()) * ss;
	m_GNode[11][1] = (E * y[0] + C * y[1] + D * y[2] - barycenter[index].getY()) * ss;

	m_GNode[12][0] = (C * x[0] + E * x[1] + D * x[2] - barycenter[index].getX()) * ss;
	m_GNode[12][1] = (C * y[0] + E * y[1] + D * y[2] - barycenter[index].getY()) * ss;


	//计算14种系数矩阵因子

	Sxishu[0] = w1 * Hanshu1(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu1(m_GNode[1][0], m_GNode[1][1]) + Hanshu1(m_GNode[2][0], m_GNode[2][1]) + Hanshu1(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu1(m_GNode[4][0], m_GNode[4][1]) + Hanshu1(m_GNode[5][0], m_GNode[5][1]) + Hanshu1(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu1(m_GNode[7][0], m_GNode[7][1]) + Hanshu1(m_GNode[8][0], m_GNode[8][1]) + Hanshu1(m_GNode[9][0], m_GNode[9][1])
		+Hanshu1(m_GNode[10][0], m_GNode[10][1]) + Hanshu1(m_GNode[11][0], m_GNode[11][1]) + Hanshu1(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[1] = w1 * Hanshu2(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu2(m_GNode[1][0], m_GNode[1][1]) + Hanshu2(m_GNode[2][0], m_GNode[2][1]) + Hanshu2(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu2(m_GNode[4][0], m_GNode[4][1]) + Hanshu2(m_GNode[5][0], m_GNode[5][1]) + Hanshu2(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu2(m_GNode[7][0], m_GNode[7][1]) + Hanshu2(m_GNode[8][0], m_GNode[8][1]) + Hanshu2(m_GNode[9][0], m_GNode[9][1])
		+Hanshu2(m_GNode[10][0], m_GNode[10][1]) + Hanshu2(m_GNode[11][0], m_GNode[11][1]) + Hanshu2(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[2] = w1 * Hanshu3(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu3(m_GNode[1][0], m_GNode[1][1]) + Hanshu3(m_GNode[2][0], m_GNode[2][1]) + Hanshu3(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu3(m_GNode[4][0], m_GNode[4][1]) + Hanshu3(m_GNode[5][0], m_GNode[5][1]) + Hanshu3(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu3(m_GNode[7][0], m_GNode[7][1]) + Hanshu3(m_GNode[8][0], m_GNode[8][1]) + Hanshu3(m_GNode[9][0], m_GNode[9][1])
		+Hanshu3(m_GNode[10][0], m_GNode[10][1]) + Hanshu3(m_GNode[11][0], m_GNode[11][1]) + Hanshu3(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[3] = w1 * Hanshu4(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu4(m_GNode[1][0], m_GNode[1][1]) + Hanshu4(m_GNode[2][0], m_GNode[2][1]) + Hanshu4(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu4(m_GNode[4][0], m_GNode[4][1]) + Hanshu4(m_GNode[5][0], m_GNode[5][1]) + Hanshu4(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu4(m_GNode[7][0], m_GNode[7][1]) + Hanshu4(m_GNode[8][0], m_GNode[8][1]) + Hanshu4(m_GNode[9][0], m_GNode[9][1])
		+Hanshu4(m_GNode[10][0], m_GNode[10][1]) + Hanshu4(m_GNode[11][0], m_GNode[11][1]) + Hanshu4(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[4] = w1 * Hanshu5(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu5(m_GNode[1][0], m_GNode[1][1]) + Hanshu5(m_GNode[2][0], m_GNode[2][1]) + Hanshu5(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu5(m_GNode[4][0], m_GNode[4][1]) + Hanshu5(m_GNode[5][0], m_GNode[5][1]) + Hanshu5(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu5(m_GNode[7][0], m_GNode[7][1]) + Hanshu5(m_GNode[8][0], m_GNode[8][1]) + Hanshu5(m_GNode[9][0], m_GNode[9][1])
		+Hanshu5(m_GNode[10][0], m_GNode[10][1]) + Hanshu5(m_GNode[11][0], m_GNode[11][1]) + Hanshu5(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[5] = w1 * Hanshu6(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu6(m_GNode[1][0], m_GNode[1][1]) + Hanshu6(m_GNode[2][0], m_GNode[2][1]) + Hanshu6(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu6(m_GNode[4][0], m_GNode[4][1]) + Hanshu6(m_GNode[5][0], m_GNode[5][1]) + Hanshu6(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu6(m_GNode[7][0], m_GNode[7][1]) + Hanshu6(m_GNode[8][0], m_GNode[8][1]) + Hanshu6(m_GNode[9][0], m_GNode[9][1])
		+Hanshu6(m_GNode[10][0], m_GNode[10][1]) + Hanshu6(m_GNode[11][0], m_GNode[11][1]) + Hanshu6(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[6] = w1 * Hanshu7(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu7(m_GNode[1][0], m_GNode[1][1]) + Hanshu7(m_GNode[2][0], m_GNode[2][1]) + Hanshu7(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu7(m_GNode[4][0], m_GNode[4][1]) + Hanshu7(m_GNode[5][0], m_GNode[5][1]) + Hanshu7(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu7(m_GNode[7][0], m_GNode[7][1]) + Hanshu7(m_GNode[8][0], m_GNode[8][1]) + Hanshu7(m_GNode[9][0], m_GNode[9][1])
		+Hanshu7(m_GNode[10][0], m_GNode[10][1]) + Hanshu7(m_GNode[11][0], m_GNode[11][1]) + Hanshu7(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[7] = w1 * Hanshu8(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu8(m_GNode[1][0], m_GNode[1][1]) + Hanshu8(m_GNode[2][0], m_GNode[2][1]) + Hanshu8(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu8(m_GNode[4][0], m_GNode[4][1]) + Hanshu8(m_GNode[5][0], m_GNode[5][1]) + Hanshu8(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu8(m_GNode[7][0], m_GNode[7][1]) + Hanshu8(m_GNode[8][0], m_GNode[8][1]) + Hanshu8(m_GNode[9][0], m_GNode[9][1])
		+Hanshu8(m_GNode[10][0], m_GNode[10][1]) + Hanshu8(m_GNode[11][0], m_GNode[11][1]) + Hanshu8(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[8] = w1 * Hanshu9(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu9(m_GNode[1][0], m_GNode[1][1]) + Hanshu9(m_GNode[2][0], m_GNode[2][1]) + Hanshu9(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu9(m_GNode[4][0], m_GNode[4][1]) + Hanshu9(m_GNode[5][0], m_GNode[5][1]) + Hanshu9(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu9(m_GNode[7][0], m_GNode[7][1]) + Hanshu9(m_GNode[8][0], m_GNode[8][1]) + Hanshu9(m_GNode[9][0], m_GNode[9][1])
		+Hanshu9(m_GNode[10][0], m_GNode[10][1]) + Hanshu9(m_GNode[11][0], m_GNode[11][1]) + Hanshu9(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[9] = w1 * Hanshu10(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu10(m_GNode[1][0], m_GNode[1][1]) + Hanshu10(m_GNode[2][0], m_GNode[2][1]) + Hanshu10(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu10(m_GNode[4][0], m_GNode[4][1]) + Hanshu10(m_GNode[5][0], m_GNode[5][1]) + Hanshu10(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu10(m_GNode[7][0], m_GNode[7][1]) + Hanshu10(m_GNode[8][0], m_GNode[8][1]) + Hanshu10(m_GNode[9][0], m_GNode[9][1])
		+Hanshu10(m_GNode[10][0], m_GNode[10][1]) + Hanshu10(m_GNode[11][0], m_GNode[11][1]) + Hanshu10(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[10] = w1 * Hanshu11(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu11(m_GNode[1][0], m_GNode[1][1]) + Hanshu11(m_GNode[2][0], m_GNode[2][1]) + Hanshu11(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu11(m_GNode[4][0], m_GNode[4][1]) + Hanshu11(m_GNode[5][0], m_GNode[5][1]) + Hanshu11(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu11(m_GNode[7][0], m_GNode[7][1]) + Hanshu11(m_GNode[8][0], m_GNode[8][1]) + Hanshu11(m_GNode[9][0], m_GNode[9][1])
		+Hanshu11(m_GNode[10][0], m_GNode[10][1]) + Hanshu11(m_GNode[11][0], m_GNode[11][1]) + Hanshu11(m_GNode[12][0], m_GNode[12][1]));


	Sxishu[11] = w1 * Hanshu12(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu12(m_GNode[1][0], m_GNode[1][1]) + Hanshu12(m_GNode[2][0], m_GNode[2][1]) + Hanshu12(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu12(m_GNode[4][0], m_GNode[4][1]) + Hanshu12(m_GNode[5][0], m_GNode[5][1]) + Hanshu12(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu12(m_GNode[7][0], m_GNode[7][1]) + Hanshu12(m_GNode[8][0], m_GNode[8][1]) + Hanshu12(m_GNode[9][0], m_GNode[9][1])
		+Hanshu12(m_GNode[10][0], m_GNode[10][1]) + Hanshu12(m_GNode[11][0], m_GNode[11][1]) + Hanshu12(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[12] = w1 * Hanshu13(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu13(m_GNode[1][0], m_GNode[1][1]) + Hanshu13(m_GNode[2][0], m_GNode[2][1]) + Hanshu13(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu13(m_GNode[4][0], m_GNode[4][1]) + Hanshu13(m_GNode[5][0], m_GNode[5][1]) + Hanshu13(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu13(m_GNode[7][0], m_GNode[7][1]) + Hanshu13(m_GNode[8][0], m_GNode[8][1]) + Hanshu13(m_GNode[9][0], m_GNode[9][1])
		+Hanshu13(m_GNode[10][0], m_GNode[10][1]) + Hanshu13(m_GNode[11][0], m_GNode[11][1]) + Hanshu13(m_GNode[12][0], m_GNode[12][1]));

	Sxishu[13] = w1 * Hanshu14(m_GNode[0][0], m_GNode[0][1])
		+ w2 * (Hanshu14(m_GNode[1][0], m_GNode[1][1]) + Hanshu14(m_GNode[2][0], m_GNode[2][1]) + Hanshu14(m_GNode[3][0], m_GNode[3][1]))
		+ w3 * (Hanshu14(m_GNode[4][0], m_GNode[4][1]) + Hanshu14(m_GNode[5][0], m_GNode[5][1]) + Hanshu14(m_GNode[6][0], m_GNode[6][1]))
		+ w4 * (Hanshu14(m_GNode[7][0], m_GNode[7][1]) + Hanshu14(m_GNode[8][0], m_GNode[8][1]) + Hanshu14(m_GNode[9][0], m_GNode[9][1])
		+Hanshu14(m_GNode[10][0], m_GNode[10][1]) + Hanshu14(m_GNode[11][0], m_GNode[11][1]) + Hanshu14(m_GNode[12][0], m_GNode[12][1]));
}


void CTriangleInfo::solveCoefficient(int index, double Sxishu[])
{
	double A2[4], B2[2], A3[9], B3[3], A4[16], B4[4], A5[25], B5[5];

	double POLYA[6], POLYB[6], POLYC[6], POLYD[6], POLYE[6], POLYF[6];


	POLYA[0] = 0;
	POLYB[0] = 0;
	POLYC[0] = 0;
	POLYD[0] = 0;
	POLYE[0] = 0;
	POLYF[0] = 1;

	POLYA[1] = 0;
	POLYB[1] = 0;
	POLYC[1] = 0;
	POLYD[1] = 1;
	POLYE[1] = 0;
	POLYF[1] = 0;

	//////////////////////////////////////////////////////////////////////////
	//系数矩阵A2，右端向量B2
	//////////////////////////////////////////////////////////////////////////

	A2[0] = Sxishu[12];
	A2[1] = 1.0;
	A2[2] = Sxishu[9];
	A2[3] = Sxishu[12];

	B2[0] = -Sxishu[13];
	B2[1] = -Sxishu[10];

	Gauss(2, &A2[0], &B2[0]);

	POLYA[2] = 0;
	POLYB[2] = 0;
	POLYC[2] = 0;
	POLYD[2] = B2[0];
	POLYE[2] = 1.0;
	POLYF[2] = B2[1];

	//////////////////////////////////////////////////////////////////////////
	//系数矩阵A3,右端向量B3
	//////////////////////////////////////////////////////////////////////////

	A3[0] = Sxishu[12]; A3[1] = Sxishu[13]; A3[2] = 1.0;
	A3[3] = Sxishu[9]; A3[4] = Sxishu[10]; A3[5] = Sxishu[12];
	A3[6] = POLYD[2]*Sxishu[9] + POLYE[2]*Sxishu[10] + POLYF[2]*Sxishu[12];
	A3[7] = POLYD[2]*Sxishu[10] + POLYE[2]*Sxishu[11] + POLYF[2]*Sxishu[12];
	A3[8] = POLYD[2]*Sxishu[12] + POLYE[2]*Sxishu[13] + POLYF[2];

	B3[0] = -Sxishu[9];
	B3[1] = -Sxishu[5];
	B3[2] = -(POLYD[2] * Sxishu[5] + POLYE[2]*Sxishu[6] + POLYF[2]*Sxishu[9]);

	Gauss(3, &A3[0], &B3[0]);


	POLYA[3] = 1;
	POLYB[3] = 0;
	POLYC[3] = 0;
	POLYD[3] = B3[0];
	POLYE[3] = B3[1];
	POLYF[3] = B3[2];

	//////////////////////////////////////////////////////////////////////////
	//系数矩阵A4,右端向量B4
	//////////////////////////////////////////////////////////////////////////

	A4[0] = Sxishu[9]; A4[1] = Sxishu[12]; A4[2] = Sxishu[13]; A4[3] = 1.0;
	A4[4] = Sxishu[5];  A4[5] = Sxishu[9]; A4[6] = Sxishu[10]; A4[7] = Sxishu[12];

	A4[8] = POLYD[2] * Sxishu[5] + POLYE[2] * Sxishu[6] + POLYF[2] * Sxishu[9];
	A4[9] = POLYD[2] * Sxishu[9] + POLYE[2] * Sxishu[10] + POLYF[2] * Sxishu[12];
	A4[10] = POLYD[2] * Sxishu[10] + POLYE[2] * Sxishu[11] + POLYF[2] * Sxishu[13];
	A4[11] = POLYD[2] * Sxishu[12] + POLYE[2] * Sxishu[13] + POLYF[2];

	A4[12] = POLYA[3] * Sxishu[0] + POLYD[3] * Sxishu[5] + POLYE[3] * Sxishu[6] + POLYF[3] * Sxishu[9];
	A4[13] = POLYA[3] * Sxishu[5] + POLYD[3] * Sxishu[9] + POLYE[3] * Sxishu[10] + POLYF[3] * Sxishu[12];
	A4[14] = POLYA[3] * Sxishu[6] + POLYD[3] * Sxishu[10] + POLYE[3] * Sxishu[11] + POLYF[3] * Sxishu[13];
	A4[15] = POLYA[3] * Sxishu[9] + POLYD[3] * Sxishu[12] + POLYE[3] * Sxishu[13] + POLYF[3];

	B4[0] = -Sxishu[10];
	B4[1] = -Sxishu[6];
	B4[2] = -(POLYD[2] * Sxishu[6] + POLYE[2] * Sxishu[7] + POLYF[2] * Sxishu[10]);
	B4[3] = -(POLYA[3] * Sxishu[1] + POLYD[3] * Sxishu[6] + POLYE[3] * Sxishu[7] + POLYF[3] * Sxishu[10]);


	Gauss(4, &A4[0], &B4[0]);

	POLYA[4] = B4[0];
	POLYB[4] = 1;
	POLYC[4] = 0;
	POLYD[4] = B4[1];
	POLYE[4] = B4[2];
	POLYF[4] = B4[3];

	//////////////////////////////////////////////////////////////////////////
	//系数矩阵A5,右端向量B5
	//////////////////////////////////////////////////////////////////////////

	A5[0] = Sxishu[9]; A5[1] = Sxishu[10]; A5[2] = Sxishu[12]; A5[3] = Sxishu[13]; A5[4] = 1.0;
	A5[5] = Sxishu[5];  A5[6] = Sxishu[6];  A5[7] = Sxishu[9]; A5[8] = Sxishu[10]; A5[9] = Sxishu[12];

	A5[10] = POLYD[2]*Sxishu[5] + POLYE[2]*Sxishu[6] + POLYF[2]*Sxishu[9];
	A5[11] = POLYD[2]*Sxishu[6] + POLYE[2]*Sxishu[7] + POLYF[2]*Sxishu[10];
	A5[12] = POLYD[2]*Sxishu[9] + POLYE[2]*Sxishu[10] + POLYF[2]*Sxishu[12];
	A5[13] = POLYD[2]*Sxishu[10] + POLYE[2]*Sxishu[11] + POLYF[2]*Sxishu[13];
	A5[14]  = POLYD[2]*Sxishu[12] + POLYE[2]*Sxishu[13] + POLYF[2];

	A5[15] = POLYA[3] * Sxishu[0] + POLYD[3] * Sxishu[5] + POLYE[3] * Sxishu[6] + POLYF[3] * Sxishu[9];
	A5[16] = POLYA[3] * Sxishu[1] + POLYD[3] * Sxishu[6] + POLYE[3] * Sxishu[7] + POLYF[3] * Sxishu[10];
	A5[17] = POLYA[3] * Sxishu[5] + POLYD[3] * Sxishu[9] + POLYE[3] * Sxishu[10] + POLYF[3] * Sxishu[12];
	A5[18] = POLYA[3] * Sxishu[6] + POLYD[3] * Sxishu[10] + POLYE[3] * Sxishu[11] + POLYF[3] * Sxishu[13];
	A5[19] = POLYA[3] * Sxishu[9] + POLYD[3] * Sxishu[12] + POLYE[3] * Sxishu[13] + POLYF[3];

	A5[20] = POLYA[4] * Sxishu[0] + POLYB[4] * Sxishu[1] + POLYD[4] * Sxishu[5] + POLYE[4] * Sxishu[6] + POLYF[4] * Sxishu[9];
	A5[21] = POLYA[4] * Sxishu[1] + POLYB[4] * Sxishu[2] + POLYD[4] * Sxishu[6] + POLYE[4] * Sxishu[7] + POLYF[4] * Sxishu[10];
	A5[22] = POLYA[4] * Sxishu[5] + POLYB[4] * Sxishu[6] + POLYD[4] * Sxishu[9] + POLYE[4] * Sxishu[10] + POLYF[4] * Sxishu[12];
	A5[23] = POLYA[4] * Sxishu[6] + POLYB[4] * Sxishu[7] + POLYD[4] * Sxishu[10] + POLYE[4] * Sxishu[11] + POLYF[4] * Sxishu[13];
	A5[24] = POLYA[4] * Sxishu[9] + POLYB[4] * Sxishu[10] + POLYD[4] * Sxishu[12] + POLYE[4] * Sxishu[13] + POLYF[4];

	B5[0] = -Sxishu[11];
	B5[1] = -Sxishu[7];
	B5[2] = -(POLYD[2] * Sxishu[7] + POLYE[2] * Sxishu[8] + POLYF[2] * Sxishu[11]);
	B5[3] = -(POLYA[3] * Sxishu[2] + POLYD[3] * Sxishu[7] + POLYE[3] * Sxishu[8] + POLYF[3] * Sxishu[11]);
	B5[4] = -(POLYA[4] * Sxishu[2] + POLYB[4] * Sxishu[3] + POLYD[4] * Sxishu[7] + POLYE[4] * Sxishu[8] + POLYF[4] * Sxishu[11]);

	Gauss(5, &A5[0], &B5[0]);

	POLYA[5] = B5[0];
	POLYB[5] = B5[1];
	POLYC[5] = 1;
	POLYD[5] = B5[2];
	POLYE[5] = B5[3];
	POLYF[5] = B5[4];


	
	int gindex = index*BASIS_FUNCTIONS*BASIS_FUNCTION_COEFFS;

	memset(&basis_fun_coeff[gindex], 0, sizeof(double)*BASIS_FUNCTIONS*BASIS_FUNCTION_COEFFS);

	basis_fun_coeff[gindex] = 1;

	basis_fun_coeff[gindex+BASIS_FUNCTION_COEFFS+1] = 1/sqrt(area[index]);

	basis_fun_coeff[gindex+2*BASIS_FUNCTION_COEFFS] = POLYF[2];
	basis_fun_coeff[gindex+2*BASIS_FUNCTION_COEFFS+1] = POLYD[2] /sqrt(area[index]);
	basis_fun_coeff[gindex+2*BASIS_FUNCTION_COEFFS+2] = 1/sqrt(area[index]);
	
	basis_fun_coeff[gindex+3*BASIS_FUNCTION_COEFFS] = POLYF[3];
	basis_fun_coeff[gindex+3*BASIS_FUNCTION_COEFFS+1] = POLYD[3] /sqrt(area[index]);
	basis_fun_coeff[gindex+3*BASIS_FUNCTION_COEFFS+2] = POLYE[3] /sqrt(area[index]);
	basis_fun_coeff[gindex+3*BASIS_FUNCTION_COEFFS+3] = 1 / area[index];

	basis_fun_coeff[gindex+4*BASIS_FUNCTION_COEFFS] = POLYF[4] ;
	basis_fun_coeff[gindex+4*BASIS_FUNCTION_COEFFS+1] = POLYD[4] /sqrt(area[index]);
	basis_fun_coeff[gindex+4*BASIS_FUNCTION_COEFFS+2] = POLYE[4] /sqrt(area[index]);
	basis_fun_coeff[gindex+4*BASIS_FUNCTION_COEFFS+3] = POLYA[4] / area[index];
	basis_fun_coeff[gindex+4*BASIS_FUNCTION_COEFFS+4] = 1 / area[index];

	basis_fun_coeff[gindex+5*BASIS_FUNCTION_COEFFS] = POLYF[5] ;
	basis_fun_coeff[gindex+5*BASIS_FUNCTION_COEFFS+1] = POLYD[5] /sqrt(area[index]);
	basis_fun_coeff[gindex+5*BASIS_FUNCTION_COEFFS+2] = POLYE[5] /sqrt(area[index]);
	basis_fun_coeff[gindex+5*BASIS_FUNCTION_COEFFS+3] = POLYA[5] / area[index];
	basis_fun_coeff[gindex+5*BASIS_FUNCTION_COEFFS+4] = POLYB[5] / area[index];
	basis_fun_coeff[gindex+5*BASIS_FUNCTION_COEFFS+5] = 1 / area[index];
}


void CTriangleInfo::initMassMatrix(int index)
{
	assert(index>=0 && index<_cell_num);

	int gindex = index*VOLUME_GPOINTS;
	double x, y;

	for(int i = 0; i < BASIS_FUNCTIONS; ++i)
	{
		mass_coeff[index+i*_cell_num] = 0;
		for ( int j=0; j<VOLUME_GPOINTS; ++j)
		{
			x = vol_gauss_vertice[gindex+j].getX();
			y = vol_gauss_vertice[gindex+j].getY();

			mass_coeff[index+i*_cell_num] += vol_gauss_weight[index+j*_cell_num]
								* basisFunction(index, i, x, y)* basisFunction(index, i, x, y);
		}
	}
}


void CTriangleInfo::allocateMemory(int num)
{
	assert(num>0);
	_cell_num = num;
	// 为三角形信息类分配空间
	barycenter.resize(num);

	edge_middle_vertice.resize(num*TRIANGLE_EDGES*EDGE_GPOINTS);

	radius_of_inscribed_circle.resize(num);

	vol_gauss_vertice.resize(VOLUME_GPOINTS*num);

	edge_gauss_vertice.resize(EDGE_GPOINTS*TRIANGLE_EDGES*num);

	basis_fun_coeff.resize(BASIS_FUNCTION_COEFFS*BASIS_FUNCTIONS*num);


	area = new double[num];

	perimeter = new double[num];

	outer_normal_vector = new double[num*TRIANGLE_EDGES*2];

	mass_coeff = new double[num*BASIS_FUNCTIONS];

	vol_bf_value = new double[num*VOLUME_GPOINTS*BASIS_FUNCTIONS];

	vol_bdf_value = new double[num*2*VOLUME_GPOINTS*BASIS_FUNCTIONS];

	edge_bf_value = new double[num*TRIANGLE_EDGES*EDGE_GPOINTS*BASIS_FUNCTIONS];

	vol_gauss_weight = new double[num*VOLUME_GPOINTS];

	edge_gauss_weight = new double[num*EDGE_GPOINTS*TRIANGLE_EDGES];

}

CTriangleInfo::~CTriangleInfo()
{
	delete []area;				
	delete []perimeter;			
	delete []outer_normal_vector;
	delete []mass_coeff;			
	delete []vol_bf_value;		
	delete []vol_bdf_value;		
    delete []edge_bf_value; 		
    delete []vol_gauss_weight;	
    delete []edge_gauss_weight;	
}