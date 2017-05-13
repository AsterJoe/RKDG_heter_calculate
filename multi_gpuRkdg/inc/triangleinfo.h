/**
 * 三角形信息类，该类保存了除三角形组成之外的额外信息。
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */
 
#pragma once

#include "cppstdheaders.h"
#include "gauss.h"
#include "vertice2d.h"
#include "defines.h"
using namespace std;

/**
 * 三角形信息类，该类保存除顶点、边之外的额外信息
 * 重心，周长，面积，外法向量等信息将都保存在该类中
 */
class CTriangleInfo{
	// 属性
	public:
		
		vector<CVertice2D> barycenter;				/**< 三角形重心 */
		
		vector<CVertice2D> edge_middle_vertice;		/**< 边中点 */
		
		vector<double> radius_of_inscribed_circle;	/**< 内切圆半径 */
		
		
		vector<CVertice2D> vol_gauss_vertice; 		/**< 体积分的高斯系数 */
		vector<CVertice2D> edge_gauss_vertice;		/**< 边上的高斯积分点 */
		
		vector<double> basis_fun_coeff;				/**< 基函数系数 */
	
	protected:
		int _cell_num;								/**< 总单元数目 */


	// 以下数据需要被复制到GPU上	
	public:
		double *area;								/**< 面积 */
		double *perimeter;							/**< 周长 */
		double *outer_normal_vector;				/**< 单元边的外法向量 */
		
		double *mass_coeff;							/**< 基函数的质量系数 */
		
		double *vol_bf_value;						/**< 基函数在体高斯积分节点上的值 */
		
		double *vol_bdf_value;						/**< 基函数的导数在体高斯积分节点的值 */
		
		double *edge_bf_value; 						/**< 基函数在边上高斯积分节点的值 */
	
		double *vol_gauss_weight;					/**< 高斯体积分的权重 */

		double *edge_gauss_weight;					/**< 高斯边积分的权重 */
		
	// 行为
	protected:
		/** 
		 * 初始化系数矩阵，本函数来自程剑
		 * @param[in] index integer 三角形编号
		 * @param[in] x double[] 三角形顶点横坐标
		 * @param[in] y double[] 三角形顶点横坐标
		 * @param[in|out] s double[] 未知
		 */
		void initXishuMatrix(int index, double x[], double y[], double s[]);


		/** 
		 * 求解基函数系数， 该函数来自程剑
		 * @param[in] index integer 三角形编号
		 * @param[in|out] s double[] 未知
		 */
		void solveCoefficient(int index, double s[]);

		/** 
		 * 初始化系数矩阵
		 * @param[in] index integer 三角形编号
		 */
		void initMassMatrix(int index);
	
	public:
		/** 构造函数 */
		CTriangleInfo();

		// getter
		/** 
		 * 获取指定三角形单元的重心
		 * @param[in] index integer 三角形编号
		 * @return CVertice2D 三角形重心
 		 */
		CVertice2D& getBarycenter(int index);

		/**
		 * 返回边中点
		 * @param[in] tindex integer 三角形编号
		 * @param[in] eindex integer 第几条边
		 * @return CVertice2D 该边中点
		 */
		CVertice2D& getEdgeMiddleVertice(int tindex, int eindex);

		/**
		 * @param[in] index integer 三角形编号
		 * @return double 三角形外切圆半径
		 */
		double getRadiusOfInscribedCircle(int index);
		
		/**
		 * 基函数
		 * @param[in] tindex integer 三角形编号
		 * @param[in] findex integer 第几个基函数
		 * @param[in] x	double 横坐标
		 * @param[in] y double 纵坐标
		 * @return double 基函数在该点的值
		 */
		double basisFunction(int tindex, int findex, double x, double y);
		
		/**
		 * 基函数的导数在指定点的值
		 * @param[in] tindex integer 三角形编号
		 * @param[in] findex integer 第几个函数（基函数的导数）
		 * @param[in] x	double 横坐标
		 * @param[in] y double 纵坐标
		 * @return double 基函数在该点的值
		 */
		double derivativeFunction(int tindex, int findex, double x, double y);
		

		/**
		 * 设置三角形内切圆半径
		 * @param[in] tindex integer 三角编号
		 * @param[in] radius double 三角形内切圆半径
		 */
		void setRadiusOfInscribedCircle(int tindex, double radius);


		/** 
		 * 根据顶点坐标初始化网格信息
		 * @param[in] index integer 三角形编号
		 * @param[in] x double[] 三个顶点的横坐标
		 * @param[in] y double[] 三个顶点的纵坐标
		 */
		void initInformation(int index, double x[3], double y[3]);

		
		/**
		 * 数组分配内存
		 * @param[in] num integer 网格总单元数目
		 */
		void allocateMemory(int num);

		/** 析构函数 */
		~CTriangleInfo();
	
};