/**
 * 定义需要在GPU上使用的数组
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */
 
#pragma once

#include "cppstdheaders.h"
#include "cuda_runtime.h"
#include "defines.h"
#include "myexception.h"
using namespace std;

/**
 * 定义需要在GPU上使用的数组，这些数组将在GPU计算的时候使用
 */
class CCUDAArrays {
	// 属性
	public:
		// 单元信息
		int *neighbour;							/**< 单元邻居信息 */
		int *sharedEdge;						/**< 单元共享边信息 */
		int *triangle_flag;						/**< 单元标记 */

		double *area;							/**< 网格单元面积 */
		double *perimeter;						/**< 周长 */
		double *outer_normal_vector;			/**< 单元边的外法向量 */
		double *mass_coeff;						/**< 基函数的质量系数 */
		double *vol_bf_value;					/**< 基函数在体高斯积分节点上的值 */
		double *vol_bdf_value;					/**< 基函数的导数在体高斯积分节点的值 */
		double *edge_bf_value; 					/**< 基函数在边上高斯积分节点的值 */
		double *vol_gauss_weight;				/**< 高斯体积分的权重 */
		double *edge_gauss_weight;				/**< 高斯边积分的权重 */

		// 求解过程量
		double *freedom_rho;			/**< rho自由度 */
		double *freedom_rhou;			/**< rhou自由度 */
		double *freedom_rhov;			/**< rhov自由度 */
		double *freedom_rhoE;			/**< rhoE自由度 */

		double *convar_rho_vol;			/**< 守恒量rho在面高斯积分节点的值 */
		double *convar_rhou_vol;		/**< 守恒量rhou在面高斯积分节点的值 */
		double *convar_rhov_vol;		/**< 守恒量rhov在面高斯积分节点的值 */
		double *convar_rhoE_vol;		/**< 守恒量rhoE在面高斯积分节点的值 */
		double *convar_rho_edge;		/**< 守恒量rho在边高斯积分节点的值 */
		double *convar_rhou_edge;		/**< 守恒量rhou在边高斯积分节点的值 */
		double *convar_rhov_edge;		/**< 守恒量rhov在边高斯积分节点的值 */
		double *convar_rhoE_edge;		/**< 守恒量rhoE在边高斯积分节点的值 */

		double *rhs_volume_rho;			/**< rho的体积分残差 */
		double *rhs_volume_rhou;		/**< rhou的体积分残差 */
		double *rhs_volume_rhov;		/**< rhov的体积分残差 */
		double *rhs_volume_rhoE;		/**< rhoE的体积分残差 */

		double *rhs_edge_rho;			/**< rho的线积分残差 */
		double *rhs_edge_rhou;			/**< rhou的线积分残差 */
		double *rhs_edge_rhov;			/**< rhov的线积分残差 */
		double *rhs_edge_rhoE;			/**< rhoe的线积分残差 */

		double *lfflux_coeff;			/**< LF通量系数 */

		double *fedge_rho;				/**< f在边上的值rho */
		double *fedge_rhou;				/**< f在边上的值rhou */
		double *fedge_rhov;				/**< f在边上的值rhov */
		double *fedge_rhoE;				/**< f在边上的值rhoE */
               
		double *gedge_rho;				/**< g在边上的值rho */
		double *gedge_rhou;				/**< g在边上的值rhou */
		double *gedge_rhov;				/**< g在边上的值rhov */
		double *gedge_rhoE;				/**< g在边上的值rhoE */

		double *lfflux_rho;				/**< lf通量rho */
		double *lfflux_rhou;			/**< lf通量rhou */
		double *lfflux_rhov;			/**< lf通量rhov */
		double *lfflux_rhoE;			/**< lf通量rhoE */

		// 旧自由度
		double *freedom_rho_old;			/**< 旧自由度 */
		double *freedom_rhou_old;			/**< 旧自由度 */
		double *freedom_rhov_old;			/**< 旧自由度 */
		double *freedom_rhoE_old;			/**< 旧自由度 */

		double *ddt;					/**< GPU上存放时间步长的变量 */
		double *residual;				/**< 残差 */

		int *elem_index;

	protected:
		size_t _int_pitch;				/**< 整型数组的填充 */
		size_t _double_pitch;			/**< 浮点型数组的填充 */

	// 方法
	public:
		/** 构造函数 */
		CCUDAArrays();
		
		/**
		 * 在GPU上分配内存
		 * @param[in] num integer 网格单元总数
		 */
		void allocateMemory(int num);

		/**
		 * 将GPU上的数组初始化
		 * @param[in] num integer 网格单元总数
		 */
		//void memsetArrays(int num);

		/** @return 整型数组的填充长度 */
		size_t getIntPitch(void) const;

		/** @return 双精度数组的填充长度 */
		size_t getDoublePitch(void) const;
		
		/** 析构函数 */
		~CCUDAArrays();
};