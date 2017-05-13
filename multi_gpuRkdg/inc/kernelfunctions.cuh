/**
 * 定义CUDA的核函数
 * 由于CUDA的核函数必须为全局函数，从而需要在类的外部单独定义
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */
 
#pragma once 
#include "cppstdheaders.h"
#include "cuda_runtime.h"
#include "defines.h"
using namespace std;

extern "C"{

	/**
	 * 在GPU上求解时间步长
	 * @param[in] tnum integer 实际三角形数目
	 * @param[in] gamma double 气体比热比
	 * @param[out] ddt double 时间步长
	 * @param[in] freedom_rho double[] rho自由度
	 * @param[in] freedom_rhou double[] rhou自由度
	 * @param[in] freedom_rhov double[] rhov自由度
	 * @param[in] freedom_rhoE double[] rhoE自由度
	 * @param[in] perimeter double[] 三角形周长
	 * @param[in] area double[] 三角形面积
	 */
	__global__ void kernel_getTimeStep(int tnum, double gamma, double cfl, double * ddt,

				const double *freedom_rho, const double *freedom_rhou, 
				const double *freedom_rhov, const double *freedom_rhoE, 

				const double *perimeter, const double *area
				);

	/**
	 * 计算守恒量
	 * @param[in] tnum integer 实际三角形数目
	 * @param[in] double_pitch integer 对齐后的浮点型数组宽度
	 * @param[in] freedom_rho double[] rho自由度
	 * @param[in] freedom_rhou double[] rhou自由度
	 * @param[in] freedom_rhov double[] rhov自由度
	 * @param[in] freedom_rhoE double[] rhoE自由度
	 * @param[in] convar_rho_vol double[] 单元体上守恒量rho的值
	 * @param[in] convar_rhou_vol double[] 单元体上守恒量rhou的值
	 * @param[in] convar_rhov_vol double[] 单元体上守恒量rhov的值
	 * @param[in] convar_rhoE_vol double[] 单元体上守恒量rhoE的值
	 * @param[in] convar_rho_edge double[] 单元边上守恒量rho的值
	 * @param[in] convar_rhou_edge double[] 单元边上守恒量rhou的值
	 * @param[in] convar_rhov_edge double[] 单元边上守恒量rhov的值
	 * @param[in] convar_rhoE_edge double[] 单元边上守恒量rhoE的值
	 * @param[in] vol_bf_value double[] 基函数在单元体内的值
	 * @param[in] edge_bf_value double[] 基函数在单元边上的值
	 */
	__global__ void kernel_calculateConVars(int tnum, int double_pitch,
		const double *freedom_rho,  const double *freedom_rhou,
		const double *freedom_rhov, const double *freedom_rhoE,

		double *convar_rho_vol,   double *convar_rhou_vol,
		double *convar_rhov_vol,  double *convar_rhoE_vol,

		double *convar_rho_edge,  double *convar_rhou_edge, 
		double *convar_rhov_edge, double *convar_rhoE_edge, 

		const double *vol_bf_value, const double *edge_bf_value
		);

	/** 
	 * 处理边界条件
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] num integer  网格中单元总数
	 * @param[in] double_pitch integer  对齐后的双精度数组宽度
	 * @param[in] rho double 远场密度
	 * @param[in] rhou double 远场守恒量rhou
	 * @param[in] rhov double 远场守恒量rhov
	 * @param[in] rhoE double 远场守恒量rhoE
	 * @param[in|out] freedom_rho double[] rho自由度
	 * @param[in|out] freedom_rhou double[] rhou自由度
	 * @param[in|out] freedom_rhov double[] rhov自由度
	 * @param[in|out] freedom_rhoE double[] rhoE自由度
	 * @param[in|out] convar_rho_edge double[] 边界上守恒量rho的值
	 * @param[in|out] convar_rhou_edge double[] 边界上守恒量rhou的值
	 * @param[in|out] convar_rhov_edge double[] 边界上守恒量rhov的值
	 * @param[in|out] convar_rhoE_edge double[] 边界上守恒量rhoE的值
	 * @param[in] neighbour int[] 单元邻居信息
	 * @param[in] sharedEdge int[] 单元共享边信息
	 * @param[in] triangle_flag int[] 单元标记
	 * @param[in] outer_normal double[] 单元边外法向量
	 */
	__global__ void kernel_boundaryCondition(int tnum, int num, int double_pitch, 
		double rho, double rhou, double rhov, double rhoE,

		double *convar_rho_edge,  double *convar_rhou_edge, 
		double *convar_rhov_edge, double *convar_rhoE_edge,

		double *freedom_rho,  double *freedom_rhou, 
		double *freedom_rhov, double *freedom_rhoE,

		const int *neighbour, const int *sharedEdge, 

		const int *triangle_flag, const double *outer_normal
		);

	/** 
	 * 计算体积分守恒残差
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] double_pitch integer  对齐后的双精度浮点型数组宽度
	 * @param[in] gamma double 气体比热比
	 * @param[in] convar_rho_vol double[] 单元内守恒量rho的值
	 * @param[in] convar_rhou_vol double[] 单元内守恒量rhou的值
	 * @param[in] convar_rhov_vol double[] 单元内守恒量rhov的值
	 * @param[in] convar_rhoE_vol double[] 单元内守恒量rhoE的值
	 * @param[out] rhs_volume_rho double[] 守恒量rho在单元体积分残差
	 * @param[out] rhs_volume_rhou double[] 守恒量rhou在单元体积分残差
	 * @param[out] rhs_volume_rhov double[] 守恒量rhov在单元体积分残差
	 * @param[out] rhs_volume_rhoE double[] 守恒量rhoE在单元体积分残差
	 * @param[in] vol_gauss_weight double[] 体积分高斯权重
	 * @param[in] vol_bf_value double[] 基函数导数在单元体上的值
	 */
	__global__ void kernel_calculateVolumeRHS(int tnum, int double_pitch, double gamma,

		const double *convar_rho_vol,  const double *convar_rhou_vol, 
		const double *convar_rhov_vol, const double *convar_rhoE_vol,

		double *rhs_volume_rho,  double *rhs_volume_rhou, 
		double *rhs_volume_rhov, double *rhs_volume_rhoE,

		const double *vol_gauss_weight, const double *vol_bdf_value
		);

	/** 
	 * 计算LF通量系数
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] int_pitch integer  对齐后的整型数组宽度
	 * @param[in] double_pitch integer  对齐后的双精度数组宽度
	 * @param[in] outer_normal_vector double [] 单元边的外法向量
	 * @param[in] neighbour int[] 单元共享边信息
	 * @param[in] freedom_rho double[] 守恒量rho自由度
	 * @param[in] freedom_rhou double[] 守恒量rhou自由度
	 * @param[in] freedom_rhov double[] 守恒量rhov自由度
	 * @param[in] freedom_rhoE double[] 守恒量rhoE自由度
	 * @param[out] lfflux_coeff double[] LF通量系数
	 */
	__global__ void kernel_calculateLFCoeff(
		int tnum, int int_pitch, int double_pitch, double gamma, 

		const double *outer_normal_vector,  const int *neighbour, 

		const double* freedom_rho, const double *freedom_rhou, 
		const double *freedom_rhov, const double *freedom_rhoE, 

		double *lfflux_coeff
		);

	/** 
	 * 计算边上FG的值
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] num integer  网格总数
	 * @param[in] double_pitch integer  对齐后的数组宽度
	 * @param[in] gamma double 气体比热比
	 * @param[in] convar_rho_edge double[] 单元边上守恒量rho的值
	 * @param[in] convar_rhou_edge double[] 单元边上守恒量rhou的值
	 * @param[in] convar_rhov_edge double[] 单元边上守恒量rhov的值
	 * @param[in] convar_rhoE_edge double[] 单元边上守恒量rhoE的值
	 * @param[out] fedge_rho double[] 单元边上f通量
	 * @param[out] fedge_rhou double[] 单元边上f通量
	 * @param[out] fedge_rhov double[] 单元边上f通量
	 * @param[out] fedge_rhoE double[] 单元边上f通量
	 * @param[out] gedge_rho double[] 单元边上f通量
	 * @param[out] gedge_rhou double[] 单元边上f通量
	 * @param[out] gedge_rhov double[] 单元边上f通量
	 * @param[out] gedge_rhoE double[] 单元边上f通量
	 */
	__global__ void kernel_calculateEdgeFG(
		 int tnum, int num, int double_pitch, double gamma, 

		 const double *convar_rho_edge, const double *convar_rhou_edge,
		 const double *convar_rhov_edge, const double *convar_rhoE_edge,

		 double *fedge_rho, double *fedge_rhou,
		 double *fedge_rhov, double *fedge_rhoE,

		 double *gedge_rho, double *gedge_rhou,
		 double *gedge_rhov, double *gedge_rhoE
		);


	/** 
	 * 计算边上通量的值
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] int_pitch integer  对齐后整型数组宽度
	 * @param[in] double_pitch integer 对齐后双精度数组宽度
	 * @param[in] neighbour int[] 单元邻居信息
	 * @param[in] sharedEdge int[] 单元共享边信息
	 * @param[in] convar_rho_edge double[] 单元边上守恒量rho的值
	 * @param[in] convar_rhou_edge double[] 单元边上守恒量rhou的值
	 * @param[in] convar_rhov_edge double[] 单元边上守恒量rhov的值
	 * @param[in] convar_rhoE_edge double[] 单元边上守恒量rhoE的值
	 * @param[in] fedge_rho double[] 单元边上f通量
	 * @param[in] fedge_rhou double[] 单元边上f通量
	 * @param[in] fedge_rhov double[] 单元边上f通量
	 * @param[in] fedge_rhoE double[] 单元边上f通量
	 * @param[in] gedge_rho double[] 单元边上f通量
	 * @param[in] gedge_rhou double[] 单元边上f通量
	 * @param[in] gedge_rhov double[] 单元边上f通量
	 * @param[in] gedge_rhoE double[] 单元边上f通量
	 * @param[in] outer_normal_vector double[] 单元边外法向量
	 * @param[in] lfflux_coeff double[] 单元边上的LF通量系数
	 * @param[out] lfflux_rho double[] 单元边上各个点处的通量
	 * @param[out] lfflux_rhou double[] 单元边上各个点处的通量
	 * @param[out] lfflux_rhov double[] 单元边上各个点处的通量
	 * @param[out] lfflux_rhoE double[] 单元边上各个点处的通量
	 */
/*	__global__ void kernel_calculateFlux(
		int tnum, int int_pitch, int double_pitch, 

		const int *neighbour, const int *sharedEdge,

		const double *convar_rho_edge, const double *convar_rhou_edge,
		const double *convar_rhov_edge, const double *convar_rhoE_edge,
		
		const double *fedge_rho, const double *fedge_rhou,
		const double *fedge_rhov, const double *fedge_RhoE,
		
		const double *gedge_rho, const double *gedge_rhou,
		const double *gedge_rhov, const double *gedge_rhoE,

		const double *outer_normal_vector, const double *lfflux_coeff,

		double *lfflux_rho, double *lfflux_rhou,
		double *lfflux_rhov, double *lfflux_rhoE
		);
*/
	__global__ void kernel_calculateFlux(
		int tnum, int int_pitch, int double_pitch, 

		const int *neighbour, const int *sharedEdge,

		const double *convar_rho_edge, const double *convar_rhou_edge,
//		const double *convar_rhov_edge, const double *convar_rhoE_edge,
		
		const double *fedge_rho, const double *fedge_rhou,
//		const double *fedge_rhov, const double *fedge_RhoE,
		
		const double *gedge_rho, const double *gedge_rhou,
//		const double *gedge_rhov, const double *gedge_rhoE,

		const double *outer_normal_vector, const double *lfflux_coeff,

		double *lfflux_rho, double *lfflux_rhou
//		double *lfflux_rhov, double *lfflux_rhoE
		);
	
	/** 
	 * 计算单元各条边上的积分残差
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] double_pitch integer  对齐后的双精度数组宽度
	 * @param[in] edge_gauss_weight double[] 单元边的积分权重
	 * @param[in] edge_bf_value double[] 基函数在边上的值
	 * @param[out] rhs_edge_rho double[] 单元边上的积分残差
	 * @param[out] rhs_edge_rhou double[] 单元边上的积分残差
	 * @param[out] rhs_edge_rhov double[] 单元边上的积分残差
	 * @param[out] rhs_edge_rhoE double[] 单元边上的积分残差
	 */
	__global__ void kernel_calculateEdgeRHS(
		int tnum, int double_pitch, 

		const double *edge_gauss_weight, const double *edge_bf_value,

		const double *lfflux_rho,  const double *lfflux_rhou, 
		const double *lfflux_rhov, const double *lfflux_rhoE, 

		double *rhs_edge_rho,  double *rhs_edge_rhou, 
		double *rhs_edge_rhov, double *rhs_edge_rhoE,

		const double *rhs_volume_rho,  const double *rhs_volume_rhou, 
		const double *rhs_volume_rhov, const double *rhs_volume_rhoE
		);

	/** 
	 * rkdg时间推进第一步
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] double_pitch integer  对齐后的双精度数组宽度
	 * @param[in] dt  double 时间步长
	 * @param[in] mass_coeff double[] 基函数质量系数
	 * @param[in|out] freedom_rho double[] 守恒量rho自由度
	 * @param[in|out] freedom_rhou double[] 守恒量rhou自由度
	 * @param[in|out] freedom_rhov double[] 守恒量rhov自由度
	 * @param[in|out] freedom_rhoE double[] 守恒量rhoE自由度
	 * @param[in] rhs_edge_rho double[] 单元边上的积分残差
	 * @param[in] rhs_edge_rhou double[] 单元边上的积分残差
	 * @param[in] rhs_edge_rhov double[] 单元边上的积分残差
	 * @param[in] rhs_edge_rhoE double[] 单元边上的积分残差
	 */
	__global__ void kernel_rkdgStepOne(
		int tnum, int double_pitch, double dt, const double *mass_coeff,

		double *freedom_rho,  double *freedom_rhou, 
		double *freedom_rhov, double *freedom_rhoE, 

		const double *rhs_edge_rho,  const double *rhs_edge_rhou,
		const double *rhs_edge_rhov, const double *rhs_edge_rhoE
		);

	/** 
	 * rkdg时间推进第二步
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] double_pitch integer  对齐后的双精度数组宽度
	 * @param[in] dt  double 时间步长
	 * @param[in] mass_coeff double[] 基函数质量系数
	 * @param[out] freedom_rho double[] 守恒量rho自由度
	 * @param[out] freedom_rhou double[] 守恒量rhou自由度
	 * @param[out] freedom_rhov double[] 守恒量rhov自由度
	 * @param[out] freedom_rhoE double[] 守恒量rhoE自由度
	 * @param[in] rhs_edge_rho double[] 单元边上的积分残差
	 * @param[in] rhs_edge_rhou double[] 单元边上的积分残差
	 * @param[in] rhs_edge_rhov double[] 单元边上的积分残差
	 * @param[in] rhs_edge_rhoE double[] 单元边上的积分残差
	 * @param[in] freedom_rho_old  double[] rho旧自由度
	 * @param[in] freedom_rhou_old double[] rhou旧自由度
	 * @param[in] freedom_rhov_old double[] rhov旧自由度
	 * @param[in] freedom_rhoE_old double[] rhoE旧自由度
	 */
	__global__ void kernel_rkdgStepTwo(
		int tnum, int double_pitch, double dt, const double *mass_coeff,

		double *freedom_rho,  double *freedom_rhou, 
		double *freedom_rhov, double *freedom_rhoE, 

		const double *rhs_edge_rho,  const double *rhs_edge_rhou,
		const double *rhs_edge_rhov, const double *rhs_edge_rhoE,

		const double *freedom_rho_old,  const double *freedom_rhou_old, 
		const double *freedom_rhov_old, const double *freedom_rhoE_old
		);

	/** 
	 * rkdg时间推进第三步
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] double_pitch integer  对齐后的双精度数组宽度
	 * @param[in] dt  double 时间步长
	 * @param[in] mass_coeff double[] 基函数质量系数
	 * @param[in|out] freedom_rho double[] 守恒量rho自由度
	 * @param[in|out] freedom_rhou double[] 守恒量rhou自由度
	 * @param[in|out] freedom_rhov double[] 守恒量rhov自由度
	 * @param[in|out] freedom_rhoE double[] 守恒量rhoE自由度
	 * @param[in] rhs_edge_rho double[] 单元边上的积分残差
	 * @param[in] rhs_edge_rhou double[] 单元边上的积分残差
	 * @param[in] rhs_edge_rhov double[] 单元边上的积分残差
	 * @param[in] rhs_edge_rhoE double[] 单元边上的积分残差
	 * @param[in] freedom_rho_old double[] rho旧自由度
	 * @param[in] freedom_rhou_old double[] rhou旧自由度
	 * @param[in] freedom_rhov_old double[] rhov旧自由度
	 * @param[in] freedom_rhoE_old double[] rhoE旧自由度
	 */
	__global__ void kernel_rkdgStepThree(
		int tnum, int double_pitch, double dt, const double *mass_coeff,

		double *freedom_rho,  double *freedom_rhou, 
		double *freedom_rhov, double *freedom_rhoE, 

		const double *rhs_edge_rho,    const double *rhs_edge_rhou,
		const double *rhs_edge_rhov,   const double *rhs_edge_rhoE,

		const double *freedom_rho_old,     const double *freedom_rhou_old, 
		const double *freedom_rhov_old,    const double *freedom_rhoE_old
		);

	/** 
	 * 计算残差
	 * @param[in] tnum integer 计算区域中真实三角形数目
	 * @param[in] freedom_rho double[] 守恒量rho的自由度
	 * @param[in] freedom_rhoE double[] 守恒量rhoE的自由度
	 * @param[in] freedom_rho_old  double[] 旧自由度
	 * @param[in] freedom_rhoE_old double[] 旧自由度
	 * @param[out] residual double[] 残差
	 */
	__global__ void kernel_calculateResidual(
		int tnum, 

		const double *freedom_rho,     const double *freedom_rhoE,
		const double *freedom_rho_old, const double *freedom_rhoE_old,

		double *residual
		);
};
