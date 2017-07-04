/**
 * CUDA实现的二维非结构网格上RKDG算法解法器
 * 解法器负责程序的所有流程：读取配置文件，调用网格接口读入和输出网格，进行RKDG方向时间推进，以及输出结果
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */
 
#pragma once

#include "cppstdheaders.h"
#include "config.h"
#include "unstructuredgrid.h"
#include "cudaarrays.h"
#include "kernelfunctions.cuh"
#include "defines.h"
#include "myexception.h"
#include "mytime.h"
using namespace std;

/**
 * 二维非结构网格上基于CUDA的RKDG解法器
 * 解法器负责程序的所有流程：读取配置文件，调用网格接口读入和输出网格，进行RKDG方向时间推进，输出结果
 * 类中成员函数有可能抛出标准异常，请在调用的时候捕捉。
 * @example 该类用法如下
 *  CCUDARkdgSolver solver;
 *  solver.config_file = config_filename
 *  solver.run();
 *  其中config_filename是程序的配置文件名，在该配置文件中请按照格式填写程序配置
 */
class CCUDARkdgSolver {
	// 属性
	public:
		int myid;
		int nprocs;
		int num_commu;           /**< 信息交换的次数>*/

		string title;				/**< 程序名称 */

		double alpha;				/**< 翼型攻角 */
		double gamma;				/**< 气体比热比 */
		double mach;				/**< 气体马赫数 */
		double cfl;					/**< cfl数 */
		double rhoref;				/**< 初始参考密度 */
		double pref;				/**< 初始参考压强 */
		
		char log_history;			/**< 是否记录收敛历史 */
		int print_interval;			/**< 屏幕信息输出间隔 */
	
		string config_file;			/**< 程序配置文件 */
		string gridconf;			/**< 网格配置输入文件 */
		string solution_file;		/**< 问题的解输出文件 */
		string log_file;			/**< 程序日志文件 */
		string residual_file;		/**< 残差记录文件 */
		
		CUnstructuredGrid grid;		/**< 网格接口 */
		
		int threads_per_block;		/**< 每个CUDA块中线程的数目 */
		int reduction_threads;		/**< 在聚合函数里（通常只有一个线程块）的线程数目 */
		CCUDAArrays _cuarrays;		/**< GPU上求解使用的数组 */
//	protected:
		double _terminal_time;		/**< 求解总时间 */

		double *_freedom_rho;		/**< rho自由度 */
		double *_freedom_rhou;		/**< rhou自由度 */
		double *_freedom_rhov;		/**< rhov自由度 */
		double *_freedom_rhoE;		/**< rhoE自由度 */

		double *_dt;				/**< 时间步长 */
		double *_residual;			/**< 残差 */
		//CCUDAArrays _cuarrays;		/**< GPU上求解使用的数组 */
		
		
	private:

	// 方法
	public:
		/** 默认构造函数 */
		CCUDARkdgSolver();
		
		/** 初始化程序配置 */
		void initConfig(void);

		/** 运行程序求解流场 */
		void run(int myid, int nprocs);
		
		/** 检测CUDA设备 */
		void detectCUDADevice(void);

		/** 传送三角形信息到GPU */
		void copyTriangleInfosToGPU(void);

		/** 初始化RKDG自由度并将数据传送到GPU上 */
		void initRKDG(void);

		/** rkdg时间推进 */
		//void rkdgAdvance(void);

		void runNext();

		void runAfter();
		/**
		 * 求解时间步长
		 * @param[in] tnum integer 三角单元数目
		 */
		void getTimeStep(int tnum);

		/**
		 *求解守恒量,
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] double_pitch integer 对齐后双精度数组的宽度
		 * @param[in] blocks integer 线程块的数目
		 */
		void calculateConVars(int tnum, int double_pitch, int blocks);

		/**
		 * 处理边界条件
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] num integer 网格单元总数
		 * @param[in] double_pitch integer 对齐后的双精度数组宽度
		 * @param[in] rho double 远场压强
		 * @param[in] rhou double 远场动量
		 * @param[in] rhov double 远场动量
		 * @param[in] rhoE double 远场能量
		 */
		void boundaryCondition(int tnum, int num, int double_pitch, double rho, double rhou, double rhov, double rhoE);

		/**
		 * 计算体积分
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] double_pitch integer 对齐后双精度数组宽度
		 * @param[in] blocks integer 线程块的数目
		 */
		void calculateVolumeRHS(int tnum, int double_pitch, int blocks);

		/**
		 * 计算LF通量系数
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] int_pitch integer 对齐后的整型数组宽度
		 * @param[in] double_pitch integer 对齐后的双精度数组宽度
		 * @param[in] blocks integer 线程块的数目
		 */
		void calculateLFCoeff(int tnum, int int_pitch, int double_pitch, int blocks);

		/**
		 * 计算边上的FG
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] num integer 网格单元总数
		 * @param[in] double_pitch integer 对齐后的双精度数组宽度
		 * @param[in] blocks integer 线程块的数目
		 */
		void calculateEdgeFG(int tnum, int num, int double_pitch, int blocks);

		/**
		 * 计算通量
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] int_pitch integer 对齐后整型数组宽度
		 * @param[in] double_pitch integer 对齐后双精度数组宽度
		 * @param[in] blocks integer 线程块的数目
		 */
		void calculateFlux(int tnum, int int_pitch, int double_pitch, int blocks);

		/**
		 * 计算边上的残差
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] double_pitch integer 对齐后双精度数组宽度
		 * @param[in] blocks integer 线程块的数目
		 */
		void calculateEdgeRHS(int tnum, int double_pitch, int blocks);

		/**
		 * rkdg 时间步第一步推进
		 * @param[in] dt double 时间步长
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] double_pitch integer 对齐后双精度数组宽度
		 * @param[in] blocks integer 线程块的数目
		 */
		void rkdgStepOne(double dt, int tnum, int double_pitch, int blocks);

		/**
		 * rkdg 时间步第一步推进
		 * @param[in] dt double 时间步长
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] double_pitch integer 对齐后双精度数组宽度
		 * @param[in] blocks integer 线程块的数目
		 */
		void rkdgStepTwo(double dt, int tnum, int double_pitch, int blocks);

		/**
		 * rkdg 时间步第一步推进
		 * @param[in] dt double 时间步长
		 * @param[in] tnum integer 三角单元数目
		 * @param[in] double_pitch integer 对齐后双精度数组的宽度
		 * @param[in] blocks integer 线程块的数目
		 */
		void rkdgStepThree(double dt, int tnum, int double_pitch, int blocks);

		/** 
		 * 计算守恒量残差
		 * @param[in] tnum integer 三角单元数目
		 */
		void calculateResidual(int tnum);

		/** 复制自由度到主机 */
		void copyFreedomToHost(void);

		/** 输出解 */
		void outputSolution(int* data_size, double* result_rho,double* result_rhou, double* result_rhov, double* result_rhoE);

		/**
		 * 输出程序配置信息
		 * @param[in] out ostream 输出流
		 */
		void printConfig(ostream& out);

		void commuInfo();   //执行数据交换

		void dealCommuData();    //处理交换的数据
	
		/** 析构函数 */
		~CCUDARkdgSolver();
	protected:
	
	private:
	
};