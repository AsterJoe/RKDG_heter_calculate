/**
 * 定义程序中的常量或者宏定义
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */
 
#pragma once 

 /** 
 * 定义三角形的常量
 */
#ifndef TRIANGLE_CONSTANTS
	#define TRIANGLE_CONSTANTS
	#define TRIANGLE_VERTICES 3  // 三角形顶点个数
	#define TRIANGLE_EDGES 3     // 三角形边数
	#define BASIS_FUNCTIONS 6    // 基函数个数
	#define VOLUME_GPOINTS  7    // 面积分高斯积分节点个数
	#define EDGE_GPOINTS  3		 // 线积分高斯积分节点个数
	#define BASIS_FUNCTION_COEFFS 6  // 基函数系数个数
#endif

// 定义守恒量个数
#ifndef CONSERVATIVE_VARS
#define CONSERVATIVE_VARS 4
#endif

// 定义龙格库塔步数
#ifndef RUNGE_KUTTA_STEPS
#define RUNGE_KUTTA_STEPS 3
#endif

// 设置缓存常数
#ifndef RESIDUAL_VARS
	#define RESIDUAL_VARS 2  // 缓存的变量个数
	#define BUFFER_PORTIONS 2 // 缓冲区有几部分
#endif

// 文件名或者缓冲区的最大长度
#ifndef MAX_LENGTH
	#define MAX_LENGTH 1024 
#endif

/** 
 * 定义单元类型
 */
enum cell_tags{
	CELL_INTERIOR = 0,  // 单元是内部单元
	CELL_FARFIELD  = 1,  // 远场边界条件
	CELL_OUTFLOW  = 2,  // 外流
	CELL_REFLECTION = 3,  // 反射边界
	CELL_INNER_BOUNDARY = 4  //分割线
};