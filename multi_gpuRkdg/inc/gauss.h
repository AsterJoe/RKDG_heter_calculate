/**
 * 高斯列主元消元法程序
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */

#pragma once
#include "cppstdheaders.h"

using namespace std;

extern "C"{

/**
 * 列主元素的Gauss消去法求解简单的线性方程组
 *
 * @param[in]  N : 左端矩阵行数/列数
 * @param[in]  A : 左端矩阵(以一维向量方式存储)
 * @param[in]  B : 右端向量,求解方程组完毕后存放解向量
 */
void Gauss(int N, double* A, double* B);

}