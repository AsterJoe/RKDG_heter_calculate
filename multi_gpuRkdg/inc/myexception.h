/**
 * 自定义的异常类
 * 该类捕捉除了标准异常类外的异常
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */
 

#pragma once

#include "cppstdheaders.h"
using namespace std;

/** 
 * 自定义的异常类
 */
class CMyException: public exception
{
	// 属性
	public:

	protected:
		string _msg;		/**< 错误信息描述 */

	// 方法
	public:
		// 构造函数
		CMyException(const string &error="Unknow exception.");
		
		/**< 重写基函数的what方法 */
		const char *what()const throw();

		/** 析构函数 */
		~CMyException() throw();
	protected:
};