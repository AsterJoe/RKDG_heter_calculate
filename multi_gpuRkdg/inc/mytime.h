/** 
 * 自定义的时间类
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */

#pragma once

#ifndef MAX_LENGTH
#define MAX_LENGTH 256
#endif

#include "cppstdheaders.h"
using namespace std;

/**
 * 自定义的时间类，可实现CPU计时，墙上时间计时以及格式化本地时间
 */
class CMyTime {
	
	// 属性
	protected:
		clock_t _cpu_start;			/**< 程序开始的处理器时间 */
		time_t _wall_start; 		/**< 程序开始的墙上时间 */
		clock_t _cpu_end;			/**< 程序结束的处理器时间 */
		time_t _wall_end;			/**< 程序结束的墙上时间 */
	
	public:
	
	// 方法
	protected:
	
	public:
		// 构造函数
		CMyTime();
		
		
		// 计时开始
		void beginTimer( void );
		
		// 计时结束
		void endTimer( void );
		
		/**
		 * 获取CPU时间
		 * @return 以秒形式返回的CPU时间
		 */
		double getCPUElapsedTime( void );
		
		/**
		 * 获取墙上时间
		 * @return 以秒形式返回的墙上时间
		 */
		double getWallElapsedTime( void );
		
		
		/**
		 * 获取当前时刻，以字符串表示
		 */
		string getCurrentTime( void );
		
		// 析构函数
		~CMyTime();
};