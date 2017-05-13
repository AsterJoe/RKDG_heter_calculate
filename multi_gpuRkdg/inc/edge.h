/**
 * 该单元定义边
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */
 
#pragma once 
#include "cppstdheaders.h"
using namespace std;


/**
 * 定义单元的边，定义方式为指定边的起点和终点编号
 * @pre 为了能够正常使用该类，请先初始化顶点数组
 */
class CEdge{
	// 属性
	public:
	
	protected:
		int _start; 		/**< 起点编号 */
		int _terminal; 		/**< 终点编号 */
	
	// 行为
	public:
		/**
		 * 默认构造函数
		 * @param[in] start integer 起点
		 * @param[in] terminal integer 终点
		 */
		CEdge(int start=0, int terminal=1);
		
		// setter
		/**
		 * 重新设置边
		 * @param[in] start integer 起点
		 * @param[in] terminal integer 终点
		 */
		void setEdge(int start, int terminal);
		
		/**
		 * 设置边的起点
		 * @param[in] start integer 起点
		 */
		void setStart(int start);
		
		/**
		 * 设置边的终点
		 * @param[in] terminal integer 终点
		 */
		void setTerminal(int terminal);
		
		// getter
		
		/**
		 * 获取起点
		 */
		int getStart(void) const;
		
		/**
		 * 获取终点
		 */
		int getTerminal(void) const;
		
		/** 析构函数 */
		~CEdge();
	private:
};