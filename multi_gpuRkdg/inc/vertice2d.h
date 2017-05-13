/**
 * 该文件定义二维顶点类
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
 * 二维顶点类，类属性为顶点的x, y坐标
 * 请通过接口访问类属性
 */
class CVertice2D
{
	// 属性
	public:
	
	protected:
		double _x; /**< 顶点的横坐标 */
		double _y; /**< 顶点的纵坐标 */
	
	// 操作
	public:
		/**
		 * 默认构造函数
		 * @param[in] x double 顶点的横坐标
		 * @param[in] y double 顶点的纵坐标
		 */
		CVertice2D(double x=0, double y=0);

		// setter
		/**
		 * 设置顶点坐标
		 * @param[in] x double 顶点的新横坐标
		 * @param[in] y double 顶点的新纵坐标
		 */
		void setVertice(double x, double y);
		
		/**
		 * 设置顶点的横坐标
		 * @param[in] x double 顶点横坐标
		 */
		void setX(double x);
		
		/**
		 * 设置顶点的纵坐标
		 * @param[in] y double 顶点纵坐标
		 */
		void setY(double y);
		
		// getter
		
		/**
		 * 获取顶点的横坐标
		 * @return double 顶点的横坐标
		 */
		double getX(void) const;
		
		/**
		 * 获取顶点的纵坐标
		 * @return double 顶点的纵坐标
		 */
		double getY(void) const;
		
		
		/** 析构函数 */
		~CVertice2D();
	private:
	
};