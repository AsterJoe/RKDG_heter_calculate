/**
 * 本文件是程序配置类的头文件
 * 通过本程序能够实现从文件中读取给定的配置
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */

#pragma once

#include "cppstdheaders.h"
#include "myexception.h"
using namespace std;


extern "C"{
	/**
	 * 配置行的类型: 无法识别的行，注释行，空行，配置项行， 不支持的配置项
	 * 其中不支持的配置行指的是有配置项，但是该配置项在程序中不使用。
	 */
	enum cl_type{
		CL_TYPE_UNKNOWN   = 0,
		CL_TYPE_COMMENT   = 1,
		CL_TYPE_EMPTY     = 2,
		CL_TYPE_CONFITEM  = 3,
		CL_TYPE_UNSUPPORT = 4
	};

}

/**
 * 文件配置类，该类定义哪些配置可以从文件中获取配置值
 * 通过传入的配置文件名称解析配置文件
 * 对于文件解析的结果保存在_message中
 * @example 该类用法如下
 * // 先定义那些项需要从文件中读取
 * string items[] = {"ia", "ib", "ic"}
 * // 使用默认构造函数定义配置项
 * CConfig conf("path/to/configfile", items, 3) // 其中3代表有3个配置项需从文件中读取
 * if ( !conf.parseConfigFile() )
 *	// 对解析不成功的处理
 *  // 输出从文件中读取的值
 * conf.printConfigItems();
 * // 获取文件解析信息
 * conf.getMessage()
 */
class CConfig{
// 属性
protected:
	/**
	 * 用于打印的信息
	 */
	string _message;

	/** 
	 * 标记数组，标记传入的配置项是否都在文件中初始化过
	 */
	bool *item_flag;

	/** 
	 * 标记符号，指示是否所有配置项都成功初始化
	 */
	bool all_items_parsed;

public:
    /**
     * 配置文件名
     */
    string config_file;

	/**
	 * 支持的配置项数组
	 * 如何初始化配置项，请查看 $默认构造函数$ 的说明
	 */
	map<string, string> config_items;

	/** 
	 * 配置项的个数
	 */
	int item_num;
	 
// 方法
protected:
	 /**
     * 构造函数
	 * 该构造函数被禁用，请勿直接调用该函数
     */
    CConfig(){};

	 /** 
	  * 配置文件读取之后检查是否给定的配置项都已经在配置文件中得到了初始化
	  */
	virtual void afterParseConfigFile();

public:
	/**
     * 默认构造函数
	 * @param string filename 配置文件名
	 * @param string[] items 需要从文件中读取的配置项
	 * @param n 配置项的个数
	 */
	CConfig(const string& filename, string items[], int n);

	/**
	 * 解析配置文件
	 * 如果文件不存在，则抛出异常
	 * 解析成功的结果可以通过调用 getMessage() 函数调用
	 * @return boolean 如果给定的配置项全部在配置文件中出现则返回true, 否则返回false
	 */
	bool parseConfigFile(void);

	/**
	 * 解析文件配置行
	 * @param integer lineno 行号
	 * @param string line 行数据
	 * @return integer 返回配置行的类型
	 */
	int parseConfigLine(int lineno, const string& line);

	/** 
	 * 某个配置项是否为类中的配置项
	 * @param string item 该配置项是否在支持的配置项中
	 * @return integer 如果找到配置项，则返回索引，否则返回-1 
	 */
	int configItemIndex(const string& item);

	/**
	 * 打印出信息,请按照需要自定义此函数
	 */
	void printMessage(void);

	/**
	 * 输出读取的值
	 */
	void printConfigItems(void);

	/**
	 * 获取文件解析结果
	 */
	string& getMessage(void){ return _message; };

	/** 
	 * 打印配置错误，并提示程序退出
	 */
	void printQuitMessage(void);

    /**
     * 析构函数
     */
    ~CConfig();

};
