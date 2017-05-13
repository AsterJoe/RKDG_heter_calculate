#include "../inc/config.h"

// 配置类的默认构造函数
CConfig::CConfig(const string& filename, string items[], int n):
_message(""),
config_file(filename),
item_flag(NULL),
item_num(n),
all_items_parsed(true)
{
	item_flag = new bool[n];
	
	for ( int i=0; i<n; ++i )
	{
		item_flag[i] = false;
		config_items[items[i]] = "";
	}
}

void CConfig::printMessage(void)
{
	cout<<_message<<endl;
}

bool CConfig::parseConfigFile(void)
{

	ifstream fin(config_file.c_str());
	if (!fin)
		throw CMyException("\nFailed to open configure file "+config_file);

	int line_number = 0;
	string tmp_line;

	_message += "\n======================================\n";
	_message += "    Parse file   "+config_file;
	_message += " \n======================================\n\n";

	while(getline(fin, tmp_line))
	{
		++ line_number;
		parseConfigLine(line_number,tmp_line);
	}

	afterParseConfigFile();

	// 判断是否给出的配置项都在文件中读取完毕
	_message += "\n======================================\n";
	_message += "   Parse file finished!\n";
	_message += "======================================\n";

	fin.close();

	return all_items_parsed;
}

int CConfig::parseConfigLine(int lineno, const string& line)
{
	char buffer[256];

	// 去掉行首空格
	string::size_type first_index = line.find_first_not_of(" \t");

	// 处理空行和注释行
	if (first_index==string::npos)
		return CL_TYPE_EMPTY;
	if (line.at(first_index)=='#')
		return CL_TYPE_COMMENT;

	// 找出配置项名称
	string separator = ",; #\t";
	string::size_type mid_index = line.find_first_of(" =\t", first_index);
	// 该行给出了配置项名称却没有给出配置值
	if (mid_index==string::npos)
		return CL_TYPE_UNKNOWN;
	string item(line, first_index, mid_index-first_index);
	int index = configItemIndex(item);
	if ( -1==index )
	{
		sprintf(buffer, "note: unsurpported item found in line %d : %s\n", lineno, item.c_str());
		_message += string(buffer);
		return CL_TYPE_UNSUPPORT;
	}

	// 去掉配置项值的前导空格或者等于
	string::size_type another_index = line.find_first_not_of(" =\t", mid_index);
	// 该行给出了配置项名称却没有给出配置值
	if (another_index==string::npos)
		return CL_TYPE_UNKNOWN;


	// 找到配置值的终止位置
	string::size_type final_index = line.find_first_of(separator, another_index);
	string value;
	if ( string::npos==final_index )
	{
		value = line.substr(another_index);
	} 
	else if ( final_index==another_index )
	{
		// 没有给出配置值
		return CL_TYPE_UNKNOWN;
	}
	else
	{
		value = line.substr(another_index,final_index-another_index);
	}

	config_items[item] = value;
	item_flag[index] = true;
	return CL_TYPE_CONFITEM;
}


int CConfig::configItemIndex(const string& item)
{
	map<string,string>::iterator iter;
	int i;
	for ( i=0,iter=config_items.begin(); iter!=config_items.end(); ++iter,++i )
	{
		if ( iter->first==item )
			return i;
	}

	return -1;
}

void CConfig::printConfigItems(void)
{
	map<string,string>::iterator iter;
	cout<<"\n======================================\n";
	cout<<"    configures in file "+config_file+":\n";
	cout<<"======================================\n";
	int i=0;
	for ( iter=config_items.begin(); iter!=config_items.end(); ++iter, ++i )
	{
		if ( item_flag[i] )
			cout<<iter->first<<" : "<<iter->second<<endl;
	}
	cout<<"\n======================================\n";
}


void CConfig::afterParseConfigFile(void)
{
	map<string,string>::iterator iter;
	int i;
	string tmp;
//	bool unfinished = false;
	tmp += "\n===============  WARNING  ==============\n\n";
	for ( i=0, iter=config_items.begin(); iter!=config_items.end(); ++iter, ++i )
	{
		if ( !item_flag[i] )
		{
//			unfinished = true;
			tmp += "configure item: "+iter->first+" is not be initialized.\n";
			all_items_parsed = false;
		}
	}

	if ( !all_items_parsed )
		_message += tmp;
}

void CConfig::printQuitMessage(void)
{
	printMessage();
	cout<<"Some configure items are not provided initial value."<<endl;
	cout<<"press ENTER to exit."<<endl;
	getchar();
	exit(-1);
}

CConfig::~CConfig()
{
	if ( item_flag!=NULL )
	{
		delete[] item_flag;
	}
}