#include "../inc/myexception.h"

CMyException::CMyException(const string &error):
_msg(error)
{

}


const char * CMyException::what() const throw()
{
	return _msg.c_str();
}

CMyException::~CMyException() throw()
{

}