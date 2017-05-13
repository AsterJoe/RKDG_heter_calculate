#include "../inc/mytime.h"

CMyTime::CMyTime():
_cpu_start(0),
_wall_start(0),
_cpu_end(0),
_wall_end(0)
{}

void CMyTime::beginTimer( void )
{
	_cpu_start = clock();
	time(&_wall_start);
}

void CMyTime::endTimer( void )
{
	_cpu_end = clock();
	time(&_wall_end);
}

double CMyTime::getCPUElapsedTime( void )
{
	return static_cast<double>(_cpu_end-_cpu_start) / CLOCKS_PER_SEC;
}

double CMyTime::getWallElapsedTime( void )
{
	return difftime( _wall_end, _wall_start );
}

string CMyTime::getCurrentTime( void )
{
	char buffer[MAX_LENGTH];
	
	time_t tms = time(NULL);
	tm *tml = localtime(&tms);
	strftime(buffer, MAX_LENGTH, "%Y-%m-%d %H:%M:%S", tml);

	return string(buffer);
}

CMyTime::~CMyTime()
{}