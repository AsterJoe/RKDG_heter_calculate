#include "../inc/vertice2d.h"


CVertice2D::CVertice2D( double x/*=0*/, double y/*=0*/ ):
_x(x),_y(y)
{

}

void CVertice2D::setVertice(double x, double y)
{
	_x = x;
	_y = y;
}

void CVertice2D::setX(double x)
{
	_x = x;
}

void CVertice2D::setY(double y)
{
	_y = y;
}

double CVertice2D::getX(void) const
{
	return _x;
}

double CVertice2D::getY(void) const
{
	return _y;
}

CVertice2D::~CVertice2D()
{

}
