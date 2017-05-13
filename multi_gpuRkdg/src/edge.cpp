#include "../inc/edge.h"

CEdge::CEdge(int start, int terminal):
_start(start),
_terminal(terminal)
{

}

int CEdge::getStart(void) const
{
	return _start;
}

int CEdge::getTerminal(void) const
{
	return _terminal;
}

void CEdge::setEdge(int start, int terminal)
{
	_start = start;
	_terminal = terminal;
}

void CEdge::setStart(int start)
{
	_start = start;
}

void CEdge::setTerminal(int terminal)
{
	_terminal = terminal;
}

CEdge::~CEdge()
{

}