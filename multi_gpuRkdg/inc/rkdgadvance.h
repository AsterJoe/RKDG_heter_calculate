#include "cudarkdgsolver.h"

class RkdgAdvance {
public:
	CCUDARkdgSolver solver;

	RkdgAdvance(CCUDARkdgSolver* rkdg_solver);

	void advance();

	void commuInfo();

	void dealCommuData();
};