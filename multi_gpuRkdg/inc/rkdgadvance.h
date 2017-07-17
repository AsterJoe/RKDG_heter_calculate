#include "cudarkdgsolver.h"

class RkdgAdvance {
public:
	CCUDARkdgSolver solver;
	bool isCommun;
	double *rho_buffer, *rhou_buffer, *rhov_buffer, *rhoE_buffer;
	double *convar_rho_edge, *convar_rhou_edge, *convar_rhov_edge, *convar_rhoE_edge;
	double *convar_rho_edge_buffer, *convar_rhou_edge_buffer, *convar_rhov_edge_buffer, *convar_rhoE_edge_buffer;

	RkdgAdvance(CCUDARkdgSolver* rkdg_solver);

	void advance();

	void commuInfo(int commu_count);

	void dealCommuData(int i);

	void dealCommuConvar();
};