#include "cudarkdgsolver.h"
#define MAX_INNER_BOUNDARY_SIZE 50

class RkdgAdvance {
public:
	CCUDARkdgSolver solver;
	//int buffer_num;
	bool isCommun;
	double *rho_buffer, *rhou_buffer, *rhov_buffer, *rhoE_buffer;
	double *convar_rho_edge, *convar_rhou_edge, *convar_rhov_edge, *convar_rhoE_edge;
	double **convar_rho_edge_buffer, **convar_rhou_edge_buffer, **convar_rhov_edge_buffer, **convar_rhoE_edge_buffer;

	RkdgAdvance(CCUDARkdgSolver* rkdg_solver);

	void advance();

	void commuInfo(int commu_count);

	void dealCommuData(int i);
	void dealCommuConvar(int buffer_num);
	//void dealCommuConvar(int buffer_num, double **convar_rho_edge_buffer, double **convar_rhou_edge_buffer, double **convar_rhov_edge_buffer, double **convar_rhoE_edge_buffer);

	void commuConvar(int size, int pitch, int buffer_num, int commu_count, int *innerBoudary_count);
};