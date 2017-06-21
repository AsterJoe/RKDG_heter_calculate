#include <mpi.h>
#include "cudarkdgsolver.h"
#include<set>

class SendBoundary {
	MPI_Status status;
	MPI_Request request;
	CCUDARkdgSolver solver;
public:
	SendBoundary(CCUDARkdgSolver* solver);

	void send();
};