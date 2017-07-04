#include <mpi.h>
#include "cudarkdgsolver.h"
#include<set>

class SendBoundary {
	int myid;
	MPI_Status status;
	MPI_Request request, req[2];
	CCUDARkdgSolver solver;
public:
	SendBoundary(CCUDARkdgSolver* solver);

	void send();
};