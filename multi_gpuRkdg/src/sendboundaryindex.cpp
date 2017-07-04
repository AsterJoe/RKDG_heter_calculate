#include "../inc/sendboundaryindex.h"

SendBoundary::SendBoundary(CCUDARkdgSolver* solver) {
	this->solver = *solver;
	myid = solver->myid;
}

void SendBoundary::send() {
	for (int i = 0; i < solver.nprocs - myid - 1; ++i) {
		int* local_exchange_elem = new int[solver.grid.local_innerBoundary_index[i + myid].size()];
		int* neigh_exchange_elem = new int[solver.grid.local_innerBoundary_index[i + myid].size()];

		for (int j = 0; j < solver.grid.local_innerBoundary_index[i + myid].size(); j++) {
			local_exchange_elem[j] = solver.grid.elem_index[solver.grid.local_innerBoundary_index[i + myid].at(j)];
			neigh_exchange_elem[j] = solver.grid.neigh_innerBoundary_index[i+ myid].at(j);
		}
		MPI_Isend(local_exchange_elem, solver.grid.local_innerBoundary_index[i + myid].size(), MPI_INT, myid + i + 1, 0, MPI_COMM_WORLD, &request);
		MPI_Isend(neigh_exchange_elem, solver.grid.local_innerBoundary_index[i + myid].size(), MPI_INT, myid + i + 1, 1, MPI_COMM_WORLD, &request);
	}
	int data_size;
	for (int j = 0; j < myid; ++j) {

		MPI_Probe(j,0,MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &data_size);
		int* local_exchange_elem = new int[data_size];
		int* neigh_exchange_elem = new int[data_size];
		MPI_Irecv(local_exchange_elem, data_size,MPI_INT, j, 1, MPI_COMM_WORLD, &req[0]);
		MPI_Irecv(neigh_exchange_elem, data_size,MPI_INT, j, 0, MPI_COMM_WORLD, &req[1]);
		
		MPI_Waitall(2, req, &status);
		for (int i = 0; i < data_size; i++) {
			vector<int>::iterator iter;
			iter = find(solver.grid.elem_index.begin(), solver.grid.elem_index.end(), local_exchange_elem[i]);
			solver.grid.local_innerBoundary_index[j].push_back(distance(solver.grid.elem_index.begin(), iter));
			solver.grid.neigh_innerBoundary_index[j].push_back(neigh_exchange_elem[i]);
		}
	}
} 