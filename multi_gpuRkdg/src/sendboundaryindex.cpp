#include "../inc/sendboundaryindex.h"

SendBoundary::SendBoundary(CCUDARkdgSolver* solver) {
	this->solver = *solver;
}

void SendBoundary::send() {
	for (int i = 0; i < solver.nprocs - solver.grid.area_index - 1; ++i) {
		int* local_exchange_elem = new int[solver.grid.local_innerBoundary_index[i].size()];
		int* neigh_exchange_elem = new int[solver.grid.local_innerBoundary_index[i].size()];
		int index = 0;
		for (int j = 0; j < solver.grid.local_innerBoundary_index[i].size(); j++) {
			local_exchange_elem[index] = solver.grid.elem_index[solver.grid.local_innerBoundary_index[i].at(j)];
			neigh_exchange_elem[index++] = solver.grid.neigh_innerBoundary_index[i].at(j);
		}
		MPI_Isend(local_exchange_elem, solver.grid.local_innerBoundary_index[i].size(), MPI_INT, i + 1, 0, MPI_COMM_WORLD, &request);
		MPI_Isend(neigh_exchange_elem, solver.grid.local_innerBoundary_index[i].size(), MPI_INT, i + 1, 1, MPI_COMM_WORLD, &request);
	}
	int data_size;
	for (int j = 0; j < solver.grid.area_index; ++j) {

		MPI_Probe(j,1,MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &data_size);
		int* local_exchange_elem = new int[data_size];
		int* neigh_exchange_elem = new int[data_size];
		MPI_Irecv(local_exchange_elem, data_size,MPI_INT, j, 1, MPI_COMM_WORLD, &request);
		MPI_Irecv(neigh_exchange_elem, data_size,MPI_INT, j, 0, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
		for (int i = 0; i < data_size; i++) {
			vector<int>::iterator iter;
			iter = find(solver.grid.elem_index.begin(), solver.grid.elem_index.end(), local_exchange_elem[i]);
			solver.grid.local_innerBoundary_index[j].push_back(distance(solver.grid.elem_index.begin(), iter));
			solver.grid.neigh_innerBoundary_index[j].push_back(neigh_exchange_elem[i]);
		}
	}
} 