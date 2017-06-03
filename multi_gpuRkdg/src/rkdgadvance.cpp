#include <mpi.h>
#include "../inc/rkdgadvance.h"

RkdgAdvance::RkdgAdvance(CCUDARkdgSolver* rkdg_solver) {
	solver = *rkdg_solver;
}

void RkdgAdvance::advance() {
	ofstream fout;
	if ( solver.log_history=='Y' )
	{
		fout.open(solver.residual_file.c_str());
		if ( !fout )
			throw CMyException("Failed to open residual log file: "+solver.residual_file);
		
		fout<<"N, rho"<<endl;
	}
			
			
	double nt(0);
	int count(0);
		
	/*int tnum = grid.getTriangleNumber();
	int num  = grid.getCellNumber();*/
	int tnum = solver.grid.getLocalTriangleNumber();
	int num = solver.grid.getLocalCellNumber();
		
		
	int blocks = (tnum%solver.threads_per_block) ? tnum/solver.threads_per_block+1 : tnum/solver.threads_per_block;
		
	double ut   = sqrt(solver.gamma*solver.pref/solver.rhoref)*solver.mach;
	double rhou = solver.rhoref*ut*cos(solver.alpha);
	double rhov = solver.rhoref*ut*sin(solver.alpha);
	double rhoE = 0.5*solver.rhoref*(ut*ut) + solver.pref/(solver.gamma-1);
		
	bool copy(false);
		
	cudaError_t error;
	size_t pitch = solver._cuarrays.getDoublePitch();
	int pitch_num = pitch / sizeof(double);
	int ipitch_num = solver._cuarrays.getIntPitch() / sizeof(int);
		
	cudaEvent_t time_start, time_stop;
		
	cudaEventCreateWithFlags(&time_start, cudaEventDisableTiming|cudaEventBlockingSync);
	cudaEventCreateWithFlags(&time_stop,  cudaEventDisableTiming|cudaEventBlockingSync);
		
	if ( solver.log_history=='Y' )
		copy = true;
		
	// ȷ��֮ǰCUDA�ĳ�ʼ���������Ѿ����
	cudaDeviceSynchronize();
		
	do 
	{
		++ count;
				
		cudaEventRecord(time_start);
		// ���㵱ǰʱ�䲽��
		solver.getTimeStep(tnum);
		
		cudaEventRecord(time_stop);
		
		// ��������ɶ�
		/*cudaMemcpy2DAsync(_cuarrays.freedom_rho_old,  pitch, _cuarrays.freedom_rho,  pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(_cuarrays.freedom_rhou_old, pitch, _cuarrays.freedom_rhou, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(_cuarrays.freedom_rhov_old, pitch, _cuarrays.freedom_rhov, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							 
		cudaMemcpy2DAsync(_cuarrays.freedom_rhoE_old, pitch, _cuarrays.freedom_rhoE, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);*/
					
		cudaMemcpy2DAsync(solver._freedom_rho, pitch, solver._cuarrays.freedom_rho, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
		cudaMemcpy2DAsync(solver._freedom_rho, pitch, solver._cuarrays.freedom_rho, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
		cudaMemcpy2DAsync(solver._freedom_rho, pitch, solver._cuarrays.freedom_rho, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
		cudaMemcpy2DAsync(solver._freedom_rho, pitch, solver._cuarrays.freedom_rho, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
		
				
		
		MPI_Request request[4], request1;
		MPI_Status status;
		int x = 1;         //֪ͨ�����ڵ���Ϣ������0��ʾ����
		if(solver.myid = 0 && count % 5 == 0) {
			for (int i = 1; i < solver.nprocs; i++) {
				MPI_Isend(&x, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request1);
			}
			commuInfo();
		} else if(solver.myid != 0) {
			int flag;
			MPI_Status status;
			if (count == 0) {
				MPI_Irecv(&x,1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request1);
			} 
			MPI_Test(&request1, &flag, &status);
			if (flag == 1 && x == 1) {
				commuInfo();
				MPI_Irecv(&x,1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request1);
			}
		}
		
				
		for ( int i=0; i<RUNGE_KUTTA_STEPS; ++i )
		{
			// �����غ�����ֵ
			solver.calculateConVars(tnum, pitch_num, blocks);
		
			// ����߽�����
			solver.boundaryCondition(tnum, num, pitch_num, solver.rhoref, rhou, rhov, rhoE);
		
			// ��������ֲв�
			solver.calculateVolumeRHS(tnum, pitch_num, blocks);
		
			// ����LFͨ��ϵ��
			solver.calculateLFCoeff(tnum, ipitch_num, pitch_num, blocks);
		
			// ����f, g�ڱ��ϵ�ֵ
			solver.calculateEdgeFG(tnum, num, pitch_num, blocks);
		
			solver.calculateFlux(tnum, ipitch_num, pitch_num, blocks);
		
			// �����߻��ֲв�
			solver.calculateEdgeRHS(tnum, pitch_num, blocks);
		
			// ʱ���ƽ�
			switch (i)
			{
			case 0:
				cudaEventSynchronize(time_stop);
		
				// ��ʱ�䲽�����͵�����
				cudaMemcpy(solver._dt, solver._cuarrays.ddt, sizeof(double), cudaMemcpyDeviceToHost);
		
				if ( 0==(count-1)%solver.print_interval )
					cout<<"Step: "<<count<<", time step: "<<solver._dt[0]<<endl;
		
				if ( (solver._terminal_time-nt)<solver._dt[0] )
				{
					solver._dt[0] = solver._terminal_time -  nt;
				}
		
				// ʱ�䲽�ƽ�
				solver.rkdgStepOne(solver._dt[0], tnum, pitch_num, blocks);
		
				break;
		
			case 1:
				solver.rkdgStepTwo(solver._dt[0], tnum, pitch_num, blocks);
				break;
		
			case 2:
				solver.rkdgStepThree(solver._dt[0], tnum, pitch_num, blocks);
				break;
		
			default:
				throw CMyException("impossible case!");
				break;
			}
		}
		
				
		if ( copy && (count-1) )
		{
			// ���Ʋв�����
			cudaMemcpy(solver._residual, solver._cuarrays.residual,
						sizeof(double)*RESIDUAL_VARS, cudaMemcpyDeviceToHost);
		
			if ( 0==(count-1)%solver.print_interval )
				cout<<"Current time: "<<nt<<"   rhomax: "<<solver._residual[0]/solver.rhoref<<"   E: "<<solver._residual[1]/rhoE<<endl;
		
			fout<<count<<"   "<<log(solver._residual[0]/solver.rhoref)/log(10.0)<<endl;
		}
		
		// ����в�
		solver.calculateResidual(tnum);
		
				
		// ��ʱ�ƽ�
		nt += solver._dt[0];
		
		error = cudaPeekAtLastError();
		if ( error!=cudaSuccess )
			throw CMyException(cudaGetErrorString(error));
		
	} while ( nt<solver._terminal_time );
		
	cudaDeviceSynchronize();
			
	if ( copy )
	{
		// ���Ʋв�����
		cudaMemcpy(solver._residual, solver._cuarrays.residual,
			sizeof(double)*RESIDUAL_VARS, cudaMemcpyDeviceToHost);
		
		if ( 0==(count-1)%solver.print_interval )
			cout<<"��ǰʱ�䣺 "<<nt-solver._dt[0]<<"   rhomax: "<<solver._residual[0]/solver.rhoref<<"   E: "<<solver._residual[1]/rhoE<<endl;
		
		fout<<count<<"   "<<log(solver._residual[0]/solver.rhoref)/log(10.0)<<endl;
	}
		
		
	cudaEventDestroy(time_start);
	cudaEventDestroy(time_stop);
		
	if ( solver.log_history=='Y' )
		fout.close();
}


void RkdgAdvance::commuInfo() {
	MPI_Request request;
			
	int *rho_buffer, *rhou_buffer, *rhov_buffer, *rhoE_buffer;
	for (int i = 0; i < solver.nprocs - 1; i++) {
		int num = solver.grid.local_innerBoundary_index[i].size();
		rho_buffer = new int[num];
		rhou_buffer = new int[num];
		rhov_buffer = new int[num];
		rhoE_buffer = new int[num];
		for (int j = 0; j < num; j++) {
			rho_buffer[j] = solver._cuarrays.freedom_rho[solver.grid.local_innerBoundary_index[i].at(j)];
			rhou_buffer[j] = solver._cuarrays.freedom_rhou[solver.grid.local_innerBoundary_index[i].at(j)];
			rhov_buffer[j] = solver._cuarrays.freedom_rhov[solver.grid.local_innerBoundary_index[i].at(j)];
			rhoE_buffer[j] = solver._cuarrays.freedom_rhoE[solver.grid.local_innerBoundary_index[i].at(j)];
		}
		int dest = i < solver.myid ? i : i + 1;
		MPI_Isend(rho_buffer, num, MPI_DOUBLE, dest, 1 + 4 * solver.num_commu, MPI_COMM_WORLD, &request);
		MPI_Isend(rhou_buffer, num, MPI_DOUBLE, dest, 2 + 4 * solver.num_commu, MPI_COMM_WORLD, &request);
		MPI_Isend(rhov_buffer, num, MPI_DOUBLE, dest, 3 + 4 * solver.num_commu, MPI_COMM_WORLD, &request);
		MPI_Isend(rhoE_buffer, num, MPI_DOUBLE, dest, 4 + 4 * solver.num_commu, MPI_COMM_WORLD, &request);
				
		MPI_Irecv(rho_buffer, num, MPI_DOUBLE, dest, 1 + 4 * solver.num_commu, MPI_COMM_WORLD, &request);
		MPI_Irecv(rhou_buffer, num, MPI_DOUBLE, dest, 2 + 4 * solver.num_commu, MPI_COMM_WORLD, &request);
		MPI_Irecv(rhov_buffer, num, MPI_DOUBLE, dest, 3 + 4 * solver.num_commu, MPI_COMM_WORLD, &request);
		MPI_Irecv(rhoE_buffer, num, MPI_DOUBLE, dest, 4 + 4 * solver.num_commu++, MPI_COMM_WORLD, &request);	
	}
	MPI_Barrier(MPI_COMM_WORLD);
	dealCommuData();
}

void RkdgAdvance::dealCommuData() {
	size_t pitch = solver._cuarrays.getDoublePitch();
	for (int i = 0; i < solver.nprocs - 1; i++) {
		for (int j = 0; j < solver.grid.local_innerBoundary_index[i].size(); j++) {
			solver._freedom_rho[solver.grid.local_innerBoundary_index[i].at(j)] =  solver.rho_buffer[j];
			solver._freedom_rhou[solver.grid.local_innerBoundary_index[i].at(j)] = solver.rhou_buffer[j];
			solver._freedom_rhov[solver.grid.local_innerBoundary_index[i].at(j)] = solver.rhov_buffer[j];
			solver._freedom_rhoE[solver.grid.local_innerBoundary_index[i].at(j)] = solver.rhoE_buffer[j];
		}
	}
	
	cudaMemcpy2DAsync(solver._cuarrays.freedom_rho, pitch, solver._freedom_rho, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
	cudaMemcpy2DAsync(solver._cuarrays.freedom_rhou, pitch, solver._freedom_rhou, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
	cudaMemcpy2DAsync(solver._cuarrays.freedom_rhov, pitch, solver._freedom_rhov, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
	cudaMemcpy2DAsync(solver._cuarrays.freedom_rhoE, pitch, solver._freedom_rhoE, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
}
