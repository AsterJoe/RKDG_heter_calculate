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
		
	// 确保之前CUDA的初始化工作都已经完成
	cudaDeviceSynchronize();
	int commu_count(0);	
	int continued(0);
	int depart_num  = solver.grid.local_innerBoundary_index[0].size();
	
	convar_rho_edge = (double*)malloc(num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
	convar_rhou_edge = (double*)malloc(num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
	convar_rhov_edge = (double*)malloc(num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
	convar_rhoE_edge = (double*)malloc(num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));

	/*convar_rho_edge_buffer = (double*)malloc(depart_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
	convar_rhou_edge_buffer = (double*)malloc(depart_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
	convar_rhov_edge_buffer = (double*)malloc(depart_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
	convar_rhoE_edge_buffer = (double*)malloc(depart_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
*/
	int buffer_num(100);
	rho_buffer = new double*[solver.nprocs - 1];
	rhou_buffer = new double*[solver.nprocs - 1];
	rhov_buffer = new double*[solver.nprocs - 1];
	rhoE_buffer = new double*[solver.nprocs - 1];

	convar_rho_edge_buffer = new double*[solver.nprocs - 1];
	convar_rhou_edge_buffer = new double*[solver.nprocs - 1];
	convar_rhov_edge_buffer = new double*[solver.nprocs - 1];
	convar_rhoE_edge_buffer = new double*[solver.nprocs - 1];

	for (int i = 0; i < solver.nprocs - 1; i++) {
		rho_buffer[i] = (double*)malloc(buffer_num * BASIS_FUNCTIONS * sizeof(double));
		rhou_buffer[i] = (double*)malloc(buffer_num * BASIS_FUNCTIONS * sizeof(double));
		rhov_buffer[i] = (double*)malloc(buffer_num * BASIS_FUNCTIONS * sizeof(double));
		rhoE_buffer[i] = (double*)malloc(buffer_num * BASIS_FUNCTIONS * sizeof(double));
		
		convar_rho_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof(double));
		convar_rhou_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof(double));
		convar_rhov_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof(double));
		convar_rhoE_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof(double));
	}

	int x = 1;
	do 
	{
		++ count;
		isCommun = false;
		cudaEventRecord(time_start);
		// 计算当前时间步长
		solver.getTimeStep(tnum);
		cout<<"7"<<endl;
		//getDt();
		cout<<"8"<<endl;
		cudaEventRecord(time_stop);
		
		// 保存旧自由度
		cudaMemcpy2DAsync(solver._cuarrays.freedom_rho_old,  pitch, solver._cuarrays.freedom_rho,  pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(solver._cuarrays.freedom_rhou_old, pitch, solver._cuarrays.freedom_rhou, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(solver._cuarrays.freedom_rhov_old, pitch, solver._cuarrays.freedom_rhov, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							 
		cudaMemcpy2DAsync(solver._cuarrays.freedom_rhoE_old, pitch, solver._cuarrays.freedom_rhoE, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);
		cudaDeviceSynchronize();
		size_t size = sizeof(double) * num;
		MPI_Request request[4], request1, request2;
		MPI_Status status;
		if (count == 1) {
			isCommun = true;
		}else if (solver.myid == 0 && count % 2 == 0) {
			for (int i = 1; i < solver.nprocs; i++) {
				MPI_Isend(&x, 1, MPI_INT, i, solver.nprocs + 13 * commu_count, MPI_COMM_WORLD, &request2);
			}
			commu_count++;
			isCommun = true;
		} else if(solver.myid != 0) {
			int flag(0);
			if (count == 1) {
				MPI_Irecv(&x, 1, MPI_INT, 0, solver.nprocs + 13 * commu_count++, MPI_COMM_WORLD, &request2);
			}
			MPI_Test(&request2, &flag, &status);
			if (flag == 1 && x == 1) {
				isCommun = true;
			} else if (flag == 1 && x == 0) {
				break;
			}
		}
		for ( int i=0; i<RUNGE_KUTTA_STEPS; ++i )
		{
			// 计算守恒量的值
			solver.calculateConVars(tnum, pitch_num, blocks);
			
			// 处理边界条件
			int local_inner_count(0);
			for (int t = 0; t < solver.nprocs - 1; t++)
			{
				local_inner_count += solver.grid.local_innerBoundary_index[t].size(); 
			}
			solver.boundaryCondition(tnum + local_inner_count, num, pitch_num, solver.rhoref, rhou, rhov, rhoE);
			if (isCommun) {
				commuInfo(commu_count);
			}
			// 计算体积分残差
			solver.calculateVolumeRHS(tnum, pitch_num, blocks);

			// 计算LF通量系数
			solver.calculateLFCoeff(tnum, ipitch_num, pitch_num, blocks);
		
			// 计算f, g在边上的值
			solver.calculateEdgeFG(tnum, num, pitch_num, blocks);
		
			solver.calculateFlux(tnum, ipitch_num, pitch_num, blocks);

			// 计算线积分残差
			solver.calculateEdgeRHS(tnum, pitch_num, blocks);

			// 时间推进
			switch (i)
			{
			case 0:
				cudaEventSynchronize(time_stop);
		
				// 将时间步长传送到本地
				cudaMemcpy(solver._dt, solver._cuarrays.ddt, sizeof(double), cudaMemcpyDeviceToHost);
		
				if ( 0==(count-1)%solver.print_interval )
					cout<<"Step: "<<count<<", time step: "<<solver._dt[0]<<endl;
		
				/*if ( (solver._terminal_time-nt)<solver._dt[0] )
				{
					solver._dt[0] = solver._terminal_time -  nt;
				}*/
		
				// 时间步推进
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
		cout<<"dt:"<<solver._dt[0]<<endl;
		if (isCommun && solver.myid != 0) {
			MPI_Irecv(&x,1, MPI_INT, 0, solver.nprocs + 13 * commu_count++, MPI_COMM_WORLD, &request2);
		}

		if ( copy && (count-1) )
		{
			// 复制残差数据
			cudaMemcpy(solver._residual, solver._cuarrays.residual,
						sizeof(double)*RESIDUAL_VARS, cudaMemcpyDeviceToHost);
		
			if ( 0==(count-1)%solver.print_interval )
				cout<<"Current time: "<<nt<<"   rhomax: "<<solver._residual[0]/solver.rhoref<<"   E: "<<solver._residual[1]/rhoE<<endl;
		
			fout<<count<<"   "<<log(solver._residual[0]/solver.rhoref)/log(10.0)<<endl;
		}
		// 计算残差
		solver.calculateResidual(tnum);
		//cout<<"calculateResidual"<<endl;
				
		// 计时推进
		nt += solver._dt[0];
		
		error = cudaPeekAtLastError();
		//cout<<"error test"<<endl;
		if ( error!=cudaSuccess ) {
			cout<<"do error"<<endl;
			throw CMyException(cudaGetErrorString(error));
		}
		
		/*int end(0);
		if (nt<solver._terminal_time)
		{
		end = 1;
		}
		MPI_Reduce(&end, &continued, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		cout<<"reduce end"<<endl;
		MPI_Bcast(&continued, 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout<<"continued:"<<continued<<endl;*/

		if (solver.myid == 0 && count > 10) {
			x = 0;
			for (int i = 1; i < solver.nprocs; i++) {
				MPI_Isend(&x, 1, MPI_INT, i, solver.nprocs + 13 * commu_count, MPI_COMM_WORLD, &request2);
			}
			break;
		}
	} while ( true );
	cout<<solver.myid<<" process end while, nt"<<nt<<", terminal:"<<solver._terminal_time<<endl;	
	cudaDeviceSynchronize();
			
	if ( copy )
	{
		// 复制残差数据
		cudaMemcpy(solver._residual, solver._cuarrays.residual,
			sizeof(double)*RESIDUAL_VARS, cudaMemcpyDeviceToHost);
		if ( 0==(count-1)%solver.print_interval )
			cout<<"当前时间： "<<nt-solver._dt[0]<<"   rhomax: "<<solver._residual[0]/solver.rhoref<<"   E: "<<solver._residual[1]/rhoE<<endl;
		
		fout<<count<<"   "<<log(solver._residual[0]/solver.rhoref)/log(10.0)<<endl;
	}	
		
	cudaEventDestroy(time_start);
	cudaEventDestroy(time_stop);
		
	if ( solver.log_history=='Y' )
		fout.close();
}

double RkdgAdvance::getDt() {
	double tmx(0);
	for (int i = 0; i < solver.grid.getLocalTriangleNumber(); i++) {
		double rho = solver._freedom_rho[i];
		double u   = solver._freedom_rhou[i] / rho;
		double v   = solver._freedom_rhov[i] / rho;
		double p   = (solver.gamma-1)*(solver._freedom_rhoE[i]-0.5*rho*(u*u+v*v));
		double a   = sqrt(solver.gamma*p/rho);

		double dx  = (sqrt(u*u+v*v)+a)*solver.grid.triangle_infos.perimeter[i]/solver.grid.triangle_infos.area[i];

		if ( dx>tmx )
		{
			tmx = dx;
		}
	}
	//solver._dt[0] = solver.cfl / tmx;
	cout<<"tmx:"<<tmx<<endl;
	cout<<"my dt:" << solver.cfl / tmx<<endl;
	return solver.cfl / tmx;
}

void RkdgAdvance::commuInfo(int commu_count) {
	int buffer_num(100);
	int num = solver.grid.getLocalCellNumber();
	size_t pitch = solver._cuarrays.getDoublePitch();
	size_t size = sizeof(double) * num;
	cudaMemcpy2DAsync(solver._freedom_rho,  size, solver._cuarrays.freedom_rho,  pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			  
	cudaMemcpy2DAsync(solver._freedom_rhou, size, solver._cuarrays.freedom_rhou, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			  
	cudaMemcpy2DAsync(solver._freedom_rhov, size, solver._cuarrays.freedom_rhov, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			 
	cudaMemcpy2DAsync(solver._freedom_rhoE, size, solver._cuarrays.freedom_rhoE, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
	
	cudaMemcpy2DAsync(convar_rho_edge, size, solver._cuarrays.convar_rho_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);			  
	cudaMemcpy2DAsync(convar_rhou_edge, size, solver._cuarrays.convar_rhou_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);			 
	cudaMemcpy2DAsync(convar_rhov_edge, size, solver._cuarrays.convar_rhov_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);
	cudaMemcpy2DAsync(convar_rhoE_edge, size, solver._cuarrays.convar_rhoE_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	
	isCommun = true;
	MPI_Request *request,req, request1;
	MPI_Status status[4];	
	request = new MPI_Request[8 * (solver.nprocs - 1)];
	for (int i = 0; i < solver.nprocs - 1; i++) {
		int num = solver.grid.local_innerBoundary_index[i].size();
		for (int t = 0; t < BASIS_FUNCTIONS; t++) {
			for (int j = 0; j < num; j++) {
				rho_buffer[i][t * buffer_num + j] = solver._freedom_rho[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
				rhou_buffer[i][t * buffer_num + j] = solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
				rhov_buffer[i][t * buffer_num + j] = solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
				rhoE_buffer[i][t * buffer_num + j] = solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
			}
		}
	
		int dest = i < solver.myid ? i : i + 1;
		
		MPI_Isend(rho_buffer[i], buffer_num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 1 + 13 * commu_count, MPI_COMM_WORLD, &req);
		MPI_Isend(rhou_buffer[i], buffer_num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 2 + 13 * commu_count, MPI_COMM_WORLD, &req);
		MPI_Isend(rhov_buffer[i], buffer_num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 3 + 13 * commu_count, MPI_COMM_WORLD, &req);
		MPI_Isend(rhoE_buffer[i], buffer_num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 4 + 13 * commu_count, MPI_COMM_WORLD, &req);

		MPI_Irecv(rho_buffer[i], buffer_num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 1 + 13 * commu_count, MPI_COMM_WORLD, &request[0]);
		MPI_Irecv(rhou_buffer[i], buffer_num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 2 + 13 * commu_count, MPI_COMM_WORLD, &request[1]);
		MPI_Irecv(rhov_buffer[i], buffer_num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 3 + 13 * commu_count, MPI_COMM_WORLD, &request[2]);
		MPI_Irecv(rhoE_buffer[i], buffer_num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 4 + 13 * commu_count, MPI_COMM_WORLD, &request[3]);			
		//MPI_Waitall(4, request, status);


		//......................................................................................
		for (int t = 0; t < solver.grid.local_innerBoundary_index[i].size(); t++)
		{
			for (int p = 0; p < TRIANGLE_EDGES * EDGE_GPOINTS; p++) {
				convar_rho_edge_buffer[i][t + p * buffer_num] = convar_rho_edge[solver.grid.local_innerBoundary_index[i].at(t) + p * num];
				convar_rhou_edge_buffer[i][t + p * buffer_num] = convar_rhou_edge[solver.grid.local_innerBoundary_index[i].at(t) + p * num];
				convar_rhov_edge_buffer[i][t + p * buffer_num] = convar_rhov_edge[solver.grid.local_innerBoundary_index[i].at(t) + p * num];
				convar_rhoE_edge_buffer[i][t + p * buffer_num] = convar_rhoE_edge[solver.grid.local_innerBoundary_index[i].at(t) + p * num];

			}
		}
		MPI_Isend(convar_rho_edge_buffer[i], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 13 * commu_count, MPI_COMM_WORLD, &request1);
		MPI_Irecv(convar_rho_edge_buffer[i], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 13 * commu_count, MPI_COMM_WORLD, &request[4]);
		MPI_Isend(convar_rhou_edge_buffer[i], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 13 * commu_count, MPI_COMM_WORLD, &request1);
		MPI_Irecv(convar_rhou_edge_buffer[i], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 13 * commu_count, MPI_COMM_WORLD, &request[5]);
		MPI_Isend(convar_rhov_edge_buffer[i], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 13 * commu_count, MPI_COMM_WORLD, &request1);
		MPI_Irecv(convar_rhov_edge_buffer[i], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 13 * commu_count, MPI_COMM_WORLD, &request[6]);
		MPI_Isend(convar_rhoE_edge_buffer[i], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 13 * commu_count, MPI_COMM_WORLD, &request1);
		MPI_Irecv(convar_rhoE_edge_buffer[i], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 13 * commu_count, MPI_COMM_WORLD, &request[7]);

		MPI_Waitall(8, request, status);
	}
	//MPI_Waitall(8 * (solver.nprocs - 1), request, status);
	dealCommuData();
}

void RkdgAdvance::dealCommuData() {
	int buffer_num(100);
	size_t pitch = solver._cuarrays.getDoublePitch();

	int tvn = solver.grid.getLocalTriangleNumber();
	int cell = solver.grid.getLocalCellNumber();
	for (int i = 0; i < solver.nprocs - 1; i++) {
		int cur_loc(solver.grid.getLocalTriangleNumber());
		for (int t = 0; t < i; t++) {
			cur_loc += solver.grid.local_innerBoundary_index[t].size();
		}
		int num = solver.grid.local_innerBoundary_index[i].size();
		for (int t = 0; t < BASIS_FUNCTIONS; t++) {
			for (int j = 0; j < num; j++) {
				solver._freedom_rho[t * solver.grid.getLocalCellNumber() + cur_loc + j] =  rho_buffer[i][t * buffer_num + j];
				solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhou_buffer[i][t * buffer_num + j];
				solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhov_buffer[i][t * buffer_num + j];
				solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhoE_buffer[i][t * buffer_num + j];
			}
		}


		for (int j = 0; j < solver.grid.local_innerBoundary_index[i].size(); j++) {
			for (int t = 0; t < TRIANGLE_EDGES * EDGE_GPOINTS; t++) {
				convar_rho_edge[tvn + t * cell] = 
					convar_rho_edge_buffer[i][j + t * buffer_num];
				convar_rhou_edge[tvn + t * cell] = 
					convar_rhou_edge_buffer[i][j + t * buffer_num];
				convar_rhov_edge[tvn + t * cell] = 
					convar_rhov_edge_buffer[i][j + t * buffer_num];
				convar_rhoE_edge[tvn + t * cell] = 
					convar_rhoE_edge_buffer[i][j + t * buffer_num];
			}
			tvn++;
		}
	}
	size_t size = sizeof(double)*solver.grid.getLocalCellNumber();
	cudaMemcpy2DAsync(solver._cuarrays.freedom_rho,  pitch, solver._freedom_rho,  size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			  
	cudaMemcpy2DAsync(solver._cuarrays.freedom_rhou, pitch, solver._freedom_rhou, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			  
	cudaMemcpy2DAsync(solver._cuarrays.freedom_rhov, pitch, solver._freedom_rhov, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			 
	cudaMemcpy2DAsync(solver._cuarrays.freedom_rhoE, pitch, solver._freedom_rhoE, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(solver._cuarrays.convar_rho_edge, pitch, convar_rho_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);			  
	cudaMemcpy2DAsync(solver._cuarrays.convar_rhou_edge, pitch, convar_rhou_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);			 
	cudaMemcpy2DAsync(solver._cuarrays.convar_rhov_edge, pitch, convar_rhov_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);
	cudaMemcpy2DAsync(solver._cuarrays.convar_rhoE_edge, pitch, convar_rhoE_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);
	cudaError_t error = cudaPeekAtLastError();
	if ( error!=cudaSuccess ) {
		throw CMyException(cudaGetErrorString(error));
	}
	cudaDeviceSynchronize();
}

RkdgAdvance::~RkdgAdvance() {
	for (int i = 0; i <solver.nprocs - 1; i++) {
		delete[] rho_buffer[i];
		delete[] rhou_buffer[i];
		delete[] rhov_buffer[i];
		delete[] rhoE_buffer[i];
		delete[] convar_rho_edge_buffer[i];
		delete[] convar_rhou_edge_buffer[i];
		delete[] convar_rhov_edge_buffer[i];
		delete[] convar_rhoE_edge_buffer[i];
	}
	delete[] rho_buffer;
	delete[] rhou_buffer;
	delete[] rhov_buffer;
	delete[] rhoE_buffer;

	delete[] convar_rho_edge;
	delete[] convar_rhou_edge;
	delete[] convar_rhov_edge;
	delete[] convar_rhoE_edge;
}