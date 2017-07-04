#include <mpi.h>
#include "../inc/rkdgadvance.h"

RkdgAdvance::RkdgAdvance(CCUDARkdgSolver* rkdg_solver) {
	solver = *rkdg_solver;
}

void RkdgAdvance::advance() {
	clock_t start, end;
	double commu_time(0);
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
		
	int tnum = solver.grid.getLocalTriangleNumber();
	int num = solver.grid.getLocalCellNumber();
	int* innerBoudary_count;
	innerBoudary_count = (int*)malloc((solver.nprocs - 1) * sizeof(int));
	for (int i = 0; i < solver.nprocs - 1; i++)
	{
		innerBoudary_count[i] = solver.grid.local_innerBoundary_index[i].size();
	}
		
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

	convar_rho_edge_buffer = new double*[solver.nprocs - 1];
	convar_rhou_edge_buffer = new double*[solver.nprocs - 1];
	convar_rhov_edge_buffer = new double*[solver.nprocs - 1];
	convar_rhoE_edge_buffer = new double*[solver.nprocs - 1];
	int buffer_num(500);
	/*for (int i = 0; i < solver.nprocs - 1; i++) {
		if (buffer_num < innerBoudary_count[i]) {
			buffer_num = innerBoudary_count[i];
		}
	}*/
	//cout<<"a1"<<endl;
	for (int i = 0; i < solver.nprocs - 1; i++) {
		convar_rho_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
		convar_rhou_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
		convar_rhov_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
		convar_rhoE_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
	}
	//cout<<"b1"<<endl;

	MPI_Request *request, request1, request2;
	MPI_Status status;
	int x = 1;    //通知其他节点信息交换，0表示结束
	do 
	{
		++ count;
		isCommun = false;
		cudaEventRecord(time_start);
		// 计算当前时间步长
		solver.getTimeStep(tnum);
		
		cudaEventRecord(time_stop);
		
		// 保存旧自由度
		cudaMemcpy2DAsync(solver._cuarrays.freedom_rho_old,  pitch, solver._cuarrays.freedom_rho,  pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(solver._cuarrays.freedom_rhou_old, pitch, solver._cuarrays.freedom_rhou, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(solver._cuarrays.freedom_rhov_old, pitch, solver._cuarrays.freedom_rhov, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							 
		cudaMemcpy2DAsync(solver._cuarrays.freedom_rhoE_old, pitch, solver._cuarrays.freedom_rhoE, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);
		
		size_t size = sizeof(double) * num;
		cudaMemcpy2DAsync(solver._freedom_rho,  size, solver._cuarrays.freedom_rho,  pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			  
		cudaMemcpy2DAsync(solver._freedom_rhou, size, solver._cuarrays.freedom_rhou, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			  
		cudaMemcpy2DAsync(solver._freedom_rhov, size, solver._cuarrays.freedom_rhov, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			 
		cudaMemcpy2DAsync(solver._freedom_rhoE, size, solver._cuarrays.freedom_rhoE, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		
		if (count == 1) {
			commuInfo(commu_count);
		}else if (solver.myid == 0 && count % 5 == 0) {
			for (int i = 1; i < solver.nprocs; i++) {
				MPI_Isend(&x, 1, MPI_INT, i, solver.nprocs + 9 * commu_count, MPI_COMM_WORLD, &request2);
			}
			commu_count++;
			start = clock();
			commuInfo(commu_count);
			end = clock();
			commu_time += (double)(end - start) / CLOCKS_PER_SEC;
			//isCommun = true;
		} else if(solver.myid != 0) {
			int flag(0);
			if (count == 1) {
				MPI_Irecv(&x, 1, MPI_INT, 0, solver.nprocs + 9 * commu_count++, MPI_COMM_WORLD, &request2);
			}
			MPI_Test(&request2, &flag, &status);
			if (flag == 1 && x == 1) {
				//MPI_Request_free(&request2);
				start = clock();
				commuInfo(commu_count);
				end = clock();
				commu_time += (double)(end - start) / CLOCKS_PER_SEC;
				//isCommun = true;
				//MPI_Irecv(&x,1, MPI_INT, 0, solver.nprocs + 9 * commu_count++, MPI_COMM_WORLD, &request2);
			} else if (flag == 1 && x == 0) {
				break;
			}
		}

		/*if (isCommun) {
			cout<<"if co:"<<commu_count<<endl;
			commuInfo(commu_count);
		}*/
		//start = clock();
		//commuInfo(commu_count++);
		//end = clock();
		//commu_time += (double)(end - start) / CLOCKS_PER_SEC;
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
				start = clock();
				commuConvar(size, pitch, buffer_num, commu_count, innerBoudary_count);
				end = clock();
				commu_time += (double)(end - start) / CLOCKS_PER_SEC;
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
		
				if ( (solver._terminal_time-nt)<solver._dt[0] )
				{
					solver._dt[0] = solver._terminal_time -  nt;
				}
		
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

		if (isCommun && solver.myid != 0) {
			MPI_Irecv(&x,1, MPI_INT, 0, solver.nprocs + 9 * commu_count++, MPI_COMM_WORLD, &request2);
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
				
		// 计时推进
		nt += solver._dt[0];
		
		error = cudaPeekAtLastError();
		if ( error!=cudaSuccess ) {
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
		if (solver.myid == 0 && count > 10000) {
			x = 0;
			for (int i = 1; i < solver.nprocs; i++) {
				MPI_Isend(&x, 1, MPI_INT, i, solver.nprocs + 9 * commu_count, MPI_COMM_WORLD, &request2);
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
	cout<<"total commu time:"<<commu_time<<endl;
	if ( solver.log_history=='Y' )
		fout.close();
}


void RkdgAdvance::commuInfo(int commu_count) {
	isCommun = true;
	MPI_Request *request,req;
	MPI_Status status[4];	
	request = (MPI_Request*)malloc(4 * (solver.nprocs - 1) * sizeof(MPI_Request));
	for (int i = 0; i < solver.nprocs - 1; i++) {
		int num = solver.grid.local_innerBoundary_index[i].size();
		rho_buffer = (double*)malloc(num * BASIS_FUNCTIONS * sizeof(double));
		rhou_buffer = (double*)malloc(num * BASIS_FUNCTIONS * sizeof(double));
		rhov_buffer = (double*)malloc(num * BASIS_FUNCTIONS * sizeof(double));
		rhoE_buffer = (double*)malloc(num * BASIS_FUNCTIONS * sizeof(double));
		for (int t = 0; t < BASIS_FUNCTIONS; t++) {
			for (int j = 0; j < num; j++) {
				rho_buffer[t * num + j] = solver._freedom_rho[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
				rhou_buffer[t * num + j] = solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
				rhov_buffer[t * num + j] = solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
				rhoE_buffer[t * num + j] = solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
			}
		}
		int dest = i < solver.myid ? i : i + 1;
		
		MPI_Isend(rho_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 1 + 9 * commu_count, MPI_COMM_WORLD, &req);
		MPI_Isend(rhou_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 2 + 9 * commu_count, MPI_COMM_WORLD, &req);
		MPI_Isend(rhov_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 3 + 9 * commu_count, MPI_COMM_WORLD, &req);
		MPI_Isend(rhoE_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 4 + 9 * commu_count, MPI_COMM_WORLD, &req);

		MPI_Irecv(rho_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 1 + 9 * commu_count, MPI_COMM_WORLD, &request[i * 4]);
		MPI_Irecv(rhou_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 2 + 9 * commu_count, MPI_COMM_WORLD, &request[i * 4 + 1]);
		MPI_Irecv(rhov_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 3 + 9 * commu_count, MPI_COMM_WORLD, &request[i * 4 + 2]);
		MPI_Irecv(rhoE_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 4 + 9 * commu_count, MPI_COMM_WORLD, &request[i * 4 + 3]);			

		MPI_Waitall(4, &request[i * 4], status);
		dealCommuData(i);
	} 
}

void RkdgAdvance::dealCommuData(int i) {
	size_t pitch = solver._cuarrays.getDoublePitch();
	int cur_loc(solver.grid.getLocalTriangleNumber());
	for (int t = 0; t < i; t++) {
		cur_loc += solver.grid.local_innerBoundary_index[t].size();
	}
	int num = solver.grid.local_innerBoundary_index[i].size();
	for (int t = 0; t < BASIS_FUNCTIONS; t++) {
		for (int j = 0; j < num; j++) {
			/*solver._freedom_rho[t * solver.grid.getLocalCellNumber() + cur_loc + j] =  rho_buffer[t * num + j];
			solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhou_buffer[t * num + j];
			solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhov_buffer[t * num + j];
			solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhoE_buffer[t * num + j];*/

			solver._freedom_rho[t * solver.grid.getLocalCellNumber() + cur_loc + j] =  rho_buffer[t * num + j];
			solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhou_buffer[t * num + j];
			solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhov_buffer[t * num + j];
			solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhoE_buffer[t * num + j];
		}
	}
	size_t size = sizeof(double)*solver.grid.getLocalCellNumber();
	cudaMemcpy2D(solver._cuarrays.freedom_rho,  pitch, solver._freedom_rho,  size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			  
	cudaMemcpy2D(solver._cuarrays.freedom_rhou, pitch, solver._freedom_rhou, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			  
	cudaMemcpy2D(solver._cuarrays.freedom_rhov, pitch, solver._freedom_rhov, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			 
	cudaMemcpy2D(solver._cuarrays.freedom_rhoE, pitch, solver._freedom_rhoE, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
	cudaError_t error = cudaPeekAtLastError();
	if ( error!=cudaSuccess ) {
		throw CMyException(cudaGetErrorString(error));
	}
}

void RkdgAdvance::commuConvar(int size, int pitch, int buffer_num, int commu_count, int *innerBoudary_count) {
		MPI_Request request, request1;
		MPI_Status status;
		int num = solver.grid.getLocalCellNumber();
		cudaMemcpy2DAsync(convar_rho_edge, size, solver._cuarrays.convar_rho_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);			  
		cudaMemcpy2DAsync(convar_rhou_edge, size, solver._cuarrays.convar_rhou_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);			 
		cudaMemcpy2DAsync(convar_rhov_edge, size, solver._cuarrays.convar_rhov_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);
		cudaMemcpy2DAsync(convar_rhoE_edge, size, solver._cuarrays.convar_rhoE_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		//request = (MPI_Request *) malloc(4 * (solver.nprocs - 1) * sizeof(MPI_Request));
		for (int j = 0; j < solver.nprocs - 1; j++) {
			//request = (MPI_Request *) malloc(4 * (solver.nprocs - 1) * sizeof(MPI_Request));
			for (int t = 0; t < innerBoudary_count[j]; t++)
			{
				for (int p = 0; p < TRIANGLE_EDGES * EDGE_GPOINTS; p++) {
					convar_rho_edge_buffer[j][t + p * buffer_num] = convar_rho_edge[solver.grid.local_innerBoundary_index[j].at(t) + p * num];
					convar_rhou_edge_buffer[j][t + p * buffer_num] = convar_rhou_edge[solver.grid.local_innerBoundary_index[j].at(t) + p * num];
					convar_rhov_edge_buffer[j][t + p * buffer_num] = convar_rhov_edge[solver.grid.local_innerBoundary_index[j].at(t) + p * num];
					convar_rhoE_edge_buffer[j][t + p * buffer_num] = convar_rhoE_edge[solver.grid.local_innerBoundary_index[j].at(t) + p * num];
							
				}
			}
			int dest = j < solver.myid ? j : j+1;
			/*cout<<solver.myid<<" count:"<<innerBoudary_count[j]<<endl;
			cout<<"nprocs:"<<solver.nprocs<<",j:"<<j<<",commu_count:"<<commu_count<<",myid:"<<solver.myid<<endl;
			MPI_Isend(convar_rho_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 9 * commu_count, MPI_COMM_WORLD, &request1);
			MPI_Isend(convar_rhou_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 9 * commu_count, MPI_COMM_WORLD, &request1);
			MPI_Isend(convar_rhov_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 9 * commu_count, MPI_COMM_WORLD, &request1);
			MPI_Isend(convar_rhoE_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 9 * commu_count, MPI_COMM_WORLD, &request1);					
			cout<<"s1"<<endl;
			MPI_Irecv(convar_rho_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 9 * commu_count, MPI_COMM_WORLD, &request[j * 4]);	
			cout<<"m1"<<endl;
			MPI_Irecv(convar_rhou_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 9 * commu_count, MPI_COMM_WORLD, &request[j * 4 + 1]);
			cout<<"n1"<<endl;
			MPI_Irecv(convar_rhov_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 9 * commu_count, MPI_COMM_WORLD, &request[j * 4 + 2]);
			cout<<"l1"<<endl;
			MPI_Irecv(convar_rhoE_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 9 * commu_count, MPI_COMM_WORLD, &request[j * 4 +3]);
			cout<<"d1"<<endl;*/
			//cout<<"myid:"<<solver.myid<<",convar 4"<<endl;


			if(dest > solver.myid) {
				MPI_Send(convar_rho_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 9 * commu_count, MPI_COMM_WORLD);
				MPI_Send(convar_rhou_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 9 * commu_count, MPI_COMM_WORLD);
				MPI_Send(convar_rhov_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 9 * commu_count, MPI_COMM_WORLD);
				MPI_Send(convar_rhoE_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 9 * commu_count, MPI_COMM_WORLD);					
				
				MPI_Recv(convar_rho_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 9 * commu_count, MPI_COMM_WORLD, &status);	
				MPI_Recv(convar_rhou_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 9 * commu_count, MPI_COMM_WORLD, &status);
				MPI_Recv(convar_rhov_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 9 * commu_count, MPI_COMM_WORLD, &status);
				MPI_Recv(convar_rhoE_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 9 * commu_count, MPI_COMM_WORLD, &status);
			} else {
				MPI_Recv(convar_rho_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 9 * commu_count, MPI_COMM_WORLD, &status);	
				MPI_Recv(convar_rhou_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 9 * commu_count, MPI_COMM_WORLD, &status);
				MPI_Recv(convar_rhov_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 9 * commu_count, MPI_COMM_WORLD, &status);
				MPI_Recv(convar_rhoE_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 9 * commu_count, MPI_COMM_WORLD, &status);

				MPI_Send(convar_rho_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 9 * commu_count, MPI_COMM_WORLD);
				MPI_Send(convar_rhou_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 9 * commu_count, MPI_COMM_WORLD);
				MPI_Send(convar_rhov_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 9 * commu_count, MPI_COMM_WORLD);
				MPI_Send(convar_rhoE_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 9 * commu_count, MPI_COMM_WORLD);					
			}
		}
		//MPI_Waitall(4 * (solver.nprocs - 1), request, &status);
		dealCommuConvar(buffer_num);
}

void RkdgAdvance::dealCommuConvar(int buffer_num) {
	int cell = solver.grid.getLocalCellNumber();
	int start_loc = solver.grid.getLocalTriangleNumber();
	for (int j = 0; j < solver.nprocs - 1; j++) {
		for (int i = 0; i < solver.grid.local_innerBoundary_index[j].size(); i++) {
			for (int t = 0; t < TRIANGLE_EDGES * EDGE_GPOINTS; t++) {
				convar_rho_edge[start_loc + i + t * cell] = 
					convar_rho_edge_buffer[j][i + t * buffer_num];
				convar_rhou_edge[start_loc + i + t * cell] = 
					convar_rhou_edge_buffer[j][i + t * buffer_num];
				convar_rhov_edge[start_loc + i + t * cell] = 
					convar_rhov_edge_buffer[j][i + t * buffer_num];
				convar_rhoE_edge[start_loc + i + t * cell] = 
					convar_rhoE_edge_buffer[j][i + t * buffer_num];
			}
		}
		start_loc += solver.grid.local_innerBoundary_index[j].size();
	}

	size_t size = sizeof(double) * solver.grid.getLocalCellNumber();
	size_t pitch = solver._cuarrays.getDoublePitch();
	cudaMemcpy2DAsync(solver._cuarrays.convar_rho_edge, pitch, convar_rho_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);			  
	cudaMemcpy2DAsync(solver._cuarrays.convar_rhou_edge, pitch, convar_rhou_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);			 
	cudaMemcpy2DAsync(solver._cuarrays.convar_rhov_edge, pitch, convar_rhov_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);
	cudaMemcpy2DAsync(solver._cuarrays.convar_rhoE_edge, pitch, convar_rhoE_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}
