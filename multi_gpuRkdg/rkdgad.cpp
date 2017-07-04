//#include <mpi.h>
//#include "../inc/rkdgadvance.h"
//
//RkdgAdvance::RkdgAdvance(CCUDARkdgSolver* rkdg_solver) {
//	solver = *rkdg_solver;
//}
//
//void RkdgAdvance::advance() {
//	ofstream fout;
//	if ( solver.log_history=='Y' )
//	{
//		fout.open(solver.residual_file.c_str());
//		if ( !fout )
//			throw CMyException("Failed to open residual log file: "+solver.residual_file);
//		
//		fout<<"N, rho"<<endl;
//	}
//			
//			
//	double nt(0);
//	int count(0);
//		
//	/*int tnum = grid.getTriangleNumber();
//	int num  = grid.getCellNumber();*/
//	int tnum = solver.grid.getLocalTriangleNumber();
//	int num = solver.grid.getLocalCellNumber();
//	int* innerBoudary_count;
//	innerBoudary_count = (int*)malloc((solver.nprocs - 1) * sizeof(int));
//	for (int i = 0; i < solver.nprocs - 1; i++)
//	{
//		innerBoudary_count[i] = solver.grid.local_innerBoundary_index[i].size();
//	}
//		
//	int blocks = (tnum%solver.threads_per_block) ? tnum/solver.threads_per_block+1 : tnum/solver.threads_per_block;
//		
//	double ut   = sqrt(solver.gamma*solver.pref/solver.rhoref)*solver.mach;
//	double rhou = solver.rhoref*ut*cos(solver.alpha);
//	double rhov = solver.rhoref*ut*sin(solver.alpha);
//	double rhoE = 0.5*solver.rhoref*(ut*ut) + solver.pref/(solver.gamma-1);
//		
//	bool copy(false);
//		
//	cudaError_t error;
//	size_t pitch = solver._cuarrays.getDoublePitch();
//	int pitch_num = pitch / sizeof(double);
//	int ipitch_num = solver._cuarrays.getIntPitch() / sizeof(int);
//		
//	cudaEvent_t time_start, time_stop;
//		
//	cudaEventCreateWithFlags(&time_start, cudaEventDisableTiming|cudaEventBlockingSync);
//	cudaEventCreateWithFlags(&time_stop,  cudaEventDisableTiming|cudaEventBlockingSync);
//		
//	if ( solver.log_history=='Y' )
//		copy = true;
//		
//	// 确保之前CUDA的初始化工作都已经完成
//	cudaDeviceSynchronize();
//	int commu_count(0);	
//	int continued(0);
//	int depart_num  = solver.grid.local_innerBoundary_index[0].size();
//	
//	convar_rho_edge = (double*)malloc(num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//	convar_rhou_edge = (double*)malloc(num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//	convar_rhov_edge = (double*)malloc(num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//	convar_rhoE_edge = (double*)malloc(num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//
//	convar_rho_edge_buffer = new double*[solver.nprocs - 1];
//	convar_rhou_edge_buffer = new double*[solver.nprocs - 1];
//	convar_rhov_edge_buffer = new double*[solver.nprocs - 1];
//	convar_rhoE_edge_buffer = new double*[solver.nprocs - 1];
//	int buffer_num(0);
//	for (int i = 0; i < solver.nprocs - 1; i++) {
//		if (buffer_num < innerBoudary_count[i]) {
//			buffer_num = innerBoudary_count[i];
//		}
//	}
//	cout<<"a1"<<endl;
//	for (int i = 0; i < solver.nprocs - 1; i++) {
//		convar_rho_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//		convar_rhou_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//		convar_rhov_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//		convar_rhoE_edge_buffer[i] = (double*)malloc(buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//	}
//	cout<<"b1"<<endl;
//	do 
//	{
//		++ count;
//		isCommun = false;
//		cudaEventRecord(time_start);
//		// 计算当前时间步长
//		solver.getTimeStep(tnum);
//		
//		cudaEventRecord(time_stop);
//		
//		// 保存旧自由度
//		cudaMemcpy2DAsync(solver._cuarrays.freedom_rho_old,  pitch, solver._cuarrays.freedom_rho,  pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
//		cudaMemcpy2DAsync(solver._cuarrays.freedom_rhou_old, pitch, solver._cuarrays.freedom_rhou, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
//		cudaMemcpy2DAsync(solver._cuarrays.freedom_rhov_old, pitch, solver._cuarrays.freedom_rhov, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							 
//		cudaMemcpy2DAsync(solver._cuarrays.freedom_rhoE_old, pitch, solver._cuarrays.freedom_rhoE, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);
//		error = cudaPeekAtLastError();
//		if ( error!=cudaSuccess ) {
//			throw CMyException(cudaGetErrorString(error));
//		}
//		size_t size = sizeof(double) * num;
//		//cudaMemcpy2DAsync(solver._freedom_rho, host_pitch, solver._cuarrays.freedom_rho, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
//		cudaMemcpy2D(solver._freedom_rho,  size, solver._cuarrays.freedom_rho,  pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			  
//		
//		error = cudaPeekAtLastError();
//		if ( error!=cudaSuccess ) {
//			throw CMyException(cudaGetErrorString(error));
//		}
//		cudaMemcpy2D(solver._freedom_rhou, size, solver._cuarrays.freedom_rhou, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			  
//		cudaMemcpy2D(solver._freedom_rhov, size, solver._cuarrays.freedom_rhov, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			 
//		cudaMemcpy2D(solver._freedom_rhoE, size, solver._cuarrays.freedom_rhoE, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
//		cudaDeviceSynchronize();
//		//for(int t = 0; t < 10; t++) {
//		//	cout<<"myid:"<<solver.myid<<",previous rhou:"<<solver._freedom_rhou[t]<<endl;
//		//	/*cout<<"myid:"<<solver.myid<<",previous rho:"<<solver._freedom_rho[t]<<endl;
//		//	cout<<"myid:"<<solver.myid<<",previous rhov:"<<solver._freedom_rhov[t]<<endl;
//		//	cout<<"myid:"<<solver.myid<<",previous rhoE:"<<solver._freedom_rhoE[t]<<endl;*/
//		//}
//		/*cudaMemcpy2DAsync(solver._freedom_rho, host_pitch, solver._cuarrays.freedom_rho, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
//		cudaMemcpy2DAsync(solver._freedom_rho, host_pitch, solver._cuarrays.freedom_rho, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
//		cudaMemcpy2DAsync(solver._freedom_rho, host_pitch, solver._cuarrays.freedom_rho, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);*/
//		MPI_Request *request, request1, request2;
//		MPI_Status status;
//		cout<<"c1"<<endl;
//		//int x = 1;         //通知其他节点信息交换，0表示结束
//		//if(count == 1 || (solver.myid == 0 && count % 5 == 0)) {
//		//	for (int i = 1; i < solver.nprocs; i++) {
//		//		MPI_Isend(&x, 1, MPI_INT, i, solver.nprocs + 9 * commu_count, MPI_COMM_WORLD, &request2);
//		//	}
//		//	commu_count++;
//		//	commuInfo(commu_count);
//		//} else if(solver.myid != 0) {
//		//	int flag(0);
//		//	MPI_Status status;
//		//	if (count == 1) {
//		//		MPI_Irecv(&x, 1, MPI_INT, 0, solver.nprocs + 9 * commu_count++, MPI_COMM_WORLD, &request2);
//		//	} 
//		//	MPI_Test(&request2, &flag, &status);
//		//	if (flag == 1 && x == 1) {
//		//		commuInfo(commu_count);
//		//		MPI_Irecv(&x,1, MPI_INT, 0, solver.nprocs + 9 * commu_count++, MPI_COMM_WORLD, &request2);
//		//	} else if (x == 0) {
//		//		break;
//		//	}
//		//}
//		commuInfo(commu_count++);
//		cout<<"v1"<<endl;
//		for ( int i=0; i<RUNGE_KUTTA_STEPS; ++i )
//		{
//			// 计算守恒量的值
//			solver.calculateConVars(tnum, pitch_num, blocks);
//			//cout<<solver.myid<<"tnum:"<<tnum<<",pitch_num:"<<pitch_num<<",blocks:"<<blocks<<endl;
//			// 处理边界条件
//			int local_inner_count(0);
//			for (int t = 0; t < solver.nprocs - 1; t++)
//			{
//				local_inner_count += solver.grid.local_innerBoundary_index[t].size(); 
//			}
//			solver.boundaryCondition(tnum + local_inner_count, num, pitch_num, solver.rhoref, rhou, rhov, rhoE);
//			cout<<"f1"<<endl;
//			if (isCommun) {
//				//cout<<"convar 1"<<endl;
//				cudaMemcpy2D(convar_rho_edge, size, solver._cuarrays.convar_rho_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);			  
//				cudaMemcpy2D(convar_rhou_edge, size, solver._cuarrays.convar_rhou_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);			 
//				cudaMemcpy2D(convar_rhov_edge, size, solver._cuarrays.convar_rhov_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);
//				cudaMemcpy2D(convar_rhoE_edge, size, solver._cuarrays.convar_rhoE_edge, pitch, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyDeviceToHost);
//				cudaDeviceSynchronize();
//				for(int i = 0; i < 10; i++) {
//					cout<<"r  rho:"<<convar_rho_edge[i]<<endl;
//				}
//				//cout<<"myid:"<<solver.myid<<", convar 2"<<endl;
//				for (int j = 0; j < solver.nprocs - 1; j++) {
//					/*convar_rho_edge_buffer = (double*)malloc(solver.grid.local_innerBoundary_index[j].size() * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//					convar_rhou_edge_buffer = (double*)malloc(solver.grid.local_innerBoundary_index[j].size() * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//					convar_rhov_edge_buffer = (double*)malloc(solver.grid.local_innerBoundary_index[j].size() * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));
//					convar_rhoE_edge_buffer = (double*)malloc(solver.grid.local_innerBoundary_index[j].size() * TRIANGLE_EDGES * EDGE_GPOINTS * sizeof (double));*/
//					//cout<<"myid:"<<solver.myid<<", convar 2.5"<<endl;
//					request = (MPI_Request *) malloc(4 * (solver.nprocs - 1) * sizeof(MPI_Request));
//					for (int t = 0; t < innerBoudary_count[j]; t++)
//					{
//						for (int p = 0; p < TRIANGLE_EDGES * EDGE_GPOINTS; p++) {
//							convar_rho_edge_buffer[j][t + p * buffer_num] = convar_rho_edge[solver.grid.local_innerBoundary_index[j].at(t) + p * num];
//							convar_rhou_edge_buffer[j][t + p * buffer_num] = convar_rhou_edge[solver.grid.local_innerBoundary_index[j].at(t) + p * num];
//							convar_rhov_edge_buffer[j][t + p * buffer_num] = convar_rhov_edge[solver.grid.local_innerBoundary_index[j].at(t) + p * num];
//							convar_rhoE_edge_buffer[j][t + p * buffer_num] = convar_rhoE_edge[solver.grid.local_innerBoundary_index[j].at(t) + p * num];
//							
//						}
//					}
//					for(int i = 0; i < 10; i++) {
//						cout<<"rhou:"<<convar_rho_edge[i]<<endl;
//					}
//					for(int i = 0; i < 10; i++) {
//						cout<<"receive rhou:"<<convar_rho_edge_buffer[0][i]<<endl;
//					}
//					//cout<<"myid:"<<solver.myid<<", convar3"<<endl;
//					int dest = j < solver.myid ? j : j+1;
//					cout<<solver.myid<<" count:"<<innerBoudary_count[j]<<endl;
//					cout<<"nprocs:"<<solver.nprocs<<",j:"<<j<<",commu_count:"<<commu_count<<",myid:"<<solver.myid<<endl;
//					MPI_Isend(convar_rho_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 9 * commu_count, MPI_COMM_WORLD, &request1);
//					MPI_Isend(convar_rhou_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 9 * commu_count, MPI_COMM_WORLD, &request1);
//					MPI_Isend(convar_rhov_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 9 * commu_count, MPI_COMM_WORLD, &request1);
//					MPI_Isend(convar_rhoE_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 9 * commu_count, MPI_COMM_WORLD, &request1);					
//					cout<<"s1"<<endl;
//					MPI_Irecv(convar_rho_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 5 + 9 * commu_count, MPI_COMM_WORLD, &request[j * 4]);	
//					cout<<"m1"<<endl;
//					MPI_Irecv(convar_rhou_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 6 + 9 * commu_count, MPI_COMM_WORLD, &request[j * 4 + 1]);
//					cout<<"n1"<<endl;
//					MPI_Irecv(convar_rhov_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 7 + 9 * commu_count, MPI_COMM_WORLD, &request[j * 4 + 2]);
//					cout<<"l1"<<endl;
//					MPI_Irecv(convar_rhoE_edge_buffer[j], buffer_num * TRIANGLE_EDGES * EDGE_GPOINTS, MPI_DOUBLE, dest, solver.nprocs + 8 + 9 * commu_count++, MPI_COMM_WORLD, &request[j * 4 +3]);
//					cout<<"d1"<<endl;
//					//cout<<"myid:"<<solver.myid<<",convar 4"<<endl;
//				}
//				cout<<"u1"<<endl;
//				MPI_Waitall(4 * (solver.nprocs - 1), request, &status);
//				cout<<"h1"<<endl;
//				for(int i = 0; i < 10; i++) {
//					cout<<"receive convar rhou:"<<convar_rho_edge_buffer[0][i]<<endl;
//				}
//				dealCommuConvar(buffer_num, convar_rho_edge_buffer, convar_rhou_edge_buffer, convar_rhov_edge_buffer, convar_rhoE_edge_buffer);
//				commu_count++;
//			}
//			cout<<"d1"<<endl;
//			// 计算体积分残差
//			solver.calculateVolumeRHS(tnum, pitch_num, blocks);
//
//			// 计算LF通量系数
//			solver.calculateLFCoeff(tnum, ipitch_num, pitch_num, blocks);
//		
//			// 计算f, g在边上的值
//			solver.calculateEdgeFG(tnum, num, pitch_num, blocks);
//		
//			solver.calculateFlux(tnum, ipitch_num, pitch_num, blocks);
//
//			// 计算线积分残差
//			solver.calculateEdgeRHS(tnum, pitch_num, blocks);
//
//			// 时间推进
//			switch (i)
//			{
//			case 0:
//				cudaEventSynchronize(time_stop);
//		
//				// 将时间步长传送到本地
//				cudaMemcpy(solver._dt, solver._cuarrays.ddt, sizeof(double), cudaMemcpyDeviceToHost);
//		
//				if ( 0==(count-1)%solver.print_interval )
//					cout<<"Step: "<<count<<", time step: "<<solver._dt[0]<<endl;
//		
//				if ( (solver._terminal_time-nt)<solver._dt[0] )
//				{
//					solver._dt[0] = solver._terminal_time -  nt;
//				}
//		
//				// 时间步推进
//				solver.rkdgStepOne(solver._dt[0], tnum, pitch_num, blocks);
//		
//				break;
//		
//			case 1:
//				solver.rkdgStepTwo(solver._dt[0], tnum, pitch_num, blocks);
//				break;
//		
//			case 2:
//				solver.rkdgStepThree(solver._dt[0], tnum, pitch_num, blocks);
//				break;
//		
//			default:
//				throw CMyException("impossible case!");
//				break;
//			}
//		}
//				
//		if ( copy && (count-1) )
//		{
//			// 复制残差数据
//			cudaMemcpy(solver._residual, solver._cuarrays.residual,
//						sizeof(double)*RESIDUAL_VARS, cudaMemcpyDeviceToHost);
//		
//			if ( 0==(count-1)%solver.print_interval )
//				cout<<"Current time: "<<nt<<"   rhomax: "<<solver._residual[0]/solver.rhoref<<"   E: "<<solver._residual[1]/rhoE<<endl;
//		
//			fout<<count<<"   "<<log(solver._residual[0]/solver.rhoref)/log(10.0)<<endl;
//		}
//		// 计算残差
//		solver.calculateResidual(tnum);
//		//cout<<"calculateResidual"<<endl;
//				
//		// 计时推进
//		nt += solver._dt[0];
//		
//		error = cudaPeekAtLastError();
//		//cout<<"error test"<<endl;
//		if ( error!=cudaSuccess ) {
//			cout<<"do error"<<endl;
//			throw CMyException(cudaGetErrorString(error));
//		}
//		//cout<<solver.myid<<"no error"<<endl;
//		/*int end(0);
//		if (nt<solver._terminal_time)
//		{
//		end = 1;
//		}
//
//		MPI_Reduce(&end, &continued, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//		cout<<"reduce end"<<endl;
//		MPI_Bcast(&continued, 1, MPI_INT, 0, MPI_COMM_WORLD);
//		cout<<"continued:"<<continued<<endl;*/
//		/*if (solver.myid == 0 && nt >= 0.02) {
//			x = 0;
//			for (int i = 1; i < solver.nprocs; i++) {
//				MPI_Isend(&x, 1, MPI_INT, i, solver.nprocs + 9 * commu_count++, MPI_COMM_WORLD, &request2);
//			}
//			break;
//		}*/
//	} while ( count < 5 );
//	cout<<solver.myid<<" process end while, nt"<<nt<<", terminal:"<<solver._terminal_time<<endl;	
//	cudaDeviceSynchronize();
//			
//	if ( copy )
//	{
//		// 复制残差数据
//		cudaMemcpy(solver._residual, solver._cuarrays.residual,
//			sizeof(double)*RESIDUAL_VARS, cudaMemcpyDeviceToHost);
//		if ( 0==(count-1)%solver.print_interval )
//			cout<<"当前时间： "<<nt-solver._dt[0]<<"   rhomax: "<<solver._residual[0]/solver.rhoref<<"   E: "<<solver._residual[1]/rhoE<<endl;
//		
//		fout<<count<<"   "<<log(solver._residual[0]/solver.rhoref)/log(10.0)<<endl;
//	}	
//		
//	cudaEventDestroy(time_start);
//	cudaEventDestroy(time_stop);
//		
//	if ( solver.log_history=='Y' )
//		fout.close();
//}
//
//
//void RkdgAdvance::commuInfo(int commu_count) {
//	isCommun = true;
//	MPI_Request request[4],req;
//	MPI_Status status[4];	
//	for (int i = 0; i < solver.nprocs - 1; i++) {
//		int num = solver.grid.local_innerBoundary_index[i].size();
//		rho_buffer = (double*)malloc(num * BASIS_FUNCTIONS * sizeof(double));
//		rhou_buffer = (double*)malloc(num * BASIS_FUNCTIONS * sizeof(double));
//		rhov_buffer = (double*)malloc(num * BASIS_FUNCTIONS * sizeof(double));
//		rhoE_buffer = (double*)malloc(num * BASIS_FUNCTIONS * sizeof(double));
//		for (int t = 0; t < BASIS_FUNCTIONS; t++) {
//			for (int j = 0; j < num; j++) {
//				rho_buffer[t * num + j] = solver._freedom_rho[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
//				rhou_buffer[t * num + j] = solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
//				rhov_buffer[t * num + j] = solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
//				rhoE_buffer[t * num + j] = solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + solver.grid.local_innerBoundary_index[i].at(j)];
//			}
//		}
//		/*for(int i = 0; i < 10; i++) {
//			cout<<"index:"<<solver.grid.local_innerBoundary_index[0].at(i)<<",freedom rhou:"<<solver._freedom_rhou[solver.grid.local_innerBoundary_index[0].at(i)]<<endl;
//		}
//		for(int t = 0; t < 10; t++) {
//			cout<<"buffer rhou:"<<rhou_buffer[t]<<endl;
//		}*/
//		int dest = i < solver.myid ? i : i + 1;
//		
//		MPI_Isend(rho_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 1 + 9 * commu_count, MPI_COMM_WORLD, &req);
//		MPI_Isend(rhou_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 2 + 9 * commu_count, MPI_COMM_WORLD, &req);
//		MPI_Isend(rhov_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 3 + 9 * commu_count, MPI_COMM_WORLD, &req);
//		MPI_Isend(rhoE_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 4 + 9 * commu_count, MPI_COMM_WORLD, &req);
//
//		MPI_Irecv(rho_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 1 + 9 * commu_count, MPI_COMM_WORLD, &request[0]);
//		MPI_Irecv(rhou_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 2 + 9 * commu_count, MPI_COMM_WORLD, &request[1]);
//		MPI_Irecv(rhov_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 3 + 9 * commu_count, MPI_COMM_WORLD, &request[2]);
//		MPI_Irecv(rhoE_buffer, num * BASIS_FUNCTIONS, MPI_DOUBLE, dest, solver.nprocs + 4 + 9 * commu_count, MPI_COMM_WORLD, &request[3]);			
//		MPI_Waitall(4, request, status);
//		/*for(int i = 0; i < 10; i++) {
//			cout<<"receive rhou:"<<rhou_buffer[i]<<endl;
//		}*/
//		dealCommuData(i);
//	} 
//	//cout<<"wait end"<<endl;
//	//cout<<"myid2:"<<solver.myid<<"  ";
//	//for(int m = 0; m < 5; m++) {
//	//	cout<<rho_buffer[m]<<"  ";
//	//}	
//	//cout<<endl;
//	
//}
//
//void RkdgAdvance::dealCommuData(int i) {
//	size_t pitch = solver._cuarrays.getDoublePitch();
//	int cur_loc(solver.grid.getLocalTriangleNumber());
////	for (int i = 0; i < solver.nprocs - 1; i++) {
//	for (int t = 0; t < i; t++) {
//		cur_loc += solver.grid.local_innerBoundary_index[t].size();
//	}
//	int num = solver.grid.local_innerBoundary_index[i].size();
//	for (int t = 0; t < BASIS_FUNCTIONS; t++) {
//		for (int j = 0; j < num; j++) {
//			/*solver._freedom_rho[t * solver.grid.getLocalCellNumber() + cur_loc + j] =  rho_buffer[t * num + j];
//			solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhou_buffer[t * num + j];
//			solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhov_buffer[t * num + j];
//			solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhoE_buffer[t * num + j];*/
//
//			solver._freedom_rho[t * solver.grid.getLocalCellNumber() + cur_loc + j] =  rho_buffer[t * num + j];
//			solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhou_buffer[t * num + j];
//			solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhov_buffer[t * num + j];
//			solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + cur_loc + j] = rhoE_buffer[t * num + j];
//		}
//	}
////	}
//	size_t size = sizeof(double)*solver.grid.getLocalCellNumber();
//	//cudaMemcpy2DAsync(solver._cuarrays.freedom_rho, pitch, solver._freedom_rho, host_pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
//	//cudaMemcpy2DAsync(solver._cuarrays.freedom_rhou, pitch, solver._freedom_rhou, host_pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
//	//cudaMemcpy2DAsync(solver._cuarrays.freedom_rhov, pitch, solver._freedom_rhov, host_pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
//	//cudaMemcpy2DAsync(solver._cuarrays.freedom_rhoE, pitch, solver._freedom_rhoE, host_pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
//	cudaMemcpy2D(solver._cuarrays.freedom_rho,  pitch, solver._freedom_rho,  size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			  
//	cudaMemcpy2D(solver._cuarrays.freedom_rhou, pitch, solver._freedom_rhou, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			  
//	cudaMemcpy2D(solver._cuarrays.freedom_rhov, pitch, solver._freedom_rhov, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			 
//	cudaMemcpy2D(solver._cuarrays.freedom_rhoE, pitch, solver._freedom_rhoE, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);
//	cudaError_t error = cudaPeekAtLastError();
//	if ( error!=cudaSuccess ) {
//		throw CMyException(cudaGetErrorString(error));
//	}
//}
//
//void RkdgAdvance::dealCommuConvar(int buffer_num, double **convar_rho_edge_buffer, double **convar_rhou_edge_buffer, double **convar_rhov_edge_buffer, double **convar_rhoE_edge_buffer) {
//	int cell = solver.grid.getLocalCellNumber();
//	int start_loc = solver.grid.getLocalTriangleNumber();
//	for (int j = 0; j < solver.nprocs - 1; j++) {
//		for (int i = 0; i < solver.grid.local_innerBoundary_index[j].size(); i++) {
//			for (int t = 0; t < TRIANGLE_EDGES * EDGE_GPOINTS; t++) {
//				convar_rho_edge[start_loc + i + t * cell] = 
//					convar_rho_edge_buffer[j][i + t * buffer_num];
//				convar_rhou_edge[start_loc + i + t * cell] = 
//					convar_rhou_edge_buffer[j][i + t * buffer_num];
//				convar_rhov_edge[start_loc + i + t * cell] = 
//					convar_rhov_edge_buffer[j][i + t * buffer_num];
//				convar_rhoE_edge[start_loc + i + t * cell] = 
//					convar_rhoE_edge_buffer[j][i + t * buffer_num];
//				//cout<<"commu convar:"<<convar_rhou_edge_buffer[i + t * TRIANGLE_EDGES * EDGE_GPOINTS]<<endl;
//			}
//		}
//		start_loc += solver.grid.local_innerBoundary_index[j].size();
//	}
//	
//	
//	/*delete []convar_rho_edge_buffer;
//	delete []convar_rho_edge_buffer;
//	delete []convar_rho_edge_buffer;
//	delete []convar_rho_edge_buffer;*/
//
//	size_t size = sizeof(double) * solver.grid.getLocalCellNumber();
//	size_t pitch = solver._cuarrays.getDoublePitch();
//	cudaMemcpy2D(solver._cuarrays.convar_rho_edge, pitch, convar_rho_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);			  
//	cudaMemcpy2D(solver._cuarrays.convar_rhou_edge, pitch, convar_rhou_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);			 
//	cudaMemcpy2D(solver._cuarrays.convar_rhov_edge, pitch, convar_rhov_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);
//	cudaMemcpy2D(solver._cuarrays.convar_rhoE_edge, pitch, convar_rhoE_edge, size, size, TRIANGLE_EDGES * EDGE_GPOINTS, cudaMemcpyHostToDevice);
//	cudaDeviceSynchronize();
//}
//*/