/**
 * 在GPU上使用RKDG方法计算二维非结构网格上NACA0012机翼扰流，特点如下：
 * 1). 在每个三角单元上使用二次正交多项式作为基函数
 * 2). 三阶龙格库塔法进行时间推进
 *
 * 感谢： cjbuaa 提供算法和参考源程序
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */

#include<mpi.h>
#include "../inc/cppstdheaders.h"
#include "../inc/cudarkdgsolver.h"
#include"../inc/unstructuredgrid.h"
#include "../inc/rkdgadvance.h"
#include "../inc/sendboundaryindex.h"
using namespace std;


int main(int argc, char* argv[])
{
	CCUDARkdgSolver solver;
//	
//
	int myid;
	int nprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Request request, req[4];
	MPI_Status status;
	clock_t start, end;
//	CCUDARkdgSolver solver;
	//CUnstructuredGrid grid;
	//grid.config_file = "input/mesh.conf";
	//grid.initializeGrid(myid,2);

//	ofstream of("output/local.dat");
//	for (int i = 0; i < grid.elem_index.size() - 1; i++) {
////		cout<<grid.vertice.at(grid.elem_index[i]).getX()<<"  "<<grid.vertice.at(grid.elem_index[i]).getX()<<endl;
////		of<<grid.vertice.at(grid.elem_index[i]).getX()<<"  "<<grid.vertice.at(grid.elem_index[i]).getY()<<endl;
//		of<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3]).getX()<<"  "<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3]).getY()<<endl;
//		of<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3 + 1]).getX()<<"  "<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3 + 1]).getY()<<endl;
//		of<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3 + 2]).getX()<<"  "<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3 + 2]).getY()<<endl;
//	}
//	of.close();



	//try {
	start = clock();
		solver.config_file = "input/main.conf";

		solver.run(myid, nprocs);
		SendBoundary sendBoundary(&solver);
		sendBoundary.send();
		solver.runNext();
		RkdgAdvance advance(&solver);
		advance.advance();
		cout<<myid<<"end advance"<<endl;
		solver.runAfter();
		cout<<myid<<"end run after"<<endl;
	/*}
	catch ( const exception& e )
	{
		cout<<endl<<"Error occured, description: "<<endl;
		cout<<e.what()<<endl;
		cout<<endl<<"Please modify the program according to the hint."<<endl;
		getchar();

		exit(-1);
	}*/
	MPI_Barrier(MPI_COMM_WORLD);
	//double *rho_result, *rhou_result, *rhov_result, *rhoE_result;
	//double *rho_buffer, *rhou_buffer, *rhov_buffer, *rhoE_buffer;
	//int num = solver.grid.getTriangleNumber();
	//rho_result = (double*)malloc(BASIS_FUNCTIONS * solver.grid.getTriangleNumber() * sizeof(double));
	//rhou_result = (double*)malloc(BASIS_FUNCTIONS * solver.grid.getTriangleNumber() * sizeof(double));
	//rhov_result = (double*)malloc(BASIS_FUNCTIONS * solver.grid.getTriangleNumber() * sizeof(double));
	//rhoE_result = (double*)malloc(BASIS_FUNCTIONS * solver.grid.getTriangleNumber() * sizeof(double));
	
	/*MPI_Gatherv(solver._freedom_rho, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
		rho_result, solver.grid.getTriangleNumber(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(solver._freedom_rhou, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
		rhou_result, solver.grid.getTriangleNumber(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(solver._freedom_rhov, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
		rhov_result, solver.grid.getTriangleNumber(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(solver._freedom_rhoE, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
		rhoE_result, solver.grid.getTriangleNumber(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
*/
	cout<<"ggg"<<endl;

	if (myid != 0) {
		MPI_Send(solver._freedom_rho, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
			0, 10000,MPI_COMM_WORLD);
		MPI_Send(solver._freedom_rhou, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
			0, 10001,MPI_COMM_WORLD);
		MPI_Send(solver._freedom_rhov, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
			0, 10002,MPI_COMM_WORLD);
		MPI_Send(solver._freedom_rhoE, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
			0, 10003,MPI_COMM_WORLD);
	}
	int* data_size;
	data_size = (int*)malloc(nprocs * sizeof(int));
	data_size[0] = solver.grid.getLocalTriangleNumber();
	if (myid == 0) {
		int start_loc = solver.grid.getLocalTriangleNumber();
		for (int n = 1; n < nprocs; n++) {
			MPI_Probe(n,10000,MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_DOUBLE, &data_size[n]);

			//nprocs不大于6 
			MPI_Recv(solver._freedom_rho + start_loc, data_size[n], MPI_DOUBLE, n, 10000, MPI_COMM_WORLD, &status);
			MPI_Recv(solver._freedom_rhou + start_loc, data_size[n], MPI_DOUBLE, n, 10001, MPI_COMM_WORLD, &status);
			MPI_Recv(solver._freedom_rhov + start_loc, data_size[n], MPI_DOUBLE, n, 10002, MPI_COMM_WORLD, &status);
			MPI_Recv(solver._freedom_rhoE + start_loc, data_size[n], MPI_DOUBLE, n, 10003, MPI_COMM_WORLD, &status);
			start_loc += data_size[n];
		}
		cout<<"hhh"<<endl;
		//int i =0;
		//for (int t = 0; t < BASIS_FUNCTIONS; t++) {
		//	for(i = 0; i < solver.grid.getLocalTriangleNumber(); i++) {
		//		rho_result[t * num + i] = solver._freedom_rho[t * solver.grid.getLocalCellNumber() + i];
		//		rhou_result[t * num + i] = solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + i];
		//		rhov_result[t * num + i] = solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + i];
		//		rhoE_result[t * num + i] = solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + i];
		//	}
		//}
		///*for(int i = 0; i < 10; i++) {
		//cout<<"rhou:"<<rhou_result[i]<<endl;
		//}*/
		//MPI_Waitall(4, req, &status);
		//cout<<"sss"<<endl;
		//for (int t = 0; t < BASIS_FUNCTIONS; t++) {
		//	for (int j = solver.grid.getLocalTriangleNumber(); j < solver.grid.getTriangleNumber(); j++) {
		//		rho_result[t * num + j] = rho_buffer[t * data_size / BASIS_FUNCTIONS + j - solver.grid.getLocalTriangleNumber()];
		//		rhou_result[t * num + j] = rhou_buffer[t * data_size / BASIS_FUNCTIONS + j - solver.grid.getLocalTriangleNumber()];
		//		rhov_result[t * num + j] = rhov_buffer[t * data_size / BASIS_FUNCTIONS + j - solver.grid.getLocalTriangleNumber()];
		//		rhoE_result[t * num + j] = rhoE_buffer[t * data_size / BASIS_FUNCTIONS + j - solver.grid.getLocalTriangleNumber()];
		//	}
		//}
		for(int i = 0; i < 10; i++) {
			cout<<"rhou:"<<solver._freedom_rho[i]<<endl;
		}
		solver.outputSolution(data_size, solver._freedom_rho, solver._freedom_rhou, solver._freedom_rhov, solver._freedom_rhoE);
		cout<<"mmm"<<endl;
	}
	cout<<"complete solving flow, press ENTER to exit."<<endl;
	end = clock();
	if (solver.myid == 0) {
		cout<<"Wall time:"<<(double)(end - start) / CLOCKS_PER_SEC<<endl;
	}
	MPI_Finalize();
	return 0;
}