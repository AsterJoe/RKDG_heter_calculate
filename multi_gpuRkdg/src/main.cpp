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
		solver.config_file = "input/main.conf";

		solver.run(myid, nprocs);
		SendBoundary sendBoundary(&solver);
		sendBoundary.send();
		solver.runNext();
		RkdgAdvance advance(&solver);
		advance.advance();
		cout<<"end advance"<<endl;
		solver.runAfter();
	/*}
	catch ( const exception& e )
	{
		cout<<endl<<"Error occured, description: "<<endl;
		cout<<e.what()<<endl;
		cout<<endl<<"Please modify the program according to the hint."<<endl;
		getchar();

		exit(-1);
	}*/
	cout<<"end run after"<<endl;
	cout<<"987"<<endl;
	MPI_Barrier(MPI_COMM_WORLD);
	cout<<"654"<<endl;
	double *rho_result, *rhou_result, *rhov_result, *rhoE_result;
	double *rho_buffer, *rhou_buffer, *rhov_buffer, *rhoE_buffer;
	int num = solver.grid.getTriangleNumber();
	rho_result = (double*)malloc(BASIS_FUNCTIONS * solver.grid.getTriangleNumber() * sizeof(double));
	rhou_result = (double*)malloc(BASIS_FUNCTIONS * solver.grid.getTriangleNumber() * sizeof(double));
	rhov_result = (double*)malloc(BASIS_FUNCTIONS * solver.grid.getTriangleNumber() * sizeof(double));
	rhoE_result = (double*)malloc(BASIS_FUNCTIONS * solver.grid.getTriangleNumber() * sizeof(double));
	/*MPI_Gather(solver._freedom_rho, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
		result_rho, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(solver._freedom_rhou, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
		result_rhou, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(solver._freedom_rhov, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
		result_rhov, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(solver._freedom_rhoE, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE,
		result_rhoE, solver.grid.getLocalTriangleNumber(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
*/
	if (myid != 0) {
		MPI_Isend(solver._freedom_rho, solver.grid.getLocalCellNumber() * BASIS_FUNCTIONS, MPI_DOUBLE,
			0, 10000,MPI_COMM_WORLD, &request);
		MPI_Isend(solver._freedom_rhou, solver.grid.getLocalCellNumber() * BASIS_FUNCTIONS, MPI_DOUBLE,
			0, 10001,MPI_COMM_WORLD, &request);
		MPI_Isend(solver._freedom_rhov, solver.grid.getLocalCellNumber() * BASIS_FUNCTIONS, MPI_DOUBLE,
			0, 10002,MPI_COMM_WORLD, &request);
		MPI_Isend(solver._freedom_rhoE, solver.grid.getLocalCellNumber() * BASIS_FUNCTIONS, MPI_DOUBLE,
			0, 10003,MPI_COMM_WORLD, &request);
		cout<<"send num:"<<solver.grid.getLocalCellNumber() * BASIS_FUNCTIONS<<endl;
	} 
	cout<<"mmm"<<endl;
	int data_size;
	if (myid == 0) {
		MPI_Probe(1,10000,MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &data_size);
		rho_buffer = (double*)malloc(data_size*sizeof(double));
		rhou_buffer = (double*)malloc(data_size*sizeof(double));
		rhov_buffer = (double*)malloc(data_size*sizeof(double));
		rhoE_buffer = (double*)malloc(data_size*sizeof(double));
		cout<<"datasize:"<<data_size<<endl;
		MPI_Irecv(rho_buffer, data_size, MPI_DOUBLE, 1, 10000, MPI_COMM_WORLD, &req[0]);
		MPI_Irecv(rhou_buffer, data_size, MPI_DOUBLE, 1, 10001, MPI_COMM_WORLD, &req[1]);
		MPI_Irecv(rhov_buffer, data_size, MPI_DOUBLE, 1, 10002, MPI_COMM_WORLD, &req[2]);
		MPI_Irecv(rhoE_buffer, data_size, MPI_DOUBLE, 1, 10003, MPI_COMM_WORLD, &req[3]);
		cout<<"nnn"<<endl;
		/*for(int i = 0; i < 10; i++) {
			cout<<"freedom rhou:"<<solver._freedom_rhou[i]<<endl;
		}*/
		int i =0;
		for (int t = 0; t < BASIS_FUNCTIONS; t++) {
			for(i = 0; i < solver.grid.getLocalTriangleNumber(); i++) {
				rho_result[t * num + i] = solver._freedom_rho[t * solver.grid.getLocalCellNumber() + i];
				rhou_result[t * num + i] = solver._freedom_rhou[t * solver.grid.getLocalCellNumber() + i];
				rhov_result[t * num + i] = solver._freedom_rhov[t * solver.grid.getLocalCellNumber() + i];
				rhoE_result[t * num + i] = solver._freedom_rhoE[t * solver.grid.getLocalCellNumber() + i];
			}
		}
		/*for(int i = 0; i < 10; i++) {
		cout<<"rhou:"<<rhou_result[i]<<endl;
		}*/
		MPI_Waitall(4, req, &status);
		cout<<"sss"<<endl;
		for (int t = 0; t < BASIS_FUNCTIONS; t++) {
			for (int j = solver.grid.getLocalTriangleNumber(); j < solver.grid.getTriangleNumber(); j++) {
				rho_result[t * num + j] = rho_buffer[t * data_size / BASIS_FUNCTIONS + j - solver.grid.getLocalTriangleNumber()];
				rhou_result[t * num + j] = rhou_buffer[t * data_size / BASIS_FUNCTIONS + j - solver.grid.getLocalTriangleNumber()];
				rhov_result[t * num + j] = rhov_buffer[t * data_size / BASIS_FUNCTIONS + j - solver.grid.getLocalTriangleNumber()];
				rhoE_result[t * num + j] = rhoE_buffer[t * data_size / BASIS_FUNCTIONS + j - solver.grid.getLocalTriangleNumber()];
			}
		}
		/*for(int i = 0; i < 10; i++) {
		cout<<"rhou:"<<rhou_result[i]<<endl;
		}*/
		cout<<"ttt"<<endl;
		solver.outputSolution(rho_result, rhou_result, rhov_result, rhoE_result);
		cout<<"output end"<<endl;
	}
	cout<<"complete solving flow, press ENTER to exit."<<endl;

	
	MPI_Finalize();
	
	//CUnstructuredGrid grid0, grid1, grid2;
	//grid0.config_file = "input/mesh.conf";
	//grid0.initializeGrid(0, 2);
	////grid0.readMeshPartitionInfo("input/depart5.dat", 0);
	//ofstream of("output/de0.dat");
	//for (int i = 0; i < grid0.elem_index.size(); i++) {
	//	of<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3]).getX()<<"  "<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3]).getY()<<endl;
	//	of<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3 + 1]).getX()<<"  "<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3 + 1]).getY()<<endl;
	//	of<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3 + 2]).getX()<<"  "<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3 + 2]).getY()<<endl;
	//}
	//ofstream of2("output/innerBoundary.dat");
	//for (int i = 0; i < grid0.innerBoundary_edge[0].size(); i++) {
	//	of2<<grid0.vertice.at(grid0.edge.at(grid0.innerBoundary_edge[0].at(i)).getStart()).getX()<<"  "<<grid0.vertice.at(grid0.edge.at(grid0.innerBoundary_edge[0].at(i)).getStart()).getY()<<endl;
	//	of2<<grid0.vertice.at(grid0.edge.at(grid0.innerBoundary_edge[0].at(i)).getTerminal()).getX()<<"  "<<grid0.vertice.at(grid0.edge.at(grid0.innerBoundary_edge[0].at(i)).getTerminal()).getY()<<endl;
	//}
	//grid1.config_file = "input/mesh.conf";
	//grid1.initializeGrid(1, 2);
	////grid1.readMeshPartitionInfo("input/depart5.dat", 1);
	//ofstream of1("output/de1.dat");
	//for (int i = 0; i < grid1.elem_index.size(); i++) {
	//	of1<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3]).getX()<<"  "<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3]).getY()<<endl;
	//	of1<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3 + 1]).getX()<<"  "<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3 + 1]).getY()<<endl;
	//	of1<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3 + 2]).getX()<<"  "<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3 + 2]).getY()<<endl;
	//}



	//CUnstructuredGrid grid;
	//grid.config_file = "input/mesh.conf";
	//grid.initializeGrid(myid,2);

	//MPI_Finalize();
	return 0;
}