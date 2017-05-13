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

#include "../inc/cppstdheaders.h"
#include "../inc/cudarkdgsolver.h"
#include"../inc/unstructuredgrid.h"
#include<mpi.h>
using namespace std;


int main(int argc, char* argv[])
{
//	CCUDARkdgSolver solver;
//	
//
	int myid;
	int nprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//
//	CUnstructuredGrid grid;
//	grid.config_file = "input/mesh.conf";
////	grid.initializeGrid(0);
////
////	ofstream of("output/local.dat");
////	for (int i = 0; i < grid.elem_index.size() - 1; i++) {
//////		cout<<grid.vertice.at(grid.elem_index[i]).getX()<<"  "<<grid.vertice.at(grid.elem_index[i]).getX()<<endl;
//////		of<<grid.vertice.at(grid.elem_index[i]).getX()<<"  "<<grid.vertice.at(grid.elem_index[i]).getY()<<endl;
////		of<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3]).getX()<<"  "<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3]).getY()<<endl;
////		of<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3 + 1]).getX()<<"  "<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3 + 1]).getY()<<endl;
////		of<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3 + 2]).getX()<<"  "<<grid.vertice.at(grid.tri_vertice[grid.elem_index[i] * 3 + 2]).getY()<<endl;
////	}
////	of.close();
//
//	try {
//		solver.config_file = "input/main.conf";
//
//		solver.run(myid, nprocs);
//	}
//	catch ( const exception& e )
//	{
//		cout<<endl<<"Error occured, description: "<<endl;
//		cout<<e.what()<<endl;
//		cout<<endl<<"Please modify the program according to the hint."<<endl;
//		getchar();
//
//		exit(-1);
//	}
//	
//	
//	
//	cout<<"complete solving flow, press ENTER to exit."<<endl;
//	MPI_Finalize();
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

	CUnstructuredGrid grid;
	grid.config_file = "input/mesh.conf";
	grid.initializeGrid(myid,2);

	MPI_Finalize();
	return 0;
}