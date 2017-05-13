#include "inc/unstructuredgrid.h"

void testLocalElemIndex() {

	CUnstructuredGrid grid0, grid1, grid2;
	grid0.config_file = "input/mesh.conf";
	grid0.initializeGrid(0, 3);
	//grid0.readMeshPartitionInfo("input/depart5.dat", 0);
	ofstream of("output/de0.dat");
	for (int i = 0; i < grid0.elem_index.size(); i++) {
		of<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3]).getX()<<"  "<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3]).getY()<<endl;
		of<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3 + 1]).getX()<<"  "<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3 + 1]).getY()<<endl;
		of<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3 + 2]).getX()<<"  "<<grid0.vertice.at(grid0.tri_vertice[grid0.elem_index[i] * 3 + 2]).getY()<<endl;
	}
	grid1.config_file = "input/mesh.conf";
	grid1.initializeGrid(1, 3);
	//grid1.readMeshPartitionInfo("input/depart5.dat", 1);
	ofstream of1("output/de1.dat");
	for (int i = 0; i < grid1.elem_index.size(); i++) {
		of1<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3]).getX()<<"  "<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3]).getY()<<endl;
		of1<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3 + 1]).getX()<<"  "<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3 + 1]).getY()<<endl;
		of1<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3 + 2]).getX()<<"  "<<grid1.vertice.at(grid1.tri_vertice[grid1.elem_index[i] * 3 + 2]).getY()<<endl;
	}


	grid2.config_file = "input/mesh.conf";
	grid2.initializeGrid(2, 3);
	//grid1.readMeshPartitionInfo("input/depart5.dat", 1);
	ofstream of2("output/de2.dat");
	for (int i = 0; i < grid2.elem_index.size(); i++) {
		of2<<grid2.vertice.at(grid2.tri_vertice[grid2.elem_index[i] * 3]).getX()<<"  "<<grid1.vertice.at(grid2.tri_vertice[grid2.elem_index[i] * 3]).getY()<<endl;
		of2<<grid2.vertice.at(grid2.tri_vertice[grid2.elem_index[i] * 3 + 1]).getX()<<"  "<<grid1.vertice.at(grid2.tri_vertice[grid2.elem_index[i] * 3 + 1]).getY()<<endl;
		of2<<grid2.vertice.at(grid2.tri_vertice[grid2.elem_index[i] * 3 + 2]).getX()<<"  "<<grid1.vertice.at(grid2.tri_vertice[grid2.elem_index[i] * 3 + 2]).getY()<<endl;
	}
}