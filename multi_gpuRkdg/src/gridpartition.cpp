#include "../inc/gridpartition.h"

GridPartition::GridPartition() {
	
}

void GridPartition::readMeshPartitionInfo(const string& partition_file, int nprocs ) {
	area_num = nprocs;
	elem_index = new vector<int>[nprocs];
	ifstream fin(partition_file.c_str());
	if(!fin)
		throw CMyException("Failed to open partition file: " + partition_file);
	int index = 0;
	string tmp_string;
	
	while(fin.good()) {
		getline(fin, tmp_string);
		for (int i = 0; i < area_num; i++) {
			if(atoi(tmp_string.c_str()) == i) {
				elem_index[i].push_back(index);
			}
		}	
		index++;
	}
}

void GridPartition::readNeigh() {
	for (int i = 0; i < area_num; i++) {
		for (int j = 0; j < area_num; j++) {
			if (i != j) {
				if(hasCommon(i, j)) {
					
				}
			}
		}
	}
}

bool GridPartition::hasCommon(int first_elem, int second_elem) {
	int tvn = _local_vertice_num;
	for (int p = 0; p < TRIANGLE_EDGES; p++) {
		for (int q = 0; q < TRIANGLE_EDGES; q++) {
			if (tri_edge[i * TRIANGLE_EDGES + p] == tri_edge[j * TRIANGLE_EDGES + q]) {
				innerBoundary_edge.push_back(tri_edge[i  * TRIANGLE_EDGES + p]);
				CEdge boundary_edge0 = edge.at(tri_edge[i * TRIANGLE_EDGES + 1]);
				CEdge boundary_edge1 = edge.at(tri_edge[i * TRIANGLE_EDGES + 1]);
				CEdge boundary_edge2 = edge.at(tri_edge[i * TRIANGLE_EDGES + 2]);
				switch (q) {
				case 0:
					if (boundary_edge1.getStart() == boundary_edge2.getTerminal()) {
						vertice[tvn].setVertice(vertice[boundary_edge1.getStart()].getX(), 
							vertice[boundary_edge1.getStart()].getY());
					} else {
						vertice[tvn].setVertice(vertice[boundary_edge1.getTerminal()].getX(), 
							vertice[boundary_edge1.getTerminal()].getY());
					}
					++ tvn;
					break;
				case 1:
					if (boundary_edge0.getStart() == boundary_edge2.getTerminal()) {
						vertice[tvn].setVertice(vertice[boundary_edge0.getStart()].getX(), 
							vertice[boundary_edge0.getStart()].getY());
					} else {
						vertice[tvn].setVertice(vertice[boundary_edge0.getTerminal()].getX(), 
							vertice[boundary_edge0.getTerminal()].getY());
					}
					++ tvn;
					break;
				case 2:
					if (boundary_edge0.getStart() == boundary_edge1.getTerminal()) {
						vertice[tvn].setVertice(vertice[boundary_edge0.getStart()].getX(), 
							vertice[boundary_edge0.getStart()].getY());
					} else {
						vertice[tvn].setVertice(vertice[boundary_edge0.getTerminal()].getX(), 
							vertice[boundary_edge0.getTerminal()].getY());
					}
					++ tvn;
					break;
				}
				return true;
			}
		}
	}
	return false;
}
