#include "cppstdheaders.h"
#include "myexception.h"
#include "defines.h"
using namespace std;

struct GridInfo {
	int triangle_num;
	int vertice_num;

}

class GridPartition {
private:
	int area_num;

public:
	vector<int>* elem_index;
	vector<int> innerBoundary_index;
	vector<int> innerBoundary_edge;
	

	GridPartition();

	void readMeshPartitionInfo(const string& partition_file, int nprocs);

	void readNeigh();

	bool hasCommon(int first_elem, int second_elem);
}