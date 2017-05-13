#include 
using namespace std;

class GridPartition {
public:
	vector<int> elem_index;
	vector<int> innerBoundary_index;
	vector<int> innerBoundary_edge;

	void readMeshPartitionInfo(const string& partition_file, int rank);
}