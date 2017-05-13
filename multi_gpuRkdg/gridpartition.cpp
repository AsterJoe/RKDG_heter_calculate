#include "../inc

void GridPartition::readMeshPartitionInfo(const string& partition_file, int rank ) {
	ifstream fin(partition_file.c_str());
	if(!fin)
		throw CMyException("Failed to open partition file: " + partition_file);
	int index = 0;
	string tmp_string;
	
	while(fin.good()) {
		getline(fin, tmp_string);
		if(atoi(tmp_string.c_str()) == rank) {
			elem_index.push_back(index);
		}
		index++;
	}
}
