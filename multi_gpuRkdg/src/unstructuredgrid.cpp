#include "../inc/unstructuredgrid.h"

CUnstructuredGrid::CUnstructuredGrid():
_vertice_num(0),
_edge_num(0),
_triangle_num(0),
_ghost_triangle_num(0),
_cell_num(0),
tri_neighbour(NULL),
tri_sharedEdge(NULL),
tri_flag(NULL)
{}

void CUnstructuredGrid::initializeGrid(int myid, int nprocs)
{
	this->nprocs = nprocs;
	area_index = myid;
	// 需要从配置文件里读取的配置项
	string conf_items[] = {
		"trianglenum",
		"verticenum", 
		"edgenum", 
		"trivertice_filename",
		"trineigh_filename",
		"triedge_filename", 
		"vertice_filename", 
		"edgevertice_filename",
		"output_filename"
	};
	CConfig grid_conf(config_file, conf_items, 9);
	if ( !grid_conf.parseConfigFile() )
		grid_conf.printQuitMessage();
	
	_triangle_num = atoi(grid_conf.config_items["trianglenum"].c_str());
	_vertice_num  = atoi(grid_conf.config_items["verticenum"].c_str());
	_edge_num = atoi(grid_conf.config_items["edgenum"].c_str());
	output_filename = grid_conf.config_items["output_filename"];	
	readMeshPartitionInfo("input/ElemDepart.dat", myid);
	parseGhostNum(grid_conf.config_items["trineigh_filename"]);	
	_cell_num = _triangle_num + _ghost_triangle_num;

	// 分配空间
	vertice.resize((_vertice_num+_ghost_triangle_num));
	edge.resize(_edge_num);
	tri_vertice.resize(3*_cell_num);
	tri_edge.resize(3*_cell_num);
	tri_neighbour  = new int[3*_cell_num];
	tri_sharedEdge = new int[3*_cell_num];
	tri_flag       = new int[_cell_num];

	local_innerBoundary_index = new vector<int>[nprocs - 1];
	neigh_innerBoundary_index = new vector<int>[nprocs - 1];
	readTriangleVertice(grid_conf.config_items["trivertice_filename"]);
	readTriangleNeighbour(grid_conf.config_items["trineigh_filename"]);
	readTriangleEdge(grid_conf.config_items["triedge_filename"]);
	readVertice(grid_conf.config_items["vertice_filename"]);
	readEdgeVertice(grid_conf.config_items["edgevertice_filename"]);
	cout<<"11"<<endl;
	markInnerBoundary();
	cout<<"12"<<endl;
	for (int i = 0; i < nprocs - area_index - 1; ++i) {
		int* local_exchange_elem = new int[local_innerBoundary_index[i].size()];
		int* neigh_exchange_elem = new int[local_innerBoundary_index[i].size()];
		int index = 0;
		for (int j = 0; j < local_innerBoundary_index[i].size(); j++) {
			local_exchange_elem[index] = local_innerBoundary_index[i].at(j);
			neigh_exchange_elem[index++] = neigh_innerBoundary_index[i].at(j);
		}
		MPI_Isend(local_exchange_elem, local_innerBoundary_index[i].size(), MPI_INT, i + 1, 0, MPI_COMM_WORLD, &request);
		MPI_Isend(neigh_exchange_elem, local_innerBoundary_index[i].size(), MPI_INT, i + 1, 1, MPI_COMM_WORLD, &request);
	}
	int data_size;
	for (int j = 0; j < area_index; ++j) {
		
		MPI_Probe(j,1,MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &data_size);
		int* local_exchange_elem = new int[data_size];
		int* neigh_exchange_elem = new int[data_size];
		MPI_Irecv(local_exchange_elem, data_size,MPI_INT, j, 1, MPI_COMM_WORLD, &request);
		MPI_Irecv(neigh_exchange_elem, data_size,MPI_INT, j, 0, MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
		for (int i = 0; i < data_size; i++) {
			local_innerBoundary_index[j].push_back(local_exchange_elem[i]);
			neigh_innerBoundary_index[j].push_back(neigh_exchange_elem[i]);
		}
	}
	cout<<"13"<<endl;
	_local_triangle_num = elem_index.size();
	_local_vertice_num = getLocalMeshInfo(MESH_VERTICE_NUM);
	_local_edge_num = getLocalMeshInfo(MESH_EDGE_NUM);
	//parseGhostNum;
	_local_cell_num = _local_triangle_num + _local_ghost_triangle_num;
	for (int i = 0; i < nprocs - 1; i++) {
		_local_cell_num += local_innerBoundary_index[i].size();
	}
	local_tri_neighbour = new int[3 * _local_cell_num];
	local_tri_sharedEdge = new int[3 * _local_cell_num];
	local_tri_flag = new int[_local_cell_num];
	cout<<"14"<<endl;
	initLocalTriNeigh();

//	initGhostGrid();
	cout<<"a"<<endl;
	initLocalGhostGrid();
	cout<<"b"<<endl;
//	initSharedEdge();
	initLocalSharedEdge();
	cout<<"c"<<endl;

}

void CUnstructuredGrid::initLocalTriNeigh() {
	vector<int>::iterator iter;
	for ( int i = 0; i < _local_triangle_num; ++i) {
		for ( int j = 0; j < TRIANGLE_EDGES; ++j) {
//			local_tri_neighbour[i + j *_cell_num] = tri_neighbour[elem_index[i]+j*_cell_num];
			iter = find(elem_index.begin(), elem_index.end(), tri_neighbour[elem_index[i]+j*_cell_num]);
			if (iter != elem_index.end()) {
				local_tri_neighbour[i + j * _local_cell_num] = distance(elem_index.begin(), iter);
			} else {
				local_tri_neighbour[i + j * _local_cell_num] = -2;
			}
		}
	}
}

int CUnstructuredGrid::getLocalMeshInfo(int info_type) {
		ofstream of2("output/index.dat");
	for (int i = 0; i < elem_index.size(); i++) {
		of2<<elem_index[i]<<endl;
	}
	of2.close();
	set<int> edges;
	for (int i = 0; i < _local_triangle_num - 1; i++) {
		for (int j = 0; j < TRIANGLE_EDGES; j++) {
			switch (info_type){
			case MESH_VERTICE_NUM:
				edges.insert(tri_vertice[elem_index[i] * TRIANGLE_EDGES + j]);
				break;
			case MESH_EDGE_NUM:
				edges.insert(tri_edge[elem_index[i] * TRIANGLE_EDGES + j]);
				break;
			}
		}
	}
	return edges.size();
}

void CUnstructuredGrid::initializeTriangleInfos(void)
{
	int i, j;
	double x[3], y[3];

	/*for ( i=0; i<_triangle_num; ++i )
	{
		for ( j=0; j<TRIANGLE_VERTICES; ++j )
		{
			x[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getX();
			y[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getY();
		}

		triangle_infos.initInformation(i, x, y);
	}*/

	for ( i=0; i<_local_triangle_num; ++i )
	{
		for ( j=0; j<TRIANGLE_VERTICES; ++j )
		{
			x[j] = vertice[tri_vertice[elem_index[i]*TRIANGLE_VERTICES+j]].getX();
			y[j] = vertice[tri_vertice[elem_index[i]*TRIANGLE_VERTICES+j]].getY();
		}

		triangle_infos.initInformation(i, x, y);
	}

	/*for ( ; i<_cell_num; ++i )
	{
		for ( j=0; j<TRIANGLE_VERTICES; ++j )
		{
			x[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getX();
			y[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getY();
		}

		triangle_infos.barycenter[i].setVertice((x[0]+x[1]+x[2])/3.0, (y[0]+y[1]+y[2])/3.0);
	}*/

	for ( ; i<_local_cell_num; ++i )
	{
		for ( j=0; j<TRIANGLE_VERTICES; ++j )
		{
			x[j] = vertice[tri_vertice[elem_index[i]*TRIANGLE_VERTICES+j]].getX();
			y[j] = vertice[tri_vertice[elem_index[i]*TRIANGLE_VERTICES+j]].getY();
		}

		triangle_infos.barycenter[i].setVertice((x[0]+x[1]+x[2])/3.0, (y[0]+y[1]+y[2])/3.0);
	}
}

void CUnstructuredGrid::parseGhostNum(const string& trineigh_filename)
{
	ifstream fin(trineigh_filename.c_str());
	if ( !fin )
		throw CMyException("Failed to open neighbour file: "+trineigh_filename);
	
	int s(0);
	string tmp_str;
	/*while ( fin>>tmp_str )
	{
		if ( tmp_str=="-1" )
			++ s;
	}
	_ghost_triangle_num = s;*/

	int t(0);
	vector<int>::iterator it;
	int i = 0;
	int j = 0;
	while ( fin>>tmp_str )
	{
		if ( tmp_str == "-1") 
			++ s;
		if ( tmp_str=="-1" && (find(elem_index.begin(), elem_index.end(), j) != elem_index.end()))
			++ t;
		++ i;
		if (i % 3 == 0) {
			++ j;
		}	
	}
	_ghost_triangle_num = s;
	_local_ghost_triangle_num = t;
	/*for (int i = 0; i < nprocs; i++) {
	 _local_ghost_triangle_num += innerBoundary_edge[i].size();
	}*/

	fin.close();
}

void CUnstructuredGrid::parseLocalGhostNum() {

}

void CUnstructuredGrid::readTriangleVertice(const string& trivertice_filename)
{
	ifstream fin(trivertice_filename.c_str());
	if ( !fin )
		throw CMyException("Failed to open triangle vertice file: "+trivertice_filename);
	
	string tmp_string;
	getline(fin, tmp_string);
	
	// 读取组成单元的顶点编号
	for ( int i=0; i<_triangle_num; ++i )
	{
		for ( int j=0; j<TRIANGLE_VERTICES; ++j )
			fin>>tri_vertice[TRIANGLE_VERTICES*i+j];
			
		tri_flag[i] = CELL_INTERIOR;
	}
	for ( int i = 0; i < _local_triangle_num; ++i) {
		local_tri_flag[i] = CELL_INTERIOR;
	}

	fin.close();
}

void CUnstructuredGrid::readTriangleNeighbour(const string& trineigh_filename)
{
	ifstream fin(trineigh_filename.c_str());
	if ( !fin )
		throw CMyException("Failed to open triangle neighbour file: "+trineigh_filename);

	string tmp_string;
	getline(fin, tmp_string);

	// 读取单元的相邻单元编号
	for ( int i=0; i<_triangle_num; ++i )
	{
		for ( int j=0; j<TRIANGLE_EDGES; ++j )
			fin>>tri_neighbour[i+j*_cell_num];
	}
//	vector<int>::iterator iter;
//	for ( int i = 0; i < _local_triangle_num; ++i) {
//		for ( int j = 0; j < TRIANGLE_EDGES; ++j) {
////			local_tri_neighbour[i + j *_cell_num] = tri_neighbour[elem_index[i]+j*_cell_num];
//			iter = find(elem_index.begin(), elem_index.end(), tri_neighbour[elem_index[i]+j*_cell_num]);
//			if (iter != elem_index.end()) {
//				local_tri_neighbour[i + j * _local_cell_num] = distance(elem_index.begin(), iter);
//			} else {
//				local_tri_neighbour[i + j * _local_cell_num] = -2;
//			}
//		}
//	}

	fin.close();
}

void CUnstructuredGrid::readTriangleEdge(const string& triedge_filename)
{
	ifstream fin(triedge_filename.c_str());
	
	if ( !fin )
		throw CMyException("Failed to open triangle edge file: "+triedge_filename);

	string tmp_string;
	getline(fin, tmp_string);

	// 读取网格边编号
	for ( int i=0; i<_triangle_num; ++i )
	{
		for ( int j=0; j<TRIANGLE_EDGES; ++j )
			fin>>tri_edge[i*TRIANGLE_EDGES+j];
	}

	fin.close();
}

void CUnstructuredGrid::readVertice(const string& vertice_filename)
{
	ifstream fin(vertice_filename.c_str());	
	if ( !fin )
		throw CMyException("Failed to open vertice file: "+vertice_filename);

	string tmp_string;
	getline(fin, tmp_string);
	
	double x, y;
	// 读取顶点坐标
	for ( int i=0; i<_vertice_num; ++i )
	{
		fin>>x>>y;
		vertice[i].setVertice(x,y);
	}

	fin.close();
}

void CUnstructuredGrid::readEdgeVertice(const string& edgevertice_filename)
{
	ifstream fin(edgevertice_filename.c_str());
	if ( !fin )
		throw CMyException("Failed to open edge vertice file: "+edgevertice_filename);

	string tmp_string;
	getline(fin, tmp_string);

	int start, terminal;
	// 读取顶点坐标
	for ( int i=0; i<_edge_num; ++i )
	{
		fin>>start>>terminal;
		edge.at(i).setEdge(start,terminal);
	}

	fin.close();
}

void CUnstructuredGrid::readMeshPartitionInfo(const string& partition_file, int rank) {
	ifstream fin(partition_file.c_str());
	if(!fin)
		throw CMyException("Failed to open partition file: " + partition_file);
	int index = 0;
	string tmp_string;
	
	while(fin.good()) {
		getline(fin, tmp_string);
		elem_location.push_back(atoi(tmp_string.c_str()));
		if(atoi(tmp_string.c_str()) == rank) {
			elem_index.push_back(index);
		}
		index++;
	}
}

void CUnstructuredGrid::markInnerBoundary() {
	ofstream ofs("output/boundary.dat");
	//ofstream ofs2("output/boundary2.dat");
	//for (int i = 0; i < elem_index.size(); i++) {
	//	for (int j = 0; j < TRIANGLE_EDGES; j++) {
	//		if (tri_neighbour[elem_index[i] + j * _cell_num] != -1) {
	//			int loc = elem_location[tri_neighbour[elem_index[i] + j * _cell_num]];
	//			
	//			if (area_index < loc) {
	//				local_innerBoundary_index[loc - 1].push_back(i);
	//				neigh_innerBoundary_index[loc - 1].push_back(j);
	//			}
	//		}
	//	}
	//}

	for (int i = 0; i < elem_index.size(); i++) {
		for (int j = 0; j < _triangle_num; j++) {
			if (area_index < elem_location[j]) {
				if (hasCommonEdge(elem_index[i], j)) {
					int loc = elem_location[j];
					int t = loc < area_index ? loc : loc - 1;
					local_innerBoundary_index[t].push_back(i);
					neigh_innerBoundary_index[t].push_back(j);
					ofs<<i<<" "<<j<<endl;
				}
				/*	vector<int>::iterator it;

				it = find(innerBoundary_index.begin(), innerBoundary_index.end(), i);
				if (it != innerBoundary_index.end()) {
					innerBoundary_index.push_back(i);
				}*/
			}
		}
	}
}

bool CUnstructuredGrid::hasCommonEdge(int i, int j) {
	int tvn = _local_vertice_num;
	int ten = _local_triangle_num;
	bool result = false;
	for (int p = 0; p < TRIANGLE_EDGES; p++) {
		for (int q = 0; q < TRIANGLE_EDGES; q++) {
			if (tri_edge[i * TRIANGLE_EDGES + p] == tri_edge[j * TRIANGLE_EDGES + q]) {
				//innerBoundary_edge[elem_location[j] - 1].push_back(tri_edge[i  * TRIANGLE_EDGES + p]);
				
				//CEdge boundary_edge0 = edge.at(tri_edge[first * TRIANGLE_EDGES]);
				//CEdge boundary_edge1 = edge.at(tri_edge[first * TRIANGLE_EDGES + 1]);
				//CEdge boundary_edge2 = edge.at(tri_edge[first * TRIANGLE_EDGES + 2]);
				//
				//switch (q) {
				//case 0:
				//	if (boundary_edge1.getStart() == boundary_edge2.getTerminal()) {
				//		//vertice[tvn].setVertice(vertice[boundary_edge1.getStart()].getX(), 
				//		//	vertice[boundary_edge1.getStart()].getY());

				//		// 添加虚拟单元的邻居关系
				//		local_tri_neighbour[ten] = i;
				//		local_tri_neighbour[ten+_cell_num] = -1;
				//		local_tri_neighbour[ten+2*_cell_num] = -1;
				//		local_tri_sharedEdge[ten] = 0;

				//		// 更新所在单元邻居列表
				//		local_tri_neighbour[i] = ten;
				//	} else {
				//		/*vertice[tvn].setVertice(vertice[boundary_edge1.getTerminal()].getX(), 
				//			vertice[boundary_edge1.getTerminal()].getY());*/
				//		// 添加虚拟单元的邻居关系
				//		local_tri_neighbour[ten] = i;
				//		local_tri_neighbour[ten+_cell_num] = -1;
				//		local_tri_neighbour[ten+2*_cell_num] = -1;
				//		local_tri_sharedEdge[ten] = 0;

				//		// 更新所在单元邻居列表
				//		local_tri_neighbour[i] = ten;
				//	}
				//	++ tvn;
				//	++ ten;
				//	break;
				//case 1:
				//	if (boundary_edge0.getStart() == boundary_edge2.getTerminal()) {
				//		/*vertice[tvn].setVertice(vertice[boundary_edge0.getStart()].getX(), 
				//			vertice[boundary_edge0.getStart()].getY());*/
				//	} else {
				//		/*vertice[tvn].setVertice(vertice[boundary_edge0.getTerminal()].getX(), 
				//			vertice[boundary_edge0.getTerminal()].getY());*/
				//	}
				//	++ tvn;
				//	break;
				//case 2:
				//	if (boundary_edge0.getStart() == boundary_edge1.getTerminal()) {
				//		/*vertice[tvn].setVertice(vertice[boundary_edge0.getStart()].getX(), 
				//			vertice[boundary_edge0.getStart()].getY());*/
				//	} else {
				//		/*vertice[tvn].setVertice(vertice[boundary_edge0.getTerminal()].getX(), 
				//			vertice[boundary_edge0.getTerminal()].getY());*/
				//	}
				//	++ tvn;
				//	break;
				//}
				result = true;
			}
		}
	}
	return result;
}

void CUnstructuredGrid::initGhostGrid(void)
{
	int i, j, k, e;
	double x1, x2, x3, y1, y2, y3, xt, yt;
	
	int ten = _triangle_num;
	int tvn = _vertice_num;
	int local_ten = _local_triangle_num;
	int local_tvn = _local_vertice_num;
	bool isLocal;

	for ( e=0; e<_triangle_num; ++e )
	{
		if (find(elem_index.begin(), elem_index.end(), e) == elem_index.end()) {
			isLocal = false;
		} else {
			isLocal = true;
		}
		if ( tri_neighbour[e]==-1 )
		{
			// 单元编号和点坐标
			i = tri_vertice[e*TRIANGLE_VERTICES];
			j = tri_vertice[e*TRIANGLE_VERTICES+1];
			k = tri_vertice[e*TRIANGLE_VERTICES+2];

			x1 = vertice[i].getX();
			x2 = vertice[j].getX();
			x3 = vertice[k].getX();
			y1 = vertice[i].getY();
			y2 = vertice[j].getY();
			y3 = vertice[k].getY();

			// 新增顶点的坐标
			xt = x2 + x3 - x1;
			yt = y2 + y3 - y1;
			vertice[tvn].setVertice(xt,yt);

			// 更新单元列表
			tri_vertice[ten*TRIANGLE_VERTICES] = tvn;
			tri_vertice[ten*TRIANGLE_VERTICES+1] = k;
			tri_vertice[ten*TRIANGLE_VERTICES+2] = j;
			
			// 添加虚拟单元的邻居关系
			tri_neighbour[ten] = e;
			tri_neighbour[ten+_cell_num] = -1;
			tri_neighbour[ten+2*_cell_num] = -1;
			tri_sharedEdge[ten] = 0;

			// 更新所在单元邻居列表
			tri_neighbour[e] = ten;


			if (isLocal) {
				local_tri_neighbour[local_ten] = e;
				local_tri_neighbour[local_ten+_local_cell_num] = -1;
				local_tri_neighbour[local_ten+2*_local_cell_num] = -1;
				local_tri_sharedEdge[local_ten] = 0;

				local_tri_neighbour[e] = ten;
			}

			++ tvn;
			++ ten;

			++local_tvn;
			++local_ten;
		}
		else if ( tri_neighbour[e+_cell_num]==-1 )
		{
			// 单元编号和点坐标
			i = tri_vertice[e*TRIANGLE_VERTICES];
			j = tri_vertice[e*TRIANGLE_VERTICES+1];
			k = tri_vertice[e*TRIANGLE_VERTICES+2];

			x1 = vertice[i].getX();
			x2 = vertice[j].getX();
			x3 = vertice[k].getX();
			y1 = vertice[i].getY();
			y2 = vertice[j].getY();
			y3 = vertice[k].getY();

			// 新增顶点坐标
			xt = x1 + x3 - x2;
			yt = y1 + y3 - y2;
			vertice[tvn].setVertice(xt,yt);

			// 更新单元列表
			tri_vertice[ten*TRIANGLE_VERTICES] = tvn;
			tri_vertice[ten*TRIANGLE_VERTICES+1] = i;
			tri_vertice[ten*TRIANGLE_VERTICES+2] = k;
			
			// 添加虚拟单元的邻居关系
			tri_neighbour[ten] = e;
			tri_neighbour[ten+_cell_num] = -1;
			tri_neighbour[ten+2*_cell_num] = -1;
			tri_sharedEdge[ten] = 1;

			// 更新所在单元邻居列表
			tri_neighbour[e+_cell_num] = ten;

			if (isLocal) {
				local_tri_neighbour[local_ten] = e;
				local_tri_neighbour[local_ten+_local_cell_num] = -1;
				local_tri_neighbour[local_ten+2*_local_cell_num] = -1;
				local_tri_sharedEdge[local_ten] = 1;

				local_tri_neighbour[e+_local_cell_num] = ten;
			}

			++ tvn;
			++ ten;

			++local_tvn;
			++local_ten;
		}
		else if ( tri_neighbour[e+2*_cell_num]==-1 )
		{
			// 单元编号和点坐标
			i = tri_vertice[e*TRIANGLE_VERTICES];
			j = tri_vertice[e*TRIANGLE_VERTICES+1];
			k = tri_vertice[e*TRIANGLE_VERTICES+2];

			x1 = vertice[i].getX();
			x2 = vertice[j].getX();
			x3 = vertice[k].getX();
			y1 = vertice[i].getY();
			y2 = vertice[j].getY();
			y3 = vertice[k].getY();

			// 新增顶点坐标
			xt = x1 + x2 - x3;
			yt = y1 + y2 - y3;
			vertice[tvn].setVertice(xt,yt);

			// 更新单元列表
			tri_vertice[ten*TRIANGLE_VERTICES] = tvn;
			tri_vertice[ten*TRIANGLE_VERTICES+1] = i;
			tri_vertice[ten*TRIANGLE_VERTICES+2] = j;

			// 添加虚拟单元的邻居关系
			tri_neighbour[ten] = e;
			tri_neighbour[ten+_cell_num] = -1;
			tri_neighbour[ten+2*_cell_num] = -1;
			tri_sharedEdge[ten] = 2;
			
			// 更新所在单元邻居列表
			tri_neighbour[e+2*_cell_num] = ten;

			if (isLocal) {
				local_tri_neighbour[local_ten] = e;
				local_tri_neighbour[local_ten+_local_cell_num] = -1;
				local_tri_neighbour[local_ten+2*_local_cell_num] = -1;
				local_tri_sharedEdge[local_ten+2*_local_cell_num] = 2;

				local_tri_neighbour[e] = ten;
			}

			++ tvn;
			++ ten;

			++local_tvn;
			++local_ten;
		}
		
	}

	for (e=0; e<_local_triangle_num; ++e){

	}
}

void CUnstructuredGrid::initLocalGhostGrid(void) {
	int i, j, k, e;
	double x1, x2, x3, y1, y2, y3, xt, yt;
	
	int ten = _local_triangle_num;
	int ten_global = _triangle_num;
	int tvn = _local_vertice_num;

	for (int i = 0; i < nprocs - 1; i++) {
		for ( int j = 0; j < local_innerBoundary_index[i].size(); j++) {
			local_tri_neighbour[ten] = local_innerBoundary_index[i].at(j);
			local_tri_neighbour[ten + _local_cell_num] = -1;
			local_tri_neighbour[ten + 2 * _local_cell_num] = -1;
			for (int t = 0; t < TRIANGLE_EDGES; t++) {
				if (local_tri_neighbour[local_innerBoundary_index[i].at(j) + t * _local_cell_num] == -2) {
					local_tri_sharedEdge[ten] = t;
				}
			}	
			local_tri_neighbour[local_innerBoundary_index[i].at(j)] = ten;

			++ ten;
		}
	}
	for ( e=0; e<_local_triangle_num; ++e )
	{
		if (tri_neighbour[elem_index[e]]==-1 )
		{
			// 单元编号和点坐标
			i = tri_vertice[elem_index[e]*TRIANGLE_VERTICES];
			j = tri_vertice[elem_index[e]*TRIANGLE_VERTICES+1];
			k = tri_vertice[elem_index[e]*TRIANGLE_VERTICES+2];

			x1 = vertice[i].getX();
			x2 = vertice[j].getX();
			x3 = vertice[k].getX();
			y1 = vertice[i].getY();
			y2 = vertice[j].getY();
			y3 = vertice[k].getY();

			// 新增顶点的坐标
			xt = x2 + x3 - x1;
			yt = y2 + y3 - y1;
			vertice[tvn].setVertice(xt,yt);

			// 更新单元列表
			tri_vertice[ten_global*TRIANGLE_VERTICES] = tvn;
			tri_vertice[ten_global*TRIANGLE_VERTICES+1] = k;
			tri_vertice[ten_global*TRIANGLE_VERTICES+2] = j;
			
			// 添加虚拟单元的邻居关系
			/*tri_neighbour[ten] = elem_index[e];
			tri_neighbour[ten+_cell_num] = -1;
			tri_neighbour[ten+2*_cell_num] = -1;
			tri_sharedEdge[ten] = 0;*/
			local_tri_neighbour[ten] = elem_index[e];
			local_tri_neighbour[ten + _local_cell_num] = -1;
			local_tri_neighbour[ten + 2 * _local_cell_num] = -1;
			local_tri_sharedEdge[ten] = 0;

			// 更新所在单元邻居列表
			//tri_neighbour[e] = ten;
			local_tri_neighbour[e] = ten;

			++ tvn;
			++ ten;
			++ ten_global;
		}
		else if ( tri_neighbour[elem_index[e] + _cell_num]==-1 )
		{
			// 单元编号和点坐标
			i = tri_vertice[e*TRIANGLE_VERTICES];
			j = tri_vertice[e*TRIANGLE_VERTICES+1];
			k = tri_vertice[e*TRIANGLE_VERTICES+2];

			x1 = vertice[i].getX();
			x2 = vertice[j].getX();
			x3 = vertice[k].getX();
			y1 = vertice[i].getY();
			y2 = vertice[j].getY();
			y3 = vertice[k].getY();

			// 新增顶点坐标
			xt = x1 + x3 - x2;
			yt = y1 + y3 - y2;
			vertice[tvn].setVertice(xt,yt);

			// 更新单元列表
			tri_vertice[ten_global*TRIANGLE_VERTICES] = tvn;
			tri_vertice[ten_global*TRIANGLE_VERTICES+1] = i;
			tri_vertice[ten_global*TRIANGLE_VERTICES+2] = k;
			
			// 添加虚拟单元的邻居关系
			/*tri_neighbour[ten] = elem_index[e];
			tri_neighbour[ten+_cell_num] = -1;
			tri_neighbour[ten+2*_cell_num] = -1;
			tri_sharedEdge[ten] = 1;*/

			local_tri_neighbour[ten] = elem_index[e];
			local_tri_neighbour[ten + _local_cell_num] = -1;
			local_tri_neighbour[ten + 2 * _local_cell_num] = -1;
			local_tri_sharedEdge[ten] = 1;

			// 更新所在单元邻居列表
			tri_neighbour[e+_local_cell_num] = ten;

			++ tvn;
			++ ten;
			++ ten_global;
		}
		else if ( tri_neighbour[elem_index[e] + 2 * _cell_num]==-1 )
		{
			// 单元编号和点坐标
			i = tri_vertice[e*TRIANGLE_VERTICES];
			j = tri_vertice[e*TRIANGLE_VERTICES+1];
			k = tri_vertice[e*TRIANGLE_VERTICES+2];

			x1 = vertice[i].getX();
			x2 = vertice[j].getX();
			x3 = vertice[k].getX();
			y1 = vertice[i].getY();
			y2 = vertice[j].getY();
			y3 = vertice[k].getY();

			// 新增顶点坐标
			xt = x1 + x2 - x3;
			yt = y1 + y2 - y3;
			vertice[tvn].setVertice(xt,yt);

			// 更新单元列表
			tri_vertice[ten_global*TRIANGLE_VERTICES] = tvn;
			tri_vertice[ten_global*TRIANGLE_VERTICES+1] = i;
			tri_vertice[ten_global*TRIANGLE_VERTICES+2] = j;

			// 添加虚拟单元的邻居关系
			/*tri_neighbour[ten] = elem_index[e];
			tri_neighbour[ten+_cell_num] = -1;
			tri_neighbour[ten+2*_cell_num] = -1;
			tri_sharedEdge[ten] = 2;*/

			local_tri_neighbour[ten] = elem_index[e];
			local_tri_neighbour[ten + _local_cell_num] = -1;
			local_tri_neighbour[ten + 2 * _local_cell_num] = -1;
			local_tri_sharedEdge[ten] = 2;
			
			// 更新所在单元邻居列表
			tri_neighbour[e+2*_local_cell_num] = ten;

			++ tvn;
			++ ten;
			++ ten_global;
		} 
		
	}
}

int CUnstructuredGrid::getGhostTriangleNumber(void) const
{
	return _ghost_triangle_num;
}

int CUnstructuredGrid::getLocalGhostTriangleNumber(void) const
{
	return _local_ghost_triangle_num;
}

int CUnstructuredGrid::getEdgeNumber() const
{
	return _edge_num;
}

int CUnstructuredGrid::getLocalEdgeNumber() const
{
	return _local_edge_num;
}

int CUnstructuredGrid::getTriangleNumber() const
{
	return _triangle_num;
}

int CUnstructuredGrid::getLocalTriangleNumber() const
{
	return _local_triangle_num;
}

int CUnstructuredGrid::getVerticeNumber() const
{
	return _vertice_num;
}

int CUnstructuredGrid::getLocalVerticeNumber() const
{
	return _local_vertice_num;
}

int CUnstructuredGrid::getCellNumber() const
{
	return _cell_num;
}

int CUnstructuredGrid::getLocalCellNumber() const
{
	return _local_cell_num;
}

void CUnstructuredGrid::markBoundaryTriangles( void )
{
	double bcx, bcy;

	for ( int i=_triangle_num; i<_cell_num; ++i )
	{
		// 重心坐标
		bcx = triangle_infos.getBarycenter(i).getX();
		bcy = triangle_infos.getBarycenter(i).getY();

		if ( bcx<-15.0 ) // 左方来流
		{
			tri_flag[i] = CELL_FARFIELD;
		} 
		else if ( bcx>15.0 ) // 右方去流
		{
		//	_triangles[i].setFlag(bnd_outflow);
			tri_flag[i] = CELL_FARFIELD;
		}
		else if ( bcy<-15.0 ) // 下方来流
		{
			tri_flag[i] = CELL_FARFIELD;
		}
		else if ( bcy>15.0 ) // 上方来流
		{
			tri_flag[i] = CELL_FARFIELD;
		}
		else
		{
			tri_flag[i] = CELL_REFLECTION; // 机翼边界
		}
	}
}

void CUnstructuredGrid::markLocalBoundaryTriangles() {
	double bcx, bcy;

	for ( int i=_local_triangle_num; i<_local_cell_num; ++i )
	{
		// 重心坐标
		bcx = triangle_infos.getBarycenter(elem_index[i]).getX();
		bcy = triangle_infos.getBarycenter(elem_index[i]).getY();

		if ( bcx<-15.0 ) // 左方来流
		{
			local_tri_flag[i] = CELL_FARFIELD;
		} 
		else if ( bcx>15.0 ) // 右方去流
		{
		//	_triangles[i].setFlag(bnd_outflow);
			local_tri_flag[i] = CELL_FARFIELD;
		}
		else if ( bcy<-15.0 ) // 下方来流
		{
			local_tri_flag[i] = CELL_FARFIELD;
		}
		else if ( bcy>15.0 ) // 上方来流
		{
			local_tri_flag[i] = CELL_FARFIELD;
		}
		else
		{
			local_tri_flag[i] = CELL_REFLECTION; // 机翼边界
		}
	}

	/*for (int i = 0; i < innerBoundary_index.size(); i++) {
		local_tri_flag[i] = CELL_INNER_BOUNDARY;
	}*/
}

void CUnstructuredGrid::testTrianglesAntiwise() const
{
	int i, j;
	double x[3], y[3];
	double delt;

	for ( i=0; i<_triangle_num; ++i )
	{
		// 坐标
		for ( j=0; j<TRIANGLE_VERTICES; ++j )
		{
			x[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getX();
			y[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getY();
		}

		// 根据行列式判断是否逆序
		delt = (x[1]-x[0])*(y[2]-y[1]) - (y[1]-y[0])*(x[2]-x[1]);

		if ( delt<0 )
		{
			cout<<"The order of vertices in  "<<i+1<<"th triangle is not antiwise."<<endl;
		}
	}
}

void CUnstructuredGrid::outputGrid(void) const
{
	// 先输出网格单元信息
	ofstream fout(output_filename.c_str());
	if ( !fout )
	{
		cout<<"\nFailed to open grid output file: "<<output_filename<<endl;
		cout<<"Grid output will be omitted."<<endl;
		return;
	}
	
	cout<<"\n========================================\nGrid information:"<<endl;
	cout<<"triangle number: "<<_triangle_num<<endl;
	cout<<"vertice number:  "<<_vertice_num<<endl;
	cout<<"edge number:     "<<_edge_num<<endl;
	cout<<"========================================"<<endl;
		
	fout<<"TITLE=RKDG"<<endl;
	fout<<"VARIABLES="<<"X"<<" , "<<"Y"<< endl;
	fout<<"ZONE N= "<<_vertice_num<<" , "<<"E= "<<_triangle_num<<" , "<<"F=FEPOINT"<<" , "<<"ET=TRIANGLE"<<endl;

	fout.precision(10);
	int i;

	// 输出点坐标
	for ( i=0; i<_vertice_num; ++i )
	{
		fout<<vertice.at(i).getX()<<"    "<<vertice.at(i).getY()<<endl;
	}

	// 输出单元顶点
	for ( i=0; i<_triangle_num; ++i )
	{
		// 由于行编号是从1开始，故而需要+1
		for ( int j=0; j<TRIANGLE_VERTICES; ++j )
			fout<<tri_vertice[TRIANGLE_VERTICES*i+j]+1<<"    ";
			
		fout<<endl;
	}

	fout.close();
}

void CUnstructuredGrid::outputGridWithGhostCells(const string& filename) const
{
	ofstream fout(filename.c_str(), ios::out);
	if ( !fout )
	{
		cout<<"Failed to open ghost mesh file: "<<filename<<" and output will be omitted."<<endl<<endl;
		return;
	}
	int vnum = _vertice_num + _ghost_triangle_num;

	fout<<"TITLE=mesh-with-ghost-cells"<<endl;
	fout<<"VARIABLES="<<"x"<<" , Y"<<endl;
	fout<<"ZONE N="<<vnum<<" , E="<<_cell_num<<" , F=FEPOINT , ET=TRIANGLE"<<endl;

	int i;
	for ( i=0; i<vnum; ++i )
	{
		fout<<vertice[i].getX()<<"    "<<vertice[i].getY()<<endl;
	}

	for ( i=0; i<_cell_num; ++i )
	{
		for ( int j=0; j<TRIANGLE_VERTICES; ++j )
			fout<<tri_vertice[i*TRIANGLE_VERTICES+j]+1<<"    ";
		
		fout<<endl;
	}

	fout.close();
}

void CUnstructuredGrid::initSharedEdge()
{
	int neigh_index, edge_index;

	for ( int e=0; e<_triangle_num; ++e )
	{
		for ( int k=0; k<TRIANGLE_EDGES; ++k )
		{
			neigh_index = tri_neighbour[e+k*_cell_num];
			// 如果该单元为虚拟单元，根据虚拟单元生成规则，必定是0号边
			if ( neigh_index>=_triangle_num )
			{
				tri_sharedEdge[e+k*_cell_num] = 0;
			} 
			else
			{
				edge_index = tri_edge[e*TRIANGLE_EDGES+k];
				if ( edge_index==tri_edge[neigh_index*TRIANGLE_EDGES] )
				{
					tri_sharedEdge[e+k*_cell_num] = 0;
				} 
				else if ( edge_index==tri_edge[neigh_index*TRIANGLE_EDGES+1] )
				{
					tri_sharedEdge[e+k*_cell_num] = 1;
				}
				else if ( edge_index==tri_edge[neigh_index*TRIANGLE_EDGES+2] )
				{
					tri_sharedEdge[e+k*_cell_num] = 2;
				}
				else
				{
					cout<<"\nedge relation error:"<<endl;
					cout<<"detail: triangle "<<e+1<<"th's  "<<k+1<<"th edge is the share edge of "<<neigh_index+1<<"th triangle, but indice of "<<neigh_index+1<<"th edge is not correct."<<endl;
					cout<<"press ENTER to exit."<<endl;
					getchar();
					exit(-1);
				}
			}
		}
	}
}

void CUnstructuredGrid::initLocalSharedEdge() {
	int neigh_index, edge_index;
	for ( int e=0; e<_local_triangle_num; ++e )
	{
		for ( int k=0; k<TRIANGLE_EDGES; ++k )
		{
			//neigh_index = tri_neighbour[elem_index[e]+k*_cell_num];
			neigh_index = local_tri_neighbour[e + k * _local_cell_num];
			// 如果该单元为虚拟单元，根据虚拟单元生成规则，必定是0号边
			if ( neigh_index>=_local_triangle_num )
			{
				local_tri_sharedEdge[e+k*_local_cell_num] = 0;
			}
			else
			{
				edge_index = tri_edge[elem_index[e]*TRIANGLE_EDGES+k];
				if ( edge_index==tri_edge[elem_index[neigh_index]*TRIANGLE_EDGES] )
				{
					local_tri_sharedEdge[e+k*_cell_num] = 0;
				} 
				else if ( edge_index==tri_edge[elem_index[neigh_index]*TRIANGLE_EDGES+1] )
				{
					local_tri_sharedEdge[e+k*_cell_num] = 1;
				}
				else if ( edge_index==tri_edge[elem_index[neigh_index]*TRIANGLE_EDGES+2] )
				{
					local_tri_sharedEdge[e+k*_cell_num] = 2;
				}
				else
				{
					cout<<"\nedge relation error:"<<endl;
					cout<<"detail: triangle "<<e+1<<"th's  "<<k+1<<"th edge is the share edge of "<<neigh_index+1<<"th triangle, but indice of "<<neigh_index+1<<"th edge is not correct."<<endl;
					cout<<"press ENTER to exit."<<endl;
					getchar();
					exit(-1);
				}

			}
		}
	}
}

CUnstructuredGrid::~CUnstructuredGrid()
{
	// 释放资源
	delete []tri_neighbour;
	delete []tri_sharedEdge;
	delete []tri_flag;
}