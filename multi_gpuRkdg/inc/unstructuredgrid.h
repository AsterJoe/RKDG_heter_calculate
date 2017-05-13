/**
 * 二维非结构网格类
 * 该类为流场解法器提供网格接口
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 */
 
#pragma once

#include "cppstdheaders.h"
#include "vertice2d.h"
#include "edge.h"
#include "config.h"
#include "triangleinfo.h"
#include "defines.h"
#include<set>
#include<mpi.h>
using namespace std;

extern "C"{

	enum cl_type2{
		MESH_VERTICE_NUM  = 0,
		MESH_EDGE_NUM     = 1

	};

}

/**
 * 二维非结构网格类
 * 该类为流场解法器提供网格接口，同时提供网格单元的详细信息
 * 该类主要包含的信息： 网格顶点，网格的边，三角形组成以及详细的三角形信息（邻居，基函数等）
 */
class CUnstructuredGrid {
	// 属性
	public:
		MPI_Status status;
		MPI_Request request;

		string config_file;				/**< 网格配置文件 */
		string output_filename;			/**< 网格输出文件 */

		vector<CVertice2D> vertice;		/**< 网格顶点 */
		vector<CEdge> edge;				/**< 网格边 */
		vector<int> tri_vertice;		/**< 三角形的顶点组成 */
		vector<int> tri_edge;			/**< 三角形的边组成 */
		vector<int> elem_index;
		vector<int> elem_location;

		vector<int>* local_innerBoundary_index;
		vector<int>* neigh_innerBoundary_index;
		vector<int>* innerBoundary_edge;
		//此信息需要发送到其他节点
		//int** innerBoundary_index;

		// 此三项需要被复制到GPU上
		int *tri_neighbour;				/**< 三角形邻居信息 */
		int *tri_sharedEdge;			/**< 三角形共享边信息 */
		int *tri_flag;					/**< 三角形标记（边界条件等信息） */

		int *local_tri_neighbour;
		int *local_tri_sharedEdge;
		int *local_tri_flag;
		
		CTriangleInfo triangle_infos;	/**< 三角形额外信息 */
		CTriangleInfo local_triangle_infos;
		
	protected:
		int nprocs;
		int area_index;
		int _vertice_num;				/**< 顶点数目 */
		int _edge_num;					/**< 边数目 */
		int _triangle_num;				/**< 三角形数目 */
		int _ghost_triangle_num;		/**< 虚拟三角单元数目 */
		int _cell_num;					/**< 网格总总单元数目（包括虚拟网格） */

		int _local_vertice_num;
		int _local_edge_num;
		int _local_triangle_num;
		int _local_ghost_triangle_num;
		int _local_cell_num;
		
	private:
	
	// 方法
	public:
		/** 默认构造函数 */
		CUnstructuredGrid();
		
		/**
		 * 解析虚拟三角形数目, 如果某个单元邻居编号为-1，则需要生成虚拟网格
		 * @param[in] trineigh_filename string 三角形邻居编号文件
		 */
		void parseGhostNum(const string& trineigh_filename);

		/** 初始化网格数据 */
		void initializeGrid(int gpu_index, int nprocs);
		
		/** 初始化三角形数据 */
		void initializeTriangleInfos(void);
		
		/** 
		 * 导入三角形顶点编号
		 * @param[in] trivertice_file string 三角形顶点编号文件
		 */
		void readTriangleVertice(const string& trivertice_file);
	
		/** 
		* 导入相邻三角形编号
		* @param[in] trineigh_file string 三角形邻居编号文件
		*/
		void readTriangleNeighbour(const string& trineigh_file);
	
		/**
		* 导入三角形边编号
		* @param triedge_file string 三角形边组成文件
		*/
		void readTriangleEdge(const string& triedge_file);
	
		/** 
		* 导入顶点坐标
		* @param[in] vertice_file string 顶点坐标文件
		*/
		void readVertice(const string& vertice_file);
	
		/** 
		* 导入边的顶点编号
		* @param[in] edge_file 边的顶点编号文件
		*/
		void readEdgeVertice(const string& edge_file);

		/**
		* 导入网格剖分后此gpu计算的网格单元 
		*/
		void readMeshPartitionInfo(const string& partition_file, int gpu_index);

		void markInnerBoundary();

		int getLocalMeshInfo(int info_type);

		/**
		* 判断单元i与单元j是否有公共边
		*/
		bool hasCommonEdge(int i, int j);
	
		/** 
		* 初始化虚拟单元
		*/
		void initGhostGrid(void);

		void initLocalGhostGrid(void);
		
		/** @return integer 三角形数目 */
		int getTriangleNumber(void) const;

		int getLocalTriangleNumber(void) const;

		/** @return integer 边数目 */
		int getEdgeNumber(void) const;

		int getLocalEdgeNumber(void) const;
		
		/** @return integer 顶点数目 */
		int getVerticeNumber(void) const;

		int getLocalVerticeNumber(void) const;
		
		/** @return integer 虚拟三角形数目 */
		int getGhostTriangleNumber(void) const;

		int getLocalGhostTriangleNumber(void) const;

		/** @return integer 返回总单元数目 */
		int getCellNumber(void) const;

		int getLocalCellNumber(void) const;

		/** 初始化单元之间共享边 */
		void initSharedEdge(void);

		void initLocalSharedEdge(void);

		/** 标记网格边界条件 */
		void markBoundaryTriangles(void);

		void markLocalBoundaryTriangles(void);
		
		/** 输出网格 */
		void outputGrid(void) const;

		/** 测试网格单元是不是都是逆时针编号 */
		void testTrianglesAntiwise(void) const;
		
		/**
		 * 输出计算网格(包括虚拟网格)
		 * @param[in] filename string 虚拟网格输出文件
		 */
		void outputGridWithGhostCells(const string& filename) const;
		
		/** 析构函数 */
		~CUnstructuredGrid();
	protected:
	
	private:
};