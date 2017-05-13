/**
 * ��ά�ǽṹ������
 * ����Ϊ�����ⷨ���ṩ����ӿ�
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
 * ��ά�ǽṹ������
 * ����Ϊ�����ⷨ���ṩ����ӿڣ�ͬʱ�ṩ����Ԫ����ϸ��Ϣ
 * ������Ҫ��������Ϣ�� ���񶥵㣬����ıߣ�����������Լ���ϸ����������Ϣ���ھӣ��������ȣ�
 */
class CUnstructuredGrid {
	// ����
	public:
		MPI_Status status;
		MPI_Request request;

		string config_file;				/**< ���������ļ� */
		string output_filename;			/**< ��������ļ� */

		vector<CVertice2D> vertice;		/**< ���񶥵� */
		vector<CEdge> edge;				/**< ����� */
		vector<int> tri_vertice;		/**< �����εĶ������ */
		vector<int> tri_edge;			/**< �����εı���� */
		vector<int> elem_index;
		vector<int> elem_location;

		vector<int>* local_innerBoundary_index;
		vector<int>* neigh_innerBoundary_index;
		vector<int>* innerBoundary_edge;
		//����Ϣ��Ҫ���͵������ڵ�
		//int** innerBoundary_index;

		// ��������Ҫ�����Ƶ�GPU��
		int *tri_neighbour;				/**< �������ھ���Ϣ */
		int *tri_sharedEdge;			/**< �����ι������Ϣ */
		int *tri_flag;					/**< �����α�ǣ��߽���������Ϣ�� */

		int *local_tri_neighbour;
		int *local_tri_sharedEdge;
		int *local_tri_flag;
		
		CTriangleInfo triangle_infos;	/**< �����ζ�����Ϣ */
		CTriangleInfo local_triangle_infos;
		
	protected:
		int nprocs;
		int area_index;
		int _vertice_num;				/**< ������Ŀ */
		int _edge_num;					/**< ����Ŀ */
		int _triangle_num;				/**< ��������Ŀ */
		int _ghost_triangle_num;		/**< �������ǵ�Ԫ��Ŀ */
		int _cell_num;					/**< �������ܵ�Ԫ��Ŀ�������������� */

		int _local_vertice_num;
		int _local_edge_num;
		int _local_triangle_num;
		int _local_ghost_triangle_num;
		int _local_cell_num;
		
	private:
	
	// ����
	public:
		/** Ĭ�Ϲ��캯�� */
		CUnstructuredGrid();
		
		/**
		 * ����������������Ŀ, ���ĳ����Ԫ�ھӱ��Ϊ-1������Ҫ������������
		 * @param[in] trineigh_filename string �������ھӱ���ļ�
		 */
		void parseGhostNum(const string& trineigh_filename);

		/** ��ʼ���������� */
		void initializeGrid(int gpu_index, int nprocs);
		
		/** ��ʼ������������ */
		void initializeTriangleInfos(void);
		
		/** 
		 * ���������ζ�����
		 * @param[in] trivertice_file string �����ζ������ļ�
		 */
		void readTriangleVertice(const string& trivertice_file);
	
		/** 
		* �������������α��
		* @param[in] trineigh_file string �������ھӱ���ļ�
		*/
		void readTriangleNeighbour(const string& trineigh_file);
	
		/**
		* ���������α߱��
		* @param triedge_file string �����α�����ļ�
		*/
		void readTriangleEdge(const string& triedge_file);
	
		/** 
		* ���붥������
		* @param[in] vertice_file string ���������ļ�
		*/
		void readVertice(const string& vertice_file);
	
		/** 
		* ����ߵĶ�����
		* @param[in] edge_file �ߵĶ������ļ�
		*/
		void readEdgeVertice(const string& edge_file);

		/**
		* ���������ʷֺ��gpu���������Ԫ 
		*/
		void readMeshPartitionInfo(const string& partition_file, int gpu_index);

		void markInnerBoundary();

		int getLocalMeshInfo(int info_type);

		/**
		* �жϵ�Ԫi�뵥Ԫj�Ƿ��й�����
		*/
		bool hasCommonEdge(int i, int j);
	
		/** 
		* ��ʼ�����ⵥԪ
		*/
		void initGhostGrid(void);

		void initLocalGhostGrid(void);
		
		/** @return integer ��������Ŀ */
		int getTriangleNumber(void) const;

		int getLocalTriangleNumber(void) const;

		/** @return integer ����Ŀ */
		int getEdgeNumber(void) const;

		int getLocalEdgeNumber(void) const;
		
		/** @return integer ������Ŀ */
		int getVerticeNumber(void) const;

		int getLocalVerticeNumber(void) const;
		
		/** @return integer ������������Ŀ */
		int getGhostTriangleNumber(void) const;

		int getLocalGhostTriangleNumber(void) const;

		/** @return integer �����ܵ�Ԫ��Ŀ */
		int getCellNumber(void) const;

		int getLocalCellNumber(void) const;

		/** ��ʼ����Ԫ֮�乲��� */
		void initSharedEdge(void);

		void initLocalSharedEdge(void);

		/** �������߽����� */
		void markBoundaryTriangles(void);

		void markLocalBoundaryTriangles(void);
		
		/** ������� */
		void outputGrid(void) const;

		/** ��������Ԫ�ǲ��Ƕ�����ʱ���� */
		void testTrianglesAntiwise(void) const;
		
		/**
		 * �����������(������������)
		 * @param[in] filename string ������������ļ�
		 */
		void outputGridWithGhostCells(const string& filename) const;
		
		/** �������� */
		~CUnstructuredGrid();
	protected:
	
	private:
};