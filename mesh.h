#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <iterator>
#include<assert.h>


#pragma once
class mesh{
public:
	mesh();
	~mesh();
	void vfh_reserve(int vn, int fn, int hn);
	int vfh_capacity(std::string name);
	int vfh_size(std::string name);
	void add_vert_property(std::vector<double> property_i, const bool has_normals, const bool has_texcoords, const bool has_colors);
	void add_face(std::vector<int> property_facei);
	bool is_boundary(int h);
	bool is_valid(int i);
	int get_face(int h);
	int out_halfedge(int v);
	int find_halfedge(int start, int end);
	int to_vertex(int i);
	int opposite_halfedge(int i);
	int next_halfedge(int i);
	int prev_halfedge(int i);
	int new_edge(int start, int end);
	int new_face(std::vector<int> face_point);
	double distance(int vv1, int vv2);
	void set_tex(int v, double ui, double vi);
	double get_tex(int v, std::string s);
	void set_locked(int v1, int v2);
	int get_face_halfedge(int f);
	int get_face_p(int i,int j);
	bool get_locked(int i);
	void set_pos2d(int v, double xi, double yi);
	double get_pos2d(int v, std::string s);
	double get_pos(int v, std::string s);
	void set_face_W(int f, int i, std::string s,double w);
	double get_face_W(int f, int i, std::string s);
	
	
private:
	struct Pos {
		double x;
		double y;
		double z;
	};
	struct Nor {
		double nx;
		double ny;
		double nz;
	};
	struct Col {
		double r;
		double g;
		double b;
	};
	struct Tex {
		double u;
		double v;
	};
	struct Pos2d {
		double x;
		double y;
	};
	struct W {
		double r;
		double i;
	};
	struct Vertex {
		int halfedge = -1;  //ÿ������Ĵ�����ߣ��߽綥��߽���
	};
	struct Halfedge {
		int toVertex = -1;  //���������ߵ���
		int face = -1;  //���������ߵ���
		int prevHalfedge = -1;  //��һ�����
		int nextHalfedge = -1;  //��һ�����
	};
	struct Face {
		int halfedge = -1;  //�������һ�����
	};
	struct VertexProperty {
		Vertex vertex;
		Pos pos;  //ÿ�����������
		Nor nor;  //ÿ����ķ�����Ϣ
		Col col;  //ÿ���������ɫ��Ϣ
		Tex tex;  //ÿ������Ĳ�������
		Pos2d pos2d;  //ÿ������ľֲ�����
		bool locked = false;
	};
	struct HalfedgeProperty {
		Halfedge halfedge;  //���������ߵ���
	};

	struct FaceProperty {
		Face face;
		int p[3]; //���������
		W W[3];  //ÿ�����w��ʵ�����ֺ���������
	};
	
	void set_next_halfedge(int h, int nh);
	void set_halfedge_to_vertex(int h, int v);
	void set_face_halfedge(int f, int h);
	void set_point_out_halfedge(int v, int h);
	void set_halfedge_face(int h, int f);
	void adjust_outgoing_halfedge(int v);

	std::vector<VertexProperty> Vertexs;      //���е�
	std::vector<HalfedgeProperty> Halfedges;  //���а��
	std::vector<FaceProperty> Faces;          //������
};

