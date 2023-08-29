#include"mesh.h"
#include<algorithm>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<Eigen/Dense>
#include<Eigen/LU>
#include <Eigen/Sparse>
#include <functional>
#include <windows.h> 
#include<boost/filesystem.hpp>




#pragma once
class parameterization{
public:
	parameterization();
	~parameterization();
	void Parameterization(mesh& mesh);
private:
	struct Triplet {
		int row;
		int col;
		double value;
		Triplet(int row, int col, double value)
			: row(row), col(col), value(value) {}
		bool operator<(const Triplet& rhs) const;
	};
	struct Key {
		int row;
		int col;
		/*size_t operator()(const Key& key) const {
			size_t h1 = std::hash<int>{}(key.row);
			size_t h2 = std::hash<int>{}(key.col);
			return h1 ^ (h2 << 1);
		}*/

		size_t operator()(const Key& key) const {
			size_t h1 = std::hash<int>{}(key.row);
			size_t h2 = std::hash<int>{}(key.col);
			return h1 ^ (h2 << 1);
		}
		bool operator==(const Key& rhs) const {
			return row == rhs.row && col == rhs.col;
		}
	};


	bool has_boundary(mesh& mesh);
	void setup_lscm_boundary(mesh& mesh, int& v1, int& v2);
	void count_allface_w(mesh& mesh);
	void new_idx_for_freepoint(mesh& mesh, std::vector<int>& free_vertexs);
	double** create_matrix_solver(mesh& mesh, std::vector<int>free_vertexs, int v1, int v2, int m, int n);
	double norm(double x, double y, double z);
	std::vector<double> normalize(double x, double y, double z);
	std::vector<double> cross(double x1, double y1, double z1, double x2, double y2, double z2);
	double dot(double x1, double y1, double z1, double x2, double y2, double z2);
	double count_area(double x1, double y1, double x2, double y2, double x3, double y3);
	double** matrix_multiplication(double** A, double** B, int a, int b, int c);
	double** matrix_transpose(double** A, int m, int n);
	double LU_decomposition_det(double** A, int n);
	double** Principal_component_selection(double** A, int n);
	double** LU_decomposition_inverse(double** A, int n);
	Eigen::MatrixXd eigen_create_matrix_solver(mesh& mesh, std::vector<int>free_vertexs, int v1, int v2, int m, int n);
	void show_matrix(Eigen::MatrixXd A, int m, int n);
	std::vector<double> sparse_matrix_solving(mesh& mesh, std::vector<int>free_vertexs, int v1, int v2, int m, int n);
	std::vector<double> sparse_matrix_solving_vector(mesh& mesh, std::vector<int>free_vertexs, int v1, int v2, int m, int n);
	void tripletToCSC(const std::vector<Triplet>& triplets, int n_rows, int n_cols, std::vector<int>& column_pointers, std::vector<int>& row_indices, std::vector<double>& nonzero_values);


	
};

