#include "parameterization.h"

parameterization::parameterization() {

}

parameterization::~parameterization(){

}

void parameterization::Parameterization(mesh& mesh){
	if (!has_boundary(mesh)) {
		auto what = "Mesh has no boundary";
		throw what;
	}
	int v1, v2;
	setup_lscm_boundary(mesh,v1,v2);
	count_allface_w(mesh);
	std::vector<int> free_vertexs;
	free_vertexs.resize(mesh.vfh_size("vertex"), 0);
	new_idx_for_freepoint(mesh, free_vertexs);
	int m = mesh.vfh_size("face");
	int n = mesh.vfh_size("vertex") - 2;
	/*Eigen::MatrixXd x = eigen_create_matrix_solver(mesh, free_vertexs, v1, v2, m, n);
	for (int i = 0; i < mesh.vfh_size("vertex"); i++) {
		if (mesh.get_locked(i) == false) {
			int i_idx = free_vertexs[i];
			mesh.set_tex(i, x(i_idx, 0), x(i_idx + n, 0));
		}
	}*/
	double** x = create_matrix_solver(mesh, free_vertexs, v1, v2, m, n);
	for (int i = 0; i < mesh.vfh_size("vertex"); i++) {
		if (mesh.get_locked(i) == false) {
			int i_idx = free_vertexs[i];
			mesh.set_tex(i, x[i_idx][0], x[i_idx + n][0]);
		}
	}
	double maxu = 1.0, maxv = 1.0, minu = 0.0, minv = 0.0;
	for (int i = 0; i < mesh.vfh_size("vertex"); i++) {
		maxu = std::max(maxu, mesh.get_tex(i, "u"));
		maxv = std::max(maxv, mesh.get_tex(i, "v"));
		minu = std::min(minu, mesh.get_tex(i, "u"));
		minv = std::min(minv, mesh.get_tex(i, "v"));
	}
	maxu = maxu - minu;
	maxv = maxv - minv;
	double s = std::max(maxu, maxv);
	for (int i = 0; i < mesh.vfh_size("vertex"); i++) {
		mesh.set_tex(i, (mesh.get_tex(i, "u") - minu) / s, (mesh.get_tex(i, "v") - minv) / s);
	}
	return;
}

bool parameterization::has_boundary(mesh& mesh){
	for (int i = 0; i < mesh.vfh_size("vertex"); i++) {
		if (mesh.is_boundary(mesh.out_halfedge(i))) {
			return true;
		}
	}
	return false;
}

void parameterization::setup_lscm_boundary(mesh& mesh, int& v1, int& v2){
	std::vector<int> boundary;
	for (int i = 0; i < mesh.vfh_size("vertex"); i++) {
		if (mesh.is_boundary(mesh.out_halfedge(i))) {
			boundary.push_back(i);
		}
	}
	double d_max = 0;
	for (int vv1 = 0; vv1 < boundary.size(); vv1++) {
		for (int vv2 = 0; vv2 < boundary.size(); vv2++) {
			double d = mesh.distance(vv1, vv2);
			if (d > d_max) {
				d_max = d;
				v1 = vv1;
				v2 = vv2;
			}
		}
	}
	mesh.set_tex(v1, 0.0, 0.0);
	mesh.set_tex(v2, 1.0, 1.0);
	mesh.set_locked(v1, v2);
}

void parameterization::count_allface_w(mesh& mesh){
	for (int i = 0; i < mesh.vfh_size("face"); i++) {
		int a, b, c;
		a = mesh.get_face_p(i, 0);
		b = mesh.get_face_p(i, 1);
		c = mesh.get_face_p(i, 2);
		mesh.set_pos2d(a, 0.0, 0.0);
		mesh.set_pos2d(b, norm((mesh.get_pos(b, "x") - mesh.get_pos(a, "x")), (mesh.get_pos(b, "y") - mesh.get_pos(a, "y")), (mesh.get_pos(b, "z") - mesh.get_pos(a, "z"))), 0.0);
		std::vector<double> pos0 = normalize((mesh.get_pos(b, "x") - mesh.get_pos(a, "x")), (mesh.get_pos(b, "y") - mesh.get_pos(a, "y")), (mesh.get_pos(b, "z") - mesh.get_pos(a, "z")));
		double c2dx = dot((mesh.get_pos(c, "x") - mesh.get_pos(a, "x")), (mesh.get_pos(c, "y") - mesh.get_pos(a, "y")), (mesh.get_pos(c, "z") - mesh.get_pos(a, "z")), pos0[0], pos0[1], pos0[2]);
		std::vector<double> pos1= normalize((mesh.get_pos(c, "x") - mesh.get_pos(b, "x")), (mesh.get_pos(c, "y") - mesh.get_pos(b, "y")), (mesh.get_pos(c, "z") - mesh.get_pos(b, "z")));
		std::vector<double> pos2 = normalize((mesh.get_pos(a, "x") - mesh.get_pos(b, "x")), (mesh.get_pos(a, "y") - mesh.get_pos(b, "y")), (mesh.get_pos(a, "z") - mesh.get_pos(b, "z")));
		pos1 = cross(pos1[0], pos1[1], pos1[2], pos2[0], pos2[1], pos2[2]);
		pos2 = cross(pos1[0], pos1[1], pos1[2], pos0[0], pos0[1], pos0[2]);
		pos1 = normalize(pos2[0], pos2[1], pos2[2]);
		double c2dy= dot((mesh.get_pos(c, "x") - mesh.get_pos(a, "x")), (mesh.get_pos(c, "y") - mesh.get_pos(a, "y")), (mesh.get_pos(c, "z") - mesh.get_pos(a, "z")), pos1[0], pos1[1], pos1[2]);
		mesh.set_pos2d(c, c2dx, c2dy);
		double area = count_area(mesh.get_pos2d(a, "x"), mesh.get_pos2d(a, "y"), mesh.get_pos2d(b, "x"), mesh.get_pos2d(b, "y"), mesh.get_pos2d(c, "x"), mesh.get_pos2d(c, "y"));
		mesh.set_face_W(i, 0, "r", (mesh.get_pos2d(c, "x") - mesh.get_pos2d(b, "x")) / area);
		mesh.set_face_W(i, 0, "i", (mesh.get_pos2d(c, "y") - mesh.get_pos2d(b, "y")) / area);
		mesh.set_face_W(i, 1, "r", (mesh.get_pos2d(a, "x") - mesh.get_pos2d(c, "x")) / area);
		mesh.set_face_W(i, 1, "i", (mesh.get_pos2d(a, "y") - mesh.get_pos2d(c, "y")) / area);
		mesh.set_face_W(i, 2, "r", (mesh.get_pos2d(b, "x") - mesh.get_pos2d(a, "x")) / area);
		mesh.set_face_W(i, 2, "i", (mesh.get_pos2d(b, "y") - mesh.get_pos2d(a, "y")) / area);
	}
}

void parameterization::new_idx_for_freepoint(mesh& mesh, std::vector<int>& free_vertexs){
	int j = 0;
	for (int i = 0; i < mesh.vfh_size("vertex"); i++) {
		if (mesh.get_locked(i) == false) {
			free_vertexs[i] = j;
			j++;
		}
	}
}

double** parameterization::create_matrix_solver(mesh& mesh, std::vector<int> free_vertexs, int v1, int v2, int m, int n){
	double** A = new double* [2 * m];
	for (int i = 0; i < 2*m; i++) {
		A[i] = new double[2 * n];
		for (int j = 0; j < 2 * n; j++) {
			A[i][j] = 0;
		}
	}
	double** B = new double* [2 * m];
	for (int i = 0; i < 2 * m; i++) {
		B[i] = new double[4];
		for (int j = 0; j < 4; j++) {
			B[i][j] = 0;
		}
	}
	double** b = new double* [2 * m];
	for (int i = 0; i < 2 * m; i++) {
		b[i] = new double[1];
		for (int j = 0; j < 1; j++) {
			b[i][j] = 0;
		}
	}
	double** Up = new double* [4];
	for (int i = 0; i < 4; i++) {
		Up[i] = new double[1];
		for (int j = 0; j < 1; j++) {
			Up[i][j] = 0;
		}
	}
	for (int i = 0; i < mesh.vfh_size("face"); i++) {
		for (int j = 0; j < 3; j++) {
			int k = mesh.get_face_p(i, j);
			int k_id = free_vertexs[k];
			if (mesh.get_locked(k) == false) {
				A[i][k_id] = mesh.get_face_W(i, j, "r");
				A[i + m][k_id] = -mesh.get_face_W(i, j, "i");
				A[i][k_id + n] = mesh.get_face_W(i, j, "i");
				A[i + m][k_id + n] = mesh.get_face_W(i, j, "r");
			}
			else if (mesh.get_locked(k) == true) {
				if (k == v1) {
					B[i][0] = mesh.get_face_W(i, j, "r");
					B[i + m][0] = -mesh.get_face_W(i, j, "i");
					B[i][2] = mesh.get_face_W(i, j, "i");
					B[i + m][2] = mesh.get_face_W(i, j, "r");
				}
				else if (k == v2) {
					B[i][1] = mesh.get_face_W(i, j, "r");
					B[i + m][1] = -mesh.get_face_W(i, j, "i");
					B[i][3] = mesh.get_face_W(i, j, "i");
					B[i + m][3] = mesh.get_face_W(i, j, "r");
				}
			}
		}
	}
	//show_matrix(A, 2 * m, 2 * n);
	//show_matrix(B, 2 * m, 4);
	Up[0][0] = mesh.get_tex(v1, "u");
	Up[1][0] = mesh.get_tex(v2, "u");
	Up[2][0] = mesh.get_tex(v1, "v");
	Up[3][0] = mesh.get_tex(v2, "v");
	b = matrix_multiplication(B, Up, 2 * m, 4, 1);
	for (int i = 0; i < 2 * m; i++) {
		for(int j = 0; j < 1; j++) {
			b[i][j] = -b[i][j];
		}
	}
	double** tA = matrix_transpose(A, 2 * m, 2 * n);
	double** tAA = matrix_multiplication(tA, A, 2 * n, 2 * m, 2 * n);
	double determinant = LU_decomposition_det(tAA, 2 * n);
	if (determinant == 0) {
		auto what = "Error: Matrix is not reversible";
		throw what;
	}
	double** P = Principal_component_selection(tAA, 2 * n);
	double** tAA_inv = LU_decomposition_inverse(tAA, 2 * n);
	double** tAb = matrix_multiplication(tA, b, 2 * n, 2 * m, 1);
	double** x = matrix_multiplication(tAA_inv, tAb, 2 * n, 2 * n, 1);
	
	x = matrix_multiplication(P, x, 2 * n, 2 * n, 1);
	//show_matrix(x, 2 * n, 1);
	return x;
}

double parameterization::norm(double x, double y, double z){
	double s;
	s = x * x + y * y + z * z;
	s = sqrt(s);
	return s;
}

std::vector<double> parameterization::normalize(double x, double y, double z){
	std::vector<double> newpos;
	newpos.resize(3);
	double n = norm(x, y, z);
	n = 1.0 / n;
	newpos[0] = x * n;
	newpos[1] = y * n;
	newpos[2] = z * n;
	return newpos;
}

std::vector<double> parameterization::cross(double x1, double y1, double z1, double x2, double y2, double z2){
	std::vector<double> newpos;
	newpos.resize(3);
	newpos[0] = y1 * z2 - z1 * y2;
	newpos[1] = z1 * x2 - x1 * z2;
	newpos[2] = x1 * y2 - y1 * x2;
	return newpos;
}

double parameterization::dot(double x1, double y1, double z1, double x2, double y2, double z2){
	double p;
	p = x1 * x2 + y1 * y2 + z1 * z2;
	return p;
}

double parameterization::count_area(double x1, double y1, double x2, double y2, double x3, double y3){
	double s;
	s = (x1 * y2 - x2 * y1) + (x2 * y3 - x3 * y2) + (y1 * x3 - y3 * x1);
	return s;
}

double** parameterization::matrix_multiplication(double** A, double** B, int a, int b, int c){
	double** M = new double* [a];
	for (int i = 0; i < a; i++) {
		M[i] = new double[c];
		for (int j = 0; j < c; j++) {
			M[i][j] = 0.0;
			for (int k = 0; k < b; k++) {
				M[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return M;
}

double** parameterization::matrix_transpose(double** A, int m, int n){
	double** tA = new double* [n];   //��A��ת�þ���
	for (int i = 0; i < n; i++) {
		tA[i] = new double[m];
		for (int j = 0; j < m; j++) {
			tA[i][j] = A[j][i];
		}
	}
	return tA;
}

double parameterization::LU_decomposition_det(double** A, int n){
	double determinant = 1.0;
	double determinantL = 1.0;
	double determinantU = 1.0;
	double** L = new double* [n];
	for (int i = 0; i < n; i++) {
		L[i] = new double[n];
		for (int j = 0; j < n; j++) {
			L[i][j] = 0.0;
		}
	}
	double** U = new double* [n];
	for (int i = 0; i < n; i++) {
		U[i] = new double[n];
		for (int j = 0; j < n; j++) {
			U[i][j] = 0.0;
		}
	}
	for (int i = 0; i < n; i++) {
		L[i][i] = 1.0;
	}
	for (int j = 0; j < n; j++) {
		U[0][j] = A[0][j];
	}
	for (int i = 1; i < n; i++) {
		L[i][0] = A[i][0] / U[0][0];
	}
	for (int k = 1; k < n; k++) {
		for (int j = k; j < n; j++) {
			double s = 0.0;
			for (int t = 0; t < k; t++) {
				s += L[k][t] * U[t][j];
			}
			U[k][j] = A[k][j] - s;
		}
		for (int i = k; i < n; i++) {
			double s = 0.0;
			for (int t = 0; t < k; t++) {
				s += L[i][t] * U[t][k];
			}
			L[i][k] = (A[i][k] - s) / U[k][k];
		}
	}
	for (int i = 0; i < n; i++) {
		determinantL = determinantL * L[i][i];
		determinantU = determinantU * U[i][i];
	}
	determinant = determinantL * determinantU;
	return determinant;
}

double** parameterization::Principal_component_selection(double** A, int n){
	double** P = new double* [n];
	for (int i = 0; i < n; i++) {
		P[i] = new double[n];
		for (int j = 0; j < n; j++) {
			if (i == j) {
				P[i][j] = 1;
			}
			else {
				P[i][j] = 0;
			}
		}
	}
	for (int k = 0; k < n; k++) {
		double max = 0.0;
		int k_ = k;
		for (int i = k; i < n; i++) {
			if (fabs(A[i][k] > max)) {
				max = fabs(A[i][k]);
				k_ = i;
			}
		}
		if (max == 0) {
			auto what = "Error: singular  matrix";
			throw what;
		}
		for (int j = 0; j < n; ++j) {
			double temp = A[k][j];
			A[k][j] = A[k_][j];
			A[k_][j] = temp;
			double p = P[k][j];
			P[k][j] = P[k_][j];
			P[k_][j] = p;
		}
	}
	for (int k = 0; k < n; k++) {
		double max = 0.0;
		int k_ = k;
		for (int j = k; j < n; j++) {
			if (fabs(A[k][j] > max)) {
				max = fabs(A[k][j]);
				k_ = j;
			}
		}
		for (int i = 0; i < n; i++) {
			double temp = A[i][k];
			A[i][k] = A[i][k_];
			A[i][k_] = temp;
			double p = P[i][k];
			P[i][k] = P[i][k_];
			P[i][k_] = p;
		}
	}
	return P;
}

double** parameterization::LU_decomposition_inverse(double** A, int n){
	double** L = new double* [n];
	for (int i = 0; i < n; i++) {
		L[i] = new double[n];
		for (int j = 0; j < n; j++) {
			L[i][j] = 0.0;
		}
	}
	double** U = new double* [n];
	for (int i = 0; i < n; i++) {
		U[i] = new double[n];
		for (int j = 0; j < n; j++) {
			U[i][j] = 0.0;
		}
	}
	double** L_inv = new double* [n];
	for (int i = 0; i < n; i++) {
		L_inv[i] = new double[n];
		for (int j = 0; j < n; j++) {
			L_inv[i][j] = 0.0;
		}
	}
	double** U_inv = new double* [n];
	for (int i = 0; i < n; i++) {
		U_inv[i] = new double[n];
		for (int j = 0; j < n; j++) {
			U_inv[i][j] = 0.0;
		}
	}
	for (int i = 0; i < n; i++) {
		L[i][i] = 1.0;
	}
	for (int j = 0; j < n; j++) {
		U[0][j] = A[0][j];
	}
	for (int i = 1; i < n; i++) {
		L[i][0] = A[i][0] / U[0][0];
	}
	for (int k = 1; k < n; k++) {
		for (int j = k; j < n; j++) {
			double s = 0.0;
			for (int t = 0; t < k; t++) {
				s += L[k][t] * U[t][j];
			}
			U[k][j] = A[k][j] - s;
		}
		for (int i = k; i < n; i++) {
			double s = 0.0;
			for (int t = 0; t < k; t++) {
				s += L[i][t] * U[t][k];
			}
			L[i][k] = (A[i][k] - s) / U[k][k];
		}
	}
	for (int j = 0; j < n; j++) {
		for (int i = j; i < n; i++) {
			if (i == j) {
				L_inv[i][j] = 1 / L[i][j];
			}
			else if (i < j) {
				L_inv[i][j] = 0.0;
			}
			else {
				double s = 0.0;
				for (int k = j; k < i; k++) {
					s += L[i][k] * L_inv[k][j];
				}
				L_inv[i][j] = -L_inv[j][j] * s;
			}
		}
	}
	for (int j = 0; j < n; j++) {
		for (int i = j; i >= 0; i--) {
			if (i == j) {
				U_inv[i][j] = 1 / U[i][j];
			}
			else if (i > j) {
				U_inv[i][j] = 0.0;
			}
			else {
				double s = 0.0;
				for (int k = i + 1; k <= j; k++) {
					s += U[i][k] * U_inv[k][j];
				}
				U_inv[i][j] = -1 / U[i][i] * s;
			}
		}
	}
	double** A_inv = matrix_multiplication(U_inv, L_inv, n, n, n);
	return A_inv;
}

Eigen::MatrixXd parameterization::eigen_create_matrix_solver(mesh& mesh, std::vector<int> free_vertexs, int v1, int v2, int m, int n){
	Eigen::MatrixXd A(2 * m, 2 * n);
	Eigen::MatrixXd B(2 * m, 4);
	Eigen::MatrixXd Up(4, 1);
	for (int i = 0; i < 2 * m; i++) {
		for (int j = 0; j < 2 * n; j++) {
			A(i, j) = 0;
		}
	}
	for (int i = 0; i < 2 * m; i++) {
		for (int j = 0; j < 4; j++) {
			B(i, j) = 0;
		}
	}
	for (int i = 0; i < mesh.vfh_size("face"); i++) {
		for (int j = 0; j < 3; j++) {
			int k = mesh.get_face_p(i, j);
			if (mesh.get_locked(k) == false) {
				int k_id = free_vertexs[k];
				A(i, k_id) = mesh.get_face_W(i, j, "r");
				A(i + m, k_id) = -mesh.get_face_W(i, j, "i");
				A(i, k_id + n) = mesh.get_face_W(i, j, "i");
				A(i + m, k_id + n) = mesh.get_face_W(i, j, "r");
			}
			else if (mesh.get_locked(k) == true) {
				if (k == v1) {
					B(i, 0) = mesh.get_face_W(i, j, "r");
					B(i + m, 0) = -mesh.get_face_W(i, j, "i");
					B(i, 2) = mesh.get_face_W(i, j, "i");
					B(i + m, 2) = mesh.get_face_W(i, j, "r");
				}
				else if (k == v2) {
					B(i, 1) = mesh.get_face_W(i, j, "r");
					B(i + m, 1) = -mesh.get_face_W(i, j, "i");
					B(i, 3) = mesh.get_face_W(i, j, "i");
					B(i + m, 3) = mesh.get_face_W(i, j, "r");
				}
			}
		}
	}
	/*for (int i = 0; i < 2 * m; i++) {
		for (int j = 0; j < 2 * n; j++) {
			std::cout << A(i, j) << " ";
		}
		std::cout << std::endl;
	}*/
	/*for (int i = 0; i < 2 * m; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << B(i, j) << " ";
		}
		std::cout << std::endl;
	}*/
	Up(0, 0) = mesh.get_tex(v1, "u");
	Up(1, 0) = mesh.get_tex(v2, "u");
	Up(2, 0) = mesh.get_tex(v1, "v");
	Up(3, 0) = mesh.get_tex(v2, "v");
	Eigen::MatrixXd b = B * Up;
	for (int i = 0; i < 2 * m; i++) {
		for (int j = 0; j < 1; j++) {
			b(i, j) = -b(i, j);
		}
	}
	Eigen::MatrixXd x = (A.transpose() * A).inverse() * (A.transpose() * b);
	/*for (int i = 0; i < 2 * m; i++) {
		for (int j = 0; j < 1; j++) {
			std::cout << x(i, j) << " ";
		}
		std::cout << std::endl;
	}*/
	return x;
}

void parameterization::show_matrix(double** A, int m, int n){
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
}


