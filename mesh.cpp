#include "mesh.h"
#pragma once


mesh::mesh(){

}

mesh::~mesh(){

}

void mesh::vfh_reserve(int nv, int nf, int nh) {
	Vertexs.reserve(nv);
	Faces.reserve(nf);
	Halfedges.reserve(nh);
}

int mesh::vfh_capacity(std::string name) {
	if (name == "vertex") {
		return Vertexs.capacity();
	}
	else if (name == "face") {
		return Faces.capacity();
	}
	else if (name == "halfedge") {
		return Halfedges.capacity();
	}
}

int mesh::vfh_size(std::string name){
	if (name == "vertex") {
		return Vertexs.size();
	}
	else if (name == "face") {
		return Faces.size();
	}
	else if (name == "halfedge") {
		return Halfedges.size();
	}
}

void mesh::add_vert_property(std::vector<double> property_i, const bool has_normals, const bool has_texcoords, const bool has_colors) {
	if (property_i.size() < 3) {
		auto what = "data error ! ";
		throw what;
	}
	VertexProperty vert;
	int c = 0;
	vert.pos.x = property_i[c];
	c++;
	vert.pos.y = property_i[c];
	c++;
	vert.pos.z = property_i[c];
	c++;
	if (has_normals) {
		vert.nor.nx = property_i[c];
		c++;
		vert.nor.ny = property_i[c];
		c++;
		vert.nor.nz = property_i[c];
		c++;
	}
	if (has_colors) {
		double R = property_i[c];
		c++;
		double G = property_i[c];
		c++;
		double B = property_i[c];
		c++;
		if (R > 1.0f || G > 1.0f || B > 1.0f){
			R /= 255.0f;
			G /= 255.0f;
			B /= 255.0f;
		}
		vert.col.r = R;
		vert.col.g = G;
		vert.col.b = B;

	}
	if (has_texcoords) {
		vert.tex.u = property_i[c];
		c++;
		vert.tex.v = property_i[c];
	}
	Vertexs.push_back(vert);
	return;
}

void mesh::add_face(std::vector<int> property_facei) {
	int c = 0;
	const int n = property_facei[c];
	if (n != 3) {
		auto what = "Not a triangular mesh";
		throw what;
	}
	c++;
	std::vector<int>face_point;
	face_point.resize(n);
	for (int j = 0; j < face_point.size(); j++) {
		face_point[j] = property_facei[c];
		c++;
	}
	int i, ii;
	std::vector<int> halfedges;
	std::vector<bool> is_new;
	std::vector<bool> needs_adjust;
	halfedges.resize(n);
	is_new.resize(n);
	needs_adjust.resize(n, false);
	for (i = 0, ii = 1; i < n; i++, ii++, ii %= n) {
		if (!is_boundary(out_halfedge(face_point[i]))) {
			auto what = "SurfaceMesh::add_face: Complex vertex";
			throw what;
		}
		halfedges[i] = find_halfedge(face_point[i], face_point[ii]);
		is_new[i] = !is_valid(halfedges[i]);
		if (!is_new[i] && !is_boundary(halfedges[i])){
			auto what = "SurfaceMesh::add_face: Complex edge.";
			throw what;
		}
	}
	int inner_next, inner_prev, outer_next, outer_prev, boundary_next, boundary_prev, patch_start, patch_end;
	for (i = 0, ii = 1; i < n; i++, ii++, ii %= n){
		if (!is_new[i] && !is_new[ii]){
			inner_prev = halfedges[i];
			inner_next = halfedges[ii];
			if (next_halfedge(inner_prev) != inner_next){
				outer_prev = opposite_halfedge(inner_next);
				outer_next = opposite_halfedge(inner_prev);
				boundary_prev = outer_prev;
				do
				{
					boundary_prev =opposite_halfedge(next_halfedge(boundary_prev));
				} while (!is_boundary(boundary_prev) ||boundary_prev == inner_prev);
				boundary_next = next_halfedge(boundary_prev);
				assert(is_boundary(boundary_prev));
				assert(is_boundary(boundary_next));
				if (boundary_next == inner_next){
					auto what ="SurfaceMesh::add_face: Patch re-linking failed.";
					throw what;
				}
				patch_start = next_halfedge(inner_prev);
				patch_end = prev_halfedge(inner_next);
				set_next_halfedge(boundary_prev, patch_start);
				set_next_halfedge(patch_end, boundary_next);
				set_next_halfedge(inner_prev, inner_next);
			}
		}
	}
	for (i = 0, ii = 1; i < n; ++i, ++ii, ii %= n){
		if (is_new[i]){
			halfedges[i] = new_edge(face_point[i], face_point[ii]);
		}
	}
	int f = new_face(face_point);
	set_face_halfedge(f,halfedges[n-1]);
	for (i = 0, ii = 1; i < n; ++i, ++ii, ii %= n){
		int v = face_point[ii];
		inner_prev = halfedges[i];
		inner_next = halfedges[ii];
		int id = 0;
		if (is_new[i])
			id |= 1;
		if (is_new[ii])
			id |= 2;

		if (id){
			outer_prev = opposite_halfedge(inner_next);
			outer_next = opposite_halfedge(inner_prev);
			switch (id){
			case 1: 
				boundary_prev = prev_halfedge(inner_next);
				set_next_halfedge(boundary_prev, outer_next);
				set_point_out_halfedge(v, outer_next);
				break;

			case 2: 
				boundary_next = next_halfedge(inner_prev);
				set_next_halfedge(outer_prev, boundary_next);
				set_point_out_halfedge(v, boundary_next);
				break;

			case 3: 
				if (!is_valid(out_halfedge(v))){
					set_point_out_halfedge(v, outer_next);
					set_next_halfedge(outer_prev, outer_next);
				}
				else {
					boundary_next = out_halfedge(v);
					boundary_prev = prev_halfedge(boundary_next);
					set_next_halfedge(boundary_prev, outer_next);
					set_next_halfedge(outer_prev, boundary_next);
				}
				break;
			}
			set_next_halfedge(inner_prev, inner_next);
		}
		else {
			needs_adjust[ii] = (out_halfedge(v) == inner_next);
		}	
		set_halfedge_face(halfedges[i], f);
	}
	for (i = 0; i < n; ++i){
		if (needs_adjust[i]){
			adjust_outgoing_halfedge(face_point[i]);
		}
	}
}

bool mesh::is_boundary(int h) {
	return (!(is_valid(h) && is_valid(get_face(h))));
}

bool mesh::is_valid(int i) {
	return (i != -1);
}

int mesh::get_face(int i) {
	return Halfedges[i].halfedge.face;
}

int mesh::out_halfedge(int i) {
	return Vertexs[i].vertex.halfedge;
}

int mesh::find_halfedge(int start, int end) {
	assert(is_valid(start) && is_valid(end));
	int h = out_halfedge(start);
	const int hh = h;
	if (is_valid(h)){
		do
		{
			if (to_vertex(h) == end)
				return h;
			h = next_halfedge(opposite_halfedge(h));
		} while (h != hh);
	}
	return -1;
}

int mesh::to_vertex(int i) {
	return Halfedges[i].halfedge.toVertex;
}

int mesh::opposite_halfedge(int i) {
	return ((i & 1) ?i - 1 : i + 1);
}

int mesh::next_halfedge(int i) {
	return Halfedges[i].halfedge.nextHalfedge;
}

int mesh::prev_halfedge(int i) {
	return Halfedges[i].halfedge.prevHalfedge;
}

void mesh::set_next_halfedge(int h, int nh){
	Halfedges[h].halfedge.nextHalfedge = nh;
	Halfedges[nh].halfedge.prevHalfedge = h;
}

int mesh::new_edge(int start, int end) {
	assert(start != end);
	HalfedgeProperty h;
	Halfedges.push_back(h);
	Halfedges.push_back(h);
	int h0 = Halfedges.size() - 2;
	int h1 = Halfedges.size() - 1;
	set_halfedge_to_vertex(h0, end);
	set_halfedge_to_vertex(h1, start);
	return h0;
}

void mesh::set_halfedge_to_vertex(int h, int v) {
	Halfedges[h].halfedge.toVertex = v;
}

int mesh::new_face(std::vector<int> face_point) {
	FaceProperty F;
	for (int i = 0; i < face_point.size(); i++) {
		F.p[i] = face_point[i];
	}
	Faces.push_back(F);
	int f = Faces.size() - 1;
	return f;
}

double mesh::distance(int vv1, int vv2){
	double d;
	double dx = (Vertexs[vv1].pos.x - Vertexs[vv2].pos.x) * (Vertexs[vv1].pos.x - Vertexs[vv2].pos.x);
	double dy = (Vertexs[vv1].pos.y - Vertexs[vv2].pos.y) * (Vertexs[vv1].pos.y - Vertexs[vv2].pos.y);
	double dz = (Vertexs[vv1].pos.z - Vertexs[vv2].pos.z) * (Vertexs[vv1].pos.z - Vertexs[vv2].pos.z);
	d = sqrt(dx + dy + dz);
	return d;
}

void mesh::set_tex(int v, double ui, double vi){
	Vertexs[v].tex.u = ui;
	Vertexs[v].tex.v = vi;
}

double mesh::get_tex(int v, std::string s){
	if (s == "u") {
		return Vertexs[v].tex.u;
	}
	else if (s == "v") {
		return Vertexs[v].tex.v;
	}
}

void mesh::set_locked(int v1, int v2){
	Vertexs[v1].locked = true;
	Vertexs[v2].locked = true;
}

int mesh::get_face_halfedge(int f){
	return Faces[f].face.halfedge;
}

int mesh::get_face_p(int i, int j){
	return Faces[i].p[j];
}

bool mesh::get_locked(int i){
	return Vertexs[i].locked;
}

void mesh::set_pos2d(int v, double xi, double yi){
	Vertexs[v].pos2d.x = xi;
	Vertexs[v].pos2d.y = yi;
}

double mesh::get_pos2d(int v, std::string s) {
	if (s == "x") {
		return Vertexs[v].pos2d.x;
	}
	else if (s == "y") {
		return Vertexs[v].pos2d.y;
	}
}

double mesh::get_pos(int v, std::string s)
{
	if (s == "x") {
		return Vertexs[v].pos.x;
	}
	else if (s == "y") {
		return Vertexs[v].pos.y;
	}
	else if (s == "z") {
		return Vertexs[v].pos.z;
	}
}

void mesh::set_face_W(int f, int i, std::string s,double w){
	if (s == "r") {
		Faces[f].W[i].r = w;
	}
	else if (s == "i") {
		Faces[f].W[i].i = w;
	}
}

double mesh::get_face_W(int f, int i, std::string s){
	if (s == "r") {
		return Faces[f].W[i].r;
	}
	else if (s == "i") {
		return Faces[f].W[i].i;
	}
}

void mesh::set_face_halfedge(int f, int h) {
	Faces[f].face.halfedge = h;
}

void mesh::set_point_out_halfedge(int v, int h) {
	Vertexs[v].vertex.halfedge = h;
}

void mesh::set_halfedge_face(int h, int f) {
	Halfedges[h].halfedge.face = f;
}

void mesh::adjust_outgoing_halfedge(int v) {
	int h = out_halfedge(v);
	const int hh = h;
	if (is_valid(h)){
		do{
			if (is_boundary(h)){
				set_point_out_halfedge(v, h);
				return;
			}
			h = next_halfedge(opposite_halfedge(h));
		} while (h != hh);
	}
}
