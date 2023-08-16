#include "io.h"
io::io(){

}
io::~io(){

}
mesh io::read(std::string file) {
    mesh mesh;
	std::string ext = file.substr(file.find_last_of('.') + 1);
	if (ext == "off") {
		read_off(file, mesh);
	}
	else {
        auto what = "Could not find reader for " + file;
		throw what;
	}
    return mesh;
}
void io::write(std::string file, mesh& mesh){
    std::fstream f;
    f.open(file, std::ios::out);
    f << "OFF" << std::endl;
    f << mesh.vfh_size("vertex") << " " << mesh.vfh_size("face") << " " << "0" << " " << std::endl;
    for (int i = 0; i < mesh.vfh_size("vertex"); i++) {
        f << mesh.get_tex(i, "u") << " " << mesh.get_tex(i, "v") << " " << "0" << " " << std::endl;
    }
    for (int i = 0; i < mesh.vfh_size("face"); i++) {
        f << "3" << " " << mesh.get_face_p(i, 0) << " " << mesh.get_face_p(i, 1) << " " << mesh.get_face_p(i, 2) << " " << std::endl;
    }
    f.close();
    return;
}
void io::read_off(std::string file, mesh& mesh) {
	std::vector<char> val;
	bool has_texcoords = false;
	bool has_normals = false;
	bool has_colors = false;
	bool has_hcoords = false;
	bool has_dim = false;
	bool is_Binary = false;
	std::ifstream in(file);
	std::string line;
	std::getline(in, line);
	std::stringstream ss(line);
	char s;
	while (ss >> s) {
		val.push_back(s);
	}
	int c = 0;
    if (val[c] == 'S' && val[c] == 'T')
    {
        has_texcoords = true;
        c += 2;
    }
    if (val[c] == 'C')
    {
        has_colors = true;
        c++;
    }
    if (val[c] == 'N')
    {
        has_normals = true;
        ++c;
    }
    if (val[c] == '4')
    {
        has_hcoords = true;
        c++;
    }
    if (val[c] == 'n')
    {
        has_dim = true;
        c++;
    }
    /*std::string is_off="";
    for (int i = 0; i < 3; i++) {
        is_off += std::to_string("val[c]");
        c++;
    }
    if (is_off.compare("OFF") != 0)
    {
        in.close();
        auto what = "Failed to parse OFF header";
        throw what;
    }
    char is_binary[6];
    for (int i = 0; i < 6; i++) {
        is_binary[i] = val[c];
        c++;
    }
    if (strcmp(is_binary, "BINARY") == 0)
        is_Binary = true;*/

    if (has_hcoords == true) {
        in.close();
        auto what = "Error: Homogeneous coordinates not supported.";
        throw what;
    }
    if (has_dim == true) {
        in.close();
        auto what = "Error: vertex dimension != 3 not supported.";
        throw what;
    }
    /*if (is_binary) {
        read_off_binary(mesh, in, has_normals, has_texcoords, has_colors);
    }
    else {
        read_off_ascii(mesh, in, has_normals, has_texcoords, has_colors);
    } */
    read_off_ascii(mesh, in, has_normals, has_texcoords, has_colors);
    in.close();
}

void io::read_off_binary(mesh& mesh, std::ifstream& in, const bool has_normals, const bool has_texcoords, const bool has_colors) {
    return;
}
void io::read_off_ascii(mesh& mesh, std::ifstream& in, const bool has_normals, const bool has_texcoords, const bool has_colors) {
    std::vector<int> vfh_size;
    std::string line;
    std::getline(in, line);
    std::stringstream ss(line);
    int number;
    while (ss >> number) {
        vfh_size.push_back(number);
    }
    mesh.vfh_reserve(vfh_size[0], vfh_size[1], 2 * std::max(3 * vfh_size[1], vfh_size[2]));
    for (int i = 0; i < mesh.vfh_capacity("vertex"); i++) {
        std::getline(in, line);
        std::stringstream ss(line);
        double p;
        std::vector<double> property_i;
        while (ss >> p) {
            property_i.push_back(p);
        }
        mesh.add_vert_property(property_i, has_normals, has_texcoords, has_colors);
    }
    for (int i = 0; i < mesh.vfh_capacity("face"); i++) {
        std::getline(in, line);
        std::stringstream ss(line);
        double p;
        std::vector<int> property_facei;
        while (ss >> p) {
            property_facei.push_back(p);
        }
        mesh.add_face(property_facei);
    }
    return;
}