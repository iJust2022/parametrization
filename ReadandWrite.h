#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <iterator>
#include"mesh.h"




class ReadandWrite{
public:
	ReadandWrite();
	~ReadandWrite();
	mesh read(std::string file);
	void write(std::string file, mesh& mesh);
	
private:
	void read_off(std::string file,mesh &mesh);
	void read_off_binary(mesh& mesh, std::ifstream& in, const bool has_normals, const bool has_texcoords, const bool has_colors);
	void read_off_ascii(mesh& mesh, std::ifstream& in, const bool has_normals, const bool has_texcoords, const bool has_colors);
};

