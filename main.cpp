#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <iterator>
#include"mesh.h"
#include"io.h"
#include"parameterization.h"

int main() {
	io io;
    mesh mesh;
    try
    {
        mesh = io.read("C:\\Users\\ÑîºÀ\\Desktop\\2000.off");
    }
    catch (const std::exception& e)
    {
        std::cout << "Error: " << e.what() << std::endl;
        throw;
    }
    parameterization parameterization;
    parameterization.Parameterization(mesh);
    io.write("C:\\Users\\ÑîºÀ\\Desktop\\nefertiti_result.off", mesh);      //nefertiti  merged_submaps_repair
    return 0;
}