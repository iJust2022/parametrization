#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <iterator>
#include "mesh.h"
#include "parameterization.h"
#include "ReadandWrite.h"

int main() {
	ReadandWrite ReadandWrite;
    mesh mesh;
    try
    {
        mesh = ReadandWrite.read("C:\\Users\\ÑîºÀ\\Desktop\\2000.off");
    }
    catch (const std::exception& e)
    {
        std::cout << "Error: " << e.what() << std::endl;
        throw;
    }
    parameterization parameterization;
    parameterization.Parameterization(mesh);
    ReadandWrite.write("C:\\Users\\ÑîºÀ\\Desktop\\2000_result.off", mesh);      //nefertiti  merged_submaps_repair
    return 0;
}