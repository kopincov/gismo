/** @file gsUV_XYZ_2_gismo_xml.cpp

    @brief Converts point cloud UV and XYZ to gismo xml file.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Bracco, C. Giannelli, J. Speh
*/

#include <iostream>
#include <fstream>
#include <vector>

#include <gismo.h>

using namespace gismo;

gsMatrix<> read_parameters(const std::string& filename)
{
    std::fstream stream;
    stream.open(filename.c_str());
    
    real_t u, v;
    std::vector<real_t> vector;
    
    while (stream >> u >> v)
    {
        vector.push_back(u);
        vector.push_back(v);
    }

    int num_cols = vector.size() / 2;
    gsMatrix<> parameters(2, num_cols);
    
    for (int col = 0; col != num_cols; col++)
    {
        parameters(0, col) = vector[2 * col + 0];
        parameters(1, col) = vector[2 * col + 1];
    }
    
    return parameters;
}

gsMatrix<> read_points(const std::string& filename)
{
    std::fstream stream;
    stream.open(filename.c_str());
    
    real_t x, y, z;
    std::vector<real_t> vector;
    
    while (stream >> x >> y >> z)
    {
        vector.push_back(x);
        vector.push_back(y);
        vector.push_back(z);
    }

    int num_cols = vector.size() / 3;
    gsMatrix<> points(3, num_cols);
    
    for (int col = 0; col != num_cols; col++)
    {
        points(0, col) = vector[3 * col + 0];
        points(1, col) = vector[3 * col + 1];
        points(2, col) = vector[3 * col + 2];
    }
    
    return points;
}

int main(int argc, char *argv[])
{
    std::string inputUV("");
    std::string inputXYZ("");
    std::string output("fitting_data");
    
    // ----------------------------------------------------------------------

    gsCmdLine cmd("Scattered fitting using hierarchical splines.");

    cmd.addString("u", "inputUV", "File with u v data", inputUV);
    cmd.addString("x", "inputXYZ", "File with x y z data", inputXYZ);
    cmd.addString("o", "output", "Output file", output);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // ----------------------------------------------------------------------
    
    gsInfo << "----------------------------------------\n"
           << "INPUT DATA:\n\n" 
           << "inputUV:   " << inputUV << "\n\n"
           << "inputXYZ:  " << inputXYZ << "\n\n"
           << "output:    " << output << "\n\n"
           << "----------------------------------------" << std::endl;

    // ----------------------------------------------------------------------
    
    if (inputUV == "" || inputXYZ == "")
    {
        gsInfo << "No input. Exiting...\n";
        return 0;
    }
    
    gsMatrix<> uv = read_parameters(inputUV);
    gsMatrix<> xyz = read_points(inputXYZ);
    
    gsFileData<> fd;
    fd << uv;
    fd << xyz;

    gsInfo << "Output: " << output << std::endl;
    fd.dump(output);
    
    return 0;
}


    
    


