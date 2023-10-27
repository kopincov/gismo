/** @file gsVTK.cpp

    @brief It allows to convert a .par file into a .vtk file.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <gsFitting/gsVTK.h>
#include <gsFitting/gsPointCloudTools.h>

#include <gismo.h>
using namespace gismo;


int main(int argc, char *argv[])
{
    std::string input("");
    std::string output("");
    
// ----------------------------------------------------------------------

    gsCmdLine cmd("File in .vtk format");
    cmd.addString("i", "input", ".par file with parameters and points", input);
    cmd.addString("o", "output", ".vtk file", output);
     
    
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    gsInfo << "----------------------------------------\n"
           << "INPUT DATA:\n\n" 
           << "input:                " << input << "\n\n"
           << "output:               " << output << "\n\n"
           << "----------------------------------------" << std::endl;
    
            
    if (input == "")
    {
        gsInfo << "No input. Exiting...\n";
        return 0;
    }
    
    
    gsMatrix<real_t> parameters;
    gsMatrix<real_t> points;
    fillTheMatrices(input, parameters, points);
    
    std::vector<real_t> info;
    for(index_t i = 0; i < parameters.cols(); i++)
        info.push_back(parameters(0,i));
    
    vtk_converter(points, info, output);
    
return 0;
}
