/** @file gsGeometryToMultipatch.cpp

    @brief Constructs a multipatch geometry (creates an xml file)
    from several single geometries stored in separate xml files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius, J. Speh
*/
// ./bin/gsGeometryToMultipatch -f "2d_target_geometry_curve_1.xml 2d_target_geometry_curve_2.xml" -s 100 -o 2d_target_geometry --plot

#include <iostream>
#include <sstream>
#include <gismo.h>
#include "gsMotorIOUtils.h"

using namespace gismo;

int main(int argc, char *argv[])
{

    // Options with default values
    bool save = true;  // If set to true, XML file is generated
    bool plot = false; // If set to true, ParaView file is generated and launched on exit
    int npts  = 100;
    std::string file("");
    // input is a string of several file names
    std::string input(MOTOR_DATA_DIR "jku/target_rectangle.xml " MOTOR_DATA_DIR "jku/target_input.xml");
    std::string outputFile("target_multipatch");

    gsCmdLine cmd("This program constructs a multipatch geometry (creates an xml file)\n from several single geometries stored in separate files.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addString("f", "input", "Names of the input files", input);
    cmd.addInt("s", "npts", "Number of points used for sampling each patch (ParaView visualization)", npts);
    cmd.addString("o", "output", "The name of output file", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input files:                    " << input << "\n"
           << "s (for ParaView visualization): " << npts << "\n"
           << "plot:                           " << plot << "\n"
           << "Output file:                    " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    std::istringstream iss(input);

    gsGeometry<>::uPtr geom;
    gsMultiPatch<> mp;

    while(iss >> file)
    {
        geom = gsGeometry<>::uPtr(readGeometry(file));
        mp.addPatch(give(geom));
    }

    // Output an xml file
    if(save)
    {
        writeMultipatch(mp, outputFile);
    }
    
    if(plot)
    {
        // Output a paraview file
        gsInfo << "Writing " << outputFile << ".pvd" << "\n";
        gsWriteParaview(mp , outputFile, npts);
        //Call paraview on exit
        if (mp.nPatches() == 1)
        {
            gsFileManager::open("paraview " + outputFile +".vts");
        }
        else
        {
            gsFileManager::open("paraview " + outputFile +".pvd");
        }
    }
    return 0;
}

