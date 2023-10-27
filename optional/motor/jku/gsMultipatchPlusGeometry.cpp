/** @file gsMultipatchPlusGeometry.cpp

    @brief This program adds geometries (given in separate xml files)
    to a multipatch (given in an xml file).
    (The code is based on gsGeometryToMultipatch.cpp)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius, J. Speh
*/


#include <iostream>
#include <sstream>
#include <gismo.h>
//#include "gsOffset.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

int main(int argc, char *argv[])
{

    // Options with default values
    bool save = true; // If set to true, XML file is generated
    bool plot = false; // If set to true, ParaView file is generated and launched on exit
    int npts  = 100;
    std::string file("");
    std::string multipatchFile(MOTOR_DATA_DIR "jku/target_segmentation.xml");
    // geoFiles is a string of several file names
    std::string geometryFiles(MOTOR_DATA_DIR "jku/target_rectangle.xml " MOTOR_DATA_DIR "jku/target_input_offset.xml");
    std::string outputFile("target_segmentation_and_geometry");
  
    gsCmdLine cmd("This program constructs a multipatch geometry (creates an xml file)\n from separate multipatch and geometry files.");
    //cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addString("m", "multipatch", "Name of the input (multipatch) file", multipatchFile);
    cmd.addString("f", "geometry", "Names of the input (geometries) files", geometryFiles);
    cmd.addInt("s", "npts", "Number of points used for sampling each patch (Paraview visualization)", npts);
    cmd.addString("o", "output", "The name of output file (prefix for all output files)", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Multipatch file:                " << multipatchFile << "\n"
           << "Geometry files:                 " << geometryFiles << "\n"
           << "s (for ParaView visualization): " << npts << "\n"
           << "plot:                           " << plot << "\n"
           << "Output file:                    " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    std::istringstream iss(geometryFiles);

    gsGeometry<>::uPtr geom;
    gsMultiPatch<> * mp = readMultipatch(multipatchFile);

    while(iss >> file)
    {
        geom = gsGeometry<>::uPtr(readGeometry(file));
        mp->addPatch(give(geom));
    }

    // Output an xml file
    if(save)
    {
        writeMultipatch(*mp, outputFile);
    }
    
    if(plot)
    {
        // Output a paraview file
        gsInfo << "Writing " << outputFile << ".pvd" << "\n";
        gsWriteParaview(*mp , outputFile, npts);
        //Call paraview on exit
        if (mp->nPatches() == 1)
        {
            std::string command = "paraview " + outputFile +".vts &";
            return system(command.c_str());
        }
        else
        {
            std::string command = "paraview " + outputFile +".pvd &";
            return system(command.c_str());
        }
    }
    else
    {
        return 0;
    }
}

