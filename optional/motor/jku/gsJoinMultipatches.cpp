/** @file gsJoinMultipatches.cpp

    @brief This program joins several multipatch files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
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
    // multipatchFiles is a string of several file names
    std::string multipatchFiles(MOTOR_DATA_DIR "jku/target_multipatch.xml " MOTOR_DATA_DIR  "jku/target_segmentation.xml " MOTOR_DATA_DIR "jku/target_extended_segmentation.xml");
    std::string outputFile("target_segmentation_and_geometry");
  
    gsCmdLine cmd("This program constructs a multipatch geometry (creates an xml file)\n from separate multipatch and geometry files.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addString("f", "multipatches", "Names of the input files", multipatchFiles);
    cmd.addInt("s", "npts", "Number of points used for sampling each patch (Paraview visualization)", npts);
    cmd.addString("o", "output", "The name of output file (prefix for all output files)", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input files:                    " << multipatchFiles << "\n"
           << "s (for ParaView visualization): " << npts << "\n"
           << "plot:                           " << plot << "\n"
           << "Output file:                    " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    std::istringstream iss(multipatchFiles);

    gsMultiPatch<>* mp     = new gsMultiPatch<>();
    gsMultiPatch<>* mp_new = NULL;

    while(iss >> file)
    {
        mp_new = readMultipatch(file);
        for (std::size_t i = 0; i != mp_new->nPatches(); i++)
        {
            const gsGeometry<>& geom = (*mp_new)[i];
            mp->addPatch(geom);
        }
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

