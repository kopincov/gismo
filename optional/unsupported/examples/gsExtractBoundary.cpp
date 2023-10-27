/** @file gsExtractBoundary.cpp

    @brief Provides implementation of input/output XML utilities struct.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Moore
*/

#include <iostream>

#include <gismo.h>



using namespace gismo;


int main(int argc, char *argv[])
{
    bool plot = false;
    std::string fn("volumes/cube.xml");
    real_t tol = 1e-4;

    gsCmdLine cmd("Extracting the boundary of B-spline geometries.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addString("g","geometry","File containing Geometry (.axl, .txt)", fn);
    cmd.addReal  ("t","tolerance","Tolerance for identifing patch interfaces", tol);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr input = gsReadFile<>(fn);
    gsInfo << "Loaded: "<< *input << "\n";

    // Compute the interfaces and boundaries of the input (in case it
    // is not already loaded from the file)
    input->computeTopology(tol);
    gsInfo << "Identified topology: "<< *input << "\n";

    gsMultiPatch<> result;

    for (gsMultiPatch<>::const_biterator it = input->bBegin();
         it != input->bEnd(); ++it)
    {
        const gsGeometry<> & g = input->patch( it->patch );
        result.addPatch( g.boundary( it->side() ) );
    }

    /*
    // Equivalent code, but inefficient
    const int d = input->parDim();
    for (index_t i = 0; i< input->nPatches(); ++i)
    {
        const gsGeometry<> & g = input->patch(i);

        for (int k = 1; k <= 2*d; ++k) // for all boundaries
        {
            if ( input->isBoundary(i,k) )
                result.addPatch( g.boundary(k) );
        }
    }
    //*/

    // Get the topology of the boundaries
    result.computeTopology(tol);

    gsWrite(result, "extracted_boundaries.xml");

    gsInfo << "Extracted "<< result.nPatches() << " boundaries in extracted_boundaries.xml.\n";

    int exitCommand(0);
    if (plot) 
    {
        gsWrite(result, "boundary_out");

        gsInfo<<"Plotting result in Paraview...\n";
        gsWriteParaview( result , "boundary_paraview", 500);

        // Run paraview on exit
        exitCommand = system("paraview boundary_paraview.pvd &");
    }

    return exitCommand;
}
