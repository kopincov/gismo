/** @file gsNonMatching_test.cpp

    @brief Testing patches with non-matching parametrization at the interface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Seiler
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{   
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 0;
    index_t numNHref    = 0;
    index_t basisDegree = 0;
    bool plot       = false;

    // Multipatch object


    int result = 0;
    std::string fn( "planar/two_squares.xml");
    //std::string fn_pde("");
    
    gsCmdLine cmd("Testing a multipatch problem.");
    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("R","NonMatch", "Number of non-matching elements steps to perform before solving",
               numNHref);
    cmd.addInt("p","degree", 
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)", fn);
    //cmd.addString("q","pde","File containing a poisson PDE (.xml)", fn_pde);
    cmd.addInt("e","degreeElevation", 
               "Number of degree elevation steps to perform on the Geometry's basis before solving", 
               numElevate);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
        
    // Read the input
    gsFileData<> fileData(fn);

    if ( !fileData.has< gsMultiPatch<> >() )
    {
        gsWarn<<"Reading "<<fn<<" failed.\n";
        return -1;
    }

    // We have the input
    gsMultiPatch<> mp;
    fileData.getFirst< gsMultiPatch<> >(mp);
    
    gsInfo <<"Read a "<< mp <<"\n";

    mp.computeTopology();

    gsInfo <<"Now: "<< mp <<"\n";
    

    return result;
}
