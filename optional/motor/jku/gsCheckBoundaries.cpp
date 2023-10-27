/** @file gsCheckBoundaries

    @brief Plots all the boundaries of multipatch with appropriate names. 

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

// This program is usefull for setting boundary conditions to various
// simulation problems.
//
// The program prints all the boundaries with appropriate names. The names help
// a programmer to determine how to set boundary conditions.



#include <iostream>
#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{ 
    std::string multipatchFile(MOTOR_DATA_DIR "jku/airPassageParameterization.xml");
    std::string output("boundaries");
    
    gsCmdLine cmd("Check the boundaries");
    cmd.addString("m", "m", "Multi patch file", multipatchFile);
    cmd.addString("o", "out", "Output file", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "---------------------------------------------------------\n\n"
           << "Input Arguments: \n\n"
           << "Multipatch file: " << multipatchFile << "\n\n"
           << "Output: " << output << "\n\n"
           << "---------------------------------------------------------\n"
           << std::endl;


    gsFileData<> fd(multipatchFile);
    gsMultiPatch<>* mp = fd.getAnyFirst< gsMultiPatch<> >().release();

    
    for (gsMultiPatch<>::const_biterator it = mp->bBegin();
         it != mp->bEnd(); ++it)
    {
        const gsGeometry<> & g = mp->patch( it->patch );
        gsGeometry<>::uPtr b = g.boundary(it->side());
        std::stringstream ss;
        ss << it->side();
        std::string s;
        ss >> s;
        
        std::string out = output + "_Patch_" + util::to_string(it->patch)
            + "_Side_" + s;

        std::cout << "Boundary: patch = " << it->patch 
                  << " side = " << it->side() << std::endl;

        std::cout << "  writing: " << out << std::endl;
        gsWriteParaview(*b, out);
    }

    return 0;
}
