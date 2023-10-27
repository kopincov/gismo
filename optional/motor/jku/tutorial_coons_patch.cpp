/** @file tutorial_coons_patch.cpp

    @brief Constructs a patch (surface, volume, ...) from a given set of
    domain boundaries using a Coon's parameterization.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#include <gismo.h>
#include <gsModeling/gsCoonsPatch.h>
#include "gsParameterLines.h"

using namespace gismo;

int main(int argc, char* argv[])
{
    // Try the following boundaries.
    // MOTOR_DATA_DIR "jku/wavy_boundary_curves.xml"
    // MOTOR_DATA_DIR "jku/skrew_boundary_surfaces.xml

    std::string boundary_filename = MOTOR_DATA_DIR "jku/wavy_boundary_curves.xml";
    std::string output = "coons_patch";
    
    gsCmdLine cmd("Constructs a patch (surface, volume, ...) from a given set of domain boundaries using a Coon's parameterization");
    cmd.addString("o", "output", "Patch computed from a given boundary data", output);
    cmd.addString("B", "boundary", "File containing boundary data", boundary_filename);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    gsInfo << "Input Arguments: \n"
           << "Boundary: " << boundary_filename << "\n"
           << "Output:   " << output << "\n"
           << "---------------------------------------------------------\n\n";
    
    
    // Boundary must contain four curves (6 surfaces)
    // All of the boundaries must have the same degree!

    gsMultiPatch<> boundary;
    gsReadFile<>(boundary_filename, boundary);

    gsCoonsPatch<real_t> coons_patch(boundary);
    const gsGeometry<>& patch = coons_patch.compute();
    gsInfo << "Created " << patch << "\n\n";
    

    gsInfo << "Saving: " << output << "\n";
    gsWrite(patch, output);
    std::string parameter_lines_filename = output + "_parameter_lines";
    writeParameterLines(patch, parameter_lines_filename, 500, 10);
    
    return 0;
}

    
