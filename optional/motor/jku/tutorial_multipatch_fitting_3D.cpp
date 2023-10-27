/** @file tutorial_multipatch_fitting_3D.cpp

    @brief Constructs a volumetric multi-patch from discrete sets of
    parameters and points sampled on the boundaries.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#include <gismo.h>
#include "gsParameterLines.h"
#include "gsFittingUtils.h"

using namespace gismo;

int main(int argc, char* argv[])
{
    std::string boundary_data = MOTOR_DATA_DIR "jku/chair_stand_boundary_data.xml";
    std::string multibasis_file = MOTOR_DATA_DIR "jku/chair_stand_multibasis.xml";
    std::string output = "parameterization";
    double lambda = 9e-14;
    
    gsCmdLine cmd("Constructs a volumetric multi-patch from discrete sets of parameters and points sampled on the boundaries");
    cmd.addString("o", "output", "Multi-patch computed from a given boundary data", output);
    cmd.addString("T", "topology", "File containing topology data", multibasis_file);
    cmd.addString("B", "boundary", "File containing boundary data", boundary_data);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input Arguments: \n"
           << "Boundary data: " << boundary_data << "\n"
           << "Multibasis:    " << multibasis_file << "\n"
           << "Output:        " << output << "\n"
           << "---------------------------------------------------------\n\n";
    
    gsFileData<> fd(multibasis_file);
    gsMultiBasis<>* basis = fd.getAnyFirst< gsMultiBasis<> >().release();
    
    std::vector< gsMatrix<> > parameters;
    std::vector< gsMatrix<> > points;
    getPoints(parameters, points, boundary_data);

    gsMultiPatch<> solution = fitting(*basis, parameters, points, lambda);
    gsInfo << "Created: " << solution << "\n\n";
    
    gsInfo << "Writing: " << output << "\n";
    gsWrite(solution, output);
    gsInfo << "solution:\n" << solution << std::endl;
    std::string parameter_lines_filename = output + "_parameter_lines";
    writeParameterLines(solution, parameter_lines_filename, 500, 7);
    
    return 0;
}


