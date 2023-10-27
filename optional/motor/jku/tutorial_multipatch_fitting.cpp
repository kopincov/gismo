/** @file tutorial_multipatch_fitting.cpp

    @brief Constructs a planar multi-patch from discrete sets of
    parameters and points sampled on the boundaries, and topology data.

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
    std::string boundary_data = MOTOR_DATA_DIR "jku/triangle_boundary_data.xml";
    std::string topology_data = MOTOR_DATA_DIR "jku/triangle_topology.xml";
    std::string output = "parameterization";
    int num_internal_knots = 5;
    int degree = 2;
    double lambda = 9e-14;
    
    gsCmdLine cmd("Constructs a planar multi-patch from discrete sets of parameters and points sampled on the boundaries, and topology data");
    cmd.addString("o", "output", "Multi-patch computed from a given boundary data", output);
    cmd.addInt("k", "numKnots", "Number of interior knots in one parameter direction", num_internal_knots);
    cmd.addInt("d", "degree", "Degree", degree);
    cmd.addString("T", "topology", "File containing topology data", topology_data);
    cmd.addString("B", "boundary", "File containing boundary data", boundary_data);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input Arguments: \n"
           << "Boundary data: " << boundary_data << "\n"
           << "Topology data: " << topology_data << "\n"
           << "Output:        " << output << "\n"
           << "numKnots:      " << num_internal_knots << "\n"
           << "degree:        " << degree << "\n"
	   << "---------------------------------------------------------\n\n";
    
    gsFileData<> fd(topology_data);
    gsMultiPatch<>* topology = fd.getAnyFirst< gsMultiPatch<> >().release();
    
    std::vector< gsMatrix<> > parameters;
    std::vector< gsMatrix<> > points;
    getPoints(parameters, points, boundary_data);

    gsMultiBasis<> basis = makeMultiBasisBSplines(*topology, degree, num_internal_knots);
    
    gsMultiPatch<> solution = fitting(basis, parameters, points, lambda);
    gsInfo << "Created: " << solution << "\n\n";

    gsInfo << "Writing: " << output << "\n";    
    gsWrite(solution, output);  
    std::string parameter_lines_filename = output + "_parameter_lines";
    writeParameterLines(solution, parameter_lines_filename, 500, 10);
    
    return 0;
}

/*
fitting with uniform refinement

const int num_iterations = 6;
const real_t threshold = 1e-4;
for (int i = 0; i != num_iterations; ++i)
{
    fitting.compute(lambda);
    fitting.computeErrors();
    real_t max_error = fitting.maxPointError();
    if (max_error < threshold)
    {
        break;
    }
    else
    {
        basis.uniformRefine();
    }
}
*/


