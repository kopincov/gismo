/** @file tutorial_parameterization_fitting.cpp

    @brief Constructs a patch (volume) from discrete sets of
    parameters and points sampled on the boundaries (surfaces)
    (see gsModelling/gsFitting).

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

// this file is very similar to the tutorial_parameterization_fitting.cpp


#include <gismo.h>
#include "gsParameterLines.h"

using namespace gismo;

int main(int argc, char* argv[])
{
    std::string boundary_data = MOTOR_DATA_DIR "jku/skrew_boundary_data.xml";
    std::string output = "parameterization";
    index_t num_internal_knots = 5;
    index_t degree = 2;
    double lambda = 9e-14;

    gsCmdLine cmd("Constructs a patch (volume) from discrete sets of parameters and points sampled on the boundaries (surfaces) (see gsModelling/gsFitting)");
    cmd.addString("o", "output", "Patch computed from a given boundary data", output);
    cmd.addInt("k", "numKnots", "Number of interior knots in one parameter direction", num_internal_knots);
    cmd.addInt("d", "degree", "Degree", degree);
    cmd.addString("B", "boundary", "File containing boundary data", boundary_data);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    gsInfo << "Input Arguments: \n"
           << "Boundary data: " << boundary_data << "\n"
           << "Output:        " << output << "\n"
           << "numKnots:      " << num_internal_knots << "\n"
           << "degree:        " << degree << "\n"
           << "---------------------------------------------------------\n\n";
    
    
    gsFileData<> fd(boundary_data);
    gsMatrix<> parameters = *fd.getId< gsMatrix<> >(0);
    gsMatrix<> points = *fd.getId< gsMatrix<> >(1);
    
    gsKnotVector<> kv1(0.0, 1.0, num_internal_knots, degree);
    gsKnotVector<> kv2(0.0, 1.0, num_internal_knots, degree);
    gsKnotVector<> kv3(0.0, 1.0, num_internal_knots, degree);
    gsTensorBSplineBasis<3> basis(kv1, kv2, kv3);
    
    gsFitting<> fitting(parameters, points, basis);
    fitting.compute(lambda);
    const gsGeometry<>* result = fitting.result();
    gsInfo << "Created: " << *result << "\n\n";

    gsInfo << "Saving: " << output << "\n";
    gsWrite(*result, output);
    std::string parameter_lines_filename = output + "_parameter_lines";
    writeParameterLines(*result, parameter_lines_filename, 500, 10);


    return 0;
}


/*
fitting with uniform refinement
for fitting with adaptive refinement look at stable/examples/fitting_example.cpp

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





