/** @file tutorial_multipatch_fitting_2.cpp

    @brief Constructs a planar multi-patch from discrete sets of
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
    std::string boundary_data = MOTOR_DATA_DIR "jku/triangle_boundary_data.xml";
    std::string output = "parameterization_2";
    int num_internal_knots = 9;
    int degree = 2;
    double lambda = 9e-14;

    gsCmdLine cmd("Constructs a planar multi-patch from discrete sets of parameters and points samples on the boundaries");
    cmd.addString("o", "output", "Multi-patch computed from a given boundary data", output);
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
    

    std::vector< gsMatrix<> > parameters;
    std::vector< gsMatrix<> > points;
    getPoints(parameters, points, boundary_data);
    
    
    gsMultiPatch<> multi_patch;
    for (std::size_t patch = 0; patch != parameters.size(); ++patch)
    {
	gsKnotVector<> kv1(0.0, 1.0, num_internal_knots, degree);
	gsKnotVector<> kv2(0.0, 1.0, num_internal_knots, degree);
	gsTensorBSplineBasis<2> basis(kv1, kv2);
	
	gsFitting<> fitting(parameters[patch], points[patch], basis);
	fitting.compute(lambda);
	gsGeometry<>* result = fitting.result();
	multi_patch.addPatch(result->clone());
	
	gsInfo << "Fitting " << patch + 1 << " complete.\n";
    }
    
    gsInfo << "\nCreated: " << multi_patch << "\n";
    std::string out = output + "_no_topology";
    gsInfo << "Saving " << out << std::endl;
    gsWrite(multi_patch, out);
    writeParameterLines(multi_patch, out, 500, 10);
    
    gsInfo << "Computing topology...\n";
    multi_patch.computeTopology();
    multi_patch.closeGaps();


    gsInfo << "\nCreated: " << multi_patch << "\n";
    out = output + "_with_topology";
    gsInfo << "Saving " << out << "\n";
    gsWrite(multi_patch, out);
    writeParameterLines(multi_patch, out, 500, 10);
    
    return 0;
}










