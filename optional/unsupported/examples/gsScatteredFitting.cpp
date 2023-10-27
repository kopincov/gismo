/** @file gsScatteredFitting.cpp

    @brief An example of scattered fitting. This is the fitting algorithm
    described in
    C. Bracco et al.: Adaptive fitting with THB-splines: Error analysis and industrial applications,
    CAGD 2018.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Bracco, C. Giannelli, J. Speh, D. Mokris
*/

#include <iostream>
#include <fstream>
#include <vector>

#include <gismo.h>
#include <gsHSplines/gsScatteredHFitting.h>
#include <gsHSplines/gsScatteredHFittingLSPlotting.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    std::string input("");
    index_t degree               = 2;
    index_t num_u_internal_knots = 3;
    index_t num_v_internal_knots = 3;
    index_t num_iterations       = 2;
    real_t threshold             = 5e-5;
    real_t sigma                 = 1e+8;
    real_t alpha                 = 0.1;
    index_t n_loc                = 15;
    bool circle                  = true;
    index_t plot_cp              = -1;
    index_t strategy             = 0;
    

    // ----------------------------------------------------------------------

    gsCmdLine cmd("Scattered fitting using hierarchical splines.");

    cmd.addString("p", "point", "G+Smo xml file having parameters and points", input);
    cmd.addInt("d", "deg", "Degree of initial basis.", degree);
    cmd.addInt("u", "numInternalKnotsU", "Number of internal knots of the initial basis in the u-direction",
               num_u_internal_knots);
    cmd.addInt("v", "numInternalKnotsV", "Number of internal knots of the initial basis in the v-direction",
               num_v_internal_knots);
    cmd.addInt("i", "iterations", "Maximum number of iterations", num_iterations);
    cmd.addReal("t", "threshold", "Maximum threshold", threshold);
    cmd.addReal("s", "sigma", "Look at the papers for the definition", sigma);
    cmd.addReal("a", "alpha", "Look at the papers for the definition", alpha);
    cmd.addInt("n", "n_loc", "Look at the papers for the definition", n_loc);
    cmd.addSwitch("c", "circle", "if true, enlarge with circles in parameter domain; else with rectangles", circle);
    cmd.addInt("k", "plot_cp", "If set to a non-negative integer, the local fitting surface corresponding to that CP is plotted in paraview. Default: false.", plot_cp);
    cmd.addInt("r", "strategy", "how to go about refinement; 0: check the support for n_loc points (as in the paper); 1: check quadrants of the support for n_loc points, 2: check support of each child for n_loc points", strategy);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // ----------------------------------------------------------------------
    
    gsInfo << "----------------------------------------\n"
           << "INPUT DATA:\n\n" 
           << "input:                " << input << "\n\n"
           << "degree:               " << degree << "\n\n"
           << "num internal knots u: " << num_u_internal_knots << "\n\n"
           << "num internal knots v: " << num_v_internal_knots << "\n\n"
           << "max. num iterations:  " << num_iterations << "\n\n"
           << "max. threshold:       " << threshold << "\n\n"
           << "sigma:                " << sigma << "\n\n"
           << "alpha:                " << alpha << "\n\n"
           << "n_loc:                " << n_loc << "\n\n"
	       << "plot_cp:              " << plot_cp << "\n\n"
           << "circle:               " << circle << "\n\n"
           << "----------------------------------------" << std::endl;
    
    // ----------------------------------------------------------------------
        
    if (input == "")
    {
        gsInfo << "No input. Exiting...\n";
        return 0;
    }
    
    
    gsFileData<> fd(input); // The data contained in this file must be parametrized between 0 and 1, otherwise the construction with the gsKnotVector does not work. 
			    // If they are not properly scaled, we can scale them with the function scalingPoints in "examples" > gsParToXML.cpp
    gsMatrix<> parameters, points;
    fd.getId< gsMatrix<> >(0, parameters);
    fd.getId< gsMatrix<> >(1, points);
    
    gsKnotVector<> kv1(0.0, 1.0, num_u_internal_knots, degree + 1);
    gsKnotVector<> kv2(0.0, 1.0, num_v_internal_knots, degree + 1);
    gsTensorBSplineBasis<2> b(kv1, kv2);
    gsTHBSplineBasis<2> basis(b);


    gsWriteParaviewPoints(parameters, "parameters");
    gsWriteParaviewPoints(points, "points");

    scatteredHFitting(basis, parameters, points, num_iterations, threshold, sigma, alpha, n_loc, circle, strategy, plot_cp);

    return 0;
}


