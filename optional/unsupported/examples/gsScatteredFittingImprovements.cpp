/** @file gsScatteredFittingImprovements.cpp

    @brief Improvements of an example of scattered fitting. This is the fitting
    algorithm described in the master thesis of S. Imperatore.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Bracco, C. Giannelli, S. Imperatore, J. Speh, D. Mokris
*/

#include <gismo.h>
#include <gsHSplines/gsScatteredHFitting.h>
#include <gsHSplines/gsScatteredHFittingImprovements.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string input("");
    std::string output("");
    index_t degree               = 2;
    index_t num_u_internal_knots = 3;
    index_t num_v_internal_knots = 3;
    index_t num_iterations       = 2;
    real_t threshold             = 5e-5;
    real_t sigma                 = 1e+8;
    real_t alpha                 = 0.1;
    index_t n_loc                = 15;
    index_t strategy             = 0;
    index_t cp_magnitude         = 0;
    real_t hole                  = 0;
    real_t lambda                = 1e-06;
    index_t circle               = 2;
    bool dirU                    = false;
    bool dirV                    = false;
    

    // ----------------------------------------------------------------------

    gsCmdLine cmd("Scattered fitting using hierarchical splines.");

    cmd.addString("p", "point", "G+Smo xml file having parameters and points", input);
    cmd.addString("f","name","Name of the file", output);
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
    cmd.addInt("r", "strategy", "how to go about refinement\n"
                                 "0: check the support for n_loc points (as in the paper);\n"
                                 "1: check quadrants of the support for n_loc points;\n"
                                 "2: check support of each child for n_loc points;", strategy);
    
    cmd.addInt("m", "cp_magnitude", "minimum number of parameters for each control point", cp_magnitude);
    cmd.addReal("o", "hole", "distribuition of points taking care of eventual holes", hole);
    cmd.addInt("c", "circle", "enlargment in parameter domain made by\n"
                              "0: circles;\n"
                              "1: rectangles;\n"
                              "2: knot-span;", circle);
    cmd.addReal("l", "lambda", "stiffness weight", lambda);
    cmd.addSwitch("x", "dirU", "if true: periodic basis in u-direction", dirU);
    cmd.addSwitch("y", "dirV", "it true: periodic basis in v-direction", dirV);
    
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
           << "strategy:             " << strategy << "\n\n"
           << "cp_magnitude:         " << cp_magnitude << "\n\n"
           << "hole:                 " << hole << "\n\n"
           << "circle:               " << circle << "\n\n"
           << "lambda:               " << lambda << "\n\n"
           << "periodicity u dir     " << dirU << "\n\n"
           << "periodicity v dir     " << dirV << "\n\n"
           << "----------------------------------------" << std::endl;
    
    // ----------------------------------------------------------------------
    
    
    if (input == "")
    {
        gsInfo << "No input. Exiting...\n";
        return 0;
    }
    
    if (dirU && dirV)
    {
        gsInfo << "Periodicity in both directions not implemented yet." << std::endl;
        return 0;
    }
    
    gsFileData<> fd(input);
    gsMatrix<> parameters, points;
    fd.getId< gsMatrix<> >(0, parameters);
    fd.getId< gsMatrix<> >(1, points);
   
    gsKnotVector<> kv1(0, 1.0, num_u_internal_knots, degree + 1);
    gsKnotVector<> kv2(0, 1.0, num_v_internal_knots, degree + 1);
    gsTensorBSplineBasis<2> b(kv1, kv2);
    gsTHBSplineBasis<2, real_t> basis(b);

    gsWriteParaviewPoints(parameters, "parametrized_points");
    gsWriteParaviewPoints(points, "data_points");
    
    scatteredHFitting_improvements(basis, parameters, points, num_iterations, threshold, lambda, alpha, cp_magnitude, hole, n_loc, dirU, dirV, circle, strategy);

    //gsInfo << "The basis has " << basis.size() << " functions, out of which " << basis.numTruncated() << " are truncated." << std::endl;
        
    return 0;
}


