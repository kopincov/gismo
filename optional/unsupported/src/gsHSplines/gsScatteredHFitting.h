/** @file gsScatteredHFitting.h

    @brief Scattered fitting using hierarchical splines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Bracco, C. Giannelli, S. Imperatore, D. Mokris, J. Speh
*/

#pragma once

#include <gsIO/gsIOUtils.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsScatteredHFittingRefinement.h>
#include <gsHSplines/gsScatteredHFittingCoefficients.h>
#include <gsHSplines/gsScatteredHFittingLSPlotting.h>
#include <gsFitting/gsVTK.h>
#include <gsCore/gsGeometry.h>

namespace gismo {

void scatteredHFitting(gsTHBSplineBasis<2>& basis,
                       const gsMatrix<>& parameters,
                       const gsMatrix<>& points,
                       const int num_iterations,
                       const real_t threshold,
                       const real_t sigma,
                       const real_t alpha,
                       const int n_loc,
                       bool circle,
                       int strategy,
                       int plot_cp = -1)
{
    for (int iter = 0; iter != num_iterations; iter++)
    {
        gsInfo << "\nIteration: " << iter + 1 << " / " << num_iterations << std::endl;

        gsMatrix<> coefficients;
        if(plot_cp >= 0)
            coefficients = get_coefficients_print_LS_surface(basis, parameters, points, sigma, iter, 1, plot_cp, num_iterations, circle);
        else
            coefficients = get_coefficients(basis, parameters, points, sigma, circle);

        gsTHBSpline<2> QI(basis, coefficients); // quasi interpolant

        if(plot_cp >= 0 && iter == num_iterations - 1)
        {
            gsMatrix<real_t> transpose = coefficients.transpose();
            std::vector<real_t> coef_num;
            for(int i = 0; i < basis.size(); i++)
                coef_num.push_back(i);

            vtk_converter(transpose, coef_num, "coefficients.vtk");
        }
        
        gsMatrix<> errors = compute_errors(QI, parameters, points);
        real_t percentage = percentage_of_points_below_threshold(errors, threshold);

        std::cout << "error percentage below threshold: " << 100 * percentage << "%\n";

        // saving mesh
        gsMesh<> mesh;
        makeMesh<>(basis, mesh, 4);
        std::string output = "iteration";
        std::string prefix = output + internal::to_string(iter);

        gsWriteParaview(mesh, prefix + "mesh");
        QI.evaluateMesh(mesh);
        gsWriteParaview(mesh, prefix + "geometry");
        gsWriteParaview(QI, prefix + "thb", 50000, false, true);
        
        if ((iter != num_iterations - 1) && (alpha < 1 - percentage))
        {
            scatteredHFittingRefine(basis, parameters, errors, threshold, n_loc, strategy);
        }
        else
        {
            break;
        }
    }
}

} //namespace gismo




