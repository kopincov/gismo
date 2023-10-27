/** @file gsScatteredHFittingImprovements.h

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

#include <gsHSplines/gsScatteredHFittingCoefficientsOpt.h>
#include <gsHSplines/gsScatteredHFittingCoefficientsPeriodic.h>
#include <gsHSplines/gsScatteredHFittingRefinement.h>
#include <gsHSplines/gsScatteredHFittingRefinementExtended.h>

#include <gsFitting/gsVTK.h>
// Uncomment for plotting in Parasolid.
//#include <gsParasolid/gsWriteParasolid.h>

namespace gismo{
/** \brief returns true if a periodic fitting in \a u or \a v direction is required.
 */
bool periodic_surface(const bool u, const bool v)
{
    if(!( u || v ))
        return false;
    return true;
}

/** \brief returns false if the surface is \a u periodic, true if the surface is \a v periodic. 
 */
bool get_periodic_dir(const bool u, const bool v)
{
    if(u)
        return false;
    return true;
}

/** \brief Hierarchical fitting of scattered data points
 *  \param[in] parameters parametrized data
 *  \param[in] points data points in the spatial domain
 *  \param[in] num_iteration maximum level of refinement
 *  \param[in] threshold maximum threshold
 *  \param[in] lambda weight for the smoothing term
 *  \param[in] alpha 1 - \em alpha percentage of points below \em threshold
 *  \param[in] cp_magnitude minimum number of points for the local fitting surface
 *  \param[in] hole minimum density of points for the local fitting inside the local domain
 *  \param[in] n_loc minimum number of points for refinement
 *  \param[in] dirU if true periodic surface in \a u direction
 *  \param[in] dirV if true periodic surface in \a v direction 
 *  \param[in] circle shape of enlargement
 *  \param[in] strategy refinement strategy
 */
void scatteredHFitting_improvements(gsTHBSplineBasis<2>& basis,
                                    const gsMatrix<>& parameters,
                                    const gsMatrix<>& points,
                                    const int num_iterations,
                                    const real_t threshold,
                                    const real_t lambda,
                                    const real_t alpha,
                                    const int cp_magnitude,
                                    const real_t hole,
                                    const int n_loc,
                                    const bool dirU,
                                    const bool dirV,
                                    const int circle,
                                    const int strategy)
{
    gsMatrix<> coefficients;
    gsMatrix<> old_coefficients;
    std::vector<std::map<index_t, index_t> > tensorIndices;  // tensorIndices[lvl] assigns the Hierarchical index to the tensor-product index of level lvl.
    int dir = -1;
    
    std::vector<index_t> enlargement;
    std::vector<index_t> old_enlargement;
 
    const bool is_periodic = periodic_surface(dirU, dirV);
    if(is_periodic)
        dir = get_periodic_dir(dirU, dirV);
    
    for (int iter = 0; iter != num_iterations; iter++)
    {
        gsInfo << "\nIteration: " << iter + 1 << " / " << num_iterations << std::endl;
        
        if(iter == 0)
        {
            gsStopwatch time;
            time.restart();
            if(is_periodic)
                coefficients = get_coefficients_periodic(basis, parameters, points, iter, cp_magnitude, hole, lambda, dir, circle, enlargement);
            else
                coefficients = get_coefficients_not_scaled_stiff(basis, parameters, points, iter, cp_magnitude, hole, lambda, circle, enlargement);
            time.stop();
            gsInfo << "Fitting time: " << time << std::endl;
        }
        else
        {
            old_coefficients = coefficients;
            old_enlargement = enlargement;
            coefficients.resize(basis.size(), points.rows());
            gsStopwatch time;
            time.restart();
            if(is_periodic)
                coefficients = get_coefficients_periodic_opt(old_coefficients, tensorIndices, basis, parameters, points, iter, cp_magnitude, hole, lambda, dir, circle, old_enlargement, enlargement);
            else
                coefficients = get_coefficients_opt(old_coefficients, tensorIndices, basis, parameters, points, iter, cp_magnitude, hole, lambda, circle, old_enlargement, enlargement);
            time.stop();
             gsInfo << "Fitting time: " << time << std::endl;
        }
         
        
        gsTHBSpline<2> QI(basis, coefficients); // quasi interpolant
        
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
        gsWriteParaview(QI, prefix + "thb", 10000, false, true);
        
        
        // tensorBspline
        gsTensorBSpline<2, real_t> tensorB;
        QI.convertToBSpline(tensorB);
        gsTensorBSplineBasis<2, real_t> tensBasis = tensorB.basis();
        gsMesh<> TensorMesh;
        makeMesh<>(tensBasis, TensorMesh, 4);
        gsWriteParaview(TensorMesh, prefix + "tensor_mesh");
        QI.evaluateMesh(TensorMesh);
        gsWriteParaview(TensorMesh, prefix + "tensor_geometry");
        
        std::cout << "Tensor DOFs: " << tensBasis.size() << std::endl;
        std::cout << "THB DOFs: " << basis.size() << std::endl;
        
         // Advanced visualisation with NX:
        /*gsTHBSpline<2> copy(QI);
        gsTensorBSpline<2> TP_spline;
        copy.convertToBSpline(TP_spline);
        extensions::gsWritePK_SHEET(TP_spline, "iteration" + std::to_string(iter));*/

        if ((iter != num_iterations - 1) && (alpha < 1 - percentage))
        {
            // Fill tensorIndices before refinement.
            tensorIndices.clear();
            tensorIndices.resize(basis.getBases().size()+1); // +1, because there will be an extra level in the next step.
            for(index_t k = 0; k != basis.size(); k++)
            {
                int lvl = basis.levelOf(k);
                tensorIndices[lvl].insert(std::pair<index_t, index_t>(basis.flatTensorIndexOf(k), k));
            }
            if(is_periodic)
                scatteredHFittingRefine_periodic(basis, parameters, errors, threshold, n_loc, dir, strategy, circle, enlargement);
            else
                scatteredHFittingRefine(basis, parameters, errors, threshold, n_loc, strategy, circle, enlargement);
            
            if(old_coefficients.rows() == coefficients.rows())
            {
                std::cout << "No further refinement: alpha stopping criteria not satisfied. " << std::endl;
                break;
            }
        }
        else
        {
            break;
        }
    }
}
}// namespace gismo
