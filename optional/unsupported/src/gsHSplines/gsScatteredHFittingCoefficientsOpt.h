/** @file gsScatteredHFittingCoefficientsOpt.h

    @brief Computation of coefficients with Stiffness, taking care of the previous iteration.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#pragma once
#include <gismo.h> // enough if only the stable (public) part of gismo is used.
#include <algorithm>
#include <cmath>

//#include <gsHSplines/gsTHBSplineBasis.h>
//#include <gsHSplines/gsTHBSpline.h>

//#include <gsHSplines/gsScatteredHFittingImprovementsGeneralUtilities.h>
//#include <gsHSplines/gsScatteredHFittingImprovementsCoefficients.h>
//#include <gsHSplines/gsScatteredHFittingImprovementsStiffness.h>
#include <gsHSplines/gsScatteredHFittingCoefficientsStiffness.h>
#include <gsHSplines/gsScatteredHFittingCoefficientsPeriodic.h>



namespace gismo{

/** \brief Computation of coefficients for gsTHBSplineBasis \em basis.
 *         coefficients that do not change with respect to tha basis functions after refinement, are not recomputed.
 */    
gsMatrix<> get_coefficients_opt(const gsMatrix<>& old_coefficients,
                                const std::vector<std::map<index_t, index_t> >& tensorIndices,
                                const gsTHBSplineBasis<2>& basis,
                                const gsMatrix<>& parameters,
                                const gsMatrix<>& points,
                                int iteration,
                                const int cp_magnitude,
                                const real_t hole,
                                const real_t lambda,
                                const int circle,
                                std::vector<index_t>& old_enlargement,
                                std::vector<index_t>& enlargement)
{
    gsMatrix<> coefficients(basis.size(), points.rows());
    
    int recomputed = 0;
     
    for (index_t h = 0; h != basis.size(); h++) // h is the hierarchical index, t is the local tensor-product index
    {
        gsMatrix<> coefficient;
        
        const std::map<index_t, index_t>& tensorMapOfThisLvl = tensorIndices[basis.levelOf(h)];
        int t = basis.flatTensorIndexOf(h);
        
        std::map<index_t, index_t>::const_iterator found = tensorMapOfThisLvl.find(t);
        if(found==tensorMapOfThisLvl.end())
        {
            coefficient = get_coefficient_not_scaled_stiff(basis, parameters, points, iteration, h, cp_magnitude, hole, lambda, circle, enlargement);
            recomputed++;
        }
        else
        {
            coefficient = old_coefficients.row(found->second);
            enlargement.push_back(old_enlargement[found->second]);
        }
                
        coefficients.row(h) = coefficient.row(0);
        /*gsMatrix<> printCoefficient = coefficient.transpose();
        gsWriteParaviewPoints(printCoefficient, "iteration" + std::to_string(iteration) + "_control point" + std::to_string(h));*/
    }
    
    //gsInfo << "Nr. of recomputed coefficients: " << recomputed << std::endl;
    return coefficients;
}

/** \brief Computation of coefficients for gsTHBSplineBasis \em basis for periodic surfacesS.
 *         coefficients that do not change with respect to tha basis functions after refinement, are not recomputed.
 */
gsMatrix<> get_coefficients_periodic_opt(const gsMatrix<>& old_coefficients,
                                         const std::vector<std::map<index_t, index_t> >& tensorIndices,
                                         const gsTHBSplineBasis<2>& basis,
                                         const gsMatrix<>& parameters,
                                         const gsMatrix<>& points,
                                         int iteration,
                                         const int cp_magnitude,
                                         const real_t hole,
                                         const real_t lambda,
                                         const int dir,
                                         const int circle,
                                         std::vector<index_t>& old_enlargement,
                                         std::vector<index_t>& enlargement)
{
    gsMatrix<> coefficients(basis.size(), points.rows());
    int recomputed = 0;
     
    for (index_t h = 0; h != basis.size(); h++) // h is the hierarchical index, t is the local tensor-product index
    {
        gsMatrix<> coefficient;
        
        const std::map<index_t, index_t>& tensorMapOfThisLvl = tensorIndices[basis.levelOf(h)];
        int t = basis.flatTensorIndexOf(h);
        
        std::map<index_t, index_t>::const_iterator found = tensorMapOfThisLvl.find(t);
        if(found==tensorMapOfThisLvl.end())
        {
            coefficient = get_coefficient_not_scaled_stiff(basis, parameters, points, iteration, h, cp_magnitude, hole, lambda, circle, enlargement);
            recomputed++;
        }
        else
        {
            coefficient = old_coefficients.row(found->second);
            enlargement.push_back(old_enlargement[found->second]);
        }
        
        coefficients.row(h) = coefficient.row(0);
        
    }
    
     // composite basis
    std::vector< gsHTensorBasis<2>* > vectorOfBasis;
    gsTHBSplineBasis<2> copy(basis);
    vectorOfBasis.push_back(&copy);
    gsCompositeHBasis<2, real_t> periodicBasis(vectorOfBasis, get_periodic_topology(dir));
    
    gsMatrix<> global_coefs;
    periodicBasis.local_coef_to_global_coef(coefficients, global_coefs);
    
    gsMatrix<> result;
    periodicBasis.global_coef_to_local_coef(global_coefs, result);

    return result;
}


    
} //namespace gismo
