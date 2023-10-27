/** @file gsScatteredHFittingCoefficientsPeriodic.h

    @brief Extension of Scattered fitting with hierarchical splines to periodic data.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore, D. Mokris
*/

#pragma once

#include <gismo.h>
#include <gsSmoothPatches/gsCompositeHBasis.h>
#include <gsHSplines/gsScatteredHFittingImprovementsGeneralUtilities.h>
#include <gsHSplines/gsScatteredHFittingCoefficientsStiffness.h>

namespace gismo
{


/** \brief Topology of the interface for periodic surfaces
 * \param[in]  dir periodic direction: 0 for \a u direction, 1 for \a v direction
 * \param[out] topology of the interface 
 */
gsBoxTopology get_periodic_topology(bool dir)
{
    gsBoxTopology topology(2, 1);
    gsVector<index_t> dirMap(2);
    dirMap << 0, 1; //the coordinate direction 0 is mapped to coordinate direction dirMap[0] and the coordinate direction 1 with dirMap[1].
    gsVector<bool> dirOrientation(2);
    dirOrientation << true, true;    
    
    //patchSide ps1, ps2;
    patchSide ps1(0, dir ? boundary::down : boundary::left );
    patchSide ps2(0, dir ? boundary::up   : boundary::right);

    boundaryInterface boundInter(ps1, ps2, true);
    
    topology.addInterface(boundInter);
    return topology;
    
}


/** \brief computation of coefficients of gsTHBSplineBasis \em thbBasis for periodic fitting with smoothing term.
 */
gsMatrix<> get_coefficients_periodic(const gsTHBSplineBasis<2>& thbBasis, 
                                     const gsMatrix<>& parameters, 
                                     const gsMatrix<>& points,
                                     int iteration,
                                     const int cp_magnitude,
                                     const real_t hole,
                                     const real_t lambda,
                                     const int dir,
                                     const int circle,
                                     std::vector<index_t>& enlargement)
{
    
    //coefficients of the THB fit.
    gsMatrix<> coefficients(thbBasis.size(), points.rows());
    for (index_t k = 0; k < thbBasis.size(); k++)
    {
        gsMatrix<> coefficient = get_coefficient_not_scaled_stiff(thbBasis, parameters, points, iteration, k, cp_magnitude, hole, lambda, circle, enlargement);        
        coefficients.row(k) = coefficient.row(0);
    }
    
    
    // composite basis
    std::vector< gsHTensorBasis<2>* > vectorOfBasis;
    gsTHBSplineBasis<2> copy(thbBasis);
    vectorOfBasis.push_back(&copy);
    gsCompositeHBasis<2, real_t> periodicBasis(vectorOfBasis, get_periodic_topology(dir));
    
    gsMatrix<> global_coefs;
    periodicBasis.local_coef_to_global_coef(coefficients, global_coefs);
    
    gsMatrix<> result;
    periodicBasis.global_coef_to_local_coef(global_coefs, result);

    return result;

}

} // namespace gismo
