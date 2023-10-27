/** @file gsPatchPreconditionersCreator2.h

    @brief Provides robust preconditioners for single patch geometries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, C. Hofreither
*/
#pragma once

#include <gsSolver/gsPatchPreconditionersCreator.h>

namespace gismo
{

/// @brief Provides robust preconditioners for single patch geometries.
///
/// This class extends gsPatchPreconditionersCreator with devel stuff
///
/// @ingroup Solver
template<typename T>
class gsPatchPreconditionersCreator2
{
    typedef typename gsLinearOperator<T>::uPtr OpUPtr;
    typedef typename gsLinearOperator<T>::Ptr  OpPtr;
public:

    /// @brief Like \a gsPatchPreconditionersCreator<T>::fastDiagonalization, but with geometry.
    static OpUPtr fastDiagonalizationWithGeoOp(gsGeometry<T>& geo, gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, T alpha = 1);

};

} // namespace gismo
