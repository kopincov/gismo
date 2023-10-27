/** @file gsExperimentalSweeps.h

    @brief Collection of some preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{
    
GISMO_EXPORT void dampedRichardsonSweepBoundary(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau, int numBdNodes);
GISMO_EXPORT void kaczmarzSweepBoundary(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau, int numBdNodes);
GISMO_EXPORT void dampedPreRichardsonSweep(const gsSparseMatrix<real_t>& A, const gsSparseMatrix<real_t>& P, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau);
GISMO_EXPORT void dampedPreJacobiSweep(const gsSparseMatrix<real_t>& A, const gsSparseMatrix<real_t>& P, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(1.0/2.0));
GISMO_EXPORT void preGaussSeidelSweep(const gsSparseMatrix<real_t>& A, const gsSparseMatrix<real_t>& P, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau, bool reverse);

} // namespace gismo
