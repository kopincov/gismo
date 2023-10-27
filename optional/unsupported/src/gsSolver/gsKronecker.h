/** @file gsKronecker.h

    @brief Provides functions and classes for working with Kronecker products of matrices and operators.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, S. Takacs
*/
#pragma once

#include <gsSolver/gsKroneckerOp.h>

namespace gismo
{

/// @brief Compute the Kronecker product of a vector of matrices as a matrix
///
/// \ingroup Solver
template <typename MatrixType>
MatrixType kroneckerProduct(const std::vector< MatrixType >& matrices)
{
    if ( matrices.size() == 0 )
    {
        MatrixType result;
        return result;
    }
    else if ( matrices.size() == 1 )
        return matrices[0];
    else
    {
        MatrixType result = matrices[0].kron(matrices[1]);
        for ( unsigned i = 2; i<matrices.size(); ++ i )
            result = result.kron(matrices[i]);
        return result;
    }
}



}
