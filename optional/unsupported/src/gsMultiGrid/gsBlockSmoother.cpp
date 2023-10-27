/** @file gsBlockSmoother.cpp

    @brief Provides multigrid block smoothers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/


#include <gsMultiGrid/gsBlockSmoother.h>
#include <gsSolver/gsSimplePreconditioners.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

void gsGaussSeidelBlockSmoother::step(const gsMatrix<real_t>& f, gsMatrix<real_t>& x) const
{
    assert( m_APtr->rows() == x.rows() && x.rows() == f.rows() );
    assert( m_APtr->cols() == m_APtr->rows() && x.cols() == 1 && f.cols() == 1);
    assert( m_blockInfo.size() > 0);

    const index_t numberOfBlocks = m_blockInfo.size();
    for (index_t k = 0; k < numberOfBlocks; k++)
        gaussSeidelSingleBlock(*m_APtr, x, f, m_blockInfo[k]);
}

void gsGaussSeidelBlockSmoother::stepT(const gsMatrix<real_t>& f, gsMatrix<real_t>& x) const
{
    assert( m_APtr->rows() == x.rows() && x.rows() == f.rows() );
    assert( m_APtr->cols() == m_APtr->rows() && x.cols() == 1 && f.cols() == 1);
    assert( m_blockInfo.size() > 0);

    const index_t numberOfBlocks = m_blockInfo.size();
    for (index_t k = 0; k < numberOfBlocks; k++)
        gaussSeidelSingleBlock(*m_APtr, x, f, m_blockInfo[numberOfBlocks - 1 - k]);
}


} // namespace gismo
