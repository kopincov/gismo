/** @file gsILUTPreconditioners.h

    @brief Collection of some simple preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsPreconditioner.h>

namespace gismo
{


/// @brief  ILUT (incomplete LU with thresholding) smoother
/// Assumes A is symmetric (not needed)!
///
/// \ingroup Solvers
template<class T = real_t>
class GISMO_EXPORT gsILUTPreconditionerOp : public gsPreconditionerOp<T>
{
public:
    /// Shared pointer for gsILUTPreconditionerOp
    typedef memory::shared_ptr< gsILUTPreconditionerOp > Ptr;

    /// Unique pointer for gsILUTPreconditionerOp
    typedef memory::unique_ptr< gsILUTPreconditionerOp > uPtr;

    /// Base class
    typedef gsPreconditionerOp<T> Base;

    /// @brief Constructor with given matrix
    gsILUTPreconditionerOp(const gsSparseMatrix<T>& A, const gsSparseMatrix<T>& K)
        : m_A(A)
    {
        m_ilu.setFillfactor(1);
        m_ilu.compute(K);
    }

    /// @brief Make function returning smart pointer
    static uPtr make(const gsSparseMatrix<T>& A, const gsSparseMatrix<T>& K)
        { return uPtr( new gsILUTPreconditionerOp( A, K ) ); }

    void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        x += m_ilu.solve( rhs - m_A * x ).eval();
    }

    // We use our own apply implementation as we can save one multiplication. This is important if the number
    // of sweeps is 1: Then we can save *all* multiplications.
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        // For the first sweep, we do not need to multiply with the matrix
        x = m_ilu.solve( input ).eval();

        for (index_t k = 1; k < m_num_of_sweeps; ++k)
            x += m_ilu.solve( input - m_A * x ).eval();

    }

    index_t rows() const {return m_A.rows();}
    index_t cols() const {return m_A.cols();}
    typename gsLinearOperator<T>::Ptr underlyingOp() const { return makeMatrixOp(m_A); }

private:
    gsSparseMatrix<T> m_A;
    Eigen::IncompleteLUT<T> m_ilu;
};


} // namespace gismo

