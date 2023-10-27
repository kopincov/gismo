/** @file gsBlockSmoother.h

    @brief Provides multigrid block smoothers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

#pragma once

#include <gsCore/gsDebug.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsSolver/gsPreconditioner.h>
#include <gsMultiGrid/gsBlockOrganizer.h>

namespace gismo
{

/// @brief Preforms a block Gauss-Seidel on the degrees of freedom in DoFs.
/// \ingroup Solver
template<typename T>
void gaussSeidelSingleBlock(const gsSparseMatrix<T>& A, gsMatrix<T>& x, const gsMatrix<T>& f, gsVector<index_t>& DoFs)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");

    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    //Sorting from lowest to highest
    DoFs.sortByColumn(0);
    const index_t size = DoFs.rows();

    GISMO_ASSERT ( DoFs(size-1)< A.cols(), "The given DoF is higher than the size of the matrix");

    gsMatrix<T> Dblock(size, size);
    gsMatrix<T> residual(size, 1);
    for (int i = 0; i< size; i++)
    {
        //Symmetry is assumed here!
        residual(i,0) = f(DoFs(i),0) - (A.innerVector(DoFs(i)).transpose()*x).value();
        for (int j = 0; j< size; j++)
            Dblock(i,j) = A.coeff(DoFs(i), DoFs(j));
    }
    //Multiply residual with the inverse of the diagonal block
    residual = Dblock.inverse()*residual;

    //Update
    for (int i = 0; i< size; i++)
       x(DoFs(i),0) += residual(i,0);
}


/// @brief Block Gauss-Seidel smoother
/// Assumes A is symmetric (not needed)!
///
/// \ingroup Solvers
class GISMO_EXPORT gsGaussSeidelBlockSmoother : public gsPreconditionerOp<real_t>
{
public:

    typedef memory::shared_ptr<gsGaussSeidelBlockSmoother> Ptr;
    typedef memory::unique_ptr<gsGaussSeidelBlockSmoother> uPtr;
    typedef gsPreconditionerOp<real_t> Base;

    /// @brief Constructor
    ///
    /// Each elemet in \a blockInfo contains a vector with the indices of the DoFs with are grouped
    /// together in one block
    gsGaussSeidelBlockSmoother(gsSparseMatrix<> A, std::vector<gsVector<index_t> > &blockInfo)
        : m_APtr(A.moveToPtr()), m_blockInfo(blockInfo)
    { }

    gsGaussSeidelBlockSmoother(const gsSparseMatrix<>::Ptr& APtr, std::vector<gsVector<index_t> > &blockInfo)
        :  m_APtr(APtr), m_blockInfo(blockInfo)
    {}

    static uPtr make(gsSparseMatrix<> A, std::vector<gsVector<index_t> > &blockInfo)
        { return uPtr( new gsGaussSeidelBlockSmoother(give(A), blockInfo ) ); }

    static uPtr make(const gsSparseMatrix<>::Ptr & APtr, std::vector<gsVector<index_t> > &blockInfo)
        { return uPtr( new gsGaussSeidelBlockSmoother( APtr , blockInfo) ); }

    void step(const gsMatrix<>& f, gsMatrix<>& x) const;
    void stepT(const gsMatrix<>& f, gsMatrix<>& x) const;

    index_t rows() const { return m_APtr->rows();}
    index_t cols() const { return m_APtr->cols();}

    gsLinearOperator<>::Ptr underlyingOp() const { return makeMatrixOp(m_APtr); }

private:
    gsSparseMatrix<>::Ptr m_APtr;
    mutable std::vector<gsVector<index_t> > m_blockInfo;
};

} // namespace gismo
