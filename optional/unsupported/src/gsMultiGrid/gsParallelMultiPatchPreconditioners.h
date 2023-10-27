/** @file gsParallelMultiPatchSmoothers.h

    @brief Provides Multigrid smoothers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/
#pragma once

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>

#include <gsSolver/gsPreconditioner.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsIETI/gsDistributedOperator.h>
#include <gsIETI/gsParallelOperator.h>

namespace gismo
{

/// @brief Generic smoother which applies an arbitrary linear operator to the residual. The linear operators
/// are assumed to be somewhat local. This class combines them.
///
/// \ingroup Solvers
template<class T=real_t>
class gsParallelAdditivePreconditionerOp : public gsPreconditionerOp<T>
{
    typedef gsPreconditionerOp<T> Base;
    typedef typename gsLinearOperator<T>::Ptr OpPtr;
    typedef typename gsDistributedOperator<T>::Ptr DOpPtr;
public:

    /// Shared pointer
    typedef memory::shared_ptr<gsParallelAdditivePreconditionerOp> Ptr;

    /// Unique pointer
    typedef memory::unique_ptr<gsParallelAdditivePreconditionerOp> uPtr;

    /// Type of container of transfer matrices
    typedef std::vector< gsSparseMatrix<T> > TransferContainer;

    /// Type of container of linear operators
    typedef std::vector<typename gsLinearOperator<T>::Ptr > OpContainer;


    /// Default Constructor
    gsParallelAdditivePreconditionerOp() : m_A(), m_restriction(), m_prolongation(), m_ops(), m_damping(1.) {}

    /// Constructor
    gsParallelAdditivePreconditionerOp(const DOpPtr& A,
                               const TransferContainer& transfers,
                               const OpContainer& ops,
                               real_t damping = 1.0)
        : m_A(A), m_restriction(transfers.size()), m_prolongation(transfers.size()), m_ops(ops), m_damping(damping)
    {
        GISMO_ASSERT( m_ops.size() == transfers.size(), "Sizes do not agree" );

        for ( size_t i=0; i<m_ops.size(); ++i )
        {
            memory::shared_ptr< gsSparseMatrix<T> > transferMatrix = const_cast<gsSparseMatrix<T>*>(&transfers[i])->moveToPtr();
            // TODO: this setup is opposite to that of gsMultiGridOp
            // TODO: if we swap restrict and prolong, we also should swap gsSparseMatrix<T> with gsSparseMatrix<T,RowMajor>
            m_restriction[i] = makeMatrixOp(transferMatrix);
            m_prolongation[i] = makeMatrixOp(transferMatrix->transpose()); // note that this works as we know that
            // m_prolong does not get destroyed before
            // m_restrict, so the shared pointer stored in m_prolong
            // will make sure that the matrix will not get destroyed
        }
    }

    /// Make function
    static uPtr make(const DOpPtr& A,
                     const TransferContainer& transfers,
                     const OpContainer& ops,
                     real_t damping = 1.0)
    { return uPtr( new gsParallelAdditivePreconditionerOp( A, transfers, ops, damping ) ); }

    /// Make function
    static uPtr make(const DOpPtr& A,
                     const std::pair<
                         TransferContainer,
                         OpContainer
                     >& transfer_and_ops,
                     real_t damping = 1.0)
    { return uPtr( new gsParallelAdditivePreconditionerOp( A, transfer_and_ops.first, transfer_and_ops.second, damping ) ); }

    virtual void step(const gsMatrix<T>& f, gsMatrix<T>& x) const;

    index_t rows() const { return m_A->rows();}
    index_t cols() const { return m_A->cols();}
    OpPtr underlyingOp() const { return m_A; }

protected:
    DOpPtr m_A;
    OpContainer m_restriction;
    OpContainer m_prolongation;
    OpContainer m_ops;
    T m_damping;
    mutable gsMatrix<T> m_res, m_corr, m_corr_global, m_local_res, m_corr_global_buffer ;

};

/*
/// @brief To be used within gsAdditiveSmoother or gsMultiplicativeSmoother
///
/// \ingroup Solvers
template <typename T>
std::vector< gsSparseMatrix<T> > getPatchwiseTransfers(
        const gsMultiBasis<T>& mb,
        const gsBoundaryConditions<T>& bc,
        const gsOptionList& opt,
        const gsSortedVector<size_t>& myPatches,
        const gsMpiComm comm
        );
*/
/// @brief To be used within gsAdditiveSmoother or gsMultiplicativeSmoother
///
/// We clearly know what this does in 2D, but we might change our mind what it does in 3D.
///
/// \ingroup Solvers
template <typename T>
std::vector< std::vector< gsSparseMatrix<T> > > getPiecewiseTransfers(
        const gsMultiPatch<T>& mp,
        const gsMultiBasis<T>& mb,
        const gsBoundaryConditions<T>& bc,
        const gsOptionList& opt,
        const gsSortedVector<size_t>& myPatches,
        const gsParallelGlobalLocalHandler& A_handler,
        const gsMpiComm comm,
        gsSortedVector<index_t>& interfaceDofs
        );


/// @brief To be used within gsAdditiveSmoother or gsMultiplicativeSmoother
///
/// We clearly know what this does in 2D, but we might change our mind what it does in 3D.
///
/// \ingroup Solvers
template <typename T>
std::pair< std::vector< gsSparseMatrix<T> >, std::vector< typename gsLinearOperator<T>::Ptr > > setupPiecewisePreconditioner(
        const gsParallelOperator<T>& A,
        const gsParallelGlobalLocalHandler& A_handler,
        const std::vector< typename gsLinearOperator<T>::Ptr >& ops,
        const gsMultiPatch<T>& mp,
        const gsMultiBasis<T>& mb,
        const gsBoundaryConditions<T>& bc,
        const gsOptionList& opt,
        const gsSortedVector<size_t>& myPatches,
        const gsMpiComm comm
        );

/*
/// @brief Set up local exact solvers
///
/// \ingroup Solvers
template <typename T>
std::vector<gsLinearOperator<T>::Ptr> getLocalExactSolvers(
        const gsSparseMatrix<T>& A,
        const std::vector< gsSparseMatrix<T> >& transfers
        );

*/

} // namespace gismo

#endif
