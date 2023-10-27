/** @file gsMultiPatchPreconditioners.h

    @brief Provides multi-patch preconditioners, particularly smoothers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsPreconditioner.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsCore/gsMultiPatch.h>
#include <gsMultiGrid/gsDomainDecomposition.h>

namespace gismo
{

// TODO: make transfer consistent with gsMultiGrid, i.e, transpose it _and_ change its type
//       from gsSparseMatrix<T> to gsSparseMatrix<T,RowMajor>

/// @brief Generic preconditioner which applies an arbitrary linear operator to the residual.
///
/// @ingroup Solvers
template<class T=real_t>
class gsAdditivePreconditionerOp : public gsPreconditionerOp<T>
{
    typedef gsPreconditionerOp<T> Base;
    typedef typename gsLinearOperator<T>::Ptr OpPtr;
public:

    /// Shared pointer
    typedef memory::shared_ptr<gsAdditivePreconditionerOp> Ptr;

    /// Unique pointer
    typedef memory::unique_ptr<gsAdditivePreconditionerOp> uPtr;

    /// Type of container of transfer matrices
    typedef std::vector< gsSparseMatrix<T> > TransferContainer;

    /// Type of container of linear operators
    typedef std::vector<OpPtr> OpContainer;

    /// Default Constructor
    gsAdditivePreconditionerOp() : m_A(), m_transfers(), m_ops(), m_damping(1) {}

    /// Constructor
    gsAdditivePreconditionerOp(const OpPtr& A,
                       TransferContainer transfers,
                       OpContainer ops,
                       T damping = 1.0)
    : m_A(A), m_transfers(give(transfers)), m_ops(give(ops)), m_damping(damping)
    {
        GISMO_ASSERT( m_ops.size() == m_transfers.size(), "Sizes do not agree" );
    }

    /// Make function
    static uPtr make(const OpPtr& A,
                     TransferContainer transfers,
                     OpContainer ops,
                     T damping = 1.0)
    { return uPtr( new gsAdditivePreconditionerOp( A, give(transfers), give(ops), damping ) ); }

    /// Make function
    static uPtr make(const OpPtr& A,
                     std::pair<
                         TransferContainer,
                         OpContainer
                     > transfer_and_ops,
                     T damping = 1.0)
    { return uPtr( new gsAdditivePreconditionerOp( A, give(transfer_and_ops.first), give(transfer_and_ops.second), damping ) ); }

    virtual void step(const gsMatrix<T>& f, gsMatrix<T>& x) const;

    index_t rows() const { return m_A->rows();}
    index_t cols() const { return m_A->cols();}
    OpPtr underlyingOp() const { return m_A; }

protected:
    OpPtr m_A;
    TransferContainer m_transfers;
    OpContainer m_ops;
    T m_damping;
    mutable gsMatrix<T> m_res, m_res_local, m_corr_local;

};

/// @brief Generic preconditioner which applies an arbitrary linear operator to the residual.
///
/// @ingroup Solvers
template<class T=real_t>
class gsMultiplicativePreconditionerOp : public gsPreconditionerOp<real_t>
{
    typedef gsPreconditionerOp<T> Base;
    typedef typename gsLinearOperator<T>::Ptr OpPtr;
public:

    /// Shared pointer
    typedef memory::shared_ptr<gsMultiplicativePreconditionerOp> Ptr;

    /// Unique pointer
    typedef memory::unique_ptr<gsMultiplicativePreconditionerOp> uPtr;

    /// Type of container of transfer matrices
    typedef std::vector< gsSparseMatrix<T> > TransferContainer;

    /// Type of container of linear operators
    typedef std::vector<OpPtr> OpContainer;

    /// Default Constructor
    gsMultiplicativePreconditionerOp() : m_A(), m_transfers(), m_ops(), m_damping(1) {}

    /// Constructor
    gsMultiplicativePreconditionerOp(const OpPtr& A,
                       TransferContainer transfers,
                       OpContainer ops,
                       T damping = 1.0)
    : m_A(A), m_transfers(give(transfers)), m_ops(give(ops)), m_damping(damping)
    {
        GISMO_ASSERT( m_ops.size() == m_transfers.size(), "Sizes do not agree" );
    }

    /// Make function
    static uPtr make(const OpPtr& A,
                     TransferContainer transfers,
                     OpContainer ops,
                     T damping = 1.0)
    { return uPtr( new gsMultiplicativePreconditionerOp( A, give(transfers), give(ops), damping ) ); }

    /// Make function
    static uPtr make(const OpPtr& A,
                     std::pair<
                         TransferContainer,
                         OpContainer
                     > transfer_and_ops,
                     T damping = 1.0)
    { return uPtr( new gsMultiplicativePreconditionerOp( A, give(transfer_and_ops.first), give(transfer_and_ops.second), damping ) ); }

    virtual void step(const gsMatrix<T>& f, gsMatrix<T>& x) const;

    index_t rows() const { return m_A->rows();}
    index_t cols() const { return m_A->cols();}
    OpPtr underlyingOp() const { return m_A; }

protected:
    OpPtr m_A;
    TransferContainer m_transfers;
    OpContainer m_ops;
    T m_damping;
    mutable gsMatrix<T> m_res, m_corr_local;

};

/// Wrapper for compatability (deprecated)
/// @ingroup Solvers
typedef gsAdditivePreconditionerOp<real_t> gsAdditiveSmoother;
/// Wrapper for compatability (deprecated)
/// @ingroup Solvers
typedef gsMultiplicativePreconditionerOp<real_t> gsMultiplicativeSmoother;


/// @brief To be used within gsAdditivePreconditionerOp or gsMultiplicativePreconditionerOp
///
/// @ingroup Solvers
template<class T>
std::vector< gsSparseMatrix<T> > getPatchwiseTransfers(
                            const gsMultiBasis<T>& mb,
                            const gsBoundaryConditions<T>& bc,
                            const gsOptionList& opt
                        );

/// @brief To be used within gsAdditivePreconditionerOp or gsMultiplicativePreconditionerOp
///
/// We clearly know what this does in 2D, but we might change our mind what it does in 3D.
///
/// @ingroup Solvers
template<class T>
std::vector< std::vector< gsSparseMatrix<T> > > getPiecewiseTransfers(
                            const gsMultiBasis<T>& mb,
                            const gsBoundaryConditions<T>& bc,
                            const gsOptionList& opt
                        );

/// @brief To be used within gsAdditivePreconditionerOp or gsMultiplicativePreconditionerOp
///
/// We clearly know what this does in 2D, but we might change our mind what it does in 3D.
///
/// @ingroup Solvers
template<class T>
std::pair< std::vector< gsSparseMatrix<T> >, std::vector< typename gsLinearOperator<T>::Ptr > > setupPiecewisePreconditioner(
                            gsSparseMatrix<T> A,
                            const std::vector< typename gsLinearOperator<T>::Ptr >& ops,
                            const gsMultiBasis<T>& mb,
                            const gsBoundaryConditions<T>& bc,
                            const gsOptionList& opt
                        );

template<class T>
GISMO_DEPRECATED std::pair< std::vector< gsSparseMatrix<T> >, std::vector< typename gsLinearOperator<T>::Ptr > > setupPiecewisePreconditioner(
                            gsSparseMatrix<T> A,
                            const std::vector< typename gsLinearOperator<T>::Ptr >& ops,
                            const gsMultiPatch<T>& mp,
                            const gsMultiBasis<T>& mb,
                            const gsBoundaryConditions<T>& bc,
                            const gsOptionList& opt
                        )
{ return setupPiecewisePreconditioner(A, ops, mb, bc, opt); }


/// @brief Set up local exact solvers
///
/// @ingroup Solvers
template<class T>
std::vector< typename gsLinearOperator<T>::Ptr> getLocalExactSolvers(
                    const gsSparseMatrix<T>& A,
                    const std::vector< gsSparseMatrix<T> >& transfers
                    );



// Highly experimental
template<class T>
std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<T> > > > constructPieces(
                                    const gsMultiBasis<T>& mb,
                                    const std::vector< gsVector<index_t> >& locals,
                                    index_t totalNumberDof,
                                    bool combineVertices = true
                                );

// Highly experimental
template<class T>
std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<T> > > > constructPieces(
                                    const gsMultiBasis<T>& mb,
                                    const gsDofMapper& dm,
                                    bool combineVertices = true
                                );

// Highly experimental
template<class T>
std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<T> > > > constructPieces(
                                    const gsMultiBasis<T>& mb,
                                    const gsBoundaryConditions<T>& bc,
                                    const gsOptionList& opt,
                                    bool combineVertices = true
                                );

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMultiPatchPreconditioners.hpp)
#endif
