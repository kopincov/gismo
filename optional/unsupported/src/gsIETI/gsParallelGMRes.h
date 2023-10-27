/** @file gsParallelGMRes.h

    @brief Parallel Implementation of the GMRes method


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on:
*/

#pragma once

#include <gsCore/gsConfig.h>
#ifdef GISMO_WITH_MPI


#include <gsIETI/gsDistributedOperator.h>
#include <gsMpi/gsMpi.h>
#include <gsMpi/gsMpiComm.h>
#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

//In contrast to the serial class, we use the maxInter to allocate memory for the orthogonal vectors

/// @brief The generalized minimal residual (GMRES) method.
///
/// \ingroup Solver

template<typename T = real_t>
class gsParallelGMRes: public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;

    typedef typename Base::VectorType VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;
    typedef typename gsDistributedOperator<T>::Ptr DistOpPtr;

    typedef memory::shared_ptr<gsParallelGMRes> Ptr;
    typedef memory::unique_ptr<gsParallelGMRes> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    gsParallelGMRes( const OperatorType& mat, const LinOpPtr & precond = LinOpPtr(), MPI_Comm comm = gsMpi::worldComm()  )
    : Base(mat, precond), m_comm(comm) {}

    /// @brief Make function using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat,
                      const LinOpPtr& precond = LinOpPtr(), gsMpiComm comm = gsMpi::worldComm() )
    { return uPtr( new gsParallelGMRes(mat, precond, comm) ); }

    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );
    void finalizeIteration( VectorType& x );

private:

    /// Solves the Upper triangular system Ry = gg
    /// and stores the solution in the private member y.
    void solveUpperTriangular(const gsMatrix<T> & R, const gsMatrix<T> & gg)
    {
       y = R.template triangularView<Eigen::Upper>().solve(gg);
    }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsParallelGMRes\n";
        return os;
    }

    /// @ brief implementation of the dot product
    T dot(const VectorType & a, const VectorType & b);

     /// @ brief implementation of the dot product for two vectors
    void dot(const VectorType &a, const VectorType &b, const VectorType &c, const VectorType &d, T* res);

private:
    void accumulate(const VectorType& input, VectorType& accumulated, const gsLinearOperator<T>& op) const
    {
        const gsDistributedOperator<T>* parOp = dynamic_cast<const gsDistributedOperator<T>*>(&op);
        GISMO_ASSERT(parOp!=NULL, "your operator is not derived from gsDistributedOperator");
        parOp->accumulate(input, accumulated);
    }

    void postAccumulate(const gsLinearOperator<T>& op)
    {
        const gsDistributedOperator<T>* parOp = dynamic_cast<const gsDistributedOperator<T>*>(&op);
        GISMO_ASSERT(parOp!=NULL, "your operator is not derived from gsDistributedOperator");
        parOp->postAccumulate();
    }

    void startAccumulate(const  VectorType & input, const gsLinearOperator<T>& op)
    {
        const gsDistributedOperator<T>* parOp = dynamic_cast<const gsDistributedOperator<T>*>(&op);
        GISMO_ASSERT(parOp!=NULL, "your operator is not derived from gsDistributedOperator");
        parOp->startAccumulate(input);
    }


    void finishAccumulate(VectorType & result, const gsLinearOperator<T>& op)
    {
        const gsDistributedOperator<T>* parOp = dynamic_cast<const gsDistributedOperator<T>*>(&op);
        GISMO_ASSERT(parOp!=NULL, "your operator is not derived from gsDistributedOperator");
        parOp->finishAccumulate(result);
    }

    void distribute(const VectorType& input, VectorType& distributed, const gsLinearOperator<T>& op) const
    {
        const gsDistributedOperator<T>* parOp = dynamic_cast<const gsDistributedOperator<T>*>(&op);
        GISMO_ASSERT(parOp!=NULL, "your operator is not derived from gsDistributedOperator");
        parOp->distribute(input, distributed);
    }


private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;


    gsMatrix<T> tmp,tmp_acc, g, g_tmp, h_tmp, y, w; // w_acc;
    gsMatrix<T> residual,residual_acc;
    gsMatrix<T> H_prev, H, Omega, Omega_prev, Omega_tmp, Omega_prev_tmp;
    std::vector<gsMatrix<T> > v_acc;
    T beta;

    gsMpiComm m_comm;
};

} // namespace gismo

#endif

#ifndef GISMO_BUILD_LIB
#include <gsIETI/gsParallelGMRes.hpp>
#endif
