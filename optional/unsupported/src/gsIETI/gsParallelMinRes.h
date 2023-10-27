/** @file gsParallelMinRes.h

    @brief Parallel implementation of the MinRes method

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/


#pragma once

#include <gsCore/gsConfig.h>
#ifdef GISMO_WITH_MPI

#include <gsSolver/gsIterativeSolver.h>
#include <gsMpi/gsMpi.h>
#include <gsMpi/gsMpiComm.h>
#include <gsIETI/gsDistributedOperator.h>

namespace gismo {

/** @brief The minimal residual (MinRes) method.
  *
  * \ingroup Solver
  */
template<class T = real_t>
class gsParallelMinRes : public gsIterativeSolver<T>
{

public:
    typedef gsIterativeSolver<T> Base;

    typedef typename  Base::VectorType VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsParallelMinRes> Ptr;
    typedef memory::unique_ptr<gsParallelMinRes> uPtr;

    /// Contructor. See gsIterativeSolver for details.
    template< typename OperatorType >
    gsParallelMinRes( const OperatorType& mat, const LinOpPtr & precond = LinOpPtr(), MPI_Comm comm = gsMpi::worldComm())
        : Base(mat, precond), m_inexact_residual(false), m_comm(comm){}

    /// @brief Make function using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat, const LinOpPtr & precond = LinOpPtr(), MPI_Comm comm = gsMpi::worldComm() )
    { return uPtr( new gsParallelMinRes(mat, precond, comm) ); }

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addSwitch( "InexactResidual", "If true, the residual is estimated, not accurately computed.", false );
        return opt;
    }

    gsParallelMinRes& setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_inexact_residual = opt.askSwitch("InexactResidual", m_inexact_residual);
        return *this;
    }


    bool initIteration(const VectorType& rhs, VectorType& x);
    void finalizeIteration( VectorType& x );

    bool step( VectorType& x );

    /// @brief If true, the residual is estimated, not accurately computed.
    void setInexactResidual( bool flag )     { m_inexact_residual = flag; }

    /// @ brief implementation of the dot product
    T dot(const VectorType & a, const VectorType & b);

     /// @ brief implementation of the dot product for two vectors
    void dot(const VectorType &a, const VectorType &b, const VectorType &c, const VectorType &d, T* res);

    /// @ brief implementation of the dot product for three
    void dot(const VectorType &a, const VectorType &b, const VectorType &c, const VectorType &d, const VectorType &e, const VectorType &f, T* res);


    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsParallelMinRes\n";
        return os;
    }

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

    gsMatrix<T> negResidual,
                     vPrev, v, vNew,v_acc,
                     wPrev_acc, w_acc, wNew_acc, AwPrev, Aw, AwNew,
                     zNew_acc, z_acc, Az;

    T eta,
           gammaPrev, gamma, gammaNew,
           sPrev, s, sNew,
           cPrev, c, cNew;

    bool m_inexact_residual;

    gsMpiComm m_comm;
};



} // namespace gismo

#endif

#ifndef GISMO_BUILD_LIB
#include <gsIETI/gsParallelMinRes.hpp>
#endif
