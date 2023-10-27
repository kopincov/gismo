/** @file gsParallelGradientMethod.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2017-11-22
*/


#pragma once
#include <gsCore/gsConfig.h>
#ifdef GISMO_WITH_MPI

#include <gsSolver/gsIterativeSolver.h>
#include <gsIETI/gsDistributedOperator.h>
#include <gsMpi/gsMpiComm.h>
#include <gsMpi/gsMpi.h>

namespace gismo {

/// @brief The Gradient iteration method.
///
/// \ingroup Solver

template<typename T>
class gsParallelGradientMethod : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;

    typedef typename  Base::VectorType VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;
    typedef typename gsDistributedOperator<T>::Ptr DistOpPtr;

    typedef memory::shared_ptr<gsParallelGradientMethod> Ptr;
    typedef memory::unique_ptr<gsParallelGradientMethod> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    /// @param damping The damping parameter (step width)
    template< typename OperatorType >
    explicit gsParallelGradientMethod( const OperatorType& mat,
                                  const LinOpPtr& precond = LinOpPtr(), T damping = 1, gsMpiComm comm = gsMpi::worldComm() )
    : Base(mat, precond), m_damping(damping), m_comm(comm) {}
    
    /// @brief Make function using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    /// @param damping The damping parameter (step width)
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat,
                      const LinOpPtr& precond = LinOpPtr(), T damping = 1, gsMpiComm comm = gsMpi::worldComm() )
    { return uPtr( new gsParallelGradientMethod(mat, precond, damping, comm) ); }

    /// @brief Returns the chosen damping parameter (step width)
    T getDamping() { return m_damping; }

    /// @brief Set the damping parameter (step width)
    void setDamping(T damping) { m_damping = damping; }

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addReal("Damping", "Damping (step width)", 1 );
        return opt;
    }

    /// @brief Set the options based on a gsOptionList
    gsParallelGradientMethod& setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_damping = opt.askReal("Damping", m_damping);
        return *this;
    }

    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsParallelGradientMethod\n";
        return os;
    }

private:
    void accumulate(const VectorType& input, VectorType& accumulated, const gsLinearOperator<T>& op) const
    {
        const gsDistributedOperator<T>* parOp = dynamic_cast<const gsDistributedOperator<T>*>(&op);
        GISMO_ASSERT(parOp!=NULL, "your operator is not derived from gsDistributedOperator");
        parOp->accumulate(input, accumulated);
    }

    void distribute(const VectorType& input, VectorType& distributed, const gsLinearOperator<T>& op) const
    {
        const gsDistributedOperator<T>* parOp = dynamic_cast<const gsDistributedOperator<T>*>(&op);
        GISMO_ASSERT(parOp!=NULL, "your operator is not derived from gsDistributedOperator");
        parOp->distribute(input, distributed);
    }

    /// @ brief implementation of the dot product
    T dot(const VectorType & a, const VectorType & b);

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;

    T m_damping;

    VectorType m_res, m_res_acc;
    VectorType m_update_acc;
    VectorType m_tmp;

    gsMpiComm m_comm;
};


} // namespace gismo

#endif

#ifndef GISMO_BUILD_LIB
#include <gsIETI/gsParallelGradientMethod.hpp>
#endif
