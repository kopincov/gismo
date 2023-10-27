/** @file gsParallelCG.h

    @brief Parallel Implementation of the GMRes method

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
//#include <gsSolver/gsSolverUtils.h>
#include <gsIETI/gsDistributedOperator.h>

namespace gismo {

/** The conjugate gradient implementation from Eigen, adapted to allow for more
 *  general preconditioners and better iteration control. Also capable of using
 *  a gsLinearOperator as matrix.
 *
 *
 *  Only implemented for single right hand side!
 */
template<class T>
class gsParallelCG : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;
    
    typedef typename Base::VectorType  VectorType;
    
    typedef typename Base::LinOpPtr LinOpPtr;

    typedef typename memory::shared_ptr<gsParallelCG<T> > Ptr;
    typedef typename memory::unique_ptr<gsParallelCG<T> > uPtr;

    /// Contructor. See gsIterativeSolver for details.
    template< typename OperatorType >
    gsParallelCG( const OperatorType& mat, const LinOpPtr & precond = LinOpPtr(), MPI_Comm comm = gsMpi::worldComm())
    : Base(mat, precond),
      m_calcEigenvals(false), m_eigsAreCalculated(false), m_comm(comm)
    { }

    /// @brief Make function using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat, const LinOpPtr & precond = LinOpPtr(), MPI_Comm comm = gsMpi::worldComm() )
    { return uPtr( new gsParallelCG<T>(mat, precond, comm) ); }
    
    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addSwitch("CalcEigenvalues", "Additionally to solving the system,"
                      " CG computes the eigenvalues of the Lanczos matrix", false );
        return opt;
    }
    
    gsParallelCG<T>& setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_calcEigenvals = opt.askSwitch("CalcEigenvalues", m_calcEigenvals);
        return *this;
    }
    

    bool initIteration(const VectorType& rhs, VectorType& x);
/*
    void solve(const VectorType& rhs, VectorType& x)
        {
      //  gsInfo<<"rank: "<<m_comm.rank()<<" before init in PCG\n ";
            if(initIteration(rhs, x)==true)
            {
              //  gsInfo<<"rank: "<<m_comm.rank()<<" after init in PCG\n ";
                m_error = 0;
                if(m_calcEigenvals)
                    m_eigsAreCalculated =  true;
                return;
            }
            while(m_num_iter < m_max_iters)
            {
               // gsInfo<<"rank: "<<m_comm.rank()<<" -> "<<x.transpose(    Created on: 2016-04-22)<<"\n\n";
                m_num_iter++;
                if (step(x))
                    break;
            }
            m_error = math::sqrt(residualNorm2 / rhsNorm2);
        }
*/
    bool step( VectorType& x );

    /// @brief specify if you want to store data for eigenvalue estimation
    /// @param flag true stores the coefficients of the lancos matrix, false not.
    void setCalcEigenvalues(bool flag) {m_calcEigenvals = flag;}

    /// @brief returns the condition number of the (preconditioned) system matrix
    T getConditionNumber();

    /// @brief returns the eigenvalues of the Lanczos matrix
    void getEigenvalues(gsMatrix<T>& eigs);

    /// @ brief implementation of the dot product
    T dot(const VectorType & a, const VectorType & b);

     /// @ brief implementation of the dot product for two vectors
    void dot(const VectorType &a, const VectorType &b, const VectorType &c, const VectorType &d, T* res);

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsParallelCG\n";
        return os;
    }

    //  /// @ brief implementation of the norm
    //  T norm(const VectorType &a, const gsLinearOperator& precond);
    
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

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;

    VectorType z, z_acc, tmp, p,p_acc;
    VectorType residual, residual_acc;
    T absNew, residualNorm2, threshold, rhsNorm2;

    bool m_calcEigenvals;
    bool m_eigsAreCalculated;

    std::vector<T> delta, gamma;

    gsMpiComm m_comm;
};


} // namespace gismo

#endif

#ifndef GISMO_BUILD_LIB
#include <gsIETI/gsParallelCG.hpp>
#endif
