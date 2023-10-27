/** @file gsIterativeSolverOp.h

    @brief Wraps an iterative solver in a gsLinearOperator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#include <gsCore/gsConfig.h>
#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

/*
template<class IterativeSolver>
class gsIterativeSolverOp :public gsLinearOperator<typename IterativeSolver::ScalarType>
{
public:
    typedef T IterativeSolver::ScalarType;

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<gsIterativeSolverOp> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<gsIterativeSolverOp> uPtr;

    gsIterativeSolverOp(const typename gsLinearOperator<T>::Ptr& op, const typename gsLinearOperator<T>::Ptr& prec, gsOptionList options)
        : m_op(op) , m_prec(prec), m_options(options)
    {
        m_solver = IterativeSolver::make(*m_op,*m_prec);
        m_solver->setOptions(m_options);
    }

    static uPtr make(const typename gsLinearOperator<T>::Ptr& op, const typename gsLinearOperator<T>::Ptr& prec, gsOptionList options)
        { return memory::make_unique( new gsIterativeSolverOp(op,prec,options) ); }

    void apply(const gsMatrix<T> &input, gsMatrix<T> &x) const
    {
        //otherwise use x as initial guess
        if(x.rows() != input.rows() || x.cols() != input.cols())
            x.setZero(input.rows(),input.cols());

        if(input.cols()==1)
            m_solver.solve(input,x);
        else //one additional copy is required for multiple right hand sides
        {
            gsMatrix<T> temp;
            for(index_t i=0; i<input.cols();++i)
            {
                m_solver.solve(input.col(i),temp);
                x.col(i) = temp;
            }
        }
    }

    /// Returns the number of rows of the operator
    virtual index_t rows() const {return m_op->rows();}

    /// Returns the number of columns of the operator
    virtual index_t cols() const {return m_op->cols();}


   const typename gsLinearOperator<T>::Ptr& m_op;
   const typename gsLinearOperator<T>::Ptr& m_prec;
   gsOptionList m_options;
   typename IterativeSolver::Ptr m_solver;
};
*/


#ifdef GISMO_WITH_MPI
#include <gsIETI/gsDistributedOperator.h>
#include <gsIETI/gsParallelOperator.h>

template<class IterativeSolver, class T = real_t>
class gsDistributedIterativeSolverOp :public gsDistributedOperator<T>
{
public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<gsDistributedIterativeSolverOp> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<gsDistributedIterativeSolverOp> uPtr;

    gsDistributedIterativeSolverOp(const typename gsConnectionHandler<T>::Ptr& connHandler, const typename gsDistributedOperator<T>::Ptr& op, const typename gsDistributedOperator<T>::Ptr& prec, gsOptionList options)
        : m_connHandler(connHandler), m_op(op) , m_prec(prec), m_options(options)
    {
        m_solver = IterativeSolver::make(m_op,m_prec,m_connHandler->getComm());
        m_solver->setOptions(m_options);
    }

    static uPtr make(const typename gsConnectionHandler<T>::Ptr& connHandler, const typename gsDistributedOperator<T>::Ptr& op, const typename gsDistributedOperator<T>::Ptr& prec, gsOptionList options)
        { return memory::make_unique( new gsDistributedIterativeSolverOp(connHandler,op,prec,options) ); }


    void apply(const gsMatrix<T> &input, gsMatrix<T> &x) const
    {
        //otherwise use x as initial guess
        if(x.rows() != input.rows() || x.cols() != input.cols())
            x.setZero(input.rows(),input.cols());

        if(input.cols()==1)
            m_solver->solve(input,x);
        else //one additional copy is required for multiple right hand sides
        {
            gsMatrix<T> temp;
            for(index_t i=0; i<input.cols();++i)
            {
                m_solver->solve(input.col(i),temp);
                x.col(i) = temp;
            }
        }
        gsInfo<<"gsDistributedIterativeSolverOp numIT: "<<m_solver->iterations()<<"  -- error: "<<m_solver->error()<<"\n";
    }

    /// Returns the number of rows of the operator
    virtual index_t rows() const {return m_op->rows();}

    /// Returns the number of columns of the operator
    virtual index_t cols() const {return m_op->cols();}

    virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
    {
        distributed = input;
        m_connHandler->distributeAccumulatedVector(distributed);
    }

    virtual void accumulate(const gsMatrix<T> &input, gsMatrix<T> &accumulated) const
    {
        accumulated = input;
        m_connHandler->accumulateDistributedVector(accumulated);
    }

    virtual void postAccumulate() const
    {
        m_connHandler->postAccumulate();
    }
    virtual void startAccumulate(const  gsMatrix<T> & input) const
    {
        m_connHandler->startAccumulate(input);
    }
    virtual void finishAccumulate(gsMatrix<T> & result) const
    {
        m_connHandler->finishAccumulate(result);
    }



   typename gsConnectionHandler<T>::Ptr m_connHandler;
   typename gsDistributedOperator<T>::Ptr m_op;
   typename gsDistributedOperator<T>::Ptr m_prec;
   gsOptionList m_options;
   typename IterativeSolver::Ptr m_solver;
};

#endif
}

