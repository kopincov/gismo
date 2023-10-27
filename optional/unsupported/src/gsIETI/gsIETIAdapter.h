/** @file gsIETIAdapter.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    Created on: 2018-06-05
*/


#pragma once


#include <gsCore/gsConfig.h>
#include <gsSolver/gsLinearOperator.h>
#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIScaledDirichlet.h>


namespace gismo {

template<typename T>
class IETIAdapter : public gsLinearOperator<T>
{

public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<IETIAdapter> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<IETIAdapter> uPtr;

    IETIAdapter(typename gsIETIAssembler<T>::Ptr IETIAss, unsigned maxIter, real_t tol,bool prec,bool GMRES=false,real_t lambda = -1);

    static uPtr make(typename gsIETIAssembler<T>::Ptr IETIAss, unsigned maxIter, real_t tol,bool prec,bool GMRES=false, real_t lambda= -1 ) {return memory::make_unique(new IETIAdapter(IETIAss,maxIter,tol,prec,GMRES,lambda));}

    virtual ~IETIAdapter() {}

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

    /// Returns the number of rows of the operator
    virtual index_t rows() const {return m_IETIAss->getInfo().origSystemSize;}

    /// Returns the number of columns of the operator
    virtual index_t cols() const {return rows();}

    real_t getMaxEigenvalue() const {return m_lambda;}

private:

    typename gsIETIAssembler<T>::Ptr m_IETIAss;
    mutable typename gsIETISolver<T>::Ptr m_IETISolv;
    mutable gsMatrix<T> m_solLambda;
    real_t m_lambda;
    mutable typename gsIterativeSolver<T>::uPtr m_iterationMethod;
};

}

#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>

#include <gsIETI/gsDistributedOperator.h>
#include <gsIETI/gsParallelOperator.h>
#include <gsIETI/gsIETIAssemblerMPI.h>
#include <gsIETI/gsIETISolverMPI.h>
#include <gsIETI/gsIETIScaledDirichletMPI.h>
#include <gsSolver/gsPowerIteration.h>

namespace gismo
{

template<typename T>
class IETIAdapterMPI : public gsDistributedOperator<T>
{

public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<IETIAdapterMPI> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<IETIAdapterMPI> uPtr;

    IETIAdapterMPI(typename gsConnectionHandler<T>::Ptr A, typename gsIETIAssemblerMPI<T>::Ptr IETIAss, gsParallelGlobalLocalHandler::Ptr parHandler, unsigned maxIter, real_t tol, bool prec, bool GMRES=false, real_t lambda = -1, bool ouptput=false);

    static uPtr make(typename gsConnectionHandler<T>::Ptr A,typename gsIETIAssemblerMPI<T>::Ptr IETIAss, gsParallelGlobalLocalHandler::Ptr parHandler,unsigned maxIter, real_t tol,bool prec,bool GMRES=false, real_t lambda= -1, bool output =false ) {return memory::make_unique(new IETIAdapterMPI(A,IETIAss,parHandler,maxIter,tol,prec,GMRES,lambda,output));}

    virtual ~IETIAdapterMPI() {}

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

    /// Returns the number of rows of the operator
    virtual index_t rows() const {return m_parHandler->localSize();}

    /// Returns the number of columns of the operator
    virtual index_t cols() const {return rows();}

    real_t getMaxEigenvalue() const {return m_lambda;}

    virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const;

    virtual void postAccumulate() const;

    virtual void startAccumulate(const  gsMatrix<T> & input) const;

    virtual void finishAccumulate(gsMatrix<T> & result) const;


private:
    typename gsConnectionHandler<T>::Ptr  m_A;
    typename gsIETIAssemblerMPI<T>::Ptr m_IETIAss;
    mutable typename gsIETISolverMPI<T>::Ptr m_IETISolv;
    mutable gsMatrix<T> m_solLambda;
    mutable gsMatrix<T> m_globVec;
    real_t m_lambda;
    mutable typename gsIterativeSolver<T>::uPtr m_iterationMethod;
    gsParallelGlobalLocalHandler::Ptr m_parHandler;
    bool m_isPreconditioner;
    bool m_output;
};

} // namespace gismo


#endif

#ifndef GISMO_BUILD_LIB
#include <gsIETI/gsIETIAdapter.hpp>
#endif
