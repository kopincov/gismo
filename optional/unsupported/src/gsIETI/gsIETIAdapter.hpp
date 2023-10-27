/**  gsIETIAdapter.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    Created on:  2018-06-05
*/


#include  <gsIETI/gsIETIAdapter.h>

#include <gsSolver/gsPowerIteration.h>
#include <gsSolver/gsConjugateGradient.h>
#include <gsSolver/gsGMRes.h>
#include <gsSolver/gsGradientMethod.h>
#include <gsSolver/gsPowerIteration.h>

namespace gismo {

template<typename T>
IETIAdapter<T>::IETIAdapter(typename gsIETIAssembler<T>::Ptr IETIAss, unsigned maxIter, real_t tol,bool prec,bool GMRES,real_t lambda)
    :m_IETIAss(IETIAss) , m_IETISolv(gsIETISolver<T>::make(*m_IETIAss))
{
    typename gsScaledDirichletPrecond<T>::Ptr precond = gsScaledDirichletPrecond<T>::make(*m_IETIAss);
    if(prec && false)
    {
        if( lambda==-1)
        {
            //typename IETIApplication::Ptr IETIApp = IETIApplication::make(IETIAss);

            m_lambda = powerIteration<T>(m_IETISolv,precond);
            gsInfo<<"maxEig: "<<m_lambda<<"\n";
        }
        else
            m_lambda = lambda;
        //gsChebyshevSemiIteration::uPtr uPtr = memory::make_unique<gsChebyshevSemiIteration>(new gsChebyshevSemiIteration(m_IETISolv,1,m_lambda,gsScaledDirichletPrecond<T>::make(*m_IETIAss)));
        m_iterationMethod = gsGradientMethod<T>::make(m_IETISolv,precond,2/(1+m_lambda));
    }
    else
    {
        if(GMRES)
            m_iterationMethod = gsGMRes<T>::make(m_IETISolv,precond);
        else
            m_iterationMethod = gsConjugateGradient<T>::make(m_IETISolv,precond);

    }

    m_iterationMethod->setMaxIterations(maxIter);
    m_iterationMethod->setTolerance(tol);
    m_solLambda.setZero(m_IETISolv->rows(),1);
}

template<typename T>
void IETIAdapter<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    m_solLambda.setZero(m_solLambda.rows(),input.cols());
    m_IETIAss->setNewRhs(input);
    m_IETISolv->init();

    /*
        gsMatrix<T> mat, col;
        mat.setZero(m_IETISolv->cols(),m_IETISolv->cols());
        for(int c=0; c<m_IETISolv->cols();++c)
        {
            m_IETISolv->apply(gsMatrix<T>::Identity(m_IETISolv->cols(),m_IETISolv->cols()).col(c),col);
            mat.col(c)=col;
        }
        gsInfo<<"inp:\n "<<input<<"\n\n F: \n" <<mat<<"\n\n rhs: "<<m_IETISolv->getRhs()<<"\n\n";
        */

    m_iterationMethod->solve(m_IETISolv->getRhs(),m_solLambda);
    gsInfo<<"\t \t IETI: it "<<m_iterationMethod->iterations()<<"  --  "<<"res "<<m_iterationMethod->error()<<"\n"<<std::flush;
    m_IETISolv->calculateSolution(m_solLambda,x);
    /*
        gsInfo<<" IETI Input: \n "<<input.transpose()<<"\n";
        gsInfo<<"x.size(): "<<x.rows() <<" | Op.rows(): "<<rows()<<" \n "<<x.transpose()<<"\n";
        gsInfo<<"m_solLambda: \n" <<m_solLambda<<"\n\n";
        */
}

}

#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>

#include <gsIETI/gsDistributedOperator.h>
#include <gsIETI/gsParallelOperator.h>
#include <gsIETI/gsIETIAssemblerMPI.h>
#include <gsIETI/gsIETISolverMPI.h>
#include <gsIETI/gsIETIScaledDirichletMPI.h>
#include <gsIETI/gsParallelGradientMethod.h>
#include <gsSolver/gsPowerIteration.h>

#include <gsIETI/gsParallelCG.h>
#include <gsIETI/gsParallelGMRes.h>
#include <gsIETI/gsParallelGradientMethod.h>

namespace gismo
{

template<typename T>
IETIAdapterMPI<T>::IETIAdapterMPI(typename gsConnectionHandler<T>::Ptr A, typename gsIETIAssemblerMPI<T>::Ptr IETIAss, gsParallelGlobalLocalHandler::Ptr parHandler, unsigned maxIter, real_t tol,bool prec,bool GMRES,real_t lambda, bool ouptput)
    : m_A(A), m_IETIAss(IETIAss) , m_IETISolv(gsIETISolverMPI<T>::make(*m_IETIAss)), m_parHandler(parHandler), m_output(ouptput)
{
    // gsInfo<<"rank: "<<gsMpi::worldRank()<<"setup IETI-Adapter\n";
    typename gsScaledDirichletPrecondMPI<T>::Ptr precond = gsScaledDirichletPrecondMPI<T>::make(*m_IETIAss);
    if(prec && false)
    {
        if( lambda==-1 )
        {
            //typename IETIApplicationMPI::Ptr IETIApp = IETIApplicationMPI::make(IETIAss);
            m_lambda = powerIteration_MPI<T>(m_IETISolv,precond,m_parHandler->getComm(),15);
            //m_lambda = parHandler->getComm().sum(m_lambda)/parHandler->getComm().size();
            if(gsMpi::worldRank() == 0) gsInfo<<"rank: "<<gsMpi::worldRank()<<"maxEig: "<<m_lambda<<"\n";
        }
        else
            m_lambda = lambda;
        //gsParallelGradientMethod::uPtr uPtr = gsParallelGradientMethod::uPtr(new gsParallelGradientMethod(m_IETISolv,gsScaledDirichletPrecondMPI<T>::make(*m_IETIAss),2/(1+m_lambda),m_parHandler->getComm()));
        //gsChebyshevSemiIteration::uPtr uPtr = memory::make_unique<gsChebyshevSemiIteration>(new gsChebyshevSemiIteration(m_IETISolv,1,m_lambda,gsScaledDirichletPrecond<T>::make(*m_IETIAss)));
        m_iterationMethod = gsParallelGradientMethod<T>::make(m_IETISolv,precond,1/(m_lambda*m_lambda),m_parHandler->getComm());
    }
    else
    {
        if(GMRES)
            m_iterationMethod = gsParallelGMRes<T>::make(m_IETISolv,precond,m_parHandler->getComm());
        else
            m_iterationMethod = gsParallelCG<T>::make(m_IETISolv,precond,m_parHandler->getComm());

    }
    m_isPreconditioner = prec;

    m_iterationMethod->setMaxIterations(maxIter);
    m_iterationMethod->setTolerance(tol);
    m_solLambda.setZero(m_IETISolv->rows(),1);
}

template<typename T>
void IETIAdapterMPI<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    //if(gsMpi::worldRank()==0) gsInfo<<"rank: "<<gsMpi::worldRank()<<"before solve A_inv:\n"<<input.transpose()<<"\n\n";

    m_solLambda.setZero(m_solLambda.rows(),input.cols());
    m_globVec.setZero(m_parHandler->globalSize(),input.cols());

    //  gsInfo<<"input before adding: "<<input.transpose()<<"\n\n";
    //m_parHandler->addLocalVectorToGlobal(input,m_globVec);

    //  glob.setZero(m_globVec.rows(),m_globVec.cols());
    //  m_parHandler->buildGlobalVector(input,glob);
    // if(gsMpi::worldRank()==0) gsInfo<<"rank: "<<gsMpi::worldRank()<<"before solve A_inv:\n"<<glob.transpose()<<"\n\n";
    gsMatrix<T> inpD = input;
    m_A->accumulateDistributedVector(inpD);
    m_A->distributeAccumulatedVector(inpD);
    //glob.setZero(m_globVec.rows(),m_globVec.cols());
    m_parHandler->addLocalVectorToGlobal(inpD,m_globVec);
    //   gsInfo<<"input after adding: "<<m_globVec.transpose()<<"\n\n";
    //  m_IETIAss->setNewRhs(m_globVec);
    m_IETIAss->setNewRhs(m_globVec);
    m_IETISolv->init();


    /*
        gsMatrix<T> mat, col;
        mat.setZero(m_IETISolv->cols(),m_IETISolv->cols());
        for(int c=0; c<m_IETISolv->cols();++c)
        {
            m_IETISolv->apply(gsMatrix<T>::Identity(m_IETISolv->cols(),m_IETISolv->cols()).col(c),col);
            mat.col(c)=col;
        }

        for(int pp = 0; pp< m_IETIAss->m_comm.size();++pp)
        {
            m_IETIAss->m_comm.barrier();
            if(pp==m_IETIAss->m_comm.rank())
            {
                gsInfo<<"inp:\n "<<input.transpose()<<"\n\n F: \n" <<mat<<"\n\n rhs: "<<m_IETISolv->getRhs().transpose()<<"\n\n";
            }
        }
*/
    m_iterationMethod->solve(m_IETISolv->getRhs(),m_solLambda);
    //  if(gsMpi::worldRank()==0 && iterCount>=7550 && iterCount<=7560)      gsInfo<<"lambda: "<<m_solLambda.transpose()<<"\n";
    if(m_output && gsMpi::worldRank()==0) gsInfo<<"rank: "<<gsMpi::worldRank()<<"\t \t IETI: it "<<m_iterationMethod->iterations()<<"  --  "<<"res "<<m_iterationMethod->error()<<"\n"<<std::flush;


    m_IETISolv->calculateSolution(m_solLambda,m_globVec);
    //   m_parHandler->extractLocalVector(m_globVec,x);
    //  m_IETIAss->combineToCommonSolution(m_globVec);
    //  if(gsMpi::worldRank()==0) gsInfo<<"globVec: \n"<<m_globVec.transpose()<<"\n";


    m_parHandler->extractLocalVector(m_globVec,x);
    //  if(gsMpi::worldRank()==0) gsInfo<<"globVec: \n"<<x.transpose()<<"\n";
    // if(gsMpi::worldRank()==0) gsInfo<<"rank: "<<gsMpi::worldRank()<<"x.size(): "<<m_globVec.rows() <<" | Op.rows(): "<<rows()<<" \n "<<globVec.transpose()<<"\n\n";


    /*
        for(int pp = 0; pp< m_IETIAss->m_comm.size();++pp)
        {
            m_IETIAss->m_comm.barrier();
            if(pp==m_IETIAss->m_comm.rank())
            {
                gsInfo<<"rank: "<<gsMpi::worldRank()<<" IETI Input: \n "<<input.transpose()<<"\n";
                gsInfo<<"rank: "<<gsMpi::worldRank()<<"x.size(): "<<x.rows() <<" | Op.rows(): "<<rows()<<" \n "<<x.transpose()<<"\n";

                gsInfo<<"m_solLambda: \n" <<m_solLambda.transpose()<<"\n\n";
            }
        }
        */
}

template<typename T>
void IETIAdapterMPI<T>::distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
{
    distributed = input;
    m_A->distributeAccumulatedVector(distributed);
}

template<typename T>
void IETIAdapterMPI<T>::postAccumulate() const
{
    m_A->postAccumulate();
}

template<typename T>
void IETIAdapterMPI<T>::startAccumulate(const  gsMatrix<T> & input) const
{
    m_A->startAccumulate(input);
}

template<typename T>
void IETIAdapterMPI<T>::finishAccumulate(gsMatrix<T> & result) const
{
    m_A->finishAccumulate(result);
}

} // namespace gismo
#endif
