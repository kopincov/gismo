/** @file gsNewtonIETI.h

    @brief A version of the Newton iterator, where as linear solver the
    IETI method is used.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C.Hofer
    Created on: 2016-06-06
*/


#pragma once
#include <gsPde/gsNewtonIterator.h>
#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIScaledDirichlet.h>
#include <gsSolver/gsConjugateGradient.h>

namespace gismo {


template<class T>
class gsNewtonIETI : public gsNewtonIterator<T>
{
    typedef gsNewtonIterator<T> Base;

public:
    gsNewtonIETI( gsAssembler<T> & assembler, const gsMultiPatch<T> & initialSolution,
                 gsOptionList opt = gsIETIAssembler<T>::defaultOptions()):
        gsNewtonIterator<T>(assembler,initialSolution)
      , m_assembler(assembler )
      ,m_PCGTol(Base::m_tolerance*1.e-1)
    {
        m_assembler.setOptions(opt);
        m_assembler.init();
    }

    gsNewtonIETI( gsAssembler<T> & assembler, gsOptionList opt = gsIETIAssembler<T>::defaultOptions()):
        gsNewtonIterator<T>(assembler)
      , m_assembler(assembler)
      ,m_PCGTol(Base::m_tolerance*1.e-1)
    {
        m_assembler.setOptions(opt);
        m_assembler.init();
    }

    void setPCGTol(T tol) {m_PCGTol=tol;}

protected:

    virtual void solveLinearProblem(gsMatrix<T> &updateVector);

    virtual void solveLinearProblem(const gsMultiPatch<T> & currentSol, gsMatrix<T> &updateVector);

    virtual T getResidue() {return m_res;}
protected:
    gsIETIAssembler<T> m_assembler;
    IETIPrecondScaling::strategy m_scaling;

    gsMatrix<T> m_lagrangeMult;

    T m_res;

    T m_PCGTol;

};


template <class T>
void gsNewtonIETI<T>::solveLinearProblem(gsMatrix<T>& updateVector)
{
    m_assembler.getOptions().opt.setSwitch("NonlinearMode",false);
    m_assembler.getOptions().opt.setSwitch("CalcRhsNorm",true);

    m_assembler.assemble();

    typename gsIETISolver<T>::Ptr solv = memory::make_shared(new gsIETISolver<T>(m_assembler));
    solv->init();
    typename gsScaledDirichletPrecond<T>::Ptr prec =
        gsScaledDirichletPrecond<T>::make(m_assembler);
    gsConjugateGradient<> PCG(solv,prec);
    PCG.setMaxIterations(100);
    PCG.setTolerance(m_PCGTol);
    PCG.setCalcEigenvalues(true);
    m_lagrangeMult.setZero(   m_assembler.systemSize() , m_assembler.numberRhs());

    PCG.solve(solv->getRhs(), m_lagrangeMult);

    solv->calculateSolution(m_lagrangeMult,updateVector);

    m_res = m_assembler.getRhsNorm();

 //   gsDebugVar(PCG.iterations());

 //   gsDebugVar(updateVector.transpose());

    gsMatrix<T> eigs;
    PCG.getEigenvalues(eigs);
  //  gsDebugVar(eigs.maxCoeff()/eigs.minCoeff());
    gsInfo<<"IETI-Iterations: "<<PCG.iterations()<<" ; Estimated ConditionNumber: "<<eigs.maxCoeff()/eigs.minCoeff()<<"\n";
}

template <class T>
void gsNewtonIETI<T>::solveLinearProblem(const gsMultiPatch<T> & currentSol, gsMatrix<T>& updateVector)
{
    m_assembler.getOptions().opt.setSwitch("NonlinearMode",true);
    m_assembler.getOptions().opt.setSwitch("CalcRhsNorm",true);


    m_assembler.assemble( currentSol);
    typename gsIETISolver<T>::Ptr solv = memory::make_shared(new gsIETISolver<T>(m_assembler));
    solv->init();
    typename gsScaledDirichletPrecond<T>::Ptr prec =
        gsScaledDirichletPrecond<T>::make(m_assembler);
    gsConjugateGradient<> PCG(solv,prec);
    PCG.setMaxIterations(100);
    PCG.setTolerance(m_PCGTol);
    PCG.setCalcEigenvalues(true);
   // m_lagrangeMult.setZero( m_assembler.systemSize() , m_assembler.numberRhs());

    PCG.solve(solv->getRhs(),m_lagrangeMult);

    solv->calculateSolution(m_lagrangeMult,updateVector);

    m_res = m_assembler.getRhsNorm();

 //   gsDebugVar(updateVector.transpose());

   // gsDebugVar(PCG.iterations());

    gsMatrix<T> eigs;
    PCG.getEigenvalues(eigs);
  //  gsDebugVar(eigs.maxCoeff()/eigs.minCoeff());
    gsNewtonIterator<T>::m_assembler.homogenizeFixedDofs(-1);

    gsInfo<<"IETI-Iterations: "<<PCG.iterations()<<" ; Estimated ConditionNumber: "<<eigs.maxCoeff()/eigs.minCoeff()<<"\n";
}







} // namespace gismo

