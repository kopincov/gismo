/** @file gsAdaptiveSolver.h

    @brief abstract implementation of an adaptive method based on the
           ( assemble, solve, estimate, mark, refine) loop

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsRecipeAssemblerAdaptive/gsErrorEstimator.h>
#include <gsRecipeAssemblerAdaptive/gsMarker.h>
#include <gsRecipeAssemblerAdaptive/gsSpaceRefiner.h>
#include <gsRecipeAssemblerAdaptive/gsStopCriteria.h>

#include <numeric>

namespace gismo {




class gsAdaptiveSolver
{
protected:
    gsRecipeAssembler                     *m_assembler;
    gsSparseSolver<real_t>                *m_solver;
    gsSparseSolver<real_t>                *m_solverEli;
    gsErrorEstimator                      *m_estimator;
    gsMarker                              *m_marker;
    gsSpaceRefiner                        *m_refiner;

    gsStopCriteria                        *m_criteria;
protected:
    gsMatrix<real_t>                       m_modifiedRhs;
    gsMatrix<real_t>                       m_sol;
    gsMatrix<real_t>                       m_eliSol;

    std::vector<gsMatrix<real_t> >         m_solCoefs;

    bool                                   m_internalErrorFlag;
protected:
    gsAdaptiveSolver()
        : m_assembler(NULL), m_solver(NULL), m_solverEli(NULL), m_estimator(NULL), m_marker(NULL), m_refiner(NULL), m_criteria(NULL),
          m_internalErrorFlag(false)
    {}

public:
    gsAdaptiveSolver(gsRecipeAssembler &assembler, gsSparseSolver<real_t> &solverSys, gsSparseSolver<real_t> &solverEli,gsErrorEstimator &estimator,gsMarker &marker, gsSpaceRefiner &refiner, gsStopCriteria &crit)
        : m_assembler(&assembler), m_solver(&solverSys), m_solverEli(&solverEli), m_estimator(&estimator), m_marker(&marker), m_refiner(&refiner), m_criteria(&crit),
          m_internalErrorFlag(false)
    {}

    virtual ~gsAdaptiveSolver(){}

    virtual void logProgress()
    {
        std::cout << "iteration: " << m_criteria->getCurStepNumber()
                  << " totalError: " << math::sqrt(m_estimator->getTotalErrorEstimate().sum())
                  << " dofs: " << m_assembler->getSysSize()
                  << std::endl;
    }

    virtual void assemble()
    {
        m_assembler->setSpace(m_refiner->getSpaces());
        m_assembler->assemble();
    }

    virtual bool solveEliminatedDofs()
    {
        const gsMatrix<real_t>       &eliRhs = m_assembler->getEliminatedRhs();
        const gsSparseMatrix<real_t> &eliSys = m_assembler->getEliminatedMatrix();
        const gsSparseMatrix<real_t> &rhsMod = m_assembler->getRhsModMatrix();
        const gsMatrix<real_t>       &rhs    = m_assembler->getSystemRhs();

        if(eliSys.rows()>0&&eliSys.cols()>0)
        {
            m_solverEli->compute(eliSys);
            m_eliSol=m_solverEli->solve(eliRhs);
            if(!m_solverEli->succeed())
            {
                gsWarn<<"Failed to compute eliminated dofs.\n";
                m_internalErrorFlag=true;
                return false;
            }
            m_modifiedRhs=rhs-rhsMod*m_eliSol;
        }
        else
            m_modifiedRhs=rhs;
        return true;
    }

    virtual bool solve()
    {
        return this->solveEliminatedDofs() && this->solveSystem();
    }

    virtual bool solveSystem ()
    {
        const gsSparseMatrix<real_t> &sysMat = m_assembler->getSystemMatrix();

        m_solCoefs.clear();
        m_solver->compute(sysMat);
        m_sol=m_solver->solve(m_modifiedRhs);
        if(!m_solver->succeed())
        {
            gsWarn<<"Failed to compute solution.\n";
            m_internalErrorFlag=true;
            return false;
        }
        for (size_t i=0; i<m_assembler->getSpace().size();++i)
            m_solCoefs.push_back(m_assembler->reconstructSolution(i,m_sol,m_eliSol));
        return true;
    }

    virtual void errorEstimate ()
    {
        m_estimator->estimate(m_refiner->getSpaces(),m_solCoefs);
    }

    virtual bool succeed() const
    {
        return m_criteria->succeed() && !m_internalErrorFlag;
    }

    virtual void mark()
    {
        m_marker->mark(m_estimator->getLocalErrorEstimate());
    }

    virtual void refine()
    {
        m_refiner->updateSpaces(m_marker->getMarked());
    }

    virtual void initStep()
    {
        m_assembler->reset();
        if(m_marker)
            m_marker->reset();
        m_estimator->reset();
    }

    virtual void adaptiveSolve ()
    {
        m_criteria->reset();
        do
        {
            this->initStep();
            this->assemble();
            if(!this->solve())
            {
                gsWarn<<"Exiting...\n";
                m_internalErrorFlag=true;
                break;
            }
            this->errorEstimate();
            this->mark();
            logProgress();
            if (m_criteria->stop())
                break;
            this->refine();
        }
        while (true);
    }
};

}

