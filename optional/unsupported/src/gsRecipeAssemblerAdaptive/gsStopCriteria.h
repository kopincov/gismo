/** @file gsStopCriteria.h

    @brief DESCRIPTION

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/



#pragma once

#include <gsRecipeAssemblerAdaptive/gsAdaptiveSolver.h>

namespace gismo {




/// @brief Abstract class for stop criteria for iterative solvers.
///
/// \ingroup Solver
class gsStopCriteria
{
protected:
    int  m_curStep;
    bool m_succeedFlag;
public:
    gsStopCriteria() : m_curStep(0),m_succeedFlag(false) { }

    virtual ~gsStopCriteria(){}

    bool stop()
    {
        m_curStep++;
        return this->stopImpl(m_curStep);
    }

    int getCurStepNumber()
    {
        return m_curStep;
    }

    virtual void reset()
    {
        m_curStep=0;
        m_succeedFlag=false;
    }

    virtual bool succeed() const
    { return m_succeedFlag; }

protected:
    virtual bool stopImpl(int iterationNr) = 0;
};

class gsStopCriteriaIterationOnly : public gsStopCriteria
{
    int  m_maxIter;
public:
    gsStopCriteriaIterationOnly(int maxSteps) : gsStopCriteria(), m_maxIter(maxSteps) { }

    ~gsStopCriteriaIterationOnly(){}
protected:
    bool stopImpl(int iterationNr)
    {
        if(iterationNr==m_maxIter)
            m_succeedFlag=true;
        return iterationNr >= m_maxIter ? true : false;
    }
};

class gsStopCriteriaIterationAndTotalError : public gsStopCriteria
{
protected:
    real_t   m_target;
    int      m_maxIter;

    gsErrorEstimator *m_estimator;

public:
    gsStopCriteriaIterationAndTotalError(gsErrorEstimator *estimator, unsigned maxIt=10, real_t target=NAN )
        : m_estimator(estimator)
    {
        m_target =  math::isnan(target) ? 1e-10 :  target;
        m_maxIter   =  maxIt;
    }
protected:
    virtual bool stopImpl(int iterationNr)
    {
        real_t tot=math::sqrt(m_estimator->getTotalErrorEstimate().sum());
        if(math::sqrt(tot)<m_target)
        {
            m_succeedFlag=true;
            return true;
        }
        else
            return iterationNr >= m_maxIter;
    }
};




}
