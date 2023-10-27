/** @file uwbINSSolverSteady.h

Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once

#include "uwbINSSolverBase.h"
#include "uwbINSAssemblerSteady.h"
#include "uwbINSAssemblerSteady_AFC.h"
//#include "uwbINSSUPGAssemblerSteady.h"
#include "uwbINSAssemblerSteadyPeriodic.h"
#include <gsTrilinos/gsTrilinos.h>

namespace gismo
{

template<class T>
class uwbINSSolverSteady : public uwbINSSolverBase<T>
{

public:
    typedef uwbINSSolverBase<T> Base;

public:
    uwbINSSolverSteady(uwbINSSolverParams<T>& params)
    {
        //create assembler
        if (params.settings().get(constantsINS::AFC) || params.settings().get(constantsINS::AFC_HO))
        {
            params.getAssemblerOptions().dirStrategy = dirichlet::none;
            m_pAssembler = new uwbINSAssemblerSteady_AFC<T>(params);
        }
        else if (params.getBCs().numPeriodic())
            m_pAssembler = new uwbINSAssemblerSteadyPeriodic<T>(params);
        else
            m_pAssembler = new uwbINSAssemblerSteady<T>(params);

        Base::initMembers();

        m_alpha_u = params.settings().get(constantsINS::alpha_u);
        m_alpha_p = params.settings().get(constantsINS::alpha_p);
    }

    virtual ~uwbINSSolverSteady()
    {
    }

public:
    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        this->updateAssembler();

        // in the first iteration the pattern will get analyzed
        if (!m_iterationNumber)
            this->initIteration();

        this->applySolver(m_solution);
        //this->applySolver(m_solution, m_alpha_u, m_alpha_p);

        m_iterationNumber++;
    }

    virtual T residualRelNorm() const
    {
        gsMatrix<T> residual = getAssembler()->rhs() - getAssembler()->matrix() * m_solution;
        T resNorm = residual.norm() / getAssembler()->rhs().norm();

        return resNorm;
    }

    virtual gsMatrix<T> getSolution_full(const gsMatrix<T>& solVector)
    {
        return getAssembler()->getSolution_full(solVector);
    }


    virtual uwbINSAssemblerSteady<T>* getAssembler() const
    {
        uwbINSAssemblerSteadyPeriodic<T>* pAssembler = dynamic_cast<uwbINSAssemblerSteadyPeriodic<T>*>(m_pAssembler);
        uwbINSAssemblerSteady_AFC<T>* pAssemblerAFC = dynamic_cast<uwbINSAssemblerSteady_AFC<T>*>(m_pAssembler);

        if (pAssembler != NULL)
            return pAssembler;
        else if (pAssemblerAFC != NULL)
            return dynamic_cast<uwbINSAssemblerSteady_AFC<T>*>(m_pAssembler);
        else
            return dynamic_cast<uwbINSAssemblerSteady<T>*>(m_pAssembler);
    }

    void evalElWiseForLocRef(const gsMatrix<T> & solVector, std::vector<gsVector<T> > & elWiseVals, bool outputInQuadPoints = false)
    {
        getAssembler()->evalElWiseForLocRef(solVector, elWiseVals, outputInQuadPoints);
    }

    void evalElWiseForLocRef(std::vector<gsVector<T> > & elWiseVals, bool outputInQuadPoints = false)
    {
        getAssembler()->evalElWiseForLocRef(m_solution, elWiseVals, outputInQuadPoints);
    }

    void plotResiduum(std::string const & fn, gsMatrix<T>& solution, unsigned npts = 10000)
    {
        getAssembler()->plotResiduum(fn, solution, npts);
    }
    void plotResiduum(std::string const & fn, unsigned npts = 10000)
    {
        getAssembler()->plotResiduum(fn, m_solution, npts);
    }


protected:
    real_t m_alpha_u, m_alpha_p;

    // members from uwbINSSolverBase
    using Base::m_pAssembler;
    using Base::m_solution;
    using Base::m_iterationNumber;

}; //uwbINSSolverSteady

} //namespace gismo
