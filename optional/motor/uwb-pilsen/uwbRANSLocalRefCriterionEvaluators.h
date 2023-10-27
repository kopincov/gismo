/** @file uwbRANSLocalRefCriterionEvaluators.h

    Author(s): E. Turnerova, K. Slaba
*/

#pragma once

#include "uwbLocalRefCriterionEvaluators.h"
#include "uwbTMEvaluators.h"
#include "uwbTMSolverBase.h"
//#include "uwbStabilizationEvaluator.h"

namespace gismo
{
template <class T>
class uwbLocRefRANSResiduumEvaluator : public uwbLocRefEvaluator<T>
{

public:
    typedef uwbLocRefEvaluator<T> Base;

public:

    uwbLocRefRANSResiduumEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_bIsUnsteady(false), m_bOldSolutionSet(false)
    {
        m_pTMsolver = NULL;
        m_sTMEvaluatorType = "";
        m_elemValue.setZero(1);
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity and pressure solution set in the RANS locRef evaluator.");
        if (m_bIsUnsteady)
            GISMO_ASSERT(m_bOldSolutionSet, "No old time step velocity and pressure solution set in the RANS locRef evaluator.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        // Evaluate turbulent viscosity at element nodes
        evalTurbViscosity(m_solUGrads, quNodes, geoEval);

        const gsMatrix<T> & basisGradsU = m_basisDataU[1];
        const gsMatrix<T> & bHessian_u = m_basisDataU[2]; //-> obsahuje i smisenou druhou derivaci

        gsMatrix<T> physLaplacian;

        m_solUValsOld = m_solUVals;
        if (m_bIsUnsteady)
            m_solUValsOld = m_oldSolU.value(quNodes, m_patchIndex);

        m_residuumMomentumEq.resize(m_dim);
        index_t nQuPoints = quNodes.cols();
        for (int var = 0; var < m_dim; var++)
            m_residuumMomentumEq[var].setZero(nQuPoints);
        m_residuumContinuityEq.setZero(nQuPoints);

        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
            gsMatrix<T> solLaplacian = m_solActUCoeffs * physLaplacian.transpose();

            for (int var = 0; var < m_dim; var++)
            {
                m_residuumMomentumEq[var](k) = m_solUGrads[k].row(var) * m_solUVals.col(k) //advection
                                              - (m_viscosity + getTurbViscosityVal(k)) * solLaplacian(var, 0) //diffusion
                                              + m_solPGrads[k](0, var); //pressure term
                                              //- f; //source term ... -2/3 grad(k)

                if (m_bIsUnsteady)
                    m_residuumMomentumEq[var](k) += 1./m_timeStep * (m_solUVals(var, k) - m_solUValsOld(var, k));

                m_residuumContinuityEq(k) += m_solUGrads[k](var, var);
            }
        }

        /* part of the residuum
         * for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            for (int var = 0; var < m_dim; var++)
            {
                m_residuumMomentumEq[var](k) = (m_solUGrads[k].row(var)).norm();
                m_residuumContinuityEq(k) = (m_solPGrads[k]).norm();
            }
        }*/
    }

    void evalTurbViscosity(std::vector<gsMatrix<T> >& solUGrads, gsMatrix<T>& quNodes, gsGeometryEvaluator<T> & geoEval)
    {
        GISMO_ENSURE(m_pTMsolver != NULL, "uwbRANSResiduumEvaluator: No turbulent model set!");
        m_sTMEvaluatorType = m_pTMsolver->getTMEvaluator();
        GISMO_ASSERT(m_sTMEvaluatorType != "", "No evaluator type set.");

        m_numTMvar = m_pTMsolver->getAssembler()->getNumVar();
        //if (m_sTMEvaluatorType == "koWilcoxLRN" || m_sTMEvaluatorType == "koSST" || m_sTMEvaluatorType == "koSSTMenter2009")
        //    numVar = 2;

        typename uwbTMEvaluator<T>::uPtr evaluatorTM = uwbTMEvaluator<T>::make(m_sTMEvaluatorType);

        evaluatorTM->initialize(m_viscosity, quNodes.cols());
        evaluatorTM->setKOmegaVariant(m_sTMEvaluatorType);

        //solution from the last time step is stored in solver
        gsField<T> solTMfield = m_pTMsolver->constructSolution();
        gsMatrix<T> solKOVals = solTMfield.value(quNodes, m_patchIndex);

        //--- KOGrads
        gsMatrix<T> solActKOCoeffs, physGradKO;//, bGradsKO;
        gsMatrix<index_t> activesKO;
        std::vector<gsMatrix<T> > basisDataKO;

        gsMultiBasis<T> basisKO = m_pTMsolver->getAssembler()->getBlockAssembler().getSolBasis();
        basisKO.back().active_into(quNodes.col(0), activesKO);
        basisKO.back().evalAllDers_into(quNodes, 2, basisDataKO);
        //basisKO.back().deriv_into(quNodes, bGradsKO);

        const index_t numActKO = activesKO.rows();

        solActKOCoeffs.setZero(m_numTMvar, numActKO);
        for (int j = 0; j < numActKO; j++)
            solActKOCoeffs.col(j) = solTMfield.coefficientVector(m_patchIndex).row(activesKO(j)).transpose();

        index_t nQuPoints = quNodes.cols();
        std::vector<gsMatrix<T> > solKOGrads(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, basisDataKO[1], physGradKO);
            solKOGrads[k].noalias() = solActKOCoeffs * physGradKO.transpose();
        }

        if (this->checkWallDistanceBasedTM())
        {
            gsField<T> solPoisson = m_pTMsolver->getAssembler()->getBlockAssembler().getPoissonSolution();
            gsMultiPatch<T> patchesPoisson = m_pTMsolver->getAssembler()->getBlockAssembler().getPoissonPatches();
            gsMultiBasis<T> basesPoisson = m_pTMsolver->getAssembler()->getBlockAssembler().getPoissonBasis();

            unsigned evFlagsPoisson = NEED_MEASURE | NEED_GRAD_TRANSFORM;
            const gsGeometry<T>& geo = patchesPoisson.patch(m_patchIndex);
            typename gsGeometryEvaluator<T>::uPtr geoEvalPoisson(getEvaluator(evFlagsPoisson, geo));

            gsMatrix<index_t> activesPoisson;
            gsMatrix<T> basisGradsPoisson, physGradPoisson;
            basesPoisson.basis(m_patchIndex).active_into(quNodes.col(0), activesPoisson);
            const index_t numActPoisson = activesPoisson.rows();
            basesPoisson.basis(m_patchIndex).deriv_into(quNodes, basisGradsPoisson);

            geoEvalPoisson->evaluateAt(quNodes);

            gsMatrix<T> solPoissonVals = solPoisson.value(quNodes, m_patchIndex);
            gsMatrix<T> solActPoissonCoeffs;
            solActPoissonCoeffs.setZero(1, numActPoisson);
            for (int j = 0; j < numActPoisson; j++)
                solActPoissonCoeffs.col(j) = solPoisson.coefficientVector(m_patchIndex).row(activesPoisson(j)).transpose();

            std::vector<gsMatrix<T> > solPoissonGrads(nQuPoints);
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEvalPoisson->transformGradients(k, basisGradsPoisson, physGradPoisson);
                solPoissonGrads[k].noalias() = solActPoissonCoeffs * physGradPoisson.transpose();
            }

            evaluatorTM->evalWallDistance(solPoissonVals, solPoissonGrads);
        }

        evaluatorTM->initAtElement(solUGrads, solKOVals, solKOGrads);
        evaluatorTM->evalQuantities_turbViscosity();
        m_turbViscosityVals = evaluatorTM->getTurbViscosity();
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t nQuPoints = quWeights.rows();

        T elResMomentum = 0.;
        T elResContinuity = 0.;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            for (index_t var = 0; var != m_dim; ++var)
                elResMomentum += weight * math::pow(m_residuumMomentumEq[var](k), 2);

            elResContinuity += weight * math::pow(m_residuumContinuityEq(k), 2);
        }
        m_elemValue(0) = (math::sqrt(elResMomentum) + math::sqrt(elResContinuity)) * m_elementLength;
    }

    void setTurbulenceSolver(uwbTMSolverBase<T> * pTMsolver) { m_pTMsolver = pTMsolver; }

    gsMatrix<T> getTurbViscosityVals()
    {
        GISMO_ASSERT(m_turbViscosityVals.rows() > 0, "Turbulent viscosity not evaluated yet.");
        gsMatrix<T> turbViscosity = m_turbViscosityVals;
        return turbViscosity;
    }

    T getTurbViscosityVal(int k)
    {
        GISMO_ASSERT(m_turbViscosityVals.rows() > 0, "Turbulent viscosity not evaluated yet.");
        return m_turbViscosityVals(k);
    }

    bool checkWallDistanceBasedTM()
    {
        return (m_sTMEvaluatorType == "koSST" || m_sTMEvaluatorType == "koSSTMenter2009");
        // || m_sTMEvaluatorType == "koSAS" || m_sTMEvaluatorType == "koSAS_SS" || m_sTMEvaluatorType == "koSAS_SO" || m_sTMEvaluatorType == "koSAS_OO");
    }

    void setOldSolutionField(bool unsteady, gsField<T>& solution, T timeStep)
    {
        m_oldSolU = solution;
        m_bIsUnsteady = unsteady;
        m_timeStep = timeStep;
        m_bOldSolutionSet = true;
    }

protected:
    bool m_bIsUnsteady;
    bool m_bOldSolutionSet;

    T m_timeStep;
    gsField<T> m_oldSolU;
    gsMatrix<T> m_solUValsOld;

    std::vector<gsVector<T> > m_residuumMomentumEq;
    gsVector<T> m_residuumContinuityEq;
    gsVector<T> m_turbViscosityVals;
    uwbTMSolverBase<T> * m_pTMsolver;
    std::string m_sTMEvaluatorType;
    int m_numTMvar;

protected:
    using Base::m_elemValue;
    using Base::m_viscosity;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_solActUCoeffs;
    using Base::m_solUGrads;
    using Base::m_solPGrads;
    using Base::m_basisDataU;
    using Base::m_solUVals;
    using Base::m_patchIndex;

    using uwbLocRefEvaluatorBase<T>::m_elementLength;
};

} // namespace gismo

