/** @file uwbTMLocalRefCriterionEvaluators.h

    Author(s): E. Turnerova
*/

#pragma once

#include "uwbLocalRefCriterionEvaluators.h"
#include "uwbTMEvaluators.h"

namespace gismo
{
// ============================================================= PARENT ============================================================= //
template <class T>
class uwbLocRefTMResiduumEvaluator : public uwbLocRefEvaluatorBase<T>
{
public:
    uwbLocRefTMResiduumEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        m_viscosity(viscosity), m_bOldSolSet(false), m_bIsUnsteady(false), m_numVar(2)//,
        //m_Umap(dofMappers.front()),
        //m_Pmap(dofMappers.back())
    {
        m_evaluatorType = "";
        m_pEvaluator = NULL;
        m_elemValue.setZero(2);
    }

    virtual void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        const gsAssemblerOptions& options,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        const gsBasis<T>& basis = basisRefs.back();

        m_dim = basisRefs.front().dim();
        m_patchIndex = patchIndex;

        // Setup Quadrature
        rule = gsGaussRule<T>(basis, options.quA, options.quB);// harmless slicing occurs here

        if (getTMEvaluator() != NULL)
            getTMEvaluator()->initialize(m_viscosity, rule.numNodes());

        GISMO_ASSERT(this->m_bSolutionSet, "No solution set in the residuum visitor.");

        initializeSpecific(evFlags);
    }

    virtual void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        const gsBasis<T>& basis = basisRefs.back();

        m_dim = basisRefs.front().dim();
        m_patchIndex = patchIndex;

        gsVector<index_t> numQuadNodes(m_dim);
        for (int i = 0; i < m_dim; ++i) // to do: improve
            numQuadNodes[i] = basis.maxDegree() + 1; // take quadrature from highest degree

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        if (getTMEvaluator() != NULL)
            getTMEvaluator()->initialize(m_viscosity, rule.numNodes());

        GISMO_ASSERT(this->m_bSolutionSet, "No solution set in the residuum visitor.");

        initializeSpecific(evFlags);
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element nodes
        basisRefs.front().active_into(quNodes.col(0), m_activesU);
        basisRefs.back().active_into(quNodes.col(0), m_activesKO);

        basisRefs.front().deriv_into(quNodes, m_basisGradsU);
        basisRefs.back().evalAllDers_into(quNodes, 2, m_basisDataKO);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate velocity solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);

        const index_t numActKOmega = m_activesKO.rows();
        const index_t numActU = m_activesU.rows();

        // Evaluate solution on element nodes
        m_solKOmegaVals = m_solTM.value(quNodes, m_patchIndex);
        gsMatrix<T> solKOmegaValsOld = m_solKOmegaVals;
        if (m_bIsUnsteady)
        {
            GISMO_ASSERT(m_bOldSolSet, "Old time step solution not set in Crosswind visitor.");
            solKOmegaValsOld = m_oldSolTM.value(quNodes, m_patchIndex);
        }

        gsMatrix<T> solActUCoeffs;
        solActUCoeffs.setZero(m_dim, numActU);
        for (int j = 0; j < numActU; j++) {
            solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();
        }
        gsMatrix<T> solActKOmegaCoeffs;
        solActKOmegaCoeffs.setZero(m_numVar, numActKOmega);
        for (int j = 0; j < numActKOmega; j++)
            solActKOmegaCoeffs.col(j) = m_solTM.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();

        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        index_t nQuPoints = quNodes.cols();
        m_solUGrads.resize(nQuPoints);
        m_solKOmegaGrads.resize(nQuPoints);
        gsMatrix<T> physGrad_u, physGradKO;;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, m_basisGradsU, physGrad_u);
            m_solUGrads[k].noalias() = solActUCoeffs * physGrad_u.transpose();
            geoEval.transformGradients(k, bGrads_komega, physGradKO);
            m_solKOmegaGrads[k].noalias() = solActKOmegaCoeffs * physGradKO.transpose();
        }

        getTMEvaluator()->initAtElement(m_solUGrads, m_solKOmegaVals, m_solKOmegaGrads);
        if (this->checkWallDistanceBasedTM())
        {
            getTMEvaluator()->setKOmegaVariant(m_evaluatorType);

            unsigned evFlagsPoisson = NEED_MEASURE | NEED_GRAD_TRANSFORM;
            const gsGeometry<T>& geo = m_patchesPoisson.patch(m_patchIndex);
            typename gsGeometryEvaluator<T>::uPtr geoEvalPoisson(getEvaluator(evFlagsPoisson, geo));

            gsMatrix<index_t> activesPoisson;
            gsMatrix<T> basisGradsPoisson, physGradPoisson;
            m_basesPoisson.basis(m_patchIndex).active_into(quNodes.col(0), activesPoisson);
            const index_t numActPoisson = activesPoisson.rows();
            m_basesPoisson.basis(m_patchIndex).deriv_into(quNodes, basisGradsPoisson);

            geoEvalPoisson->evaluateAt(quNodes);

            m_solPoissonVals = m_solPoisson.value(quNodes, m_patchIndex);
            gsMatrix<T> solActPoissonCoeffs;
            solActPoissonCoeffs.setZero(1, numActPoisson);
            for (int j = 0; j < numActPoisson; j++)
                solActPoissonCoeffs.col(j) = m_solPoisson.coefficientVector(m_patchIndex).row(activesPoisson(j)).transpose();

            std::vector<gsMatrix<T> > solPoissonGrads(nQuPoints);
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEvalPoisson->transformGradients(k, basisGradsPoisson, physGradPoisson);
                solPoissonGrads[k].noalias() = solActPoissonCoeffs * physGradPoisson.transpose();
            }

            getTMEvaluator()->evalWallDistance(m_solPoissonVals, solPoissonGrads);
        }

        getTMEvaluator()->evalQuantities_nonlinearBlocksPart();
        getTMEvaluator()->evalQuantities_rhsPart();

        gsMatrix<T> physLaplacian;

        m_residuumTM.resize(m_numVar);
        for (int var = 0; var < m_numVar; var++)
            m_residuumTM[var].setZero(nQuPoints);

        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformLaplaceHgrad(k, m_basisDataKO[1], m_basisDataKO[2], physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
            gsMatrix<T> solLaplacian = solActKOmegaCoeffs * physLaplacian.transpose();

            m_residuumTM[0](k) = m_solKOmegaGrads[k].row(0) * m_solUVals.col(k) //advection
                           - getTMEvaluator()->getDiffusionCoefficient(k, 0) * solLaplacian(0, 0) //diffusion
                           + getTMEvaluator()->getReactionCoeff(k, 0) * m_solKOmegaVals(1, k) * m_solKOmegaVals(0, k) //reaction
                           - getTMEvaluator()->getRhs(k, 0); //source

            m_residuumTM[1](k) = m_solKOmegaGrads[k].row(1) * m_solUVals.col(k) //advection
                           - getTMEvaluator()->getDiffusionCoefficient(k, 1) * solLaplacian(1, 0) //diffusion
                           + getTMEvaluator()->getReactionCoeff(k, 1) * m_solKOmegaVals(1, k) * m_solKOmegaVals(1, k) // reaction
                           - getTMEvaluator()->getRhs(k, 1) // source
                           - getTMEvaluator()->getBlendCoeff(k, 1) / math::max(m_solKOmegaVals(1, k), math::pow(10, -15)) *
                             m_solKOmegaGrads[k].row(0).dot(m_solKOmegaGrads[k].row(1)); //blend

            if (m_bIsUnsteady)
            {
                m_residuumTM[0](k) += 1./m_timeStep * (m_solKOmegaVals(0, k) - solKOmegaValsOld(0, k));
                m_residuumTM[1](k) += 1./m_timeStep * (m_solKOmegaVals(1, k) - solKOmegaValsOld(1, k));
            }
        }

        /* part of the residuum
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            m_residuumTM[0](k) = (m_solKOmegaGrads[k].row(0)).norm();
            m_residuumTM[1](k) = (m_solKOmegaGrads[k].row(1)).norm();
        }*/

    }

    virtual inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    { GISMO_NO_IMPLEMENTATION }

    void setCurrentSolution(std::vector<gsField<T> >& solutions)
    {
        m_solU = solutions.front();
        m_solTM = solutions.back();
        this->m_bSolutionSet = true;
    }

    void setOldSolutionField(bool unsteady, gsField<T>& solution, T timeStep)
    {
        m_oldSolTM = solution;
        m_bIsUnsteady = unsteady;
        m_timeStep = timeStep;
        m_bOldSolSet = true;
    }

    gsVector<T> getElementValue() const { return m_elemValue; }

    gsVector<T> getElemValsInQuadPoints () const { return m_elemValsQP; }


    virtual uwbTMEvaluatorKOmega<T>* getTMEvaluator()
    {
        GISMO_ASSERT(m_evaluatorType != "", "No evaluator type set in the visitor.");

        if (m_evaluatorType == "koWilcoxLRN")
            return dynamic_cast<uwbTMEvaluatorKOmegaWilcoxLRN<T>*>(m_pEvaluator);
        else // (m_evaluatorType == "koSST") or its variants
            return dynamic_cast<uwbTMEvaluatorKOmegaSST<T>*>(m_pEvaluator);
    }

    void createTMEvaluator(std::string evaluatorType)
    {
        m_evaluatorType = evaluatorType;
        if (m_evaluatorType == "koWilcoxLRN")
            m_pEvaluator = new uwbTMEvaluatorKOmegaWilcoxLRN<T>();
        else
            m_pEvaluator = new uwbTMEvaluatorKOmegaSST<T>();
    }

    bool checkWallDistanceBasedTM()
    {
        return (m_evaluatorType == "koSST" || m_evaluatorType == "koSSTMenter2009" || m_evaluatorType == "koSAS"
                || m_evaluatorType == "koSAS_SS" || m_evaluatorType == "koSAS_SO" || m_evaluatorType == "koSAS_OO");
    }

    void setPoissonSolution(gsField<T>& solPoisson, const gsMultiPatch<T>& patchesPoisson, const gsMultiBasis<T>& basesPoisson)
    { m_solPoisson = solPoisson; m_patchesPoisson = patchesPoisson; m_basesPoisson = basesPoisson; }

protected:
    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;
    }

protected:
    const T m_viscosity;
    bool m_bOldSolSet;
    bool m_bIsUnsteady;
    index_t m_numVar;

    std::string m_evaluatorType;
    uwbTMEvaluator<T>* m_pEvaluator;

    gsVector<T> m_elemValue, m_elemValsQP;
    T m_timeStep;
    gsField<T> m_solU, m_solTM;
    index_t m_dim; // Velocity vector dimension
    int m_patchIndex;
    gsField<T> m_oldSolTM;
    gsMatrix<index_t> m_activesU;
    gsMatrix<T> m_basisGradsU;
    std::vector<gsMatrix<T> > m_solUGrads;

    gsMatrix<index_t> m_activesKO;
    std::vector<gsMatrix<T> > m_basisDataKO;
    gsMatrix<T> m_solKOmegaVals;
    std::vector<gsMatrix<T> > m_solKOmegaGrads;

    gsMatrix<T> m_solUVals;

    gsField<T> m_solPoisson;
    gsMultiPatch<T> m_patchesPoisson;
    gsMultiBasis<T> m_basesPoisson;
    gsMatrix<T> m_solPoissonVals;

    std::vector<gsVector<T> > m_residuumTM;
};

// ============================================================= Residuum TM L2norm ===================================================== //
template <class T>
class uwbLocRefTMResiduumEvaluatorKOL2norm : public uwbLocRefTMResiduumEvaluator<T>
{

public:
    typedef uwbLocRefTMResiduumEvaluator<T> Base;

public:

    uwbLocRefTMResiduumEvaluatorKOL2norm(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t nQuPoints = quWeights.rows();

        gsVector<T> elResTM;
        elResTM.setZero(m_numVar);
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            for (index_t var = 0; var < m_numVar; ++var)
                elResTM(var) += weight * math::pow(m_residuumTM[var](k), 2);
        }
        for (index_t var = 0; var < m_numVar; ++var)
            m_elemValue(var) = math::sqrt(elResTM(var)) * m_elementLength;
    }

protected:
    using Base::m_elemValue;
    using Base::m_elemValsQP;
    using Base::m_residuumTM;
    using Base::m_numVar;

    using uwbLocRefEvaluatorBase<T>::m_elementLength;

};

} // namespace gismo

