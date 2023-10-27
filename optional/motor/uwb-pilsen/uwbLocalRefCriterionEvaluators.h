/** @file uwbLocalRefCriterionEvaluators.h

    Author(s): E. Turnerova, K. Slaba
*/

#pragma once

namespace gismo
{
// ========================================================== SUPER CLASS ========================================================== //
template <class T>
class uwbLocRefEvaluatorBase
{
public:
    uwbLocRefEvaluatorBase()
    {
        m_bSolutionSet = false;
        m_elementLength = 0.;
    }

    virtual void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { GISMO_NO_IMPLEMENTATION }
    virtual void setOldSolutionField(bool unsteady, gsField<T>& solution, T timeStep)
    { GISMO_NO_IMPLEMENTATION }

    void setElementLength(T elementLength) { m_elementLength = elementLength; }

protected: 
    bool m_bSolutionSet;
    T m_elementLength;
};

// ============================================================= PARENT ============================================================= //
template <class T>
class uwbLocRefEvaluator : public uwbLocRefEvaluatorBase<T>
{
public:
    uwbLocRefEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        m_viscosity(viscosity)
        //,m_Umap(dofMappers.front()),
        //m_Pmap(dofMappers.back())
    {
        m_areaVals.setZero(2);
    }

    virtual void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        const gsAssemblerOptions& options,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        const gsBasis<T>& basis = basisRefs.front();

        m_dim = basis.dim();
        m_patchIndex = patchIndex;

        // Setup Quadrature
        rule = gsGaussRule<T>(basis, options.quA, options.quB);// harmless slicing occurs here

        initializeSpecific(evFlags);
    }

    virtual void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        const gsBasis<T>& basis = basisRefs.front();

        m_dim = basis.dim();
        m_patchIndex = patchIndex;

        gsVector<index_t> numQuadNodes(m_dim);
        for (int i = 0; i < m_dim; ++i) // to do: improve
            numQuadNodes[i] = basis.maxDegree() + 1; // take quadrature from highest degree

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        initializeSpecific(evFlags);
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        gsMatrix<index_t> activesP;
        basisRefs.front().active_into(quNodes.col(0), m_activesU);
        basisRefs.back().active_into(quNodes.col(0), activesP);

        basisRefs.front().evalAllDers_into(quNodes, 1, m_basisDataU);

        gsMatrix<T> basisGradsP;
        basisRefs.back().deriv_into(quNodes, basisGradsP);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate velocity gradients
        index_t numActU = m_activesU.rows();
        const index_t numActP = activesP.rows();
        m_nQuPoints = quNodes.cols();

        m_solActUCoeffs.resize(m_dim, numActU);
        for (int j = 0; j < numActU; j++)
            m_solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();

        gsMatrix<T> solActPCoeffs;
        solActPCoeffs.setZero(1, numActP);
        for (int j = 0; j < numActP; j++)
            solActPCoeffs.col(j) = m_solP.coefficientVector(m_patchIndex).row(activesP(j)).transpose();

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);
        m_solPVals = m_solP.value(quNodes, m_patchIndex);

        const gsMatrix<T> & basisGradsU = m_basisDataU[1];

        m_solUGrads.resize(m_nQuPoints);
        m_solPGrads.resize(m_nQuPoints);
        gsMatrix<T> physGradU, physGradP;
        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, basisGradsU, physGradU);
            m_solUGrads[k].noalias() = m_solActUCoeffs * physGradU.transpose();
            geoEval.transformGradients(k, basisGradsP, physGradP);
            m_solPGrads[k].noalias() = solActPCoeffs * physGradP.transpose();
        }
    }

    virtual inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    { GISMO_NO_IMPLEMENTATION }

    void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { 
        m_solU = solutions.front();
        m_solP = solutions.back();
        this->m_bSolutionSet = true;
    }

    gsVector<T> getElementValue() const { return m_elemValue; }

    gsVector<T> getElemValsInQuadPoints () const { return m_elemValsQP; }

    gsVector<T> getAreaVals() const { return m_areaVals; }

protected:
    virtual inline void initializeSpecific(unsigned & evFlags)
    { GISMO_NO_IMPLEMENTATION }

protected:

    gsVector<T> m_elemValue, m_elemValsQP, m_areaVals;
    const T m_viscosity;
    gsField<T> m_solU, m_solP;
    index_t m_dim; // Velocity vector dimension
    int m_patchIndex;
    index_t m_nQuPoints;

    gsMatrix<index_t> m_activesU;
    std::vector<gsMatrix<T> > m_basisDataU;
    std::vector<gsMatrix<T> > m_solUGrads, m_solPGrads;

    gsMatrix<T> m_solUVals, m_solPVals;

    gsMatrix<T> m_solActUCoeffs;
};

// ============================================================= Residuum NS ===================================================== //
template <class T>
class uwbLocRefINSResiduumEvaluator : public uwbLocRefEvaluator<T>
{

public:
    typedef uwbLocRefEvaluator<T> Base;

public:

    uwbLocRefINSResiduumEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_bIsUnsteady(false), m_bOldSolutionSet(false)
    {
        m_elemValue.setZero(1);
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity and pressure solution set in the locRef evaluator.");
        if (m_bIsUnsteady)
            GISMO_ASSERT(m_bOldSolutionSet, "No old time step velocity and pressure solution set in the locRef evaluator.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        gsMatrix<T> deriv2_u;
        basisRefs.front().deriv2_into(quNodes, deriv2_u);

        m_basisDataU.push_back(deriv2_u);

        const gsMatrix<T> & basisGradsU = m_basisDataU[1];
        const gsMatrix<T> & bHessian_u = m_basisDataU[2]; //-> obsahuje i smisenou druhou derivaci

        gsMatrix<T> physLaplacian;

        m_solUValsOld = m_solUVals;
        if (m_bIsUnsteady)
            m_solUValsOld = m_oldSolU.value(quNodes, m_patchIndex);

        m_residuumMomentumEq.resize(m_dim);
        for (int var = 0; var < m_dim; var++)
            m_residuumMomentumEq[var].setZero(m_nQuPoints);
        m_residuumContinuityEq.setZero(m_nQuPoints);

        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
            gsMatrix<T> solLaplacian = m_solActUCoeffs * physLaplacian.transpose();

            for (int var = 0; var < m_dim; var++)
            {
                m_residuumMomentumEq[var](k) = m_solUGrads[k].row(var) * m_solUVals.col(k) //advection
                                              - m_viscosity * solLaplacian(var, 0) //diffusion
                                              + m_solPGrads[k](0, var); //pressure term
                                              //- f; //source term

                if (m_bIsUnsteady)
                    m_residuumMomentumEq[var](k) += 1./m_timeStep * (m_solUVals(var, k) - m_solUValsOld(var, k));

                m_residuumContinuityEq(k) += m_solUGrads[k](var, var);
            }
        }
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        T elResMomentum = 0.;
        T elResContinuity = 0.;
        m_elemValsQP.setZero(m_nQuPoints);
        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);
            //T tmp = 0.;
            for (index_t var = 0; var != m_dim; ++var)
                m_elemValsQP(k) += math::pow(m_residuumMomentumEq[var](k), 2);
            elResMomentum += weight * m_elemValsQP(k);

            m_elemValsQP(k) = math::sqrt(m_elemValsQP(k));
            //tmp = math::pow(m_residuumContinuityEq(k), 2);
            m_elemValsQP(k) += math::abs(m_residuumContinuityEq(k));//math::sqrt(tmp);
            m_elemValsQP(k) *= m_elementLength;

            elResContinuity += weight * math::pow(m_residuumContinuityEq(k), 2);
        }
        m_elemValue(0) = (math::sqrt(elResMomentum) + math::sqrt(elResContinuity)) * m_elementLength;
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

protected:
    using Base::m_elemValue;
    using Base::m_elemValsQP;
    using Base::m_viscosity;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_solActUCoeffs;
    using Base::m_solUGrads;
    using Base::m_solPGrads;
    using Base::m_basisDataU;
    using Base::m_solUVals;
    using Base::m_patchIndex;
    using Base::m_nQuPoints;

    using uwbLocRefEvaluatorBase<T>::m_elementLength;
};


// ============================================================= Residuum momentum NS ===================================================== //
template <class T>
class uwbLocRefINSResiduumMomentumEvaluator : public uwbLocRefEvaluator<T>
{

public:
    typedef uwbLocRefEvaluator<T> Base;

public:

    uwbLocRefINSResiduumMomentumEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_bIsUnsteady(false), m_bOldSolutionSet(false)
    {
        m_elemValue.setZero(1);
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity and pressure solution set in the locRef evaluator.");
        if (m_bIsUnsteady)
            GISMO_ASSERT(m_bOldSolutionSet, "No old time step velocity and pressure solution set in the locRef evaluator.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        gsMatrix<T> deriv2_u;
        basisRefs.front().deriv2_into(quNodes, deriv2_u);

        m_basisDataU.push_back(deriv2_u);

        const gsMatrix<T> & basisGradsU = m_basisDataU[1];
        const gsMatrix<T> & bHessian_u = m_basisDataU[2]; //-> obsahuje i smisenou druhou derivaci

        gsMatrix<T> physLaplacian;

        m_solUValsOld = m_solUVals;
        if (m_bIsUnsteady)
            m_solUValsOld = m_oldSolU.value(quNodes, m_patchIndex);

        m_residuumMomentumEq.resize(m_dim);
        for (int var = 0; var < m_dim; var++)
            m_residuumMomentumEq[var].setZero(m_nQuPoints);

        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
            gsMatrix<T> solLaplacian = m_solActUCoeffs * physLaplacian.transpose();

            for (int var = 0; var < m_dim; var++)
            {
                m_residuumMomentumEq[var](k) = m_solUGrads[k].row(var) * m_solUVals.col(k) //advection
                                              - m_viscosity * solLaplacian(var, 0) //diffusion
                                              + m_solPGrads[k](0, var); //pressure term
                                              //- f; //source term

                if (m_bIsUnsteady)
                    m_residuumMomentumEq[var](k) += 1./m_timeStep * (m_solUVals(var, k) - m_solUValsOld(var, k));
            }
        }
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        T elResMomentum = 0.;
        m_elemValsQP.setZero(m_nQuPoints);
        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);
            for (index_t var = 0; var != m_dim; ++var)
                m_elemValsQP(k) += math::pow(m_residuumMomentumEq[var](k), 2);
            elResMomentum += weight * m_elemValsQP(k);

            m_elemValsQP(k) = math::sqrt(m_elemValsQP(k));
            m_elemValsQP(k) *= m_elementLength;
        }
        m_elemValue(0) = math::sqrt(elResMomentum) * m_elementLength;
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

protected:
    using Base::m_elemValue;
    using Base::m_elemValsQP;
    using Base::m_viscosity;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_solActUCoeffs;
    using Base::m_solUGrads;
    using Base::m_solPGrads;
    using Base::m_basisDataU;
    using Base::m_solUVals;
    using Base::m_patchIndex;
    using Base::m_nQuPoints;

    using uwbLocRefEvaluatorBase<T>::m_elementLength;
};

// ============================================================= residuum continuity ===================================================== //
template <class T>
class uwbLocRefResiduumContinuityEvaluator : public uwbLocRefEvaluator<T>
{

public:
    typedef uwbLocRefEvaluator<T> Base;

public:

    uwbLocRefResiduumContinuityEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    {
        m_elemValue.setZero(1);
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity and pressure solution set in the locRef evaluator.");
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        m_elemValue(0) = 0.;
        m_elemValsQP.setZero(m_nQuPoints);
        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            for (index_t var = 0; var != m_dim; ++var)
                m_elemValsQP(k) += m_solUGrads[k](var, var);
            m_elemValue(0) += weight * m_elemValsQP(k);
        }
        m_elemValue(0) *= m_elementLength;
    }

protected:
    using Base::m_elemValue;
    using Base::m_elemValsQP;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_solUGrads;
    using Base::m_nQuPoints;

    using uwbLocRefEvaluatorBase<T>::m_elementLength;
};

// ============================================================= Sistek ===================================================== //
template <class T>
class uwbLocRefSistekEvaluator : public uwbLocRefEvaluator<T>
{

public:
    typedef uwbLocRefEvaluator<T> Base;

public:

    uwbLocRefSistekEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_bIsUnsteady(false), m_bOldSolutionSet(false)
    {
        m_elemValue.setZero(1);//numerator
        m_areaVals.setZero(2);//denominator, element area
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity and pressure solution set in the locRef evaluator.");
        if (m_bIsUnsteady)
            GISMO_ASSERT(m_bOldSolutionSet, "No old time step velocity and pressure solution set in the locRef evaluator.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        gsMatrix<T> deriv2_u;
        basisRefs.front().deriv2_into(quNodes, deriv2_u);

        m_basisDataU.push_back(deriv2_u);

        const gsMatrix<T> & basisGradsU = m_basisDataU[1];
        const gsMatrix<T> & bHessian_u = m_basisDataU[2]; //-> obsahuje i smisenou druhou derivaci

        gsMatrix<T> physLaplacian;

        m_solUValsOld = m_solUVals;
        if (m_bIsUnsteady)
            m_solUValsOld = m_oldSolU.value(quNodes, m_patchIndex);

        m_residuumMomentumEq.resize(m_dim);
        for (int var = 0; var < m_dim; var++)
            m_residuumMomentumEq[var].setZero(m_nQuPoints);
        m_residuumContinuityEq.setZero(m_nQuPoints);

        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
            gsMatrix<T> solLaplacian = m_solActUCoeffs * physLaplacian.transpose();

            for (int var = 0; var < m_dim; var++)
            {
                m_residuumMomentumEq[var](k) = m_solUGrads[k].row(var) * m_solUVals.col(k) //advection
                                              - m_viscosity * solLaplacian(var, 0) //diffusion
                                              + m_solPGrads[k](0, var); //pressure term
                                              //- f; //source term

                if (m_bIsUnsteady)
                    m_residuumMomentumEq[var](k) += 1./m_timeStep * (m_solUVals(var, k) - m_solUValsOld(var, k));

                m_residuumContinuityEq(k) += m_solUGrads[k](var, var);
            }
        }
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        T elResMomentum = 0.;
        T elResContinuity = 0.;
        T elP = 0.;
        T elVel = 0.;
        m_elemValsQP.setZero(m_nQuPoints);
        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);
            //----- numerator -----
            for (index_t var = 0; var != m_dim; ++var)
                m_elemValsQP(k) += math::pow(m_residuumMomentumEq[var](k), 2);
            elResMomentum += weight * m_elemValsQP(k);

            m_elemValsQP(k) = math::pow(m_elementLength, 2) * m_elemValsQP(k);
            m_elemValsQP(k) += math::pow(m_residuumContinuityEq(k), 2);

            elResContinuity += weight * math::pow(m_residuumContinuityEq(k), 2);
            //---- denominator ----
            T tmp = 0.;
            for (index_t var = 0; var != m_dim; ++var)
            {
                for (index_t col = 0; col != m_dim; col++)
                    tmp += math::pow(m_solUGrads[k](var, col), 2);

                tmp += math::pow(m_solUVals(var, k), 2);
            }
            elVel += weight * tmp;

            elP += weight * math::pow(m_solPVals(0, k), 2);
            //---- element area ----
            m_areaVals(1) += weight;
        }
        m_elemValue(0) = (math::pow(m_elementLength, 2) * elResMomentum + elResContinuity) / m_areaVals(1); //numerator / elemArea
        m_areaVals(0) = elVel + elP; //denominator
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

protected:
    using Base::m_elemValue;
    using Base::m_elemValsQP;
    using Base::m_areaVals;
    using Base::m_viscosity;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_solActUCoeffs;
    using Base::m_solPVals;
    using Base::m_solUGrads;
    using Base::m_solPGrads;
    using Base::m_basisDataU;
    using Base::m_solUVals;
    using Base::m_patchIndex;
    using Base::m_nQuPoints;

    using uwbLocRefEvaluatorBase<T>::m_elementLength;

};



// ============================================================= grad u ===================================================== //
template <class T>
class uwbLocRefGradientUEvaluator : public uwbLocRefEvaluator<T>
{

public:
    typedef uwbLocRefEvaluator<T> Base;

public:

    uwbLocRefGradientUEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    {
        m_elemValue.setZero(1);
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity and pressure solution set in the locRef evaluator.");
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        m_elemValue(0) = 0.;
        m_elemValsQP.setZero(m_nQuPoints);
        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            for (index_t var = 0; var != m_dim; ++var)
                m_elemValsQP(k) += (m_solUGrads[k].row(var)).norm();
            m_elemValue(0) += weight * m_elemValsQP(k);
        }
        m_elemValue(0) *= m_elementLength;
    }

protected:
    using Base::m_elemValue;
    using Base::m_elemValsQP;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_solUGrads;
    using Base::m_nQuPoints;

    using uwbLocRefEvaluatorBase<T>::m_elementLength;
};


// ============================================================= grad p ===================================================== //
template <class T>
class uwbLocRefGradientPEvaluator : public uwbLocRefEvaluator<T>
{

public:
    typedef uwbLocRefEvaluator<T> Base;

public:

    uwbLocRefGradientPEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    {
        m_elemValue.setZero(1);
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity and pressure solution set in the locRef evaluator.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
            gsGeometryEvaluator<T> & geoEval,
            gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        gsMatrix<index_t> activesP;
        basisRefs.back().active_into(quNodes.col(0), activesP);

        gsMatrix<T> basisGradsP;
        basisRefs.back().deriv_into(quNodes, basisGradsP);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        const index_t numActP = activesP.rows();
        m_nQuPoints = quNodes.cols();

        gsMatrix<T> solActPCoeffs;
        solActPCoeffs.setZero(1, numActP);
        for (int j = 0; j < numActP; j++)
            solActPCoeffs.col(j) = m_solP.coefficientVector(m_patchIndex).row(activesP(j)).transpose();

        m_solPGrads.resize(m_nQuPoints);
        gsMatrix<T> physGradP;
        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, basisGradsP, physGradP);
            m_solPGrads[k].noalias() = solActPCoeffs * physGradP.transpose();
        }
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        m_elemValue(0) = 0.;
        m_elemValsQP.setZero(m_nQuPoints);
        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);
            m_elemValsQP(k) = m_solPGrads[k].norm();
            m_elemValue(0) += weight * m_elemValsQP(k);
        }
        m_elemValue(0) *= m_elementLength;
    }

protected:
    using Base::m_elemValue;
    using Base::m_elemValsQP;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_solP;
    using Base::m_solPGrads;
    using Base::m_nQuPoints;
    using Base::m_patchIndex;

    using uwbLocRefEvaluatorBase<T>::m_elementLength;
};



// ============================================================= vorticity ===================================================== //
template <class T>
class uwbLocRefVorticityEvaluator : public uwbLocRefEvaluator<T>
{

public:
    typedef uwbLocRefEvaluator<T> Base;

public:

    uwbLocRefVorticityEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    {
        m_elemValue.setZero(1);
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity and pressure solution set in the locRef evaluator.");
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        if (m_dim != 2)
            GISMO_ERROR("Vorticity is implemented only for 2D so far.");

        m_elemValue(0) = 0.;
        m_elemValsQP.setZero(m_nQuPoints);
        for (index_t k = 0; k < m_nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);
            m_elemValsQP(k) = m_solUGrads[k](1, 0) - m_solUGrads[k](0, 1);
            m_elemValue(0) += weight * m_elemValsQP(k);
        }
    }

protected:
    using Base::m_elemValue;
    using Base::m_elemValsQP;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_nQuPoints;
    using Base::m_solUGrads;
    //using Base::m_solPGrads;

};


} // namespace gismo

