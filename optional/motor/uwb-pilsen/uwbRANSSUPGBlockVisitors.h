/** @file uwbRANSSUPGBlockVisitors.h

Author(s): E. Turnerova
*/

#pragma once
#include "uwbRANSBlockVisitors.h"
//#include "uwbINSSUPGBlockVisitors.h"

namespace gismo
{

// ============================================================= PARENT ============================================================= //
template <class T>
class uwbRANSSUPGBlockVisitor : public uwbRANSBlockVisitor<T>
{

public:
    typedef uwbRANSBlockVisitor<T> Base;

public:
    uwbRANSSUPGBlockVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_degU(0), m_tauStabType(2), m_timeStep(0.)
    {
        m_pStabEvaluator = new uwbStabilizationEvaluator<T>();
        m_bAverageTau = false;
        m_dirElemLength = false;
    }
    ~uwbRANSSUPGBlockVisitor()
    {
        /*if (m_pTMsolver)
        {
            delete m_pTMsolver;
            m_pTMsolver = NULL;
        }*/
        if (m_pStabEvaluator)
        {
            delete m_pStabEvaluator;
            m_pStabEvaluator = NULL;
        }
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        const gsAssemblerOptions& options,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        Base::initialize(basisRefs, patchIndex, options, rule, evFlags);
        if (m_pStabEvaluator != NULL)
            m_pStabEvaluator->initialize(rule.numNodes(), m_dim);
    }

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        Base::initialize(basisRefs, patchIndex, rule, evFlags);
        if (m_pStabEvaluator != NULL)
            m_pStabEvaluator->initialize(rule.numNodes(), m_dim);
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        m_solUVals = m_solU.value(quNodes, m_patchIndex);

        m_degU = basisRefs.front().maxDegree();

        //m_turbViscosityField = m_pTMsolver->constructTurbulentViscositySol();

        index_t nQuPoints = quNodes.cols();

        const index_t numActU = m_activesU.rows();
        gsMatrix<T> solActUCoeffs(m_dim, numActU);
        for (int j = 0; j < numActU; j++)
            solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();

        //-----
        gsMatrix<T> physGradU, solUGrad;

        std::vector<gsVector<T>> h_advection, h_diffusion;
        h_advection.resize(m_dim);
        for (int s = 0; s < m_dim; s++)
            h_advection[s].setZero(nQuPoints);
        h_diffusion = h_advection;
        if (m_bDirElemLength)
        {
            if (m_hDirType == 0) //h_RGN
            {
                gsMatrix<T> solNormGrad;
                T gradProduct;
                T solNormGradNorm;
                for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
                {
                    geoEval.transformGradients(k, m_parGradsU, physGradU);
                    solUGrad.noalias() = solActUCoeffs * physGradU.transpose();

                    solNormGrad.setZero(1, m_dim);
                    gsMatrix<T> product;
                    for (index_t var = 0; var < m_dim; var++)
                    {
                        product = ((m_solUVals.col(k)).transpose() * solUGrad.col(var)) / (m_solUVals.col(k)).norm();
                        solNormGrad(0, var) = product(0, 0);
                    }
                    solNormGradNorm = (solNormGrad.row(0)).norm();

                    gradProduct = 0.;
                    for (index_t a = 0; a < numActU; a++)
                        if (solNormGradNorm != 0.)
                            gradProduct += math::abs(solNormGrad.row(0) * physGradU.col(a)) / solNormGradNorm;

                    if (gradProduct == 0.)
                        h_advection[0](k) = m_elementLength;
                    else
                        h_advection[0](k) = 2. / gradProduct;
                    if (math::isnan(h_advection[0](k)))
                        h_advection[0](k) = m_elementLength;
                }
                for (int var = 1; var < m_dim; var++)
                    h_advection[var] = h_advection[0];

                /*gsInfo << "h_RGN RANS\n";
                gsInfo << "m_elementLength" << m_elementLength << "\n";
                gsInfo << "h_advection" << h_advection[0] << "\n";*/
            }//-------------------------------------------------------------------------------
            else if (m_hDirType == 1) //h_UGN
            {
                T gradProduct;
                for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
                {
                    geoEval.transformGradients(k, m_parGradsU, physGradU);

                    gsVector<T> advection = m_solUVals.col(k);

                    gradProduct = 0.;
                    for (index_t a = 0; a < numActU; a++)
                        gradProduct += math::abs(advection.transpose() * physGradU.col(a));

                    if (gradProduct == 0.)
                        h_advection[0](k) = m_elementLength;
                    else
                        h_advection[0](k) = 2. * advection.norm() / gradProduct;
                    if (math::isnan(h_advection[0](k)))
                        h_advection[0](k) = m_elementLength;
                }
                for (int var = 1; var < m_dim; var++)
                    h_advection[var] = h_advection[0];
                /*gsInfo << "h_RANS\n";
                gsInfo << "m_elementLength" << m_elementLength << "\n";
                gsInfo << "h_advection u" << h_advection[0] << "\n";*/
            }//------------------------------------------------------------------------------
            else
                GISMO_ERROR("Wrong hDirType selected.");
        }
        else
        {
            for (index_t k = 0; k < nQuPoints; ++k)
                for (index_t var = 0; var < m_dim; var++)
                    h_advection[var](k) = m_elementLength;
        }
        h_diffusion = h_advection;
        //-----

        GISMO_ENSURE(m_bEffectiveViscSet, "uwbRANSSUPGBlockVisitor: effective viscosity not set!");
        m_pStabEvaluator->initAtElement(m_solUVals, this->m_diffusionCoeff, m_bTauDeg);
        m_pStabEvaluator->setSUPGvars(m_tauStabType, m_degU, m_timeStep);
        m_pStabEvaluator->setElemLength(m_elementLength, h_advection, h_diffusion);
    }

    void setTauStabType(const int tauStabType, const T timeStep = 0)
    { m_tauStabType = tauStabType; m_timeStep = timeStep; }

protected:
    int m_degU;
    int m_tauStabType;
    T m_timeStep;

    bool m_bAverageTau;
    bool m_dirElemLength;

    gsMatrix<T> m_solUVals;

    uwbStabilizationEvaluator<T>* m_pStabEvaluator;

    //gsField<T> m_turbViscosityField;

    // Members from Base (uwbRANSBlockVisitor)
    //using uwbRANSBlockVisitor<T>::m_pTMsolver;
    using Base::m_parGradsU;
    using Base::m_activesU;
    using Base::m_bEffectiveViscSet;

    // Members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_solU;
    using uwbINSBlockVisitor<T>::m_patchIndex;
    using uwbINSBlockVisitor<T>::m_viscosity;
    using uwbINSBlockVisitor<T>::m_dim; // Velocity vector dimension

    // Members from uwbVisitorBase
    using uwbVisitorBase<T>::m_elementLength;
    using uwbVisitorBase<T>::m_bDirElemLength;
    using uwbVisitorBase<T>::m_hDirType;
    using uwbVisitorBase<T>::m_bTauDeg;
};

// ============================================================= BLOCK A_RANS SUPG ============================================================= //

template <class T>
class uwbRANSSUPGBlockAVisitor : public uwbRANSSUPGBlockVisitor<T>
{

public:
    typedef uwbRANSSUPGBlockVisitor<T> Base;

public:

    uwbRANSSUPGBlockAVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        gsMatrix<T> deriv2_u;
        basisRefs.front().deriv2_into(quNodes, deriv2_u);

        m_basisDataU.push_back(m_parGradsU);
        m_basisDataU.push_back(deriv2_u);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localMat.setZero(numActU, numActU);

        const gsMatrix<T> & basisGradsU = m_basisDataU[0]; //first partial derivatives
        const gsMatrix<T> & bHessian_u = m_basisDataU[1]; //second partial derivatives

        const index_t nQuPoints = quWeights.rows();

        gsVector<T> tau_s;
        tau_s.setZero(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k)
            tau_s(k) = m_pStabEvaluator->getTauS(k);

        if (Base::m_bAverageTau)
        {
            T sumTau = 0.;
            for (index_t k = 0; k < nQuPoints; ++k)
                sumTau += tau_s(k);
            tau_s.setConstant(nQuPoints, sumTau / nQuPoints);
        }

        gsMatrix<T> physLaplacianU;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical second partial derivatives at k
            geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacianU); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n

            // Local block
            //T tau_s = m_pStabEvaluator->getTauS(k);
            m_localMat.noalias() -= weight * getTurbViscosityVal(k) * tau_s(k) * ((m_physGradsU[k].transpose() * m_solUVals.col(k)) * physLaplacianU);
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {
        const index_t usz = m_Umap.freeSize();

        // Local Dofs to global dofs
        m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
        const index_t numActU = m_activesU.rows();

        for (index_t i = 0; i < numActU; ++i)
        {
            const int ii = m_activesU(i);
            if (m_Umap.is_free_index(ii))
            {
                for (index_t j = 0; j < numActU; ++j)
                {
                    const int jj = m_activesU(j);
                    if (m_Umap.is_free_index(jj))
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                    else // m_Umap.is_boundary_index(jj)
                    {
                        const int bb = m_Umap.global_to_bindex(jj);
                        for (index_t s = 0; s != m_dim; ++s)
                            rhs(s*usz + ii, 0) -= m_localMat(i, j) * eliminatedDofs[0](bb, s); // assuming single rhs
                    }
                }
            }
        }
    }

protected:

    gsMatrix<T> m_localMat; // Local matrix
    std::vector<gsMatrix<T> > m_basisDataU;

    // Members from Base (uwbRANSSUPGBlockVisitor)
    using Base::m_solUVals;
    using Base::m_pStabEvaluator;

    // Members and functions from uwbRANSBlockVisitor
    using uwbRANSBlockVisitor<T>::m_parGradsU;
    using uwbRANSBlockVisitor<T>::m_activesU;
    using uwbRANSBlockVisitor<T>::m_physGradsU;
    using uwbRANSBlockVisitor<T>::getTurbViscosityVal;

    // Members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_Umap;
    using uwbINSBlockVisitor<T>::m_dim; // Velocity vector dimension
    using uwbINSBlockVisitor<T>::m_patchIndex;
};

} // namespace gismo
