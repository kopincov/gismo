/** @file uwbRANSAssembler.h

Author(s): E. Turnerova, H. Hornikova
*/

#pragma once

#include "uwbINSAssemblerUnsteady.h"
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

template<class T>
class uwbRANSAssembler : public uwbINSAssemblerUnsteady<T>
{

public:
    typedef uwbINSAssemblerUnsteady<T> Base;

public:
    uwbRANSAssembler(uwbINSSolverParams<T>& params, uwbTMSolverBase<T>* pTMsolver) : 
        Base(params), m_pTurbulenceSolver(pTMsolver)
    {
        m_blockAssembler.setRANS();
    }

    virtual ~uwbRANSAssembler()
    {
    }

public:

    virtual void updateAssembly()
    {
        Base::updateAssembly();

        m_blockAssembler.assembleRANSPart(m_solution);

        fillNMatrix();

        if (m_blockAssembler.isRANScrosswind())
        {
            m_blockAssembler.assembleRANScrosswindPart(m_solution);

            fillRANScrosswindMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsCrosswind();
        }
        if (m_blockAssembler.isSRBAV())
        {
            m_blockAssembler.assembleRANSsrbavPart(m_solution);
            fillRANSsrbavMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsSRBAV();
        }
        if (m_blockAssembler.isRANSisoAD())
        {
            m_blockAssembler.assembleRANSisoADpart(m_solution);
            fillRANSisoADmatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsIsoAD();
        }
        if (m_blockAssembler.isRANSad())
        {
            m_blockAssembler.assembleBlockRANSad(m_solution);

            fillRANSadMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhs_RANS_AD();
        }
    }

    virtual void updatePicardAssembly(const gsMatrix<T> & solVector)
    {
        Base::updatePicardAssembly(solVector);

        m_baseMatrix = m_NMatrix;
        m_baseRhs = m_NRhs;

        if (m_blockAssembler.isRANScrosswind())
        {
            m_blockAssembler.assembleRANScrosswindPart(solVector, false);

            fillRANScrosswindMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsCrosswind();
        }
        if (m_blockAssembler.isSRBAV())
        {
            m_blockAssembler.assembleRANSsrbavPart(solVector, false);
            fillRANSsrbavMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsSRBAV();
        }
        if (m_blockAssembler.isRANSisoAD())
        {
            m_blockAssembler.assembleRANSisoADpart(solVector, false);
            fillRANSisoADmatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsIsoAD();
        }
        if (m_blockAssembler.isRANSad())
        {
            m_blockAssembler.assembleBlockRANSad(solVector, false);

            fillRANSadMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhs_RANS_AD();
        }
    }

protected:

    virtual void fillBase()
    {
        Base::fillBase();

        m_StokesMatrix = m_baseMatrix;
        m_StokesRhs = m_baseRhs;
    }
    
    virtual void fillNMatrix()
    {
        m_baseMatrix = m_StokesMatrix;
        m_baseRhs = m_StokesRhs;

        int udofs = m_blockAssembler.getUdofs();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < udofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockA_RANS(), col); it; ++it)
                for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
                    m_baseMatrix.coeffRef(s*udofs + it.row(), s*udofs + it.col()) += it.value();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < udofs; ++col)
            for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockEdiag_RANS(s), col); it; ++it)
                    m_baseMatrix.coeffRef(s*udofs + it.row(), s*udofs + it.col()) += it.value();

        //m_baseRhs += m_blockAssembler.getRhsA_RANS();
        m_baseRhs += m_blockAssembler.getRhsA_RANS() + m_blockAssembler.getRhsEdiag_RANS() + m_blockAssembler.getRhsEnondiag_RANS();

/*        T kKoef = -2 / 3;
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
            m_baseRhs.middleRows(s * udofs, udofs).noalias() += kKoef * m_blockAssembler.getBlockMinusBT(s) * m_pTurbulenceSolver->getSolutionK_full(m_blockAssembler.getMappers().back());
*/

        if (m_blockAssembler.isSUPG())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < udofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockA_RANS_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
                        m_baseMatrix.coeffRef(s*udofs + it.row(), s*udofs + it.col()) += it.value();

            m_baseRhs += m_blockAssembler.getRhsA_RANS_SUPG();

            //#pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            //for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
            //    m_baseRhs.middleRows(s * udofs, udofs).noalias() += kKoef * m_blockAssembler.getBlockB_SUPG(s) * m_pTurbulenceSolver->getSolutionK_full(m_blockAssembler.getMappers().back());

        }

        m_NMatrix = m_baseMatrix;
        m_NRhs = m_baseRhs;
    }

    void fillRANScrosswindMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockCrosswind(s), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

    void fillRANSsrbavMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockSRBAV(s), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

    void fillRANSisoADmatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockIsoAD(s), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

    void fillRANSadMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlock_RANS_AD(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

public:
    void plotResiduum(std::string const & fn, gsMatrix<T> & solution, gsMatrix<T>& solutionOld, gsMatrix<T>& TMsolution, unsigned npts,
                      bool tanh = false, bool scale = false, bool tanhScale = false, bool print = false)
    {
        std::string fnM = fn + "MomentumEq";
        std::string fnC = fn + "ContinuityEq";

        gsParaviewCollection collectionM(fnM);
        gsParaviewCollection collectionC(fnC);
        std::string fileNameM, fileNameC;

        gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);
        gsField<T> pSolField = m_blockAssembler.constructSolution(solution, 1);
        gsField<T> uSolFieldOld = m_blockAssembler.constructSolution(solutionOld, 0);

        m_averageResidual.setConstant(m_blockAssembler.getPatches().patch(0).targetDim(), 0.);
        m_averageAbsResidual.setConstant(m_blockAssembler.getPatches().patch(0).targetDim(), 0.);
        T numPoints = 0.;

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            fileNameM = fnM + util::to_string(p);
            fileNameC = fnC + util::to_string(p);

            gsMatrix<T> geoVals, residuumMomentumVals, residuumContinuityVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
            const int parDim = patch->domainDim();
            const int tarDim = patch->targetDim();

            gsMatrix<T> ab = patch->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);
            numPoints += pts.cols();

            evalRANSresiduum_singlePatch(p, pts, uSolField, pSolField, uSolFieldOld, TMsolution, residuumMomentumVals, residuumContinuityVals,
                                         tanh, scale, tanhScale);

            // Evaluate geometry at pts
            geoVals = patch->eval(pts);

            if (3 - parDim > 0)
            {
                np.conservativeResize(3);
                np.bottomRows(3 - parDim).setOnes();
            }
            if (3 - tarDim > 0)
            {
                geoVals.conservativeResize(3, geoVals.cols());
                geoVals.bottomRows(3 - tarDim).setZero();
            }

            gsWriteParaviewTPgrid(geoVals, residuumMomentumVals, np.template cast<index_t>(), fileNameM);
            collectionM.addPart(fileNameM, ".vts");
            if (!tanh && !scale && !tanhScale)
            {
                gsWriteParaviewTPgrid(geoVals, residuumContinuityVals, np.template cast<index_t>(), fileNameC);
                collectionC.addPart(fileNameC, ".vts");
            }
        }

        if (print)
        {
            //gsInfo << "numPoints = " << numPoints << "\n";
            gsInfo << "\naveraged RANS residuum = \n" << m_averageResidual / numPoints << "\n";
            gsInfo << "averaged RANS abs(residuum) = \n" << m_averageAbsResidual / numPoints << "\n\n";
        }

        collectionM.save();
        if (!tanh && !scale && !tanhScale)
            collectionC.save();
    }

protected:
    void evalRANSresiduum_singlePatch(index_t patchIndex, gsMatrix<T>& pts, gsField<T>& uSolField, gsField<T>& pSolField,
                                      gsField<T>& uSolFieldOld, gsMatrix<T>& TMsolution, gsMatrix<T>& residuumMvals,
                                      gsMatrix<T>& residuumCvals, bool tanh = false, bool scale = false, bool tanhScale = false)
    {
        gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
        gsBasisRefs<T> bases(m_blockAssembler.getBases(), patchIndex);

        const T viscosity = m_blockAssembler.getViscosity();

        const int tarDim = patch->targetDim();

        unsigned evFlags = NEED_VALUE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *patch));

        gsMatrix<index_t> activesU, activesP;

        gsMatrix<T> basisGradsU, physGradU;
        bases.front().deriv_into(pts, basisGradsU);
        gsMatrix<T> bHessianU;
        bases.front().deriv2_into(pts, bHessianU);

        gsMatrix<T> basisGradsP, physGradP;
        bases.back().deriv_into(pts, basisGradsP);

        geoEval->evaluateAt(pts);

        gsMatrix<T> solUVals = uSolField.value(pts, patchIndex);
        gsMatrix<T> solUValsOld = uSolFieldOld.value(pts, patchIndex);

        index_t nQuPoints = pts.cols();

        gsMatrix<T> actCoeffsU, actCoeffsP;
        std::vector<gsMatrix<T> > uGrads, pGrads, solLaplacian;
        uGrads.resize(nQuPoints);
        pGrads.resize(nQuPoints);
        solLaplacian.resize(nQuPoints);
        gsMatrix<T> physLaplacian;//, physDeriv2U;
        gsMatrix<T> resMvals;
        resMvals.setZero(nQuPoints, tarDim);
        gsVector<T> resCvals;
        resCvals.setZero(nQuPoints);
        gsVector<T> nuT(nQuPoints);
        // Evaluate vorticity at pts
        for (int k = 0; k < nQuPoints; k++)
        {
            // eval uGrads at all pts
            bases.front().active_into(pts.col(k), activesU);
            int numActU = activesU.rows();
            actCoeffsU.setZero(tarDim, numActU);
            bases.back().active_into(pts.col(k), activesP);
            int numActP = activesP.rows();
            actCoeffsP.setZero(1, numActP);

            for (int j = 0; j < numActU; j++)
                actCoeffsU.col(j) = uSolField.coefficientVector(patchIndex).row(activesU(j)).transpose();
            for (int j = 0; j < numActP; j++)
                actCoeffsP.col(j) = pSolField.coefficientVector(patchIndex).row(activesP(j)).transpose();

            geoEval->transformGradients(k, basisGradsU, physGradU);
            uGrads[k].noalias() = actCoeffsU * physGradU.transpose();
            geoEval->transformGradients(k, basisGradsP, physGradP);
            pGrads[k].noalias() = actCoeffsP * physGradP.transpose();
            geoEval->transformLaplaceHgrad(k, basisGradsU, bHessianU, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
            solLaplacian[k].noalias() = actCoeffsU * physLaplacian.transpose();
            //geoEval->transformDeriv2Hgrad(k, basisGradsU, bHessianU, physDeriv2U);
        }

        m_pTurbulenceSolver->evalTurbViscosity_into(patchIndex, pts, TMsolution, uGrads, nuT);

        for (int k = 0; k < nQuPoints; k++)
        {
            for (int var = 0; var < tarDim; var++)
            {
                resMvals(k, var) = 1./m_timeStepSize * (solUVals(var, k) - solUValsOld(var, k)) //time derivative
                                 + uGrads[k].row(var) * solUVals.col(k) //advection
                                 - (viscosity + nuT(k)) * solLaplacian[k](var, 0) //diffusion
                                 + pGrads[k](0, var); //pressure term

                m_averageResidual(var) += resMvals(k, var);
                m_averageAbsResidual(var) += math::abs(resMvals(k, var));

                /*if (tanh)
                    resMvals(k, var) = math::pow(math::tanh(resMvals(k, var)), 2);
                else
                {
                    T L_step = 0.0127;//8. * 0.0127;
                    T U = 44.2;
                    T scaleFactor = L_step / math::pow(U, 2);
                    T K = 1.;//100.;
                    if (scale)
                        resMvals(k, var) = K * scaleFactor * resMvals(k, var);
                    else if (tanhScale)
                        resMvals(k, var) = math::pow(math::tanh(K * scaleFactor * resMvals(k, var)), 2);
                }*/

                resCvals(k) += uGrads[k](var, var);
            }
        }

        residuumMvals = resMvals.transpose();
        residuumCvals = resCvals.transpose();
    }

    void evalElWiseForLocRef(const gsMatrix<T>& solVector, const gsMatrix<T>& oldSolVector, std::vector<gsVector<T> > & elWiseVals, bool outputInQuadPoints = false)
    {
        gsInfo << "Evaluating elWiseVals of RANS problem for local refinement...\n";
        m_blockAssembler.evaluateRANSLocRefCritElWiseVal(solVector, oldSolVector, elWiseVals, outputInQuadPoints);
        gsInfo << "Evaluating done.\n";
    }

    void evalElWiseForLocRef(const gsMatrix<T>& solVector, std::vector<gsVector<T> > & elWiseVals, bool outputInQuadPoints = false)
    {
        gsInfo << "Evaluating elWiseVals of RANS problem for local refinement...\n";
        m_blockAssembler.evaluateRANSLocRefCritElWiseVal(solVector, this->getSolution(), elWiseVals, outputInQuadPoints);
        gsInfo << "Evaluating done.\n";
    }

protected:
    uwbTMSolverBase<T>* m_pTurbulenceSolver;

    gsSparseMatrix<T> m_StokesMatrix;
    gsMatrix<T> m_StokesRhs;

    gsSparseMatrix<T> m_NMatrix;
    gsMatrix<T> m_NRhs;

    gsVector<T> m_averageResidual;
    gsVector<T> m_averageAbsResidual;

    // members from uwbINSAssemblerBase
    using uwbINSAssemblerBase<T>::m_blockAssembler;
    using uwbINSAssemblerBase<T>::m_solution;
    using uwbINSAssemblerBase<T>::m_bInitialized;

    // members from uwbINSAssemblerUnsteady
    using Base::m_baseMatrix;
    using Base::m_baseRhs;
    using Base::m_timeStepSize;
    

}; // class uwbRANSAssembler

} // namespace gismo
