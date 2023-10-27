/** @file uwbINSAssemblerUnsteady.h

Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once

#include "uwbINSAssemblerBase.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerUnsteady : public uwbINSAssemblerBase<T>
{

public:
    typedef uwbINSAssemblerBase<T> Base;

public:
    uwbINSAssemblerUnsteady(uwbINSSolverParams<T>& params) : Base(params)
    {
        m_timeStepSize = params.settings().get(constantsINS::timeStep);
        initMembers();
    }

    virtual ~uwbINSAssemblerUnsteady()
    {
    }

protected:

    void initMembers()
    {
        int dofs = Base::numDofs();

        m_baseMatrix.resize(dofs, dofs);
        m_baseRhs.setZero(dofs, 1);

        m_matrix.resize(dofs, dofs);
        m_rhs.setZero(dofs, 1);

        m_solution.setZero(dofs, 1);

        m_bInitialized = false;
        m_bMatrixReady = false;
        m_bRhsReady = false;

        m_blockAssembler.setUnsteady(true);
    }

    virtual void reinitMembers() { initMembers(); }

public:

    virtual void initialize()
    {
        m_blockAssembler.assembleLinearStokesPart();
        m_blockAssembler.assembleBlockNpattern();

        fillBase();

        if (m_blockAssembler.isCROSSWIND())
        {
            m_blockAssembler.assembleBlockCROSSWINDpattern();
            fillCROSSWINDBase();
        }

        if (m_blockAssembler.isSUPG())
        {
            m_blockAssembler.assembleBlockNonlinearSUPGpattern();
            fillSupgCsdBase();
        }

        if (m_blockAssembler.isTCSD())
        {
            m_blockAssembler.assembleBlockNonlinearCSDpattern();
            fillSupgCsdBase();
        }

        if (!m_baseMatrix.isCompressed())
            m_baseMatrix.makeCompressed();

        m_bInitialized = true;
    }

    virtual void initialize(gsMatrix<T> & initialConditionCoeffVector)
    {
        initialize();
        Base::update(initialConditionCoeffVector);
    }

    virtual void updateAssembly()
    {
        Base::updateAssembly();

        if (m_blockAssembler.isSUPG())
            m_blockAssembler.assembleNonlinearSUPGPart(m_solution);

        if (m_blockAssembler.isTCSD())
            m_blockAssembler.assembleNonlinearCSDPart(m_solution);

        if (m_blockAssembler.isCROSSWIND())
            m_blockAssembler.assembleCrosswindPart(m_solution);
        if (m_blockAssembler.isSRBAV() && !m_blockAssembler.isRANS())
            m_blockAssembler.assembleINSsrbavPart(m_solution);
    }

    virtual void updatePicardAssembly(const gsMatrix<T> & solVector)
    {
        m_blockAssembler.assembleNonlinearPart(solVector, false);

        if (m_blockAssembler.isSUPG())
            m_blockAssembler.assembleNonlinearSUPGPart(solVector, false);

        if (m_blockAssembler.isTCSD())
            m_blockAssembler.assembleNonlinearCSDPart(solVector, false);

        if (m_blockAssembler.isCROSSWIND())
            m_blockAssembler.assembleCrosswindPart(solVector, false);
        if (m_blockAssembler.isSRBAV() && !m_blockAssembler.isRANS())
            m_blockAssembler.assembleINSsrbavPart(solVector, false);
    }

    virtual void changeTimeStep(const T timeStepSize)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_timeStepSize = timeStepSize;

        m_matrix.resize(Base::numDofs(), Base::numDofs());
        m_rhs.setZero();

        fillBase();
    }

protected:

    virtual void fillBase()
    {
        const T invTimeStep = 1. / m_timeStepSize;
        int numDofs = Base::numDofs();
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        gsSparseMatrix<T> stokesMatrix(numDofs, numDofs);

        m_blockAssembler.fillStokesSystem_into(stokesMatrix, m_baseRhs);

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(numDofs);
        for (int s = 0; s < tarDim; ++s)
            for (int i = 0; i < uDofs; i++)
                nonZerosPerColumnVector(i + s * uDofs) = m_blockAssembler.getBlockNpattern().col(i).nonZeros();

        gsSparseMatrix<T> blockNpatternMatrix(numDofs, numDofs);
        blockNpatternMatrix.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockNpattern(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    blockNpatternMatrix.insert(it.row() + s * uDofs, it.col() + s * uDofs) = 0.;

        m_baseMatrix = stokesMatrix + blockNpatternMatrix;

        if (m_blockAssembler.isTimeDerTerm())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockM(), col); it; ++it)
                    for (index_t s = 0; s < tarDim; ++s)
                        m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += invTimeStep * it.value();
        }

        m_bMatrixReady = false;
        m_bRhsReady = false;
    }

    virtual void fillMatrix()
    {
        int uDofs = Base::getUdofs();

        m_matrix = m_baseMatrix;
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN(), col); it; ++it)
                for (index_t s = 0; s < Base::getTarDim(); ++s)
                    m_matrix.coeffRef(it.row() + s*uDofs, it.col() + s*uDofs) += it.value();

        if (m_blockAssembler.isSUPG() || m_blockAssembler.isTCSD())
            fillSupgCsdMatrix();
        if (m_blockAssembler.isCROSSWIND())
            fillCROSSWINDMatrix();
        if (m_blockAssembler.isSRBAV() && !m_blockAssembler.isRANS())
            fillINSsrbavMatrix();

        m_bMatrixReady = true;

        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();
    }

    void fillRhs()
    {
        const T invTimeStep = 1. / m_timeStepSize;
        int uDofs = Base::getUdofs();

        m_rhs.noalias() = m_baseRhs + m_blockAssembler.getRhsN();

        if (m_blockAssembler.isTimeDerTerm())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t s = 0; s < Base::getTarDim(); ++s)
                m_rhs.middleRows(s * uDofs, uDofs).noalias() += invTimeStep * m_blockAssembler.getBlockM() * m_blockAssembler.getSolution().middleRows(s * uDofs, uDofs);
        }

        if (m_blockAssembler.isSUPG() || m_blockAssembler.isTCSD())
            fillSupgCsdRhs();
        if (m_blockAssembler.isCROSSWIND())
            m_rhs.noalias() += m_blockAssembler.getRhsCrosswind();
        if (m_blockAssembler.isSRBAV() && !m_blockAssembler.isRANS())
            m_rhs.noalias() += m_blockAssembler.getRhsSRBAV();

        m_bRhsReady = true;
    }

    //--------------------------------------------------------------------------------------------

    void fillSupgCsdBase()
    {
        int numDofs = Base::numDofs();
        int uDofs = Base::getUdofs();
        int pDofs = Base::getPdofs();
        int tarDim = Base::getTarDim();
        int pShift = Base::getPshift();

        gsVector<int> nonZerosPerColumnVector;
        gsSparseMatrix<T> blockApatternMatrix_SUPG(numDofs, numDofs);
        gsSparseMatrix<T> blockBpatternMatrix_SUPG(numDofs, numDofs);
        if (m_blockAssembler.isSUPG())
        {
            //------ A_SUPG -----------
            nonZerosPerColumnVector.setZero(numDofs);
            for (int s = 0; s < tarDim; ++s)
                for (int i = 0; i < uDofs; i++)
                    nonZerosPerColumnVector(s * uDofs + i) = m_blockAssembler.getBlockApattern_SUPG().col(i).nonZeros();

            blockApatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockApattern_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < tarDim; ++s)
                        blockApatternMatrix_SUPG.insert(s * uDofs + it.row(), s * uDofs + it.col()) = 0.;

            //------ B_SUPG -----------
            nonZerosPerColumnVector.setZero(numDofs);
            for (int i = 0; i < pDofs; i++)
                for (int s = 0; s < tarDim; ++s)
                    nonZerosPerColumnVector(i + pShift) += m_blockAssembler.getBlockBpattern_SUPG(s).col(i).nonZeros();

            blockBpatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < pDofs; ++col)
                for (index_t s = 0; s < tarDim; ++s)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockBpattern_SUPG(s), col); it; ++it)
                        blockBpatternMatrix_SUPG.insert(it.row() + s*uDofs, it.col() + pShift) = 0.;
        }

        //------ N_SUPG -----------
        nonZerosPerColumnVector.setZero(numDofs);
        for (int s = 0; s < tarDim; ++s)
            for (int i = 0; i < uDofs; i++)
                nonZerosPerColumnVector(i + s *  uDofs) = m_blockAssembler.getBlockNpattern_SUPG().col(i).nonZeros();

        gsSparseMatrix<T> blockNpatternMatrix_SUPG(numDofs, numDofs);
        blockNpatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockNpattern_SUPG(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    blockNpatternMatrix_SUPG.insert(it.row() + s* uDofs, it.col() + s* uDofs) = 0.;

        if (m_blockAssembler.isTimeDerTerm())
        {
            //------ M_SUPG unsteady-----------
            nonZerosPerColumnVector.setZero(numDofs);
            for (int s = 0; s < tarDim; ++s)
                for (int i = 0; i < uDofs; i++)
                    nonZerosPerColumnVector(i + s *  uDofs) = m_blockAssembler.getBlockMpattern_SUPG().col(i).nonZeros();

            gsSparseMatrix<T> blockMpatternMatrix_SUPG(numDofs, numDofs);
            blockMpatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockMpattern_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < tarDim; ++s)
                        blockMpatternMatrix_SUPG.insert(it.row() + s* uDofs, it.col() + s* uDofs) = 0.;

            m_baseMatrix += blockMpatternMatrix_SUPG;
        }

        if (m_blockAssembler.isTCSD())
            m_baseMatrix += blockNpatternMatrix_SUPG;
        else
            m_baseMatrix += blockNpatternMatrix_SUPG + blockBpatternMatrix_SUPG + blockApatternMatrix_SUPG;

        if (m_blockAssembler.isRotation())
        {
            //------ M_SUPG rotation -----------
            nonZerosPerColumnVector.setZero(numDofs);
            for (int s = 0; s < tarDim; ++s)
                for (int i = 0; i < uDofs; i++)
                    nonZerosPerColumnVector(i + s *  uDofs) = m_blockAssembler.getBlockMpattern_SUPG().col(i).nonZeros();

            gsSparseMatrix<T> blockMpatternMatrix_SUPG_rot(numDofs, numDofs);
            blockMpatternMatrix_SUPG_rot.reserve(nonZerosPerColumnVector);

            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockMpattern_SUPG(), col); it; ++it)
                {
                    blockMpatternMatrix_SUPG_rot.insert(it.row() + ((tarDim - 2) * uDofs), it.col() + ((tarDim - 1) * uDofs)) = 0.;
                    blockMpatternMatrix_SUPG_rot.insert(it.row() + ((tarDim - 1) * uDofs), it.col() + ((tarDim - 2) * uDofs)) = 0.;
                }

            m_baseMatrix += blockMpatternMatrix_SUPG_rot;
        }

        m_bMatrixReady = false;
        m_bRhsReady = false;
    }

    void fillSupgCsdMatrix()
    {
        const T invTimeStep = 1. / m_timeStepSize;
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        if (m_blockAssembler.isSUPG())
        {
            //------ A_SUPG -----------
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockA_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < tarDim; ++s)
                        m_matrix.coeffRef(s * uDofs + it.row(), s * uDofs + it.col()) += it.value();

            //------ B_SUPG -----------
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < Base::getPdofs(); ++col)
                for (index_t s = 0; s < tarDim; ++s)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockB_SUPG(s), col); it; ++it)
                        m_matrix.coeffRef(it.row() + s*uDofs, it.col() + Base::getPshift()) += it.value();
        }

        //------ N_SUPG -----------
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN_SUPG(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    m_matrix.coeffRef(it.row() + s* uDofs, it.col() + s* uDofs) += it.value();

        if (m_blockAssembler.isTimeDerTerm())
        {
            //------ M_SUPG unsteady-----------
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockM_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < tarDim; ++s)
                        m_matrix.coeffRef(it.row() + s* uDofs, it.col() + s* uDofs) += invTimeStep * it.value();
        }

        if (m_blockAssembler.isRotation())
        {
            //------ M_SUPG rotation-----------
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockM_SUPG(), col); it; ++it)
                {
                    m_matrix.coeffRef(it.row() + ((tarDim - 2) * uDofs), it.col() + ((tarDim - 1) * uDofs)) = -m_blockAssembler.getOmega() * it.value();
                    m_matrix.coeffRef(it.row() + ((tarDim - 1) * uDofs), it.col() + ((tarDim - 2) * uDofs)) = m_blockAssembler.getOmega() * it.value();
                }

        }
    }

    virtual void fillSupgCsdRhs()
    {
        const T invTimeStep = 1. / m_timeStepSize;
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        if (m_blockAssembler.isSUPG())
            m_rhs.noalias() += m_blockAssembler.getRhsN_SUPG() + m_blockAssembler.getRhsB_SUPG() + m_blockAssembler.getRhsA_SUPG();
        else
            m_rhs.noalias() += m_blockAssembler.getRhsN_SUPG();

        if (m_blockAssembler.isTimeDerTerm())
        {
            //m_rhs.noalias() += invTimeStep * m_blockAssembler.getRhsM_SUPG();

            // unsteady term on the RHS
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t s = 0; s < tarDim; ++s)
                m_rhs.middleRows(s * uDofs, uDofs).noalias() += invTimeStep * m_blockAssembler.getBlockM_SUPG() * m_blockAssembler.getSolution().middleRows(s * uDofs, uDofs);
        }

        // BC for rotation term
        if (m_blockAssembler.isRotation())
        {
            m_rhs.middleRows((tarDim - 2) * uDofs, uDofs) -= m_blockAssembler.getOmega() * m_blockAssembler.getRhsM_SUPG().middleRows((tarDim - 1) * uDofs, uDofs);
            m_rhs.middleRows((tarDim - 1) * uDofs, uDofs) += m_blockAssembler.getOmega() * m_blockAssembler.getRhsM_SUPG().middleRows((tarDim - 2) * uDofs, uDofs);
        }
    }

    void fillCROSSWINDBase()
    {
        int numDofs = Base::numDofs();
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(numDofs);
        for (int i = 0; i < uDofs; i++)
            for (int s = 0; s < tarDim; ++s)
                nonZerosPerColumnVector(i + s * uDofs) += m_blockAssembler.getBlockCrosswindPattern(s).col(i).nonZeros();

        gsSparseMatrix<T> blockCROSSWINDMatrix(numDofs, numDofs);
        blockCROSSWINDMatrix.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockCrosswindPattern(s), col); it; ++it)
                    blockCROSSWINDMatrix.insert(it.row() + s * uDofs, it.col() + s * uDofs) = 0.;


        m_baseMatrix += blockCROSSWINDMatrix;

        m_bMatrixReady = false;
        m_bRhsReady = false;
    }

    void fillCROSSWINDMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockCrosswind(s), col); it; ++it)
                    m_matrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

    void fillINSsrbavMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockSRBAV(s), col); it; ++it)
                    m_matrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

    //---------------------------------------------------------------------------------------------------

public:
    void plotResiduum(std::string const & fn, gsMatrix<T> & solution, gsMatrix<T>& solutionOld, unsigned npts)
    {
        std::string fnM = fn + "MomentumEq";
        std::string fnC = fn + "ContinuityEq";

        gsParaviewCollection collectionM(fnM);
        gsParaviewCollection collectionC(fnC);
        std::string fileNameM, fileNameC;

        gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);
        gsField<T> pSolField = m_blockAssembler.constructSolution(solution, 1);
        gsField<T> uSolFieldOld = m_blockAssembler.constructSolution(solutionOld, 0);

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

            evalUnsteadyINSresiduum_singlePatch(p, pts, uSolField, pSolField, uSolFieldOld, residuumMomentumVals, residuumContinuityVals);

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
            gsWriteParaviewTPgrid(geoVals, residuumContinuityVals, np.template cast<index_t>(), fileNameC);

            collectionM.addPart(fileNameM + ".vts");
            collectionC.addPart(fileNameC + ".vts");
        }
        collectionM.save();
        collectionC.save();
    }

protected:
    void evalUnsteadyINSresiduum_singlePatch(index_t patchIndex, gsMatrix<T>& pts, gsField<T>& uSolField, gsField<T>& pSolField,
                                      gsField<T>& uSolFieldOld, gsMatrix<T>& residuumMvals, gsMatrix<T>& residuumCvals)
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
        gsMatrix<T> uGrads, pGrads, solLaplacian;
        gsMatrix<T> physLaplacian;//, physDeriv2U;
        gsMatrix<T> resMvals;
        resMvals.setZero(nQuPoints, tarDim);
        gsVector<T> resCvals;
        resCvals.setZero(nQuPoints);
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
            uGrads.noalias() = actCoeffsU * physGradU.transpose();
            geoEval->transformGradients(k, basisGradsP, physGradP);
            pGrads.noalias() = actCoeffsP * physGradP.transpose();
            geoEval->transformLaplaceHgrad(k, basisGradsU, bHessianU, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
            solLaplacian.noalias() = actCoeffsU * physLaplacian.transpose();
            //geoEval->transformDeriv2Hgrad(k, basisGradsU, bHessianU, physDeriv2U);

            for (int var = 0; var < tarDim; var++)
            {
                resMvals(k, var) = 1./m_timeStepSize * (solUVals(var, k) - solUValsOld(var, k)) //time derivative
                                 + uGrads.row(var) * solUVals.col(k) //advection
                                 - viscosity * solLaplacian(var, 0) //diffusion
                                 + pGrads(0, var); //pressure term

                resCvals(k) += uGrads(var, var);
            }
        }

        residuumMvals = resMvals.transpose();
        residuumCvals = resCvals.transpose();
    }

public:

    virtual const gsSparseMatrix<T> & matrix() const
    {
        GISMO_ASSERT(m_bMatrixReady, "Matrix not ready, update() must be called first");
        return m_matrix;
    }

    virtual const gsMatrix<T> & rhs() const
    {
        GISMO_ASSERT(m_bRhsReady, "Rhs not ready, update() must be called first");
        return m_rhs;
    }

    virtual void fillPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType, bool assembAp, bool assembFp, bool lumping)
    {
        Fp = (1. / m_timeStepSize) * m_blockAssembler.getBlockMp();

     
        Base::fillPCDblocks(Ap, Fp, bcType, assembAp, assembFp, lumping);
    }

    const T getTimeStepSize() const { return m_timeStepSize; }

    virtual void evalElWiseForLocRef(const gsMatrix<T>& solVector, const gsMatrix<T>& oldSolVector, std::vector<gsVector<T> > & elWiseVals, bool outputInQuadPoints = false)
    {
        gsInfo << "Evaluating elWiseVals of unsteady N-S problem for local refinement...\n";
        m_blockAssembler.evaluateINSunsteadyLocRefCritElWiseVal(solVector, oldSolVector, elWiseVals, outputInQuadPoints);
        gsInfo << "Evaluating done.\n";
    }

    void evalElWiseForLocRef(const gsMatrix<T>& solVector, std::vector<gsVector<T> > & elWiseVals, bool outputInQuadPoints = false)
    {
        gsInfo << "Evaluating elWiseVals of unsteady N-S problem for local refinement...\n";
        m_blockAssembler.evaluateINSunsteadyLocRefCritElWiseVal(solVector, this->getSolution(), elWiseVals, outputInQuadPoints);
        gsInfo << "Evaluating done.\n";
    }

protected:

    gsSparseMatrix<T> m_baseMatrix;
    gsMatrix<T> m_baseRhs;

    bool m_bMatrixReady;
    bool m_bRhsReady;

    T m_timeStepSize;

    // members from uwbINSAssemblerBase
    using Base::m_blockAssembler;
    using Base::m_solution;
    using Base::m_bInitialized;

    using Base::m_matrix;
    using Base::m_rhs;

}; // class uwbINSAssemblerUnsteady

} // namespace gismo
