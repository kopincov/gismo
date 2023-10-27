/** @file uwbRANSAssembler_AFC.h

Author(s): E. Turnerova
*/

#pragma once

//#include "uwbINSAssemblerUnsteady.h"
#include "uwbRANSAssembler.h"
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

template<class T>
class uwbRANSAssembler_AFC : public uwbRANSAssembler<T>
{

public:
    typedef uwbRANSAssembler<T> Base;

public:
    uwbRANSAssembler_AFC(uwbINSSolverParams<T>& params, uwbTMSolverBase<T>* pTMsolver) :
        Base(params, pTMsolver)
    {
        GISMO_ASSERT(params.settings().get(constantsINS::timeDerTerm) == true, "AFC method implemented only for unsteady case.");

        m_bAFC_HO = params.settings().get(constantsINS::AFC_HO);

        uwbINSSolverParams<T> params1(params);
        params1.getAssemblerOptions().dirStrategy = dirichlet::elimination;
        m_pBlockAssemblerBC = new uwbINSBlockAssembler<T>(params1);
    }

    virtual ~uwbRANSAssembler_AFC()
    {
        if (m_pBlockAssemblerBC)
        {
            delete m_pBlockAssemblerBC;
            m_pBlockAssemblerBC = NULL;
        }
    }

public:

    virtual void updatePicard(const gsMatrix<T> & solVector)
    {
        m_solPicard = solVector;
        Base::updatePicard(solVector);
    }

    //=========================================================================================================
protected:

    virtual void fillBase()
    {
        const T invTimeStep = 1. / m_timeStepSize;
        int numDofs = this->numDofs();
        int uDofs = this->getUdofs();
        int tarDim = this->getTarDim();

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
            // lumped mass matrix M
            m_lumpedM.resize(numDofs, numDofs);
            fillLumpedMassMatrix();
            m_baseMatrix += invTimeStep * m_lumpedM;
        }

        m_bMatrixReady = false;
        m_bRhsReady = false;

        m_StokesMatrix = m_baseMatrix;
        m_StokesRhs = m_baseRhs;
    }

    void fillLumpedMassMatrix()
    {
        int uDofs = this->getUdofs();
        int tarDim = this->getTarDim();

        //int varDofs = this->numVarDofs();
        //int numVar = this->getNumVar();

        const gsSparseMatrix<T>& blockMref = m_blockAssembler.getBlockM();

        m_lumpedM.reserve(gsVector<int>::Constant(m_lumpedM.cols(), 1));

        for (int j = 0; j < uDofs; j++)
            for (int s = 0; s < tarDim; s++)
                m_lumpedM.insert(j + s*uDofs, j + s*uDofs) = blockMref.at(j, j);

        for (int j = 0; j < uDofs; j++)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(blockMref, j); it; ++it)
            {
                int k = it.row();
                int l = it.col();

                if ((blockMref.at(k, l) != 0) && k != l)
                {
                    for (int s = 0; s < tarDim; s++)
                        m_lumpedM.coeffRef(k + s*uDofs, k + s*uDofs) += blockMref.at(k, l);
                }
            }
        }
    }

    //====================================================================================================================
    /*virtual void fillNMatrix()
    {
        m_baseMatrix = m_StokesMatrix;
        m_baseRhs = m_StokesRhs;

        //int udofs = m_blockAssembler.getUdofs();

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

        m_baseRhs += m_blockAssembler.getRhsA_RANS() + m_blockAssembler.getRhsEdiag_RANS() + m_blockAssembler.getRhsEnondiag_RANS();
        m_baseRhs += m_blockAssembler.getRhsEdiag_RANS();

        m_NMatrix = m_baseMatrix;
        m_NRhs = m_baseRhs;
    }*/

    //====================================================================================================================

    virtual void fillMatrix()
    {
        //gsInfo << "fillMatrix RANS_AFC\n";
        int uDofs = this->getUdofs();
        int tarDim = this->getTarDim();

        m_matrix = m_baseMatrix;

        //--------------------------------- block N AFC modified -------------------------------------
        //const gsSparseMatrix<T>& blockNref = m_blockAssembler.getBlockN();
        const gsSparseMatrix<T>& blockNrefTemp = m_blockAssembler.getBlockN();
        //const gsSparseMatrix<T>& blockNref = m_blockAssembler.getBlockN() + m_blockAssembler.getBlockA() + m_blockAssembler.getBlockA_RANS();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
        {
            for (index_t s = 0; s < tarDim; ++s)
            {
                //const gsSparseMatrix<T>& blockNrefTemp = blockNref + m_blockAssembler.getBlockEdiag_RANS(s);
                m_matrix.coeffRef(col + s*uDofs, col + s*uDofs) += blockNrefTemp.at(col, col);
            }

            for (typename gsSparseMatrix<T>::InnerIterator it(blockNrefTemp, col); it; ++it)
            {
                int row = it.row();

                if (col > row)
                {
                    T temp1 = it.value(); // N(row, col)
                    T temp2 = blockNrefTemp.at(col, row);
                    T delta = math::max(0.0, math::max(temp1, temp2));

                    for (index_t s = 0; s < tarDim; ++s)
                    {
                        m_matrix.coeffRef(row + s*uDofs, col + s*uDofs) += (temp1 - delta);
                        m_matrix.coeffRef(col + s*uDofs, row + s*uDofs) += (temp2 - delta);
                        m_matrix.coeffRef(row + s*uDofs, row + s*uDofs) += delta;
                        m_matrix.coeffRef(col + s*uDofs, col + s*uDofs) += delta;
                    }
                }
            }
            //}
        }
        //---------------------------------------------------------------------------------

        m_bMatrixReady = true;

        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();
    }

    void fillRhs()
    {
        //gsInfo << "fillRhs RANS_AFC\n";
        const T invTimeStep = 1. / m_timeStepSize;

        m_rhs.noalias() = m_baseRhs + m_blockAssembler.getRhsN();// + m_blockAssembler.getRhsA() + m_blockAssembler.getRhsA_RANS();

        if (m_blockAssembler.isTimeDerTerm())
            m_rhs.noalias() += invTimeStep * m_lumpedM * m_solution;
        /*{
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t s = 0; s < this->getTarDim(); ++s)
                m_rhs.middleRows(s * uDofs, uDofs).noalias() += invTimeStep * m_blockAssembler.getBlockM() * m_blockAssembler.getSolution().middleRows(s * uDofs, uDofs);
        }*/

        GISMO_ASSERT(m_bAFC_HO == false, "AFC correction implementation not finished for RANS.");
        /*if (m_bAFC_HO)
        {
            gsMatrix<T> afcFlux(this->getUdofs(), this->getNumVar());
            fillAfcFluxes(afcFlux);
            m_rhs += afcFlux;
        }*/

        applyDirichletConditionsAFC();

        m_bRhsReady = true;
    }

    void applyDirichletConditionsAFC()
    {
        int unk = 0; // velocity
        //gsInfo << "applyDirichletConditionsAFC RANS_AFC\n\n\n";
        int uDofs = this->getUdofs();
        int tarDim = this->getTarDim();

        const std::vector<gsMatrix<T> > & ddofs = m_pBlockAssemblerBC->getDirichletDofs();
        const gsMultiBasis<T>& basis = m_pBlockAssemblerBC->getBases().at(0);

        for (typename gsBoundaryConditions<>::const_iterator it = m_pBlockAssemblerBC->getBCs().dirichletBegin(); it != m_pBlockAssemblerBC->getBCs().dirichletEnd(); ++it)
        {
            gsMatrix<index_t> boundary = basis.piece(it->patch()).boundary(it->side());

            for (index_t k = 0; k != boundary.rows(); ++k)
            {
                // dof number with dirichlet::none
                index_t ii = m_blockAssembler.getMappers().front().index(boundary.at(k), it->patch());

                for (int col = 0; col < uDofs; ++col)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_matrix, col); it; ++it)
                    {
                        if (it.row() == ii)
                        {
                            for (index_t s = 0; s < tarDim; ++s)
                                m_matrix.coeffRef(ii + s*uDofs, col + s*uDofs) = 0.0;
                        }
                    }

                index_t bi = m_pBlockAssemblerBC->getMappers().front().bindex(boundary.at(k), it->patch());

                for (index_t s = 0; s < tarDim; ++s)
                {
                    m_matrix.coeffRef(ii + s*uDofs, ii + s*uDofs) = 1.0;
                    m_rhs(ii + s*uDofs, 0) = ddofs[unk](bi, s);
                }
            }
        }
    }

public:
    gsMatrix<T> getSolution_full(const gsMatrix<T>& solVector) const
    {
        //only BC for velocity is expected here!!!
        int unk = 0; // velocity

        const gsDofMapper& mapper = m_blockAssembler.getMappers().front();
        const gsDofMapper& mapperBC = m_pBlockAssemblerBC->getMappers().front();
        const std::vector<gsMatrix<T> > & ddofs = m_pBlockAssemblerBC->getDirichletDofs();

        //index_t totalDofs = mapper.size();
        gsMatrix<T> solFull = m_blockAssembler.getSolution();
        //solFull.setZero(totalDofs, 1);

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            const int sz = m_blockAssembler.getBases().at(0).piece(p).size();
            for (index_t i = 0; i < sz; ++i)
            {
                index_t ii = mapperBC.index(i, p);
                index_t jj = mapper.index(i, p);

                if (mapper.is_free_index(ii))
                {
                    solFull(jj, 0) = solVector(ii, 0);
                }
                else
                {
                    solFull(jj, 0) = ddofs[unk](mapperBC.global_to_bindex(ii), 0);
                }
            }
        }

        return solFull;
    }

    //================================================================================================================================

protected:
    uwbINSBlockAssembler<T>* m_pBlockAssemblerBC;
    bool m_bAFC_HO; //high order ... corrected LO
    gsMatrix<T> m_solPicard;
    gsSparseMatrix<T> m_lumpedM;

    // members from uwbRANSAssembler
    using Base::m_StokesMatrix;
    using Base::m_StokesRhs;
    using Base::m_NMatrix;
    using Base::m_NRhs;

    // members from uwbINSAssemblerBase
    using uwbINSAssemblerBase<T>::m_blockAssembler;
    using uwbINSAssemblerBase<T>::m_solution;
    using uwbINSAssemblerBase<T>::m_bInitialized;
    using uwbINSAssemblerBase<T>::m_matrix;
    using uwbINSAssemblerBase<T>::m_rhs;

    // members from uwbINSAssemblerUnsteady
    using uwbINSAssemblerUnsteady<T>::m_baseMatrix;
    using uwbINSAssemblerUnsteady<T>::m_baseRhs;
    using uwbINSAssemblerUnsteady<T>::m_timeStepSize;
    using uwbINSAssemblerUnsteady<T>::m_bMatrixReady;
    using uwbINSAssemblerUnsteady<T>::m_bRhsReady;
    

}; // class uwbRANSAssembler_AFC

} // namespace gismo
