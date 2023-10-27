/** @file uwbINSAssemblerSteady_AFC.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbINSAssemblerBase.h"
#include "uwbINSAssemblerSteady.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerSteady_AFC : public uwbINSAssemblerSteady<T>
{

public:
    typedef uwbINSAssemblerSteady<T> Base;

public:
    uwbINSAssemblerSteady_AFC(uwbINSSolverParams<T>& params) : Base(params)
    {
        m_bAFC_HO = params.settings().get(constantsINS::AFC_HO);

        uwbINSSolverParams<T> params1(params);
        params1.getAssemblerOptions().dirStrategy = dirichlet::elimination;
        m_pBlockAssemblerBC = new uwbINSBlockAssembler<T>(params1);
    }

    virtual ~uwbINSAssemblerSteady_AFC()
    {
        if (m_pBlockAssemblerBC)
        {
            delete m_pBlockAssemblerBC;
            m_pBlockAssemblerBC = NULL;
        }
    }

public:

    virtual void updateAssembly()
    {
        m_solPicard = m_solution;
        Base::updateAssembly();
    }

protected:

    virtual void fillBase()
    {
        int numDofs = this->numDofs();
        int uDofs = this->getUdofs();
        int tarDim = this->getTarDim();

        gsSparseMatrix<T> stokesMatrix(numDofs, numDofs);

        m_blockAssembler.fillStokesSystem_into(stokesMatrix, m_baseRhs);

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(numDofs);
        //-------------------- AFC LO applied to N block -----------------------------------
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

        //--------------------- AFC LO applied to N and A blocks ---------------------------
        /*for (int s = 0; s < tarDim; ++s)
            for (int i = 0; i < uDofs; i++)
                nonZerosPerColumnVector(i + s * uDofs) = m_blockAssembler.getBlockNpattern().col(i).nonZeros() + m_blockAssembler.getBlockA().col(i).nonZeros();

        gsSparseMatrix<T> blockNpatternMatrix(numDofs, numDofs);
        blockNpatternMatrix.reserve(nonZerosPerColumnVector);

        gsSparseMatrix<T> blockAN = m_blockAssembler.getBlockNpattern()+m_blockAssembler.getBlockA();
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(blockAN, col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    blockNpatternMatrix.insert(it.row() + s * uDofs, it.col() + s * uDofs) = 0.;*/

        //------------------------------------------------

        m_baseMatrix = stokesMatrix + blockNpatternMatrix;

        m_bMatrixReady = false;
        m_bRhsReady = false;
    }

    //====================================================================================================================

    virtual void fillMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = this->getTarDim();

        m_matrix = m_baseMatrix;

        //--------------------------------- block N AFC modified -------------------------------------
        const gsSparseMatrix<T>& blockNrefTemp = m_blockAssembler.getBlockN();
        //blockNrefTemp = blockNrefTemp + m_blockAssembler.getBlockA();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
        {
            for (index_t s = 0; s < tarDim; ++s)
                m_matrix.coeffRef(col + s*uDofs, col + s*uDofs) += blockNrefTemp.at(col, col);

            for (typename gsSparseMatrix<T>::InnerIterator it(blockNrefTemp, col); it; ++it)
            {
                int row = it.row();

                if (col > row)
                {
                    T temp1 = blockNrefTemp.at(row, col);//it.value(); // N(row, col)
                    T temp2 = blockNrefTemp.at(col, row);
                    T delta = math::max(0.0, math::max(temp1, temp2));

                    for (index_t s = 0; s < tarDim; ++s)
                    {
                        m_matrix.coeffRef(row + s*uDofs, col + s*uDofs) += (temp1 - delta);
                        m_matrix.coeffRef(col + s*uDofs, row + s*uDofs) += (temp2 - delta);
                        m_matrix.coeffRef(row + s*uDofs, row + s*uDofs) += delta;//blockNrefTemp.at(row, row) + delta;
                        m_matrix.coeffRef(col + s*uDofs, col + s*uDofs) += delta;//blockNrefTemp.at(col, col) + delta;
                    }
                }
            }
        }

        m_bMatrixReady = true;

        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();
    }

    void fillRhs()
    {
        //int numDofs = this->numDofs();

        //m_rhs.setZero(numDofs, 1);
        m_rhs.noalias() = m_baseRhs + m_blockAssembler.getRhsN();
        //m_rhs.noalias() = m_baseRhs + m_blockAssembler.getRhsN() + m_blockAssembler.getRhsA();

        GISMO_ASSERT(m_bAFC_HO == false, "AFC correction implementation not finished for steady NS.");
        /*if (m_bAFC_HO)
        {
            gsMatrix<T> afcFlux(numDofs, 1);
            fillAfcFluxes(afcFlux);
            m_rhs += afcFlux;
        }*/

        applyDirichletConditionsAFC();

        m_bRhsReady = true;
    }

protected:
    void fillAfcFluxes(gsMatrix<T>& ff)
    {
        int uDofs = this->getUdofs();
        int tarDim = this->getTarDim();
        int pShift = m_blockAssembler.getPshift(); //= tarDim * uDofs
        int dofs = this->numDofs();

        const gsSparseMatrix<T>& blockMref = m_blockAssembler.getBlockM();
        const gsSparseMatrix<T>& blockNref = m_blockAssembler.getBlockN();

        gsSparseMatrix<T> flux;
        flux.resize(uDofs, pShift);
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t s = 0; s < tarDim; ++s)
            flux.middleCols(s * uDofs, uDofs) = blockMref + blockNref;

        gsMatrix<T> Pp, Pm, Qp, Qm;
        Pp.setZero(pShift, 1);
        Pm.setZero(pShift, 1);
        Qp.setZero(pShift, 1);
        Qm.setZero(pShift, 1);

        for (int j = 0; j < uDofs; j++)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(flux, j); it; ++it)
            {
                int k = it.row();
                int l = it.col();

                T temp1 = blockNref.at(k, l);
                T temp2 = blockNref.at(l, k);//puvodne toto: blockNref.at(k, l);
                T delta = math::max(0.0, math::max(temp1, temp2));

                for (int s = 0; s < tarDim; s++)
                {
                    flux.coeffRef(k, l + s*uDofs) = delta * (m_solPicard(k + s*uDofs, 0) - m_solPicard(l + s*uDofs, 0));

                    if (k != l)
                    {
                        Pp.coeffRef(k + s*uDofs, 0) += math::max(0.0, flux.at(k, l + s*uDofs));
                        Pm.coeffRef(k + s*uDofs, 0) += math::min(0.0, flux.at(k, l + s*uDofs));

                        Qp.coeffRef(k + s*uDofs, 0) += math::max(0.0, -flux.at(k, l + s*uDofs));
                        Qm.coeffRef(k + s*uDofs, 0) += math::min(0.0, -flux.at(k, l + s*uDofs));
                        Qp.coeffRef(l + s*uDofs, 0) += math::max(0.0, flux.at(k, l + s*uDofs));
                        Qm.coeffRef(l + s*uDofs, 0) += math::min(0.0, flux.at(k, l + s*uDofs));
                        //Qp.coeffRef(k + s*uDofs, 0) = math::max(Qp.coeffRef(k + s*uDofs, 0), m_solPicard(l + s*uDofs, 0) - m_solPicard(k + s*uDofs, 0));
                        //Qm.coeffRef(k + s*uDofs, 0) = math::min(Qm.coeffRef(k + s*uDofs, 0), m_solPicard(l + s*uDofs, 0) - m_solPicard(k + s*uDofs, 0));
                    }
                }
            }
        }

        gsMatrix<T> Rp, Rm;
        Rp.setZero(pShift, 1);
        Rm.setZero(pShift, 1);

        for (int k = 0; k < uDofs; k++)
        {
            for (int s = 0; s < tarDim; s++)
            {
                if (Pp(k + s*uDofs, 0) == 0)
                    Rp.coeffRef(k + s*uDofs, 0) = 1.0;
                else
                    Rp.coeffRef(k + s*uDofs, 0) = math::min(1.0, Qp(k + s*uDofs, 0) / Pp(k + s*uDofs, 0));
                if (Pm(k + s*uDofs, 0) == 0)
                    Rm.coeffRef(k + s*uDofs, 0) = 1.0;
                else
                    Rm.coeffRef(k + s*uDofs, 0) = math::min(1.0, Qm(k + s*uDofs, 0) / Pm(k + s*uDofs, 0));
            }
        }

        ff.setZero(dofs, 1);
        for (int j = 0; j < uDofs; j++)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(flux, j); it; ++it)
            {
                int k = it.row();
                int l = it.col();

                for (int s = 0; s < tarDim; s++)
                {
                    if (flux.at(k, l + s*uDofs) > 0.0)
                        ff.coeffRef(k + s*uDofs, 0) += Rp(k + s*uDofs, 0) * flux.at(k, l + s*uDofs);
                    //else if (flux.at(k, l + s*uDofs) < 0.0)
                    else
                        ff.coeffRef(k + s*uDofs, 0) += Rm(k + s*uDofs, 0) * flux.at(k, l + s*uDofs);
                }
            }
        }
    }

    void applyDirichletConditionsAFC()
    {
        int unk = 0; // velocity
        //gsInfo << "applyDirichletConditionsAFC INS_AFC\n\n\n";
        int uDofs = this->getUdofs();
        //int numDofs = this->numDofs();
        int tarDim = this->getTarDim();

        const std::vector<gsMatrix<T> > & ddofs = m_pBlockAssemblerBC->getDirichletDofs();
        const gsMultiBasis<T>& basis = m_pBlockAssemblerBC->getBases().at(0);

        for (typename gsBoundaryConditions<>::const_iterator it = m_pBlockAssemblerBC->getBCs().dirichletBegin(); it != m_pBlockAssemblerBC->getBCs().dirichletEnd(); ++it)
        {
            if (it->unknown() != unk)
                continue;

            gsMatrix<index_t> boundary = basis.piece(it->patch()).boundary(it->side());

            gsSparseMatrix<T> mat(m_matrix.rows(), m_matrix.cols());
            mat.setIdentity();

            for (index_t k = 0; k != boundary.rows(); ++k)
            {
                if (m_blockAssembler.getMappers().front().is_free(boundary.at(k), it->patch())) // DoF value is in the solVector
                {
                    index_t ii = m_blockAssembler.getMappers().front().index(boundary.at(k), it->patch());
                    for (index_t s = 0; s < tarDim; ++s)
                        mat.coeffRef(ii + s*uDofs, ii + s*uDofs) = 0.;
                }
             }

            mat.prune(0,0);

            m_matrix = mat * m_matrix;

            for (index_t k = 0; k != boundary.rows(); ++k)
            {
                if (m_blockAssembler.getMappers().front().is_free(boundary.at(k), it->patch())) // DoF value is in the solVector
                {
                    index_t ii = m_blockAssembler.getMappers().front().index(boundary.at(k), it->patch());
                    index_t bi = m_pBlockAssemblerBC->getMappers().front().bindex(boundary.at(k), it->patch());

                    for (index_t s = 0; s < tarDim; ++s)
                    {
                        m_matrix.coeffRef(ii + s*uDofs, ii + s*uDofs) = 1.0;
                        m_rhs(ii + s*uDofs, 0) = ddofs[unk](bi, s);
                    }
                }
            }
        }
    }

public:
    gsMatrix<T> getSolution_full(const gsMatrix<T>& solVector) const
    {
        //only BC for velocity is expected here!!!
        int unk = 0; // velocity

        int tarDim = this->getTarDim();

        const gsDofMapper& mapper = m_blockAssembler.getMappers().front();
        const gsDofMapper& mapperBC = m_pBlockAssemblerBC->getMappers().front();
        const std::vector<gsMatrix<T> > & ddofs = m_pBlockAssemblerBC->getDirichletDofs();

        const index_t usz = mapper.freeSize();
        const index_t uszBC = mapperBC.freeSize();

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

                for (index_t s = 0; s < tarDim; ++s)
                {
                    if (mapperBC.is_free_index(ii))
                    {
                        solFull(jj + s*usz, 0) = solVector(ii + s*uszBC, 0);
                    }
                    else
                    {
                        solFull(jj + s*usz, 0) = ddofs[unk](mapperBC.global_to_bindex(ii), s);
                    }
                }
            }

            const gsDofMapper& mapperp = m_blockAssembler.getMappers().back();
            const int szp = m_blockAssembler.getBases().at(1).piece(p).size();
            for (index_t i = 0; i < szp; ++i)
            {
                index_t jj = mapperp.index(i, p);
                solFull(jj + tarDim*usz, 0) = solVector(jj + tarDim*uszBC, 0);
            }
        }

        return solFull;
    }
    //---------------------------------------------------------------------------------------------------

protected:
    uwbINSBlockAssembler<T>* m_pBlockAssemblerBC;
    bool m_bAFC_HO; //high order ... corrected LO
    gsMatrix<T> m_solPicard;

    gsSparseMatrix<T> matPom;
    gsMatrix<T> rhsPom;

    // members from uwbINSAssemblerUnsteady
    using Base::m_baseMatrix;
    using Base::m_baseRhs;

    using Base::m_bMatrixReady;
    using Base::m_bRhsReady;

    // members from uwbINSAssemblerBase
    using uwbINSAssemblerBase<T>::m_blockAssembler;
    using uwbINSAssemblerBase<T>::m_solution;

    using uwbINSAssemblerBase<T>::m_matrix;
    using uwbINSAssemblerBase<T>::m_rhs;

}; // class uwbINSAssemblerSteady_AFC

} // namespace gismo
