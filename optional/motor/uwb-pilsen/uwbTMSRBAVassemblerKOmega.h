/** @file uwbTMSRBAVassemblerKOmega.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbTMAssemblerBase.h"
#include "uwbTMAssemblerKOmega.h"

namespace gismo
{

template<class T>
class uwbTMSRBAVassemblerKOmega : public uwbTMAssemblerKOmega<T>
{

public:
    typedef uwbTMAssemblerKOmega<T> Base;

public:
    uwbTMSRBAVassemblerKOmega(uwbINSSolverParams<T>& params) : Base(params)
    {
    }

    virtual ~uwbTMSRBAVassemblerKOmega()
    { }

protected:

    virtual void initAssembly(const gsField<T> & uSolField)
    {
        Base::initAssembly(uSolField);
        m_blockAssembler.assembleSRBAVpatternBlock_kOmega();
    }

    virtual void updateAssembly(const gsField<T> & uSolField)
    {
        Base::updateAssembly(uSolField);
        m_blockAssembler.assembleSRBAVBlock_kOmega(m_solution);
    }

    virtual void updatePicardAssembly(const gsMatrix<T> & solVector)
    {
        Base::updatePicardAssembly(solVector);
        m_blockAssembler.assembleSRBAVBlock_kOmega(solVector, false);
    }

public:
    //===================================================== fillExplicit ================================================================
    virtual void fillExplicitPartMatrix()
    {
        Base::fillExplicitPartMatrix();

        int dofs = this->numDofs();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockSRBAV(), col); it; ++it)
                m_nExplicitPartMatrix.coeffRef(it.row(), it.col()) += it.value();
    }

    //===================================================== fillBase ================================================================

    virtual void fillBase()  //+pattern
    {
        Base::fillBase();
        
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();


        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t i = 0; i < dofs; ++i)
            nonZerosPerColumnVector(i) = m_blockAssembler.getBlockSRBAVpattern().col(i).nonZeros();

        gsSparseMatrix<T> matrixSRBAVpattern(varDofs, dofs);
        matrixSRBAVpattern.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockSRBAVpattern(), col); it; ++it)
                matrixSRBAVpattern.insert(it.row(), it.col()) = 0.;

        m_baseMatrix += matrixSRBAVpattern;
        //m_nExplicitPartMatrix += matrixSRBAVpattern;

        if (!m_baseMatrix.isCompressed())
            m_baseMatrix.makeCompressed();

        m_bSystemReady = false;

    } //end fillBase

public:
    //===================================================== fillSystem ==============================================================

    virtual void fillSystem()
    {
        Base::fillSystem();

        int dofs = this->numDofs();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockSRBAV(), col); it; ++it)
                m_matrix.coeffRef(it.row(), it.col()) += m_theta * it.value();

        m_rhs += m_blockAssembler.getRhsSRBAV();

        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();

        m_bSystemReady = true;

    } //end fillSystem

protected:
    //uwbTMStabilizationBlockAssembler<T>* m_pBlockAssembler;

    // members from uwbTMAssemblerBaseUnsteady
    using uwbTMAssemblerBaseUnsteady<T>::m_nBaseMatrix;
    using uwbTMAssemblerBaseUnsteady<T>::m_nBaseRhs;
    using uwbTMAssemblerBaseUnsteady<T>::m_nExplicitPartRhs;
    using uwbTMAssemblerBaseUnsteady<T>::m_nExplicitPartMatrix;
    using uwbTMAssemblerBaseUnsteady<T>::m_timeStepSize;
    using uwbTMAssemblerBaseUnsteady<T>::m_theta;

    // members from uwbTMAssemblerBase
    using uwbTMAssemblerBase<T>::m_blockAssembler;
    using uwbTMAssemblerBase<T>::m_baseMatrix;
    using uwbTMAssemblerBase<T>::m_matrix;
    using uwbTMAssemblerBase<T>::m_rhs;
    using uwbTMAssemblerBase<T>::m_solution;
    using uwbTMAssemblerBase<T>::m_bSystemReady;

}; //uwbTMSRBAVassemblerKOmega

} //namespace gismo
