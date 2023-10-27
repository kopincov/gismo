/** @file gsBramblePascialCG.hpp

    @brief General Conjugate gradient solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#include <gsSolver/gsBramblePasciakCG.h>
#include <gsSolver/gsLanczosMatrix.h>

namespace gismo
{

namespace internal
{

template<typename VectorType>
void gsBPCG_Preconditioner<gsBPCG_Types::BP,VectorType>::apply( VectorType& input, VectorType& result, VectorType& hResult, gsBramblePasciakCG<gsBPCG_Types::BP,typename VectorType::Scalar>& bpcg )
{
    VectorType tempA_(bpcg.m_nA,1);
    VectorType tempC_(bpcg.m_nC,1);

    bpcg.m_precondA->apply(input.topRows(bpcg.m_nA),tempA_);
    tempA_/=bpcg.m_precond_scaling[0];
    result.topRows(bpcg.m_nA) = tempA_;
    bpcg.m_matB->apply(tempA_,tempC_);
    result.bottomRows(bpcg.m_nC)= tempC_ +input.bottomRows(bpcg.m_nC)*(1-2*bpcg.m_scaling[0]);

    bpcg.m_matA->apply(result.topRows(bpcg.m_nA),tempA_);
    hResult.topRows(bpcg.m_nA) = tempA_ +(1-2*bpcg.m_scaling[0])*input.topRows(bpcg.m_nA);
    hResult.bottomRows(bpcg.m_nC) = result.bottomRows(bpcg.m_nC);
}

template<typename VectorType>
void gsBPCG_Preconditioner<gsBPCG_Types::BP_Schur,VectorType>::apply( VectorType& input, VectorType& result, VectorType& hResult, gsBramblePasciakCG<gsBPCG_Types::BP_Schur,typename VectorType::Scalar>& bpcg )
{
    VectorType tempA_(bpcg.m_nA,1);
    VectorType tempC_(bpcg.m_nC,1);

    bpcg.m_precondA->apply(input.topRows(bpcg.m_nA),tempA_);
    tempA_/=bpcg.m_precond_scaling[0];
    result.topRows(bpcg.m_nA) = tempA_;
    bpcg.m_matB->apply(tempA_,tempC_);

    bpcg.m_matA->apply(result.topRows(bpcg.m_nA),tempA_);
    hResult.topRows(bpcg.m_nA) = tempA_ +(1-2*bpcg.m_scaling[0])*input.topRows(bpcg.m_nA);
    hResult.bottomRows(bpcg.m_nC) = tempC_ +input.bottomRows(bpcg.m_nC)*(1-2*bpcg.m_scaling[0]);

    bpcg.m_precondC->apply(tempC_,tempA_);
    result.bottomRows(bpcg.m_nC)= tempA_/bpcg.m_precond_scaling[1];
    bpcg.m_precondC->apply(input.bottomRows(bpcg.m_nC),tempA_);
    result.bottomRows(bpcg.m_nC)+= (1-2*bpcg.m_scaling[0])*tempA_/bpcg.m_precond_scaling[1];
}

template<typename VectorType>
void gsBPCG_Preconditioner<gsBPCG_Types::SchoeberlZulehner,VectorType>::apply( VectorType& input, VectorType& result, VectorType& hResult, gsBramblePasciakCG<gsBPCG_Types::SchoeberlZulehner,typename VectorType::Scalar>& bpcg )
{
    VectorType tempA_, tempA_2(bpcg.m_nA,1);
    VectorType tempC_, tempC_2(bpcg.m_nC,1);

    bpcg.m_precondA->apply(input.topRows(bpcg.m_nA),tempA_);
    tempA_/=bpcg.m_precond_scaling[0];

    bpcg.m_matB->apply(tempA_,tempC_);
    bpcg.m_precondC->apply(tempC_ -input.bottomRows(bpcg.m_nC),tempC_2);
    result.bottomRows(bpcg.m_nC)=tempC_2 / bpcg.m_precond_scaling[1];
    bpcg.m_matBT->apply(result.bottomRows(bpcg.m_nC),tempA_);
    bpcg.m_precondA->apply(input.topRows(bpcg.m_nA)-tempA_,tempA_2);
    result.topRows(bpcg.m_nA)=tempA_2/bpcg.m_precond_scaling[0];

    bpcg.m_mat->apply(-result,hResult);
    hResult+=input;
}

}// namespace internal

template<gsBPCG_Types::Type type, typename T>
bool gsBramblePasciakCG<type,T>::initIteration( const VectorType& rhs, VectorType& x )
{
    if (m_calcEigenvals)
    {
        m_delta.clear();
        m_delta.resize(1,0);
        m_delta.reserve(m_max_iters / 3);

        m_gamma.clear();
        m_gamma.reserve(m_max_iters / 3);
    }

    if (Base::initIteration(rhs,x))
        return true;

    int n = m_mat->cols();
    int m = 1;                                                          // == rhs.cols();
    m_tmp.resize(n,m);
    m_update.resize(n,m);
    m_z.resize(n,m);
    m_hz.resize(n,m);
    m_precTemp.resize(n,m);
    m_hPtmp.resize(n,m);

    m_mat->apply(x,m_tmp);                                              // apply the system matrix
    //VectorType res = rhs - m_tmp;
    m_residual = rhs - m_tmp;                                                  // initial residual

    m_error = m_residual.norm() / m_rhs_norm;                                  // ?????? TODO: check, whats the right error
    if (m_error < m_tol)
        return true;

    applyPreconditioner(m_residual,m_z,m_hz);

  //  m_error = m_hz.norm()/m_rhs_norm;

   // m_error = 1;
   // m_rhs_norm = sqrt(std::abs<T>(res.col(0).dot(m_hz.col(0))));

    m_update = m_z;

    m_abs_new = m_z.col(0).dot(m_hz.col(0));                             // the square of the absolute value of r scaled by invM

    return false;
}

template<gsBPCG_Types::Type type, typename T>
bool gsBramblePasciakCG<type,T>::step( VectorType& x )
{
    m_mat->apply(m_update,m_tmp);                                      // apply system matrix

    applyPreconditioner(m_tmp,m_precTemp,m_hPtmp);

  //  T gamma = m_tmp.col(0).dot(m_hPtmp.col(0));
    T alpha = m_abs_new / m_update.col(0).dot(m_hPtmp.col(0));         // the amount we travel on dir
    if (m_calcEigenvals)
        m_delta.back()+=(1./alpha);

    x += alpha * m_update;                                             // update solution
    m_z -= alpha * m_precTemp;                                         // update residual
    m_residual -= alpha*m_tmp;
    m_error = m_residual.norm()/m_rhs_norm;

  //  m_error = m_z.norm() / m_rhs_norm;

   // m_error = sqrt(std::abs<T>(gamma))/m_rhs_norm;
    if (m_error < m_tol)
        return true;

    m_hz-= alpha*m_hPtmp;

    T abs_old = m_abs_new;

    m_abs_new = m_z.col(0).dot(m_hz.col(0));                          // update the absolute value of r
    T beta = m_abs_new / abs_old;                                     // calculate the Gram-Schmidt value used to create the new search direction
    m_update = m_z + beta * m_update;                                 // update search direction

    if (m_calcEigenvals)
    {
        m_gamma.push_back(-math::sqrt(beta)/alpha);
        m_delta.push_back(beta/alpha);
    }
    return false;
}


template<gsBPCG_Types::Type type, typename T>
T gsBramblePasciakCG<type,T>::getConditionNumber()
{
    if ( m_delta.empty() )
    {
        gsWarn<< "Condition number needs eigenvalues set setCalcEigenvalues(true)"
                 " and call solve with an arbitrary right hand side";
        return -1;
    }

    gsLanczosMatrix<T> L(m_gamma,m_delta);
    return L.maxEigenvalue()/L.minEigenvalue();
}

template<gsBPCG_Types::Type type, typename T>
void gsBramblePasciakCG<type,T>::getEigenvalues( gsMatrix<T>& eigs )
{
    if ( m_delta.empty() )
    {
        gsWarn<< "Eigenvalues were not computed, set setCalcEigenvalues(true)"
                 " and call solve with an arbitrary right hand side";
        eigs.clear();
        return;
    }

   gsLanczosMatrix<T> LM(m_gamma,m_delta);
   gsSparseMatrix<T> L = LM.matrix();
   // there is probably a better option...
   typename gsMatrix<T>::SelfAdjEigenSolver eigensolver(L);
   eigs = eigensolver.eigenvalues();
}

} // end namespace gismo


