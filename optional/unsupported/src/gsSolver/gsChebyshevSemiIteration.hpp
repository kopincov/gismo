/** @file gsChebyshevSemiIteration.hpp

    @brief Chebyshev semiteration solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/
#include <gsSolver/gsChebyshevSemiIteration.h>
#include <gsSolver/gsSolverUtils.h>

namespace gismo
{

template<class T>    
bool gsChebyshevSemiIteration<T>::initIteration( const typename gsChebyshevSemiIteration<T>::VectorType& rhs,
                                                 typename gsChebyshevSemiIteration<T>::VectorType& x )
{
    if (Base::initIteration(rhs,x))
        return true;

    m_d = (m_lambda_max+m_lambda_min) / (T)2;
    m_c = (m_lambda_max-m_lambda_min) / (T)2;    
    
    int n = m_mat->cols();
    int m = 1;                                                          // == rhs.cols();
    m_tmp.resize(n,m);
    m_update.resize(n,m);

    m_mat->apply(x,m_tmp);                                              // apply the system matrix
    m_res = rhs - m_tmp;                                                // initial residual

    m_error = m_res.norm() / m_rhs_norm;
    if (m_error < m_tol)
        return true;

    m_precond->apply(m_res, m_update);
    m_alpha=1/m_d;
    
    return false;
}

template<class T>
bool gsChebyshevSemiIteration<T>::step( typename gsChebyshevSemiIteration<T>::VectorType& x )
{
    // Update iterate
    x += m_alpha * m_update;
    
    // Update residual
    m_mat->apply(m_update,m_tmp);
    m_res -= m_alpha * m_tmp;

    m_error = m_res.norm() / m_rhs_norm;
    if (m_error < m_tol)
        return true;
          
    // Compute the preconditioned residual
    m_precond->apply(m_res, m_tmp);
    
    // Update search direction
    T beta = (m_c*m_c*m_alpha*m_alpha) / (T)4;
    m_alpha = (T)1 / (m_d - beta/m_alpha);
    m_update = m_tmp + beta * m_update;

    return false;
}


} // end namespace gismo


