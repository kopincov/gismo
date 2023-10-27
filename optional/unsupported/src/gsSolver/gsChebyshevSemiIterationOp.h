/** @file gsChebyshevSemiIteration.h

    @brief Chebyshev semiteration solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Chebyshev semi iteration
///
/// \ingroup Solver
template <typename T>
class gsChebyshevOp : public gsLinearOperator<T>
{
public:
    /// Shared pointer for gsRichardsonOp
    typedef memory::shared_ptr< gsChebyshevOp > Ptr;

    /// Unique pointer for gsRichardsonOp
    typedef memory::unique_ptr< gsChebyshevOp > uPtr;
 
    /// Base class
    typedef gsLinearOperator<T> Base;

    typedef typename gsLinearOperator<T>::Ptr BasePtr;

    /// @brief Constructor with given matrix
    explicit gsChebyshevOp(const BasePtr& mat, const BasePtr& precond, index_t iter, T lambda)
    : m_mat(mat), m_precond(precond), m_iter(iter), m_lambda(lambda)
    {
        GISMP_ASSERT( mat->rows() == mat->cols() && precond->rows() == precond->cols() && precond->rows() == mat->rows(),
          "Dimensions do not agree.");
    }
    
    static Ptr make(const BasePtr& mat, const BasePtr& precond, index_t iter, T lambda)
        { return Ptr(new gsChebyshevOp<T>(mat,precond,iter,lambda)); }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT( m_mat->rows() == input.rows(), "Dimensions do not match.");

        int n = m_mat->cols();
        int m = input.cols();

        gsMatrix<T> tmp(n,m);
        gsMatrix<T> update(n,m);
    
        m_mat->apply(x,tmp);
        gsMatrix<T> res = input - tmp;
    
        m_precond->apply(res, update);
        T alpha = ((T)2)/m_lambda;
        
        for (int i=0; i<m_iter; ++i)
        {
            // Update iterate
            x += alpha * update;
            
            // Update residual
            m_mat->apply(update,tmp);
            res -= alpha * tmp;

            if ( i == m_iter-1 ) return;        
                  
            // Compute the preconditioned residual
            m_precond->apply(res, tmp);
            
            // Update search direction
            real_t beta = (m_lambda*m_lambda*alpha*alpha);
            alpha = (T)1 / (m_lambda / (T)2 - beta/alpha);
            update = tmp + beta * update;
    
        }

    }

    index_t rows() const {return m_mat->rows();}
    index_t cols() const {return m_mat->cols();}


private:
    const BasePtr m_mat;
    const BasePtr m_precond;
    index_t m_iter;
    T m_lambda;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsChebyshevSemiIteration.hpp)
#endif
