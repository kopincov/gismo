/** @file gsChebyshevSemiIteration.h

    @brief Chebyshev semiteration solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

/// @brief The Chebyshev semiteration method.
///
/// This is similar to conjugate gradient, but you have to know the eigenvalues
/// of the (preconditioned) system.
///
/// \ingroup Solver
template< class T = real_t >
class GISMO_EXPORT gsChebyshevSemiIteration : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;
    
    typedef gsMatrix<T>  VectorType;
    
    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsChebyshevSemiIteration> Ptr;
    typedef memory::unique_ptr<gsChebyshevSemiIteration> uPtr;
    
    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat        The operator to be solved for, see gsIterativeSolver for details
    /// @param lambda_min The smallest eigenvalue of precond*mat
    /// @param lambda_max The largest eigenvalue of precond*mat
    /// @param precond    The preconditioner, defaulted to the identity
    template< typename OperatorType >
    explicit gsChebyshevSemiIteration( const OperatorType& mat, T lambda_min, T lambda_max,
                                  const LinOpPtr & precond = LinOpPtr() )
    : Base(mat, precond), m_lambda_min(lambda_min), m_lambda_max(lambda_max) {}

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat        The operator to be solved for, see gsIterativeSolver for details
    /// @param lambda_min The smallest eigenvalue of precond*mat
    /// @param lambda_max The largest eigenvalue of precond*mat
    /// @param precond    The preconditioner, defaulted to the identity
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat, T lambda_min, T lambda_max,
                      const LinOpPtr& precond = LinOpPtr() )
    { return uPtr( new gsChebyshevSemiIteration(mat, lambda_min, lambda_max, precond) ); }
    
    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsChebyshevSemiIteration\n";
        return os;
    }

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;

    T m_lambda_min;
    T m_lambda_max;
    
    T m_c;
    T m_d;
    
    T m_alpha;

    VectorType m_res;
    VectorType m_update;
    VectorType m_tmp;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsChebyshevSemiIteration.hpp)
#endif
