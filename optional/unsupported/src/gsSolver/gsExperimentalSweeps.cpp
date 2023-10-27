/** @file gsExperimentalSweeps.cpp

    @brief Collection of some preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

#include <gsSolver/gsExperimentalSweeps.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{
    
void dampedRichardsonSweepBoundary(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau, int numBdNodes)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    gsMatrix<real_t> temp = f - A * x;
    for (int i = 0; i < A.outerSize(); ++i)
    {
        if (i < numBdNodes || i >= A.outerSize() - numBdNodes)
            x(i) += tau * temp(i);
    }
}

// assumes symmetric matrix
void kaczmarzSweepBoundary(const gsSparseMatrix<real_t> & A,
                           gsMatrix<real_t>             & x,
                           const gsMatrix<real_t>       & f,
                           real_t                         /*tau*/,
                           int                            numBdNodes)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    for (int eq = 0; eq < A.outerSize(); ++eq)
    {
        if (!(eq < numBdNodes || eq >= A.outerSize() - numBdNodes))
            continue;

        real_t Asq = 0;
        real_t Aeqx = 0;

        for (gsSparseMatrix<real_t>::InnerIterator it(A,eq); it; ++it)
        {
            Asq += it.value() * it.value();         // compute squared norm of A_eq
            Aeqx += it.value() * x(it.index());     // compute A_eq . x
        }

        const real_t beta = (f(eq) - Aeqx) / Asq;

        for (gsSparseMatrix<real_t>::InnerIterator it(A,eq); it; ++it)
        {
            x(it.index()) += beta * it.value();     // add beta * A_eq
        }
    }
}

void dampedPreRichardsonSweep(const gsSparseMatrix<real_t>& A, const gsSparseMatrix<real_t>& P, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    gsSparseSolver<real_t>::CGDiagonal solver;
    gsMatrix<real_t> corr;
    corr = f - A * x;
    corr = solver.compute( P ).setTolerance(1e-3).solve(corr).eval();
    //corr.array() *= P.diagonal().array();
    x += tau * corr;
}

void dampedPreJacobiSweep(const gsSparseMatrix<real_t>& A, const gsSparseMatrix<real_t>& P, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    //gsMatrix<real_t> temp = f - A * x;
    //temp.array() /= A.diagonal().array();
    //x = tau * temp;

#if 0
    // original version
    gsMatrix<real_t> corr;
    corr = f - A * x;
    corr.array() /= A.diagonal().array();

    corr = gsSparseSolver<>::CGDiagonal( P ).setTolerance(1e-3).solve(corr).eval();
    corr.array() *= P.diagonal().array();
#else
    // symmetricized version
    gsMatrix<real_t> corr;
    corr = f - A * x;

    corr.array() /= A.diagonal().array().sqrt();
    corr.array() *= P.diagonal().array().sqrt();

    corr = gsSparseSolver<real_t>::CGDiagonal( P ).setTolerance(1e-3).solve(corr).eval();

    corr.array() *= P.diagonal().array().sqrt();
    corr.array() /= A.diagonal().array().sqrt();
#endif

    x += tau * corr;
}

void preGaussSeidelSweep(const gsSparseMatrix<real_t> & A,
                         const gsSparseMatrix<real_t> & P,
                         gsMatrix<real_t>             & x,
                         const gsMatrix<real_t>       & f,
                         real_t                         tau,
                         bool                           ) // reverse)
{
    GISMO_ASSERT( A.rows() == x.rows() && x.rows() == f.rows() && A.cols() == A.rows() && x.cols() == f.cols(),
        "Dimensions do not match.");
    
    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    gsMatrix<real_t> corr;
    corr.setZero( x.rows(), 1 );

    //const int start_i = reverse ? A.outerSize() - 1 : 0                 ;
    //const int end_i   = reverse ? 0                 : A.outerSize() - 1 ;
    //const int inc_i   = reverse ? -1                : 1                 ;

    // A is supposed to be symmetric, so it doesn't matter if it's stored in row- or column-major order
    for (int i = 0; i < A.outerSize(); ++i)
    {
        real_t diag = 0.0;
        real_t sum  = 0.0;

        for (gsSparseMatrix<real_t>::InnerIterator it(A,i); it; ++it)
        {
            sum += it.value() * (x( it.index() ) + corr( it.index() ));        // compute A.x
            if (it.index() == i)
                diag = it.value();
        }

        corr(i) = (f(i) - sum) / diag;
    }


    corr = gsSparseSolver<real_t>::CGDiagonal( P )
            .setTolerance(1e-3).solve(corr).eval();
    x.array() += tau * /*P.diagonal().array() * */ corr.array();
}

} // namespace gismo

