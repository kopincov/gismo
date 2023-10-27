/** @file gsExprNewtonIterator.h

    @brief Performs Newton iterations to solve a nonlinear equation system using gsExprAssembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Goyal, A. Mantzaflaris, O. Weeger (Adapted by F.  Beiser)
*/


#pragma once


namespace gismo
{

template <class T, class E1, class E2> 
class gsExprNewtonIterator
{

// Constructor giving access to the gsAssemblerBase object to create a linear system per iteration 
public:

    gsExprNewtonIterator(
                        gsMultiPatch<T> _MultiPatch,
                        gsExprAssembler<T> &_ExprAssembler,
                        const typename gsExprAssembler<T>::geometryMap &_G,
                        const typename gsExprAssembler<T>::space &_H,
                        //typename gsExprAssembler<T>::nonConstVariable _currentSolution, 
                        E1 &_exprJacobian,
                        E2 &_exprResidual,
                        gsMultiPatch<T> &_initialSolution
                        )
    :   m_A(_ExprAssembler)
    ,   m_G( _G )
    ,   m_H( _H )
    ,   m_initSolution( _initialSolution )
    ,   m_curSolution( _initialSolution )
    ,   exprJacobian(_exprJacobian)
    ,   exprResidual(_exprResidual)
    ,   m_numIterations(0)
    ,   m_maxIterations(100)
    ,   m_tolerance(1e-12)
    ,   m_converged(false) 
    ,   m_mp( _MultiPatch )
    { }


// member functions 
public:

    /// \brief Applies Newton method and Performs Newton iterations until convergence or maximum iterations.
    void solve();

    /// \brief Solves linear system in each iteration based on last solution and computes residual
    void firstIteration();

    /// \brief Solves linear system in each iteration based on last solution and computes residual
    void nextIteration();

    /// \brief Returns the latest configuration
    const gsMultiPatch<T> & solution() const { return m_curSolution; }

    /// \brief Tells whether the Newton method converged
    bool converged() const {return m_converged;}

    /// \brief Returns the number of Newton iterations performed
    index_t numIterations() const { return m_numIterations; }

    /// \brief Returns the tolerance value used
    T tolerance() const {return m_tolerance;}

    /// \brief Set the maximum number of Newton iterations allowed
    void setMaxIterations(index_t nIter) {m_maxIterations = nIter;}

    /// \brief Set the tolerance for convergence
    void setTolerance(T tol) {m_tolerance = tol;}

protected:

    void solveLinearProblem();

    T getNorm();

    void updateSolution();


// member variables
public:

    // gsExprAssembler object to generate the linear system for each iteration
    gsExprAssembler<T> & m_A;
    const typename gsExprAssembler<T>::geometryMap & m_G;
    const typename gsExprAssembler<T>::space & m_H;
    gsMultiPatch<T> & m_initSolution;
    gsMultiPatch<T> & m_curSolution;
    E1 & exprJacobian;
    E2 & exprResidual;
    // Solution of the linear system in each iteration
    gsMultiPatch<T>         m_updateMultiPatch;
    // Linear solver employed
    typename gsSparseSolver<T>::LU m_solver;
    // Number of Newton iterations performed
    index_t m_numIterations;
    // Maximum number of Newton iterations allowed
    index_t m_maxIterations;
    // Tolerance value to decide convergence
    T       m_tolerance;
    // Convergence result
    bool m_converged;
    // Norm of the current Newton update vector
	T m_updnorm;

    // tmp
    gsMultiPatch<T> m_mp;
    
};

} // namespace gismo




namespace gismo
{

template <class T, class E1, class E2> 
void gsExprNewtonIterator<T, E1, E2>::solve()
{

    firstIteration();

	const T initUpdate = m_updnorm;

    // ----- Iterations start -----
    for (m_numIterations = 1; m_numIterations < m_maxIterations; ++m_numIterations)
    {
        gsDebug << "Newton iteration " << m_numIterations << ".\n";

        nextIteration();
        
        // termination criteria
        if ( math::abs(m_updnorm / initUpdate) < m_tolerance || 
             math::abs(m_updnorm) < m_tolerance )
        {
            m_converged = true;
            break;
        }
    }

}

template <class T, class E1, class E2> 
void gsExprNewtonIterator<T, E1, E2>::firstIteration()
{
    // Solve the linaer system of the current iteration
    solveLinearProblem();

    // Update the deformed solution
    updateSolution();
    
    // Compute residue
    m_updnorm = getNorm();
    
    gsDebug << "Iteration: "   << m_numIterations << ", "
            << "update norm: " << m_updnorm
            << "\n";

}

template <class T, class E1, class E2> 
void gsExprNewtonIterator<T, E1, E2>::nextIteration()
{
    // Solve the linaer system of the current iteration
    solveLinearProblem();

    // Update the deformed solution
    updateSolution();
    
    // Compute residue
    m_updnorm = getNorm();
    
    gsDebug << "Iteration: "   << m_numIterations << ", "
            << "update norm: " << m_updnorm
            << "\n";

}

template <class T, class E1, class E2> 
void gsExprNewtonIterator<T, E1, E2>::solveLinearProblem()
{

    m_initSolution = m_curSolution; 

    ////////////////////////////////////////
    // Assembling
    ////////////////////////////////////////

    // Initialize system
    m_A.initSystem();
    // Assemble stiffness matrix
    m_A.assemble( exprJacobian, exprResidual );

    ////////////////////////////////////////
    // Solving
    ////////////////////////////////////////

    // Initialize solver
    gsSparseSolver<>::CGDiagonal solverCGD;
    solverCGD.compute( m_A.matrix() );

    // Solve for right hand side 
    gsMatrix<T> updateVector = solverCGD.solve( m_A.rhs() );
    
    typename gsExprAssembler<T>::solution updateSolution = m_A.getSolution( m_H, updateVector );
    updateSolution.extract( m_updateMultiPatch );

    gsWriteParaview( m_updateMultiPatch, "updateMultiPatch" );

}

template <class T, class E1, class E2>
void gsExprNewtonIterator<T, E1, E2>::updateSolution()
{

    for(size_t i = 0; i < m_initSolution.nPatches(); ++i)
    {
        gsMatrix<T> & currentCoefs = m_curSolution.patch(i).coefs();
        gsMatrix<T> & updateCoefs = m_updateMultiPatch.patch(i).coefs();

        for(index_t j = 0; j < currentCoefs.size(); ++j )
        {
            currentCoefs(j) += updateCoefs(j);  
        }
    }

}

template <class T, class E1, class E2> 
T gsExprNewtonIterator<T, E1, E2>::getNorm()
{

    gsExprEvaluator<T> ExprEvaluator( m_A );
    typename gsExprAssembler<T>::variable updateVariable = m_A.getCoeff( m_updateMultiPatch );
    
    return math::sqrt( ExprEvaluator.integral( updateVariable.sqNorm() * meas(m_G) ) );

}

} // namespace gismo

