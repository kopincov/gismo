/** @file gsNewtonIterator.h

    @brief Performs Newton iterations to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Goyal, A. Mantzaflaris
*/


#pragma once


namespace gismo
{

/** 
    @brief Performs Newton iterations to solve a nonlinear equation system.
    
    \tparam T coefficient type
    
    \ingroup ThinShell
*/
template <class T> 
class gsNewtonIterator
{
public:

    /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
    gsNewtonIterator(  gsShellAssembler<T> & assembler, 
                       const gsMultiPatch<T> & initialSolution)
    : m_assembler(assembler),
      m_curSolution(initialSolution),
      m_numIterations(0),
      m_maxIterations(100),
      m_tolerance(1e-12),
      m_converged(false)
    { 

    }


public:

    /// Applies Newton method and Performs Newton iterations until convergence or maximum iterations.
    void solve();

    /// Solves linear system obtained using linear elasticity in first step and computes residual
    void firstIteration();

    /// Solves linear system in each iteration based on last solution and computes residual
    void nextIteration();

public:

    /// Returns the latest configuration
    const gsMultiPatch<T> & solution() const {return  m_curSolution;}

    /// Tells if the Newton method converged
    bool converged() const {return m_converged;}

    /// Returns the number of Newton iterations performed
    index_t numIterations() const { return m_numIterations;}

    /// Returns the tolerance value used
    T tolerance() const {return m_tolerance;}

    /// Returns the error after solving the nonlinear system
    T residue()   const {return m_residue;}

    /// Set the maximum number of Newton iterations allowed
    void setMaxIterations(index_t nIter) {m_maxIterations = nIter;}

    /// Set the tolerance for convergence
    void setTolerance(T tol) {m_tolerance = tol;}

protected:

    /// gsShellAssembler object to generate the linear system for each iteration when solving shells
    gsShellAssembler<T> & m_assembler;

    /// The latest/current solution
    gsMultiPatch<T>     m_curSolution;

    /// Solution of the linear system in each iteration
    gsMatrix<T>         m_updateVector;

    /// Linear solver employed
    gsSparseSolver<>::LU  m_solver;
    //gsSparseSolver<>::BiCGSTABDiagonal solver;
    //gsSparseSolver<>::QR  solver;

protected:

    /// Number of Newton iterations performed
    index_t m_numIterations;

    /// Maximum number of Newton iterations allowed
    index_t m_maxIterations;

    /// Tolerance value to decide convergence
    T       m_tolerance;

protected:

    /// Convergence result
    bool m_converged;

    /// Final error
    T m_residue;
};


} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

namespace gismo
{


template <class T> 
void gsNewtonIterator<T>::solve()
{
    firstIteration();

    const T initResidue = m_residue;
    
    // ----- Iterations start -----
    for (m_numIterations = 1; m_numIterations < m_maxIterations; ++m_numIterations)
    {
        nextIteration();
        
        // termination criteria
        if (  m_updateVector.norm()   < m_tolerance || 
             (m_residue / initResidue)< m_tolerance )
        {
            m_converged = true;
            break;
        }
    }
}


template <class T> 
void gsNewtonIterator<T>::firstIteration()
{
    // ----- First iteration -----
    m_converged = false;

    // Construct the linear system
    m_assembler.assemble();

    // Compute the newton update
    m_solver.compute( m_assembler.matrix() );
    m_updateVector = m_solver.solve( m_assembler.rhs() );

    // Update the deformed solution
    m_assembler.updateSolution(m_updateVector, m_curSolution);

    // Compute initial residue
    m_residue = m_assembler.rhs().norm();
}

template <class T> 
void gsNewtonIterator<T>::nextIteration()
{
        // Construct linear system for next iteration
        m_assembler.assemble( m_curSolution );

        // Compute the newton update
        m_solver.compute( m_assembler.matrix() );
        m_updateVector = m_solver.solve( m_assembler.rhs() );
        
        // Update the deformed solution
        m_assembler.updateSolution(m_updateVector, m_curSolution);
        
        // Compute residue
        m_residue = m_assembler.rhs().norm();
        
        gsDebug<<"Iteration: "<< m_numIterations
               <<", residue: "<< m_residue
               <<", update norm: "<<m_updateVector.norm()
               <<"\n";
}



} // namespace gismo

