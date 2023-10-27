/** @file gsOptElectricMotor.h

    @brief Provides declaration of an optimization problem.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s):  R. Schneckenleitner, A. Mantzaflaris
*/

# include <gsCore/gsLinearAlgebra.h>
# include <gsIETI/gsIETIUtils.h>
# include <gsIETI/gsIETIAssembler.h>
# include <gsIETI/gsIETISolver.h>
# include <gsIETI/gsIETIScaledDirichlet.h>
# include <gsAssembler/gsMagnetostaticAssembler.h>
# include <gsAssembler/gsMagnetostaticAdjointAssembler.h>
# include <iostream>
# include <gsAssembler/gsMagnetostaticShapeDerivAssembler_decoupled.h>

#include <gsIpopt/gsOptProblem.h>

#pragma once

namespace gismo
{

/**
   \brief Class defining an optimization problem
*/

template <typename T>
class gsOptElectricMotor : public gsOptProblem<T>
{

    //friend class gsIpOptTNLP<T>;

public:

    gsOptElectricMotor() {};
    //*****
    /** constructor */
    gsOptElectricMotor(gsMultiPatch<T> domain, gsMultiBasis<T> mbasis, gsPiecewiseFunction<T> magRel, gsPiecewiseFunction<T> magnetization, gsFunctionExpr<T> reference, gsPiecewiseFunction<T> penalization,
                       gsOptionList opt, gsBoundaryConditions<> bcInfo, std::vector<size_t> desDomain, std::vector<size_t> desNeighbours, unsigned maxIter = 1);
    // stored as members

public:

    virtual T evalObj( const gsAsConstVector<T> & u ) const;

    virtual void gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const;
    
    virtual void evalCon_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const;

    virtual void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const;

public:
    
    void solve ();
    
    const gsMultiPatch<T> & currentPatches() const
    {
        gsAsConstVector<T> curDes (m_curDesign.data(), m_numDesignVars);
        updateTempPatches(curDes);
        return m_tmpPatches;
    }

    /// Method to set the tolerance from outside, 1.0e-8 by default
    void setTolerance(const T eps);

protected:

    typedef gsExprAssembler<real_t>::space space;
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::variable variable;
    typedef gsExprAssembler<real_t>::solution solution;

protected: // Inherited members from gsOptProblem

    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;

    using gsOptProblem<T>::numIterations;
    using gsOptProblem<T>::finalObjective;

protected: // Own members

    /// domain (object of interest )
    mutable gsMultiPatch<T> m_tmpPatches;

    /// Basis of the (perturbed) domain
    gsMultiBasis<> m_basis;

    /// Function which represents the magnetic reluctivity
    gsPiecewiseFunction<T> m_reluctivity;

    /// The reference function which belongs to the ideal property of the motor
    gsFunctionExpr<T> m_reference;

    /// Function which represents the rhs
    gsPiecewiseFunction<T> m_rhs;

    /// Penalization function for the calculation of the shape gradient
    gsPiecewiseFunction<T> m_penalization;

    /// The options for the Assembler
    gsOptionList m_opt;

    /// The boundary conditions for solving the PDEs
    gsBoundaryConditions<> m_bc;

    /// Vector which specifies the design patches
    std::vector<size_t> m_desDomain;

    /// Vector which specifies the neighbours of the design patches
    std::vector<size_t> m_desNeighbours;

    /// Collocation points for Jacobian determinant
    std::vector<gsMatrix<T> > m_colPoints;

    /// Basis for the jacobian determinant
    gsMultiBasis<T> m_jacBasis;

    /// Keeps the indices of the corner coefficients
    std::vector<gsVector<unsigned> > m_corner;

    /// Vector which contains dof mappers for the design areas
    gsDofMapper m_mapper;

    /// Helper for the sparse implementation
    std::vector<std::vector<std::pair<index_t, index_t> > > m_preImages;

    /// Maximum number of iterations
    unsigned m_maxIterations;

    /// Tolerance for the stopping criterion, default value is 1e-4
    T m_tolerance;

    /// Sparse factorization of the collocation matrix
    //std::vector<typename gsSparseSolver<T>::LU > m_jacSolver;
    mutable typename gsSparseSolver<T>::LU  m_jacSolver;

    /// Matrix which stores the global indices of the design variables
    std::vector<index_t> m_patchOrientation;

private:

    /**@name Methods to block default compiler methods.
     * The compiler automatically generates the following three methods.
     *  Since the default compiler implementation is generally not what
     *  you want (for all but the most simple classes), we usually 
     *  put the declarations of these methods in the private section
     *  and never implement them. This prevents the compiler from
     *  implementing an incorrect "default" behavior without us
     *  knowing. (See e.g. Scott Meyers book, "Effective C++")
     *  
     */
    //@{
    gsOptElectricMotor(const gsOptElectricMotor & );
    gsOptElectricMotor& operator=(const gsOptElectricMotor & );
    //@}

private:
    /// Method to compute the solution u of the state equation
    gsMultiPatch<T> solveStateForMultipatch(gsMagnetostaticPde<T> magPDE, gsOptionList opt, gsPiecewiseFunction<T> f) const;

    gsMatrix<T> solveState(gsMagnetostaticAssembler <T> & magAss) const;

    /// Method to compute the solution p of the state equation
    gsMultiPatch<T> solveAdjoint(gsMagnetostaticAdjointPde<T> magPDE, gsMultiPatch<T> mp, gsOptionList opt) const;

    /// Method to compute the gradient of the shape derivative
    gsMultiPatch<T> calcShapeDerivative(gsMagnetostaticShapeDerivPde<T> shapeD, gsOptionList opt, gsBoundaryConditions<> bcInfo) const;

    /// Method to create a multipatch from a vector
    // Probably not needed
    //void vectorToMultipatch(const gsAsConstVector<T> & u, gsMultiPatch<T> & mp) const;

    /// Method to create a gsAsVector object from a multipatch
    gsMatrix<T> multipatchToVector(const gsMultiPatch<T> & mp) const;

    /// Method to update the domain *****************
    void updateTempPatches(const gsAsConstVector<T> & u) const;

    /// Method to update the inner control points
    void patchFromBoundary(const gsAsConstVector<T> & u) const;

    //updateCurrentDesign

    /// Method to get the coefficients of the jacobian determinant of the design domain
    void jacobianDetCoefs(const gsAsConstVector<T> & u, gsAsVector<T> & result) const;

    /// Compute coefficients of the jacobian in the initial state
    bool computeInitialjacobianDetCoefs(const size_t designPatch);

    /// Compute the non-zero entries in the jacobian
    void computeJacStructure();


};


} // end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsOptElectricMotor.hpp)
#endif

