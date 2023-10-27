/** @file gsParameterDomainAssembler.h

    @brief Provides very easy assemblers for mass and stiffness matrix, living
    on the parameter domain. Used in multigrid testing

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, S. Takacs
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsSolver/gsKronecker.h>
#include <gsSolver/gsPatchPreconditionersCreator.h>

namespace gismo
{

// Wrappers for stuff that is now available in stable
    
template <typename T>
GISMO_DEPRECATED void assembleParameterMass(const gsBasis<T>& basis, gsSparseMatrix<T>& M)
{
    M = gsPatchPreconditionersCreator<T>::massMatrix(basis);
}

template <typename T>
GISMO_DEPRECATED void assembleParameterStiffness(const gsBasis<T>& basis, gsSparseMatrix<T>& K)
{
    K = gsPatchPreconditionersCreator<T>::stiffnessMatrix(basis);
}

template <typename T>
GISMO_DEPRECATED void assembleGeneralizedParameterStiffnessForTensorProductSpace(const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, T alpha, T beta, gsSparseMatrix<T>& matrix)
{
    // Assembles alpha*stiffness + beta*mas
    matrix = alpha * gsPatchPreconditionersCreator<T>::stiffnessMatrix(basis, bc, gsAssembler<T>::defaultOptions(), beta/alpha);
}

template <typename T>
GISMO_DEPRECATED void assembleParameterStiffnessForTensorProductSpace(const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& matrix) {
    matrix = gsPatchPreconditionersCreator<T>::stiffnessMatrix(basis, bc);
}

template <typename T>
GISMO_DEPRECATED void assembleParameterMassForTensorProductSpace(const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& matrix)
{
    matrix = gsPatchPreconditionersCreator<T>::massMatrix(basis, bc); 
}

template <typename T>
GISMO_DEPRECATED void assembleParameterMassInverseForTensorProductSpace(const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, typename gsKroneckerOp<T>::Ptr& matrix)
{
    
    matrix = memory::make_shared( dynamic_cast<gsKroneckerOp<T>*>( gsPatchPreconditionersCreator<T>::massMatrixInvOp(basis, bc).release() ) );
}

// Other stuff that is not available in stable

/// @brief ParameterDomainAssembler for the mass matrix.
///
/// Assembles the mass matrix on the parametric domain (no geometry transformation)
///
/// \param[in]  basis1   the basis for the ansatz functions
/// \param[in]  basis2   the basis for the test functions
/// \param[out] M        the assembled matrix
///
/// \ingroup Solver 
template <typename T>
void assembleParameterMass(const gsBasis<T>& basis1, const gsBasis<T>& basis2, gsSparseMatrix<T>& M);

/// @brief ParameterDomainAssembler for the stiffness matrix.
///
/// Assembles the stiffness matrix on the parametric domain (no geometry transformation)
///
/// \param[in]  basis1   the basis for the ansatz functions
/// \param[in]  basis2   the basis for the test functions
/// \param[out] K        the assembled matrix
///
/// \ingroup Solver 
template <typename T>
void assembleParameterStiffness(const gsBasis<T>& basis1, const gsBasis<T>& basis2, gsSparseMatrix<T>& K);

/// @brief ParameterDomainAssembler for the 1 dimensional 2nd derivativ matrix.
///
/// Assembles the 2nd derivativ matrix on the parametric domain for 1D (no geometry transformation)
///
/// \param[in]  basis    the basis
/// \param[out] B        the assembled matrix
///
/// \ingroup Solver
template <typename T>
void assembleParameter2ndDer1D(const gsBasis<T>& basis, gsSparseMatrix<T>& B);

/// @brief ParameterDomainAssembler module for eliminating the Dirichlet conditions for a 1D matrix.
///
/// NB! Only work for Dirichlet BC of first kind (Standard Poisson).
/// \param[in,out]  matrix    the matrix
/// \param[in]      bc        the boundary conditions object
/// \param[in]      west      specify which boxSide corresponds to the first dof
/// \param[in]      east      specify which boxSide corresponds to the last dof
///
/// \ingroup Solver 
template <typename T>
void handleDirichletConditions(gsSparseMatrix<T>& matrix, const gsBoundaryConditions<T>& bc, const boxSide& west, const boxSide& east);

/// @brief ParameterDomainAssembler for the stiffness and/or stiffness matrix.
///
/// Assembles the generalized stiffness matrix alpha*stiffness + beta*mass assuming basis to be a tensor-product basis (checked only at run-time)
/// and eliminates Dirichlet dofs if needed
///
/// The computation uses the tenspr-product structure and is (for d>1) very fast therefore
///
/// \param[in]  basis1   the basis for the ansatz functions
/// \param[in]  basis2   the basis for the test functions
/// \param[in]  bc       the boundary conditions object
/// \param[in]  alpha    the scaling factor for the stiffness matrix
/// \param[in]  beta     the scaling factor for the mass matrix
/// \param[out] matrix   the assembled matrix
///
/// \ingroup Solver 
template <typename T>
void assembleGeneralizedParameterStiffnessForTensorProductSpace(const gsBasis<T>& basis1, const gsBasis<T>& basis2, const gsBoundaryConditions<T>& bc, T alpha, T beta, gsSparseMatrix<T>& matrix);


/// @brief ParameterDomainAssembler for the stiffness matrix.
///
/// Assembles the stiffnes matrix assuming basis to be a tensor-product basis (checked only at run-time)
/// and eliminates Dirichlet dofs if needed
///
/// The computation uses the tenspr-product structure and is (for d>1) very fast therefore
///
/// \param[in]  basis1   the basis for the ansatz functions
/// \param[in]  basis2   the basis for the test functions
/// \param[in]  bc       the boundary conditions object
/// \param[out] matrix   the assembled matrix
///
/// \ingroup Solver 
template <typename T>
void assembleParameterStiffnessForTensorProductSpace(const gsBasis<T>& basis1, const gsBasis<T>& basis2, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& matrix) {
    assembleGeneralizedParameterStiffnessForTensorProductSpace(basis1, basis2, bc, (T)1, (T)0, matrix);
}

/// @brief ParameterDomainAssembler for the mass matrix.
///
/// Assembles the mass matrix assuming basis to be a tensor-product basis (checked only at run-time)
/// and eliminates Dirichlet dofs if needed
///
/// The computation uses the tenspr-product structure and is (for d>1) very fast therefore
///
/// \param[in]  basis1   the basis for the ansatz functions
/// \param[in]  basis2   the basis for the test functions
/// \param[in]  bc       the boundary conditions object
/// \param[out] matrix   the assembled matrix
///
/// \ingroup Solver 
template <typename T>
void assembleParameterMassForTensorProductSpace(const gsBasis<T>& basis1, const gsBasis<T>& basis2, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& matrix);

/// @brief ParameterDomainAssembler for the moments.
///
/// Assembles the moments vector on the parametric domain (no geometry transformation)
///
/// \param[in]  basis    the basis
/// \param[in]  f        the function representing the moments
/// \param[out] fh       the assembled moments vector
/// \param[in]  a
/// \param[in]  b
///
/// \ingroup Solver 
template <typename T>
void assembleParameterMoments(const gsBasis<T>& basis, gsFunction<T>& f, gsMatrix<T>& fh, index_t a=0, index_t b=0);

/// @brief ParameterDomainAssembler module for eliminating the Dirichlet conditions for a 1D moments vector.
///
/// \param[in,out]  fh        the moments vector
/// \param[in]      bc        the boundary conditions object
/// \param[in]      west      specify which boxSide corresponds to the first dof
/// \param[in]      east      specify which boxSide corresponds to the last dof
///
/// \ingroup Solver 
template <typename T>
void handleDirichletConditionsForMoments(gsMatrix<T>& fh, const gsBoundaryConditions<T>& bc, const boxSide& west, const boxSide& east);

/// @brief ParameterDomainAssembler for the moments.
///
/// Assembles the moments assuming basis to be a tensor-product basis (checked only at run-time) and eliminates
/// Dirichlet dofs if needed (only homogenous Dirichlet dofs!)
///
/// @note This uses that function f has a tensor-product structure.
///
/// \param[in]  basis    the basis
/// \param[in]  bc       the boundary conditions object
/// \param[in]  f        the moments function
/// \param[out] fh       the assembled moments vector
///
/// \ingroup Solver 
template <typename T>
void assembleParameterMomentsForTensorProduct( const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsFunction<T>& f, gsMatrix<T>& fh);

}
