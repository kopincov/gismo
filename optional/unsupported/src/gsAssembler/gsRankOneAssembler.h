/** @file gsRankOneAssembler.h

    @brief Provides very easy assemblers for mass and stiffness matrix, living
    on the parameter domain. Used in multigrid testing

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#pragma once

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsSolver/gsKronecker.h>

namespace gismo
{

/// @brief RankOneAssembler for the stiffness and/or stiffness matrix.
///
/// Assembles the generalized stiffness matrix alpha*stiffness + beta*mass assuming basis to be a tensor-product basis (checked only at run-time)
/// and eliminates Dirichlet dofs if needed
///
/// The computation uses the tensor-product structure and is (for d>1) very fast therefore
///
/// \param[in]  domain   the domain (multipatch)
/// \param[in]  basis    the basis
/// \param[in]  bc       the boundary conditions object
/// \param[in]  alpha    the scaling factor for the stiffness matrix
/// \param[in]  beta     the scaling factor for the mass matrix
/// \param[out] matrix   the assembled matrix
///
/// \ingroup Solver
template <typename T>
void assembleRankOneGeneralizedStiffness(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, T alpha, T beta, gsSparseMatrix<T>& matrix);

/// @brief RankOneAssembler for the stiffness and/or mass matrix.
///
/// Assembles the stiffnes matrix assuming basis to be a tensor-product basis (checked only at run-time)
/// and eliminates Dirichlet dofs if needed
///
/// The computation uses the tensor-product structure and is (for d>1) very fast therefore
///
/// \param[in]  domain   the domain (multipatch)
/// \param[in]  basis    the basis
/// \param[in]  bc       the boundary conditions object
/// \param[out] matrix   the assembled matrix
///
/// \ingroup Solver
template <typename T>
void assembleRankOneStiffness(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& matrix) {
    assembleRankOneGeneralizedStiffness(domain, basis, bc, (T)1, (T)0, matrix);
}

/// @brief RankOneAssembler for the simplified biharmonic.
///
/// Assembles the simplified biharmonic, i.e., not the cross derivatives.
/// Assuming basis to be a tensor-product basis (checked only at run-time)
/// and eliminates Dirichlet (of first kind) dofs if needed
///
/// The computation uses the tensor-product structure and is (for d>1) very fast therefore
///
/// \param[in]  domain   the domain (multipatch)
/// \param[in]  basis    the basis
/// \param[in]  bc       the boundary conditions object
/// \param[in]  alpha    the scaling factor for the stiffness matrix
/// \param[in]  beta     the scaling factor for the mass matrix
/// \param[out] matrix   the assembled matrix
///
/// \ingroup Solver
template <typename T>
void assembleRankOneSimpleBiharmonic(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& B);

/// @brief RankOneAssembler for the mass matrix.
///
/// Assembles the mass matrix assuming basis to be a tensor-product basis (checked only at run-time)
/// and eliminates Dirichlet dofs if needed
///
/// The computation uses the tensor-product structure and is (for d>1) very fast therefore
///
/// \param[in]  domain   the domain (multipatch)
/// \param[in]  basis    the basis
/// \param[in]  bc       the boundary conditions object
/// \param[out] matrix   the assembled matrix
///
/// \ingroup Solver
template <typename T>
void assembleRankOneMass(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& matrix);

/// @brief RankOneAssembler for the inverse of the mass matrix.
///
/// Assembles the mass matrix assuming basis to be a tensor-product basis (checked only at run-time)
/// and eliminates Dirichlet dofs if needed
///
/// The computation uses the tensor-product structure and is (for d>1) very fast therefore
///
/// \param[in]  domain   the domain (multipatch)
/// \param[in]  basis    the basis
/// \param[in]  bc       the boundary conditions object
/// \param[out] matrix   the assembled matrix
///
/// \ingroup Solver
template <typename T>
void assembleRankOneMassInverse(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, typename gsKroneckerOp<T>::Ptr& matrix);

template <typename T>
void assemble1DMass(const gsGeometry<T>& domain, const gsMatrix<T>& qp, index_t direction, const gsBasis<T>& basis, gsSparseMatrix<T>& M);

template <typename T>
void assemble1DStiffness(const gsGeometry<T>& domain, const gsMatrix<T>& qp, index_t direction, const gsBasis<T>& basis, gsSparseMatrix<T>& K);

template <typename T>
void assemble1D2ndDer(const gsGeometry<T>& domain, const gsMatrix<T>& qp, index_t direction, const gsBasis<T>& basis, gsSparseMatrix<T>& B);

}
