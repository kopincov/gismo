/** @file gsParameterDomainAssembler_.cpp

    @brief Provides very easy assemblers for mass and stiffness matrix, living
    on the parameter domain. Used in multigrid testing

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, S. Takacs
*/

#include <gsAssembler/gsParameterDomainAssembler.hpp>

namespace gismo
{

    TEMPLATE_INST void assembleParameterMass(const gsBasis<real_t>& basis1, const gsBasis<real_t>& basis2, gsSparseMatrix<real_t>& M);

    TEMPLATE_INST void assembleParameterStiffness(const gsBasis<real_t>& basis1, const gsBasis<real_t>& basis2, gsSparseMatrix<real_t>& K);

    TEMPLATE_INST void assembleParameter2ndDer1D(const gsBasis<real_t>& basis, gsSparseMatrix<real_t>& B);

    TEMPLATE_INST void handleDirichletConditions(gsSparseMatrix<real_t>& matrix, const gsBoundaryConditions<real_t>& bc, const boxSide& west, const boxSide& east);

    TEMPLATE_INST void assembleGeneralizedParameterStiffnessForTensorProductSpace(const gsBasis<real_t>& basis1, const gsBasis<real_t>& basis2, const gsBoundaryConditions<real_t>& bc, real_t alpha, real_t beta, gsSparseMatrix<real_t>& K);

    TEMPLATE_INST void assembleParameterMassForTensorProductSpace(const gsBasis<real_t>& basis1, const gsBasis<real_t>& basis2, const gsBoundaryConditions<real_t>& bc, gsSparseMatrix<real_t>& K);

    TEMPLATE_INST void assembleParameterMoments(const gsBasis<real_t>& basis, gsFunction<real_t>& f, gsMatrix<real_t>& fh, index_t a, index_t b);

    TEMPLATE_INST void handleDirichletConditionsForMoments(gsMatrix<real_t>& matrix, const gsBoundaryConditions<real_t>& bc, const boxSide& west, const boxSide& east);

    TEMPLATE_INST void assembleParameterMomentsForTensorProduct(const gsBasis<real_t>& basis, const gsBoundaryConditions<real_t>& bc, gsFunction<real_t>& f, gsMatrix<real_t>& fh);

}

