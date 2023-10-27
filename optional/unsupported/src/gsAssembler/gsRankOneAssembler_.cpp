/** @file gsRankOneAssembler_.cpp

    @brief Provides very easy assemblers for mass and stiffness matrix, living
    on the parameter domain. Used in multigrid testing

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <gsAssembler/gsRankOneAssembler.hpp>

namespace gismo
{

    TEMPLATE_INST void assembleRankOneGeneralizedStiffness(const gsGeometry<real_t>& domain, const gsBasis<real_t>& basis, const gsBoundaryConditions<real_t>& bc, real_t alpha, real_t beta, gsSparseMatrix<real_t>& K);

    TEMPLATE_INST void assembleRankOneSimpleBiharmonic(const gsGeometry<real_t>& domain, const gsBasis<real_t>& basis, const gsBoundaryConditions<real_t>& bc, gsSparseMatrix<real_t>& B);

    TEMPLATE_INST void assembleRankOneMass(const gsGeometry<real_t>& domain, const gsBasis<real_t>& basis, const gsBoundaryConditions<real_t>& bc, gsSparseMatrix<real_t>& M);

    TEMPLATE_INST void assembleRankOneMassInverse(const gsGeometry<real_t>& domain, const gsBasis<real_t>& basis, const gsBoundaryConditions<real_t>& bc, gsKroneckerOp<real_t>::Ptr& matrix);

    TEMPLATE_INST void assemble1DMass(const gsGeometry<real_t>& domain, const gsMatrix<real_t>& qp, index_t direction, const gsBasis<real_t>& basis, gsSparseMatrix<real_t>& M);

    TEMPLATE_INST void assemble1DStiffness(const gsGeometry<real_t>& domain, const gsMatrix<real_t>& qp, index_t direction, const gsBasis<real_t>& basis, gsSparseMatrix<real_t>& K);

    TEMPLATE_INST void assemble1D2ndDer(const gsGeometry<real_t>& domain, const gsMatrix<real_t>& qp, index_t direction, const gsBasis<real_t>& basis, gsSparseMatrix<real_t>& B);


}
