/** @file gsRecipeAssemblerFitting.h

    @brief provides an assembler for fitting to a given function,
    inherited from the gsRecipeAssembler. It assumes a two-dimensional
    domain and will fit the third coordinate to a function.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsRecipeAssembler/gsPhysicalSpace.h>

namespace gismo{

class GISMO_EXPORT gsRecipeAssemblerFitting : public gsRecipeAssembler
{
public:
    const gsFunction<real_t>&        m_f;
    const std::vector<boundary_condition<real_t> > m_emptyVec;

    gsRecipeAssemblerFitting(const gsMultiPatch<real_t>                &geo,
                             const gsFunction<real_t>                  &f,
                             const std::vector<gsBasis<real_t> *>      &bases,
                             const gsWeightMapper<real_t>              &mapper);

    gsRecipeAssemblerFitting(const gsMultiPatch<real_t>                &geo,
                             const gsFunction<real_t>                  &f,
                             std::vector<gsPhysicalSpace* >             space);

    ~gsRecipeAssemblerFitting()
    { freeAll(m_space); }

    gsRecipe<real_t>    getPatchRecipe     (index_t patch);

    gsRecipe<real_t>    getBoundaryRecipe  (patchSide ) // ps)
    { return gsRecipe<real_t>(); }

};
}
