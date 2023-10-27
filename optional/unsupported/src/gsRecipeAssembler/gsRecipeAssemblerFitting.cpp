/** @file gsRecipeAssemblerFitting.cpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssembler/gsRecipeAssemblerFitting.h>

namespace gismo{

gsRecipeAssemblerFitting::gsRecipeAssemblerFitting(const gsMultiPatch<real_t>                &geo,
                                                   const gsFunction<real_t>                  &f,
                                                   const std::vector<gsBasis<real_t>*>       &bases,
                                                   const gsWeightMapper<real_t>              &mapper)
    : gsRecipeAssembler(geo), m_f(f)
{
    gsPhysicalSpace* physicalSpace = new gsPhysicalSpaceScalar(bases,geo,INVERSE_COMPOSITION,mapper);
    std::vector<gsPhysicalSpace* > phySpace;
    phySpace.push_back(physicalSpace);

    setSpace(phySpace);
}

gsRecipeAssemblerFitting::gsRecipeAssemblerFitting(const gsMultiPatch<real_t>                &geo,
                                                   const gsFunction<real_t>                  &f,
                                                   std::vector<gsPhysicalSpace* >             space)
    : gsRecipeAssembler(geo), m_f(f)
{
    setSpace(space);
}

gsRecipe<real_t>    gsRecipeAssemblerFitting::getPatchRecipe     (index_t ) // patch)
{
    const unsigned limit   = getFreeLimit();

    gsRecipe<real_t>      patch_recipe;
    patch_recipe.resize(2);

    typedef gsMatAndRhsModWriter<gsSparseMatrix<real_t>,gsSparseMatrix<real_t> > Writer;
    typedef gsMatAndRhsModWriter<gsMatrix<real_t>,gsNullWriter<real_t> > WriterRhs;

    patch_recipe[0].setOperator(new gsL2ScalarOp<real_t>());
    patch_recipe[0].setTestSpace(0);
    patch_recipe[0].setUnknownSpace(0);
    patch_recipe[0].setRule(new gsL2GMapper<Writer>(*m_map,*m_map,Writer(limit,limit,m_sysM,m_rhsMod)));

    patch_recipe[1].setOperator(new gsL2TestOp<real_t>(m_f));
    patch_recipe[1].setTestSpace(0);
    patch_recipe[1].setUnknownSpace(0);
    patch_recipe[1].setRule(new gsL2GMapperRhs<WriterRhs>(*m_map,WriterRhs(limit,INT_MAX,m_rhs,gsNullWriter<real_t>::unique_instance)));

    return patch_recipe;
}
}

