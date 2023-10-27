/** @file gsRecipeAssemblerLinElast.cpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <iterator>

#include <gsPde/gsLinearElasticityPde.h>
#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsRecipeAssembler/gsRecipeAssemblerLinElast.h>

// In most parts copied from gsRecipeAssemblerStokes.
namespace gismo {

// copied from gsRecipeAssemblerStokes, not checked yet.
void gsRecipeAssemblerLinElast::collectEliminatedDofs  ()
{
    const gsBoundaryConditions<real_t>::bcContainer container =
        m_pde.boundaryConditions().allConditions();
    typedef gsBoundaryConditions<real_t>::const_iterator iter;
    for (iter bc=container.begin(); bc!=container.end(); ++bc)
    {
        if (bc->type()==condition_type::dirichlet )
        {
            std::vector<index_t> patchDofs = m_space[0]->boundaryDofs(bc->ps,0);
            eliminateDofs(patchDofs,0);
        }
    }
}

gsRecipe<real_t>    gsRecipeAssemblerLinElast::getPatchRecipe     (index_t ) // patch)
{
    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    // new PDEOperator used here:
    ingr.setOperator(new gsLinElastOp<real_t>(m_YoungsModulus,m_PoissonsRatio));    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getSysWriter());
    result.add(ingr);

    // copied from gsRecipeAssemblerStokes
    ingr.setOperator(new gsL2TestVecOp<real_t>(*m_pde.rhs()));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getRhsWriter());
    result.add(ingr);


    return result;
}

// copied from gsRecipeAssemblerStokes, not checked yet.
gsRecipe<real_t>    gsRecipeAssemblerLinElast::getBoundaryRecipe  (patchSide ps)
{
    const boundary_condition<real_t> *bc = m_pde.boundaryConditions().getConditionFromSide(ps);
    if (!bc)
        return gsRecipe<real_t>();

    switch(bc->type())
    {
    case condition_type::dirichlet:
        return getDirichletRecipe(*bc);
    case condition_type::neumann:
        return getNeumannRecipe(*bc);
    default:
        GISMO_ERROR("condition type not implemented");
    }
}

// copied from gsRecipeAssemblerStokes, not checked yet.
gsRecipe<real_t> gsRecipeAssemblerLinElast::getDirichletRecipe  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t> result;
    gsRecipeIngredient<real_t> ingr;

    {
        ingr.setOperator    ( new gsBoundaryL2ScalarOp<real_t>());
        ingr.setTestSpace   ( 0);
        ingr.setUnknownSpace( 0);
        ingr.setRule        ( getEliSysWriter() );
    }
    result.add(ingr);

    if( !bc.isHomogeneous() )
    {
        ingr.setOperator    (new gsBoundaryL2TestVecOp<real_t>(f) );
        ingr.setTestSpace   ( 0);
        ingr.setUnknownSpace( 0);
        ingr.setRule        ( getEliRhsWriter());

        result.add(ingr);
    }

    return result;
}

// copied from gsRecipeAssemblerStokes, not checked yet.
gsRecipe<real_t>    gsRecipeAssemblerLinElast::getNeumannRecipe  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryL2TestVecOp<real_t>(f));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    return result;

}

// overrides function from gsRecipeAssembler for over-integration
//gsIntegrationRule gsRecipeAssemblerLinElast::getPatchIntegration     (index_t patch )
//{
//    gsIntegrationRule result;
//    result.rule=gsGaussRule<real_t>(getPatchMaxDegree(patch).array()+2);
//    result.subdomains=m_space[0]->patchPartition(patch);
//    return result;
//}


}
