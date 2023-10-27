/** @file gsRecipeAssemblerKKTPervertedPoisson.cpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, J. Sogn
**/

#include <iterator>

#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsRecipeAssembler/gsRecipeAssemblerKKTPervertedPossion.h>


namespace gismo {


gsRecipeAssemblerKKTPervertedPoisson::gsRecipeAssemblerKKTPervertedPoisson(const gsPoissonPde<real_t> &pde, const gsBoundaryConditions<real_t> &desired, real_t alpha)
    : gsRecipeAssembler(pde.domain()), m_pde(pde), m_desired(desired), m_alpha(alpha)
{
    m_space.resize(1);
    m_space[0]=NULL;
}

void gsRecipeAssemblerKKTPervertedPoisson::collectEliminatedDofs  ()
{
    const gsBoundaryConditions<real_t>::bcContainer container =
        m_pde.boundaryConditions().allConditions();
    typedef gsBoundaryConditions<real_t>::const_iterator iter;
    for (iter bc=container.begin(); bc!=container.end(); ++bc)
    {
        if (bc->type()==condition_type::dirichlet )
        {
            std::vector<index_t> patchDofs = m_space[stateSpace]->boundaryDofs(bc->ps,0);
            eliminateDofs(patchDofs,stateSpace);
        }
    }
}

gsRecipe<real_t>    gsRecipeAssemblerKKTPervertedPoisson::getPatchRecipe     (index_t ) // patch)
{
    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;
    ingr.setOperator(new gsL2ScalarOp<real_t>());
    ingr.setTestSpace(controlSpace);
    ingr.setUnknownSpace(controlSpace);
    ingr.setRule(getSysWriter(m_alpha));
    result.add(ingr);

    ingr.setOperator(new gsL2ScalarOp<real_t>());
    ingr.setTestSpace(controlSpace);
    ingr.setUnknownSpace(multiplierSpace);
    ingr.setRule(getLagrangeMultiplierWriter());
    result.add(ingr);

    //Operatior not tested!
    ingr.setOperator(new gsLaplaceOp<real_t>());
    ingr.setTestSpace(multiplierSpace);
    ingr.setUnknownSpace(stateSpace);
    //ingr.setRule(getSysWriter(-1));
    ingr.setRule(getLagrangeMultiplierWriter(-1.0));
    result.add(ingr);

    return result;
}

gsRecipe<real_t>    gsRecipeAssemblerKKTPervertedPoisson::getBoundaryRecipe  (patchSide ps)
{
    gsRecipe<real_t> result;
    gsRecipeIngredient<real_t> ingr;

    //The surface integral
    {
    //The rhs
    const boundary_condition<real_t> *bc_d = m_desired.getConditionFromSide(ps);
    ingr.setOperator(new gsBoundaryNormalDerTestNormalOp<real_t>(*bc_d->function().get()) );
    ingr.setTestSpace(stateSpace);
    ingr.setUnknownSpace(stateSpace);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    ingr.setOperator(new gsBoundaryNormalDerNormalDerOp<real_t>());
    ingr.setTestSpace(stateSpace);
    ingr.setUnknownSpace(stateSpace);
    ingr.setRule(getSysWriter());
    result.add(ingr);
    }

    //The boundary condition
    const boundary_condition<real_t> *bc = m_pde.boundaryConditions().getConditionFromSide(ps);
    if (!bc)
        return result;

    switch(bc->type())
    {
    case condition_type::dirichlet:
    {
        gsFunction<real_t> &f=*bc->function().get();

        ingr.setOperator    ( new gsBoundaryL2ScalarOp<real_t>());
        ingr.setTestSpace   ( stateSpace);
        ingr.setUnknownSpace( stateSpace);
        ingr.setRule        ( getEliSysWriter() );
        result.add(ingr);

        if( !bc->isHomogeneous() )
        {
            ingr.setOperator    (new gsBoundaryL2TestVecOp<real_t>(f) );
            ingr.setTestSpace   ( stateSpace);
            ingr.setUnknownSpace( stateSpace);
            ingr.setRule        ( getEliRhsWriter());

            result.add(ingr);
        }
        return result;
    }
    case condition_type::neumann:
        GISMO_ERROR("condition type not implemented");
    default:
        GISMO_ERROR("condition type not implemented");
    }
}


gsIntegrationRule gsRecipeAssemblerKKTPervertedPoisson::getPatchIntegration     (index_t patch )
{
    gsIntegrationRule result;
    result.rule=gsGaussRule<real_t>(getPatchMaxDegree(patch).array()+2);
    result.subdomains=m_space[0]->patchPartition(patch);
    return result;
}


}
