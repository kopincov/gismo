/** @file gsRecipeAssemblerBiharmonicSogn.cpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
**/

#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsPde/gsPde.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsRecipeAssembler/gsRecipeAssemblerBiharmonicSogn.h>

namespace gismo {


gsRecipeAssemblerBiharmonicSogn::gsRecipeAssemblerBiharmonicSogn( const gsBiharmonicPde<real_t> &pde)
    : gsRecipeAssembler(pde.domain()), m_pde(pde), m_zeroAverage(true), m_zeroAvgSet(false),m_strategy(dirichlet::elimination)
{
    m_space.resize(1);
    m_space[0]=NULL;
}

void gsRecipeAssemblerBiharmonicSogn::collectEliminatedDofs  ()
{
    const gsBoundaryConditions<real_t>::bcContainer container =
        m_pde.boundaryConditions().allConditions();
    typedef gsBoundaryConditions<real_t>::const_iterator iter;
    for (iter bc=container.begin(); bc!=container.end(); ++bc)
    {
        if (bc->type()==condition_type::dirichlet && m_strategy == dirichlet::elimination)
        {
            std::vector<index_t> patchDofs = m_space[0]->boundaryDofs(bc->ps,0);
            eliminateDofs(patchDofs,0);
        }
    }
}

gsRecipe<real_t>    gsRecipeAssemblerBiharmonicSogn::getPatchRecipe     (index_t ) // patch)
{
    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsLaplaceLaplaceOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getSysWriter());
    result.add(ingr);

    ingr.setOperator(new gsL2TestOp<real_t>(*m_pde.rhs()));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    if (m_zeroAverage)
    {
        static const gsConstantFunction<real_t> my_one_func(1.0,m_domain.dim());

        ingr.setOperator(new gsL2TestOp<real_t>(my_one_func));
        ingr.setTestSpace(0);
        ingr.setUnknownSpace(0);
        ingr.setRule(new gsL2GMapperRhs<MulW>(*m_map,
                                              MulW(SysW(getFreeLimit(),getFreeLimit(),m_sysM,m_rhsMod),0,getFreeLimit()))
                     );
        result.add(ingr);
    }
    return result;
}

gsRecipe<real_t>    gsRecipeAssemblerBiharmonicSogn::getBoundaryRecipe  (patchSide ps)
{
    gsRecipe<real_t>            result;

    const boundary_condition<real_t> *bc1=m_pde.boundaryConditions().getConditionFromSide(ps);
    const boundary_condition<real_t> *bc2=m_pde.boundaryConditionsSecond().getConditionFromSide(ps);

    if (bc1)
    {
        switch(bc1->type())
        {
        case condition_type::dirichlet:
            if (m_strategy==dirichlet::elimination)
            {
                result.add(getDirichletFirstSysIngredient(*bc1));
                result.add(getDirichletFirstRhsIngredient(*bc1));
            }
            else
                GISMO_ERROR("strategy is not implemented");
            break;
        case condition_type::neumann:
        {
            result.add(getNeumannFirstIngredient(*bc1));
            break;
        }
        default:
        {
            GISMO_ERROR("condition type not implemented");
        }
        }
    }
    if (bc2)
    {
        switch(bc2->type())
        {
        case condition_type::dirichlet:
        {
            GISMO_ERROR("condition type not implemented");
            break;
        }
        case condition_type::neumann:
        {
            result.add(getNeumannSecondIngredient(*bc2));
            break;
        }
        default:
            GISMO_ERROR("condition type not implemented");
        }
    }
    return result;
}


gsRecipeIngredient<real_t> gsRecipeAssemblerBiharmonicSogn::getDirichletFirstSysIngredient(const boundary_condition<real_t> & ) // bc)
{

    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryL2ScalarOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getEliSysWriter());

    return ingr;
}

gsRecipeIngredient<real_t> gsRecipeAssemblerBiharmonicSogn::getDirichletFirstRhsIngredient(const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryL2TestOp<real_t>(f));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getEliRhsWriter());

    return ingr;
}

gsRecipeIngredient<real_t>    gsRecipeAssemblerBiharmonicSogn::getNeumannFirstIngredient  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryL2TestVecOp<real_t>(f));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getRhsWriter());

    return ingr;
}

gsRecipeIngredient<real_t>    gsRecipeAssemblerBiharmonicSogn::getNeumannSecondIngredient  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryNormalDerTestVecOp<real_t>(f));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getRhsWriter());

    return ingr;
}

}

