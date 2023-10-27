/** @file gsRecipeAssemblerPoisson.cpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsPde/gsPde.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>
#include <gsCore/gsTransformedFuncSet.h>

namespace gismo {


class gsRecipeAssemblerPoisson2 : public gsRecipeAssemblerPoisson
{
public:
    gsRecipeAssemblerPoisson2(const gsPoissonPde<real_t> &pde)
        : gsRecipeAssemblerPoisson(pde)
    {
    }
gsRecipeAssembler::SpaceList getPatchSpaces (index_t patch )
{
    SpaceList result;
    result.push_back(m_space[0]->getPatchSpace(patch));
    static_cast<gsTransformedFuncSet<real_t> *>(result.back().get())->activeShift+=m_shiftsSource[0];
    result.push_back(gsRecipe<real_t>::spacePtr(new gsRestrictTFS(*m_pde.rhs())));
    result.push_back(gsRecipe<real_t>::spacePtr(new gsConstantFunction<real_t>(1.0,m_domain.dim())));
    return result;
}

gsRecipe<real_t>    getPatchRecipe     (index_t ) // patch)
{
    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsGradGradOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getSysWriter());
    result.add(ingr);

    ingr.setOperator(new gsL2ScalarOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(1);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    if (m_zeroAverage)
    {
        ingr.setOperator(new gsL2ScalarOp<real_t>());
        ingr.setTestSpace(0);
        ingr.setUnknownSpace(2);
        ingr.setRule(new gsL2GMapperRhs<MulW>(*m_map,
                                              MulW(SysW(getFreeLimit(),getFreeLimit(),m_sysM,m_rhsMod),0,getFreeLimit()))
                     );
        result.add(ingr);
    }
    return result;
}
};




gsRecipeAssemblerPoisson::gsRecipeAssemblerPoisson( const gsPoissonPde<real_t> &pde)
    : gsRecipeAssembler(pde.domain()), m_pde(pde), m_zeroAverage(true), m_zeroAvgSet(false),m_strategy(dirichlet::elimination)
{
    m_space.resize(1);
    m_space[0]=NULL;
}

void gsRecipeAssemblerPoisson::collectEliminatedDofs  ()
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

gsRecipe<real_t>    gsRecipeAssemblerPoisson::getPatchRecipe     (index_t ) // patch)
{
    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsGradGradOp<real_t>());
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

gsRecipe<real_t>    gsRecipeAssemblerPoisson::getBoundaryRecipe  (patchSide ps)
{
    const boundary_condition<real_t> *bc=m_pde.boundaryConditions().getConditionFromSide(ps);
    if (!bc)
        return gsRecipe<real_t>();
    switch(bc->type())
    {
    case condition_type::dirichlet:
        if (m_strategy==dirichlet::elimination)
            return getDirichletRecipe(*bc);
        if (m_strategy==dirichlet::nitsche)
            return getNitscheRecipe(*bc);
        else
            GISMO_ERROR("strategy is not implemented");
    case condition_type::neumann:
        return getNeumannRecipe(*bc);
    case condition_type::robin:
        return getRobinRecipe(*bc);
    default:
        GISMO_ERROR("condition type not implemented");
    }
}


gsRecipe<real_t> gsRecipeAssemblerPoisson::getDirichletRecipe(const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryL2ScalarOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getEliSysWriter());
    result.add(ingr);


    ingr.setOperator(new gsBoundaryL2TestOp<real_t>(f));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getEliRhsWriter());
    result.add(ingr);

    return result;
}

gsRecipe<real_t>    gsRecipeAssemblerPoisson::getNitscheRecipe  (const boundary_condition<real_t> &bc)
{
    typedef gsCoeffWriter<SysW> CoefSysW;
    typedef gsCoeffWriter<RhsW> CoefRhsW;
    typedef gsCoeffWriter<MulW> CoefMulW;

    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    const real_t nit = nitche(bc.patch());
    const real_t hpar  = hparam(bc.patch());

    ingr.setOperator(new gsBoundaryL2ScalarOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(new gsL2GMapper<CoefSysW>(*m_map,*m_map,
                                           CoefSysW(SysW(getFreeLimit(),getFreeLimit(),m_sysM,m_rhsMod),nit/hpar)));
    result.add(ingr);

    ingr.setOperator(new gsBoundaryL2TestOp<real_t>(f));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(new gsL2GMapperRhs<CoefRhsW>(*m_map,
                                              CoefRhsW(
                                                  RhsW(getFreeLimit(),m_rhs,gsNullWriter<real_t>::unique_instance),
                                                  nit/hpar)));
    result.add(ingr);

    ingr.setOperator(new gsBoundaryNormalDerTestOp<real_t>(f));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(new gsL2GMapperRhs<CoefRhsW>(*m_map,
                                              CoefRhsW(
                                                  RhsW(getFreeLimit(),m_rhs,gsNullWriter<real_t>::unique_instance),-1)));
    result.add(ingr);

    ingr.setOperator(new gsBoundaryNormalDerValueOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(new gsL2GMapper<CoefMulW>(*m_map,*m_map,
                                           CoefMulW(MulW(SysW(getFreeLimit(),m_sysM,m_rhsMod),0,0),-1)));
    result.add(ingr);

    return result;
}

gsRecipe<real_t>    gsRecipeAssemblerPoisson::getNeumannRecipe  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryL2TestOp<real_t>(f));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    return result;

}

gsRecipe<real_t>    gsRecipeAssemblerPoisson::getRobinRecipe  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t>            result;//=getNeumannRecipe(bc); to avoid repetitions?
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryL2ScalarOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getSysWriter());
    result.add(ingr);


    ingr.setOperator(new gsBoundaryL2TestOp<real_t>(f));
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    return result;
}

real_t gsRecipeAssemblerPoisson::nitche (index_t patch)
{
    gsVector<index_t> deg=getPatchMaxDegree(patch);
    index_t mdeg=deg.array().maxCoeff();
    index_t dim =m_domain.dim();
    return ( (mdeg+dim)* (mdeg+1) * 2.0 );
}

real_t gsRecipeAssemblerPoisson::hparam (index_t patch)
{
    SpaceList  eval = getPatchSpaces(patch);
    real_t   dim  = m_domain.dim();
    real_t   size = eval[0]->size();
    return math::pow( size, -1.0 / dim );
}
}

