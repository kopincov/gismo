
/** @file gsRecipeAssemblerStokes.cpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <iterator>

#include <gsPde/gsStokesPde.h>
#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsRecipeAssembler/gsRecipeAssemblerStokes.h>

namespace gismo {


template <typename T>
T* Cloner(const T* a)
{
    return (T*)a->clone().release();
}

template <typename T>
T* Constifier(T* a)
{
    return a; //const_cast<const T*>(a);
}

GISMO_EXPORT void DegreeElevate( gsBasis<real_t>* a) {a->degreeElevate();}



std::vector<gsPhysicalSpace*> constructTHSpaces (
        const std::vector<gsBasis<real_t> *>       &bases,
        const gsMultiPatch<real_t>                 &geo,
        //    const gsMapFactory                         &factory,
        std::vector<std::vector<gsBasis<>*> >      *outBasis
        )
{
    const index_t dim=bases[0]->domainDim();
    std::vector<gsPhysicalSpace*>        result(2);

    std::vector<gsBasis<real_t>*>  temp(bases.size()), patchBasesP, patchBasisV;
    {
        cloneAll(bases.begin(),bases.end(),temp.begin());
        // build pressure space
        gsMultiBasis<real_t>                 multipatch (temp,geo);
        gsMapFactoryMultiBasis               mapFactory (multipatch);
        std::vector<gsBasis<real_t>*>        argVec;

        std::transform(bases.begin(),bases.end(),std::back_inserter(patchBasesP), Cloner<gsBasis<real_t> >);
        std::transform(patchBasesP.begin(),patchBasesP.end(),std::back_inserter(argVec), Constifier<gsBasis<real_t> >);
        result[pressure]=new gsPhysicalSpaceScalar(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
    }

    {   // build velocity space
        std::transform(bases.begin(),bases.end(),std::back_inserter(patchBasisV), Cloner<gsBasis<real_t> >);
        std::for_each(patchBasisV.begin(),patchBasisV.end(), DegreeElevate);

        temp.resize(patchBasisV.size());
        cloneAll(patchBasisV.begin(),patchBasisV.end(),temp.begin());
        gsMultiBasis<real_t>                 multipatch (temp,geo);
        gsMapFactoryMultiBasis               mapFactory (multipatch);
        std::vector<gsBasis<real_t>*>        argVec;

        std::transform(patchBasisV.begin(),patchBasisV.end(),std::back_inserter(argVec), Constifier<gsBasis<real_t> >);
        gsPhysicalSpaceScalar vel(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
        std::vector<gsPhysicalSpaceScalar*> tmp(dim,&vel);
        result[velocity]=(new gsPhysicalSpaceVector(tmp));
    }

    if(outBasis)
    {
        *outBasis= std::vector<std::vector<gsBasis<>*> >(1,patchBasisV);
        outBasis->push_back(patchBasesP);
    }
    return result;
}

gsRecipeAssemblerStokes::gsRecipeAssemblerStokes( const gsStokesPde<real_t> &pde)
    : gsRecipeAssembler(pde.domain()), m_pde(pde), m_zeroAverage(false)
{
    m_space.resize(1);
    m_space[0]=NULL;
    if (m_pde.boundaryConditions().dirichletSides().size()==static_cast<size_t>(m_domain.nBoundary()))
        m_zeroAverage=true;
}

void gsRecipeAssemblerStokes::collectEliminatedDofs  ()
{
    const gsBoundaryConditions<real_t>::bcContainer container =
        m_pde.boundaryConditions().allConditions();
    typedef gsBoundaryConditions<real_t>::const_iterator iter;
    for (iter bc=container.begin(); bc!=container.end(); ++bc)
    {
        if (bc->type()==condition_type::dirichlet )
        {
            std::vector<index_t> patchDofs = m_space[velocity]->boundaryDofs(bc->ps,0);
            eliminateDofs(patchDofs,0);
        }
    }
}

gsRecipe<real_t>    gsRecipeAssemblerStokes::getPatchRecipe     (index_t ) // patch)
{
    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsGradGradOp<real_t>());
    ingr.setTestSpace(velocity);
    ingr.setUnknownSpace(velocity);
    ingr.setRule(getSysWriter());
    result.add(ingr);

    ingr.setOperator(new gsDivergenceOp<real_t>());
    ingr.setTestSpace( pressure);
    ingr.setUnknownSpace( velocity);
    ingr.setRule(getLagrangeMultiplierWriter());
    result.add(ingr);


    ingr.setOperator(new gsL2TestVecOp<real_t>(*m_pde.rhs()));
    ingr.setTestSpace(velocity);
    ingr.setUnknownSpace(velocity);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    if (m_zeroAverage)
    {
        static const gsConstantFunction<real_t> my_one_func(1.0,m_domain.dim());

        ingr.setOperator(new gsL2TestOp<real_t>(my_one_func));
        ingr.setTestSpace(pressure);
        ingr.setUnknownSpace(pressure);
        ingr.setRule(new gsL2GMapperRhs<MulW>(*m_map,
                                              MulW(SysW(getSysSize(),getSysSize(),m_sysM,m_rhsMod),0,getFreeLimit()))
                     );
        result.add(ingr);
    }
    return result;
}

gsRecipe<real_t>    gsRecipeAssemblerStokes::getBoundaryRecipe  (patchSide ps)
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

gsRecipe<real_t> gsRecipeAssemblerStokes::getDirichletRecipe  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t> result;
    gsRecipeIngredient<real_t> ingr;

    {
        ingr.setOperator    ( new gsBoundaryL2ScalarOp<real_t>());
        ingr.setTestSpace   ( velocity);
        ingr.setUnknownSpace( velocity);
        ingr.setRule        ( getEliSysWriter() );
        result.add(ingr);
    }

    if( !bc.isHomogeneous() )
    {
        ingr.setOperator    (new gsBoundaryL2TestVecOp<real_t>(f) );
        ingr.setTestSpace   ( velocity);
        ingr.setUnknownSpace( velocity);
        ingr.setRule        ( getEliRhsWriter());
        result.add(ingr);
    }

    if (m_zeroAverage && !bc.isHomogeneous())
    {
        typedef gsAccumulatorWriter<FMatT> AccuW;
        ingr.setOperator(new gsBoundaryNormalDerTestVecOp<real_t>(f));
        ingr.setTestSpace(velocity);
        ingr.setUnknownSpace(pressure);
        ingr.setRule(new gsL2GMapperActive<AccuW>(
                         AccuW( getSysSize()-1,0,m_rhs))
                     );
        result.add(ingr);
    }

    return result;
}

gsRecipe<real_t>    gsRecipeAssemblerStokes::getNeumannRecipe  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryL2TestVecOp<real_t>(f));
    ingr.setTestSpace(velocity);
    ingr.setUnknownSpace(velocity);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    return result;

}

gsIntegrationRule gsRecipeAssemblerStokes::getPatchIntegration     (index_t patch )
{
    gsIntegrationRule result;
    result.rule=gsGaussRule<real_t>(gsVector<index_t>::Constant(m_domain.parDim(),getPatchMaxDegree(patch).maxCoeff()+1));
    result.subdomains=m_space[0]->patchPartition(patch);
    return result;
}


}
