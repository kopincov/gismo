/** @file gsRecipeAssemblerPevertedStokes.cpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, J. Sogn
**/

#include <iterator>

#include <gsPde/gsStokesPde.h>
#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsRecipeAssembler/gsRecipeAssemblerPervertedStokes.h>
//#include <gsRecipeAssembler/gsRecipeAssemblerStokes.h>


namespace gismo {


template <typename T>
T* Cloner(const T* a)
{
    return (T*)a->clone().release();    // TODO: return uPtr
}

template <typename T>
T* Constifier(T* a)
{
    return a;//const_cast<const T*>(a);
}

void DegreeElevate( gsBasis<real_t>* a);// {a->degreeElevate();} //Defined in gsRecipeAssemblerStokes.cpp

void UniformRefine( gsBasis<real_t>* a) {a->uniformRefine();}

void DegreeReduce( gsBasis<real_t>* a) {a->degreeReduce();}

void ReduceContinuity( gsBasis<real_t>* a) {a->reduceContinuity();}

std::vector<gsPhysicalSpace*> constructPervertedTHSpaces (
        const std::vector<gsBasis<real_t> *>       &bases,
        const gsMultiPatch<real_t>                 &geo,
        std::vector<std::vector<gsBasis<>*> >      *outBasis,
        index_t                              numRefine,
        index_t                              increasePolydegree
        )
{
    const index_t dim=bases[0]->domainDim();
    std::vector<gsPhysicalSpace*>        result(4);

    std::vector<gsBasis<real_t>*>        temp, patchBasisVS, patchBasisPS, patchBasisVM, patchBasisPM;

    {   // build velocity state space
        std::transform(bases.begin(),bases.end(),std::back_inserter(patchBasisVS), Cloner<gsBasis<real_t> >);

        for (index_t k = 0; k < numRefine; ++k)
            std::for_each(patchBasisVS.begin(),patchBasisVS.end(), UniformRefine);

        for (index_t k = 0; k < increasePolydegree; ++k)
            std::for_each(patchBasisVS.begin(),patchBasisVS.end(), DegreeElevate);

        temp.resize(patchBasisVS.size());
        cloneAll(patchBasisVS.begin(),patchBasisVS.end(),temp.begin());
        gsMultiBasis<real_t>                 multipatch (temp,geo);
        gsMapFactoryMultiBasis               mapFactory (multipatch);
        std::vector<gsBasis<real_t>*>        argVec;

        std::transform(patchBasisVS.begin(),patchBasisVS.end(),std::back_inserter(argVec), Constifier<gsBasis<real_t> >);
        gsPhysicalSpaceScalar vel(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
        std::vector<gsPhysicalSpaceScalar*> tmp(dim,&vel);
        result[velocityState]=(new gsPhysicalSpaceVector(tmp));
    }

    {   // build pressure state space
        std::transform(bases.begin(),bases.end(),std::back_inserter(patchBasisPS), Cloner<gsBasis<real_t> >);

        if (increasePolydegree == 0)
            std::for_each(patchBasisPS.begin(),patchBasisPS.end(), DegreeReduce);
        else
        {
            for (index_t k = 0; k < increasePolydegree-1; ++k)
                std::for_each(patchBasisPS.begin(),patchBasisPS.end(), DegreeElevate);
        }

        for (index_t k = 0; k < numRefine; ++k)
            std::for_each(patchBasisPS.begin(),patchBasisPS.end(), UniformRefine);


        temp.resize(patchBasisPS.size());
        cloneAll(patchBasisPS.begin(),patchBasisPS.end(),temp.begin());
        gsMultiBasis<real_t>                 multipatch (temp,geo);
        gsMapFactoryMultiBasis               mapFactory (multipatch);
        std::vector<gsBasis<real_t>*>  argVec;

        std::transform(patchBasisPS.begin(),patchBasisPS.end(),std::back_inserter(argVec), Constifier<gsBasis<real_t> >);
        result[pressureState]=new gsPhysicalSpaceScalar(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
    }

    {   // build velocity multiplier space
        std::transform(bases.begin(),bases.end(),std::back_inserter(patchBasisVM), Cloner<gsBasis<real_t> >);

        for (index_t k = 0; k < numRefine; ++k)
            std::for_each(patchBasisVM.begin(),patchBasisVM.end(), UniformRefine);

        for (index_t k = 0; k < increasePolydegree; ++k)
            std::for_each(patchBasisVM.begin(),patchBasisVM.end(), DegreeElevate);



        //Reduce Continuity
        std::for_each(patchBasisVM.begin(),patchBasisVM.end(), ReduceContinuity);
        std::for_each(patchBasisVM.begin(),patchBasisVM.end(), ReduceContinuity);

        temp.resize(patchBasisVM.size());
        cloneAll(patchBasisVM.begin(),patchBasisVM.end(),temp.begin());
        gsMultiBasis<real_t>                 multipatch (temp,geo);
        gsMapFactoryMultiBasis               mapFactory (multipatch);
        std::vector<gsBasis<real_t>*>  argVec;

        std::transform(patchBasisVM.begin(),patchBasisVM.end(),std::back_inserter(argVec), Constifier<gsBasis<real_t> >);
        gsPhysicalSpaceScalar vel(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
        std::vector<gsPhysicalSpaceScalar*> tmp(dim,&vel);
        result[velocityMultiplier]=(new gsPhysicalSpaceVector(tmp));
    }

    {   // build pressure multiplier space
        std::transform(bases.begin(),bases.end(),std::back_inserter(patchBasisPM), Cloner<gsBasis<real_t> >);

        //Degree refine
        if (increasePolydegree == 0)
            std::for_each(patchBasisPM.begin(),patchBasisPM.end(), DegreeReduce);
        else
        {
            for (index_t k = 0; k < increasePolydegree-1; ++k)
                std::for_each(patchBasisPM.begin(),patchBasisPM.end(), DegreeElevate);
        }

        //h-refine
        for (index_t k = 0; k < numRefine; ++k)
            std::for_each(patchBasisPM.begin(),patchBasisPM.end(), UniformRefine);


        temp.resize(patchBasisPM.size());
        cloneAll(patchBasisPM.begin(),patchBasisPM.end(),temp.begin());
        gsMultiBasis<real_t>                 multipatch (temp,geo);
        gsMapFactoryMultiBasis               mapFactory (multipatch);
        std::vector<gsBasis<real_t>*>  argVec;

        std::transform(patchBasisPM.begin(),patchBasisPM.end(),std::back_inserter(argVec), Constifier<gsBasis<real_t> >);
        result[pressureMultiplier]=new gsPhysicalSpaceScalar(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
    }

    if(outBasis)
    {
        *outBasis= std::vector<std::vector<gsBasis<>*> >(1,patchBasisVS);
        outBasis->push_back(patchBasisPS);
        outBasis->push_back(patchBasisVM);
        outBasis->push_back(patchBasisPM);
    }
    return result;
}

gsRecipeAssemblerPervertedStokes::gsRecipeAssemblerPervertedStokes(const gsPervertedStokesPde<real_t> &pde)
    : gsRecipeAssembler(pde.domain()), m_pde(pde)
{
    m_space.resize(1);
    m_space[0]=NULL;
}

void gsRecipeAssemblerPervertedStokes::collectEliminatedDofs  ()
{
    const gsBoundaryConditions<real_t>::bcContainer container =
        m_pde.boundaryConditions().allConditions();
    typedef gsBoundaryConditions<real_t>::const_iterator iter;
    for (iter bc=container.begin(); bc!=container.end(); ++bc)
    {
        if (bc->type()==condition_type::dirichlet )
        {
            std::vector<index_t> patchDofs = m_space[velocityState]->boundaryDofs(bc->ps,0);
            eliminateDofs(patchDofs,velocityState);
        }
    }
}

gsRecipe<real_t>    gsRecipeAssemblerPervertedStokes::getPatchRecipe     (index_t ) // patch)
{
    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsL2ScalarOp<real_t>());
    ingr.setTestSpace(velocityState);
    ingr.setUnknownSpace(velocityState);
    ingr.setRule(getSysWriter());
    result.add(ingr);

    ingr.setOperator(new gsL2ScalarOp<real_t>());
    ingr.setTestSpace(pressureState);
    ingr.setUnknownSpace(pressureState);
    ingr.setRule(getSysWriter());
    result.add(ingr);

    //Operatior not tested!
    ingr.setOperator(new gsLaplaceOp<real_t>());
    ingr.setTestSpace(velocityMultiplier);
    ingr.setUnknownSpace(velocityState);
    ingr.setRule(getLagrangeMultiplierWriter());
    result.add(ingr);

    ingr.setOperator(new gsDivergenceOp<real_t>());
    ingr.setTestSpace( pressureMultiplier);
    ingr.setUnknownSpace( velocityState ) ;
    ingr.setRule(getLagrangeMultiplierWriter());//B^T and B
    result.add(ingr);

    //Operator not tested!
    ingr.setOperator(new gsGradientOp<real_t>());
    ingr.setTestSpace(velocityMultiplier);
    ingr.setUnknownSpace(pressureState);
    ingr.setRule(getLagrangeMultiplierWriter());//B^T and B
    result.add(ingr);


    //MUST BE ZERO
    //ingr.setOperator(new gsL2TestVecOp<real_t>(*m_pde.rhs()));
    //ingr.setTestSpace(velocityMultiplier);
    //ingr.setUnknownSpace(velocityMultiplier);
    //ingr.setRule(getRhsWriter());
    //result.add(ingr);

    /*if (m_zeroAverage)
    {
        static const gsConstantFunction<real_t> my_one_func(real_t(1),m_domain.dim());

        ingr.setOperator(new gsL2TestOp<real_t>(my_one_func));
        ingr.setTestSpace(pressure);
        ingr.setUnknownSpace(pressure);
        ingr.setRule(new gsL2GMapperRhs<MulW>(*m_map,
                                              MulW(SysW(getSysSize(),getSysSize(),m_sysM,m_rhsMod),0,getFreeLimit()))
                     );
        result.add(ingr);
    }*/
    return result;
}

gsRecipe<real_t>    gsRecipeAssemblerPervertedStokes::getBoundaryRecipe  (patchSide ps)
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

//JS2: not thought out
gsRecipe<real_t> gsRecipeAssemblerPervertedStokes::getDirichletRecipe  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t> result;
    gsRecipeIngredient<real_t> ingr;

    {
        ingr.setOperator    ( new gsBoundaryL2ScalarOp<real_t>());
        ingr.setTestSpace   ( velocityState);
        ingr.setUnknownSpace( velocityState);
        ingr.setRule        ( getEliSysWriter() );
    }
    result.add(ingr);

    if( !bc.isHomogeneous() )
    {
        ingr.setOperator    (new gsBoundaryL2TestVecOp<real_t>(f) );
        ingr.setTestSpace   ( velocityState);
        ingr.setUnknownSpace( velocityState);
        ingr.setRule        ( getEliRhsWriter());

        result.add(ingr);
    }

    /*if (m_zeroAverage && !bc.isHomogeneous())
    {
        typedef gsAccumulatorWriter<FMatT> AccuW;
        ingr.setOperator(new gsBoundaryNormalDerTestVecOp<real_t>(f,s));
        ingr.setTestSpace(velocity);
        ingr.setUnknownSpace(pressure);
        ingr.setRule(new gsL2GMapperActive<AccuW>(
                         AccuW( getSysSize()-1,0,m_rhs))
                     );
        result.add(ingr);
    }*/

    return result;
}

//JS2: not thought out
gsRecipe<real_t>    gsRecipeAssemblerPervertedStokes::getNeumannRecipe  (const boundary_condition<real_t> &bc)
{
    gsFunction<real_t> &f=*bc.function().get();

    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsBoundaryL2TestVecOp<real_t>(f));
    ingr.setTestSpace(velocityState);
    ingr.setUnknownSpace(velocityState);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    return result;

}

gsIntegrationRule gsRecipeAssemblerPervertedStokes::getPatchIntegration     (index_t patch )
{
    gsIntegrationRule result;
    result.rule=gsGaussRule<real_t>(getPatchMaxDegree(patch).array()+2);
    result.subdomains=m_space[0]->patchPartition(patch);
    return result;
}


}
