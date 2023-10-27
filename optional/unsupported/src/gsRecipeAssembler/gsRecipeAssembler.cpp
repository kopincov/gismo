/** @file gsRecipeAssembler.cpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, F. Buchegger
*/

#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsCore/gsTransformedFuncSet.h>
#include <iterator> // needed for std::back_inserter on mvc compilers

namespace gismo{


gsRecipeAssembler::~gsRecipeAssembler()
{
    delete m_map;
    freeAll(std::unique(m_space.begin(), m_space.end()),m_space.end());
}

const char checkSpaceIdErrorMsg[]="space unknown: either you are setting eliminated dofs for a non existent space or you tried to eliminate dofs before setting the spaces.";

void gsRecipeAssembler::eliminateDofs(const std::vector<index_t> &dofs, int spaceId)
{
    GISMO_ASSERT(spaceId<static_cast<index_t>(m_eliminatedDofs.size()),checkSpaceIdErrorMsg);
    freeDofs(dofs,spaceId);
    std::vector<index_t> &dest = spaceId==target ? m_eliminatedTarget : m_eliminatedDofs[spaceId];
    dest.insert(dest.end(),dofs.begin(),dofs.end());
}

void gsRecipeAssembler::freeDofs(const std::vector<index_t> &other, int spaceId)
{
    GISMO_ASSERT(spaceId<static_cast<index_t>(m_eliminatedDofs.size()),checkSpaceIdErrorMsg);
    std::vector<index_t> &dest = spaceId==target ? m_eliminatedTarget : m_eliminatedDofs[spaceId];
    for (size_t i=0; i<other.size();++i )
        dest.erase(std::remove(dest.begin(),dest.end(),other[i]),dest.end());
}

void gsRecipeAssembler::assemble ()
{
    this->init                   ();
    this->assembleVolumes        ();
    this->assembleBoundaries     ();
    this->postProcess            ();
}

gsMatrix<real_t>  gsRecipeAssembler::reconstructSolution(index_t space, const gsMatrix<real_t> &sysSol, const gsMatrix<real_t> &elimDof ) const
{
    gsMatrix<real_t>recSol(m_map->getNrOfTargets(),sysSol.cols());
    recSol.topRows(getFreeLimit())=sysSol.topRows(getFreeLimit());
    if(getElimSize()>0)
        recSol.bottomRows(getElimSize())=elimDof;
    return (m_permutation.inverse()*recSol).middleRows(m_shiftsTarget[space],m_shiftsTarget[space+1]-m_shiftsTarget[space]);
}

void gsRecipeAssembler::setSpace (const Discretization &space)
{
    freeAll(m_space);
    m_space=space;
    m_eliminatedDofs.clear();
    m_eliminatedDofs.resize(m_space.size());
    m_eliminatedTarget.clear();
}


void gsRecipeAssembler::reset()
{
    m_eliminatedDofs.clear();
    m_eliminatedDofs.resize(m_space.size());
    m_eliminatedTarget.clear();
}

void gsRecipeAssembler::init()
{
    this->initMappers            ();
    this->remapEliminatedDOFS    ();
    this->initSystemMatrices     ();
    m_map->optimize();
}

void gsRecipeAssembler::initMappers()
{
    if(m_map) delete m_map;

    m_shiftsTarget.clear();
    index_t targerShift=0;
    std::vector<gsWeightMapper<real_t>*> mappers;
    for(size_t i=0;i<m_space.size();++i)
    {
        m_shiftsTarget.push_back(targerShift);
        mappers.push_back(m_space[i]->getMapper());
        targerShift+=mappers.back()->getNrOfTargets();
    }
    m_shiftsTarget.push_back(targerShift);
    m_map = combineMappers(mappers,m_shiftsSource);
    freeAll(mappers);
}

// TODO use a lambda expression and c++11
struct myPrivateSum
{
    index_t amount;
    myPrivateSum(index_t howmuch)
        : amount(howmuch)
    {}
    index_t operator() (index_t n) {return n+amount;}
};

void gsRecipeAssembler::remapEliminatedDOFS()
{
    m_map->optimize(gsWeightMapper<real_t>::optSourceToTarget);
    std::vector<index_t> shifted;
    std::vector<index_t> eliminated;
    for (size_t space=0; space<m_space.size();++space)
    {
        std::sort(m_eliminatedDofs[space].begin(),m_eliminatedDofs[space].end());
        std::transform(m_eliminatedDofs[space].begin(),m_eliminatedDofs[space].end(), std::back_inserter(shifted),myPrivateSum(m_shiftsSource[space]));
        m_map->sourceToTarget(shifted, eliminated);
        eliminateDofs(eliminated,target);
        shifted.clear();
        eliminated.clear();
    }
    std::sort(m_eliminatedTarget.begin(),m_eliminatedTarget.end());
    reorderMapperTarget(*m_map,m_eliminatedTarget,&m_permutation);
}

void gsRecipeAssembler::initSystemMatrices()
{
    const index_t eliSize   = getElimSize();
    const index_t size      = getSysSize();
    const index_t sizeRhs   = getRhsDim();


    // init m_sysMatrix
    m_sysM.resize(size,size);
    if(size>0)
        m_sysM.reserve(getSysEstimatedEntries());

    // init m_sysMatrix
    m_rhsMod.resize(size,eliSize);
    if(eliSize>0)
        m_rhsMod.reserve(getSysRhsModEstimatedEntries());

    // init m_elimM
    m_elimM.resize(eliSize,eliSize);
    if(eliSize>0)
        m_elimM.reserve(getElimEstimatedEntries());

    // init m_elimRhs
    m_elimRhs.setZero(eliSize,sizeRhs);

    // init m_rhs
    m_rhs.setZero(size,sizeRhs);
}

gsVector<index_t> gsRecipeAssembler::getSysEstimatedEntries()
{
    index_t estimatedOverlap=0;
    for (size_t patch=0; patch<m_domain.nPatches();++patch)
        estimatedOverlap=std::max<index_t>(getPatchMaxOverlap(patch),estimatedOverlap);
    estimatedOverlap*=m_space.size();
    estimatedOverlap=math::min(estimatedOverlap,getSysSize());
    gsVector<index_t> reserve;
    reserve.setConstant(getSysSize(),estimatedOverlap);
    return reserve;
}

gsVector<index_t> gsRecipeAssembler::getSysRhsModEstimatedEntries()
{
    index_t estimatedOverlap=0;
    for (size_t patch=0; patch<m_domain.nPatches();++patch)
        estimatedOverlap=std::max<index_t>(getPatchMaxOverlap(patch),estimatedOverlap);
    estimatedOverlap*=m_space.size();
    estimatedOverlap=math::min(estimatedOverlap,getSysSize());
    gsVector<index_t> reserve;
    reserve.setConstant(getElimSize(),estimatedOverlap);
    return reserve;
}

gsVector<index_t> gsRecipeAssembler::getElimEstimatedEntries()
{
    index_t estimatedOverlap=0;
    for (size_t patch=0; patch<m_domain.nPatches();++patch)
        estimatedOverlap=std::max<index_t>(getPatchMaxOverlap(patch),estimatedOverlap);
    estimatedOverlap*=m_space.size();
    estimatedOverlap=math::min(estimatedOverlap,getElimSize());
    gsVector<index_t> reserve;
    reserve.setConstant(getElimSize(),estimatedOverlap);
    return reserve;
}

void gsRecipeAssembler::postProcess()
{
    m_elimM.makeCompressed();
    m_sysM.makeCompressed();
    m_rhsMod.makeCompressed();
}

gsIntegrationRule gsRecipeAssembler::getPatchIntegration     (index_t patch )
{
    gsIntegrationRule result;
    result.rule=gsGaussRule<real_t>(getPatchMaxDegree(patch).array()+1);
    result.subdomains=m_space[0]->patchPartition(patch);
    return result;
}

gsIntegrationRule gsRecipeAssembler::getBoundaryIntegration  (patchSide ps)
{
    gsIntegrationRule result;
    index_t patch=ps.patch;
    gsVector<index_t> numP=getPatchMaxDegree(patch).array()+1;
    numP((ps.side()-1)/2)=1;
    result.rule=gsGaussRule<real_t>(numP);
    result.subdomains=m_space[0]->patchPartition(patch,ps.side());
    return result;
}


gsRecipeAssembler::SpaceList gsRecipeAssembler::getPatchSpaces (index_t patch )
{
    SpaceList result;
    for(unsigned i =0;i<m_space.size();++i)
    {
        result.push_back(m_space[i]->getPatchSpace(patch));
        dynamic_cast<gsTransformedFuncSet<real_t> *>(result.back().get())->activeShift+=m_shiftsSource[i];
    }
    return result;
}

gsVector<index_t> gsRecipeAssembler::getPatchMaxDegree(index_t ) // patch)
{
    gsVector<index_t> temp;
    gsVector<index_t> degree;
    degree.setZero(m_domain.dim());
    for (size_t space=0; space<m_space.size();++space)
    {
        temp=m_space[space]->getSpaceDegree();
        // stupid EIGEN this does not work
        //(degree.array().max(temp.array().cast<index_t>());
        for (index_t i=0; i<temp.size(); ++i)
            degree(i)=math::max<index_t>(temp(i),degree(i));
    }
    return degree;
}


gsLocalToGlobalMapper<real_t> *gsRecipeAssembler::getSysWriter(real_t coef)
{
    return new gsL2GMapper<SysWC>(
                *m_map,
                *m_map,
                SysWC(
                    SysW(getFreeLimit(),getFreeLimit(),m_sysM,m_rhsMod),
                    coef)
                );
}

gsLocalToGlobalMapper<real_t> *gsRecipeAssembler::getLagrangeMultiplierWriter(real_t coef)
{
    return new gsL2GMapper<SymWC>(
                *m_map,
                *m_map,
                SymWC(
                    SymW(
                        SysW(getFreeLimit(),getFreeLimit(),m_sysM,m_rhsMod)),
                    coef)
                );
}


gsLocalToGlobalMapper<real_t> *gsRecipeAssembler::getSysWriter()
{
    return new gsL2GMapper<SysW>(
                *m_map,
                *m_map,
                SysW(getFreeLimit(),getFreeLimit(),m_sysM,m_rhsMod)
                );
}

gsLocalToGlobalMapper<real_t> *gsRecipeAssembler::getLagrangeMultiplierWriter()
{
    return new gsL2GMapper<SymW>(
                *m_map,
                *m_map,
                SymW(SysW(getFreeLimit(),getFreeLimit(),m_sysM,m_rhsMod))
                );
}

gsLocalToGlobalMapper<real_t> *gsRecipeAssembler::getRhsWriter()
{
    return new gsL2GMapperRhs<RhsW>(
                *m_map,
                RhsW(getFreeLimit(),INT_MAX,m_rhs,gsNullWriter<real_t>::unique_instance)
                );
}
gsLocalToGlobalMapper<real_t> *gsRecipeAssembler::getEliSysWriter()
{
    return new  gsL2GMapper<EliSysW>(
                *m_map,
                *m_map,
                EliSysW(m_elimM,getFreeLimit(),getFreeLimit())
                );
}
gsLocalToGlobalMapper<real_t> *gsRecipeAssembler::getEliRhsWriter()
{
    return new gsL2GMapperRhs<EliRhsW>(
                *m_map,
                EliRhsW(m_elimRhs,getFreeLimit(),0)
                );
}


void gsRecipeAssembler::assembleVolumes ()
{
    for (size_t patch=0; patch<m_domain.nPatches(); ++patch)
    {
        SpaceList spList;
        spList=getPatchSpaces(patch);
        gsIntegrationRule intRule = this->getPatchIntegration(patch);
        gsRecipe<real_t>   recipe = this->getPatchRecipe(patch);
        recipe.assemble(intRule.rule,*intRule.subdomains,m_domain[patch],spList);
    }
}


void gsRecipeAssembler::assembleBoundaries()
{
    gsMultiPatch<real_t>::const_biterator bs=m_domain.bBegin();
    gsMultiPatch<real_t>::const_biterator be=m_domain.bEnd();
    for (; bs!=be; ++bs)
    {
        index_t patch=bs->patch;
        SpaceList           spList  = this->getPatchSpaces(patch);
        gsIntegrationRule intRule = this->getBoundaryIntegration(*bs);
        gsRecipe<real_t>  recipe  = this->getBoundaryRecipe(*bs);
        recipe.assemble(intRule.rule,*intRule.subdomains,m_domain[patch],spList, *bs);
    }
}


}
