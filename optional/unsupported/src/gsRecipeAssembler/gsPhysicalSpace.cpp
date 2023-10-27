/**  gsPhysicalSpace.cpp

    Implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, F. Buchegger
    Created on:  2014-10-28
*/

#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsCore/gsTransformedFuncSet.h>
#include <gsCore/gsVectorValuedFunctionSet.h>

namespace gismo{

gsVector<index_t> gsPhysicalSpace::getSpaceDegree() const
{
    size_t i = 0;
    gsBasisEvaluator<real_t> *eval=getEvaluator(i);
    gsVector<unsigned>deg=eval->getDegree();
    delete eval;
    gsVector<index_t> result(deg.rows());
    for(index_t j=0;j<deg.rows();++j)
        result(j)=deg(j);
    i++;
    while(i<nPatches())
    {
        eval=getEvaluator(i);
        deg=eval->getDegree();
        for(index_t j=0;j<deg.rows();++j)
            result(j)=math::max<index_t>(deg(j),result(j));
        delete eval;
        i++;
    }
    return result;
}

gsPhysicalSpaceScalar::gsPhysicalSpaceScalar( const gsMultiBasis<real_t>  &basis, const gsMultiPatch<real_t> &geo, ValueTransformationType transform, const gsMapFactory *factory)
    : m_geo(geo), m_transform(transform)
{
    for(unsigned i = 0;i<basis.nBases();++i)
        m_bases.push_back(const_cast<gsBasis<real_t>*>(&basis.basis(i)));
    const gsMapFactory* tempFactory= factory ? factory : new gsMapFactoryMultiBasis(basis);
    m_mapper=tempFactory->makeMapper();
    initSpace();
    if(!factory)
        delete tempFactory;
}

gsPhysicalSpaceScalar::gsPhysicalSpaceScalar( const BasisContainer &bases, const gsMultiPatch<real_t> &geo, ValueTransformationType transform,const gsWeightMapper<real_t>& mapper)
    : m_bases(bases), m_geo(geo), m_transform(transform)
{
    initSpace();
    m_mapper = new gsWeightMapper<real_t>(mapper);
}

void gsPhysicalSpaceScalar::initSpace()
{
    targetDim=1;
    domainDim=m_bases[0]->domainDim();
    for(unsigned i = 0;i<m_bases.size();++i)
    {
        if (m_transform==INVERSE_COMPOSITION)
            m_space.push_back(spacePtr(new gsGradConformingTFS(*m_bases[i])));
        else if (m_transform==DIV_CONFORMING)
            m_space.push_back(spacePtr(new gsDivConformingTFS(*m_bases[i])));
        else if (m_transform==NO_TRANSFORMATION)
            m_space.push_back(spacePtr(new gsIdentityTFS(*m_bases[i])));
        else m_space.push_back(spacePtr());
    }
}

index_t gsPhysicalSpaceScalar::patchShift(index_t patch) const
{
    index_t shift=0;
    for(index_t i=0;i<patch;++i)
        shift+=patchSize(i);
    return shift;
}

gsBasisEvaluator<real_t> *gsPhysicalSpaceScalar::getEvaluator(index_t patch) const
{
    gsBasisEvaluator<real_t>* evaluator=makeBasisEvaluator(*m_bases[patch],0,&m_geo.patch(patch),m_transform);
    evaluator->addActiveShift(patchShift(patch));
    return evaluator;
}

gsPhysicalSpaceScalar::spacePtr gsPhysicalSpaceScalar::getPatchSpace(index_t patch) const
{
    gsTransformedFuncSet<real_t> *physSpace=dynamic_cast<gsTransformedFuncSet<real_t> *>(m_space[patch].get());
    physSpace->activeShift=patchShift(patch);
    return m_space[patch];
}

gsVector<index_t> gsPhysicalSpaceScalar::getSpaceDegree() const
{
    index_t dim=m_bases[0]->domainDim();
    gsVector<index_t> deg(dim);
    for (int i=0; i<dim;++i)
        deg(i)=m_bases[0]->degree(i);
    return deg;
}

std::vector<index_t> gsPhysicalSpaceScalar::boundaryDofs(patchSide ps, index_t offset) const
{
    const index_t pShift=patchShift(ps.patch);
    std::vector<index_t> result;
    gsMatrix<index_t> matrix;
    matrix=m_bases[ps.patch]->boundaryOffset(ps.side(),offset);
    for (int i=0; i< matrix.rows();++i)
        result.push_back(matrix.at(i)+pShift);
    return result;
}







std::vector<index_t> getPositions(index_t patchId, index_t targetDim, index_t numPatches)
{
    std::vector<index_t> result;
    for(index_t i = 0; i<targetDim; ++i)
        result.push_back(patchId+i*numPatches);
    return result;
}


gsPhysicalSpaceVector::gsPhysicalSpaceVector( std::vector<gsPhysicalSpaceScalar*> spaces )
    : m_geo(spaces[0]->m_geo),m_transform(spaces[0]->m_transform),m_dim(spaces.size())
{
    std::vector<gsWeightMapper<real_t>*> mappers;
    std::vector<index_t> shifts;
    for(unsigned i = 0;i<spaces.size();++i)
        mappers.push_back(spaces[i]->m_mapper);
    m_mapper = combineMappers(mappers,shifts,true);
    for(unsigned i = 0;i<spaces.size();++i)
        for(unsigned j = 0;j<spaces[i]->m_bases.size();++j)
            m_bases.push_back(spaces[i]->m_bases[j]);
    initSpace();
}

gsBasisEvaluator<real_t> *gsPhysicalSpaceVector::getEvaluator(index_t patch) const
{
    std::vector<index_t>  pos = getPositions(patch,m_dim,nPatches());
    std::vector<index_t> shifts= computeComponentShifts(patch);
    BasisContainer tempBases;
    for(size_t i = 0;i<pos.size();++i)
        tempBases.push_back(m_bases[pos[i]]);
    gsBasisEvaluator<real_t>* evaluator=makeBasisEvaluator(tempBases,shifts,0,&m_geo.patch(patch),m_transform);
    return evaluator;
}

std::vector<index_t> gsPhysicalSpaceVector::boundaryDofs(patchSide ps, index_t offset) const
{
    gsBasisEvaluator<real_t>* evaluator=getEvaluator(ps.patch);
    std::vector<index_t> result = evaluator->getBoundary(ps.side(),offset);
    delete evaluator;
    return result;
}

gsPhysicalSpace::spacePtr gsPhysicalSpaceVector::getPatchSpace(index_t patch) const
{
    gsTransformedFuncSet<real_t> *physSpace=dynamic_cast<gsTransformedFuncSet<real_t> *>(m_space[patch].get());
    physSpace->activeShift=0;
    static_cast<gsVBasis<real_t>*>(m_vecBas[patch].get())->setShifts(computeComponentShifts(patch));
    return m_space[patch];
}


std::vector<index_t> gsPhysicalSpaceVector::computeComponentShifts(index_t patch) const
{
    std::vector<index_t> pos = getPositions(patch,m_dim,nPatches());
    for(index_t i = 0;i<m_dim;++i)
        pos.push_back(patch+i*nPatches());
    size_t index=0,cur_shift=0;
    std::vector<index_t> shifts;
    for(size_t i = 0;i<m_bases.size();++i)
    {
        if( index<pos.size() && i==static_cast<size_t>(pos[index]) )
        {
            shifts.push_back(cur_shift);
            index++;
        }
        cur_shift+=m_bases[i]->size();
    }
    return shifts;
}

void gsPhysicalSpaceVector::initSpace()
{
    targetDim=m_bases.size()/nPatches();
    domainDim=m_bases[0]->domainDim();

    m_vecBas.resize(nPatches());
    for (size_t p=0; p<nPatches();++p)
    {
        std::vector<index_t> pos = getPositions(p,m_dim,nPatches());
        BasisContainer tempBases;
        for(size_t i = 0;i<pos.size();++i)
            tempBases.push_back(m_bases[pos[i]]);
        m_vecBas[p]=spacePtr(new gsVBasis<real_t>(tempBases,computeComponentShifts(p)));
        switch(m_transform)
        {
        case INVERSE_COMPOSITION:
            m_space.push_back(spacePtr(new gsGradConformingTFS(*m_vecBas[p])));
            break;
        case DIV_CONFORMING:
            m_space.push_back(spacePtr(new gsDivConformingTFS(*m_vecBas[p])));
            break;
        case NO_TRANSFORMATION:
            m_space.push_back(spacePtr(new gsIdentityTFS(*m_vecBas[p])));
            break;
        default:
            m_space.push_back(spacePtr());
        }
    }
}

}
