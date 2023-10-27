/**  gsRecipe.hpp

    Implementation of the gsRecipeAssembler framework

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
    Created on:  2014-10-28
*/


#include  <gsRecipeAssembler/gsRecipeAssembler.h>
#include  <gsRecipeAssembler/gsPDEOperators.h>
#include  <gsCore/gsGeometry.h>
#include  <gsCore/gsDomainIterator.h>
#include  <gsAssembler/gsQuadRule.h>
#include  <gsCore/gsFunctionSet.h>
#include  <gsCore/gsTransformedFuncSet.h>

namespace gismo {


template<typename T>
void gsRecipeIngredient<T>::operator= (gsRecipeIngredient other)
{
    std::swap(m_op  , other.m_op  );
    std::swap(m_rule, other.m_rule);

    this->m_tSpace=other.m_tSpace;
    this->m_uSpace=other.m_uSpace;
}

template<typename T>
gsRecipeIngredient<T>::gsRecipeIngredient()
    : m_tSpace(0), m_uSpace(0)
{}

template<typename T>
void gsRecipeIngredient<T>::setOperator (gsBilinearOp<T> *op)
{
    GISMO_ASSERT( op != NULL, "Received NULL operator.");
    m_op = opHandle(op);
}

template<typename T>
void gsRecipeIngredient<T>::setOperator (const opHandle & op)
{
    GISMO_ASSERT( op , "Received NULL operator.");
    m_op = op;
}

template<typename T>
void gsRecipeIngredient<T>::setRule (gsLocalToGlobalMapper<T>   *rule)
{
    GISMO_ASSERT( rule != NULL, "Received NULL rule.");
    m_rule =  ruleHandle(rule);
}

template<typename T>
void gsRecipeIngredient<T>::setRule (const ruleHandle & rule)
{
    GISMO_ASSERT( rule, "Received NULL rule.");
    m_rule = rule;
}

template<typename T>
void gsRecipeIngredient<T>::setTestSpace (unsigned tSpace)
{
    m_tSpace=tSpace;
}

template<typename T>
void gsRecipeIngredient<T>::setUnknownSpace (unsigned uSpace)
{
    m_uSpace=uSpace;
}

template <typename T>
gsRecipe<T>::gsRecipe() { }

template <typename T>
gsRecipe<T>::gsRecipe(size_t size) : m_data(size) { }

template <typename T>
gsRecipeIngredient<T>& gsRecipe<T>::operator[](size_t i)
{
    return m_data[i];
}

template <typename T>
const gsRecipeIngredient<T>& gsRecipe<T>::operator[](size_t i) const
{
    return m_data[i];
}

template <typename T>
size_t gsRecipe<T>::size() const {return m_data.size();}

template <typename T>
void gsRecipe<T>::operator= (gsRecipe<T> other)
{
    std::swap(m_data, other.m_data );
}

template <typename T>
void gsRecipe<T>::append (const gsRecipe &other)
{
    m_data.insert( m_data.end(), other.m_data.begin(), other.m_data.end());
}

template <typename T>
void gsRecipe<T>::add (const gsRecipeIngredient<T> & ingredient)
{
    m_data.push_back(ingredient);
}


template <typename T>
void gsRecipe<T>::add(opHandle op, ruleHandle writer,unsigned trial, unsigned test)
{
    m_data.push_back(gsRecipeIngredient<T>(op,writer,trial,test));
}


template <typename T>
void gsRecipe<T>::resize (size_t new_size)
{
    m_data.resize(new_size);
}


template <typename T>
void gsRecipe<T>::assemble(const gsQuadRule<T>               & quadRule,
                           gsDomainIterator<T>               & element,
                           const gsGeometry<T>               & geometry,
                           spaceList                           spaces,
                           patchSide                           pside) const
{
    if (size()==0)
        return; // avoid useless loop on the elements

    // data per element
    std::vector<gsFuncData<T> >   basVal(spaces.size());
    gsMapData<T>                  geoVal;
    geoVal.side=pside.side();
    gsPointMapData<T>             ptGeoVal(geoVal);

    // initialize
    for (unsigned step=0; step<size(); ++step )
    {
        opHandle  op       = m_data[step].getOp();
        unsigned  tSpaceId = m_data[step].tSpaceId();
        unsigned  uSpaceId = m_data[step].uSpaceId();

        basVal[tSpaceId].addFlags(op->testSpaceNeeds()    | SAME_ELEMENT| NEED_ACTIVE);
        basVal[uSpaceId].addFlags(op->unknownSpaceNeeds() | SAME_ELEMENT| NEED_ACTIVE);
        basVal[uSpaceId].patchId=pside.patch;
        basVal[tSpaceId].patchId=pside.patch;

        geoVal.addFlags(op->geometryNeeds());

        const gsTransformedFuncSet<T> * transSpace;
        transSpace = dynamic_cast<gsTransformedFuncSet<T> *>(spaces[tSpaceId].get());
        if (transSpace) geoVal.addFlags(transSpace->getGeoFlags(basVal[tSpaceId].flags) );
        transSpace = dynamic_cast<gsTransformedFuncSet<T> *>(spaces[uSpaceId].get());
        if (transSpace) geoVal.addFlags(transSpace->getGeoFlags(basVal[uSpaceId].flags) );
    }
    geoVal.addFlags(SAME_ELEMENT | (pside.side()==boundary::none ? NEED_MEASURE: NEED_OUTER_NORMAL) );

    gsVector<T> quWeights;

    // loop over the elements
    for (element.reset(); element.good(); element.next ())
    {
        // map quadrature points
        quadRule.mapTo (element.lowerCorner() ,element.upperCorner(), geoVal.points,quWeights );
        geometry.computeMap(geoVal);

        if (pside.side()==boundary::none)
            quWeights.array()*= geoVal.measures.transpose().array();
        else
        {
            quWeights.array()*=geoVal.outNormals.colwise().norm().transpose().array();
            geoVal.outNormals.colwise().normalize();
        }

        for(unsigned spId=0; spId< spaces.size (); ++spId)
            dynamic_cast<gsTransformedFuncSetImpl<T> *>(spaces[spId].get())->compute(geoVal,basVal[spId]);

        gsMatrix<T> locMat;
        // compute each step
        for (unsigned step=0; step<size(); ++step )
        {
            opHandle    op          = m_data[step].getOp();
            ruleHandle  rule        = m_data[step].getRule();
            unsigned    tSpaceId    = m_data[step].tSpaceId();
            unsigned    uSpaceId    = m_data[step].uSpaceId();
            gsPointFuncData<T> tSpaceData(basVal[tSpaceId]);
            gsPointFuncData<T> uSpaceData(basVal[uSpaceId]);

            unsigned rows  = tSpaceData.actives().rows();
            unsigned cols  = uSpaceData.actives().size();
            op->outputSize(rows, cols);
            locMat.setZero(rows, cols);

            for (int nodeId=0; nodeId<quadRule.numNodes (); ++nodeId)
            {
                tSpaceData.setPoint(nodeId);
                uSpaceData.setPoint(nodeId);
                ptGeoVal.setPoint(nodeId);

                op->pointEval( tSpaceData, uSpaceData, ptGeoVal, gsRecipeAccumulator<T>(locMat,quWeights[nodeId]) );
            } // quadrature points loop
            rule->store( tSpaceData.actives(), uSpaceData.actives(), locMat);
        } // ingredient loop
    } // element loop
}




template <typename T>
void gsRecipe<T>::assemble(const gsMultiPatch<T>              &domain,
                           spaceList                          spaces,
                           gsVector<index_t>                      numGaussPoints,
                           const gsElementList<T>             mesh
                           ) const
{
    int numQuad=0;
    int dir=0;

    for (size_t mp=0; mp<mesh.size();++mp)
    {
        size_t  patch = mesh[mp].second;
        boxSide side  = mesh[mp].first->side();
        if (side!=boundary::none)
        {
            dir = side.direction();
            numQuad = numGaussPoints(dir);
            numGaussPoints(dir) = 1;
        }
        gsGaussRule<T> quad(numGaussPoints);
        assemble(quad, *mesh[mp].first, domain.patch(patch),spaces,patchSide(patch,side));
        if (side!=boundary::none)
            numGaussPoints(dir) = numQuad; // restore original quadrature
    }
}


} // namespace gismo
