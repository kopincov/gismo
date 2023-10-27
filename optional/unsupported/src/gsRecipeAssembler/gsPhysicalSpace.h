/**  gsPhysicalSpace.h

    Abstraction of a space on the physical domain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, F. Buchegger
    Created on:  2014-10-28
*/


#pragma once

#include <gsCore/gsBoundary.h>
#include <gsMapUtils/gsMapFactoryMultiBasis.h>
#include <gsCore/gsBasisEvaluator.h>
#include <gsCore/gsMultiPatch.h>
#include <gsMapUtils/gsMapperUtils.h>
#include <gsCore/gsDomainIterator.h>

namespace gismo{

/**
 * @brief The gsPhysicalSpace interface
 */
class GISMO_EXPORT gsPhysicalSpace
{
public:
    typedef memory::shared_ptr<gsFunctionSet<real_t> > spacePtr;
public:
    /**
     * @brief  getMapper
     * @return a pointer to a gsWeightMapper that maps patch local degrees of freedom
     *         to space global degrees of freedom
     */
    virtual gsWeightMapper<real_t>   *getMapper() const = 0;
    /**
     * @brief  getEvaluator
     * @param  patch
     * @return a basis evaluator that evaluates the basis functions on a patch
     */
    virtual gsBasisEvaluator<real_t> *getEvaluator(index_t patch) const = 0;
    /**
     * @brief  getPatchSpace
     * @param  patch
     * @returna shared pointer to transformed function set
     */
    virtual spacePtr                  getPatchSpace(index_t patch) const = 0;
    /**
     * @brief   boundaryDofs
     * @param   ps
     * @param   offset
     * @return  a vector containing the degrees of freedom that are active on the
     *          boundary side specified by ps with derivatives up to order offset-1
     */
    virtual std::vector<index_t>      boundaryDofs(patchSide ps, index_t offset) const = 0;
    /**
     * @brief  patchPartition
     * @param  patch
     * @param  side
     * @return a pointer to a domain iterator for the given patch or patch side
     */
    virtual gsDomainIterator<real_t> *patchPartition(index_t patch,boxSide side=0) const = 0;
    /**
     * @brief  getSpaceOverlap
     * @return the maximum number of functions active on the same point or a guess based on degree
     */
    virtual index_t                   getSpaceOverlap() const = 0;
    /**
     * @brief getSpaceDegree
     * @return the suggested direction degree for using with tensor product Gauss quadrature
     */
    virtual gsVector<index_t>         getSpaceDegree() const;

    /// gives back the number of patches
    virtual size_t nPatches() const = 0;
public:
    virtual ~gsPhysicalSpace() {}
    short_t targetDim;
    short_t domainDim;

};


/**
 * @brief The gsPhysicalSpaceScalar class
 *        an implementation of gsPhysicalSpace based on MultiBasis and a gsMapFactory
 */
class GISMO_EXPORT gsPhysicalSpaceScalar : public gsPhysicalSpace
{
    typedef std::vector<gsBasis<real_t> *> BasisContainer;
    typedef std::vector<spacePtr>                SpaceContainer;
public:
          BasisContainer          m_bases;
          SpaceContainer          m_space;
    const gsMultiPatch<real_t>   &m_geo;
          gsWeightMapper<real_t> *m_mapper;
    const ValueTransformationType m_transform;
public:
    gsPhysicalSpaceScalar( const gsMultiBasis<real_t>  &basis, const gsMultiPatch<real_t> &geo, ValueTransformationType transform, const gsMapFactory *factory=NULL);

    gsPhysicalSpaceScalar( const BasisContainer &bases, const gsMultiPatch<real_t> &geo, ValueTransformationType transform,const gsWeightMapper<real_t>& mapper);

    ~gsPhysicalSpaceScalar() { delete m_mapper; }

private:
    index_t patchSize(index_t patch) const
    { return m_bases[patch]->size(); }

    index_t patchShift(index_t patch) const;

    void initSpace();
public:
    virtual size_t                    nPatches() const {return m_bases.size();}
    virtual gsWeightMapper<real_t>   *getMapper() const                                {return new gsWeightMapper<real_t>(*m_mapper);}

    virtual gsBasisEvaluator<real_t> *getEvaluator(index_t patch) const;
    virtual spacePtr                  getPatchSpace(index_t patch) const;
    virtual gsDomainIterator<real_t> *patchPartition(index_t patch,boxSide side) const {return m_bases[patch]->makeDomainIterator(side).release();}
    virtual index_t                   getSpaceOverlap() const                          {return (getSpaceDegree().array()*2+1).prod();}
    virtual gsVector<index_t>         getSpaceDegree() const;
    virtual std::vector<index_t>      boundaryDofs(patchSide ps, index_t offset) const;
};

/**
 * @brief The gsSimpleSpace class
 *        an implementation of gsPhysicalSpace based on MultiBasis and a gsMapFactory
 */
class GISMO_EXPORT gsPhysicalSpaceVector : public gsPhysicalSpace
{
    typedef std::vector<gsBasis<real_t> *> BasisContainer;
    typedef std::vector<spacePtr>          SpaceContainer;
    void initSpace();
public:
          BasisContainer          m_bases;
          SpaceContainer          m_vecBas;
          SpaceContainer          m_space;
    const gsMultiPatch<real_t>   &m_geo;
          gsWeightMapper<real_t>         *m_mapper;
    const ValueTransformationType m_transform;
    const index_t                 m_dim;
public:
    gsPhysicalSpaceVector( std::vector<gsPhysicalSpaceScalar*> spaces );

    ~gsPhysicalSpaceVector() { delete m_mapper; }
public:
    virtual size_t                    nPatches() const {return m_bases.size()/m_dim;}
    virtual gsWeightMapper<real_t>   *getMapper() const                                {return new gsWeightMapper<real_t>(*m_mapper);}
    virtual gsBasisEvaluator<real_t> *getEvaluator(index_t patch) const;
    virtual spacePtr                  getPatchSpace(index_t patch) const;
    virtual gsDomainIterator<real_t> *patchPartition(index_t patch,boxSide side) const {return m_bases[patch]->makeDomainIterator(side).release();}
    virtual index_t                   getSpaceOverlap() const                          {return (getSpaceDegree().array()*2+1).prod();}
    virtual std::vector<index_t>      boundaryDofs(patchSide ps, index_t offset) const;
private:
    std::vector<index_t> computeComponentShifts(index_t patch) const;

};

}
