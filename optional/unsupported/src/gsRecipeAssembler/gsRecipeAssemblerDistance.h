/** @file gsRecipeAssemblerDistance.h

    @brief DESCRIPTION

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsRecipeAssembler/gsDistanceOperators.h>

namespace gismo {

template<typename Writer1=gsMatrix<real_t>, typename Writer2=gsMatrix<real_t> >
class gsNormRule : public gsL2GMapperRhs<Writer2>
{
protected:
    typedef typename WriterDestType<Writer1>::Destination Wr1;
    typedef typename WriterDestType<Writer1>::Argument    ArgT1;
    typedef typename WriterDestType<Writer2>::Destination Wr2;
    typedef typename WriterDestType<Writer2>::Argument    ArgT2;

    Wr1 m_dest;
protected:
    bool m_writePerCell;
    bool m_writePerActive;

    index_t              m_row;
    std::vector<real_t> *m_perCell;
public:
    gsNormRule( const gsWeightMapper<real_t> &map, index_t row, bool perCell, bool perActive, ArgT1 norm, ArgT2 normPerActive, std::vector<real_t> *perCellV)
        :
          gsL2GMapperRhs<Writer2>(map,normPerActive),
          m_dest(norm),
          m_writePerCell(perCell),
          m_writePerActive(perActive),
          m_row(row),
          m_perCell(perCellV)
    {}

    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<real_t>   &locMat
            )
    {
        if (m_writePerActive)
            gsL2GMapperRhs<gsMatrix<real_t> >::store(activeTest,activeUnknown,locMat);
        if (m_writePerCell)
            m_perCell->push_back(locMat.sum());
        for (int c=0; c < locMat.cols (); ++c)
            m_dest.add(m_row,c,locMat(0,c));
    }
};


/**
 * @brief The gsRecipeDistance class
 * use the recipe assembler to compute the distances between a
 * function specified by coefficients and a discrete space and another function
 *
 * By default it computes the squared L2 distance between each pair
 * (Fi,Gi) where Fi is specified by a physical space Si and coefficients Ci
 * while Gi is provided and defaults to 0.
 *
 * Different norm can be computed by using the method setSeminorms that allows
 * to choose between L2, H1, H2 and W^{k,p} seminorms with k=0,1 and positive p.
 *
 *
 */
class gsRecipeDistance : public gsRecipeAssembler
{
protected:
    std::vector<gsFunction<real_t>*>  m_ref;
    std::vector<gsMatrix<real_t> >    m_coefs;
    std::vector<gsMatrix<real_t> >    m_norms;

    bool                              m_storePerCell;
    bool                              m_storePerFunc;

    std::vector<std::vector<real_t> > m_perCell;
    gsMatrix<real_t>     m_NormM;
    gsMatrix<real_t>     m_perActive;
protected:
    gsRecipeDistance(const gsMultiPatch<real_t>& domain)
        : gsRecipeAssembler(domain)
    {}

   gsRecipeDistance(const gsMultiPatch<real_t>& domain,const std::vector<gsFunction<real_t>*> &ref)
        : gsRecipeAssembler(domain), m_ref(ref)
    {
        // avoid uninitialized errors in this class due to subclasses
        // forgetting to set variables
        m_storePerCell=false;
        m_storePerFunc=false;
        // here the number of spaces is unknown, so it is responsibility of
        // the derived classes to initialize all the data afterwards
        // see initOnce
    }

public:
    /**
     * @brief gsRecipeDistance
     * @param domain is the integration domain
     * @param space  is the list of discrete spaces
     * @param coefs  is the corresponding list of coefficients
     * @param ref    is a vector of functions
     */
    gsRecipeDistance(const gsMultiPatch<real_t> &domain,const Discretization &space, const std::vector<gsMatrix<real_t> > coefs, const std::vector<gsFunction<real_t>*> &ref)
        : gsRecipeAssembler(domain)
    {
        initOnce (space,coefs,&ref);
    }


    gsRecipeDistance(const gsMultiPatch<real_t> &domain,const Discretization &space, const std::vector<gsMatrix<real_t> > coefs)
        : gsRecipeAssembler(domain)
    {
        initOnce (space,coefs,NULL);
    }


    virtual void setOptions(bool perCell, bool perFunc)
    {
        m_storePerCell=perCell;
        m_storePerFunc=perFunc;
    }

    /// each column in the matris is one element
    virtual gsMatrix<real_t> getPerCell() const
    {
        gsMatrix<real_t> result(m_NormM.rows(),m_perCell[0].size());
        for(index_t i=0;i<m_NormM.rows();++i)
        {
            result.row(i)=gsAsConstMatrix<real_t>(m_perCell[i],1,m_perCell[0].size());
        }
        return result;
    }

    virtual gsMatrix<real_t> getDistance (index_t space) const
    {
        index_t firstRow=0;
        for (index_t s=0; s<space;++s)
            firstRow+=m_norms[s].rows();
        return m_NormM.middleRows(firstRow,m_norms[space].rows());
    }

    virtual gsMatrix<real_t> getAllDistances() const
    {
        return m_NormM;
    }

    virtual void setSpace (const Discretization &space)
    {
        gsRecipeAssembler::setSpace(space);
        m_norms.resize(m_space.size());
        m_ref.resize(m_space.size(), NULL);
    }

    virtual void setSeminorms (index_t space, gsMatrix<real_t> norms) {m_norms[space]=norms;}
    virtual void setAllSeminorms (gsMatrix<real_t> norms)
    {
        for(size_t i = 0;i<m_norms.size();++i)
            setSeminorms(i,norms);
    }
    virtual gsMatrix<real_t> getSeminorms (index_t space) const {return m_norms[space];}

    virtual void describeOutput ()
    {
        size_t space=0;
        index_t seminorm=0, row=0;
        while ( space<m_norms.size() )
        {
            if ( !(seminorm<m_norms[space].rows()) )
            {
                seminorm=0;
                ++space;
                continue;
            }
            std::cout<<"row "<<row<<" is the W^{"
                    <<m_norms[space](seminorm,0)<<","<<m_norms[space](seminorm,1)
                    <<"} seminorm of space "<<space<<"\n"<<std::flush;
            ++row;
            ++seminorm;
        }
    }

    virtual void setReference (std::vector<gsFunction<real_t>*> references)
    {
        m_ref=references;
        for(size_t i=0; i < m_space.size(); ++i)
            m_ref.push_back(NULL);
    }

protected:

    /// TODO this breaks computing error on the actives
    virtual SpaceList getPatchSpaces (index_t patch )
    {
        SpaceList result;
        for(unsigned i =0;i<m_space.size();++i)
            result.push_back(m_space[i]->getPatchSpace(patch));
        return result;
    }

    virtual void initSystemMatrices()
    {
        index_t numNorms=0;
        index_t maxCols=1;
        for (size_t sp=0;sp<m_norms.size();++sp)
        {
            numNorms+=m_norms[sp].rows();
            maxCols=math::max(maxCols,m_coefs[sp].cols());
        }
        m_NormM.setZero(numNorms,maxCols);
        if (m_storePerFunc)
            m_perActive.setZero(m_map->getNrOfTargets(),maxCols);
        m_perCell.clear();
        m_perCell.resize(numNorms);
    }

    gsIntegrationRule getPatchIntegration     (index_t patch )
    {
        gsIntegrationRule result;
        result.rule=gsGaussRule<real_t>(getPatchMaxDegree(patch).array()+3);
        result.subdomains=m_space[0]->patchPartition(patch);
        return result;
    }

    gsLocalToGlobalMapper<real_t> *getRule(index_t row)
    {
        return new gsNormRule<>(*m_map, row, m_storePerCell, m_storePerFunc, m_NormM, m_perActive, &m_perCell[row]);
    }

    virtual gsRecipe<real_t> getPatchRecipe (index_t ) // patch)
    {
        gsRecipe<real_t> result;
        index_t pos=0;
        for (size_t sp=0; sp<m_space.size(); ++sp)
        {
            index_t spM=findFirst(sp); // avoid to compute the same space through different evaluators
            for(index_t norm=0; norm<m_norms[sp].rows(); ++norm )
            {
                gsRecipeIngredient<real_t> ingr;

                if ( m_norms[sp](norm,0)==0 )
                {
                    if(m_norms[sp](norm,1)==2)
                        ingr.setOperator(new gsDistL2(m_coefs[sp],m_ref[sp]));
                    else
                        ingr.setOperator(new gsDistW0p(m_norms[sp](norm,1),m_coefs[sp],m_ref[sp]));
                }
                else if (m_norms[sp](norm,0)==1)
                {
                    if(m_norms[sp](norm,1)==2)
                        ingr.setOperator(new gsDistH1(m_coefs[sp],m_ref[sp]));
                    else
                        ingr.setOperator(new gsDistW1p(m_norms[sp](norm,1),m_coefs[sp],m_ref[sp]));
                }
                else if (m_norms[sp](norm,0)==2 && m_norms[sp](norm,1)==2)
                {
                    ingr.setOperator(new gsDistH2(m_coefs[sp],m_ref[sp]));
                }
                else
                    GISMO_ERROR("norm not implemented");
                ingr.setUnknownSpace(spM);
                ingr.setTestSpace(spM);
                ingr.setRule(getRule(pos));
                result.add(ingr);
                ++pos;
            }
        }
        return result;
    }

    /// TODO move to gsRecipeAssembler so everything uses this optimization
    index_t findFirst (index_t sp)
    {
        return std::find(m_space.begin(),m_space.end(),m_space[sp])-m_space.begin();
    }

    virtual gsRecipe<real_t> getBoundaryRecipe (patchSide ) // ps)
    {
        return gsRecipe<real_t>();
    }

protected:

    void initOnce(const Discretization                  & space,
                  const std::vector<gsMatrix<real_t> >    coefs,
                  const std::vector<gsFunction<real_t>*>* ref)
    {
        initInternalData(space,coefs);
        m_storePerCell=false;
        m_storePerFunc=false;
        m_norms.clear();
        gsMatrix<> norms(1,2);
        norms<<0,2;
        for (size_t sp=0; sp<m_space.size(); ++sp)
            m_norms.push_back(norms); // default to L2 norm
        if (ref)
            m_ref=*ref;
        while(m_ref.size()<m_space.size())
            m_ref.push_back(NULL);
    }

    void initInternalData(const Discretization                 & space,
                          const std::vector<gsMatrix<real_t> >   coefs)
    {
        setSpace(space);
        m_coefs.clear();
        for (size_t sp = 0; sp < m_space.size(); ++sp)
        {
            gsWeightMapper<real_t>* map = m_space[sp]->getMapper();
            m_coefs.push_back(map->asMatrix() * coefs[sp]);
            delete map;
        }
    }

};


}
