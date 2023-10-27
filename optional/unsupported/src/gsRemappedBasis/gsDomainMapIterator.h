/** @file gsDomainMapIterator.h

    @brief Iterator for gsDomainMap.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
*/

#pragma once

#include <gsRemappedBasis/gsTensorBasesUtils.h>
#include <gsRemappedBasis/gsDomainMap.h>
#include <gsCore/gsDomainIterator.h>

namespace gismo
{

/** Iterator for the gsDomainMap class.*/
template< int DIM>
class gsDomainMapIterator : public gsDomainIterator<real_t>
{
    typedef typename gsRemTypes<DIM>::tensorBasisT                 tensorBasisT;
    typedef gsFunctionSet<real_t>::Ptr basisPtr;
    typedef typename gsRemTypes<DIM>::knotVectorT                  knotVectorT;
    typedef typename gsRemTypes<DIM>::knotVectorT::const_uiterator breakIter;
public:
    gsDomainMapIterator( const gsDomainMap& map, const std::vector<basisPtr>& _basis )
    :   m_map(map), m_tensBases(_basis)
    {
        m_lower.setZero(DIM);
        m_upper.setZero(DIM);
        first();
    }

    bool next()
    {
        GISMO_ASSERT(m_isGood,"called next on a bad domain iterator");
        ++m_currElem;
        if(m_currElem==m_subDomEnd)
        {
            ++m_currLeaf;
            if( m_currLeaf == m_map.end<gsDomainMap::leaf_FS>() ) // We are at the end.
                m_isGood = false;
            else
                initElemIter();
        }
        if (m_isGood)
            updateCorners();
        return m_isGood;
    }

    bool next(index_t increment)
    {
        for (index_t i = 0; i < increment; i++)
        {
            GISMO_ASSERT(m_isGood,"called next on a bad domain iterator");
            ++m_currElem;
            if(m_currElem==m_subDomEnd)
                {
                    ++m_currLeaf;
                    if( m_currLeaf == m_map.end<gsDomainMap::leaf_FS>() ) // We are at the end.
                        m_isGood = false;
                    else
                        initElemIter();
                }
        }
        if (m_isGood)
            updateCorners();
        return m_isGood;
    }

    void first()
    {
        m_currLeaf = m_map.begin<gsDomainMap::leaf_FS>();
        if( m_currLeaf == m_map.end<gsDomainMap::leaf_FS>() ) // Empty tree.
            m_isGood = false;
        else
        {
            m_isGood = true;
            initElemIter();
            updateCorners();
        }
    }

    void reset () {first();}

    const gsVector<real_t>& lowerCorner() const
    {
        return m_lower;
    }

    const gsVector<real_t>& upperCorner() const
    {
        return m_upper;
    }

    void computeQuadratureRule(const gsVector<int>& numIntNodes)
    {
        GISMO_NO_IMPLEMENTATION;
    }

    const gsMatrix<unsigned>& computeActiveFunctions()
    {
        GISMO_NO_IMPLEMENTATION;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:    
    void updateCorners()
    {
        for( int i =0; i < DIM; ++i )
        {
            m_lower[i] = math::max( *(m_knotBeg[i]+(*m_currElem)(i)  ), m_leafBoundingBox(i,0));
            m_upper[i] = math::min( *(m_knotBeg[i]+(*m_currElem)(i)+1), m_leafBoundingBox(i,1));
        }
    }

    /// Sets m_currElem to the first element of m_currLeaf.
    void initElemIter()
    {
        const index_t lvl = m_currLeaf->data.space;
        tensorBasisT *lvlBasis=asBasis<DIM>(m_tensBases[lvl]);
        m_leafBoundingBox  = m_map.getBoundingBoxOptimized<gsMatrix<real_t,DIM,2> >(m_currLeaf);

        gsMatrix<index_t, DIM,2> limits;
        for( int i = 0; i < DIM; ++i )
        {
            breakIter beg=lvlBasis->knots(i).ubegin();
            breakIter end=lvlBasis->knots(i).uend();
            limits(i,0)=std::upper_bound(beg,end,m_leafBoundingBox(i,0))-beg-1;
            limits(i,1)=std::lower_bound(beg,end,m_leafBoundingBox(i,1))-beg;
            m_knotBeg[i]=beg;
        }
        m_subDomain = gsTensorIndexDomain<DIM>(limits);
        m_currElem  = m_subDomain.begin();
        m_subDomEnd = m_subDomain.end();
    }

protected:
    gsDomainMap                      m_map;
    breakIter                        m_knotBeg[DIM];
    gsTensorIndexDomain<DIM>           m_subDomain;
    gsMatrix<real_t,DIM,2>             m_leafBoundingBox; // maps subdomains are not necessarily
                                                        // aligned with knots so we need this

    std::vector< basisPtr >          m_tensBases;
    gsDomainMap::leaf_iterator_FS    m_currLeaf;  // iterate over the subdomains in the map
    typename gsTensorIndexDomain<DIM>::const_iterator
                                     m_currElem,  // iterate over the elements in the subdomain
                                     m_subDomEnd;  // iterate over the elements in the subdomain


// data needed from parentclass
    gsVector<real_t> m_lower, m_upper;
};

} // namespace gismo
