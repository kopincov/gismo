/** @file gsVBasis.h

    @brief A vector-valued basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsEvalUtils.hpp>

namespace gismo {

template <typename T> //, int
class gsVBasis : public gsBasis<T>
{
public:
    /// Shared pointer for gsVBasis
    typedef memory::shared_ptr< gsVBasis > Ptr;

    /// Unique pointer for gsVBasis
    typedef memory::unique_ptr< gsVBasis > uPtr;

protected:
    std::vector<gsBasis<T>*> m_comp;   // gsBasis<T>*
    std::vector<size_t>      m_map;    // map components to basis to avoid recomputing
    std::vector<index_t>     m_shifts; // active shift for each component
    index_t                  m_size;   // sum(m_comp[i].size())
    int                      m_domDim; // m_comp.front().targetDim()
    int                      m_tarDim; // m_comp.size()

    typedef typename std::vector<gsBasis<T>*>::const_iterator const_cIterator;
    typedef typename std::vector<gsBasis<T>*>::iterator       cIterator;

    void initSize()
    {
        index_t s=0;
        for (int i=0; i<m_tarDim;++i)
        {
            m_shifts.push_back(s);
            s+=m_comp[m_map[i]]->size();
        }
        m_size=s;
    }

public:

    gsVBasis(const gsBasis<T> & func, int copies)
    : m_comp(1, func.clone().release()),
          m_map(copies,0),
          m_domDim(func.domainDim()),
          m_tarDim(copies)
    {
        GISMO_ASSERT(func.targetDim()==1, "Cannot create vector from vector valued functions");
        initSize();
    }

    gsVBasis(std::vector<gsBasis<T>*> basis)
        : m_map(basis.size()),
          m_tarDim(basis.size())
    {
        m_domDim=basis[0]->domainDim();
        // build a map to avoid to recompute the same basis more then once
        const_cIterator beg = basis.begin();
        const_cIterator end = basis.end();

        for (const_cIterator cur = beg; cur!=end; ++cur)
        {
            const_cIterator tmp = std::find(beg,cur,*cur);
            if (tmp!=cur)
                m_map[cur-beg]=tmp-beg;
            else
            {
                m_map[cur-beg] = m_comp.size();
                m_comp.push_back(*cur);
            }
        }
        initSize();
        basis.clear();
    }


    gsVBasis(std::vector<gsBasis<T>*> basis, std::vector<index_t> shifts)
        : m_map(basis.size()),
          m_tarDim(basis.size())
    {
        m_domDim=basis[0]->domainDim();
        // build a map to avoid to recompute the same basis more then once
        const_cIterator beg = basis.begin();
        const_cIterator end = basis.end();

        for (const_cIterator cur = beg; cur!=end; ++cur)
        {
            const_cIterator tmp = std::find(beg,cur,*cur);
            if (tmp!=cur)
                m_map[cur-beg]=tmp-beg;
            else
            {
                m_map[cur-beg] = m_comp.size();
                m_comp.push_back(*cur);
            }
        }
        initSize();
        m_shifts=shifts; // shifts are initialized by initSize
        basis.clear();
    }

    ~gsVBasis()
    {
    //    freeAll(m_comp);
    }

    const    gsBasis<T> & component(int i)const {return *m_comp[m_map[i]];}
    gsBasis<T> & component(int i) {return *m_comp[m_map[i]];}

    GISMO_CLONE_FUNCTION(gsVBasis)

    memory::unique_ptr<gsGeometry<T> > makeGeometry(gsMatrix<T> coefs) const
    {
        return memory::unique_ptr<gsGeometry<T> >();
    }

    std::ostream &print(std::ostream &os) const 
    {
        os<<"gsVBasis\n";
        return os;
    }

    virtual short_t domainDim() const {return m_domDim;}
    int dim() const {return m_domDim;}
    virtual short_t targetDim() const {return m_tarDim;}

public:
    virtual unsigned getGeoFlags (unsigned flags) const {return 0;}
    virtual index_t size () const {return m_size;}

    static unsigned extendFlags (unsigned flag)
    {
        if (flag & (NEED_CURL|NEED_DIV) )
            flag |= NEED_DERIV;
        return flag;
    }
    static unsigned maskFlags (unsigned flag)
    {
        const unsigned toZero=NEED_CURL|NEED_DIV;
        return extendFlags(flag) & ~toZero;
    }

    virtual void compute        (const gsMatrix <T> &points, gsFuncData<T> &result) const
    {
        result.flags=extendFlags(result.flags);
        std::vector<gsFuncData<T> > computed(m_comp.size());
        for (size_t c=0; c<computed.size(); ++c)
        {
            computed[c].flags=maskFlags(result.flags);
            computed[c].patchId=result.patchId;
            m_comp[c]->compute(points,computed[c]);
        }
        construct(computed,result);
    }

    virtual void compute        (const gsMapData<T> &geoData, gsFuncData<T> &result) const
    {
        std::vector<gsFuncData<T> > computed(m_comp.size()); //
        for (size_t c=0; c<computed.size(); ++c)
        {
            computed[c].flags=maskFlags(result.flags);
            computed[c].patchId=result.patchId;
            m_comp[c]->compute(geoData.points,computed[c]);
        }
        construct(computed,result);
    }

    void setShifts(std::vector<index_t> newShifts)
    {
        GISMO_ASSERT(newShifts.size()==m_shifts.size(),"dimension mismatch");
        m_shifts=newShifts;
    }

private:
    void construct (const std::vector<gsFuncData<T> > &computed, gsFuncData<T> &result) const
    {
        initResult(computed,result);
        index_t written=0;
        for (int c=0; c<m_tarDim;++c)
        {
            written+=writeComponent(c,computed[m_map[c]],written, result);
        }
        computeDerived(result);
    }

    void initResult(const std::vector<gsFuncData<T> > &computed, gsFuncData<T> &result) const
    {
        const unsigned flags = result.flags;
        result.dim = std::make_pair(m_domDim, m_tarDim);

        index_t activeNum=0;
        index_t activeC=computed[0].actives.cols();

        if (flags & NEED_ACTIVE )
            for (int c=0; c<m_tarDim;++c)
                activeNum+=computed[m_map[c]].actives.rows();
        else if (flags & NEED_VALUE )
            for (int c=0; c<m_tarDim;++c)
                activeNum+=computed[m_map[c]].values[0].rows();
        else if (flags & NEED_DERIV )
            for (int c=0; c<m_tarDim;++c)
                activeNum+=computed[m_map[c]].values[1].rows()/m_domDim;
        else if (flags & NEED_DERIV2 )
        {
            const int dsz = m_domDim*(m_domDim+1) / 2;
            for (int c=0; c<m_tarDim;++c)
                activeNum+=computed[m_map[c]].values[2].rows() / dsz;
        }
        else if (flags & NEED_LAPLACIAN )
            for (int c=0; c<m_tarDim;++c)
                activeNum+=computed[m_map[c]].laplacians.rows();

        if (flags & NEED_ACTIVE )
            result.actives.resize(activeNum,activeC);

        result.values.resize(result.maxDeriv()+1);

        if (flags & NEED_VALUE)
            result.values[0].resize(activeNum*result.dim.second,computed[0].values[0].cols());
        if ( flags & NEED_DERIV )
            result.values[1].resize(activeNum*result.derivSize(),computed[0].values[1].cols());
        if ( flags & NEED_DERIV2 )
            result.values[2].resize(activeNum*result.deriv2Size(),computed[0].values[2].cols());

        // assumes first derivative is set because it is needed to compute the values
        if (flags & NEED_DIV)
            result.divs.resize(activeNum*result.divSize(),result.values[1].cols() );
        if (flags & NEED_CURL)
            result.curls.resize(activeNum*result.derivSize(),result.values[1].cols() );

        // does not assume anything as laplacian is always component wise
        if (flags & NEED_LAPLACIAN)
            result.laplacians.resize(activeNum*result.dim.second, computed[0].laplacians.cols());
    }

    index_t writeComponent (int comp, const gsFuncData<T>& data,
                            index_t written, gsFuncData<T> &result) const
    {
        const unsigned flags = result.flags;
        index_t activeNum = 0;
        result.dim = this->dimensions();

        if (flags & NEED_ACTIVE )
            activeNum=writeActives ( comp, written, data.actives,   result.actives);
        if (flags & NEED_VALUE)
            activeNum=writeBlock (1, comp, written, data.values[0], result.values[0]);
        if ( flags & NEED_DERIV )
            activeNum=writeBlock (result.dim.first, comp, written, data.values[1], result.values[1]);
        if ( flags & NEED_DERIV2 )
            activeNum=writeBlock (result.dim.first*(result.dim.first+1) / 2, comp,written, data.values[2], result.values[2]);
        if (flags & NEED_LAPLACIAN)
            activeNum=writeBlock (1, comp, written, data.laplacians, result.laplacians);
        return activeNum;
    }

    // needed for both unsigned and real_t
    index_t writeBlock (index_t dataUnitSize, index_t comp, index_t startingPos, const gsMatrix<T> &origin, gsMatrix<T> &destination) const
    {
        std::pair<short_t, short_t> info = this->dimensions();
        index_t blocksNum = origin.rows()/dataUnitSize;
        index_t startRow  = startingPos*dataUnitSize*info.second;

        index_t blockTempSize=info.second*dataUnitSize;

        gsMatrix<T> blockTemp;
        blockTemp.setZero(blockTempSize, origin.cols());
        for (index_t b=0; b<blocksNum;++b)
        {
            blockTemp.middleRows(comp*dataUnitSize,dataUnitSize)=origin.middleRows(b*dataUnitSize,dataUnitSize);
            destination.middleRows(startRow+b*blockTempSize,blockTempSize)=blockTemp;
        }
        return blocksNum;
    }



    index_t writeActives(index_t comp, index_t startingPos, const gsMatrix<index_t> &origin, gsMatrix<index_t> &destination) const
    {
        destination.middleRows(startingPos,origin.rows())=origin.array()+m_shifts[comp];
        return origin.rows();
    }

    void computeDerived (gsFuncData<T> &result) const
    {
        result.dim = this->dimensions();
        const unsigned flags=result.flags;
        if ( flags & NEED_DIV )
            convertValue<T>::derivToDiv(result.values[1], result.dim, result.divs);
        if ( flags & NEED_CURL )
            convertValue<T>::derivToCurl(result.values[1], result.dim, result.curls);
    }
};


} // namespace gismo
