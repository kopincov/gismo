/** @file gsBoxList.cpp

    @brief Bodies of non-templated functions from of gsBoxList.h.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
*/

#include <gsRemappedBasis/gsBoxList.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

namespace gismo
{
// Getters
gsBoxList::basisIdT gsBoxList::maxId (const gsMatrix<real_t> &bbox) const
{
    basisIdT M=0;
    for (size_t b=0;b<size();++b)
        if (intersect(box(b),bbox))
            M=math::max(M, basisId(b));
    return M;
}

// Modifiers
void gsBoxList::append(const gsMatrix<real_t> &box, basisIdT basisId )
{
    GISMO_ASSERT(box.rows()==m_dim, "dimension mismatch");
    m_basisId.push_back(basisId);
    m_boxes.resize(m_dim*2+m_boxes.size());
    gsAsMatrix<>(&m_boxes[(size()-1)*2*m_dim],m_dim,2)=box;
}


void gsBoxList::append(const gsBoxList &boxes )
{
    GISMO_ASSERT( m_dim==boxes.m_dim, "dimension mismatch");
    GISMO_ASSERT( !m_haveBasisId || boxes.m_haveBasisId, "dimension mismatch");
    if (m_haveBasisId)
        m_basisId.insert(m_basisId.end(), boxes.m_basisId.begin(), boxes.m_basisId.end());
    m_boxes.insert(m_boxes.end(), boxes.m_boxes.begin(), boxes.m_boxes.end());
}

void gsBoxList::remove (size_t index)
{
    GISMO_ASSERT(index<size() && index>=0,"index out of container");
    if (index!=size()-1)
    {
        box(index)=box(size()-1);
        if (m_haveBasisId)
            basisId(index)=basisId(size()-1);
    }
    if (m_haveBasisId)
        m_basisId.resize(size()-1);
    m_boxes.resize(2*m_dim*m_basisId.size());
}

gsBoxList & gsBoxList::operator-=(const gsMatrix<real_t> &other)
{
    const size_t oldSize=size();
    gsBoxList temp(m_dim, m_haveBasisId);
    for (size_t b=0; b<oldSize;++b)
        temp.append(difference(box(b),other,basisId(b)));
    std::swap(m_basisId,temp.m_basisId);
    std::swap(m_boxes,temp.m_boxes);
    return *this;
}

gsBoxList & gsBoxList::operator-=(const gsBoxList &other)
{
    gsBoxList result(m_dim,m_haveBasisId);
    for (size_t b1=0;b1<size();++b1)
    {
        gsBoxList temp(m_dim, m_haveBasisId);
        temp.append(box(b1),basisId(b1));
        for (size_t b2=0;b2<other.size() && temp.size(); ++b2)
            temp-=other.box(b2);
        result.append(temp);
    }
    std::swap(m_basisId,result.m_basisId);
    std::swap(m_boxes,result.m_boxes);
    return *this;
}

gsBoxList gsBoxList::difference( const gsMatrix<real_t> &box, const gsMatrix<real_t> &toRemove, basisIdT id)
{
    gsBoxList result(m_dim, m_haveBasisId);
    gsMatrix<real_t> newBox=box;
    for(int i=0; i<m_dim;++i)
    {
        const real_t bottom=box(i,0);
        const real_t interface1=math::min(toRemove(i,0),box(i,1));
        const real_t interface2=math::min(toRemove(i,1),box(i,1));
        const real_t top=box(i,1);

        if ( bottom<interface1 )
        {
            newBox(i,1)=interface1;
            result.append(newBox,id);
            newBox(i,1)=top;
        }
        if ( interface2<top)
        {
            newBox(i,0)=interface2;
            result.append(newBox,id);
            newBox(i,0)=bottom;
        }
        if (interface1<interface2)
        {
            newBox(i,0)=interface1;
            newBox(i,1)=interface2;
        }
        else
            break;
    }
    return result;
}

template<int DIM>
std::vector<real_t> getKnotsImpl(const gsFunctionSet<real_t> *basis, gsBoxList::directionT dir)
{
    typedef typename gsBSplineTraits<DIM,real_t>::Basis  KVBasisT;
    const KVBasisT* kv= dynamic_cast<const KVBasisT*>(basis);
    if (kv!=NULL) return kv->knots(dir).unique();
    GISMO_ERROR("unknown basis type");
}

std::vector<real_t> getKnots(const gsFunctionSet<real_t> *basis,
                             gsBoxList::directionT dir, gsBoxList::directionT domainDim)
{
    switch (domainDim)
    {
    case 1:
        return getKnotsImpl<1>(basis,dir);
    case 2:
        return getKnotsImpl<2>(basis,dir);
    case 3:
        return getKnotsImpl<3>(basis,dir);
    case 4:
        return getKnotsImpl<4>(basis,dir);
    default:
        return std::vector<real_t>();
    }
}

gsBoxList::gsBoxList( const std::vector<basisPtr> &bases,
                      directionT domainDim,
                      const std::vector<index_t> &boxes )
    : m_dim(domainDim), m_haveBasisId(true)
{
    const size_t blockSize=2 * domainDim + 1;
    GISMO_ASSERT( boxes.size() % blockSize == 0, "We need one or more full boxes." );

    std::vector<index_t>::const_iterator curBox=boxes.begin();
    std::vector<index_t>::const_iterator endBox=boxes.end();

    gsMatrix<real_t> box;
    index_t          lvl;

    for( ; curBox!=endBox; curBox += blockSize )
    {
        boxFromRefineElementFormat( bases,domainDim, curBox, box, lvl );
        append( box, static_cast<basisIdT>(lvl) );
    }
}

void gsBoxList::toRefineElementFormat( const std::vector<basisPtr> &bases,
                                       std::vector<index_t> &result
                                       ) const
{
    result.clear();
    for( size_t p = 0; p < size(); ++p )
        toRefineElementFormat( bases[basisId(p)].get(), basisId(p),box(p), result );
}


void gsBoxList::toRefineElementFormat(const gsFunctionSet<real_t> *basis,
                               basisIdT            lvl,
                               gsAsConstMatrix<real_t>  box,
                               std::vector<index_t> &boxes)
{
    // We shall be appending to boxes, thus the origSize.
    directionT domainDim = static_cast<directionT>(box.rows());
    size_t origSize = boxes.size();
    boxes.resize(origSize + 2*domainDim + 1);
    boxes[origSize] = lvl;

    const size_t begOffset = origSize + 1;
    const size_t endOffset = begOffset + domainDim;
    for( directionT dir = 0; dir < domainDim; ++dir ) // For each direction.
    {
        std::vector<real_t> uKnots=getKnots(basis,dir,domainDim);
        boxes[begOffset+dir] = indexOfLastLessOrEqual(    uKnots, box(dir,0) );
        boxes[endOffset+dir] = indexOfFirstGreaterOrEqual(uKnots, box(dir,1) );
    }
}


template <class T>
class gsKnotData
{
    typedef const gsKnotVector<T> & KVref;
    typedef const gsKnotVector<T> * KVptr;
public:

    void reserve(size_t n) { m_kv.reserve(n); }
    
    index_t size() { return static_cast<index_t>(m_kv.size()); }

    KVref at(const index_t i)
    { return *m_kv[i]; }

    void push_back(KVref kv) { m_kv.push_back(&kv); }
    
private:
    std::vector<KVptr> m_kv;

private:
    gsKnotData(const gsKnotData &);
};


void gsBoxList::boxFromRefineElementFormat(const std::vector<basisPtr> &bases,
                                            directionT domainDim,
                                            std::vector<index_t>::const_iterator inputBegin,
                                            gsMatrix<real_t> &result,
                                            index_t &level ) const
{
    GISMO_ASSERT(domainDim == m_dim, "The number of variables of hTensorBasis needs to be the same as the dimension of the boxes.");
    result.resize(domainDim, 2);
    level = *inputBegin;

    for( directionT dir = 0; dir < domainDim; ++dir )
    {
        std::vector<real_t> uKnots=getKnots(bases[level].get(),dir,domainDim);
        result( dir, 0 ) = uKnots[ *(inputBegin + dir + 1) ];
        result( dir, 1 ) = uKnots[ *(inputBegin + dir + domainDim + 1) ];
    }
}



std::vector<gsBoxList::basisIdT> gsBoxList::getBasisAt(const gsMatrix<real_t> point ) const
{
    std::vector<basisIdT> result;
    for (size_t b=0; b<size();++b)
    {
        if (contains(box(b),point))
            result.push_back(basisId(b));
    }
    return result;
}

bool gsBoxList::check() const
{
    GISMO_ASSERT(m_boxes.size()%2*m_dim==0,"corrupted boxes list");
    GISMO_ASSERT( (!m_haveBasisId) || (2*m_dim*size()==m_boxes.size()),"corrupted boxes list");
    return true;
}

} // namespace gismo
