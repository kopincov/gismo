/** @file gsTensorBasesUtils

    @brief code that deals with gsTensorBsplineBasis in the gsRemappedBasis folde.
    Mostly conversion utilities

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsCore/gsFunctionSet.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

// common typedef to shorten tensor BSpline names
template <int DIM>
struct gsRemTypes
{
    typedef gsKnotVector<real_t>                        knotVectorT;
    typedef typename gsBSplineTraits<DIM,real_t>::Basis tensorBasisT;
    // TODO C++11 move to sharedPtr
    typedef tensorBasisT*     tensorBasisPtr;
    typedef gsFunctionSet<real_t>::Ptr basisPtr;
    typedef int directionT;
    typedef unsigned basisIdT;
};

// utility function to convert to the tensor basis type
// it is used in tests
template <int DIM>
inline typename gsRemTypes<DIM>::tensorBasisT* asBasis(typename gsRemTypes<DIM>::basisPtr obn)
{
    return static_cast<typename gsRemTypes<DIM>::tensorBasisT*>(obn.get());
}

template <int DIM>
class gsTensorIndexSubDomain;

template <int DIM>
class gsTensorIndexDomain
{
public:
    typedef gsVector<index_t,DIM>    MIndexT;
    typedef gsMatrix<index_t,DIM,2>  DomainT;
    typedef index_t                  FIndexT;
    friend class gsTensorIndexSubDomain<DIM>;
protected:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    FIndexT m_schi; // scalar shift
    FIndexT m_size;
    DomainT m_domm; // domain
    MIndexT m_strd; // strides: (1, dom(0,1)-dom(0,0)),
public:
    gsTensorIndexDomain(){}

    gsTensorIndexDomain(const DomainT &dom)
        : m_domm(dom)
    {
        m_strd(0)=1;
        for(int d=1; d<DIM; ++d)
            m_strd(d)=m_strd(d-1)*(m_domm(d-1,1)-m_domm(d-1,0));
        m_size = m_strd(DIM-1)*(m_domm(DIM-1,1)-m_domm(DIM-1,0));
        m_schi = -m_domm.col(0).dot(m_strd);
    }
    gsTensorIndexDomain(const MIndexT &max)
    {
        m_domm.col(0).setZero();
        m_domm.col(1)=max;
        m_strd(0)=1;
        for(int d=1; d<DIM; ++d)
            m_strd(d)=m_strd(d-1)*(m_domm(d-1,1)-m_domm(d-1,0));
        m_size = m_strd(DIM-1)*(m_domm(DIM-1,1)-m_domm(DIM-1,0));
        m_schi = -m_domm.col(0).dot(m_strd);
    }

    operator DomainT() {return m_domm;}

    bool contains(const MIndexT &point) const
    {
        return (m_domm.col(0).array()<=point.array()).all() && (m_domm.col(1).array()>point.array()).all();
    }
    bool contains(const FIndexT &point) const
    {
        return contains(toMulti(point));
    }

    typename DomainT::constColumn lower() const {return m_domm.col(0);}
    typename DomainT::constColumn upper() const {return m_domm.col(1);}

    FIndexT minFlat() const {return toFlat(lower()); }
    FIndexT endFlat() const {return toFlat(upper().array()-1)+1; }


    index_t size()    const {return m_size;}
    MIndexT sizeM()   const {return m_domm.col(1)-m_domm.col(0);}

    MIndexT toMulti (FIndexT i) const
    {
        i-=m_schi;
        return m_domm.col(0)+toIncr(i);
    }
    FIndexT toFlat  (const MIndexT &i) const
    {
        return i.dot(m_strd)+m_schi;
    }

    // iterator helpers
    MIndexT toIncr (index_t i) const
    {
        MIndexT r;
        for (int d=DIM-1; d>0; --d)
        {
            r(d)=i/m_strd(d);
            i%=m_strd(d);
        }
        r(0)=i;
        return r;
    }
    index_t toPos (const MIndexT &incr) const
    {
        return (incr-lower() ).dot(m_strd);
    }

    index_t toPos (const FIndexT &incr) const
    {
        return toPos(toMulti(incr));
    }

    class gsTensorMIndexIter
    {
    public:
        typedef       MIndexT  value_type;
        typedef const MIndexT& reference;
        typedef const MIndexT* pointer;
        typedef ptrdiff_t difference_type;
        typedef std::random_access_iterator_tag iterator_category; //or another tag
    protected:
        const gsTensorIndexDomain<DIM> *m_dom;
        MIndexT                         m_pos;
    public:

        FIndexT flatIndex() const
        {
            return m_dom->toFlat(m_pos);
        }
        MIndexT multiIndex() const
        {
            return m_pos;
        }
        FIndexT posIndex() const
        {
            return m_dom->toPos();
        }
        // construction
        gsTensorMIndexIter(const gsTensorIndexDomain &dom, const MIndexT &pos)
            : m_dom(&dom),m_pos(pos)
        {}
        gsTensorMIndexIter()
            : m_dom(NULL)
        {}

        // iterator interface
        gsTensorMIndexIter& operator++()
        {
            int i;
            for (i=0; i<DIM && m_pos(i)==m_dom->m_domm(i,1)-1; ++i)
                m_pos(i)=m_dom->m_domm(i,0);
            if(i!=DIM)
                m_pos(i)+=1;
            else
                m_pos(DIM-1)=m_dom->m_domm(DIM-1,1);
            return *this;
        }
        gsTensorMIndexIter& operator--()
        {
            int i;
            for (i=0; i<DIM && m_pos(i)==m_dom->m_domm(i,0)-1; ++i)
                m_pos(i)=m_dom->m_domm(i,1)-1;
            if(i!=DIM)
                m_pos(i)-=1;
            return *this;
        }
        gsTensorMIndexIter operator++(int) {gsTensorMIndexIter other(*this); ++(*this); return other;}
        gsTensorMIndexIter operator--(int) {gsTensorMIndexIter other(*this); --(*this); return other;}

        gsTensorMIndexIter& operator+=(index_t a)             {m_pos+=m_dom->toIncr(a); return (*this);}
        gsTensorMIndexIter  operator+ (index_t a)   const     {gsTensorMIndexIter other(*this); return other+=a;}
        gsTensorMIndexIter& operator-=(index_t a)             {m_pos-=m_dom->toIncr(a); return (*this);}
        gsTensorMIndexIter  operator- (index_t a)   const     {gsTensorMIndexIter other(*this); return other-=a;}

        // comparison
        bool operator< (const gsTensorMIndexIter& other) const
        {
            for (int i=DIM-1;i>=0;--i)
                if (m_pos(i)<other.m_pos(i))
                    return true;
            return false;
        }
        bool operator> (const gsTensorMIndexIter& other) const
        {
            for (int i=DIM-1;i>=0;--i)
                if (m_pos(i)>other.m_pos(i))
                    return true;
            return false;
        }
        bool operator== (const gsTensorMIndexIter& other) const {return m_pos==other.m_pos;}
        bool operator!= (const gsTensorMIndexIter& other) const {return !(m_pos==other.m_pos);}
        bool operator>= (const gsTensorMIndexIter& other) const {return (*this)==other || (*this)>other;}
        bool operator<= (const gsTensorMIndexIter& other) const {return (*this)==other || (*this)<other;}

        // difference
        difference_type   operator-  (const gsTensorMIndexIter& other) const {return m_dom->toPos(m_pos-other.m_pos+m_dom->lower());}

        // access
        reference   operator*  ()           const  {return m_pos;}
        pointer     operator-> ()           const  {return &(m_pos);}
        value_type  operator[] (index_t a)  const  {gsTensorFIndexIter other(*this); other+=a; return *other;}
    };
    typedef gsTensorMIndexIter const_iterator;
    const_iterator begin() const {if (size()>0) return gsTensorMIndexIter(*this,m_domm.col(0)); else return end();}
    const_iterator end()   const {return gsTensorMIndexIter(*this,endPos());}
    MIndexT endPos() const
    {
        MIndexT endPos=m_domm.col(0);
        endPos(DIM-1)=m_domm(DIM-1,1);
        return endPos;
    }


    class gsTensorFIndexIter : public gsTensorMIndexIter
    {
    protected:
        FIndexT  m_flat;
    public:
        typedef       FIndexT  value_type;
        typedef const FIndexT& reference;
        typedef const FIndexT* pointer;
        typedef std::random_access_iterator_tag iterator_category; //or another tag
    public:
        // construction
        gsTensorFIndexIter(const gsTensorIndexDomain &dom, const MIndexT &pos)
            : gsTensorMIndexIter(dom,pos)
        {}
        gsTensorFIndexIter() {}

        value_type  operator*  ()           const  {return gsTensorMIndexIter::flatIndex();}
        pointer     operator-> ()           const  {return m_flat=gsTensorMIndexIter::flatIndex(); return &m_flat;}
        value_type  operator[] (index_t a)  const  {gsTensorMIndexIter other(*this); other+=a; return *other;}
    };
    typedef gsTensorFIndexIter flat_iterator;
    flat_iterator flat_begin() const {if (size()>0) return gsTensorFIndexIter(*this,m_domm.col(0)); else return flat_end();}
    flat_iterator flat_end()   const {return gsTensorFIndexIter(*this,endPos());}
};

template <int DIM>
class gsTensorIndexSubDomain : public gsTensorIndexDomain<DIM>
{
public:
    typedef gsVector<index_t,DIM>    MIndexT;
    typedef gsMatrix<index_t,DIM,2>  DomainT;
    typedef index_t                  FIndexT;

    using gsTensorIndexDomain<DIM>::lower;
    using gsTensorIndexDomain<DIM>::upper;
    using gsTensorIndexDomain<DIM>::size;
    using gsTensorIndexDomain<DIM>::toPos;
    using gsTensorIndexDomain<DIM>::contains;
protected:
    MIndexT m_conv;
    using gsTensorIndexDomain<DIM>::m_schi;
    using gsTensorIndexDomain<DIM>::m_domm;
    using gsTensorIndexDomain<DIM>::endPos;
public:
    gsTensorIndexSubDomain() {}
    gsTensorIndexSubDomain( const gsTensorIndexDomain<DIM> &dom, const DomainT &subDom)
        : gsTensorIndexDomain<DIM>(subDom)
    {
        m_schi = dom.m_schi;
        m_conv = dom.m_strd;
    }


    MIndexT toMulti (FIndexT i) const
    {
        MIndexT r;
        i-=m_schi;
        for (int d=DIM-1; d>=0; --d)
        {
            r(d)=i/m_conv(d);
            i%=m_conv(d);
        }
        return r;
    }
    FIndexT toFlat  (const MIndexT &i) const
    {
        return i.dot(m_conv)+m_schi;
    }

    index_t toPos (const FIndexT &incr) const
    {
        return toPos(toMulti(incr));
    }

    bool contains(const FIndexT &point) const
    {
        return contains(toMulti(point));
    }

    FIndexT minFlat() const {return toFlat(lower()); }
    FIndexT endFlat() const {return toFlat(upper().array()-1)+1; }

    class gsSubTensorMIndexIter : public gsTensorIndexDomain<DIM>::gsTensorMIndexIter
    {
        using gsTensorIndexDomain<DIM>::gsTensorMIndexIter::m_dom;
        using gsTensorIndexDomain<DIM>::gsTensorMIndexIter::m_pos;
    public:


        FIndexT flatIndex() const
        {
            return static_cast<const gsTensorIndexSubDomain*>(m_dom)->toFlat(m_pos);
        }
        FIndexT posIndex() const
        {
            return static_cast<const gsTensorIndexSubDomain*>(m_dom)->toPos(m_pos);
        }
        // construction
        gsSubTensorMIndexIter(const gsTensorIndexSubDomain &dom, const MIndexT &pos)
            : gsTensorIndexDomain<DIM>::gsTensorMIndexIter(dom,pos)
        {}
        gsSubTensorMIndexIter() {}
    };
    typedef gsSubTensorMIndexIter const_iterator;
    const_iterator begin() const {if (size()>0) return gsSubTensorMIndexIter(*this,m_domm.col(0)); else return end();}
    const_iterator end()   const {return gsSubTensorMIndexIter(*this,endPos());}

    class gsSubTensorFIndexIter : public gsSubTensorMIndexIter
    {
    protected:
        FIndexT  m_flat;
    public:
        typedef       FIndexT  value_type;
        typedef const FIndexT& reference;
        typedef const FIndexT* pointer;
        typedef std::random_access_iterator_tag iterator_category; //or another tag
    public:
        // construction
        gsSubTensorFIndexIter(const gsTensorIndexSubDomain &dom, const MIndexT &pos)
            : gsSubTensorMIndexIter(dom,pos)
        {}
        gsSubTensorFIndexIter() {}

        value_type  operator*  ()           const  {return gsSubTensorMIndexIter::flatIndex();}
        pointer     operator-> ()           const  {return m_flat=gsSubTensorMIndexIter::flatIndex(); return &m_flat;}
        value_type  operator[] (index_t a)  const  {gsSubTensorFIndexIter other(*this); other+=a; return *other;}
    };
    typedef gsSubTensorFIndexIter flat_iterator;
    flat_iterator flat_begin() const {if (size()>0) return gsSubTensorFIndexIter(*this,m_domm.col(0)); else return flat_end();}
    flat_iterator flat_end()   const {return gsSubTensorFIndexIter(*this,endPos());}
};

// TODO move this function to be a member of gsTensorBSplineBasis
/**
  Adds the indices of the basis functions of \a tensorBasis that are active on \a box to \a set
**/
template <int DIM>
static void addInfluencingForBox(const gsMatrix<real_t> &box, const typename gsRemTypes<DIM>::tensorBasisT &tensorBasis, std::set<index_t> &set)
{
    typedef typename gsRemTypes<DIM>::knotVectorT knotVectorT;
    gsVector<index_t, DIM> maxIter;
    gsMatrix<index_t, DIM,2> limit; // As in the makeSubGridIterator below.
    for( int i = 0; i < DIM; ++i )
    {
        const knotVectorT &knots = tensorBasis.knots(i);
        maxIter[i] = knots.size() - tensorBasis.degree(i)-1; // knots are numbered from 0
        limit.row(i)<< (std::upper_bound(knots.begin(),knots.end(),box(i,0))-knots.begin()- tensorBasis.degree(i)-1),
                       (std::lower_bound(knots.begin(),knots.end(),box(i,1))-knots.begin());
    }
    gsTensorIndexSubDomain<DIM> other(gsTensorIndexDomain<DIM>(maxIter),limit);
    typename gsTensorIndexSubDomain<DIM>::const_iterator it=other.begin(), end=other.end();
    while(it!=end)
    {
        set.insert(it.flatIndex());
        ++it;
    }
}


} // namespace gismo

