/** @file gsSimplexIterator.h

    @brief Provides iteration over integer or numeric points in a simplex

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsUtils/gsCompositionIterator.h>

namespace gismo
{

/**
   Iterator over points in a simplex
 */
template<typename T, int d = -1>
class gsSimplexIterator
{
public:
    static const bool isIntegral = std::numeric_limits<T>::is_integer;
    
    typedef typename util::conditional<isIntegral,
                                       gsVector<T,d>,gsMatrix<T,d> >::type point;
    
    typedef typename util::conditional<isIntegral,T,index_t>::type integer_t;
    
    typedef gsCompositionIterator<integer_t,(d==-1?-1:d+1)> CompositionIterator;
    
    typedef typename CompositionIterator::point baryIndex;
    
    // only used for non-integral types
    typedef typename util::conditional<isIntegral,gsMatrix<T,0,0>,
                                       gsMatrix<T,d,(d==-1?-1:d+1)> >::type vertex_matrix;
public:
    
    template<typename U> // enabled for integral types
    explicit gsSimplexIterator(const U side, const int dim,//U=int
                               typename util::enable_if<isIntegral,U>::type * = NULL)
    {
        reset(side, dim);
    }

    template<typename U> // enabled for integral types and static dimension
    explicit gsSimplexIterator(const U side,
                               typename util::enable_if<isIntegral && d!=-1,U>::type * = NULL)
    {
        reset(side, d);
    }

    template<typename U> // enabled for non-integral types
    gsSimplexIterator(const vertex_matrix & vertices, U pointsPerSide,//U=int
                      typename util::enable_if<!isIntegral,U>::type * = NULL)
    : m_vert(vertices)
    {
        reset(pointsPerSide, vertices.rows());
    }
    
    inline void reset(const int side, int dim = -1)
    {
        if (d!=-1)
        {
            GISMO_ASSERT(dim==-1 || dim==d, "gsSimplexIterator: Invalid or unknown dimension");
            dim = d;
        }
        
        m_comp.reset(side, dim+1);
        initCurr(m_cur,dim);
        /*
        if (isIntegral)
        {
            m_cur.derived().resize(dim);
            m_cur.setZero();
        }
        else
            m_cur = m_vert.col(0);
        */
    }

    void reset() { reset(m_comp.sum(), m_cur.rows()); }

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(point)%16)==0 );

public:

    operator bool() const {return m_comp;}
    
    const point & operator*() const {return m_cur;}
    
    const point * operator->() const {return &m_cur;}
    
    inline gsSimplexIterator & operator++()
    {
        if (++m_comp)
        {
            nextCurr(m_cur);
            /*
            if (isIntegral)
                m_cur = m_comp->bottomRows(m_comp->size()-1);
            else
                m_cur.noalias() = m_vert * m_comp->template cast<T>() / m_comp.sum();
            */
        }
        return *this;
    }

    index_t numPoints() const {return m_comp.numPoints();}

    baryIndex barycentricIndex() const {return *m_comp;}

    bool isBoundary() const {return (m_comp->array()==0).any();}
    
    //int sideLength() const {return m_comp.sum();}

private:
    
    inline void initCurr(gsVector<T,d> & cur, int vd) // isIntegral
    { cur.derived().resize(vd); cur.setZero(); }

    inline void initCurr(gsMatrix<T,d> & cur, int vd)
    { cur = m_vert.col(0); }

    inline void nextCurr(gsVector<T,d> & cur) // isIntegral
    { cur = m_comp->bottomRows(m_comp->size()-1); }

    inline void nextCurr(gsMatrix<T,d> & cur)
    { cur.noalias() = m_vert * m_comp->template cast<T>() / m_comp.sum(); }

private:

    CompositionIterator m_comp;
    
    vertex_matrix m_vert;
    
    point m_cur;
};


} // namespace gismo

