/** @file gsCompositionIterator.h

    @brief Provides iteration over integer compositions

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once


namespace gismo
{

/**
   Iterator over integer compositions
 */
template<class Z, int d = -1>
class gsCompositionIterator
{
public:
    typedef gsVector<Z,d> point;

public:

    gsCompositionIterator() : m_sum(-1)
    {
        GISMO_STATIC_ASSERT( std::numeric_limits<Z>::is_integer,
            "The template parameter needs to be an integer type." );
    }
    
    explicit gsCompositionIterator(const int sum, const int dim = -1)
    { reset(sum,dim); }

    inline void reset(const int sum, int dim = -1)
    {
        if (d!=-1)
        {
            GISMO_ASSERT(dim==-1 || dim==d, "gsCompositionIterator: Invalid or unknown dimension");
            dim = d;
        }
        GISMO_ASSERT(dim>0, "gsCompositionIterator: Invalid dimension");

        m_sum = sum;
        m_cur.derived().resize(dim);
        m_cur.setZero();
        m_cur[0] = sum;
    }

    void reset() { reset(m_cur.sum(), m_cur.rows()); }

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(point)%16)==0 );

public:

    operator bool() const {return -1 != m_sum;}
    
    const point & operator*() const {return m_cur;}
    
    const point * operator->() const {return &m_cur;}
    
    inline gsCompositionIterator & operator++()
    {
        if (m_cur[m_cur.size()-1] != static_cast<Z>(m_sum))
        {
            for (index_t i = 0; i != m_cur.size(); ++i)
            {
                if ( 0 != m_cur[i] )
                {
                    const Z t = m_cur[i];
                    m_cur[i]    = 0  ;
                    m_cur[0]    = t-1;
                    m_cur[i+1] += 1  ;
                    return *this;
                }
            }
        }
        m_sum = -1; //done
        return *this;
    }

    index_t numPoints() const
    { return binomial(m_cur.sum()+m_cur.rows()-1,m_cur.rows()-1); }

    int sum() const {return m_sum;}
    
private:
    
    point m_cur;
    
    int m_sum;
};

//firstComposition
//nextComposition
//numCompositions

} // namespace gismo

