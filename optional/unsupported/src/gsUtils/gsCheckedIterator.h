/** @file gsCheckedIterator.h

    @brief A pointer wrapper with range checking for easy debugging.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo {

namespace internal {

//#define make_checked_iterator(arr,n)
//#define make_unchecked_iterator(arr,n)

#ifdef _MSC_VER

#if  _SECURE_SCL != 0 // msvc debug mode ?
template <typename U> class gsCheckedIterator
    : public stdext::checked_array_iterator<U*>
{
public:
    typedef stdext::checked_array_iterator<U*> type;
    inline gsCheckedIterator(U* _arr, size_t _n)
    : type(_arr, _n) { }
};
#else // not msvc debug mode ?
template <typename U> class gsCheckedIterator
    : public stdext::unchecked_array_iterator<U*>
{
public:
    typedef stdext::unchecked_array_iterator<U*> type;
    inline gsCheckedIterator(U* _arr) : type(_arr) { }
};
#endif

#elif !defined(_GLIBCXX_DEBUG) // not gcc debug mode ?

template <typename U> class gsCheckedIterator
{
public:
    typedef U* type;
    inline gsCheckedIterator(U* _arr):arr(_arr){ }
    inline operator type() return arr;
private:
    type arr;
};

#else // debug mode

/**
   \brief A bi-directional knot iterator which provides extended
   information for the iterated knot

   This should be used in extra debug mode only, to avoid inefficient code
   Therefore, use it only protected by debug flags
*/
template <typename U>
class gsCheckedIterator
{
public:
    typedef std::random_access_iterator_tag iterator_category; 
    typedef std::ptrdiff_t                  difference_type;
    typedef U   value_type;
    typedef U&  reference;
    typedef U*  pointer;
    typedef U*  type;

private:
    pointer    m_raw ; ///< pointer which is being checked
    pointer    m_beg ; ///< first position of the pointer
    pointer    m_end ; ///< past-the-end position 
public:

    gsCheckedIterator()
    : m_raw(NULL), m_beg(NULL), m_end(NULL)
    { }

    /**
       \brief Constructs an iterator for the knot-vector \a KV.

       Optionally the iteration starts from from the knot with unique
       index (i.e. without repetitions) equal to \a upos
     */
    gsCheckedIterator(pointer ptr, size_t n)
    : m_raw(ptr), m_beg(ptr), m_end(ptr+n)
    { }

public:

    /// \brief Dereferences the pointer
    reference operator*  () const 
    {
        GISMO_ENSURE(m_raw >= m_beg && m_raw <= m_end, 
                     "Tried to access invalid memory position "<<m_raw<<".");
        return *m_raw;
    }

    pointer   operator-> () const {return m_raw;}

    gsCheckedIterator& operator++() 
    {
        GISMO_ENSURE(++m_raw <= m_end, "Invalid increment.");
        return *this;
    }

    gsCheckedIterator& operator--() 
    {
        GISMO_ENSURE(--m_raw >= 0, "Invalid decrement");
        return *this;
    }

    gsCheckedIterator operator++(int) { gsCheckedIterator tmp(*this); ++(*this); return tmp; }
    gsCheckedIterator operator--(int) { gsCheckedIterator tmp(*this); --(*this); return tmp; }
    bool operator == (const gsCheckedIterator& other) const { return m_raw == other.m_raw;}
    bool operator != (const gsCheckedIterator& other) const {return m_raw != other.m_raw;}
    bool operator <  (const gsCheckedIterator& other) const {return m_raw <  other.m_raw;}
    bool operator >  (const gsCheckedIterator& other) const {return m_raw >  other.m_raw;}
    bool operator <= (const gsCheckedIterator& other) const {return m_raw <= other.m_raw;}
    bool operator >= (const gsCheckedIterator& other) const {return m_raw >= other.m_raw;}

    reference operator [] (ptrdiff_t a) const
    {
        GISMO_ENSURE(m_raw+a>=m_beg && m_raw+a < m_end, 
                     "Invalid memory access "<< m_raw+a);
        return m_raw[a];
    }

    gsCheckedIterator& operator+=(const difference_type & a)
    {
        // Intentionally, we allow jumping to possibly invalid
        // positions in memory.  However, an error will be thrown in
        // the event of calling other operations at invalid state.
        m_raw += a;
        return *this;
    }

    gsCheckedIterator operator+(const difference_type & a) const
    {
        gsCheckedIterator tmp(*this);
        return tmp+=a;
    }

    gsCheckedIterator& operator-=(const difference_type & a) {return operator+=(-a);}

    gsCheckedIterator operator-(const difference_type & a) const
    {
        gsCheckedIterator tmp(*this);
        return tmp-=a;
    }

    friend difference_type operator-(const gsCheckedIterator & l, const gsCheckedIterator & r)
    {return l.m_raw - r.m_raw; }

};

#endif

} // namespace internal



}// namespace gismo

