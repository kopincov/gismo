/** @file gsCloningPtr.h

    @brief Keeps a pointer like unique_ptr, but implement copying by making
    a deep copy.

    This class can be used to store references to polymorphic memebers to
    avoid the need of writing the copy constructors and assignement operators.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsCore/gsTemplateTools.h>

#if __cplusplus >= 201103L
#define NOEXCEPT noexcept
#else
#define NOEXCEPT throw()
#endif

namespace gismo {
namespace memory {

namespace internal {
// This are the possible return values of T::clone()
//  1. T*
template<class T> inline unique_ptr<T> make_unique_if_needed(                 T*  p) { return unique_ptr<T>(p);  }
//  2. unique_ptr<T>
template<class T> inline unique_ptr<T> make_unique_if_needed(const unique_ptr<T>& p) { return p;                 }
} // namespace internal

template<class T> inline unique_ptr<T> clone(const            T & o) { return internal::make_unique_if_needed<T>( o.clone()  ); }
template<class T> inline unique_ptr<T> clone(                 T * p) { return internal::make_unique_if_needed<T>( p->clone() ); }
template<class T> inline unique_ptr<T> clone(const shared_ptr<T>& p) { return internal::make_unique_if_needed<T>( p->clone() ); }
template<class T> inline unique_ptr<T> clone(const unique_ptr<T>& p) { return internal::make_unique_if_needed<T>( p->clone() ); }

/** @class cloning_ptr

    @brief Keeps a pointer like unique_ptr, but implement copying by making
    a deep copy.

    This class can be used to store references to polymorphic memebers to
    avoid the need of writing the copy constructors and assignement operators.

    @ingroup Core
*/

template<class T>
class cloning_ptr
{
private:

    template<class S> friend class cloning_ptr;

    T * m_ptr;

public:

    explicit cloning_ptr( T * p = NULL ) NOEXCEPT : m_ptr( p ) {}

    cloning_ptr(unique_ptr<T> o) NOEXCEPT : m_ptr( o.release() ) {}         //note: no reference in the signature

    template<class S>
    cloning_ptr(unique_ptr<S> o) NOEXCEPT : m_ptr( o.release() ) {}         //note: no reference in the signature

    cloning_ptr(const cloning_ptr & other) NOEXCEPT
    {
        if (other.m_ptr)
            m_ptr = clone(other.m_ptr).release();
        else
            m_ptr = NULL;
    }

    template<class S>
    cloning_ptr(const cloning_ptr<S> & other) NOEXCEPT
    {
        if (other.m_ptr)
            m_ptr = clone(other.m_ptr).release();
        else
            m_ptr = NULL;
    }

#if EIGEN_HAS_RVALUE_REFERENCES
    cloning_ptr(cloning_ptr && other) NOEXCEPT
    {
        m_ptr = other.m_ptr;
        other.m_ptr = NULL;
    }

    template<class S>
    cloning_ptr(cloning_ptr<S> && other) NOEXCEPT
    {
        m_ptr = other.m_ptr;
        other.m_ptr = NULL;
    }
#endif

    ~cloning_ptr() NOEXCEPT
    {
        if (m_ptr) delete m_ptr;
    }

    cloning_ptr& operator=(const cloning_ptr & other) NOEXCEPT
    {
        if (m_ptr) delete m_ptr;
        m_ptr = clone(other.m_ptr).release();
        return *this;
    }

    template<class S>
    cloning_ptr& operator=(const cloning_ptr<S> & other) NOEXCEPT
    {
        if (m_ptr) delete m_ptr;
        m_ptr = clone(other.m_ptr).release();
        return *this;
    }

#if EIGEN_HAS_RVALUE_REFERENCES
    cloning_ptr& operator=(cloning_ptr && other) NOEXCEPT
    {
        if (m_ptr) delete m_ptr;
        m_ptr = other.m_ptr;
        other.m_ptr = NULL;
        return *this;
    }

    template<class S>
    cloning_ptr& operator=(cloning_ptr<S> && other) NOEXCEPT
    {
        if (m_ptr) delete m_ptr;
        m_ptr = other.m_ptr;
        other.m_ptr = NULL;
        return *this;
    }

#endif

    //void reset(T* p = NULL)
    //{
    //    if (m_ptr) delete m_ptr;
    //    m_ptr = p;
    //}

    T & operator*() const
    {
        GISMO_ASSERT( m_ptr, "Access to null pointer." );
        return *m_ptr;
    }

    T * operator->() const
    {
        GISMO_ASSERT( m_ptr, "Access to null pointer." );
        return m_ptr;
    }

    T * get() const 
    {
        return m_ptr;
    }

    T * release()
    {
        T* result = m_ptr;
        m_ptr = NULL;
        return result;
    }

    inline void swap(cloning_ptr & b)
    {
        std::swap(m_ptr,b.m_ptr);
    }

    inline void swap(unique_ptr<T>& b)
    {
        T* tmp = m_ptr;
        m_ptr = b.release();
        b.reset(tmp);
    }

    bool operator!() const
    {
        return m_ptr == NULL;
    }

#if EIGEN_HAS_RVALUE_REFERENCES
    explicit operator bool() const
    {
        return m_ptr != NULL;
    }
#else    
private:

    struct SafeBool
    {
        SafeBool(int) {}
        void dummy() {}
    };

    typedef void (SafeBool::*bool_cast_type)();

public:

    operator bool_cast_type() const
    {
        return !m_ptr ? 0 : &SafeBool::dummy;
    }
#endif

};

template<class T> bool operator==(const cloning_ptr<T> & p1, const cloning_ptr<T> & p2) { return p1.get()==p2.get(); }
template<class T> bool operator!=(const cloning_ptr<T> & p1, const cloning_ptr<T> & p2) { return p1.get()!=p2.get(); }
template<class T> bool operator< (const cloning_ptr<T> & p1, const cloning_ptr<T> & p2) { return p1.get()< p2.get(); }
template<class T> bool operator> (const cloning_ptr<T> & p1, const cloning_ptr<T> & p2) { return p1.get()> p2.get(); }
template<class T> bool operator<=(const cloning_ptr<T> & p1, const cloning_ptr<T> & p2) { return p1.get()<=p2.get(); }
template<class T> bool operator>=(const cloning_ptr<T> & p1, const cloning_ptr<T> & p2) { return p1.get()>=p2.get(); }

template<class T> inline void swap(cloning_ptr<T>&  a, cloning_ptr<T>&  b) { a.swap(b); }
template<class T> inline void swap(unique_ptr<T>& b, cloning_ptr<T>&  a) { a.swap(b); }
template<class T> inline void swap(cloning_ptr<T>&  a, unique_ptr<T>& b) { a.swap(b); }

template<class T> inline unique_ptr<T> clone(const cloning_ptr<T>& p) { return internal::make_unique_if_needed<T>( p->clone() ); }
template<class T> inline unique_ptr<T> give(cloning_ptr<T>& p) { return unique_ptr<T>(p.release()); }


template<class T> inline cloning_ptr<T> make_cloning(T* p) { return cloning_ptr<T>(p); }

template<class T> inline std::vector< cloning_ptr<T> > make_cloning(std::vector<T*>& cont)
{
    std::vector< cloning_ptr<T> > result;
    for (typename std::vector<T*>::iterator it = cont.begin(); it != cont.end(); ++it)
        result.push_back( cloning_ptr<T>(*it) );
    cont.clear();
    return result;
}

template <class T> inline std::vector<T*> get_raw(std::vector< cloning_ptr<T> >& cont)
{
    std::vector<T*> result;
    for (typename std::vector< cloning_ptr<T> >::iterator it = cont.begin(); it != cont.end(); ++it)
        result.push_back( (*it).get() );
    cont.clear();
    return result;
}


} // namespace memory
} // namespace gismo
