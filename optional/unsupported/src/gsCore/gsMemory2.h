/** @file gsMemory2.h
    @brief Provides utility function related to memory management.
    This file is part of the G+Smo library. 
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/

#pragma once

#include<gsCore/gsMemory.h>


namespace gismo {
namespace memory {


#if __cplusplus >= 201103L
using std::static_pointer_cast;
using std::dynamic_pointer_cast;
using std::const_pointer_cast;

#if  __cplusplus >= 201703L
using std::reinterpret_pointer_cast;
#else // __cplusplus < 201703L
// as in http://en.cppreference.com/w/cpp/memory/shared_ptr/pointer_cast
template<class T, class S>
shared_ptr<T> reinterpret_pointer_cast( const shared_ptr<S>& s ) noexcept
{
    auto p = reinterpret_cast<typename std::shared_ptr<T>::element_type*>(s.get());
    return shared_ptr<T>( s, p );
}
#endif

#else // __cplusplus < 201103L

// Implemented as in http://en.cppreference.com/w/cpp/memory/shared_ptr/pointer_cast
// but replace aliasing constructor by our own detail::aliasing_deleter

namespace detail {
    template<class T, class S>
    struct aliasing_deleter {
        aliasing_deleter( shared_ptr<S> s ) : m_s(give(s)) {}
        void operator()(T* ptr) {}
    private:
        shared_ptr<S> m_s;
    };
}

template<class T, class S>
shared_ptr<T> static_pointer_cast( const shared_ptr<S>& s )
{
    typedef std::shared_ptr<T>::element_type P;
    P* p = static_cast<P*>(s.get());
    return shared_ptr<T>( p, detail::aliasing_deleter<T,S>(s)  );
}


template<class T, class S>
shared_ptr<T> dynamic_pointer_cast( const shared_ptr<S>& s )
{
    typedef std::shared_ptr<T>::element_type P;
    P* p = dynamic_cast<P*>(s.get());
    if (p)
        return shared_ptr<T>( p, detail::aliasing_deleter<T,S>(s)  );
    else
        return shared_ptr<T>();
}

template<class T, class S>
shared_ptr<T> const_pointer_cast( const shared_ptr<S>& s )
{
    typedef std::shared_ptr<T>::element_type P;
    P* p = const_cast<P*>(s.get());
    return shared_ptr<T>( p, detail::aliasing_deleter<T,S>(s)  );
}

template<class T, class S>
shared_ptr<T> reinterpret_pointer_cast( const shared_ptr<S>& s )
{
    typedef std::shared_ptr<T>::element_type P;
    P* p = reinterpret_cast<P*>(s.get());
    return shared_ptr<T>( p, detail::aliasing_deleter<T,S>(s)  );
}


#endif


} // namespace memory
} // namespace gismo
