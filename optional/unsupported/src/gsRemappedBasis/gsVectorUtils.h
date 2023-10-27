/** @file gsKnotVectorUtils.h

    @brief Utility functions for knot vectors.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
**/

#pragma once

#include <algorithm>
#include <vector>
#include <iterator> // needed for std::back_inserter on mvc compilers

namespace gismo
{

template<typename T, typename  KVType>
size_t indexOfLastLessOrEqual( const KVType &knots, const T value )
{
    typename KVType::const_iterator  it;
    it=std::upper_bound(knots.begin(),knots.end(),value);
    if (it==knots.begin()) GISMO_ERROR("point outside knot vector");
    return it-knots.begin()-1;
}

template<typename T, typename KVType>
size_t indexOfFirstGreaterOrEqual( const KVType &knots, const T value )
{
    typename KVType::const_iterator  it;
    it=std::lower_bound(knots.begin(),knots.end(),value);
    if (it==knots.end()) GISMO_ERROR("point outside knot vector");
    return it-knots.begin();
}

template<typename T, template<typename > class KVType1, template<typename > class KVType2>
std::vector<T> vectorDiff( const KVType1<T> &knots1, const KVType2<T> &knots2 )
{
    std::vector<T> result;
    std::set_difference(knots1.begin(), knots1.end(), knots2.begin(), knots2.end(),std::back_inserter(result));
    return result;
}

// Since std::vector has two templates parameters, we need a separate version.
template<typename T, typename Alloc, template<typename,typename > class KVType >
KVType<T,Alloc> vectorDiff( const KVType<T,Alloc> &knots1, const KVType<T,Alloc> &knots2 )
{
    KVType<T,Alloc> result;
    std::set_difference(knots1.begin(), knots1.end(), knots2.begin(), knots2.end(),std::back_inserter(result));
    return result;
}

template<typename  KVType1, typename  KVType2, typename  KVType3>
void vectorUnion( const KVType1 &knots1, const KVType2 &knots2, KVType3 &result )
{
    std::set_union(knots1.begin(), knots1.end(), knots2.begin(), knots2.end(),std::back_inserter(result));
}



template<typename  KVType>
void vectorUnionMany (typename std::vector<KVType>::const_iterator inputBeg, typename std::vector<KVType>::const_iterator inputEnd, KVType &output)
{
    size_t numEl=inputEnd-inputBeg;
    switch (numEl)
    {
    case 0:
        break;
    case 1:
        output=*inputBeg;
        break;
    case 2:
        vectorUnion(*inputBeg,*(inputBeg+1),output);
        break;
    default:
        KVType temp1, temp2;
        vectorUnionMany(inputBeg,inputBeg+numEl/2+1,temp1);
        vectorUnionMany(inputBeg+numEl/2+1,inputEnd ,temp2);
        vectorUnion(temp1,temp2,output);
    }
}


template<typename T,typename KVType>
size_t addUnique(T val, KVType &vec)
{
    typename KVType::iterator it = std::lower_bound(vec.begin(),vec.end(),val);
    if (it !=vec.end() && *it==val)
        return it-vec.begin();
    vec.insert(it,val);
    return vec.size();
}


} // namespace gismo
