/** @file test_interfaces.h

    @brief utils to test with all possible interfaces bewtwee patches

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/



#pragma once

#include  <gsUtils/gsCombinatorics.h>
#include  <gsCore/gsBoundary.h>


namespace gismo {




void setOrientation (int s, gsVector<bool> &o)
{
    for (index_t i=0; i<o.rows();++i)
    {
        o(i)= !(s& 1<<i);
    }
}

int getOrientationIndex (gsVector<bool> &o)
{
    int s=0;
    for (index_t i=0; i<o.rows();++i)
    {
        s|= o(i) ?  0 : 1<<i;
    }
    return s;
}

void firstOrientation (gsVector<bool> &o, boxSide s1, boxSide s2)
{
    o.setConstant(true);
    o( s1.direction() ) = !(s1.parameter() == s2.parameter());
}

bool nextOrientation (gsVector<bool> &o, boxSide s1, boxSide s2)
{
    const int dim=o.rows();
    o(s1.direction() ) = false;
    int s=getOrientationIndex(o);
    ++s;
    setOrientation(s,o);
    o(s1.direction() ) = !(s1.parameter() == s2.parameter());
    return s < 1<<dim;
}

void firstCompatiblePermutation(gsVector<index_t> &per,boxSide s1, boxSide s2)
{
    firstPermutation(per);
    do {
        if ( per(s1.direction()) == s2.direction() )
        {
            break;
        }
    } while (nextPermutation(per));
}

void firstInterface (boundaryInterface &current, int dim)
{
    patchSide p1(0, boxSide::getFirst(dim));
    patchSide p2(1, boxSide::getFirst(dim));
    gsVector<index_t> permutation(dim);
    firstCompatiblePermutation(permutation,p1,p2);
    gsVector<bool>    orientation(dim);
    firstOrientation(orientation,p1,p2);
    current = boundaryInterface(p1,p2,permutation,orientation);
}


bool nextInterface (boundaryInterface &current, int dim)
{
    patchSide         p1          = current.first();
    patchSide         p2          = current.second();
    gsVector<index_t> permutation = current.dirMap(current.first());
    gsVector<bool>    orientation = current.dirOrientation(current.first());

    while ( nextOrientation(orientation,p1,p2) )
    {
        if ( orientation(p1.direction()) == !(p1.parameter()==p2.parameter()) )
            goto construct;
    }
    while (nextPermutation(permutation))
    {
        if ( permutation(p1.direction()) == p2.direction() )
        {
            firstOrientation(orientation,p1,p2);
            goto construct;
        }
    }
    if ( ++p2 < p2.getEnd(dim) )
    {
        firstCompatiblePermutation(permutation,p1,p2);
        firstOrientation(orientation,p1,p2);
        goto construct;
    }
    if ( ++p1 < p1.getEnd(dim) )
    {
        p2=patchSide(1,boxSide::getFirst(dim));
        firstCompatiblePermutation(permutation,p1,p2);
        firstOrientation(orientation,p1,p2);
        goto construct;
    }
    return false;
construct:
    current= boundaryInterface(p1,p2,permutation,orientation);
    return true;
}



gsMatrix<> getCorners(const gsMatrix<> &box)
{
    const index_t dim=box.rows();
    gsMatrix<> result(dim,1<<dim);
    for (boxCorner c=boxCorner::getFirst(dim); c<boxCorner::getEnd(dim);++c)
    {
        const gsVector<bool> param=c.parameters(dim);
        for(int i=0;i<dim;++i)
            result(i,c-1) = param(i) ? box(i,1):box(i,0);
    }
    return result;
}

}
