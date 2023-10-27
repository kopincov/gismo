/** @file gsRemappedTHB.cpp

    @brief Tests gsRemappedTHB, which is a part of the ambitious gsRemappedBasis project.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
**/

#include "gismo_unittest.h"
#include <gsRemappedBasis/gsTensorBasesUtils.h>
#include <gsRemappedBasis/gsDomainMapIterator.h>
#include <gsRemappedBasis/gsRemappedBasis.h>

typedef gsTensorBSplineBasis< 2, real_t> Tbasis;

gsBoxList::basisPtr Converter(gsBasis<>* a)
{
    return gsBoxList::basisPtr(new Tbasis(*dynamic_cast<Tbasis*>(a)));
}

template <int DIM> 
gsBoxList::basisPtr Converter(gsBasis<real_t>* a) 
{ 
    typedef typename gsRemTypes<DIM>::tensorBasisT tensorBasisT; 
    tensorBasisT *comp_basis=dynamic_cast<tensorBasisT*>(a); 
    if (comp_basis) 
        return gsBoxList::basisPtr(comp_basis->clone()); 
    GISMO_ERROR("The provided basis has an unknown type of knot vector or it is not a gsTensorBSplineBasis of the proper dimension."); 
} 

bool lexicoLess( gsVector<> a, gsVector<> b )
{
    for( index_t i = 0; i < a.rows(); ++i )
    {
        if( a(i) < b(i) )
            return true;
        else if ( a(i) > b(i) )
            return false;
    }
    return false;
}

SUITE(gsDomainMapIterator)
{
    TEST(constructor)
    {
        gsKnotVector<real_t> kv(0, 1, 4, 4, 1, 3);
        gsTensorBSplineBasis<2,real_t> tensorBasis(kv,kv);
        gsTHBSplineBasis<2> thbBasis(tensorBasis);

        gsMatrix<real_t> refinementBoxes(2,4);
        refinementBoxes << kv.uValue(0), kv.uValue(1), kv.uValue(4), kv.uValue(5), kv.uValue(2), kv.uValue(3), kv.uValue(1), kv.uValue(2);
        thbBasis.refine(refinementBoxes);

        gsMatrix<real_t> boundingBox(2,2);
        boundingBox<<0,1,0,1;

        gsBoxList bb(2);
        bb.append(refinementBoxes.block(0,0,2,2),1);
        bb.append(refinementBoxes.block(0,2,2,2),1);

        gsDomainMap domainMap;
        domainMap.initFromBoxesMax(bb,boundingBox);

        std::vector< gsBoxList::basisPtr > bPtrs;
        std::vector< Tbasis*> tbases=thbBasis.getBases();
        std::transform(tbases.begin(),tbases.end(),std::back_inserter(bPtrs),Converter<2>);

        gsDomainMapIterator<2> rmpIt( domainMap, bPtrs );
        gsHDomainIterator<real_t, 2> thbIt( thbBasis );

        std::vector< gsVector<real_t> > rmpAnswers;
        std::vector< gsVector<real_t> > thbAnswers;

        for( rmpIt.first(); rmpIt.good(); rmpIt.next() )
            rmpAnswers.push_back( rmpIt.upperCorner() );

        for( ; thbIt.good(); thbIt.next() )
            thbAnswers.push_back( thbIt.upperCorner() );

        std::sort( rmpAnswers.begin(), rmpAnswers.end(), lexicoLess );
        std::sort( thbAnswers.begin(), thbAnswers.end(), lexicoLess );


        CHECK( rmpAnswers == thbAnswers );
    }
}
