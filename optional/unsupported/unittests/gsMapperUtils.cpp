/** @file gsMapperUtils.cpp

    @brief unittest for utils

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include "gismo_unittest.h"
#include <gsUtils/gsCombinatorics.h>
#include <gsMapUtils/gsMapperUtils.h>

SUITE(gsMapperUtils)
{


    struct Identity
    {
        enum { size=6 };
        gsSparseMatrix<real_t> matrix;
        gsWeightMapper<real_t> map;
        gsVector<real_t>       vector;

        Identity()
        {
            matrix.resize(size,size);
            vector.resize(size);
            for (index_t i=0; i< size;++i)
            {
                matrix(i,i)=1;
                vector(i)=i;
            }
            map=matrix;
            map.optimize();
        }

    };


    TEST_FIXTURE(Identity,reorderFullPerm)
    {
        std::vector<index_t> perm(Identity::size);
        gsVector<real_t>     vec2(Identity::size);


        for (size_t i=0; i<perm.size();++i)
            perm[i]=i;
        for (; std::next_permutation(perm.begin(),perm.end());)
        {
            gsWeightMapper<real_t> map2(map);
            for (int i=0; i<size;++i)
                vec2(i)=perm[i];

            reorderMapperTarget(map2,perm);

            CHECK(vector==map2.asMatrix()*vec2);
        }

    }


    TEST_FIXTURE(Identity,reorderPartialPerm)
    {
        std::vector<index_t> perm;
        perm.push_back(3);
        perm.push_back(5);
        perm.push_back(1);

        index_t start=Identity::size-perm.size();

        index_t curPerm=0;
        index_t maxPerm=factorial(perm.size());

        gsVector<real_t>   vec2(Identity::size);
        vec2<<0,2,4,3,5,1;

        for (;curPerm<maxPerm ; std::next_permutation(perm.begin(),perm.end()))
        {
            gsWeightMapper<real_t> map2(map);
            reorderMapperTarget(map2,perm);
            for (int i=start; i<size;++i)
                vec2(i)=perm[i-start];
            ++curPerm;

            CHECK(vector==map2.asMatrix()*vec2);
        }

    }


}
