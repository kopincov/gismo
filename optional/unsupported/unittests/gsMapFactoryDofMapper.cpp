/** @file gsMapFactoryDofMapper.cpp

    @brief test gsMapFactoryDofMapper class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
**/

#include "gismo_unittest.h"
#include <gsMSplines/gsWeightMapper.h>
#include <gsMapUtils/gsMapFactoryDofMapper.h>

SUITE(gsMapper_Factory)
{
    struct TestDofMapper
    {
        gsDofMapper dofMapper;

        TestDofMapper()
        {
            gsVector<index_t> patchDofSizes(1);
            patchDofSizes(0)=10;
            dofMapper=gsDofMapper(patchDofSizes);
            gsMatrix<index_t> boundaryDofs(2,1);
            boundaryDofs(0,0)=0;
            boundaryDofs(1,0)=9;
            dofMapper.markBoundary(0,boundaryDofs);
            dofMapper.finalize();
        }
    };

    TEST_FIXTURE(TestDofMapper,gsMapFactory_DofMapperKeepBoundary)
    {
        gsMapFactoryDofMapper mapFactory(dofMapper,true);
        gsWeightMapper<real_t>* mapper = mapFactory.makeMapper();
        index_t globalIndex;
        std::vector<index_t> source;

        globalIndex = dofMapper.index(0,0);
        source.clear();
        mapper->sourceToTarget(0,source);
        CHECK(source.size()==1 && source[0]==globalIndex);

        globalIndex = dofMapper.index(9,0);
        source.clear();
        mapper->sourceToTarget(9,source);
        CHECK(source.size()==1 && source[0]==globalIndex);

        globalIndex = dofMapper.index(5,0);
        source.clear();
        mapper->sourceToTarget(5,source);
        CHECK(source.size()==1 && source[0]==globalIndex);

        delete mapper;
    }

    TEST_FIXTURE(TestDofMapper,gsMapFactory_DofMapperForgetBoundary)
    {
        gsMapFactoryDofMapper mapFactory(dofMapper,false);
        gsWeightMapper<real_t>* mapper = mapFactory.makeMapper();
        index_t globalIndex;
        std::vector<index_t> source;

        source.clear();
        mapper->sourceToTarget(0,source);
        CHECK(source.size()==0);

        source.clear();
        mapper->sourceToTarget(9,source);
        CHECK(source.size()==0);

        globalIndex = dofMapper.index(5,0);
        source.clear();
        mapper->sourceToTarget(5,source);
        CHECK(source.size()==1 && source[0]==globalIndex);

        delete mapper;
    }
}

