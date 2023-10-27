/** @file gsWeightMapper.cpp

    @brief test gsWeightMapper class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, F. Buchegger
**/

#include "gismo_unittest.h"
#include <gsMSplines/gsWeightMapper.h>

SUITE(gsWeightMapper)
{

    struct TestMapper
    {
        gsWeightMapper<real_t> map;
        gsSparseMatrix<real_t,Eigen::RowMajor,index_t> matrix;

        TestMapper()
        {
            matrix.resize(10,10);
            for(unsigned i = 0;i<10;++i)
                matrix.at(i,i)=1;
            matrix.at(1,1)=2;
            matrix.at(1,0)=0.5;
            map=matrix;
        }
    };

    //////////////////////////////////////////////////
    // Test: operators and constructors
    //////////////////////////////////////////////////

    TEST_FIXTURE(TestMapper,product)
    {

        gsSparseMatrix<real_t,Eigen::RowMajor, index_t> check1;
        gsWeightMapper<real_t> check2;
        check2=map*matrix;
        check1=matrix*matrix;

        CHECK(gsAllCloseAbsolute(check1.toDense(),check2.asMatrix().toDense(),std::numeric_limits<real_t>::epsilon()));
    }

    //////////////////////////////////////////////////
    // Test: functions for working with the mapper
    //////////////////////////////////////////////////

    TEST_FIXTURE(TestMapper,checksize)
    {
        // checks for getNrOfSources and getNrOfTargets
        CHECK(map.getNrOfSources()==10);
        CHECK(map.getNrOfTargets()==10);
    }

    TEST_FIXTURE(TestMapper,setEntry)
    {
        map.setEntry(0,0,2);
        CHECK(map.asMatrix().at(0,0)==2);
    }

    TEST_FIXTURE(TestMapper,getWeight)
    {
        CHECK(map.getWeight(1,1)==2);
    }

    TEST_FIXTURE(TestMapper,isId)
    {
        // checks for sourceIsId and targetIsId
        CHECK(!map.targetIsId(0));
        CHECK(!map.targetIsId(1));
        CHECK(map.targetIsId(2));
        CHECK(map.sourceIsId(0));
        CHECK(!map.sourceIsId(1));
        CHECK(map.sourceIsId(2));
    }

    //////////////////////////////////////////////////
    // Test: functions for transforming the coefficients
    //////////////////////////////////////////////////

    TEST_FIXTURE(TestMapper,mapToSourceCoefs)
    {
        gsMatrix<real_t> targetCoefs(matrix.rows(),1),sourceCoefs,checkCoefs(matrix.rows(),1);
        for(int i =0;i<targetCoefs.rows();++i)
        {
            targetCoefs(i,0)=1;
            checkCoefs(i,0)= (i==1) ? 2.5 : 1;
        }
        map.mapToSourceCoefs(targetCoefs,sourceCoefs);
        CHECK(gsAllCloseAbsolute(sourceCoefs,checkCoefs,std::numeric_limits<real_t>::epsilon()));
    }

    TEST_FIXTURE(TestMapper,mapToTargetCoefs)
    {
        gsMatrix<real_t> targetCoefs,sourceCoefs(matrix.cols(),1),checkCoefs(matrix.rows(),1);
        for(int i =0;i<sourceCoefs.rows();++i)
        {
            sourceCoefs(i,0)=1;
            checkCoefs(i,0)= (i==0) ? 1 :
                                      (i==1) ? 0.25 : 1;
        }
        map.mapToTargetCoefs(sourceCoefs,targetCoefs);
        CHECK(gsAllCloseAbsolute(targetCoefs,checkCoefs,100*std::numeric_limits<real_t>::epsilon()));
    }


    //////////////////////////////////////////////////
    // Test: functions for applying the map between target and source
    //////////////////////////////////////////////////

    TEST_FIXTURE(TestMapper,sourceToTarget)
    {
        std::vector<real_t> weights;
        std::vector<index_t> target;

        map.sourceToTarget(1,target,weights);
        CHECK(target.size()==2 && weights.size()==2);
        CHECK(target[0]==0);
        CHECK(weights[0]==0.5);
        CHECK(target[1]==1);
        CHECK(weights[1]==2);

        map.sourceToTarget(3,target,weights);
        CHECK(target.size()==1 && weights.size()==1);
        CHECK(target[0]==3);
        CHECK(weights[0]==1);

        std::vector<index_t> sources;
        sources.push_back(1);
        sources.push_back(3);
        map.sourceToTarget(sources,target);
        CHECK(target.size()==3);
        CHECK(target[0]==0);
        CHECK(target[1]==1);
        CHECK(target[2]==3);
    }

    TEST_FIXTURE(TestMapper,targetToSource)
    {
        std::vector<real_t> weights;
        std::vector<index_t> source;

        map.targetToSource(0,source,weights);
        CHECK(source.size()==2 && weights.size()==2);
        CHECK(source[0]==0);
        CHECK(weights[0]==1);
        CHECK(source[1]==1);
        CHECK(weights[1]==0.5);

        map.targetToSource(3,source,weights);
        CHECK(source.size()==1 && weights.size()==1);
        CHECK(source[0]==3);
        CHECK(weights[0]==1);

        std::vector<index_t> targets;
        targets.push_back(0);
        targets.push_back(3);
        map.targetToSource(targets,source);
        CHECK(source.size()==3);
        CHECK(source[0]==0);
        CHECK(source[1]==1);
        CHECK(source[2]==3);
    }

    //////////////////////////////////////////////////
    // Test: functions for fast access to the mapping data
    //////////////////////////////////////////////////

    TEST_FIXTURE(TestMapper,fast_sourceToTarget)
    {
        map.optimize(gsWeightMapper<real_t>::optSourceToTarget);

        std::vector<real_t> weights;
        std::vector<index_t> target;

        gsWeightMapper<real_t>::Iterator it;
        unsigned pos=0;

        for (index_t i=0; i< map.getNrOfSources(); ++i)
        {
            map.sourceToTarget(i,target,weights);
            pos=0;
            for (it=map.fastSourceToTarget(i); it; ++it)
            {
                CHECK(pos<target.size());
                CHECK(target[pos]==it.index());
                CHECK(weights[pos]==it.weight());
                ++pos;
            }
            CHECK(pos==target.size());
        }
    }

    TEST_FIXTURE(TestMapper,fast_targetToSource)
    {
        map.optimize(gsWeightMapper<real_t>::optTargetToSource);

        std::vector<real_t> weights;
        std::vector<index_t> source;

        gsWeightMapper<real_t>::Iterator it;
        unsigned pos=0;

        for (index_t i=0; i< map.getNrOfTargets(); ++i)
        {
            map.targetToSource(i,source,weights);
            pos=0;
            for (it=map.fastTargetToSource(i); it; ++it)
            {
                CHECK(pos<source.size());
                CHECK(source[pos]==it.index());
                CHECK(weights[pos]==it.weight());
                ++pos;
            }
            CHECK(pos==source.size());
        }
    }


}
