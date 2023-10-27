/** @file gsMapperL2G.cpp

    @brief test for gsL2G based on gsWeightMapper

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include "gismo_unittest.h"
#include <gsMapUtils/gsL2GMapper.h>

SUITE(gsMapperL2G)
{
    typedef gsSparseMatrix<real_t>                   SMatT;
    typedef gsMatrix<real_t>                         FMatT;
    // standard writers
    typedef gsMatAndRhsModWriter<SMatT,SMatT>        SysW;
    typedef gsMatAndRhsModWriter<FMatT,gsNullWriter<real_t> > RhsW;
    // multiplier writer
    typedef gsMultiplierWriter<SysW>                 SymW;
    typedef gsShiftWriter<gsMultiplierWriter<SysW> > MulW;
    // writers for the eliminated system
    typedef gsBoundaryWriter<SMatT>                  EliSysW;
    typedef gsBoundaryWriter<FMatT>                  EliRhsW;


    class data
    {
    public:
        real_t                   eps;
        int                      size;
        gsWeightMapper<real_t>   mapper;
        gsSparseMatrix<real_t>   matS;
        gsSparseMatrix<real_t>   matR;
        gsMatrix<real_t>         locM;

        data()
            : eps(1e-10), size(10)
        {
            gsSparseMatrix<real_t,Eigen::RowMajor,index_t> matrix;
            matrix.resize(size,size);
            for(int i = 0;i<size;++i)
                matrix.at(i,i)=1;
            matrix.at(1,1)=2;
            matrix.at(1,0)=0.5;
            mapper=matrix;
            mapper.optimize();

            locM.resize(2,2);
            locM<<1,2,3,4;

            matS.resize(size,size);
            matR.resize(size,size);
        }

        void map (gsL2GMapper<SysW> &l2g, gsVector<unsigned>  rows, gsVector<unsigned> cols)
        {
            matS.setZero();
            matR.setZero();
            l2g.store(rows,cols,locM);
//            std::cout<<"----------------------------------\n";
//            std::cout<<"\n"<<matS.toDense()<<"\n\n"<<matR.toDense()<<"\n"<<std::endl;
        }

    };



    TEST_FIXTURE(data, singleWriting)
    {
        gsL2GMapper<SysW>        l2g(mapper,mapper,SysW(size,size,matS,matR));
        gsVector<unsigned>       actC(2);
        gsVector<unsigned>       actR(2);
        actC<<8,9;
        actR<<8,9;
        map(l2g,actR,actC);
        CHECK( (matS.block(8,8,2,2).toDense()-locM).norm()<eps);

        actC<<4,5;
        actR<<8,9;
        map(l2g,actR,actC);
        CHECK( (matS.block(8,4,2,2).toDense()-locM).norm()<eps);
    }





}
