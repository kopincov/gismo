/** @file gsExportMatrix_test.cpp

    @brief Tests export and import of matrices between Eigen formats and ASCII files suitable for Matlab.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Authors: A. Bressan, D. Mokris
**/

#include "gismo_unittest.h"
#include <stdlib.h>
#include <time.h>
#include <gsUtils/gsExportMatrix.h>


template <typename ST>
ST frand(ST fMin, ST fMax)
{
    ST f = ST(rand()) / ST(RAND_MAX);
    return fMin + f * (fMax - fMin);
}

SUITE( gsExportMatrix )
{

    template <typename MatrixT>
    void testImpl()
    {
        typedef typename MatrixT::Scalar ST;
        typedef typename MatrixT::Index  IT;

        srand(static_cast<unsigned>(time(NULL)));

        IT cols=3;
        IT rows=4;

        MatrixT exported(rows,cols), imported;
        for (IT i = 0; i < rows; ++i )
            for (IT j = 0; j < cols; ++j )
                exported(i,j)= frand<ST>(-1,1);

        exportMatrixToASCII("/tmp/m.txt", exported );
        importMatrixFromASCII("/tmp/m.txt", imported );

        bool passed=true;
        for (IT i = 0; i < rows; ++i )
            for (IT j = 0; j < cols; ++j )
                passed = passed && math::abs(exported(i,j)-imported(i,j)) < std::numeric_limits<ST>::epsilon();
        //std::cout << "\n\nexported:\n" << exported << "\n\nimported:\n" << imported << "\n" << std::endl;
        CHECK( passed );
    }

    TEST( exportWithgsSparseMatrix )
    {
        testImpl< gsSparseMatrix<real_t,Eigen::RowMajor,index_t> >();
        // existing test
    }
    TEST( exportWithgsSparseMatrixColMajor )
    {
        testImpl< gsSparseMatrix<real_t, Eigen::ColMajor, index_t> >();
    }

    TEST( exportWithDenseMatrices )
    {
        testImpl< gsMatrix<real_t> >();
    }

    TEST( inlineConstruction )
    {
        gsMatrix<> m1,m2(2,3);
        m2<<0,1,2,2,3,4;
        std::stringstream s("0 1 2\n2 3 4\n");
        importMatrixFromASCII(s,m1);
        CHECK( gsAllCloseAbsolute(m1,m2,1e-15) );
    }
}
