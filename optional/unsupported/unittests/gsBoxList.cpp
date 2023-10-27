/** @file gsBoxList.cpp

    @brief Tests gsBoxList, which is a part of the ambitious
    gsRemappedBasis project.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
**/

#include "gismo_unittest.h"
#include <gsRemappedBasis/gsBoxList.h>



SUITE(gsRemappedBasis_gsBoxList)
{
    TEST( asGismoBox)
    {
        // Tests whether asGismoBox (called through asGismoBoxAll) returns the correct indices.

        // Prepare the gsTHBSplineBasis for comparison.
        gsKnotVector<real_t> kv(0, 1, 3, 4, 1, 3);
        gsTensorBSplineBasis<2,real_t> tensorBasis(kv,kv);
        gsTHBSplineBasis<2> thbBasis(tensorBasis);

        // We need to refine anyhow to have three levels;
        // however, it is more illustrative if we take the same boxes as later in the remapped case.
        index_t boxesA[] = { 1, 2, 2, 6, 6,
                              1, 4, 6, 6, 8,
                              1, 6, 4, 8, 8,
                              2, 6, 4, 8, 12,
                              2, 8, 4, 10, 16};

        std::vector<index_t> boxesTHB( boxesA, boxesA + sizeof(boxesA) / sizeof(index_t) );
        thbBasis.refineElements( boxesTHB );

        // Now create the boxes for the box list.
        gsBoxList boxList(2);
        gsMatrix<real_t> box0REM(2,2);
        box0REM << 0.25, 0.75, 0.25, 0.75;
        boxList.append(box0REM, 1);

        gsMatrix<real_t> box1REM(2,2);
        box1REM << 0.5, 0.75, 0.75, 1;
        boxList.append(box1REM, 1);

        gsMatrix<real_t> box2REM(2,2);
        box2REM << 0.75, 1, 0.5, 1;
        boxList.append(box2REM, 1);

        gsMatrix<real_t> box3REM(2,2);
        box3REM << 0.375, 0.5, 0.25, 0.75;
        boxList.append(box3REM, 2);

        gsMatrix<real_t> box4REM(2,2);
        box4REM << 0.5, 0.625, 0.25, 1;
        boxList.append(box4REM, 2);

        // This is the conversion we are checking.
        std::vector<index_t> result;
        boxList.toRefineElementFormat(thbBasis.getBases(), result);

        // Check whether the size and contents are the same (cheers to the STD Library).
        CHECK( result == boxesTHB );

        /* Las but not least, a story happening at the next table after a long debugging.
         * Jarle: "Andrea, but you haven't pulled. This is the old code."
         * Andrea: "But it is better." */
    }

    TEST(asBoxListBox)
    {
        // Tests whether asBoxListBox returns the correct values.

        // Prepare the gsTHBSplineBasis for comparison.
        gsKnotVector<real_t> kv(0, 1, 3, 4, 1, 3);
        gsTensorBSplineBasis<2,real_t> tensorBasis(kv,kv);
        gsTHBSplineBasis<2> thbBasis(tensorBasis);

        // We need to refine anyhow to have three levels;
        // however, it is more illustrative if we take the same boxes as later in the remapped case.
        unsigned boxesA[] = { 1, 2, 2, 6, 6,
                              1, 4, 6, 6, 8,
                              1, 6, 4, 8, 8,
                              2, 6, 4, 8, 12,
                              2, 8, 4, 10, 16};

        std::vector<index_t> boxesTHB( boxesA, boxesA + sizeof(boxesA) / sizeof(unsigned) );
        thbBasis.refineElements( boxesTHB );

        gsBoxList boxList(thbBasis.getBases(),thbBasis.dim(), boxesTHB);

        gsMatrix<real_t> box(2,2);
        box << 0.25, 0.75, 0.25, 0.75;
        CHECK( boxList.box(0) == box && boxList.basisId(0) == 1);

        box << 0.5, 0.75, 0.75, 1;
        CHECK( boxList.box(1) == box && boxList.basisId(1) == 1);

        box << 0.75, 1, 0.5, 1;
        CHECK( boxList.box(2) == box && boxList.basisId(2) == 1);

        box << 0.375, 0.5, 0.25, 0.75;
        CHECK( boxList.box(3) == box && boxList.basisId(3) == 2);

        box << 0.5, 0.625, 0.25, 1;
        CHECK( boxList.box(4) == box && boxList.basisId(4) == 2);
    }
}
