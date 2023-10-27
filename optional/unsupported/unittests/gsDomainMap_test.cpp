/** @file gsDomainMap_test.cpp

    @brief Tests gsDomainMap.cpp.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
**/

#include "gismo_unittest.h"
#include <gsRemappedBasis/gsDomainMap.h>
#include <gsRemappedBasis/gsSelector.h>

SUITE(gsDomainMap)
{
    TEST(polytopeUnion)
    {
        gsSelector sel;
        gsMatrix<real_t> boundingBox(2,2);
        boundingBox<<0,1,0,1;

        gsBoxList bb0(2);
        gsMatrix<real_t> dom(2,2);

        dom<<0,0.4,0.3,1;
        bb0.append(dom,1);
        dom<<0.4,1,0,0.7;
        bb0.append(dom,1);

        sel.initFromBoxesMax(bb0,boundingBox,0);

        gsBoxList bb1(2);
        dom<<0.2,0.6,0.2,0.8;
        bb1.append(dom,1);

        sel.initFromBoxesMax(bb1,boundingBox,1);

        sel.patch(2) = polytopeUnion(sel.patch(0),sel.patch(1));
        sel.patch(3) = polytopeIntersect(sel.patch(0),sel.patch(1));
        sel.patch(4) = polytopeDifference(sel.patch(0),sel.patch(1));

        gsMatrix<real_t> points(2,20);
        std::vector<gsBoxList::basisIdT> testLevelsAfterUnion;
        std::vector<gsBoxList::basisIdT> testLevelsAfterIntersect;
        std::vector<gsBoxList::basisIdT> testLevelsAfterDifference;
        points << .1, .1, .1, .1, .1, .3, .3, .3, .3, .3, .5, .5, .5, .5, .5, .7, .7, .7, .7, .7,
                  .15,.25,.55,.75,.95,.15,.25,.55,.75,.95,.15,.25,.55,.75,.95,.15,.25,.55,.75,.95;
        sel.getBasisAt(points,testLevelsAfterUnion,2);
        sel.getBasisAt(points,testLevelsAfterIntersect,3);
        sel.getBasisAt(points,testLevelsAfterDifference,4);
        real_t corrLevelsAfterUnion[] =
        {
            0, 0, 1, 1, 1,
            0, 1, 1, 1, 1,
            1, 1, 1, 1, 0,
            1, 1, 1, 0, 0
        };
        real_t corrLevelsAfterIntersect[] =
        {
            0, 0, 0, 0, 0,
            0, 0, 1, 1, 0,
            0, 1, 1, 0, 0,
            0, 0, 0, 0, 0
        };
        real_t corrLevelsAfterDifference[] =
        {
            0, 0, 1, 1, 1,
            0, 0, 0, 0, 1,
            1, 0, 0, 0, 0,
            1, 1, 1, 0, 0
        };
        CHECK( testLevelsAfterUnion.size()      == 20 );
        CHECK( testLevelsAfterIntersect.size()  == 20 );
        CHECK( testLevelsAfterDifference.size() == 20 );
        for( size_t i = 0; i < 20; ++i )
        {
            //std::cout << testLevelsAfterUnion[i]      << " ";
            //std::cout << testLevelsAfterIntersect[i]  << " ";
            //std::cout << testLevelsAfterDifference[i] << std::endl;
            CHECK( testLevelsAfterUnion[i]      == corrLevelsAfterUnion[i]      );
            CHECK( testLevelsAfterIntersect[i]  == corrLevelsAfterIntersect[i]  );
            CHECK( testLevelsAfterDifference[i] == corrLevelsAfterDifference[i] );
        }

        //sel.exportToTex("/tmp/draw");
        //system("cd /tmp;pdflatex draw.tex && okular draw.pdf &");
    }
}
