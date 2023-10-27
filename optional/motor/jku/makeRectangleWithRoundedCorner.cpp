/** @file makeRectangle.cpp

    @brief ...

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/

#include <gismo.h>
#include "gsMotorIOUtils.h"

using namespace gismo;

int main(int argc, char *argv[])
{
    // This code should be revised (problems with cdash)
    return 0;

//    real_t x1 = 0.0;
//    real_t y1 = 0.0;
//    real_t x2 = 1.0;
//    real_t y2 = 1.0;
//    std::string outputFile("target_input_rectangle_rounded");
//    bool multipatch = true;

//    int num = 7;

//    gsCmdLine cmd("...");
//    cmd.addReal("A", "x1", "...", x1);
//    cmd.addReal("B", "y1", "...", y2);
//    cmd.addReal("C", "x2", "...", x2);
//    cmd.addReal("D", "y2", "...", y2);
//    cmd.addString("o", "output", "Name of output file", outputFile);
//    cmd.addSwitch("multipatch", "Returns a multipatch, if it is set to true", multipatch);

//    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

//    gsInfo << "Input arguments:\n"
//           << "x1:     " << x1 << "\n"
//           << "y1:     " << y1 << "\n"
//           << "x2:     " << x2 << "\n"
//           << "y2:     " << y2 << "\n"
//           << "output:     " << outputFile << "\n"
//           << "multipatch: " << multipatch << "\n"
//           << "--------------------------------------------------\n" << std::endl;

//    gsMatrix<> ctrlPts(num, 2);
//    ctrlPts << x1, y1,
//               0.5*(x1+x2), y1,
//               x2, y1,
//               x2, 0.5*(y1+y2),
//               x2, y2,
//               x1, y2,
//               x1, y1;

//    int deg = 2;
//    gsKnotVector<> KV(0.0, 1.0, 5-deg-1, deg+1);

//    gsBSpline<> * bsc = new gsBSpline<>(KV, ctrlPts.block(0, 0, 5, 2));

//    gsMultiPatch<> mp;
//    mp.addPatch(bsc);
//    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(4).transpose(), ctrlPts.row(5).transpose()));
//    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(5).transpose(), ctrlPts.row(6).transpose()));
//    if (multipatch)
//    {
//        writeMultipatch(mp, outputFile);
//    }
//    else
//    {
//        gsInfo << "This case of separate geometries is not implemented (use multipatch switch)." << std::endl;
//        return -1;
// //        gsGeometry<> curve_1 = dynamic_cast< gsGeometry<> >(mp.patch(0));
// //        gsGeometry<> curve_2 = dynamic_cast< gsGeometry<> >(mp.patch(1));
// //        gsBSpline<>& curve_3 = dynamic_cast< gsBSpline<>& >(mp.patch(2));

// //        gsGeometry<>&  curve_1 = mp.patch(0);
// //        gsGeometry<>&  curve_2 = mp.patch(1);
// //        gsBSpline<> geom;
// //        writeGeometry(geom, outputFile);
//    }

//    writeControlPoints(ctrlPts, outputFile + "CtrlPts");

//    return 0;
}
