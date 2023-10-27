/** @file makeTriangle.cpp

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
    real_t x1 = 0.0;
    real_t y1 = 0.0;
    real_t x2 = 1.0;
    real_t y2 = 0.0;
    real_t x3 = 0.0;
    real_t y3 = 1.0;

    std::string outputFile("target_input_triangle");
    bool multipatch = false;

    int num = 4;

    gsCmdLine cmd("...");
    cmd.addReal("A", "x1", "...", x1);
    cmd.addReal("B", "y1", "...", y1);
    cmd.addReal("C", "x2", "...", x2);
    cmd.addReal("D", "y2", "...", y2);
    cmd.addReal("E", "x3", "...", x3);
    cmd.addReal("F", "y3", "...", y3);
    cmd.addString("o", "output", "Name of output file", outputFile);
    cmd.addSwitch("multipatch", "Returns a multipatch, if it is set to true", multipatch);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "x1:         " << x1 << "\n"
           << "y1:         " << y1 << "\n"
           << "x2:         " << x2 << "\n"
           << "y2:         " << y2 << "\n"
           << "x3:         " << x3 << "\n"
           << "y3:         " << y3 << "\n"
           << "output:     " << outputFile << "\n"
           << "multipatch: " << multipatch << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsMatrix<> ctrlPts(num, 2);
    ctrlPts.setZero();
    ctrlPts << x1, y1,
               x2, y2,
               x3, y3,
               x1, y1;

    if (multipatch)
    {
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(0).transpose(), ctrlPts.row(1).transpose()));
        mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(1).transpose(), ctrlPts.row(2).transpose()));
        mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(2).transpose(), ctrlPts.row(3).transpose()));
        writeMultipatch(mp, outputFile);
    }
    else
    {
        gsKnotVector<> KV (0.0, 1.0, num-2, 2);
        gsBSpline<> bsp(KV, ctrlPts);
        writeGeometry(bsp, outputFile);
    }

    writeControlPoints(ctrlPts, outputFile + "_CtrlPts");

    return 0;
}
