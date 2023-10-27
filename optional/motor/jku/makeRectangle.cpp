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
    real_t x1 = 0.0;
    real_t y1 = 0.0;
    real_t x2 = 1.0;
    real_t y2 = 1.0;
    std::string outputFile("target_input_rectangle");
    bool multipatch = false;

    int num = 9;

    gsCmdLine cmd("...");
    cmd.addReal("A", "x1", "...", x1);
    cmd.addReal("B", "y1", "...", y1);
    cmd.addReal("C", "x2", "...", x2);
    cmd.addReal("D", "y2", "...", y2);
    cmd.addString("o", "output", "Name of output file", outputFile);
    cmd.addSwitch("multipatch", "Returns a multipatch, if it is set to true", multipatch);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "x1:         " << x1 << "\n"
           << "y1:         " << y1 << "\n"
           << "x2:         " << x2 << "\n"
           << "y2:         " << y2 << "\n"
           << "output:     " << outputFile << "\n"
           << "multipatch: " << multipatch << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsMatrix<> ctrlPts(num, 2);
    ctrlPts.setZero();
    // Such configuration of control points is necessary for template segmentation
    ctrlPts << x2, (y2+y1)/2.0,
               x2, y2,
               (x2+x1)/2.0, y2,
               x1, y2,
               x1, (y2+y1)/2.0,
               x1, y1,
               (x2+x1)/2.0, y1,
               x2, y1,
               x2, (y2+y1)/2.0;

    if (multipatch)
    {
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(7).transpose(), ctrlPts.row(1).transpose()));
        mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(1).transpose(), ctrlPts.row(3).transpose()));
        mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(3).transpose(), ctrlPts.row(5).transpose()));
        mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(5).transpose(), ctrlPts.row(7).transpose()));
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
