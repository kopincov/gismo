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
//    int numURef = 0;
    std::string outputFile("target_input_LShaped_domain");

    int num = 5;

    gsCmdLine cmd("...");
    cmd.addReal("A", "x1", "...", x1);
    cmd.addReal("B", "y1", "...", y2);
    cmd.addReal("C", "x2", "...", x2);
    cmd.addReal("D", "y2", "...", y2);
//    cmd.addInt("F", "numURef", "...", numURef);
    cmd.addString("o", "output", "Name of output file", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;

    gsInfo << "Input arguments:\n"
           << "x1:     " << x1 << "\n"
           << "y1:     " << y1 << "\n"
           << "x2:     " << x2 << "\n"
           << "y2:     " << y2 << "\n"
           //<< "numURef: " << numURef << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsMatrix<> ctrlPts(num, 2);
    ctrlPts.setZero();
    ctrlPts << x1, y1,
               x2, y1,
               x2, y2,
               x1, y2,
               x1, y1;

    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(0).transpose(), ctrlPts.row(1).transpose()));
    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(1).transpose(), ctrlPts.row(2).transpose()));
    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(2).transpose(), ctrlPts.row(3).transpose()));
    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(ctrlPts.row(3).transpose(), ctrlPts.row(4).transpose()));

    gsFileData<> fd;
    gsInfo << "Writing " << outputFile + ".xml" << "\n";
    fd << mp;
    fd.dump(outputFile);

    gsInfo << "Writing to ParaView..." << std::endl;
    gsWriteParaview(mp, outputFile, 100);
    writeControlPoints(ctrlPts, outputFile + "CtrlPts");

    return 0;
}
