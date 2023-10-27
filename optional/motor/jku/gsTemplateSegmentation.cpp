/** @file gsTemplateSegmentation.cpp

    @brief This program constructs segmentation of the template geometry
    which consists of two boundary curves.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/


#include <iostream>
#include <sstream>
#include <gismo.h>
#include "gsMotorIOUtils.h"

using namespace gismo;

int main(int argc, char *argv[])
{

    // Options with default values
    int l = 4;
    std::string inputFile(MOTOR_DATA_DIR "jku/template_multipatch.xml");
    std::string outputPrefix("template");

    gsCmdLine cmd("...");
    cmd.addInt("l", "lines", "Number of segmentation lines", l);
    cmd.addString("f", "input", "Name of the input file (multipatch)", inputFile);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>* mp = readMultipatch(inputFile);
    const gsGeometry<>& curve1 = (*mp)[0];
    const gsGeometry<>& curve2 = (*mp)[1];
    gsMultiPatch<> mp_segm;

    gsMatrix<> param(1, l);
    const real_t angle = 1.0/l; // "Angle" between segmentation lines
    real_t mult  = 0.0;         // Auxilary variable
    if (l % 2 == 0)
    {
        mult = 0.5;
    }
    else
    {
        mult = 0.0;
    }
    for (int i = 0; i != l; i++)
    {
        param(i) = angle*(mult + i);
    }

    gsMatrix<>  points1 = curve1.eval(param);
    gsMatrix<>  points2 = curve2.eval(param);

    gsMatrix<> ctrlPts(2,2);
    ctrlPts.setZero();
    for (int i = 0; i != l; i++)
    {
        mp->addPatch(gsNurbsCreator<>::BSplineLineSegment(points1.col(i), points2.col(i)));
        mp_segm.addPatch(gsNurbsCreator<>::BSplineLineSegment(points1.col(i), points2.col(i)));
    }

    writeMultipatch(*mp, outputPrefix + "_segmentation_and_geometry");
    writeMultipatch(mp_segm, outputPrefix + "_segmentation");

    // Writing coordinates of segmentation points
    gsFileData<> fd;
    fd << points1;
    fd << points2;
    fd.dump(outputPrefix + "_segmentation_points");

    return 0;
}

