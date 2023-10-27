/** @file gsTemplateSegmentationTriangles.cpp

    @brief This program constructs segmentation of triangle template.

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
    std::string inputFile(MOTOR_DATA_DIR "jku/template_multipatch_triangles.xml");
    std::string outputPrefix("template_segmentation");

    gsCmdLine cmd("...");
    cmd.addString("f", "input", "Name of the input file (multipatch)", inputFile);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

//    gsMultiPatch<>* mp = readMultipatch(inputFile);
//    const gsGeometry<>& curve_1 = (*mp)[0];
//    const gsGeometry<>& curve_2 = (*mp)[1];
//    gsMultiPatch<> mp_segm;

//    // Points on the inner triangle
//    const gsVector<> P1 = curve_1.coef(0);
//    const gsVector<> P2 = curve_1.coef(1);
//    const gsVector<> P3 = curve_1.coef(2);

//    gsVector<> P1_mid(2);
//    P1_mid(0) = 0.5*(P1(0)+P2(0));
//    P1_mid(1) = 0.5*(P1(1)+P2(1));

//    gsVector<> P2_mid(2);
//    P2_mid(0) = 0.5*(P2(0)+P3(0));
//    P2_mid(1) = 0.5*(P2(1)+P3(1));

//    gsVector<> P3_mid(2);
//    P3_mid(0) = 0.5*(P3(0)+P1(0));
//    P3_mid(1) = 0.5*(P3(1)+P1(1));

//    const real_t alpha_1 = 1.0/3.0;
//    const real_t alpha_2 = 1.0/3.0;
//    const real_t alpha_3 = 1.0/3.0;

//    gsVector<> P_center(2);
//    P_center << alpha_1*P1 + alpha_2*P2 + alpha_3*P3;

//    mp->addPatch(gsNurbsCreator<>::BSplineLineSegment(P_center, P1_mid));
//    mp->addPatch(gsNurbsCreator<>::BSplineLineSegment(P_center, P2_mid));
//    mp->addPatch(gsNurbsCreator<>::BSplineLineSegment(P_center, P3_mid));

//    mp_segm.addPatch(gsNurbsCreator<>::BSplineLineSegment(P_center, P1_mid));
//    mp_segm.addPatch(gsNurbsCreator<>::BSplineLineSegment(P_center, P2_mid));
//    mp_segm.addPatch(gsNurbsCreator<>::BSplineLineSegment(P_center, P3_mid));

//    // Points on the outer triangle
//    const gsVector<> P1_ = curve_2.coef(0);
//    const gsVector<> P2_ = curve_2.coef(1);
//    const gsVector<> P3_ = curve_2.coef(2);

//    mp->addPatch(gsNurbsCreator<>::BSplineLineSegment(P1, P1_));
//    mp->addPatch(gsNurbsCreator<>::BSplineLineSegment(P2, P2_));
//    mp->addPatch(gsNurbsCreator<>::BSplineLineSegment(P3, P3_));

//    mp_segm.addPatch(gsNurbsCreator<>::BSplineLineSegment(P1, P1_));
//    mp_segm.addPatch(gsNurbsCreator<>::BSplineLineSegment(P2, P2_));
//    mp_segm.addPatch(gsNurbsCreator<>::BSplineLineSegment(P3, P3_));

//    writeMultipatch(*mp, outputPrefix + "_segmentation_and_geometry");
//    writeMultipatch(mp_segm, outputPrefix + "_segmentation");

    return 0;
}

