/** @file gsTemplate_inner.cpp

    @brief This program constructs a template which consists of circle and segmentation lines.
    Segmentation lines are sides of rectangle and segments joining rectangle corners and circle.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/


#include <iostream>
#include <gismo.h>
#include <gsCore/gsMath.h>
#include <math.h>
#include "gsMotorUtils.h"
//#include "gsMotorIOUtils.h"

using namespace gismo;

int main(int argc, char *argv[])
{

    // Options with default values

    // Data for ellipsis 
    real_t r1 = 1.0;
    real_t r2 = 1.0;
    real_t x_c = 0.0;
    real_t y_c = 0.0;
    int num = 10;
    const int deg = 2;
    
    // Data for inner rectangle
    real_t x1 = -0.5;
    real_t y1 = -0.5;
    real_t x2 = 0.5;
    real_t y2 = 0.5;
    real_t rot = 0.0;

    int npts  = 100;


    
    std::string outputPrefix("template");
    gsCmdLine cmd("This program constructs a template which consists of circle and segmentation lines");
    cmd.addReal("x", "x", "The x-coordinate of the center", x_c);
    cmd.addReal("y", "y", "The y-coordinate of the center", y_c);
    cmd.addReal("A", "r1", "r1", r1);
    cmd.addReal("B", "r2", "r2", r2);
    cmd.addReal("C", "x1", "x1", x1);
    cmd.addReal("D", "y1", "y1", y1);
    cmd.addReal("E", "x2", "x2", x2);
    cmd.addReal("F", "y2", "y2", y2);
    cmd.addReal("G", "rot", "Rotation angle", rot);
    cmd.addInt("N", "num", "Number of control points", num);
    cmd.addInt("s", "npts", "Number of points used for sampling each patch (used for ParaView visualization) ", npts);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "x_c:         " << x_c << "\n"
           << "y_c:         " << y_c << "\n"
           << "r1:          " << r1 << "\n"
           << "r2:          " << r2 << "\n"
           << "x1:          " << x1 << "\n"
           << "y1:          " << y1 << "\n"
           << "x2:          " << x2 << "\n"
           << "y2:          " << y2 << "\n"
           << "rot:         " << rot << "\n"
           << "npts:        " << npts << "\n"
           << "Output file: " << outputPrefix << "\n"
           << "num:         " << num << "\n"
           << "--------------------------------------------------\n" << std::endl;

    //gsMultiPatch<> mp1; // Multi-patch for template goemetry (circle)
    gsMultiPatch<> mpSegmentation; // Multi-patch for segmentation lines

    // Ellipsis
    /*
    gsNurbs<>* p_ellipsis = gsNurbsCreator<>::NurbsEllipsis(r, 0.0, 0.0).release();
    gsNurbs<>& ellipsis = *p_ellipsis;
    */



    gsBSpline<> ellipsis = constructEllipse(x_c, y_c, r1, r2, num, deg);

       
    // Rotating (optional)
    if (rot != 0)
    {
        gsMatrix<> coefs = (ellipsis.coefs()).transpose();
        coefs = rotateVectors(coefs, rot).transpose();
        // Update control points
        ellipsis.setCoefs(coefs);
    }

    //mp1.addPatch(circle);

    // Segmentation lines

    // Rectangle vertices
    //gsKnotVector<> KV2 (0.0, 1.0, n-2, 2);
    gsMatrix<> ptsRectangle(2, 4);
    ptsRectangle.setZero();
    ptsRectangle << x2, x1, x1, x2,
                    y2, y2, y1, y1;
    // Rotating (optional)
    if (rot != 0)
    {
        ptsRectangle = rotateVectors(ptsRectangle, rot);
    }

    // Lines from rectangle vertices to ellipsis
    gsMatrix<> param(1, 4);
    param << 0., 0.+0.25, 0.+2*0.25, 0.+3*0.25;
    //param << 0.125, 0.125+0.25, 0.125+2*0.25, 0.125+3*0.25;
    //param << 0.125+3*0.25, 0.125, 0.125+0.25, 0.125+2*0.25;

    gsMatrix<> ptsEllipsis = ellipsis.eval(param);
    mpSegmentation.addPatch(gsNurbsCreator<>::BSplineLineSegment(ptsEllipsis.col(0), ptsRectangle.col(0)));
    mpSegmentation.addPatch(gsNurbsCreator<>::BSplineLineSegment(ptsEllipsis.col(1), ptsRectangle.col(1)));
    mpSegmentation.addPatch(gsNurbsCreator<>::BSplineLineSegment(ptsEllipsis.col(2), ptsRectangle.col(2)));
    mpSegmentation.addPatch(gsNurbsCreator<>::BSplineLineSegment(ptsEllipsis.col(3), ptsRectangle.col(3)));

    // Rectangle edges
    mpSegmentation.addPatch(gsNurbsCreator<>::BSplineLineSegment(ptsRectangle.col(3), ptsRectangle.col(0)));
    mpSegmentation.addPatch(gsNurbsCreator<>::BSplineLineSegment(ptsRectangle.col(0), ptsRectangle.col(1)));
    mpSegmentation.addPatch(gsNurbsCreator<>::BSplineLineSegment(ptsRectangle.col(1), ptsRectangle.col(2)));
    mpSegmentation.addPatch(gsNurbsCreator<>::BSplineLineSegment(ptsRectangle.col(2), ptsRectangle.col(3)));

    // Output xml files
    writeGeometry(ellipsis, outputPrefix + "_geometry");
    writeMultipatch(mpSegmentation, outputPrefix + "_segmentation");

    // Writing coordinates of segmentation points
    gsFileData<> fd;
    fd << ptsEllipsis;
    fd.dump(outputPrefix + "_segmentation_points");

    return 0;
}



gsNurbs<>* NurbsEllipsis( const real_t r1, const real_t r2, const real_t x, const real_t y)
{
    gsNurbs<>* ellipsis = gsNurbsCreator<>::NurbsCircle(1, 0, 0).release();

    // Scaling control points
    if ( (r1 != 1.0) && (r1 != 1.0))
    {
        gsMatrix<> C = ellipsis->coefs();

        // x direction
        // C0, C4, C8
        C(0, 0) *= r1;
        C(4, 0) *= r1;
        C(8, 0) *= r1;
        // y direction
        // C2, C6
        C(2, 1) *= r2;
        C(6, 1) *= r2;
        // both directions
        // C1, C3, C5, C7
        C(1, 0) *= r1;
        C(1, 1) *= r2;
        C(3, 0) *= r1;
        C(3, 1) *= r2;
        C(5, 0) *= r1;
        C(5, 1) *= r2;
        C(7, 0) *= r1;
        C(7, 1) *= r2;

        C.col(0).array() += x;
        C.col(1).array() += y;

        // Update control points
        ellipsis->setCoefs(C);
    }

    return ellipsis;
}

