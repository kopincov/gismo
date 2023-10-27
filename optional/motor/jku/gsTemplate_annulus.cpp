/** @file gsTemplate_annulus.cpp

    @brief This program constructs B-Spline representation of annulus
    with four segmentation lines and saves the data in XML and ParaView
    formats.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius, J. Speh
*/
// ./bin/gsTemplate_annulus -r 1.0 -R 2.0 -l 4 -n 30 -d 2 -s 100 -o 2d_template --save --plot

#include <iostream>
#include <gismo.h>
#include <gsCore/gsMath.h>
#include <math.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

int main(int argc, char *argv[])
{

    // Options with default values
    bool save = false; // If set to true, XML file is generated
    bool plot = false; // If set to true, ParaView file is generated and launched on exit
    real_t r1 = 1.0;
    real_t r2 = 2.0;
    int l     = 4;
    //const double angle = 2*EIGEN_PI/l; // Angle between segmentation lines
    int n     = 10;
    int deg   = 2;
    int npts  = 100;
    std::string outputPrefix("template");
    gsCmdLine cmd("This program constructs B-Spline representation of annulus\n and saves it in XML and ParaView formats.");
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addReal("r", "radius1", "Radius of the first circle", r1);
    cmd.addReal("R", "radius2", "Radius of the second circle", r2);
    cmd.addInt("l", "lines", "Number of segmentation lines", l);
    cmd.addInt("n", "num", "Number of control points representing each circle", n);
    cmd.addInt("d", "degree", "Degree of B-Spline", deg);
    cmd.addInt("s", "npts", "Number of points used for sampling each patch (used for ParaView visualization) ", npts);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Two circles
    gsMultiPatch<> mp1;

    gsKnotVector<> KV1 (0.0, 1.0, n-deg-1, deg+1);
    gsMatrix<> C(n,2);
    C.setZero();
    gsMatrix<> C_1(n,2);
    C_1.setZero();
    gsMatrix<> C_2(n,2);
    C_2.setZero();

    for(int i = 0; i < n; i++)
    {
        C(i, 0) = math::cos((2.0*EIGEN_PI*i)/(n-1));
        C(i, 1) = math::sin((2.0*EIGEN_PI*i)/(n-1));
    }

    // Control points of 1st circle
    C_1 = r1*C;
    make_g1(C_1);
    // Control points of 2nd circle
    C_2 = r2*C;
    make_g1(C_2);

    gsBSpline<>::uPtr circle_1(new gsBSpline<>(KV1, C_1));
    gsBSpline<>::uPtr circle_2(new gsBSpline<>(KV1, C_2));

    mp1.addPatch(circle_1->clone());
    mp1.addPatch(circle_2->clone());

    // Segmentation lines
    gsMultiPatch<> mp2;

    if(l > n-deg-1)
    {
        gsInfo << "Not enough control points for " << l << " segmentation lines" << "\n";
        return 1;
    }

    gsMatrix<> param(1, 4);
    param << 0.125, 0.125+0.25, 0.125+2*0.25, 0.125+3*0.25;

    gsMatrix<>  values_1 = circle_1->eval(param);
    gsMatrix<>  values_2 = circle_2->eval(param);

    gsMatrix<> C_l(2,2);
    C_l.setZero();
    for (int i = 0; i < l; i++)
    {
        mp2.addPatch(gsNurbsCreator<>::BSplineLineSegment(values_1.col(i), values_2.col(i)));
    }

    // Output an xml file
    if(save)
    {
        writeMultipatch(mp1, outputPrefix + "_geometry");
        writeMultipatch(mp2, outputPrefix + "_segmentation");
    }

    if(plot)
    {
        //Call paraview on exit
        if (mp1.nPatches() == 1)
        {
            std::string command = "paraview " + outputPrefix +"_geometry.vts &";
            return system(command.c_str());
        }
        else
        {
            std::string command = "paraview " + outputPrefix +"_geometry.pvd &";
            return system(command.c_str());
        }
    }
    else
    {
        return 0;
    }
}
