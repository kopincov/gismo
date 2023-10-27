/** @file makeEllipse.cpp

    @brief ...

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

int main(int argc, char *argv[])
{
    real_t r1 = 1.0;
    real_t r2 = 1.0;
    real_t x_c = 0.0;
    real_t y_c = 0.0;
    int num = 10;
    std::string outputFile("target_ellipse");
    int deg = 2;
    int nref = 0;
    
    gsCmdLine cmd("Construction of offset for smooth curve");
    cmd.addReal("r", "r1", "...", r1);
    cmd.addReal("R", "r2", "...", r2);
    cmd.addReal("x", "x", "The x-coordinate of the center", x_c);
    cmd.addReal("y", "y", "The y-coordinate of the center", y_c);
    cmd.addInt("N", "num", "Number of control points", num);
    cmd.addInt("d", "deg", "Degree", deg);
    cmd.addInt("n", "nref", "Number of refinement iterations", nref);    
    cmd.addString("o", "output", "Name of output file", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "r1:  " << r1 << "\n"
           << "r2:  " << r2 << "\n"
           << "x_c: " << x_c << "\n"
           << "y_c: " << y_c << "\n"
           << "N:   " << num << "\n"
           << "deg: " << deg << "\n"
           << "nref:" << nref << "\n"           
           << "--------------------------------------------------\n" << std::endl;

    gsBSpline<> ellipsis = constructEllipse(x_c, y_c, r1, r2, num, deg);

    for (int i = 0; i < nref; i++)
    {
        ellipsis.basis().uniformRefine_withCoefs(ellipsis.coefs());

    }
    
    gsInfo << "Writing to ParaView..." << std::endl;
    gsWrite(ellipsis, outputFile);
    gsWriteParaview(ellipsis, outputFile);
    writeControlPoints(ellipsis.coefs(), outputFile + "_CtrlPts");

    return 0;
}
