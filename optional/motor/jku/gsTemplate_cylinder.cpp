/** @file gsTemplate_cylinder.cpp

    @brief This program constructs B-Spline representation of cylinder
    and saves the data in XML and ParaView formats.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius, J. Speh
*/
// ./bin/gsTemplate_cylinder -f 2d_template_geometry.xml -s 2d_template_segmentation.xml -o 3d_template --save --plot
// ./bin/gsTemplate_cylinder -f 2d_target_geometry.xml -o 3d_target --save --plot

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
    bool save = true; // If set to true, XML file is generated
    bool plot = false; // If set to true, ParaView file is generated and launched on exit
    std::string inputFile(MOTOR_DATA_DIR "jku/template_geometry.xml");
    std::string segmFile(MOTOR_DATA_DIR "jku/template_segmentation.xml");
    std::string outputPrefix("3d_template");

    gsVector<real_t> levels(2);
    levels << 0.0, 1.0;

    gsCmdLine cmd("This program constructs B-Spline representation of cylinder\n and saves it in XML and ParaView formats.");
    //cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addString("f", "input", "Name of the template file (input)", inputFile);
    cmd.addString("s", "segmFile", "Name of the segmentation file (input)", segmFile);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>* mp = readMultipatch(inputFile);
    gsMultiPatch<>* mp_new = extrudeMultipatch(*mp, 1.0);

    /* **************************************************** */

    gsBSpline<>& curve_1 = dynamic_cast< gsBSpline<>& >(mp->patch(0));
    gsBSpline<>& curve_2 = dynamic_cast< gsBSpline<>& >(mp->patch(1));
    gsKnotVector<> knots_u = curve_1.knots();
    gsKnotVector<> knots_v(0.0, 1.0, 0, 2);
    int N_1 = curve_1.coefsSize();
    int N_2 = curve_2.coefsSize();
    for (int i = 0; i<2; i++)
    {
        gsMatrix<> C(N_1+N_2, 3);
        C.setConstant(levels[i]);
        C.block(0, 0, N_1, 2)   = curve_1.coefs();
        C.block(N_1, 0, N_2, 2) = curve_2.coefs();

        gsTensorBSplineBasis<2> BSplineBasis( knots_u, knots_v );
        gsTensorBSpline<2> BSpline( BSplineBasis, C );

        mp_new->addPatch(gsTensorBSpline<2>(BSpline));
    }

    // Segmentation planes
    if (segmFile != "")
    {
        gsMultiPatch<>* mp_lines = readMultipatch(segmFile);
        gsMultiPatch<>* mp_new2 = extrudeMultipatch(*mp_lines, 1.0);
        if(save)
        {
            writeMultipatch(*mp_new2, outputPrefix + "_segmentation");
        }
    }

    if(save)
    {
        writeMultipatch(*mp_new, outputPrefix + "_geometry");
    }

    delete mp;
    return 0;
}

