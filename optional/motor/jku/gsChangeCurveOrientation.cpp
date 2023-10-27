/** @file gsChangeCurveOrientation.cpp

    @brief This program changes curve orientation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/

#include <gismo.h>
#include "gsMotorUtils.h"

using namespace gismo;


int main(int argc, char *argv[])
{
    // Options with defaul t values
    std::string inputFile(MOTOR_DATA_DIR "jku/target_input.xml");
    std::string outputFile("target_input_changed_orientation");

    gsCmdLine cmd("Changes curve orientation");
    cmd.addString("f", "input", "Name of input file", inputFile);
    cmd.addString("o", "output", "Name of output file", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input file:  " << inputFile << "\n"
           << "Output file: " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsGeometry<>* p_geometry = readGeometry(inputFile);
    gsGeometry<>& geometry = *p_geometry;

    gsBSpline<> inputGeometry = dynamic_cast< const gsBSpline<>& >(geometry);

    gsMatrix<> ctrlPts = inputGeometry.coefs();

    gsMatrix<> reversedCtrlPts = reverseMatrixRows(ctrlPts);

    gsKnotVector<> KV = inputGeometry.knots(); // (0.0, 1.0,    num-deg-1, deg+1);

    // Constructing ellipse as a B-Spline
    gsBSpline<> bsc(KV, reversedCtrlPts);

    gsInfo << "Writing to ParaView..." << std::endl;
    gsWrite(bsc, outputFile);
    gsWriteParaview(bsc, outputFile);

    return 0;
}
