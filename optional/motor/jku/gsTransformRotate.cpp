/** @file gsTransformRotate.cpp

    @brief ...

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/

#include <gismo.h>
#include "gsMotorUtils.h"
//#include "gsMotorIOUtils.h"
#include "gsUtils/gsExportMatrix.h"

using namespace gismo;

int main(int argc, char *argv[])
{

    // Options with defaul t values
    std::string inputFile(MOTOR_DATA_DIR "jku/target_input.xml");
    real_t rot   = EIGEN_PI/4.0; // Default rotation degree value is -Pi/4
    int numPoints = 1000;
    std::string outputFile("rotated_curve");

    gsCmdLine cmd("Transforming (rotating ) a B-Spline curve given in a file");
    cmd.addString("f", "input", "Name of input file", inputFile);
    cmd.addReal("r", "rot", "Rotation degree", rot);
    cmd.addInt("n", "npts", "Number of points", numPoints);
    cmd.addString("o", "output", "Name of output file", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input file:  " << inputFile << "\n"
           << "rot:         " << rot << "\n"
           << "Output file: " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    if (inputFile == "")
    {
        gsInfo << "Didn't specify the input file." << std::endl;
        return 0;
    }

    gsGeometry<>* p_geometry = readGeometry(inputFile);
    gsGeometry<>& geometry = *p_geometry;

    gsBSpline<> inputGeometry = dynamic_cast< const gsBSpline<>& >(geometry);

    gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints);
    // Points on the curve
    gsMatrix<> points = inputGeometry.eval(parameters);

    // Rotating the points on the curve
    gsMatrix<> transformedPoints = rotateVectors(points, rot);

    gsKnotVector<> kv = inputGeometry.knots();

    gsBSplineBasis<> basis( kv );
    //basis.uniformRefine( (1<<numURef)-1 );

    gsFitting<> fitting(parameters, transformedPoints, basis);
    fitting.compute();

    writeGeometry(*fitting.result(), outputFile);

    return 0;

}

