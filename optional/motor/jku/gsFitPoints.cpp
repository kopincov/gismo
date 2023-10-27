/** @file gsFitPoints.cpp

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
    // This code should be revised (problems with cdash)
    return 0;

//    // Options with defaul t values
//    std::string inputFile(MOTOR_DATA_DIR "jku/Austria.txt");
//    int numCtrlPts = 100;
//    int deg = 3;
//    real_t lambda = 1e-6;
//    int numIterations = 7;
//    std::string outputFile("Austria");

//    gsCmdLine cmd("Construction of minimal length (\"rubber band\") offset for smooth curve");
//    cmd.addString("f", "input", "Name of input file", inputFile);
//    cmd.addInt("N", "ncpts", "Number of control points", numCtrlPts);
//    cmd.addInt("p", "degree", "...", deg);
//    cmd.addReal("l", "lambda", "Smoothing parameter", lambda);
//    cmd.addInt("i", "numIter", "Number of refinement iterations.", numIterations);
//    cmd.addString("o", "output", "Name of output file", outputFile);
//    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

//    /*
//    if (numPoints < numCtrlPts)
//    {
//        gsWarn << "Number of sampling points must be not less than number of control points. Number of sampling points (numPoints) is set to " << numCtrlPts << "." << std::endl;
//    }
//    */
//    gsInfo << "Input arguments:\n"
//           << "Input file:  " << inputFile << "\n"
//           << "numCtrlPts:  " << numCtrlPts << "\n"
//           << "degree:      " << deg << "\n"
//           << "lambda:      " << lambda << "\n"
//           << "numIter:     " << numIterations << "\n"
//           << "Output file: " << outputFile << "\n"
//           << "--------------------------------------------------\n" << std::endl;

//    /*
//    if (inputFile == "")
//    {
//        gsInfo << "Didn't specify the input file." << std::endl;
//        return 0;
//    }
//    */

//    gsMatrix<> pts_tmp;
//    importMatrixFromASCII (inputFile, pts_tmp);
//    gsMatrix<> pts = pts_tmp.transpose();
//    const int numPoints = pts.cols();

//    gsMatrix<> params(1, numPoints);

//    real_t length = 0.0;
//    params(0, 0) = 0.0;
//    for(index_t i = 1; i != numPoints; i++)
//    {
//        gsMatrix<> tmp(2, 1);
//        tmp(0, 0) = pts(0, i) - pts(0, i-1);
//        tmp(1, 0) = pts(0, i) - pts(1, i-1);
//        length += tmp.squaredNorm();
//        params(0, i) = length;
//    }

//    // Scaling
//    for(index_t i = 1; i != numPoints; i++)
//    {
//        params(0, i) /= length;
//    }

//    gsWriteParaviewPoints(pts, "Austria_points");

//    gsKnotVector<> kv ( 0, 1,  numPoints-deg-1, deg+1 );
//    gsBSplineBasis<> basis( kv );

// /*    for(int i = 0; i != numIterations; i++)
//    {
//        basis.uniformRefine();
//    }
// */
//    gsFitting<> fitting(params, pts, basis);
//    fitting.compute(lambda);

//    gsGeometry<>* outputGeometry = fitting.result();

//    writeGeometry(*outputGeometry, outputFile);

//    return 0;
}

