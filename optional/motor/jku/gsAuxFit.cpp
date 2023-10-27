/** @file gsAuxFit.cpp

    @brief This is auxiliary util for point clouds fitting

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

    // Options with default values
    std::string inputFile("SRM4+6.xml");
    int id = 0;
    std::string outputFile("output");

    gsCmdLine cmd("Reading a point cloud ...");
    cmd.addString("f", "input", "Name of input file with points clouds given in gsMatrix format", inputFile);
    cmd.addInt("i", "id", "id of the point cloud (matrix) in the input file", id);
    cmd.addString("o", "output", "Name of output file", outputFile);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input file:  " << inputFile << "\n"
           << "id:          " << id << "\n"
           << "Output file: " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsFileData <> fdIn(inputFile);
    gsMatrix<> points = *fdIn.getId<gsMatrix<> >(id);

    int numPoints = points.cols();

    gsMatrix<> parameters(1, numPoints);
    parameters.setZero();

    parameters(0, 0) = 0.0;
    // Chord length parameterization
    for (int i = 1; i != numPoints; i++)
    {
        parameters(0, i) = parameters(0, i-1) + (points.col(i) - points.col(i-1)).norm();
    }
    for (int i = 1; i != numPoints; i++)
    {
        parameters(0, i) /= parameters(0, numPoints-1);
    }

    // Constructing periodic fitting curve
    int deg = 2;
    int numKnots = 50;

    gsKnotVector<> kv = getPeriodicKnotVector(deg, numKnots);
    gsCurveFitting<> fitting(parameters.transpose(), points.transpose(), kv, true);
    fitting.compute_periodic();

    gsGeometry<>::uPtr geom = fitting.curve().clone();

    writeGeometry(*geom, outputFile);
    gsWriteParaview(*geom, outputFile);

    return 0;
}


