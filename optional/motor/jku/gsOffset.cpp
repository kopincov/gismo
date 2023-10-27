/** @file gsOffset.cpp

    @brief Construction of inward and outward offsets for smooth curve
    (WORK IN PROGRESS...)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius, J. Speh
*/

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

real_t replacement_f(const real_t arg,
                     const real_t c,
                     const bool inwardOffset);
gsMatrix<> modifyAnglesSvajunas(const gsMatrix<>& angles);
gsMatrix<> modifyAnglesJaka(const gsMatrix<>& angles);
gsMatrix<> modifyAngles3(const gsMatrix<>& angles);
gsMatrix<> modifyAngles_empty(const gsMatrix<>& angles);

int main(int argc, char *argv[])
{
    // Options with defaul t values
    std::string inputFile(MOTOR_DATA_DIR "jku/target_input.xml");
    int numPoints  = 1000;
    int numCtrlPts = 100;
    real_t c = 0.1;
    std::string outputFile("target_input_offset");
    bool inwardOffset = false;

    gsCmdLine cmd("Construction of offset for smooth curve");
    cmd.addString("f", "input", "Name of input file", inputFile);
    cmd.addInt("n", "npts", "Number of sampling points (for sampling)", numPoints);
    cmd.addInt("N", "ncpts", "Number of control points", numCtrlPts);
    cmd.addReal("c", "const", "...", c);
    cmd.addString("o", "output", "Name of output file", outputFile);
    cmd.addSwitch("inwardOffset", "...", inwardOffset);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (numPoints < numCtrlPts)
    {
        gsWarn << "Number of sampling points must be not less than number of control points. Number of sampling points (numPoints) is set to " << numCtrlPts << "." << std::endl;
    }
    gsInfo << "Input arguments:\n"
           << "Input file:  " << inputFile << "\n"
           << "numPoints:   " << numPoints << "\n"
           << "numCtrlPts:  " << numCtrlPts << "\n"
           << "Output file: " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsGeometry<>* p_geometry = readGeometry(inputFile);
    gsGeometry<>& geometry = *p_geometry;

    gsBSpline<> inputGeometry = dynamic_cast< const gsBSpline<>& >(geometry);
/*
    if (0 < numURef)
    {
        for (int r = 0; r != numURef; r++)
        {
            gsKnotVector<> KV = inputGeometry.knots();

            gsKnotVector<>::knotContainer newKnots;
            KV.getUniformRefinementKnots(1, newKnots);

            inputGeometry.insertKnots(newKnots.begin(), newKnots.end());
        }
    }
*/
    gsMatrix<> nrmlVects(numPoints, 2);  // Normal vectors
    nrmlVects.setZero();

    gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints); // different
    gsMatrix<> sampledPts = inputGeometry.eval(parameters);
    gsWriteParaviewPoints(sampledPts, "sampledPts");
    gsMatrix<> derivs = inputGeometry.deriv(parameters);
    gsMatrix<> derivs2 = inputGeometry.deriv2(parameters);

    if(inwardOffset)
        nrmlVects  = rotateVectors(derivs,  EIGEN_PI/2.0);
    else // Outward offset
        nrmlVects = rotateVectors(derivs, -EIGEN_PI/2.0);
    normalize(nrmlVects);

    /**************************************************/
    // Generating new normals
    gsMatrix<> nrmlVectsAngles = getAngles(nrmlVects);
    gsMatrix<> newNrmlVectsAngles = modifyAngles3(nrmlVectsAngles);
    gsMatrix<> newNrmlVects = getNormals(newNrmlVectsAngles);
    writeToTxtFile(parameters, newNrmlVectsAngles, "angles.txt");

    std::string out = outputFile + "_oldNormals";
    plotNormals(sampledPts, nrmlVects, out, 0.1);

    out = outputFile + "_newNormals";
    plotNormals(sampledPts, newNrmlVects, out, 0.1);
    
    /**************************************************/

    // Replacing sampled points
    gsMatrix<> curvs = getCurvatures(derivs, derivs2);
    //gsInfo << "Curvatures:\n" << curvs << std::endl;

    gsMatrix<> replacedPts (2, numPoints);
    for (int i = 0; i != numPoints; i++)
    {
        replacedPts(0, i) = sampledPts(0, i) + replacement_f(curvs(0, i), c, inwardOffset)*newNrmlVects(0, i);
        replacedPts(1, i) = sampledPts(1, i) + replacement_f(curvs(0, i), c, inwardOffset)*newNrmlVects(1, i);
    }
    gsWriteParaviewPoints(replacedPts, "replacedPts");
    writeDisplacements(sampledPts.transpose(), replacedPts.transpose(), "displacements");

    int deg = 2;
    gsKnotVector<> kv ( 0, 1,  numCtrlPts-deg-1, deg+1 );

    gsBSplineBasis<> basis( kv );

    gsFitting<> fitting(parameters, replacedPts, basis);
    fitting.compute();

    writeGeometry(*fitting.result(), outputFile);

    return 0;
}

real_t replacement_f(const real_t arg,
                     const real_t c,
                     const bool inwardOffset)
{
    return c;
}

gsMatrix<> modifyAnglesSvajunas(const gsMatrix<>& angles)
{
    const index_t numCols = angles.cols();
    gsMatrix<> newAngles(1, numCols);
    newAngles.setZero();

    real_t angle_max = angles(0, 0);
    real_t angle_tmp = angles(0, 0);
    const real_t angle_diff = EIGEN_PI/100000.0;

    newAngles(0, 0) = angles(0, 0);
    for (index_t i = 1; i != numCols; i++)
    {
        if (angles(0, i) < angles(0, i-1))
        { // Decreasing is detected
            if (angle_max - angles(0, i) > angle_diff)
            { // The angle decreased too much
                 newAngles(0, i) = angle_tmp;
            }
            else
            { // The angle decreased NOT too much
                angle_tmp = angles(0, i);
                newAngles(0, i) = angle_tmp;
            }
        }
        else
        { // Decreasing is NOT detected
            if(angles(0, i) > angle_tmp)
            {
                angle_max = angles(0, i);
                angle_tmp = angles(0, i);
                newAngles(0, i) = angles(0, i);
            }
            else
            {
                newAngles(0, i) = angle_tmp;
            }
        }
    }

    return newAngles;
}

gsMatrix<> modifyAnglesJaka(const gsMatrix<>& angles)
{
    const index_t numCols = angles.cols();
    const real_t start = angles(0, 0);
    const real_t end = angles(0, numCols - 1);

    gsMatrix<> newAngles(1, numCols);
    newAngles.setZero();

    for (index_t i = 0; i != numCols; i++)
    {
        const real_t frac = i * 1.0 / (numCols - 1);
        newAngles(0, i) = start * (1 - frac) + end * frac;
    }

    return newAngles;
}

gsMatrix<> modifyAngles3(const gsMatrix<>& angles)
{
    const index_t numCols = angles.cols();
    gsMatrix<> newAngles(1, numCols);
    newAngles.setZero();

    const gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numCols); // different

    gsGeometry<>::uPtr geom = fitAngles(parameters, angles, 100, 7);
    gsBSpline<> bsc = dynamic_cast< const gsBSpline<>& >(*geom);

    newAngles = bsc.eval(parameters);

    return newAngles;
}

gsMatrix<> modifyAngles_empty(const gsMatrix<>& angles)
{
    return angles;
}
