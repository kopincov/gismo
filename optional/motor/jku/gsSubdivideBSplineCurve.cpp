/** @file gSubdivideBSplineCurve.cpp

    @brief ....

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/

#include <gismo.h>
#include "gsMotorIOUtils.h"
#include "gsMotorUtils.h"
#include "gsUtils/gsExportMatrix.h"

using namespace gismo;

int main(int argc, char *argv[])
{

    // Options with default values
    std::string inputFile(MOTOR_DATA_DIR "jku/target_input.xml");
    std::string inputPoints(MOTOR_DATA_DIR "jku/points_target.txt");
    std::string outputFile("subd");

    gsCmdLine cmd("...");
    cmd.addString("f", "input", "Name of input file with geometry", inputFile);
    cmd.addString("p", "points", "Name of input file with point coordinates", inputPoints);
    cmd.addString("o", "output", "Name of output file", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input file (geometry): " << inputFile << "\n"
           << "Input file (points):   " << inputPoints << "\n"
           << "Output file:           " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    if (inputFile == "")
    {
        gsInfo << "Didn't specify the input file." << std::endl;
        return 1;
    }

    gsFileData<> fd(inputFile);
    gsGeometry<>* p_geometry = fd.getAnyFirst< gsGeometry<> >().release();
    gsGeometry<>& geometry = *p_geometry;

    gsBSpline<> inputGeometry = dynamic_cast< const gsBSpline<>& >(geometry);
    gsMatrix<> ctrlPts = inputGeometry.coefs();
    int deg = inputGeometry.degree();
    gsKnotVector<> kv = inputGeometry.knots();

/*
    // Making geometry
    const int num = 10;
    const int deg = 2;
    gsMatrix<> ctrlPts(num, 2);
    ctrlPts.setZero();
    for (int i = 0; i != num; i++)
    {
        ctrlPts(i, 0) = math::cos((2.0*EIGEN_PI*i)/(num-1));
        ctrlPts(i, 1) = math::sin((2.0*EIGEN_PI*i)/(num-1));
    }
    make_g1(ctrlPts);
    gsKnotVector<> kv(0.0, 1.0, num-deg-1, deg+1);
    gsBSpline<> inputGeometry(kv, ctrlPts);
*/
    // /////////////////////////////

    // Selecting point and finding corresponding parameter values
    gsMatrix<> pts;
    importMatrixFromASCII (inputPoints, pts);
    gsMatrix<> pars(1, pts.cols());

    // For testing //
//    pars(0, 0) = 0.25/2.0;
//    pars(0, 1) = 0.25/2.0 + 0.25;
//    pars(0, 2) = 0.25/2.0 + 2*0.25;
//    pars(0, 3) = 0.25/2.0 + 3*0.25;

//    gsMatrix<> pts2 = inputGeometry.eval(pars);
//    exportMatrixToASCII ("points_circle.txt", pts2);
//    return 0;

    for (int i = 0; i < pts.cols(); i++)
    {
        pars.col(i) = findParameter(inputGeometry, pts.col(i), 1.0e-12);
    }

    // Adding knots
    for (int i = 0; i < pars.cols(); i++)
    {
        gsBoehm(kv, ctrlPts, pars(i), kv.degree());
    }

    // Extracting curves
    int i_start = 0;
    int i_end = 0;

    gsMultiPatch<> mp;
    for (int i = 0; i < pars.cols()-1; i++)
    {
        index_t ind1 = 0;
        for(ind1 = 0; kv.at(ind1) < pars(i); ind1++);
        ind1--;
        if (i == 0)
            i_start = ind1;

        index_t ind2 = 0;
        for(ind2 = 0; kv.at(ind2) < pars(i+1); ind2++);
        ind2--;
        if (i == pars.cols()-2)
            i_end = ind2;

        gsMatrix<> ctrlPts2 = ctrlPts.block(ind1, 0, ind2-ind1+1, 2);
        gsKnotVector<> kv2(0.0, 1.0, ctrlPts2.rows()-deg-1, deg+1);

        gsBSpline<> extractedGeom(kv2, ctrlPts2);
        mp.addPatch(extractedGeom);
    }

    // Extracting a curve to which the start point belongs
    gsMatrix<> ctrlPts3(ctrlPts.rows()-i_end-2 + i_start+2, 2);
    ctrlPts3.setZero();
    ctrlPts3.block(0, 0, ctrlPts.rows()-i_end-1, 2) = ctrlPts.block(i_end, 0, ctrlPts.rows()-i_end-1, 2);
    ctrlPts3.block(ctrlPts.rows()-i_end-1, 0, i_start+1, 2) = ctrlPts.block(0, 0, i_start+1, 2);
    gsKnotVector<> kv3(0.0, 1.0, ctrlPts3.rows()-deg-1, deg+1);

    gsBSpline<> extractedGeom2(kv3, ctrlPts3);
    mp.addPatch(extractedGeom2);

    writeMultipatch(mp, outputFile);

    return 0;
}

