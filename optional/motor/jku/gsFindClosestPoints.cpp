/** @file gsFindClosestPoints.cpp

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
    std::string inputFile1(MOTOR_DATA_DIR "jku/target_rectangle.xml");
    std::string inputFile2(MOTOR_DATA_DIR "jku/target_input_offset.xml");
    std::string inputPoints(MOTOR_DATA_DIR "jku/points_target_rectangle.txt");
    std::string outputFile("points_target_input_offset");

    gsCmdLine cmd("...");
    cmd.addString("A", "input1", "Name of input file with geometry (on which points are given)", inputFile1);
    cmd.addString("B", "input2", "Name of input file with geometry (on which points will be find)", inputFile2);
    cmd.addString("p", "points", "Name of input file with point coordinates", inputPoints);
    cmd.addString("o", "output", "Name of output file with point coordinates", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input file 1:        " << inputFile1 << "\n"
           << "Input file 2:        " << inputFile2 << "\n"
           << "Input file (points): " << inputPoints << "\n"
           << "Output file:         " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsFileData<> fd2(inputFile2);
    gsGeometry<>* p_geometry2 = fd2.getAnyFirst< gsGeometry<> >().release();
    gsGeometry<>& geometry2 = *p_geometry2;
    gsBSpline<> inputGeometry2 = dynamic_cast< const gsBSpline<>& >(geometry2);

    gsMatrix<> pts;
    importMatrixFromASCII (inputPoints, pts);

    gsMatrix<> pars(1, pts.cols());

    for (int i = 0; i < pts.cols(); i++)
    {
        pars.col(i) = findParameter(inputGeometry2, pts.col(i), 1.0e-12);
    }

    gsMatrix<> pts2 = inputGeometry2.eval(pars);

    gsInfo << "Writing " << outputFile << ".txt" << std::endl;
    exportMatrixToASCII (outputFile, pts2);

    return 0;
}

