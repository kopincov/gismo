/** @file gsPlotNormalAngles 

    @brief Plots normal angles of the curve. 

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

int main(int argc, char *argv[])
{

    std::string inputFile(MOTOR_DATA_DIR "jku/target_input.xml");
    std::string outputFile("angles.txt");
    int numParams = 100;

    gsCmdLine cmd("Plotting the normal angles of the curve.");
    cmd.addString("f", "file", "Name of input file", inputFile);
    cmd.addString("o", "output", "Name of output file", outputFile);
    cmd.addInt("n", "n", "Number of parameters.", numParams);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input file:        " << inputFile << "\n"
           << "Output file:       " << outputFile << "\n"
           << "N:                 " << numParams << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsFileData<> fd(inputFile);
    
    gsGeometry<>* geom = fd.getAnyFirst< gsGeometry<> >().release();
    gsMatrix<> parameters = uniformParameters<2> (0.0, 1.0, numParams);
    gsMatrix<> derivs = geom->deriv(parameters);
    gsMatrix<> normals = rotateVectors(derivs, -EIGEN_PI/2.0);
    normalize(normals);
    gsMatrix<> angles = getAngles(normals);

    writeToTxtFile(parameters, angles, outputFile);
    
    return 0;
}
