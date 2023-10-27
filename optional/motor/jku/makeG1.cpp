/** @file gsMakeG1.h

    @brief This program makes G1 geometry.

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

    // Options with default values
    std::string inputFile(MOTOR_DATA_DIR "jku/target_input.xml");
    std::string outputFile("target_input_G1");

    gsCmdLine cmd("This program makes G1 geometry");
    cmd.addString("f", "input", "Name of the input files", inputFile);
    cmd.addString("o", "output", "Name of the output files", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsGeometry<> * geom = NULL;
    gsFileData<> input(inputFile);
    if (input.has< gsGeometry<> >())
    {
        geom = input.getFirst< gsGeometry<> >().release();
    }
    else
    {
        gsInfo << "Input file " << inputFile << " doesn't have a geometry inside." << "\n";
        return -1;
    }
    gsInfo << "Control points (before):\n" << geom->coefs() << std::endl;
    make_g1(*geom);
    gsInfo << "Control points (after):\n" << geom->coefs() << std::endl;

    writeGeometry(*geom, outputFile);

    return 0;
}
