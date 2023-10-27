/** @file gsMap.cpp

    @brief This program applies a mapping (given in XML file).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsUtils/gsExportMatrix.h"

using namespace gismo;

int main(int argc, char *argv[])
{

    // Options with default values
    std::string dataFile(MOTOR_DATA_DIR "jku/template_segmentation.xml");
    std::string mapFile(MOTOR_DATA_DIR "jku/thb_map.xml");
    std::string outputFile("target_segmentation_optimized");
    int npts  = 100;

    gsCmdLine cmd("This program applies a mapping (given in XML file)");
    cmd.addString("f", "dataFile", "Name of the input (data) file", dataFile);
    cmd.addString("m", "mapFile", "Name of the input (map) file", mapFile);
    cmd.addInt("s", "npts", "Number of points used for sampling each patch (if data is given as multipatch)", npts);
    cmd.addString("o", "output", "The name of output file (prefix for all output files)", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input (data) file:              " << dataFile << "\n"
           << "Input(map) files:               " << mapFile << "\n"
           << "s:                              " << npts << "\n"
           << "Output file:                    " << outputFile << "\n"
           << "--------------------------------------------------\n" << std::endl;

    // Map
    gsFileData<> fd(mapFile);
    gsGeometry<>* map = NULL;
    if (fd.has< gsGeometry<> >())
    {
        map = fd.getFirst< gsGeometry<> >().release();
    }

    if (map == NULL)
    {
        gsInfo << "Didn't get the input geometry. Aborting..." << "\n";
        return 1;
    }

    gsMatrix<> pts;
    importMatrixFromASCII (dataFile, pts);
    mapPoints(pts.transpose(), *map, outputFile); // Check, if transponation is necessary

    return 0;
}

