/** @file gsMapMultipatch.cpp

    @brief This program applies a mapping for a multipatch given in an xml file.

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
    std::string multipatchFile(MOTOR_DATA_DIR "jku/template_segmentation.xml");
    std::string mapFile(MOTOR_DATA_DIR "jku/thb_map.xml");
    std::string outputFile("target_segmentation_optimized");
    int npts  = 100;

    gsCmdLine cmd("This program applies a mapping for a multipatch given in an xml file");
    cmd.addString("f", "multipatchFile", "Name of the input (multipatch) file", multipatchFile);
    cmd.addString("m", "mapFile", "Name of the input (map) file", mapFile);
    cmd.addInt("s", "npts", "Number of points used for sampling each patch", npts);
    cmd.addString("o", "output", "The name of output file (prefix for all output files)", outputFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input (data) file:              " << multipatchFile << "\n"
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

    gsMultiPatch<> * mp = readMultipatch(multipatchFile);
    const gsGeometry<>& geometry = (*mp)[0];
    if (geometry.parDim() == 1)
    {
        mapMultipatch<2>(multipatchFile, *map, npts, outputFile);
    }
    else if (geometry.parDim() == 2)
    {
        mapMultipatch<3>(multipatchFile, *map, npts, outputFile);
    }

    return 0;
}

