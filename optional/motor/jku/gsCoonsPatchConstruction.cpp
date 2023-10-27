/** @file gsCoonsPatchConstruction.cpp

    @brief ...

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
    int tempCase = 1;

    std::string tempMultipatchFile(MOTOR_DATA_DIR "jku/template_multipatch.xml");
    std::string tempSegmFile(MOTOR_DATA_DIR "jku/template_segmentation.xml");
    std::string tempSegmPointsFile(MOTOR_DATA_DIR "jku/template_segmentation_points.xml");
    std::string mapFile(MOTOR_DATA_DIR "jku/thb_map.xml");
    std::string outputPrefix("result");

    gsCmdLine cmd("Coons` patch construction");
    cmd.addString("M", "tempMultipatch", "Name of file with template multipatch", tempMultipatchFile);
    cmd.addString("S", "tempSegm",       "Name of template segmentation file", tempSegmFile);
    cmd.addString("P", "tempSegmPoints", "Name of template segmentation points file", tempSegmPointsFile);
    cmd.addInt("t", "tempCase", "Templete case (expected values: 1 or 2)", tempCase);
    cmd.addString("F", "mapFile", "Name of map file", mapFile);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if ( (tempCase != 1) && (tempCase != 2) )
    {
        gsInfo << "Unknown template case.\n";
        return 1;
    }

    gsInfo << "Input Arguments: \n"
           << "Template multipatch:           " << tempMultipatchFile << "\n"
           << "Template segmentation:         " << tempSegmFile << "\n"
           << "Template segmentation points:  " << tempSegmPointsFile << "\n"
           << "Template case:                 " << tempCase << "\n"
           << "Map:                           " << mapFile << "\n"
           << "output:                        " << outputPrefix << "\n"
           << "---------------------------------------------------------\n"
           << std::endl;

    gsFileData<> fd(tempSegmPointsFile);
    gsMatrix<>* points = fd.getId< gsMatrix<> >(0).release();

    gsMultiPatch<>* mp = readMultipatch(tempMultipatchFile);
    gsMultiPatch<>* mpSegmentation = readMultipatch(tempSegmFile);
    std::vector< gsMultiPatch<> > boundaries;


    if ( tempCase == 1 )
    {
        // Splitting the template of the type 1
        std::vector< gsMultiPatch<> > splittedTemp = splittedTemplate1(*mp, *points);
        gsMultiPatch<> mpRectangle = splittedTemp.at(0);
        gsMultiPatch<> mpEllipsis = splittedTemp.at(1);
        boundaries = templateSkeletonBoundaries1(mpRectangle, mpEllipsis, *mpSegmentation);
    }
    else if ( tempCase == 2 )
    {
        // Splitting the template of the type 2
        std::vector< gsMultiPatch<> > splittedTemp = splittedTemplate2(*mp, *points);
        gsMultiPatch<> mpEllipsis = splittedTemp.at(0);
        boundaries = templateSkeletonBoundaries2(mpEllipsis, *mpSegmentation);
    }
    boundaries = unifyKnotsBoundaries(boundaries);
    std::vector< gsTensorBSpline<2, real_t> > cpatches = coonsPatches<2>(boundaries);

    // Making multipatch from Coons` patches
    // Making file with parameters and point
    gsMultiPatch<> mpCoonsPatches;
    gsMatrix<> parameters = uniformParameters<3>(0.0, 1.0, 100);
    gsFileData<> fdPointsParameters;

    gsGeometry<>* thb_map = readGeometry(mapFile);

    for (std::size_t i = 0; i != cpatches.size(); i++)
    {
        gsTensorBSpline<2, real_t> cp = cpatches.at(i);
        gsMatrix<> cpValues = cp.eval(parameters);
        gsMatrix<> cpMappedValues = thb_map->eval(cpValues);

        mpCoonsPatches.addPatch(cp);
        fdPointsParameters << parameters;
        fdPointsParameters << cpMappedValues;
    }
    mpCoonsPatches.computeTopology(1e-9, true);
    writeMultipatch(mpCoonsPatches, outputPrefix + "_multipatch");
    fdPointsParameters.dump(outputPrefix + "_points_and_parameters");

    delete thb_map;
    return 0;
}

