/** @file gsAirfoilTest.cpp

    @brief This file deals with airfoil geometries examples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/

#include <iostream>
#include <gismo.h>
#include <gsCore/gsMath.h>
#include <math.h>
#include "gsMotorUtils.h"
#include "gsAirfoilUtils.h"

using namespace gismo;

/******************************************************************/
int main(int argc, char *argv[])
{
    int mode = -1;

    std::string tempCase("In5patches");

    std::string inputFile("");
    std::string inputFile1("");
    std::string inputFile2("");
    std::string inputFile3("");
    std::string inputFile4("");

    std::string outputFile("");

    gsCmdLine cmd("Routines for working with airfoil demonstrator");
    cmd.addString("", "outputFile",  "", outputFile);
    cmd.addString("", "inputFile4",  "", inputFile4);
    cmd.addString("", "inputFile3",  "", inputFile3);
    cmd.addString("", "inputFile2",  "", inputFile2);
    cmd.addString("", "inputFile1",  "", inputFile1);
    cmd.addString("", "inputFile",   "", inputFile);
    cmd.addString("", "tempCase", "Template case", tempCase);
    cmd.addInt("m", "mode", "Mode:"
                            "  0 - Prepare file with template segmentation & splitting points and connectivity matrix\n"
                            "  1 - Prepare file with template segmentation lines\n"
                            "  2 - Construct and map planar Coons' patches", mode);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }


    // Template case (string code converted to int)
    const int tempCaseInt = tempCaseAsInt(tempCase);

    // Functionality depends on the value of "mode"
    switch(mode)
    {
        case 0: // Preparing file with template segmentation & splitting points and connectivity matrix
        {
            gsInfo << "* Preparing file with template segmentation & splitting points and connectivity matrix" << std::endl;
            gsMultiPatch<>* mpTemp = readMultipatch(inputFile);
            const gsMatrix<> splitPts = splittingPoints(*mpTemp, tempCaseInt);
            const gsMatrix<> segmPts  = segmentationPoints(splitPts, tempCaseInt);
            const index_t numSplitPts = splitPts.cols();
            const index_t numSegmPts  = segmPts.cols();

            gsMatrix<> pts(2, numSplitPts+numSegmPts);
            pts.setZero();
            pts.block(0, 0,           2, numSplitPts) = splitPts;
            pts.block(0, numSplitPts, 2, numSegmPts)  = segmPts;

            gsMatrix<> connectMatrix = connectivityMatrix(tempCaseInt);
            gsFileData<> fd;
            fd << pts;
            fd << connectMatrix;
            fd.dump(outputFile);
            break;
        }
        case 1: // Preparing file with template segmentation lines
        {
            gsInfo << "* Preparing file with template segmentation lines" << std::endl;
            gsFileData<> fd(inputFile);
            gsMatrix<> ptsMatrix =*fd.getId< gsMatrix<> >(0);
            gsMatrix<> connectMatrix =*fd.getId< gsMatrix<> >(1);

            gsMultiPatch<>* mpSegm = segmentationLines(ptsMatrix, connectMatrix);
            writeMultipatch(*mpSegm, outputFile);
            break;
        }
        case 2: // Constructing and mapping planar Coons' patches
        {
            gsInfo << "* Constructing and mapping planar Coons' patches" << std::endl;
            gsMultiPatch<>* mpTarg = readMultipatch(inputFile1);
            gsMultiPatch<>* mpTemp = readMultipatch(inputFile2);
            gsMultiPatch<>* mpSegm = readMultipatch(inputFile3);
            gsMultiPatch<>* mpTempSplitted = splittedTemplate(*mpTemp, tempCaseInt);
            gsMultiPatch<>* mpTargSplitted = splittedTemplate(*mpTarg, tempCaseInt);
            writeMultipatch(*mpTargSplitted, "target_multipatch_splitted");
            writeMultipatch(*mpTempSplitted, "template_multipatch_splitted");
            std::vector< gsMultiPatch<> > tempBound = templateBoundaries(*mpTempSplitted, *mpSegm, tempCaseInt);
            std::vector< gsTensorBSpline<2, real_t> > cpatches = getCoonsPatchesPlanar(tempBound);
            gsGeometry<>* map = readGeometry(inputFile4);
            mapCoonsPatches<2>(map, cpatches, 1000);
            break;
        }
        default:
        {
            gsInfo << "Invalid mode value" << std::endl;
            return 0;
        }

    };

    return 0;
}


