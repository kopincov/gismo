/** @file gsSRM46RotorTest.cpp

    @brief This file deals with the SRM4+6 Rotor geometry

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
#include "gsSRM46RotorUtils.h"

using namespace gismo;


/******************************************************************/

int main(int argc, char *argv[])
{
    int mode = -1;
    // Distance from the center of the template to the rectangle
    real_t dist = 15.0;
    // Lengths of the rectangle sides
    real_t l1 = 15.0;
    real_t l2 = 15.0;
    // Parameter used to identify the lobe centers
    real_t parDelta = -0.01;
    // Parameters used to control the deviation of segmentation lines
    // joining the rectangle vertices and template boundary
    real_t parDelta1 = 0.02;
    real_t parDelta2 = 0.04;

    int numPoints = 1000;

    gsCmdLine cmd("Routines for working with ... demonstrator");
    cmd.addInt("m", "mode", "Mode:\n"
                            "  0 - prepare file with template segmentation points\n"
                            "  1 - ???\n"
                            "  9 - testing mode", mode);
    cmd.addReal("D", "dist", "Distance from the center of the template to the rectangle", dist);
    cmd.addReal("l", "l1", "Length of the rectangle side (1)", l1);
    cmd.addReal("L", "l2", "Length of the rectangle side (2)", l2);
    cmd.addReal("", "parDelta", "Parameter used to identify the lobe centers", parDelta);
    cmd.addReal("", "parDelta1", "Parameter used to control the deviation of segmentation lines (1)", parDelta1);
    cmd.addReal("", "parDelta2", "Parameter used to control the deviation of segmentation lines (2)", parDelta2);
    cmd.addInt("n", "numPoints", "Number of points", numPoints);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    switch(mode)
    {
        case 0: // Preparing file with template segmentation points (Female rotor, Template 1)
        {
            gsInfo << "* Preparing file with template segmentation points (Female rotor, Template 1)" << std::endl;
            const std::string inputFile("SRM4+6_FemaleRotor_temp1.xml");
            const std::string outputFile("SRM4+6_FemaleRotor_temp1_segm_pts.xml");

            const gsGeometry<>* boundary = readGeometry(inputFile);
            const gsMatrix<> c = FemaleRotCenter(boundary, parDelta);
            const gsMatrix<> params = FemaleRotParams(parDelta);

            gsMatrix<> rectanglePts(2, 24);
            const gsMatrix<> ptsBound = boundary->eval(params);
            for (int i = 0; i < 6; i++)
            {
                gsMatrix <> ptsRect = FemaleRotTemp1RectanglePoints(c, ptsBound.col(i), dist, l1, l2);
                rectanglePts.block(0, 4*i, 2, 4) = ptsRect;
            }


            gsMatrix<> boundaryPts = FemaleRotTemp1BoundaryPoints(boundary, parDelta, parDelta1, parDelta2);

            gsMatrix<> segmPts(2, 48);
            segmPts.block(0, 0, 2, 24) = rectanglePts;
            segmPts.block(0, 24, 2, 24) = boundaryPts;

            gsFileData<> fd;
            fd << segmPts;
            fd << FemaleRotTemp1ConnectivityMatrix();
            fd.dump(outputFile);
            break;
        }
        case 1: // Splitting boundaries of planar template and target geometries (Female rotor, Template 1)
        {
            gsInfo << "* Splitting boundaries of planar template and target geometries (Female rotor, Template 1)" << std::endl;
            const std::string inputTempFile("SRM4+6_FemaleRotor_temp1.xml");
            const std::string outputTempFile("SRM4+6_FemaleRotor_temp1_splitted");

            const std::string inputTargFile("SRM4+6_FemaleRotor.xml");
            const std::string outputTargFile("SRM4+6_FemaleRotor_splitted");

            const gsMatrix<> params = FemaleRotTemp1BoundaryParameters(parDelta, parDelta1, parDelta2);

            const gsMultiPatch<>* mpTempSplitted = FemaleRotSplittedBoundary(inputTempFile, outputTempFile, params);
            const gsMultiPatch<>* mpTargSplitted = FemaleRotSplittedBoundary(inputTargFile, outputTargFile, params);

            /*
            const gsMultiPatch<>* mpTemp = joinMultipatches("", "");
            writeMultipatch(*mpTempSplitted, outputFile);
            */
            //gsInfo << "mpTempSplitted:\n" << *mpTempSplitted << std::endl;
            //writeMultipatch(*mpTempSplitted, outputFile);
            break;
        }
        case 2: // Constructing and mapping planar Coons' patches (Female rotor, Template 1)
        {
            gsInfo << "* Constructing planar Coon's patches (Female rotor, Template 1)" << std::endl;
            gsInfo << "  ... Reading multipatches" << std::endl;
            gsMultiPatch<>* mpTemp = readMultipatch("SRM4+6_FemaleRotor_temp1_splitted.xml");
            gsMultiPatch<>* mpSegm = readMultipatch("SRM4+6_FemaleRotor_temp1_segm.xml");
            gsInfo << "  ... FemaleRotTemp1Boundaries" << std::endl;
            std::vector< gsMultiPatch<> > boundaries = FemaleRotTemp1Boundaries(*mpTemp, *mpSegm);
            gsInfo << "  ... getCoonsPatchesPlanar" << std::endl;
            std::vector< gsTensorBSpline<2, real_t> > cpatches = getCoonsPatchesPlanar(boundaries);
            gsInfo << "* Mapping planar Coon's patches" << std::endl;
            gsGeometry<>* map = readGeometry("planar_param_final_THBMap_0.xml");
            mapCoonsPatches<2>(map, cpatches, numPoints);
            break;
        }
        case 3: // Preparing file with template geometry (Male rotor, Template 1 - rectangle)
        {
            gsInfo << "* Preparing file with template geometry (Male rotor, Template 1 - rectangle)" << std::endl;
            const std::string outputFile("SRM4+6_MaleRotor_temp1");

            //gsBSpline<> boundary = MaleRotTemp1Boundary();
            //writeGeometry(boundary, outputFile);
            break;
        }
        case 4: // Preparing file with template segmentation points (Male rotor, Template 1)
        {
            gsInfo << "* Preparing file with template segmentation points (Male rotor, Template 1)" << std::endl;
            const std::string inputFile("SRM4+6_MaleRotor_temp1.xml");
            const std::string outputFile("SRM4+6_MaleRotor_temp1_segm_pts");

            const gsGeometry<>* boundary = readGeometry(inputFile);
            //const gsMatrix<> params = FemaleRotTemp1BoundaryParameters();

            gsMatrix<> rectanglePts(4, 2);
            rectanglePts << 15, 0,
                            0, 15,
                            -15, 0,
                            0, -15;

            gsMatrix<> boundaryPts = MaleRotTemp1BoundaryPoints(boundary);
            //gsWriteParaviewPoints(boundaryPts, "test_boundaryPts");

            gsMatrix<> segmPts(2, 8) ;
            gsMatrix<> rectanglePtsT = rectanglePts.transpose();
            segmPts.block(0, 0, 2, 4) = rectanglePtsT;
            segmPts.block(0, 4, 2, 4) = boundaryPts;
            gsFileData<> fd;
            fd << segmPts;
            fd << MaleRotTemp1ConnectivityMatrix();
            fd.dump(outputFile);
            break;
        }
        case 5: // Splitting boundaries of planar template geometry (Male rotor, Template 1)
        {
            gsInfo << "* Splitting boundaries of planar template geometry (Male rotor, Template 1)" << std::endl;
            const std::string inputTempFile("SRM4+6_MaleRotor_temp1.xml");
            const std::string outputTempFile("SRM4+6_MaleRotor_temp1_splitted");
            const std::string inputTargFile("SRM4+6_MaleRotor.xml");
            const std::string outputTargFile("SRM4+6_MaleRotor_splitted");

            const gsMatrix<> params = MaleRotTemp1BoundaryParameters();

            const gsMultiPatch<>* mpTempSplitted = MaleRotTemp1SplittedBoundary(inputTempFile, outputTempFile, params);
            const gsMultiPatch<>* mpTargSplitted = MaleRotTemp1SplittedBoundary(inputTargFile, outputTargFile, params);

            break;
        }
        case 6: // Constructing and mapping planar Coons' patches (Male rotor, Template 1)
        {
            gsInfo << "* Constructing planar Coon's patches (Male rotor, Template 1)" << std::endl;
            gsInfo << "  ... Reading multipatches" << std::endl;
            gsMultiPatch<>* mpTemp = readMultipatch("SRM4+6_MaleRotor_temp1_splitted.xml");
            gsMultiPatch<>* mpSegm = readMultipatch("SRM4+6_MaleRotor_temp1_segm.xml");
            gsInfo << "  ... MaleRotTemp1Boundaries" << std::endl;
            std::vector< gsMultiPatch<> > boundaries = MaleRotTemp1Boundaries(*mpTemp, *mpSegm);

            gsInfo << "  ... getCoonsPatchesPlanar" << std::endl;
            std::vector< gsTensorBSpline<2, real_t> > cpatches = getCoonsPatchesPlanar(boundaries);
            gsInfo << "* Mapping planar Coon's patches" << std::endl;
            gsGeometry<>* map = readGeometry("planar_param_final_THBMap_0.xml");
            mapCoonsPatches<2>(map, cpatches, numPoints);
            break;
        }
        default:
            gsInfo << "Invalid value of mode (should be integer between 0 and ...)" << std::endl;
            return 0;
    };

    return 0;
}
