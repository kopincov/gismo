/** @file gsTemplateMappingBarycentric.cpp

    @brief Template mapping (work in progress...)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh, S. Sajavicius
*/
// ./bin/gsTemplateMapping -l 1.0e-6 -i 7 -e 1.0e-4 -n 100 -T 2d_template_geometry.xml -G 2d_target_geometry.xml -S 2d_template_segmentation.xml -o 2d_res
// ./bin/gsTemplateMapping -l 1.0e-6 -i 7 -e 1.0e-4 -n 100 -T 3d_template_geometry.xml -G 3d_target_geometry.xml -S 3d_template_segmentation.xml -o 3d_res

#include <iostream>
#include <gismo.h>
//#include <gsIO/gsIOUtils.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

template <unsigned d>
int demo(const std::string& tempFile, const std::string& geomFile,
         const std::string& segmFile, const std::string& outputPrefix,
         const real_t lambda, const int numIterations,
         const real_t threshold, const int extension,
         const int degree, const int interiorKnots, const int N)
{

    gsMultiPatch<>* mp = NULL;

    mp = readMultipatch(tempFile);
    // Identifying id of rectangle
    int idxTemplateRectangle = rectangleIndex(tempFile);
    const gsGeometry<>& templateRectangle = (*mp)[idxTemplateRectangle];
    // Removing rectangle from multipatch
    const gsGeometry<>& templateInterior = (*mp)[1];
    //const gsMultiPatch<>* templateInterior = removePatch(*mp, idxTemplateRectangle);

    mp = readMultipatch(geomFile);
    // Identifying id of rectangle
    int idxTargetRectangle   = rectangleIndex(geomFile);
    gsInfo << "idxTargetRectangle: " << idxTargetRectangle << std::endl;
    //const gsGeometry<>& targetRectangle   = (*mp)[1];
    // Removing rectangle from multipatch
    const gsGeometry<>& targetInterior   = (*mp)[1];
    //const gsMultiPatch<>* targetInterior = removePatch(*mp, idxTargetRectangle);

    // SAMPLING POINTS
    // Detecting vertices of bounding rectangles
    gsMatrix<> templateAnchors = detectVertices(tempFile);
    gsMatrix<> targetAnchors   = detectVertices(geomFile);

    // Sampling bounding rectangle in template geometry
    gsMatrix<> uvRectangle = samplePoints<2>(templateRectangle, N);
    // Sampling circle in template geometry
    gsMatrix<> uvTarget    = samplePoints<2>(templateInterior, N);

    // Mapping rectangle points
    gsMatrix<> xyRectangle = mapPointsWachspress(templateAnchors, targetAnchors, uvRectangle);
    // Mapping target points
    gsMatrix<> xyTarget    = mapPointsWachspress(templateAnchors, targetAnchors, uvTarget);

    // Correcting mapped target points
    // Old version (to be changed):
    gsMatrix<> pars (1, xyTarget.cols());
    for (int i = 0; i < xyTarget.cols(); i++)
    {
         //gsInfo << "par: " << findParameter(inputGeometry, xy.col(i), 1.0e-12) << std::endl;
        pars.col(i) = findParameter(targetInterior, xyTarget.col(i), 1.0e-12);
    }
    xyTarget = targetInterior.eval(pars);
    // New version (to be done):
    /*
    */

    writeGeometry(targetInterior, "targetInterior");
    gsWriteParaviewPoints(uvRectangle, outputPrefix + "_uvRectangle");
    gsWriteParaviewPoints(uvTarget, outputPrefix + "_uvTarget");
    gsWriteParaviewPoints(xyRectangle, outputPrefix + "_xyRectangle");
    gsWriteParaviewPoints(xyTarget, outputPrefix + "_xyTarget");
    gsFileData<> fd_xyTarget;
    fd_xyTarget << xyTarget;
    fd_xyTarget.dump(outputPrefix + "_xyTarget");

    gsMatrix<> uv(uvRectangle.rows(), uvRectangle.cols()+uvTarget.cols());
    gsMatrix<> xy(xyRectangle.rows(), xyRectangle.cols()+xyTarget.cols());

    uv.block(0, 0, uvRectangle.rows(), uvRectangle.cols()) = uvRectangle;
    uv.block(0, uvRectangle.cols(), uvTarget.rows(), uvTarget.cols()) = uvTarget;

    xy.block(0, 0, xyRectangle.rows(), xyRectangle.cols()) = xyRectangle;
    xy.block(0, xyRectangle.cols(), xyTarget.rows(), xyTarget.cols()) = xyTarget;

    gsWriteParaviewPoints(uv, outputPrefix + "_Template_Geometry_Points");
    gsWriteParaviewPoints(xy, outputPrefix + "_Target_Geometry_Points");

    gsFileData<> fd;
    fd << uv;
    fd << xy;
    fd.dump(outputPrefix + "_Parameters_and_points");

    gsInfo << "Fitting..." << std::endl;
    gsMatrix<> bbox = boundingBox<d>(uv, 0.1);

    std::vector< gsKnotVector<> > KV;
    for(unsigned i = 0; i != d; i++)
    {
        KV.push_back(gsKnotVector<> ( bbox(i, 0), bbox(i, 1), interiorKnots, degree + 1 ));
    }

    gsTensorBSplineBasis<d> tbasis( KV );
    gsTHBSplineBasis<d> thb(tbasis);
    std::vector<unsigned> ext(d, extension);
    gsHFitting<d, real_t> fitting( uv, xy, thb, 1.0, ext, lambda);

    gsGeometry<> * thb_map_refined = NULL;
    for (int i = 0; i != numIterations; i++)
    {
        gsInfo << "Iteration: " << i + 1 << " / " << numIterations << std::endl;
        fitting.nextIteration(threshold, threshold);

        gsGeometry<> * thb_map = fitting.result();

        fitting.computeErrors();
        const std::vector<real_t>& errors = fitting.pointWiseErrors();
        real_t maxError = fitting.maxPointError();
        real_t percent = percentPointBelowTreshold(errors, threshold);

        gsInfo << "maxError: " << maxError << "\n"
               << "percent of data below threshold: " << percent << "\n";

        // saving
        writeMap(*thb_map, i, outputPrefix);
        writeKnots(*thb_map, i, outputPrefix);
        thb_map_refined = thb_map;

        if (maxError < threshold)
        {
            gsInfo << "Max error is below threshold: "
                   << maxError << " < " << threshold << std::endl;
            break;
        }

        // Giving different names for files with (possibly) final results
        writeMap(*thb_map, 0, outputPrefix + "_final");
        writeKnots(*thb_map, 0, outputPrefix + "_final");
    }

    mapMultipatch<d>(segmFile, *thb_map_refined, 100, outputPrefix+"_Segmentation_Mapped");

    return 0;
}


int main(int argc, char *argv[])
{
    // Options with default values
    std::string tempFile(MOTOR_DATA_DIR "jku/template_multipatch.xml");
    std::string geomFile(MOTOR_DATA_DIR "jku/target_multipatch.xml");
    std::string segmFile(MOTOR_DATA_DIR "jku/template_segmentation.xml");
    std::string outputPrefix("result");
    real_t lambda = 1e-6;
    int numIterations = 7;
    real_t threshold = 1e-4;
    int extension = 1;
    int degree = 3;
    int interiorKnots = 2;
    int N = 100;

    gsCmdLine cmd("Template mapping");
    cmd.addString("T", "templateFile", "Name of template geometry multipatch file (input)", tempFile);
    cmd.addString("G", "geometryFile", "Name of target geometry multipatch file (input)", geomFile);
    cmd.addString("S", "segmentationFile", "Name of template geometry segmentation multipatch file (input)", segmFile);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);
    cmd.addReal("l", "lambda", "Smoothing parameter", lambda);
    cmd.addInt("i", "numIter", "Number of refinement iterations.", numIterations);
    cmd.addReal("e", "error", "Error threshold", threshold);
    cmd.addInt("x", "extension", "The extension of the refinement", extension);
    cmd.addInt("d", "degree", "The degree of the spline", degree);
    cmd.addInt("k", "interKnots", "The number of interior knots", interiorKnots);
    cmd.addInt("n", "npts", "Number of sampling points (for fitting)", N);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments: \n"
           << "Template file:     " << tempFile << "\n"
           << "Geometry file:     " << geomFile << "\n"
           << "Segmentation file: " << segmFile << "\n"
           << "output:            " << outputPrefix << "\n"
           << "lambda:            " << lambda << "\n"
           << "numIterations:     " << numIterations << "\n"
           << "threshold:         " << threshold << "\n"
           << "extension:         " << extension << "\n"
           << "degree:            " << degree << "\n"
           << "interior knots:    " << interiorKnots << "\n"
           << "N:                 " << N << "\n"
           << "--------------------------------------------------" << std::endl;

/*
    if (tempFile == "" || geomFile == "")
    {
        gsInfo << "Wrong name of input file. Aborting..." << std::endl;
        return 0;
    }
*/

    gsMultiPatch<> * mp = readMultipatch(tempFile);
    const gsGeometry<>& geometry = (*mp)[0];
    if (geometry.parDim() == 1)
        return demo<2>(tempFile, geomFile,
                       segmFile, outputPrefix,
                       lambda, numIterations,
                       threshold, extension,
                       degree, interiorKnots, N);
    else if (geometry.parDim() == 2)
        return demo<3>(tempFile, geomFile,
                       segmFile, outputPrefix,
                       lambda, numIterations,
                       threshold, extension,
                       degree, interiorKnots, N);

}

