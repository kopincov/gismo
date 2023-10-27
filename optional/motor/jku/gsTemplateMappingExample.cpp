/** @file gsTemplateMappingExample.cpp

    @brief This file is used to generate example of simple THB-splines map.
    The map is later used to test TangentDistanceMinimization method
    (see gsTangentDistanceMinimizationTest.cpp file).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh, S. Sajavicius
*/

#include <iostream>
#include <gismo.h>
//#include <gsIO/gsIOUtils.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

template <unsigned d>
int demo(const std::string& tempFile, const std::string& geomFile,
         const std::string& outputPrefix,
         const real_t lambda, const int degree,
         const int interiorKnots, const int N,
         const bool debug)
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

    if (debug)
    {
        writeGeometry(targetInterior, outputPrefix + "targetInterior");
        gsWriteParaviewPoints(uvRectangle, outputPrefix + "_uvRectangle");
        gsWriteParaviewPoints(uvTarget, outputPrefix + "_uvTarget");
        gsWriteParaviewPoints(xyRectangle, outputPrefix + "_xyRectangle");
        gsWriteParaviewPoints(xyTarget, outputPrefix + "_xyTarget");
        gsFileData<> fd_xyTarget;
        fd_xyTarget << xyTarget;
        fd_xyTarget.dump(outputPrefix + "_xyTarget");
    }

    gsMatrix<> uv(uvRectangle.rows(), uvRectangle.cols()+uvTarget.cols());
    gsMatrix<> xy(xyRectangle.rows(), xyRectangle.cols()+xyTarget.cols());

    uv.block(0, 0, uvRectangle.rows(), uvRectangle.cols()) = uvRectangle;
    uv.block(0, uvRectangle.cols(), uvTarget.rows(), uvTarget.cols()) = uvTarget;

    xy.block(0, 0, xyRectangle.rows(), xyRectangle.cols()) = xyRectangle;
    xy.block(0, xyRectangle.cols(), xyTarget.rows(), xyTarget.cols()) = xyTarget;

    if (debug)
    {
        gsWriteParaviewPoints(uv, outputPrefix + "_Template_Geometry_Points");
        gsWriteParaviewPoints(xy, outputPrefix + "_Target_Geometry_Points");
    }

    gsInfo << "Fitting..." << std::endl;
    gsMatrix<> bbox = boundingBox<d>(uv, 0.1);
    std::vector< gsKnotVector<> > KV;
    for(unsigned i = 0; i != d; i++)
    {
        KV.push_back(gsKnotVector<> ( bbox(i, 0), bbox(i, 1), interiorKnots, degree + 1 ));
    }

    gsTensorBSplineBasis<d> tbasis( KV );
    gsTHBSplineBasis<d> thb(tbasis);

    gsFitting<> fitting( uv, xy, thb);

    fitting.compute(lambda);
    gsGeometry<> * thb_map = fitting.result();

    fitting.computeErrors();

    real_t maxError = fitting.maxPointError();
    gsInfo << "maxError: " << maxError << std::endl;

    // saving
    writeMap(*thb_map, outputPrefix);
    writeKnots(*thb_map, outputPrefix);

    return 0;
}

int main(int argc, char *argv[])
{
    // Options with default values
    std::string tempFile(MOTOR_DATA_DIR "jku/template_multipatch.xml");
    std::string geomFile(MOTOR_DATA_DIR "jku/target_multipatch.xml");
    std::string outputPrefix("result");
    real_t lambda = 1e-6;
    int degree = 3;
    int interiorKnots = 2;
    int N = 100;
    bool debug = false;

    gsCmdLine cmd("Template mapping");
    cmd.addString("T", "templateFile", "Name of template geometry multipatch file (input)", tempFile);
    cmd.addString("G", "geometryFile", "Name of target geometry multipatch file (input)", geomFile);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);
    cmd.addReal("l", "lambda", "Smoothing parameter", lambda);
    cmd.addInt("d", "degree", "The degree of the spline", degree);
    cmd.addInt("k", "interKnots", "The number of interior knots", interiorKnots);
    cmd.addInt("n", "npts", "Number of sampling points (for fitting)", N);
    cmd.addSwitch("debug", "If debug, it outputs more data.", debug);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments: \n"
           << "Template file:     " << tempFile << "\n"
           << "Geometry file:     " << geomFile << "\n"
           << "lambda:            " << lambda << "\n"
           << "degree:            " << degree << "\n"
           << "interior knots:    " << interiorKnots << "\n"
           << "N:                 " << N << "\n"
           << "Debug:             " << debug << "\n"
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
                       outputPrefix,
                       lambda, degree,
                       interiorKnots, N, debug);
    else if (geometry.parDim() == 2)
        return demo<3>(tempFile, geomFile,
                       outputPrefix,
                       lambda, degree,
                       interiorKnots, N, debug);

}

