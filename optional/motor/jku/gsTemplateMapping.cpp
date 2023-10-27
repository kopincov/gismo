/** @file gsTemplateMapping.cpp

    @brief Constructs a map which transforms (maps) the template shape
    to the target shape. The map is constructed using hierarchical fitting
    of parametrized point clouds (see gsHSplines/gsHFitting)

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
    gsMatrix<> uv;
    gsMatrix<> xy;

    gsFileData<> fdTemp(tempFile);
    gsFileData<> fdGeom(geomFile);

    gsInfo << "Sampling points..." << std::endl;
    uv = samplePoints<d>(tempFile, N);
    xy = samplePoints<d>(geomFile, N);

    // Writing points to Paraview format
    gsWriteParaviewPoints(uv, outputPrefix + "_Template_Geometry_Points");
    gsWriteParaviewPoints(xy, outputPrefix + "_Target_Geometry_Points");

    gsFileData<> fd;
    fd << uv;
    fd << xy;
    fd.dump(outputPrefix + "_Parameters_and_points");

    gsInfo << "Fitting..." << std::endl;
    gsMatrix<> bbox = boundingBox<d>(uv, 0.1); // Standard value: 0.1
    /*
    gsMatrix<> bbox(3, 2);
    bbox(0, 0) = -0.5; bbox(0, 1) = 1.5;
    bbox(1, 0) = -0.5; bbox(1, 1) = 1.5;
    bbox(2, 0) = -0.5; bbox(2, 1) = 1.5;
    */

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

        // Checking the Jacobian
        //checkJacobianDeterminant(*thb_map_refined, 100000, false, "jac_pts", i);

        // Giving different names for files with (possibly) final results
        writeMap(*thb_map, 0, outputPrefix + "_final");
        writeKnots(*thb_map, 0, outputPrefix + "_final");
    }

    // Template points mapping
    mapPoints(uv, *thb_map_refined, outputPrefix+"_Points_Mapped");
    // Template geometry mapping
    mapMultipatch<d>(tempFile, *thb_map_refined, 100, outputPrefix+"_Template_Geometry_Mapped");
    // Segmentation points mapping
    gsMatrix<> uvSegm = samplePoints<d>(segmFile, N);
    mapPoints(uvSegm, *thb_map_refined, outputPrefix+"_Segmentation_Points_Mapped");
    // Segmentation lines/surfaces mapping (if segmentation lines are given)
    if (segmFile != ""){
        mapMultipatch<d>(segmFile, *thb_map_refined, 100, outputPrefix+"_Segmentation_Mapped");
    }

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

    gsCmdLine cmd("Constructs a map which transforms (maps) the template shape to the target shape. The map is constructed using hierarchical fitting of parametrized point clouds (see gsHSplines/gsHFitting)");
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);
    cmd.addReal("t", "threshold", "Error threshold", threshold);
    cmd.addInt("x", "extension", "Extension of the refinement", extension);
    cmd.addInt("r", "numRefIter", "Number of adaptive (local) refinement procedure iterations", numIterations);
    cmd.addReal("l", "lambda", "Smoothing parameter", lambda);
    cmd.addInt("n", "numPts", "Number of sampling points (for fitting)", N);
    cmd.addInt("k", "numKnots", "Number of interior knots", interiorKnots);
    cmd.addInt("d", "degree", "Degree", degree);
    cmd.addString("S", "segmentationFile", "File containing template domain skeleton (segmentation) multi-patch", segmFile);
    cmd.addString("G", "geometryFile", "File containing target shape multi-patch", geomFile);
    cmd.addString("T", "templateFile", "File containing template shape multi-patch or point samples on the template shape", tempFile);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments: \n"
           << "Template file:     " << tempFile << "\n"
           << "Geometry file:     " << geomFile << "\n"
           << "Segmentation file: " << segmFile << "\n"
           << "output:            " << outputPrefix << "\n"
           << "lambda:            " << lambda << "\n"
           << "numRefIter:        " << numIterations << "\n"
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

