/** @file gsCrescendoUtils.h

    @brief This file contains various auxilary routines used with
    crescendo-tmtf demonstrator geometry from the MOTOR project
    (http://motor-project.eu)

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Sajavicius
*/

#pragma once

//#include "gsMotorIOUtils.h"

namespace gismo
{


// Declarations

// Takes crescendo-tmtf geometry (template OR target),
// extracts Patch 4 and does some preparation for 2-D parametrization:
// - identifies Patch 4 boundaries
// - save boundaries in gsMultiPatch format
gsMultiPatch<>* getCrescendoPatch4Boundaries(const std::string& inputFile);

// Constructs ...
gsMultiPatch<>* getCrescendo25D(const std::string& inputFile,
                                const real_t height);

// Takes crescendo-tmtf-para-mp Patch 4 boundaries
// (prepared with getCrescendoPatch4Boundaries)
// and splits them at certain points determined by segmentation lines
gsMultiPatch<>* CrescendoPatch4PlanarSplittedBoundary(const std::string& inputFile);

// ...
gsMultiPatch<>* splitCrescendoParaPatch4VolumetricBoundaries(const std::string& inputParaFile,
                                                             const std::string& inputSegmFile);
// Returns ...
std::vector< gsMultiPatch<> > getCrescendoParaSkeletonPlanarBoundaries(const gsMultiPatch<>& mpTemp,
                                                                       const gsMultiPatch<>& mpSegm);
// Returns ...
std::vector< gsMultiPatch<> > getCrescendoParaSkeletonVolumetricBoundaries(const gsMultiPatch<>& mpTemp,
                                                                           const gsMultiPatch<>& mpSegm);
gsMultiPatch<>* getCrescendoPlanarTemplate(real_t x1, real_t x2,
                                           real_t ctrlPtTemp1x, real_t ctrlPtTemp1y,
                                           real_t ctrlPtTemp2x, real_t ctrlPtTemp2y);
gsMultiPatch<>* getCrescendoPlanarSegmentation(real_t x1, real_t x2,
                                               real_t ctrlPtSegm1x, real_t ctrlPtSegm1y,
                                               real_t ctrlPtSegm2x, real_t ctrlPtSegm2y);
// Constructs a multipatch for splitted bottom and top parts of
// the volumetric template geometry
gsMultiPatch<>* getCrescendoParaSplittedTopBottom(const std::vector< gsTensorBSpline<2, real_t> > cpatches);
//
gsMatrix<> boundParameters(const int numPoints, const int id);
//
void getCrescendoPlanarBoundaryConstraints(const std::string& inputBoundaryFile,
                                           const std::string& inputCoonsFile,
                                           const std::string& inputParamsPointsFile,
                                           const std::string& outputFile);
/**************************************************/
gsMultiPatch<>* getCrescendoPatch4Boundaries(const std::string& inputFile)
{
    // Reading multipatch:
    // The template geometry
    const gsMultiPatch<> * mp = readMultipatch(inputFile);

    // Patch 4 from the geometry
    const gsGeometry<>& patch4 = (*mp)[4];
    const gsTensorBSpline<2, real_t> geom4 = dynamic_cast< const gsTensorBSpline<2, real_t>& >(patch4);

    // Control points of the geometry
    gsMatrix<> ctrlPts = geom4.coefs();

    // Initialize control point matrices
    gsMatrix<> ctrlPtsBound0(39, 2);
    gsMatrix<> ctrlPtsBound1(39, 2);
    gsMatrix<> ctrlPtsBound2(59, 2);
    gsMatrix<> ctrlPtsBound3(59, 2);
    // Filling the control points matrices
    for (int i = 0; i < 39; i++)
    {
        ctrlPtsBound0.block(i, 0, 1, 2) = ctrlPts.block(58 + 59*i, 0, 1, 2);
        ctrlPtsBound1.block(i, 0, 1, 2) = ctrlPts.block(0  + 59*i, 0, 1, 2);
    }
    for (int i = 0; i < 59; i++)
    {
        ctrlPtsBound2.block(i, 0, 1, 2) = ctrlPts.block(2242 + i, 0, 1, 2);
        ctrlPtsBound3.block(i, 0, 1, 2) = ctrlPts.block(0 + i, 0, 1, 2);
    }

    // Template geometry boundaries
    const gsGeometry<>& bound0 = *geom4.boundary(0);
    const gsGeometry<>& bound1 = *geom4.boundary(1);
    const gsGeometry<>& bound2 = *geom4.boundary(2);
    const gsGeometry<>& bound3 = *geom4.boundary(3);

    const gsBSpline<> geomBound0 = dynamic_cast< const gsBSpline<>& >(bound0);
    const gsBSpline<> geomBound1 = dynamic_cast< const gsBSpline<>& >(bound1);
    const gsBSpline<> geomBound2 = dynamic_cast< const gsBSpline<>& >(bound2);
    const gsBSpline<> geomBound3 = dynamic_cast< const gsBSpline<>& >(bound3);

    const gsBSpline<> paraBound0(geomBound0.knots(), ctrlPtsBound0);
    const gsBSpline<> paraBound1(geomBound1.knots(), ctrlPtsBound1);
    const gsBSpline<> paraBound2(geomBound3.knots(), ctrlPtsBound2); // Don't worry this is not bug - this should be like this
    const gsBSpline<> paraBound3(geomBound3.knots(), ctrlPtsBound3); // Constructing and mapping volumetric Coons' patchesknots(), ctrlPtsBound3);

    // Multi-patch for boundaries
    gsMultiPatch<>* mpBound = new gsMultiPatch<>();
    mpBound->addPatch(paraBound0);
    mpBound->addPatch(paraBound1);
    mpBound->addPatch(paraBound2);
    mpBound->addPatch(paraBound3);

    return mpBound;
}


/**************************************************/
gsMultiPatch<>* getCrescendo25D(const std::string& inputFile,
                                const real_t height)
{
    const gsMultiPatch<>* mp = readMultipatch(inputFile);
    const gsGeometry<>& patch0 = (*mp)[0];
    const gsGeometry<>& patch1 = (*mp)[1];
    const gsGeometry<>& patch2 = (*mp)[2];
    const gsGeometry<>& patch3 = (*mp)[3];
    const gsBSpline<> geom0 = dynamic_cast< const gsBSpline<>& >(patch0);
    const gsBSpline<> geom1 = dynamic_cast< const gsBSpline<>& >(patch1);
    const gsBSpline<> geom2 = dynamic_cast< const gsBSpline<>& >(patch2);
    const gsBSpline<> geom3 = dynamic_cast< const gsBSpline<>& >(patch3);

    // Extruding the boundary
    gsMultiPatch<>* mpVol = extrudeMultipatch(*mp, height);

    // Bottom
    gsMultiPatch<>* boundary = new gsMultiPatch<>();
    boundary->addPatch(geom1);
    boundary->addPatch(geom0);
    boundary->addPatch(geom3);
    boundary->addPatch(geom2);

    gsCoonsPatch<real_t> cpatch = coonsPatch(*boundary);
    cpatch.compute();
    const gsTensorBSpline<2>& cpatchInit = dynamic_cast< const gsTensorBSpline<2>& >(cpatch.result());

    // Control points of the bottom
    const gsMatrix<> coefsInit =cpatchInit.coefs();
    const index_t nCoefs = coefsInit.rows();

    gsMatrix<> coefs0(nCoefs, 3);
    coefs0.setZero();
    coefs0.block(0, 0, nCoefs, 2) = coefsInit.block(0, 0, nCoefs, 2);
    gsKnotVector<> kv1 = cpatchInit.knots(0);
    gsKnotVector<> kv2 = cpatchInit.knots(1);
    gsTensorBSpline<2> cp0(kv1, kv2, coefs0);
    mpVol->addPatch(cp0);

    // Top
    gsMatrix<> coefs1(nCoefs, 3);
    coefs1.block(0, 0, nCoefs, 3) = coefs0.block(0, 0, nCoefs, 3);
    for (int i = 0; i < nCoefs; i++)
    {
        coefs1(i, 2) = height;
    }
    gsTensorBSpline<2> cp1(kv1, kv2, coefs1);
    mpVol->addPatch(cp1);

    return mpVol;
}

/**************************************************/
gsMultiPatch<>* CrescendoPatch4PlanarSplittedBoundary(const std::string& inputFile)
{

    // Reading multipatch:
    const gsMultiPatch<>* mpPara = readMultipatch(inputFile);

    // Template geometry boundaries
    const gsGeometry<>& paraPatch0 = (*mpPara)[0];
    const gsGeometry<>& paraPatch1 = (*mpPara)[1];
    const gsGeometry<>& paraPatch2 = (*mpPara)[2];
    const gsGeometry<>& paraPatch3 = (*mpPara)[3];

    const gsBSpline<> paraGeom0 = dynamic_cast< const gsBSpline<>& >(paraPatch0);
    const gsBSpline<> paraGeom1 = dynamic_cast< const gsBSpline<>& >(paraPatch1);
    const gsBSpline<> paraGeom2 = dynamic_cast< const gsBSpline<>& >(paraPatch2);
    const gsBSpline<> paraGeom3 = dynamic_cast< const gsBSpline<>& >(paraPatch3);

    // Parameters corresponding to the template splitting points
    gsMatrix<> paramsParaGeom2(1, 2);
    gsMatrix<> paramsParaGeom3(1, 2);
    paramsParaGeom2(0, 0) = 0.3; paramsParaGeom2(0, 1) = 0.7;
    paramsParaGeom3(0, 0) = 0.3; paramsParaGeom3(0, 1) = 0.7;

    gsMultiPatch<>* mpParaGeom = new gsMultiPatch<>();

    // Right vertical side:
    mpParaGeom->addPatch(paraGeom0);
    // Left vertical side:
    mpParaGeom->addPatch(paraGeom1);
    // Splitted top horizontal side:
    gsMultiPatch<>* mpParaGeom2 = splitOpenBSpline(paraGeom2, paramsParaGeom2);
    mpParaGeom = joinMultipatches(*mpParaGeom, *mpParaGeom2);
    // Splitted bottom horizontal side:
    gsMultiPatch<>* mpParaGeom3 = splitOpenBSpline(paraGeom3, paramsParaGeom3);
    mpParaGeom = joinMultipatches(*mpParaGeom, *mpParaGeom3);

    return mpParaGeom;
}

/**************************************************/
gsMultiPatch<>* splitCrescendoParaPatch4VolumetricBoundaries(const std::string& inputParaFile,
                                                             const std::string& inputSegmFile)
{
    const gsMultiPatch<>* mpSegm = readMultipatch(inputSegmFile);

    const gsMultiPatch<>* mpParaPlanarSplitted = CrescendoPatch4PlanarSplittedBoundary(inputParaFile);
    gsMultiPatch<>* mpParaVol = extrudeMultipatch(*mpParaPlanarSplitted, 1.0);

    const std::vector< gsMultiPatch<> > boundaries = getCrescendoParaSkeletonPlanarBoundaries(*mpParaPlanarSplitted, *mpSegm);
    const std::vector< gsTensorBSpline<2, real_t> > cpatches = getCoonsPatchesPlanar(boundaries);
    const gsMultiPatch<>* crescendoTopBottom = getCrescendoParaSplittedTopBottom(cpatches);

    mpParaVol = joinMultipatches(*mpParaVol, *crescendoTopBottom);

    return mpParaVol;
}

/**************************************************/
std::vector< gsMultiPatch<> > getCrescendoParaSkeletonPlanarBoundaries(const gsMultiPatch<>& mpTemp,
                                                                       const gsMultiPatch<>& mpSegm)
{
    std::vector< gsMultiPatch<> > boundaries;
    gsMultiPatch<> boundary0;
    gsMultiPatch<> boundary1;
    gsMultiPatch<> boundary2;

    boundary0.addPatch(mpTemp[0]);
    boundary0.addPatch(mpSegm[1]);
    boundary0.addPatch(mpTemp[4]);
    boundary0.addPatch(mpTemp[7]);
    boundaries.push_back(boundary0);

    boundary1.addPatch(mpSegm[1]);
    boundary1.addPatch(mpSegm[0]);
    boundary1.addPatch(mpTemp[3]);
    boundary1.addPatch(mpTemp[6]);
    boundaries.push_back(boundary1);

    boundary2.addPatch(mpSegm[0]);
    boundary2.addPatch(mpTemp[1]);
    boundary2.addPatch(mpTemp[2]);
    boundary2.addPatch(mpTemp[5]);
    boundaries.push_back(boundary2);

    return boundaries;
}

/**************************************************/
std::vector< gsMultiPatch<> > getCrescendoParaSkeletonVolumetricBoundaries(const gsMultiPatch<>& mpTemp,
                                                                           const gsMultiPatch<>& mpSegm)
{
    std::vector< gsMultiPatch<> > boundaries;
    gsMultiPatch<> boundary0;
    gsMultiPatch<> boundary1;
    gsMultiPatch<> boundary2;

    boundary0.addPatch(mpTemp[0]);
    boundary0.addPatch(mpTemp[4]);
    boundary0.addPatch(mpSegm[1]);
    boundary0.addPatch(mpTemp[7]);
    boundary0.addPatch(mpTemp[8]);
    boundary0.addPatch(mpTemp[11]);
    boundaries.push_back(boundary0);

    boundary1.addPatch(mpSegm[1]);
    boundary1.addPatch(mpTemp[3]);
    boundary1.addPatch(mpSegm[0]);
    boundary1.addPatch(mpTemp[6]);
    boundary1.addPatch(mpTemp[9]);
    boundary1.addPatch(mpTemp[12]);
    boundaries.push_back(boundary1);

    boundary2.addPatch(mpSegm[0]);
    boundary2.addPatch(mpTemp[2]);
    boundary2.addPatch(mpTemp[1]);
    boundary2.addPatch(mpTemp[5]);
    boundary2.addPatch(mpTemp[10]);
    boundary2.addPatch(mpTemp[13]);
    boundaries.push_back(boundary2);

    return boundaries;
}

/**************************************************/
gsMultiPatch<>* getCrescendoPlanarTemplate(real_t x1, real_t x2,
                                           real_t ctrlPtTemp1x, real_t ctrlPtTemp1y,
                                           real_t ctrlPtTemp2x, real_t ctrlPtTemp2y)
{
    gsMultiPatch<>* mpTemp = new gsMultiPatch<>();

    //    +-----+-----+-----+
    //    |     \_____/     |
    //    |                 |
    //    |                 |
    //    |                 |
    //    |      _____      |
    //    |     /     \     |
    //    +-----+-----+-----+
    //          x1    x2

    // Knot vectors
    const gsKnotVector<> kv0 (0, 1, 0, 3); // {0, 0, 0, 1, 1, 1};
    const gsKnotVector<> kv1 (0, 1, 0, 3);
    const real_t tmp2[] = {0, 0, 0, .3, .3, .7, .7, 1, 1, 1};
    const std::vector<real_t> knots2(tmp2, tmp2 + sizeof(tmp2) / sizeof(tmp2[0]));
    const gsKnotVector<> kv2 (knots2, 2);
    const gsKnotVector<> kv3 (knots2, 2);

    // Control points
    gsMatrix<> ctrlPtsTemp0(3, 2);
    gsMatrix<> ctrlPtsTemp1(3, 2);
    gsMatrix<> ctrlPtsTemp2(7, 2);
    gsMatrix<> ctrlPtsTemp3(7, 2);

    ctrlPtsTemp0 << 1, 0,
                    1, 0.5,
                    1, 1;

    ctrlPtsTemp1 << 0, 0,
                    0, 0.5,
                    0, 1;

    ctrlPtsTemp2 << 0, 1,
                    x1/2, 1,
                    x1, 1,
                    ctrlPtTemp2x, ctrlPtTemp2y,
                    x2, 1,
                    (x2+1)/2, 1,
                    1, 1;

    ctrlPtsTemp3 << 0, 0,
                    x1/2, 0,
                    x1, 0,
                    ctrlPtTemp1x, ctrlPtTemp1y,
                    x2, 0,
                    (x2+1)/2, 0,
                    1, 0;

    const gsBSplineBasis<> basis0( kv0 );
    const gsBSplineBasis<> basis1( kv1 );
    const gsBSplineBasis<> basis2( kv2 );
    const gsBSplineBasis<> basis3( kv3 );

    const gsBSpline<> curveTemp0( basis0, ctrlPtsTemp0 );
    const gsBSpline<> curveTemp1( basis1, ctrlPtsTemp1 );
    const gsBSpline<> curveTemp2( basis2, ctrlPtsTemp2 );
    const gsBSpline<> curveTemp3( basis3, ctrlPtsTemp3 );

    mpTemp->addPatch(gsBSpline<>(curveTemp0));
    mpTemp->addPatch(gsBSpline<>(curveTemp1));
    mpTemp->addPatch(gsBSpline<>(curveTemp2));
    mpTemp->addPatch(gsBSpline<>(curveTemp3));

    return mpTemp;

}

/**************************************************/
gsMultiPatch<>* getCrescendoPlanarSegmentation(real_t x1, real_t x2,
                                               real_t ctrlPtSegm1x, real_t ctrlPtSegm1y,
                                               real_t ctrlPtSegm2x, real_t ctrlPtSegm2y)
{
    gsMultiPatch<>* mpSegm = new gsMultiPatch<>();

    // Knot vectors
    const gsKnotVector<> kv0 (0, 1, 0, 3); // {0, 0, 0, 1, 1, 1};
    const gsKnotVector<> kv1 (0, 1, 0, 3);

    // Control points
    gsMatrix<> ctrlPtsSegm0(3, 2);
    gsMatrix<> ctrlPtsSegm1(3, 2);

    ctrlPtsSegm0 << x1, 0,
                    ctrlPtSegm1x, ctrlPtSegm1y,
                    x1, 1;

    ctrlPtsSegm1 << x2, 0,
                    ctrlPtSegm2x, ctrlPtSegm2y,
                    x2, 1;


    const gsBSplineBasis<> basis0( kv0 );
    const gsBSplineBasis<> basis1( kv1 );

    const gsBSpline<> curveSegm0( basis0, ctrlPtsSegm0 );
    const gsBSpline<> curveSegm1( basis1, ctrlPtsSegm1 );

    mpSegm->addPatch(gsBSpline<>(curveSegm0));
    mpSegm->addPatch(gsBSpline<>(curveSegm1));

    return mpSegm;
}

/**************************************************/
gsMultiPatch<>* getCrescendoParaSplittedTopBottom(const std::vector< gsTensorBSpline<2, real_t> > cpatches)
{
    const gsTensorBSpline<2, real_t> cpatch0Init = cpatches.at(0);
    const gsTensorBSpline<2, real_t> cpatch1Init = cpatches.at(1);
    const gsTensorBSpline<2, real_t> cpatch2Init = cpatches.at(2);

    const gsMatrix<> coefs0Init = cpatch0Init.coefs();
    const gsMatrix<> coefs1Init = cpatch1Init.coefs();
    const gsMatrix<> coefs2Init = cpatch2Init.coefs();

    const index_t nCoefs0 = coefs0Init.rows();
    const index_t nCoefs1 = coefs1Init.rows();
    const index_t nCoefs2 = coefs2Init.rows();

    gsMatrix<> coefs0(nCoefs0, 3);
    coefs0.setZero();
    coefs0.block(0, 0, nCoefs0, 2) = coefs0Init.block(0, 0, nCoefs0, 2);

    gsMatrix<> coefs1(nCoefs1, 3);
    coefs1.setZero();
    coefs1.block(0, 0, nCoefs1, 2) = coefs1Init.block(0, 0, nCoefs1, 2);

    gsMatrix<> coefs2(nCoefs2, 3);
    coefs2.setZero();
    coefs2.block(0, 0, nCoefs2, 2) = coefs2Init.block(0, 0, nCoefs2, 2);

    const gsKnotVector<> kv01 = cpatch0Init.knots(0);
    const gsKnotVector<> kv02 = cpatch0Init.knots(1);
    const gsTensorBSpline<2> cpBottom0(kv01, kv02, coefs0);

    const gsKnotVector<> kv11 = cpatch1Init.knots(0);
    const gsKnotVector<> kv12 = cpatch1Init.knots(1);
    const gsTensorBSpline<2> cpBottom1(kv11, kv12, coefs1);

    const gsKnotVector<> kv21 = cpatch2Init.knots(0);
    const gsKnotVector<> kv22 = cpatch2Init.knots(1);
    const gsTensorBSpline<2> cpBottom2(kv21, kv22, coefs2);

    gsMultiPatch<>* mp = new gsMultiPatch<>();
    mp->addPatch(cpBottom0);
    mp->addPatch(cpBottom1);
    mp->addPatch(cpBottom2);

    // Modifying control points
    for(int i = 0; i < nCoefs0; i++)
    {
        coefs0(i, 2) = 1.0;
    }
    for(int i = 0; i < nCoefs1; i++)
    {
        coefs1(i, 2) = 1.0;
    }
    for(int i = 0; i < nCoefs2; i++)
    {
        coefs2(i, 2) = 1.0;
    }
    const gsTensorBSpline<2> cpTop0(kv01, kv02, coefs0);
    const gsTensorBSpline<2> cpTop1(kv11, kv12, coefs1);
    const gsTensorBSpline<2> cpTop2(kv21, kv22, coefs2);
    mp->addPatch(cpTop0);
    mp->addPatch(cpTop1);
    mp->addPatch(cpTop2);

    writeMultipatch(*mp, "crescendo-tmtf-para-mp_patch4_bottom_top");

    return mp;
}

/**************************************************/
gsMatrix<> boundParameters(const int numPoints, const int id)
{

    gsMatrix<> matZeros(1, numPoints);
    matZeros.setZero();

    gsMatrix<> matOnes(1, numPoints);
    matOnes.setOnes();

    gsVector<> start(1);
    start.setConstant(0.0);
    gsVector<> end(1);
    end.setConstant(1.0);

    gsMatrix<> matUniform = uniformPointGrid(start, end, numPoints);

    gsMatrix<> matParams(2, numPoints);
    matParams.setZero();

    switch(id)
    {
        case 1:
        {
//            gsMatrix<> matParams(2, numPoints);
//            matParams.setZero();
            matParams.row(0) = matZeros;
            matParams.row(1) = matUniform;
            break;
        }
        case 2:
        {
//            gsMatrix<> matParams(2, numPoints);
//            matParams.setZero();
            matParams.row(0) = matOnes;
            matParams.row(1) = matUniform;
            /*
            gsMatrix<> matParams(2, numPoints-2);
            matParams.setZero();
            matParams.row(0) = matOnes.block(0, 1, matOnes.rows(), matOnes.cols()-2);
            matParams.row(1) = matUniform.block(0, 1, matUniform.rows(), matUniform.cols()-2);
            */
            break;
        }
        case 3:
        {
//            gsMatrix<> matParams(2, numPoints);
//            matParams.setZero();
            matParams.row(0) = matUniform;
            matParams.row(1) = matZeros;
            break;
        }
        case 4:
        {
//            gsMatrix<> matParams(2, numPoints);
//            matParams.setZero();
            matParams.row(0) = matUniform;
            matParams.row(1) = matOnes;
            break;
        }
    }

    return matParams;

}

/**************************************************/
void getCrescendoPlanarBoundaryConstraints(const std::string& inputBoundaryFile,
                                           const std::string& inputCoonsFile,
                                           const std::string& inputParamsPointsFile,
                                           const std::string& outputFile)
{

    // Splitted boundary file
    gsMultiPatch<> mpBoundary;
    gsReadFile<>(inputBoundaryFile, mpBoundary);
    mpBoundary.computeTopology(1e-4);
    // Coon's parameterization of the template
    gsMultiPatch<> mpTempCoons;
    gsReadFile<>(inputCoonsFile, mpTempCoons);
    mpTempCoons.computeTopology(1e-4);

    gsFileData<> fd1(inputParamsPointsFile);

    /*
    Numbering of patches:
    +-----+-----+-----+
    |     |     |     |
    |  2  |  1  |  0  |
    |     |     |     |
    +-----+-----+-----+
    */

    gsFileData<> fdPatch;

    gsMatrix<> parameters;
    size_t numPoints;
    gsMatrix<> paramsBoundTarg;
    gsMatrix<> paramsBoundTargTrimmed;

    std::vector< gsMatrix <> > paramsBoundCoons;
    std::vector< gsMatrix <> > ptsBound;
    gsMatrix<> params;
    gsMatrix<> pts;
    gsGeometry<>& gBound = mpBoundary.patch(4);
    // Patch 0
    paramsBoundCoons.clear();
    ptsBound.clear();

    fd1.getId< gsMatrix<> >(0, parameters);
    numPoints = parameters.cols();
    paramsBoundTarg = uniformParameters<2>(0.0, 1.0, numPoints);
    paramsBoundTargTrimmed = paramsBoundTarg.block(0, 1, paramsBoundTarg.rows(), paramsBoundTarg.cols()-1);
    // 0
    /*
    params = boundParameters(numPoints, 2);
    paramsBoundCoons.push_back(params);
    gBound = mpBoundary->patch(0);
    pts = gBound.eval(paramsBoundTarg);
    gsInfo << "params:\n" << params << std::endl;
    gsInfo << "pts:\n" << pts << std::endl;
    ptsBound.push_back(pts);
    */
     // 4
    params = boundParameters(numPoints, 4);
    paramsBoundCoons.push_back(reverseMatrixCols(params));
    gBound = mpBoundary.patch(4);
    pts = gBound.eval(paramsBoundTarg);
    ptsBound.push_back(pts);
    // 7
    params = boundParameters(numPoints, 3);
    paramsBoundCoons.push_back(reverseMatrixCols(params));
    gBound = mpBoundary.patch(7);
    pts = gBound.eval(paramsBoundTarg);
    ptsBound.push_back(pts);

    fdPatch << stdVectorToMatrix(paramsBoundCoons);
    fdPatch << stdVectorToMatrix(ptsBound);
    gsWriteParaviewPoints(stdVectorToMatrix(paramsBoundCoons), "test_par_0");
    gsWriteParaviewPoints(stdVectorToMatrix(ptsBound), "test_pts_0");

    // Patch 1
    paramsBoundCoons.clear();
    ptsBound.clear();

    fd1.getId< gsMatrix<> >(2, parameters);
    numPoints = parameters.cols();
    paramsBoundTarg = uniformParameters<2>(0.0, 1.0, numPoints);
    // 3
    params = boundParameters(numPoints, 4);
    paramsBoundCoons.push_back(reverseMatrixCols(params));
    gBound = mpBoundary.patch(3);
    pts = gBound.eval(paramsBoundTarg);
    ptsBound.push_back(pts);
    // 6
    params = boundParameters(numPoints, 3);
    paramsBoundCoons.push_back(reverseMatrixCols(params));
    gBound = mpBoundary.patch(6);
    pts = gBound.eval(paramsBoundTarg);
    ptsBound.push_back(pts);

    fdPatch << stdVectorToMatrix(paramsBoundCoons);
    fdPatch << stdVectorToMatrix(ptsBound);
    gsWriteParaviewPoints(stdVectorToMatrix(paramsBoundCoons), "test_par_1");
    gsWriteParaviewPoints(stdVectorToMatrix(ptsBound), "test_pts_1");

    // Patch 2
    paramsBoundCoons.clear();
    ptsBound.clear();

    fd1.getId< gsMatrix<> >(4, parameters);
    numPoints = parameters.cols();
    paramsBoundTarg = uniformParameters<2>(0.0, 1.0, numPoints);
    // 2
    params = boundParameters(numPoints, 4);
    paramsBoundCoons.push_back(reverseMatrixCols(params));
    gBound = mpBoundary.patch(2);
    pts = gBound.eval(paramsBoundTarg);
    ptsBound.push_back(pts);
    // 1
    /*
    params = boundParameters(numPoints, 4);
    gsInfo << "params:\n" << params << std::endl;
    paramsBoundCoons.push_back(params);
    gBound = mpBoundary->patch(1);
    pts = gBound.eval(paramsBoundTarg);
    gsInfo << "pts:\n" << pts << std::endl;
    gsWriteParaviewPoints(pts, "test_pts");
    ptsBound.push_back(pts);
    */
    // 5
    params = boundParameters(numPoints, 3);
    paramsBoundCoons.push_back(reverseMatrixCols(params));
    gBound = mpBoundary.patch(5);
    pts = gBound.eval(paramsBoundTarg);
    ptsBound.push_back(pts);

    fdPatch << stdVectorToMatrix(paramsBoundCoons);
    fdPatch << stdVectorToMatrix(ptsBound);
    gsWriteParaviewPoints(stdVectorToMatrix(paramsBoundCoons), "test_par_2");
    gsWriteParaviewPoints(stdVectorToMatrix(ptsBound), "test_pts_2");

    fdPatch.dump(outputFile);

}

} // namespace gismo


/*
// This function changes the middle (curved) part of the Patch 3
// in the planar version of the crescendo geometry
gsMultiPatch<>* getSomething(std::string inputParaFile,
                             std::string inputGeomFile)
{
    // The multipatch file with the planar template
    gsMultiPatch<>* mpPara = readMultipatch(inputParaFile);
    gsGeometry<>& paraPatch0 = mpPara->patch(0);
    gsGeometry<>& paraPatch1 = mpPara->patch(1);
    gsGeometry<>& paraPatch2 = mpPara->patch(2);
    gsGeometry<>& paraPatch3 = mpPara->patch(3);

    // Modified middle (curved) part of Patch 3 in planar template
    gsGeometry<>* modifiedGeom = readGeometry(inputGeomFile);

    // Number of fitting point
    const int numPoints = 1e8;

    // Uniform parameters
    gsMatrix<> params = uniformParameters<2>(0.0, 1.0, numPoints);

    // Fitting points
    gsMatrix<> fittingPts(2, numPoints);
    fittingPts.setZero();
    for(int i = 0; i < numPoints; i++)
    {
        if( (params(i) <= 0.3) || (params(i) >= 0.7) )
        {
            // Point is taken from the original geometry
            gsMatrix<> pt = paraPatch3.eval(params.col(i));
            fittingPts(0, i) = pt(0, 0);
            fittingPts(1, i) = pt(1, 0);
        }
        else if((params(i) > 0.3) && (params(i) < 0.7))
        {
            // Point is taken from the modified geometry
            gsMatrix<> pt = modifiedGeom->eval(params.col(i));
            fittingPts(0, i) = pt(0, 0);
            fittingPts(1, i) = pt(1, 0);
        }
    }
    //gsWriteParaviewPoints(fittingPts, "");

    // Fitting is done in the original basis of the Patch 3
    gsFitting<> fitting(params, fittingPts, paraPatch3.basis());
    fitting.compute(1e-6);
    gsBSpline<>* fittedCurve = static_cast<gsBSpline<>*>(fitting.result()->clone());

    //writeGeometry(*fittedCurve, "test_fittedCurve");

    // Constructing the updated template geometry (multipatch) file
    gsMultiPatch<>* mpNewPara = new gsMultiPatch<>();
    mpNewPara->addPatch(paraPatch0); // Patch 0: keep unchanged
    mpNewPara->addPatch(paraPatch1); // Patch 1: keep unchanged

    gsMatrix<> coefsParaPatch3 = paraPatch3.coefs();
    // Updating the coefficients of Patch 3
    gsMatrix<> tmpCoefs = (fittedCurve->coefs()).block(13, 0, 33, 2);
    coefsParaPatch3.block(13, 0, 33, 2) = tmpCoefs;
    paraPatch3.setCoefs(coefsParaPatch3);

    // Updating the coefficients of Patch 2
    gsMatrix<> coefsParaPatch2 = paraPatch2.coefs();
    coefsParaPatch2.block(13, 0, 33, 2) = transformCtrlPts(tmpCoefs);
    paraPatch2.setCoefs(coefsParaPatch2);

    mpNewPara->addPatch(paraPatch2);
    mpNewPara->addPatch(paraPatch3);

    return mpNewPara;

}
*/


/**************************************************/
/*
// Creates segmentation lines for the planar template
gsMultiPatch<>* getSegmentationPlanar(void)
{
    std::string inputSegm1("/home/cirion/svajunas/motor/geometries/crescendo-1p/SegmentationTemplate/segmTemplate1Pro.xml");
    std::string inputSegm2("/home/cirion/svajunas/motor/geometries/crescendo-1p/SegmentationTemplate/segmTemplate2Pro.xml");

    gsMultiPatch<>* mp = new gsMultiPatch<>();

    const gsGeometry<>* segm1 = readGeometry(inputSegm1);
    const gsGeometry<>* segm2 = readGeometry(inputSegm2);
    const gsBSpline<>& volSegm1 = static_cast< const gsBSpline<>& >(*segm1);
    const gsBSpline<>& volSegm2 = static_cast< const gsBSpline<>& >(*segm2);

    int deg = 2;

    // Project segmentation surfaces into xy-plane
    gsMatrix<> ctrlPts1 = volSegm1.coefs();
    gsMatrix<> ctrlPts2 = volSegm2.coefs();

    // Getting control points
    gsMatrix<> ctrlPtsPlanar1 = ctrlPts1.block(10506, 0, 103, 2);
    gsMatrix<> ctrlPtsPlanar2 = ctrlPts2.block(10506, 0, 103, 2);

    // Knot vectors
    gsKnotVector<> kv1(0.0, 1.0, ctrlPtsPlanar1.rows()-deg-1, deg+1);
    gsKnotVector<> kv2(0.0, 1.0, ctrlPtsPlanar2.rows()-deg-1, deg+1);

    // Constructing segmentation lines
    gsBSpline<> planarSegm1(kv1, ctrlPtsPlanar1);
    gsBSpline<> planarSegm2(kv2, ctrlPtsPlanar2);

    mp->addPatch(planarSegm1);
    mp->addPatch(planarSegm2);

    return mp;
}
*/

/**************************************************/
/*
gsMultiPatch<>* getModifiedPlanarTemplate(std::string inputParaFile)
{
    // Reading multipatch:
    gsMultiPatch<>* mpPara = readMultipatch(inputParaFile);

    // Template geometry boundaries
    const gsGeometry<>& paraPatch0 = mpPara->patch(0); // Vertical line, right
    const gsGeometry<>& paraPatch1 = mpPara->patch(1); // Vertical line, left
    gsGeometry<>& paraPatch2 = mpPara->patch(2); // Horizontal curve, upper
    gsGeometry<>& paraPatch3 = mpPara->patch(3); // Horizontal curve, lower
    const gsBSpline<>& geom0 = dynamic_cast< const gsBSpline<>& >(paraPatch0);
    const gsBSpline<>& geom1 = dynamic_cast< const gsBSpline<>& >(paraPatch1);
    gsBSpline<>& geom2 = dynamic_cast< gsBSpline<>& >(paraPatch2);
    gsBSpline<>& geom3 = dynamic_cast< gsBSpline<>& >(paraPatch3);

    // Parameters corresponding to C^0-cont. points
    gsMatrix<> paramsParaGeom2(1, 2);
    gsMatrix<> paramsParaGeom3(1, 2);
    paramsParaGeom2(0, 0) = 0.3; paramsParaGeom2(0, 1) = 0.7;
    paramsParaGeom3(0, 0) = 0.3; paramsParaGeom3(0, 1) = 0.7;

    gsMatrix<> ctrlPts2 = geom2.coefs();
    gsMatrix<> ctrlPts3 = geom3.coefs();

    gsKnotVector<> kv2 = geom2.knots();
    gsKnotVector<> kv3 = geom3.knots();

    // ******************************************************** //
    index_t ind21 = 0;
    for(ind21 = 0; kv2.at(ind21) < paramsParaGeom2(0); ind21++);
    ind21--;

    index_t ind22 = 0;
    for(ind22 = 0; kv2.at(ind22) < paramsParaGeom2(1); ind22++);
    ind22--;

    gsInfo << "ctrlPts2.row(ind21): " << ctrlPts2.row(ind21) << std::endl;
    gsInfo << "ctrlPts2.row(ind22): " << ctrlPts2.row(ind22) << std::endl;

    for(int i = ind21; i < ind22+1; i++)
    {
        ctrlPts2(i, 0) = ctrlPts2(i, 0) - 0.1;
        ctrlPts2(i, 1) = ctrlPts2(i, 1);
    }

    geom2.setCoefs(ctrlPts2);

    // ******************************************************** //
    index_t ind31 = 0;
    for(ind31 = 0; kv3.at(ind31) < paramsParaGeom3(0); ind31++);
    ind31--;

    index_t ind32 = 0;
    for(ind32 = 0; kv3.at(ind32) < paramsParaGeom3(1); ind32++);
    ind32--;

    gsInfo << "ctrlPts3.row(ind31): " << ctrlPts3.row(ind31) << std::endl;
    gsInfo << "ctrlPts3.row(ind32): " << ctrlPts3.row(ind32) << std::endl;

    for(int i = ind31; i < ind32+1; i++)
    {
        ctrlPts3(i, 0) = ctrlPts3(i, 0) - 0.1;
        ctrlPts3(i, 1) = ctrlPts3(i, 1);
    }

    geom3.setCoefs(ctrlPts3);

    gsMultiPatch<>* mpParaGeom = new gsMultiPatch<>();
    mpParaGeom->addPatch(geom0);
    mpParaGeom->addPatch(geom1);
    mpParaGeom->addPatch(geom2);
    mpParaGeom->addPatch(geom3);

    return mpParaGeom;
}
*/


/**************************************************/
// UNDER CONSTRUCTION!!!
// NOT FINISHED!
/*
void unifyKnotsVolumetric(gsTensorBSpline<2, real_t>& geom1,
                          gsTensorBSpline<2, real_t>& geom2)
{
    //gsTensorBSpline<2, real_t> inputGeometry1 = dynamic_cast< const gsTensorBSpline<2, real_t>& >(geom1);
    //gsTensorBSpline<2, real_t> inputGeometry2 = dynamic_cast< const gsTensorBSpline<2, real_t>& >(geom2);

    // Degree unification (if necessary)
    int deg11 = geom1.knots(0).degree(); // Geometry 1, Direction 1
    int deg12 = geom1.knots(1).degree(); // Geometry 1, Direction 2
    int deg21 = geom2.knots(0).degree(); // Geometry 2, Direction 1
    int deg22 = geom2.knots(1).degree(); // Geometry 2, Direction 2

    if (deg11 < deg21)
    {
        geom1.degreeElevate(deg21-deg11, 0);
    }
    else if (deg21 < deg11)
    {
        geom2.degreeElevate(deg11-deg21, 0);
    }

    if (deg12 < deg22)
    {
        geom1.degreeElevate(deg22-deg12, 1);
    }
    else if (deg22 < deg12)
    {
        geom2.degreeElevate(deg12-deg22, 1);
    }

    gsKnotVector<> kv11 = geom1.knots(0); // Geometry 1, Direction 1
    gsKnotVector<> kv12 = geom1.knots(1); // Geometry 1, Direction 2
    gsKnotVector<> kv21 = geom2.knots(0); // Geometry 2, Direction 1
    gsKnotVector<> kv22 = geom2.knots(1); // Geometry 2, Direction 2

    gsInfo << "Before refinement:\n" << std::endl;
    gsInfo << "geom1.knots(0):\n" << kv11 << std::endl;
    gsInfo << "geom1.knots(1):\n" << kv12 << std::endl;
    gsInfo << "geom2.knots(0):\n" << kv21 << std::endl;
    gsInfo << "geom2.knots(1):\n" << kv22 << std::endl;

    if ((kv11.degree() != kv21.degree()) || (kv12.degree() != kv22.degree()))
    {
        gsWarn << "Degrees mismatch" << std::endl;
    }

    std::vector<real_t> kv_temp1; // Direction 1
    std::vector<real_t> kv_temp2; // Direction 2

    // Direction 1
    std::set_union(kv11.begin(), kv11.end(),
                   kv21.begin(), kv21.end(),
                   std::back_inserter(kv_temp1));
    gsKnotVector<> kv1(kv_temp1, kv11.degree());

    // Direction 2
    std::set_union(kv12.begin(), kv12.end(),
                   kv22.begin(), kv22.end(),
                   std::back_inserter(kv_temp2));
    gsKnotVector<> kv2(kv_temp2, kv12.degree());

    // Geometry 1, Direction 1
    std::vector<real_t> diff11;
    std::set_difference(kv1.begin(), kv1.end(), kv11.begin(), kv11.end(),
                        std::inserter(diff11, diff11.begin()));

     // Geometry 1, Direction 2
    std::vector<real_t> diff12;
    std::set_difference(kv2.begin(), kv2.end(), kv12.begin(), kv12.end(),
                        std::inserter(diff12, diff12.begin()));

    // Geometry 2, Direction 1
    std::vector<real_t> diff21;
    std::set_difference(kv1.begin(), kv1.end(), kv21.begin(), kv21.end(),
                        std::inserter(diff21, diff21.begin()));

    // Geometry 2, Direction 2
    std::vector<real_t> diff22;
    std::set_difference(kv2.begin(), kv2.end(), kv22.begin(), kv22.end(),
                        std::inserter(diff22, diff22.begin()));

    std::vector< std::vector<real_t> > diff1;
    diff1.push_back(diff11);
    diff1.push_back(diff12);
    std::vector< std::vector<real_t> > diff2;
    diff2.push_back(diff21);
    diff2.push_back(diff22);

    if (diff1.at(0).size() != 0 || diff1.at(1).size() != 0)
    {
        gsInfo << "diff1.at(0)" << std::endl;
        printVector(diff1.at(0));
        gsInfo << "diff2.at(0)" << std::endl;
        printVector(diff1.at(1));
        geom1.basis().refine_withCoefs(geom1.coefs(), diff1);
    }
    if (diff2.at(0).size() != 0 || diff2.at(1).size() != 0)
    {
        gsInfo << "diff2.at(0)" << std::endl;
        printVector(diff2.at(0));
        gsInfo << "diff2.at(1)" << std::endl;
        printVector(diff2.at(1));
        geom2.basis().refine_withCoefs(geom2.coefs(), diff2);
    }

    kv11 = geom1.knots(0); // Geometry 1, Direction 1
    kv12 = geom1.knots(1); // Geometry 1, Direction 2
    kv21 = geom2.knots(0); // Geometry 2, Direction 1
    kv22 = geom2.knots(1); // Geometry 2, Direction 2

    gsInfo << "After refinement:\n" << std::endl;
    gsInfo << "geom1.knots(0):\n" << kv11 << std::endl;
    gsInfo << "geom1.knots(1):\n" << kv12 << std::endl;
    gsInfo << "geom2.knots(0):\n" << kv21 << std::endl;
    gsInfo << "geom2.knots(1):\n" << kv22 << std::endl;

}
*/
/**************************************************/
/*
// This function does axisymmetrical transformation
// of certain control points from Patch 3 in planar template
gsMatrix<> transformCtrlPts(gsMatrix<>& ctrlPts)
{
    index_t numPts = ctrlPts.rows();
    // Resulting matrix of coefficients
    gsMatrix<> transformedCtrlPts(numPts, ctrlPts.cols());
    transformedCtrlPts.setZero();

    // Vertical symmetry axis
    real_t axis_x = 0.5;
    // Horizontal symmetry axis
    real_t axis_y = 0.5;
    for (int i = 0; i < numPts; i++)
    {
        transformedCtrlPts(i, 0) = 2*axis_x - ctrlPts(numPts-1-i, 0);
        transformedCtrlPts(i, 1) = 2*axis_y - ctrlPts(numPts-1-i, 1);
    }

    return transformedCtrlPts;
}
*/
