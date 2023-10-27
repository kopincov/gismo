/** @file gsCrescendoTest.cpp

    @brief This file contains various tests with crescendo-tmtf geometry

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
//#include "gsMotorIOUtils.h"

using namespace gismo;

// Takes crescendo-tmtf geometry (template OR target),
// extracts Patch 4 and does some preparation for 2-D parametrization:
// - identifies Patch 4 boundaries
// - save boundaries in gsMultiPatch format
gsMultiPatch<>* getCrescendoPatch4Boundaries(const std::string& inputFile);

// Creates segmentation line for the planar template
gsMultiPatch<>* getSegmentationPlanar( void );

// Takes crescendo-tmtfc-para-mp Patch 4 boundaries
// (prepared with getCrescendoPatch4Boundaries)
// and splits them at certain points determined by segmentation lines
gsMultiPatch<>* splitCrescendoParaPatch4PlanarBoundaries(const std::string& inputParaFile);

// ...
gsMultiPatch<>* splitCrescendoParaPatch4VolumetricBoundaries(const std::string& inputParaFile,
                                                             const std::string& inputSegmFile);

// Splits B-Spline curve at points corresponding
// to given parameteric values
//gsMultiPatch<>* splitOpenBSpline(const gsBSpline<>& geom, const gsMatrix<>& pars);

// Takes two multipatches and joins them
//gsMultiPatch<>* joinMultipatches(const gsMultiPatch<>& mp1, const gsMultiPatch<>& mp2);

// Returns ...
std::vector< gsMultiPatch<> > getCrescendoParaSkeletonPlanarBoundaries(const gsMultiPatch<>& mpTemp,
                                                                       const gsMultiPatch<>& mpSegm);

// Returns ...
std::vector< gsMultiPatch<> > getCrescendoParaSkeletonVolumetricBoundaries(const gsMultiPatch<>& mpTemp,
                                                                           const gsMultiPatch<>& mpSegm);

// Constructs an std::vector of Coon's patches defined by boundaries stored in std::vector
//std::vector< gsTensorBSpline<2, real_t> > getCoonsPatchesPlanar( const std::vector< gsMultiPatch<> > boundaries);
//std::vector< gsTensorBSpline<3, real_t> > getCoonsPatchesVolumetric(const std::vector< gsMultiPatch<> > boundaries);
// Applies map for Coon's patches stored in std::vector
//void mapCoonsPatches(const gsGeometry<>* map,
//                     const std::vector< gsTensorBSpline<2, real_t> > cpatches);

// ...
gsMultiPatch<>* getModifiedPlanarTemplate(std::string inputParaFile);

void unifyKnotsVolumetric(gsTensorBSpline<2, real_t>& geom1,
                          gsTensorBSpline<2, real_t>& geom2);

/**************************************************/
template<unsigned d>
void mapCoonsPatches(const gsGeometry<>* map,
                     const std::vector< gsTensorBSpline<d, real_t> > cpatches)
{
    gsMultiPatch<> mpCoonsPatches;
    gsMatrix<> parameters = uniformParameters<d+1>(0.0, 1.0, 100);
    gsFileData<> fdPointsParameters;

    for (std::size_t i = 0; i != cpatches.size(); i++)
    {
        gsTensorBSpline<d, real_t> cp = cpatches.at(i);
        gsMatrix<> cpValues = cp.eval(parameters);
        gsMatrix<> cpMappedValues = map->eval(cpValues);

        mpCoonsPatches.addPatch(cp);
        fdPointsParameters << parameters;
        fdPointsParameters << cpMappedValues;
    }
    mpCoonsPatches.computeTopology(1e-9, true);
    writeMultipatch(mpCoonsPatches, "coons_multipatch");
    fdPointsParameters.dump("coons_points_and_parameters");
}

// Constructs ...
gsMultiPatch<>* getCrescendo25D(const std::string& inputFile,
                                const real_t height);

// Constructs a multipatch for splitted bottom and top parts of
// the volumetric template geometry
gsMultiPatch<>* getCrescendoParaSplittedTopBottom(const std::vector< gsTensorBSpline<2, real_t> > cpatches);

/******************************************************************/

int main(int argc, char *argv[])
{

    int mode = 0;

    gsCmdLine cmd("Routines for working with crescendo-tmtf demonstrator");
    cmd.addInt("m", "mode", "Mode:\n"
                            " 0 - extract Patch 4 from the template and target geometries\n"
                            " 1 - prepare planar segmentation lines\n"
                            " 2 - make 2.5-D template and target geometries\n"
                            " 3 - split boundaries of planar template geometry\n"
                            " 4 - split boundaries of volumetric template geometry\n"
                            " 5 - construct and map planar Coons' patches\n"
                            " 6 - construct and map volumetric Coons' patches\n"
                            " 9 - testing mode", mode);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    switch(mode)
    {
        case 0: // Extracting Patch 4 from the template and target geometries
        {
            gsInfo << "* Extracting Patch 4 from the template geometry" << std::endl;
            std::string inputParaBoundFile("//home/turin/svajunas/motor_/geometries/crescendo-1p/crescendo-tmtf-para-mp.xml");
            std::string outputParaBoundFile("crescendo-tmtf-para-mp_patch4_bound_planar");
            gsMultiPatch<>* mpParaBound = getCrescendoPatch4Boundaries(inputParaBoundFile);
            writeMultipatch(*mpParaBound, outputParaBoundFile);

            gsInfo << "* Extracting Patch 4 from the target geometry" << std::endl;
            std::string inputGeoBoundFile("//home/turin/svajunas/motor_/geometries/crescendo-1p/crescendo-tmtf-geo-mp.xml");
            std::string outputGeoBoundFile("crescendo-tmtf-geo-mp_patch4_bound_planar");
            gsMultiPatch<>* mpGeoBound = getCrescendoPatch4Boundaries(inputGeoBoundFile);
            writeMultipatch(*mpGeoBound, outputGeoBoundFile);break;
            break;
        }
        case 1: // Prepearing planar segmentation lines
        {
            gsInfo << "* Prepearing planar segmentation lines" << std::endl;
            gsMultiPatch<>* mpSegmPlanar = getSegmentationPlanar();
            writeMultipatch(*mpSegmPlanar, "crescendo-tmtf-para-mp_segm_planar");
            break;
        }
        case 2: // Making 2.5-D template and target geometries
        {
            gsInfo << "* Making 2.5-D template geometry" << std::endl;
            std::string inputParaVolFile("crescendo-tmtf-para-mp_patch4_bound_planar.xml");
            std::string outputParaVolFile("crescendo-tmtf-para-mp_patch4_bound_vol");
            gsMultiPatch<>* mpParaVol = getCrescendo25D(inputParaVolFile, 1.0);
            writeMultipatch(*mpParaVol, outputParaVolFile);

            gsInfo << "* Making 2.5-D target geometry" << std::endl;
            std::string inputGeoVolFile("crescendo-tmtf-geo-mp_patch4_bound_planar.xml");
            std::string outputGeoVolFile("crescendo-tmtf-geo-mp_patch4_bound_vol");
            gsMultiPatch<>* mpGeoVol = getCrescendo25D(inputGeoVolFile, 0.1);
            writeMultipatch(*mpGeoVol, outputGeoVolFile);
            break;
        }
        case 3: // Splitting boundaries of planar template geometry
        {
            gsInfo << "* Splitting boundaries of planar template geometry" << std::endl;
            std::string inputParaBoundFile("crescendo-tmtf-para-mp_patch4_bound_planar.xml");
            std::string outputParaBoundFile("crescendo-tmtf-para-mp_patch4_bound_planar_splitted");
            gsMultiPatch<>* mpParaPlanarSplitted = splitCrescendoParaPatch4PlanarBoundaries(inputParaBoundFile);
            writeMultipatch(*mpParaPlanarSplitted, outputParaBoundFile);
            break;
        }
        case 4: // Splitting boundaries of volumetric template geometry
        {
            gsInfo << "* Splitting boundaries of volumetric template geometry" << std::endl;
            std::string inputParaBoundFile("crescendo-tmtf-para-mp_patch4_bound_planar.xml");
            std::string inputParaSegmFile("crescendo-tmtf-para-mp_segm_planar.xml");
            std::string outputParaVolBoundFile("crescendo-tmtf-para-mp_patch4_bound_vol_splitted");
            gsMultiPatch<>* mpParaVolSplitted = splitCrescendoParaPatch4VolumetricBoundaries(inputParaBoundFile, inputParaSegmFile);
            writeMultipatch(*mpParaVolSplitted, outputParaVolBoundFile);
            break;
        }
        case 5: // Constructing and mapping planar Coons' patches
        {
            gsInfo << "* Constructing planar Coon's patches" << std::endl;
            gsMultiPatch<>* mpTemp = readMultipatch("crescendo-tmtf-para-mp_patch4_bound_planar_splitted.xml");
            gsMultiPatch<>* mpSegm = readMultipatch("crescendo-tmtf-para-mp_segm_planar.xml");
            gsInfo << "  ... getCrescendoParaSkeletonPlanarBoundaries" << std::endl;
            std::vector< gsMultiPatch<> > boundaries = getCrescendoParaSkeletonPlanarBoundaries(*mpTemp, *mpSegm);
            gsInfo << "  ... getCoonsPatchesPlanar" << std::endl;
            std::vector< gsTensorBSpline<2, real_t> > cpatches = getCoonsPatchesPlanar(boundaries);

            gsInfo << "* Mapping planar Coon's patches" << std::endl;
            gsGeometry<>* map = readGeometry("planar_param_final_THBMap_0.xml");
            mapCoonsPatches<2>(map, cpatches);
            break;
        }
        case 6: // Constructing and mapping volumetric Coons' patches
        {
            gsInfo << "* Constructing volumetric Coon's patches" << std::endl;
            gsInfo << "  ... Reading multipatches" << std::endl;
            gsMultiPatch<>* mpVolTemp = readMultipatch("crescendo-tmtf-para-mp_patch4_bound_vol_splitted.xml");
            gsMultiPatch<>* mpVolSegm = readMultipatch("//home/turin/svajunas/motor_/geometries/crescendo-1p/SegmentationTemplate/segmTemplatePro.xml");
            gsInfo << "  ... getCrescendoParaSkeletonVolumetricBoundaries" << std::endl;
            std::vector< gsMultiPatch<> > boundariesVol = getCrescendoParaSkeletonVolumetricBoundaries(*mpVolTemp, *mpVolSegm);
            gsInfo << "  ... getCoonsPatchesVolumetric" << std::endl;
            std::vector< gsTensorBSpline<3, real_t> > cpatchesVol = getCoonsPatchesVolumetric(boundariesVol);

            gsInfo << "* Mapping volumetric Coon's patches" << std::endl;
            gsGeometry<>* map = readGeometry("vol_param_final_THBMap_0.xml");
            mapCoonsPatches<3>(map, cpatchesVol);
            break;
        }
        case 7: // Modifying planar template
        {
            gsInfo << "* Modifying planar template geometry" << std::endl;
            std::string inputParaBoundFile("crescendo-tmtf-para-mp_patch4_bound_planar.xml");
            std::string outputParaBoundFile("crescendo-tmtf-para-mp_patch4_bound_planar_modified");
            gsMultiPatch<>* mpTempModified = getModifiedPlanarTemplate(inputParaBoundFile);
            writeMultipatch(*mpTempModified, outputParaBoundFile);
            break;
        }
        case 9: // Tests
        {

            gsInfo << "* Constructing volumetric Coon's patches" << std::endl;
            gsInfo << "  ... Reading multipatches" << std::endl;
            gsMultiPatch<>* mpVolTemp = readMultipatch("crescendo-tmtf-para-mp_patch4_bound_vol_splitted.xml");
            gsMultiPatch<>* mpVolSegm = readMultipatch("//home/turin/svajunas/motor_/geometries/crescendo-1p/SegmentationTemplate/segmTemplatePro.xml");
            gsInfo << "  ... getCrescendoParaSkeletonVolumetricBoundaries" << std::endl;
            std::vector< gsMultiPatch<> > boundariesVol = getCrescendoParaSkeletonVolumetricBoundaries(*mpVolTemp, *mpVolSegm);

            //gsMultiPatch<> test_boundary = boundariesVol.at(0);
            //gsInfo << "test_boundary:\n" << test_boundary << std::endl;
            //gsInfo << "test_boundary[0]:\n" << test_boundary[1] << std::endl;

            writeMultipatch(boundariesVol.at(0), "test_boundary0");
            writeMultipatch(boundariesVol.at(1), "test_boundary1");
            writeMultipatch(boundariesVol.at(2), "test_boundary2");

            gsInfo << "* Constructing volumetric Coon's patches" << std::endl;
            gsInfo << "  ... Reading multipatches" << std::endl;
            gsMultiPatch<>* mpBound0 = readMultipatch("test_boundary0.xml");
            std::vector< gsMultiPatch<> > testBoundariesVol;
            testBoundariesVol.push_back(*mpBound0);
            std::vector< gsTensorBSpline<3, real_t> > cpatchesVol = getCoonsPatchesVolumetric(testBoundariesVol);
            /*
            gsInfo << "* Mapping volumetric Coon's patches" << std::endl;
            gsGeometry<>* map = readGeometry("vol_param_final_THBMap_0.xml");
            mapCoonsPatches<3>(map, cpatchesVol);
            */

            /*
            gsTensorBSpline<2, real_t> bound0 = static_cast< const gsTensorBSpline<2, real_t>& >(test_boundary[0]);
            gsTensorBSpline<2, real_t> bound1 = static_cast< const gsTensorBSpline<2, real_t>& >(test_boundary[1]);
            gsTensorBSpline<2, real_t> bound2 = static_cast< const gsTensorBSpline<2, real_t>& >(test_boundary[2]);
            gsTensorBSpline<2, real_t> bound3 = static_cast< const gsTensorBSpline<2, real_t>& >(test_boundary[3]);
            gsTensorBSpline<2, real_t> bound4 = static_cast< const gsTensorBSpline<2, real_t>& >(test_boundary[4]);
            gsTensorBSpline<2, real_t> bound5 = static_cast< const gsTensorBSpline<2, real_t>& >(test_boundary[5]);

            gsInfo << "bound0:\n" << bound0 << std::endl;
            gsInfo << "bound1:\n" << bound1 << std::endl;
            gsInfo << "bound2:\n" << bound2 << std::endl;
            gsInfo << "bound3:\n" << bound3 << std::endl;
            gsInfo << "bound4:\n" << bound4 << std::endl;
            gsInfo << "bound5:\n" << bound5 << std::endl;

            gsInfo << "Unifying knots (0, 1):\n" << std::endl;
            unifyKnotsVolumetric(bound0, bound1);
            gsInfo << "Unifying knots (0, 2):\n" << std::endl;
            unifyKnotsVolumetric(bound0, bound2);
            gsInfo << "Unifying knots (0, 3):\n" << std::endl;
            unifyKnotsVolumetric(bound0, bound3);
            gsInfo << "Unifying knots (0, 4):\n" << std::endl;
            unifyKnotsVolumetric(bound0, bound3);
            gsInfo << "Unifying knots (0, 5):\n" << std::endl;
            unifyKnotsVolumetric(bound0, bound3);

            gsInfo << "Unifying knots (1, 2):\n" << std::endl;
            unifyKnotsVolumetric(bound1, bound2);
            gsInfo << "Unifying knots (1, 3):\n" << std::endl;
            unifyKnotsVolumetric(bound1, bound3);
            gsInfo << "Unifying knots (1, 4):\n" << std::endl;
            unifyKnotsVolumetric(bound1, bound3);
            gsInfo << "Unifying knots (1, 5):\n" << std::endl;
            unifyKnotsVolumetric(bound1, bound3);

            gsInfo << "Unifying knots (2, 3):\n" << std::endl;
            unifyKnotsVolumetric(bound2, bound3);
            gsInfo << "Unifying knots (2, 4):\n" << std::endl;
            unifyKnotsVolumetric(bound2, bound3);
            gsInfo << "Unifying knots (2, 5):\n" << std::endl;
            unifyKnotsVolumetric(bound2, bound3);

            gsInfo << "Unifying knots (3, 4):\n" << std::endl;
            unifyKnotsVolumetric(bound2, bound3);
            gsInfo << "Unifying knots (3, 5):\n" << std::endl;
            unifyKnotsVolumetric(bound2, bound3);

            gsInfo << "Unifying knots (4, 5):\n" << std::endl;
            unifyKnotsVolumetric(bound2, bound3);

            */

            break;
        }
        default:
            gsInfo << "Invalid value of mode (should be integer between 0 and ...)" << std::endl;
            return 1;
    }

    return 0;
}

/**************************************************/
gsMultiPatch<>* getCrescendoPatch4Boundaries(const std::string& inputFile)
{
    // Reading multipatch:
    // The template geometry
    gsMultiPatch<> * mp = readMultipatch(inputFile);

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

    gsBSpline<> paraBound0(geomBound0.knots(), ctrlPtsBound0);
    gsBSpline<> paraBound1(geomBound1.knots(), ctrlPtsBound1);
    gsBSpline<> paraBound2(geomBound3.knots(), ctrlPtsBound2); // Don't worry this is not bug - this should be like this
    gsBSpline<> paraBound3(geomBound3.knots(), ctrlPtsBound3); // Constructing and mapping volumetric Coons' patchesknots(), ctrlPtsBound3);

    // Multi-patch for boundaries
    gsMultiPatch<>* mpBound = new gsMultiPatch<>();
    mpBound->addPatch(paraBound0);
    mpBound->addPatch(paraBound1);
    mpBound->addPatch(paraBound2);
    mpBound->addPatch(paraBound3);

    return mpBound;
}

/**************************************************/
gsMultiPatch<>* getSegmentationPlanar(void)
{
    std::string inputSegm1("//home/turin/svajunas/motor_/geometries/crescendo-1p/SegmentationTemplate/segmTemplate1Pro.xml");
    std::string inputSegm2("//home/turin/svajunas/motor_/geometries/crescendo-1p/SegmentationTemplate/segmTemplate2Pro.xml");

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

/**************************************************/
gsMultiPatch<>* splitCrescendoParaPatch4PlanarBoundaries(const std::string& inputParaFile)
{

    // Reading multipatch:
    gsMultiPatch<>* mpPara = readMultipatch(inputParaFile);

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
    // Only paraGeom2 and paraGeom3 are to be splitted
    // TO DO: Fix this place
    mpParaGeom->addPatch(paraGeom0);
    mpParaGeom->addPatch(paraGeom1);

    gsMultiPatch<>* mpParaGeom2 = splitOpenBSpline(paraGeom2, paramsParaGeom2);
    gsMultiPatch<>* mpParaGeom3 = splitOpenBSpline(paraGeom3, paramsParaGeom3);
    mpParaGeom = joinMultipatches(*mpParaGeom, *mpParaGeom2);
    mpParaGeom = joinMultipatches(*mpParaGeom, *mpParaGeom3);

    return mpParaGeom;
}

/**************************************************/
gsMultiPatch<>* splitCrescendoParaPatch4VolumetricBoundaries(const std::string& inputParaFile,
                                                             const std::string& inputSegmFile)
{
    gsMultiPatch<>* mpSegm = readMultipatch(inputSegmFile);

    gsMultiPatch<>* mpParaPlanarSplitted = splitCrescendoParaPatch4PlanarBoundaries(inputParaFile);
    gsMultiPatch<>* mpParaVol = extrudeMultipatch(*mpParaPlanarSplitted, 1.0);

    std::vector< gsMultiPatch<> > boundaries = getCrescendoParaSkeletonPlanarBoundaries(*mpParaPlanarSplitted, *mpSegm);
    std::vector< gsTensorBSpline<2, real_t> > cpatches = getCoonsPatchesPlanar(boundaries);
    gsMultiPatch<>* crescendoTopBottom = getCrescendoParaSplittedTopBottom(cpatches);

    mpParaVol = joinMultipatches(*mpParaVol, *crescendoTopBottom);

    return mpParaVol;
}


/**************************************************/
/*
gsMultiPatch<>* splitOpenBSpline(const gsBSpline<>& geom, const gsMatrix<>& pars)
{
    gsMatrix<> ctrlPts = geom.coefs();
    gsKnotVector<> kv = geom.knots();
    int deg = geom.degree();

    gsMultiPatch<>* mp = new gsMultiPatch<>();

    // Beginning of the B-Spline curve
    index_t ind1 = 0;
    index_t ind2 = 0;
    for(ind2 = 0; kv.at(ind2) < pars(0); ind2++);
    ind2--;

    gsMatrix<> ctrlPtsTmp1 = ctrlPts.block(ind1, 0, ind2-ind1+1, 2);
    gsKnotVector<> kvTmp1(0.0, 1.0, ctrlPtsTmp1.rows()-deg-1, deg+1);

    gsBSpline<> patchCurve1(kvTmp1, ctrlPtsTmp1);
    mp->addPatch(patchCurve1);

    // Inner B-Spline curves
    for (int i = 0; i < pars.cols()-1; i++)
    {
        ind1 = 0;
        for(ind1 = 0; kv.at(ind1) < pars(i); ind1++);
        ind1--;

        ind2 = 0;
        for(ind2 = 0; kv.at(ind2) < pars(i+1); ind2++);
        ind2--;

        gsMatrix<> ctrlPtsTmp = ctrlPts.block(ind1, 0, ind2-ind1+1, 2);
        gsKnotVector<> kvTmp(0.0, 1.0, ctrlPtsTmp.rows()-deg-1, deg+1);
        gsBSpline<> patchCurve(kvTmp, ctrlPtsTmp);
        mp->addPatch(patchCurve);
    }

    // End of the B-Spline curve
    ind1 = ind2;
    ind2 = kv.size()-deg-2;

    gsMatrix<> ctrlPtsTmp2 = ctrlPts.block(ind1, 0, ind2-ind1+1, 2);
    gsKnotVector<> kvTmp2(0.0, 1.0, ctrlPtsTmp2.rows()-deg-1, deg+1);

    gsBSpline<> patchCurve2(kvTmp2, ctrlPtsTmp2);
    mp->addPatch(patchCurve2);

    return mp;
}
*/
/**************************************************/
/*
gsMultiPatch<>* joinMultipatches(const gsMultiPatch<>& mp1, const gsMultiPatch<>& mp2)
{
    gsMultiPatch<>* mp = new gsMultiPatch<>();
    // Adding mp1 patches to mp
    for (index_t i = 0; i != mp1.nPatches(); i++)
    {
        const gsGeometry<>& patch = mp1[i];
        mp->addPatch(patch);
    }
    // Adding mp2 patches to mp
    for (index_t i = 0; i != mp2.nPatches(); i++)
    {
        const gsGeometry<>& patch = mp2[i];
        mp->addPatch(patch);
    }
    return mp;
}
*/

/**************************************************/
std::vector< gsMultiPatch<> > getCrescendoParaSkeletonPlanarBoundaries(const gsMultiPatch<>& mpTemp,
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
    boundaries.push_back(boundary0);

    boundary1.addPatch(mpSegm[1]);
    boundary1.addPatch(mpTemp[3]);
    boundary1.addPatch(mpSegm[0]);
    boundary1.addPatch(mpTemp[6]);
    boundaries.push_back(boundary1);

    boundary2.addPatch(mpSegm[0]);
    boundary2.addPatch(mpTemp[2]);
    boundary2.addPatch(mpTemp[1]);
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
/*
std::vector< gsTensorBSpline<2, real_t> > getCoonsPatchesPlanar(const std::vector< gsMultiPatch<> > boundaries)
{
    std::vector< gsMultiPatch<> > bound = unifyKnotsBoundaries(boundaries);
    std::vector< gsTensorBSpline<2, real_t> > cpatches = coonsPatches<2>(bound);

    return cpatches;
}
*/

/**************************************************/
std::vector< gsTensorBSpline<3, real_t> > getCoonsPatchesVolumetric(const std::vector< gsMultiPatch<> > boundaries)
{
    std::vector< gsMultiPatch<> > bound = boundaries;
    /*
    gsInfo << "getCoonsPatchesVolumetric: bound.at(0): " << bound.at(0).nBoundary() << std::endl;
    gsInfo << "getCoonsPatchesVolumetric: bound.at(1): " << bound.at(1).nBoundary() << std::endl;
    gsInfo << "getCoonsPatchesVolumetric: bound.at(2): " << bound.at(2).nBoundary() << std::endl;
    */
    std::vector< gsTensorBSpline<3, real_t> > cpatches = coonsPatches<3>(bound);

    return cpatches;

}

/**************************************************/
gsMultiPatch<>* getCrescendo25D(const std::string& inputFile,
                                const real_t height)
{
    gsMultiPatch<>* mp = readMultipatch(inputFile);
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
    boundary->addPatch(geom0);
    boundary->addPatch(geom1);
    boundary->addPatch(geom2);
    boundary->addPatch(geom3);

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
gsMultiPatch<>* getCrescendoParaSplittedTopBottom(const std::vector< gsTensorBSpline<2, real_t> > cpatches)
{
    gsTensorBSpline<2, real_t> cpatch0Init = cpatches.at(0);
    gsTensorBSpline<2, real_t> cpatch1Init = cpatches.at(1);
    gsTensorBSpline<2, real_t> cpatch2Init = cpatches.at(2);

    gsMatrix<> coefs0Init = cpatch0Init.coefs();
    gsMatrix<> coefs1Init = cpatch1Init.coefs();
    gsMatrix<> coefs2Init = cpatch2Init.coefs();

    index_t nCoefs0 = coefs0Init.rows();
    index_t nCoefs1 = coefs1Init.rows();
    index_t nCoefs2 = coefs2Init.rows();

    gsMatrix<> coefs0(nCoefs0, 3);
    coefs0.setZero();
    coefs0.block(0, 0, nCoefs0, 2) = coefs0Init.block(0, 0, nCoefs0, 2);

    gsMatrix<> coefs1(nCoefs1, 3);
    coefs1.setZero();
    coefs1.block(0, 0, nCoefs1, 2) = coefs1Init.block(0, 0, nCoefs1, 2);

    gsMatrix<> coefs2(nCoefs2, 3);
    coefs2.setZero();
    coefs2.block(0, 0, nCoefs2, 2) = coefs2Init.block(0, 0, nCoefs2, 2);

    gsKnotVector<> kv01 = cpatch0Init.knots(0);
    gsKnotVector<> kv02 = cpatch0Init.knots(1);
    gsTensorBSpline<2> cpBottom0(kv01, kv02, coefs0);

    gsKnotVector<> kv11 = cpatch1Init.knots(0);
    gsKnotVector<> kv12 = cpatch1Init.knots(1);
    gsTensorBSpline<2> cpBottom1(kv11, kv12, coefs1);

    gsKnotVector<> kv21 = cpatch2Init.knots(0);
    gsKnotVector<> kv22 = cpatch2Init.knots(1);
    gsTensorBSpline<2> cpBottom2(kv21, kv22, coefs2);

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
    gsTensorBSpline<2> cpTop0(kv01, kv02, coefs0);
    gsTensorBSpline<2> cpTop1(kv11, kv12, coefs1);
    gsTensorBSpline<2> cpTop2(kv21, kv22, coefs2);
    mp->addPatch(cpTop0);
    mp->addPatch(cpTop1);
    mp->addPatch(cpTop2);

    writeMultipatch(*mp, "crescendo-tmtf-para-mp_patch4_bottom_top");

    return mp;
}

/**************************************************/
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
    
    /********************************************************/
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
    
    /********************************************************/
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

// This function prints elements of std::vector
void printVector(std::vector<real_t> vec)
{
    for(unsigned i = 0;i<vec.size();++i)
        std::cout << vec[i] <<" ";
    std::cout << std::endl;
}

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
