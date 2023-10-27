/** @file gsSRM46RotorUtils.h

    @brief This file contains various auxilary routines used with
    SRM4+6 Rotor demonstrator geometry from the MOTOR project
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

// Returns prescribed parameters of the centers of
// lobes in the female rotor geometry or corresponding
// points in the template geometry
gsMatrix<> FemaleRotParams(const real_t parDelta);

// Returns the connectivity matrix for the segmentation
// of Template 1 (female rotor)
gsMatrix<> FemaleRotTemp1ConnectivityMatrix();

// Returns std::vector of boundaries which are used for Coon's patches construction
// on the segmented template geometry
std::vector< gsMultiPatch<> > FemaleRotTemp1Boundaries(const gsMultiPatch<>& mpTemp,
                                                       const gsMultiPatch<>& mpSegm);

// Splits the boundary of the female rotor target or template geometries at the points
// corresponding to given parameters
gsMultiPatch<>* FemaleRotSplittedBoundary(const std::string& inputFile,
                                          const std::string& outputFile,
                                          const gsMatrix<>& params);

// Returns coordinates of four rectangle vertices
// ptBound - point on the boundary (typically, the center of a lobe)
// dist - distance from the template center ptCenter to the rectangle
// l1, l2 - lengths of rectangle sides
gsMatrix<> FemaleRotTemp1RectanglePoints(const gsMatrix<>& ptCenter,
                                         const gsMatrix<>& ptBound,
                                         const real_t dist,
                                         const real_t l1,
                                         const real_t l2);

// Returns the coordinates of the female rotor geometry or template center
gsMatrix<> FemaleRotCenter(const gsGeometry<>* geom,
                           const real_t parDelta);

// Returns the parameters of segmentation points on the boundary
gsMatrix<> FemaleRotTemp1BoundaryParameters(const real_t parDelta,
                                            const real_t parDelta1,
                                            const real_t parDelta2);

// Returns the segmentation points on the boundary
gsMatrix<> FemaleRotTemp1BoundaryPoints(const gsGeometry<>* geom,
                                        const double parDelta,
                                        const double parDelta1,
                                        const double parDelta2);

//
void getFemaleRotBoundaryConstraints(const std::string& inputBoundaryFile,
                                   const std::string& inputCoonsFile,
                                   const std::string& inputParamsPointsFile,
                                   const std::string& outputFile);

// Splits the boundary of the Template 1 (male rotor) at the points
// corresponding to given parameters
gsMultiPatch<>* MaleRotTemp1SplittedBoundary(const std::string& inputFile,
                                             const std::string& outputFile,
                                             const gsMatrix<>& params);

// ...
gsMatrix<> MaleRotParams();

// ...
gsMatrix<> MaleRotTemp1ConnectivityMatrix();

// Returns the parameters of segmentation points on the boundary
gsMatrix<> MaleRotTemp1BoundaryParameters();

// Returns the segmentation points on the boundary
gsMatrix<> MaleRotTemp1BoundaryPoints(const gsGeometry<>* geom);

// ...
gsMultiPatch<>* MaleRotTemp1SplittedBoundary(const std::string& tempFile,
                                             const gsMatrix<>& params);

//
void getMaleRotBoundaryConstraints(const std::string& inputBoundaryFile,
                                   const std::string& inputCoonsFile,
                                   const std::string& inputParamsPointsFile,
                                   const std::string& outputFile);

// Returns the coordinates of intersection point of two
// lines which go through points with coordinates stored
// in matrices p1 and p2
gsMatrix<> intersectionOfLines(const gsMatrix<>& p1,
                               const gsMatrix<>& p2);

// On the line y = slope*x + b finds those two points which are
// located within the distance dist from the point pt
gsMatrix<> equidistantPoints(const gsMatrix<>& pt,
                             const real_t slope,
                             const real_t b,
                             const real_t dist);

/**************************************************/
gsMatrix<> FemaleRotParams(const real_t parDelta)
{
    gsMatrix<> params(1, 6);
    params << 1 + parDelta, 1.0/6.0 + parDelta, 2.0/6.0 + parDelta,
              3.0/6.0 + parDelta, 4.0/6.0 + parDelta, 5.0/6.0 + parDelta;

    return params;
}

/**************************************************/
gsMatrix<> FemaleRotTemp1ConnectivityMatrix()
{
    const index_t numLines = 58;
    gsMatrix<> connect(numLines, 2);
    // x, y = point No. x should be connected with point No. y
    connect <<  0,  1, // Line 0
                2,  3, // Line 1
                0,  2, // Line 2
                1,  3, // Line 3
                4,  5, // Line 4
                6,  7, // Line 5
                4,  6, // Line 6
                5,  7, // Line 7
                8,  9, // Line 8
               10, 11, // Line 9
                8, 10, // Line 10
                9, 11, // Line 11
               12, 13, // Line 12
               14, 15, // Line 13
               12, 14, // Line 14
               13, 15, // Line 15
               16, 17, // Line 16
               18, 19, // Line 17
               16, 18, // Line 18
               17, 19, // Line 19
               20, 21, // Line 20
               22, 23, // Line 21
               20, 22, // Line 22
               21, 23, // Line 23
                0, 24, // Line 24
                2, 25, // Line 25
                3, 26, // Line 26
                1, 27, // Line 27
                5, 28, // Line 28
                7, 29, // Line 29
                6, 30, // Line 30
                4, 31, // Line 31
                9, 32, // Line 32
               11, 33, // Line 33
               10, 34, // Line 34
                8, 35, // Line 35
               13, 36, // Line 36
               15, 37, // Line 37
               14, 38, // Line 38
               12, 39, // Line 39
               16, 40, // Line 40
               18, 41, // Line 41
               19, 42, // Line 42
               17, 43, // Line 43
               20, 44, // Line 44
               22, 45, // Line 45
               23, 46, // Line 46
               21, 47, // Line 47
                1,  5, // Line 48
                4,  9, // Line 49
                8, 13, // Line 50
               12, 16, // Line 51
               17, 20, // Line 52
               21,  0, // Line 53
                0,  4, // Line 54
               21,  9, // Line 55
               20,  8, // Line 56
               17, 13; // Line 57

    return connect.transpose();

}

/**************************************************/
std::vector< gsMultiPatch<> > FemaleRotTemp1Boundaries(const gsMultiPatch<>& mpTemp,
                                                       const gsMultiPatch<>& mpSegm)
{
    std::vector< gsMultiPatch<> > boundaries;
    gsMultiPatch<> boundary0;
    gsMultiPatch<> boundary1;
    gsMultiPatch<> boundary2;
    gsMultiPatch<> boundary3;
    gsMultiPatch<> boundary4;
    gsMultiPatch<> boundary5;
    gsMultiPatch<> boundary6;
    gsMultiPatch<> boundary7;
    gsMultiPatch<> boundary8;
    gsMultiPatch<> boundary9;
    gsMultiPatch<> boundary10;
    gsMultiPatch<> boundary11;
    gsMultiPatch<> boundary12;
    gsMultiPatch<> boundary13;
    gsMultiPatch<> boundary14;
    gsMultiPatch<> boundary15;
    gsMultiPatch<> boundary16;
    gsMultiPatch<> boundary17;
    gsMultiPatch<> boundary18;
    gsMultiPatch<> boundary19;
    gsMultiPatch<> boundary20;
    gsMultiPatch<> boundary21;
    gsMultiPatch<> boundary22;
    gsMultiPatch<> boundary23;
    gsMultiPatch<> boundary24;
    gsMultiPatch<> boundary25;
    gsMultiPatch<> boundary26;
    gsMultiPatch<> boundary27;
    gsMultiPatch<> boundary28;
    gsMultiPatch<> boundary29;
    gsMultiPatch<> boundary30;
    gsMultiPatch<> boundary31;
    gsMultiPatch<> boundary32;
    gsMultiPatch<> boundary33;
    gsMultiPatch<> boundary34;

    // The lines ordering in boundary patches is important
    // for unifyKnotsBoundary function.
    // Please pay attention to this or improve unifyKnotsBoundary
    // Example of allowed ordering:
    //        1
    //     +-----+
    //    ¦       ¦
    //  0 ¦       ¦ 2
    //    ¦       ¦
    //     +-----+
    //        3
    // Example of not allowed ordering:
    //        2
    //     +-----+
    //    ¦       ¦
    //  0 ¦       ¦ 1
    //    ¦       ¦
    //     +-----+
    //        3


    // Patches which include boundary edges should be
    // added in proper order

    boundary0.addPatch(mpSegm[25]);
    boundary0.addPatch(mpSegm[1]);
    boundary0.addPatch(mpSegm[26]);
    boundary0.addPatch(mpTemp[0]);
    boundaries.push_back(boundary0);

    boundary1.addPatch(mpSegm[26]);
    boundary1.addPatch(mpSegm[3]);
    boundary1.addPatch(mpSegm[27]);
    boundary1.addPatch(mpTemp[1]);
    boundaries.push_back(boundary1);

    boundary2.addPatch(mpSegm[27]);
    boundary2.addPatch(mpSegm[48]);
    boundary2.addPatch(mpSegm[28]);
    boundary2.addPatch(mpTemp[2]);
    boundaries.push_back(boundary2);

    boundary3.addPatch(mpSegm[28]);
    boundary3.addPatch(mpSegm[7]);
    boundary3.addPatch(mpSegm[29]);
    boundary3.addPatch(mpTemp[3]);
    boundaries.push_back(boundary3);

    boundary4.addPatch(mpSegm[29]);
    boundary4.addPatch(mpSegm[5]);
    boundary4.addPatch(mpSegm[30]);
    boundary4.addPatch(mpTemp[4]);
    boundaries.push_back(boundary4);

    boundary5.addPatch(mpSegm[30]);
    boundary5.addPatch(mpSegm[6]);
    boundary5.addPatch(mpSegm[31]);
    boundary5.addPatch(mpTemp[5]);
    boundaries.push_back(boundary5);

    boundary6.addPatch(mpSegm[31]);
    boundary6.addPatch(mpSegm[49]);
    boundary6.addPatch(mpSegm[32]);
    boundary6.addPatch(mpTemp[6]);
    boundaries.push_back(boundary6);

    boundary7.addPatch(mpSegm[32]);
    boundary7.addPatch(mpSegm[11]);
    boundary7.addPatch(mpSegm[33]);
    boundary7.addPatch(mpTemp[7]);
    boundaries.push_back(boundary7);

    boundary8.addPatch(mpSegm[33]);
    boundary8.addPatch(mpSegm[9]);
    boundary8.addPatch(mpSegm[34]);
    boundary8.addPatch(mpTemp[8]);
    boundaries.push_back(boundary8);

    boundary9.addPatch(mpSegm[34]);
    boundary9.addPatch(mpSegm[10]);
    boundary9.addPatch(mpSegm[35]);
    boundary9.addPatch(mpTemp[9]);
    boundaries.push_back(boundary9);

    boundary10.addPatch(mpSegm[35]);
    boundary10.addPatch(mpSegm[50]);
    boundary10.addPatch(mpSegm[36]);
    boundary10.addPatch(mpTemp[10]);
    boundaries.push_back(boundary10);

    boundary11.addPatch(mpSegm[36]);
    boundary11.addPatch(mpSegm[15]);
    boundary11.addPatch(mpSegm[37]);
    boundary11.addPatch(mpTemp[11]);
    boundaries.push_back(boundary11);

    boundary12.addPatch(mpSegm[37]);
    boundary12.addPatch(mpSegm[13]);
    boundary12.addPatch(mpSegm[38]);
    boundary12.addPatch(mpTemp[12]);
    boundaries.push_back(boundary12);

    boundary13.addPatch(mpSegm[38]);
    boundary13.addPatch(mpSegm[14]);
    boundary13.addPatch(mpSegm[39]);
    boundary13.addPatch(mpTemp[13]);
    boundaries.push_back(boundary13);

    boundary14.addPatch(mpSegm[39]);
    boundary14.addPatch(mpSegm[51]);
    boundary14.addPatch(mpSegm[40]);
    boundary14.addPatch(mpTemp[14]);
    boundaries.push_back(boundary14);

    boundary15.addPatch(mpSegm[40]);
    boundary15.addPatch(mpSegm[18]);
    boundary15.addPatch(mpSegm[41]);
    boundary15.addPatch(mpTemp[15]);
    boundaries.push_back(boundary15);

    boundary16.addPatch(mpSegm[41]);
    boundary16.addPatch(mpSegm[17]);
    boundary16.addPatch(mpSegm[42]);
    boundary16.addPatch(mpTemp[16]);
    boundaries.push_back(boundary16);

    boundary17.addPatch(mpSegm[42]);
    boundary17.addPatch(mpSegm[19]);
    boundary17.addPatch(mpSegm[43]);
    boundary17.addPatch(mpTemp[17]);
    boundaries.push_back(boundary17);

    boundary18.addPatch(mpSegm[43]);
    boundary18.addPatch(mpSegm[52]);
    boundary18.addPatch(mpSegm[44]);
    boundary18.addPatch(mpTemp[18]);
    boundaries.push_back(boundary18);

    boundary19.addPatch(mpSegm[44]);
    boundary19.addPatch(mpSegm[22]);
    boundary19.addPatch(mpSegm[45]);
    boundary19.addPatch(mpTemp[19]);
    boundaries.push_back(boundary19);

    boundary20.addPatch(mpSegm[45]);
    boundary20.addPatch(mpSegm[21]);
    boundary20.addPatch(mpSegm[46]);
    boundary20.addPatch(mpTemp[20]);
    boundaries.push_back(boundary20);

    boundary21.addPatch(mpSegm[46]);
    boundary21.addPatch(mpSegm[23]);
    boundary21.addPatch(mpSegm[47]);
    boundary21.addPatch(mpTemp[21]);
    boundaries.push_back(boundary21);

    boundary22.addPatch(mpSegm[47]);
    boundary22.addPatch(mpSegm[53]);
    boundary22.addPatch(mpSegm[24]);
    boundary22.addPatch(mpTemp[22]);
    boundaries.push_back(boundary22);

    boundary23.addPatch(mpSegm[24]);
    boundary23.addPatch(mpSegm[2]);
    boundary23.addPatch(mpSegm[25]);
    boundary23.addPatch(mpTemp[23]);
    boundaries.push_back(boundary23);

    // Patches without boundary edges
    boundary24.addPatch(mpSegm[0]);
    boundary24.addPatch(mpSegm[1]);
    boundary24.addPatch(mpSegm[2]);
    boundary24.addPatch(mpSegm[3]);
    boundaries.push_back(boundary24);

    boundary25.addPatch(mpSegm[4]);
    boundary25.addPatch(mpSegm[5]);
    boundary25.addPatch(mpSegm[6]);
    boundary25.addPatch(mpSegm[7]);
    boundaries.push_back(boundary25);

    boundary26.addPatch(mpSegm[8]);
    boundary26.addPatch(mpSegm[9]);
    boundary26.addPatch(mpSegm[10]);
    boundary26.addPatch(mpSegm[11]);
    boundaries.push_back(boundary26);

    boundary27.addPatch(mpSegm[12]);
    boundary27.addPatch(mpSegm[13]);
    boundary27.addPatch(mpSegm[15]);
    boundary27.addPatch(mpSegm[14]);
    boundaries.push_back(boundary27);

    boundary28.addPatch(mpSegm[16]);
    boundary28.addPatch(mpSegm[17]);
    boundary28.addPatch(mpSegm[18]);
    boundary28.addPatch(mpSegm[19]);
    boundaries.push_back(boundary28);

    boundary29.addPatch(mpSegm[20]);
    boundary29.addPatch(mpSegm[21]);
    boundary29.addPatch(mpSegm[22]);
    boundary29.addPatch(mpSegm[23]);
    boundaries.push_back(boundary29);

    boundary30.addPatch(mpSegm[4]);
    boundary30.addPatch(mpSegm[0]);
    boundary30.addPatch(mpSegm[48]);
    boundary30.addPatch(mpSegm[54]);
    boundaries.push_back(boundary30);

    boundary31.addPatch(mpSegm[49]);
    boundary31.addPatch(mpSegm[53]);
    boundary31.addPatch(mpSegm[54]);
    boundary31.addPatch(mpSegm[55]);
    boundaries.push_back(boundary31);

    boundary32.addPatch(mpSegm[8]);
    boundary32.addPatch(mpSegm[20]);
    boundary32.addPatch(mpSegm[55]);
    boundary32.addPatch(mpSegm[56]);
    boundaries.push_back(boundary32);

    boundary33.addPatch(mpSegm[50]);
    boundary33.addPatch(mpSegm[52]);
    boundary33.addPatch(mpSegm[56]);
    boundary33.addPatch(mpSegm[57]);
    boundaries.push_back(boundary33);

    boundary34.addPatch(mpSegm[12]);
    boundary34.addPatch(mpSegm[16]);
    boundary34.addPatch(mpSegm[57]);
    boundary34.addPatch(mpSegm[51]);
    boundaries.push_back(boundary34);

    return boundaries;

}

/**************************************************/
gsMultiPatch<>* FemaleRotSplittedBoundary(const std::string& inputFile,
                                          const std::string& outputFile,
                                          const gsMatrix<>& params)
{

    gsGeometry<>* geom = readGeometry(inputFile);

    // Because of splittedCurveParameters function we need to rearrange points
    gsMatrix<> params2(params.rows(), params.cols());
    params2.setZero();
    params2.block(0, 0, params.rows(), params.cols()-2) = params.block(0, 2, params.rows(), params.cols()-2);
    params2.block(0, params.cols()-2, params.rows(), 2) = params.block(0, 0, params.rows(), 2);

    gsMultiPatch<> mpTmp = splittedCurveParameters(*geom, params2);
    writeMultipatch(mpTmp, outputFile);

    gsMultiPatch<>* mp = &mpTmp;

    return mp;
}

/**************************************************/
gsMatrix<> FemaleRotTemp1RectanglePoints(const gsMatrix<> &ptCenter,
                                         const gsMatrix<> &ptBound,
                                         const real_t dist,
                                         const real_t l1,
                                         const real_t l2)
{
    gsMatrix<> pts(2, 4);
    pts.setZero();

    // Length of the segment between the center and boundary point
    const real_t lineLen = (ptBound-ptCenter).norm();

    // Slope of the normal to the line going through the center and boundary point
    const real_t nrmlSlope = -(ptBound(0, 0)-ptCenter(0, 0))/(ptBound(1, 0)-ptCenter(1, 0));

    real_t parRectPt = dist/lineLen;
    gsMatrix<> rectPt = parRectPt*ptBound +  (1-parRectPt)*ptCenter;
    real_t nrmlB = -nrmlSlope*rectPt(0, 0)+rectPt(1, 0);
    gsMatrix<> ptsPair = equidistantPoints(rectPt, nrmlSlope, nrmlB, l1/2.0);
    pts.block(0, 0, 2, 2) = ptsPair;

    // Do the same with different parameter value
    parRectPt = (dist+l2)/lineLen;
    rectPt = parRectPt*ptBound +  (1-parRectPt)*ptCenter;
    nrmlB = -nrmlSlope*rectPt(0, 0)+rectPt(1, 0);
    ptsPair = equidistantPoints(rectPt, nrmlSlope, nrmlB, l1/2.0);
    pts.block(0, 2, 2, 2) = ptsPair;

    return pts;
}

/**************************************************/
gsMatrix<> FemaleRotCenter(const gsGeometry<>* geom,
                           const real_t parDelta)
{

   gsMatrix<> params = FemaleRotParams(parDelta);

   // Centers of the lobes
   const gsMatrix<> pts = geom->eval(params);

   gsMatrix<> pts1(2, 2);
   pts1.col(0) = pts.col(0);
   pts1.col(1) = pts.col(3);
   gsMatrix<> pts2(2, 2);
   pts2.col(0) = pts.col(1);
   pts2.col(1) = pts.col(4);

   // The center is assumed to be the intersection point of
   // two lines going through centers of lobes 0, 3 and 1, 4
   gsMatrix<> c = intersectionOfLines(pts1, pts2);
   gsInfo << "Female rotor center:\n" << c << std::endl;

   return c;
}

/**************************************************/
gsMatrix<> FemaleRotTemp1BoundaryParameters(const real_t parDelta,
                                            const real_t parDelta1,
                                            const real_t parDelta2)
{
    gsMatrix<> params = FemaleRotParams(parDelta);

    gsMatrix<> paramsBound(1, 24);
    paramsBound.setZero();

    paramsBound << params(0, 0)-parDelta2, params(0, 0)-parDelta1, params(0, 0)+parDelta1-1.0, params(0, 0)+parDelta2-1.0,
                   params(0, 1)-parDelta2, params(0, 1)-parDelta1, params(0, 1)+parDelta1, params(0, 1)+parDelta2,
                   params(0, 2)-parDelta2, params(0, 2)-parDelta1, params(0, 2)+parDelta1, params(0, 2)+parDelta2,
                   params(0, 3)-parDelta2, params(0, 3)-parDelta1, params(0, 3)+parDelta1, params(0, 3)+parDelta2,
                   params(0, 4)-parDelta2, params(0, 4)-parDelta1, params(0, 4)+parDelta1, params(0, 4)+parDelta2,
                   params(0, 5)-parDelta2, params(0, 5)-parDelta1, params(0, 5)+parDelta1, params(0, 5)+parDelta2;

    return paramsBound;
}

/**************************************************/
gsMatrix<> FemaleRotTemp1BoundaryPoints(const gsGeometry<>* geom,
                                        const real_t parDelta,
                                        const real_t parDelta1,
                                        const real_t parDelta2)
{

    gsMatrix<> paramsBound = FemaleRotTemp1BoundaryParameters(parDelta, parDelta1, parDelta2);
    const gsMatrix<> pts = geom->eval(paramsBound);

    return pts;
}

/**************************************************/
/**************************************************/
/**************************************************/

/**************************************************/
gsMatrix<> MaleRotParams()
{
    gsMatrix<> params(1, 4);
    const real_t parDelta = -0.001;
    params << 1.0 + parDelta, 1.0/4.0 + parDelta, 2.0/4.0 + parDelta, 3.0/4.0 + parDelta;

    return params;
}

/**************************************************/
gsMatrix<> MaleRotTemp1ConnectivityMatrix()
{
    const index_t numLines = 8;
    gsMatrix<> connect(numLines, 2);
    // x, y = point No. x should be connected with point No. y
    connect <<  0,  7,
                1,  4,
                2,  5,
                3,  6,
                0,  1,
                1,  2,
                3,  2,
                0,  3;

    return connect.transpose();

}

/**************************************************/
std::vector< gsMultiPatch<> > MaleRotTemp1Boundaries(const gsMultiPatch<>& mpTemp,
                                                     const gsMultiPatch<>& mpSegm)
{
    std::vector< gsMultiPatch<> > boundaries;
    gsMultiPatch<> boundary0;
    gsMultiPatch<> boundary1;
    gsMultiPatch<> boundary2;
    gsMultiPatch<> boundary3;
    gsMultiPatch<> boundary4;

    // The lines ordering in boundary patches is important
    // for unifyKnotsBoundary function.
    // Please pay attention to this or improve unifyKnotsBoundary
    // Example of allowed ordering:
    //        1
    //     +-----+
    //    ¦       ¦
    //  0 ¦       ¦ 2
    //    ¦       ¦
    //     +-----+
    //        3
    // Example of not allowed ordering:
    //        2
    //     +-----+
    //    ¦       ¦
    //  0 ¦       ¦ 1
    //    ¦       ¦
    //     +-----+
    //        3

    boundary0.addPatch(mpSegm[0]);
    boundary0.addPatch(mpSegm[4]);
    boundary0.addPatch(mpSegm[1]);
    boundary0.addPatch(mpTemp[0]);
    boundaries.push_back(boundary0);

    boundary1.addPatch(mpSegm[1]);
    boundary1.addPatch(mpSegm[5]);
    boundary1.addPatch(mpSegm[2]);
    boundary1.addPatch(mpTemp[1]);
    boundaries.push_back(boundary1);

    boundary2.addPatch(mpSegm[2]);
    boundary2.addPatch(mpSegm[6]);
    boundary2.addPatch(mpSegm[3]);
    boundary2.addPatch(mpTemp[2]);
    boundaries.push_back(boundary2);

    boundary3.addPatch(mpSegm[3]);
    boundary3.addPatch(mpSegm[7]);
    boundary3.addPatch(mpSegm[0]);
    boundary3.addPatch(mpTemp[3]);
    boundaries.push_back(boundary3);

    boundary4.addPatch(mpSegm[4]);
    boundary4.addPatch(mpSegm[5]);
    boundary4.addPatch(mpSegm[6]);
    boundary4.addPatch(mpSegm[7]);
    boundaries.push_back(boundary4);

    writeMultipatch(boundary0, "test_boundary0");
    writeMultipatch(boundary1, "test_boundary1");
    writeMultipatch(boundary2, "test_boundary2");
    writeMultipatch(boundary3, "test_boundary3");

    return boundaries;

}

/**************************************************/
gsMultiPatch<>* MaleRotTemp1SplittedBoundary(const std::string& inputFile,
                                             const std::string& outputFile,
                                             const gsMatrix<>& params)
{

    gsGeometry<>* geom = readGeometry(inputFile);

    gsMultiPatch<> mpTmp = splittedCurveParameters(*geom, params);
    writeMultipatch(mpTmp, outputFile);

    gsMultiPatch<>* mp = &mpTmp;

    return mp;
}

/**************************************************/
gsMatrix<> MaleRotTemp1BoundaryParameters()
{
    gsMatrix<> params = MaleRotParams();

    gsMatrix<> paramsBound(1, 4);
    paramsBound.setZero();
    paramsBound << params(0, 1), params(0, 2), params(0, 3), params(0, 0);

    return paramsBound;
}

/**************************************************/
gsMatrix<> MaleRotTemp1BoundaryPoints(const gsGeometry<>* geom)
{

    gsMatrix<> paramsBound = MaleRotTemp1BoundaryParameters();
    const gsMatrix<> pts = geom->eval(paramsBound);

    return pts;
}


/**************************************************/
gsMatrix<> intersectionOfLines(const gsMatrix<>& p1,
                               const gsMatrix<>& p2)
{
    gsMatrix<> p(2, 1);

    const real_t x11 = p1(0, 0);
    const real_t y11 = p1(1, 0);
    const real_t x12 = p1(0, 1);
    const real_t y12 = p1(1, 1);

    const real_t x21 = p2(0, 0);
    const real_t y21 = p2(1, 0);
    const real_t x22 = p2(0, 1);
    const real_t y22 = p2(1, 1);

    // Expressions have been obtained using Mathematica
    p(0, 0) = (x11*(x22*(-y12 + y21) + x21*(y12 - y22)) + x12*(x22*(y11 - y21) + x21*(-y11 + y22))) / (x22*(y11 - y12) + x21*(-y11 + y12) + (x11 - x12)*(y21 - y22));
    p(1, 0) = (x22*(y11 - y12)*y21 + x11*y12*y21 - x21*y11*y22 - x11*y12*y22 + x21*y12*y22 + x12*y11*(-y21 + y22)) / (x22*(y11 - y12) + x21*(-y11 + y12) + (x11 - x12)*(y21 - y22));

    return p;
}

/**************************************************/
gsMatrix<> equidistantPoints(const gsMatrix<> &pt,
                             const real_t slope,
                             const real_t b,
                             const real_t dist)
{
    gsMatrix<> pts(2, 2);
    pts.setZero();

    const real_t x = pt(0, 0);
    const real_t y = pt(1, 0);

    // Expressions have been obtained using Mathematica
    // as result of solving equations
    // (x0 - x)^2 + (y0 - y)^2 = dist^2
    //  y0 = slope*x0 + b
    const real_t D = math::sqrt(-b*b+dist*dist*(1+slope*slope)-(-slope*x + y)*(-slope*x + y)+b*(-2*slope*x+2*y));

    real_t x0 = (-b*slope+x+slope*y - D)/(1 + slope*slope);
    real_t y0 = slope*x0 + b;
    pts(0, 0) = x0;
    pts(1, 0) = y0;

    x0 = (-b*slope+x+slope*y + D)/(1 + slope*slope);
    y0 = slope*x0 + b;
    pts(0, 1) = x0;
    pts(1, 1) = y0;

    return pts;
}

} // namespace gismo
