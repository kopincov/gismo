/** @file gsAirfoilUtils.h

    @brief This file contains functions used in gsAirfoilTest.cpp

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Sajavicius
*/

#pragma once

#include "gsMotorIOUtils.h"

namespace gismo
{

/**************************************************/
// Converts string code of the template to integer.
// Function needs to be updated each time a new
// template is introduced.
int tempCaseAsInt(const std::string tempCase)
{
    int tempCaseInt = 0;
    if (tempCase == "In5patches")
    {
        tempCaseInt = 1;
    }
    else if (tempCase == "In11patches")
    {
        tempCaseInt = 2;
    }
    else if (tempCase == "Out4patches")
    {
        tempCaseInt = 3;
    }
    else if (tempCase == "Out6patches")
    {
        tempCaseInt = 4;
    }
    else if (tempCase == "Out7patches")
    {
        tempCaseInt = 5;
    }
    else if (tempCase == "Out8patches")
    {
        tempCaseInt = 6;
    }
    else if (tempCase == "UniKl22patches")
    {
        tempCaseInt = 10;
    }
    return tempCaseInt;
}

/**************************************************/
// Returns number of splitting points
int numSplitPoints(const int tempCase)
{
    int numSplitPts = 0;
    switch(tempCase)
    {
        case 1: // "In5patches"
        {
            numSplitPts  = 4;
            break;
        }
        case 2: // "In11patches"
        {
            numSplitPts  = 8;
            break;
        }
        case 3: // "Out4patches"
        {
            numSplitPts  = 8;
            break;
        }
        case 4: // "Out6patches"
        {
            numSplitPts  = 12;
            break;
        }
        case 5: // "Out7patches"
        {
            numSplitPts  = 14;
            break;
        }
        case 6: // "Out8patches"
        {
            numSplitPts  = 16;
            break;
        }
        case 10: // "UniKl22patches"
        {
            numSplitPts  = 24;
            break;
        }
        default:
        {
            gsInfo << "Invalid template case" << std::endl;
        }
    };

    return numSplitPts;
}

/**************************************************/
// Returns parameters corresponding to the points
// splitting the boundary curve (ellipse)
gsMatrix<> splittingParameters(const int tempCase)
{
    gsMatrix<> result;

    // PARAMETER SHOULD NOT COINCIDE WITH KNOT VALUES
    // Otherwise, error appears
    switch(tempCase)
    {
        case 1: // "In5patches"
        {
            gsMatrix<> splittingParams(1, 4);
            splittingParams << 0.0, 0.24875, 0.5, 0.75125;
            result = splittingParams;
            break;
        }
        case 2: // "In11patches"
        {
            gsMatrix<> splittingParams(1, 8);
            splittingParams << 0.0, 0.1125, 0.1575, 0.34, 0.5, 0.66, 0.8425, 0.887;
            result = splittingParams;
            break;
        }
        case 3: // "Out4patches"
        {
            gsMatrix<> splittingParams(1, 3);
            splittingParams << 0.0, 0.4, 0.6;
            result = splittingParams;
            break;
        }
        case 4: // "Out6patches"
        {
            gsMatrix<> splittingParams(1, 2);
            splittingParams << 0.0, 0.5;
            result = splittingParams;
            break;
        }
        case 5: // "Out7patches"
        {
            gsMatrix<> splittingParams(1, 3);
            splittingParams << 0.0, 0.4, 0.6;
            result = splittingParams;
            break;
        }
        case 6: // "Out8patches"
        {
            gsMatrix<> splittingParams(1, 4);
            splittingParams << 0.0, 0.25, 0.5, 0.75;
            result = splittingParams;
            break;
        }
        case 10: // "UniKl22patches"
        {
            gsMatrix<> splittingParams(1, 10);
            splittingParams << 0.0, 0.15, 0.24, 0.30, 0.35, 0.49, 0.65, 0.70, 0.76, 0.85;
            result = splittingParams;
            break;
        }
        default:
        {

            gsInfo << "Invalid template case" << std::endl;
        }
    };

    return result;
}

/**************************************************/
// Returns matrix of splitting points
gsMatrix<> splittingPoints(const gsMultiPatch<>& temp,
                           const int tempCase)
{
    const int numSplitPts = numSplitPoints(tempCase);

    gsMatrix<> splitPts(2, numSplitPts);
    splitPts.setZero();

    switch(tempCase)
    {
        case 1: // "In5patches"
        {
            break;
        }
        case 2: // "In11patches"
        {
            break;
        }
        case 3: // "Out4patches"
        {
            break;
        }
        case 4: // "Out6patches"
        {
            break;
        }
        case 5: // "Out7patches"
        {
            break;
        }
        case 6: // "Out8patches"
        {
            break;
        }
        case 10: // "UniKl22patches"
        {
            const gsGeometry<>& tempBound1 = temp[1];
            gsMatrix<> splitParamsEllipse = splittingParameters(tempCase);
            gsMatrix<> splitPtsEllipse = tempBound1.eval(splitParamsEllipse);

            gsMatrix<> splitPtsRectangle(2, numSplitPts - splitPtsEllipse.cols());

            splitPtsRectangle(0, 0) = 2.0;
            splitPtsRectangle(1, 0) = splitPtsEllipse(1, 0);

            splitPtsRectangle(0, 1) = 2.0;
            splitPtsRectangle(1, 1) = 0.8;

            splitPtsRectangle(0, 2) = 0.5*(splitPtsEllipse(0, 0)+2.0);
            splitPtsRectangle(1, 2) = 0.8;

            splitPtsRectangle(0, 3) = splitPtsEllipse(0, 1);
            splitPtsRectangle(1, 3) = 0.8;

            splitPtsRectangle(0, 4) = splitPtsEllipse(0, 2);
            splitPtsRectangle(1, 4) = 0.8;

            splitPtsRectangle(0, 5) = splitPtsEllipse(0, 3);
            splitPtsRectangle(1, 5) = 0.8;

            splitPtsRectangle(0, 6) = -1.0;
            splitPtsRectangle(1, 6) =  0.8;

            splitPtsRectangle(0, 7) = -1.0;
            splitPtsRectangle(1, 7) = splitPtsEllipse(1, 5);

            splitPtsRectangle(0, 8) = -1.0;
            splitPtsRectangle(1, 8) = -0.8;

            splitPtsRectangle(0, 9) = splitPtsEllipse(0, 7);
            splitPtsRectangle(1, 9) = -0.8;

            splitPtsRectangle(0, 10) = splitPtsEllipse(0, 8);
            splitPtsRectangle(1, 10) = -0.8;

            splitPtsRectangle(0, 11) = splitPtsEllipse(0, 9);
            splitPtsRectangle(1, 11) = -0.8;

            splitPtsRectangle(0, 12) = 0.5*(splitPtsEllipse(0, 0)+2.0);
            splitPtsRectangle(1, 12) = -0.8;

            splitPtsRectangle(0, 13) = 2.0;
            splitPtsRectangle(1, 13) = -0.8;

            splitPts.block(0, 0, 2, splitPtsRectangle.cols()) = splitPtsRectangle;
            splitPts.block(0, splitPtsRectangle.cols(), 2, splitPtsEllipse.cols()) = splitPtsEllipse;

            break;
        }
        default:
        {
            gsInfo << "Invalid template case" << std::endl;
        }
    };

    return splitPts;
}

/**************************************************/
// Returns number of segmentation points inside domain
int numSegmPoints(const int tempCase)
{
    int numSegmPts = 0;
    switch(tempCase)
    {
        case 1: // "In5patches"
        {
            numSegmPts  = 4;
            break;
        }
        case 2: // "In11patches"
        {
            numSegmPts  = 8;
            break;
        }
        case 3: // "Out4patches"
        {
            numSegmPts  = 0;
            break;
        }
        case 4: // "Out6patches"
        {
            numSegmPts  = 0;
            break;
        }
        case 5: // "Out7patches"
        {
            numSegmPts  = 0;
            break;
        }
        case 6: // "Out8patches"
        {
            numSegmPts  = 0;
            break;
        }
        case 10: // "UniKl22patches"
        {
            numSegmPts  = 6;
            break;
        }
        default:
        {
            gsInfo << "Invalid template case" << std::endl;
        }
    };

    return numSegmPts;
}

/**************************************************/
// Returns matrix of segmentation points
gsMatrix<> segmentationPoints(const gsMatrix<>& splitPoints,
                              const int tempCase)
{
    const int numSegmPts = numSegmPoints(tempCase);

    gsMatrix<> segmPts(2, numSegmPts);
    segmPts.setZero();

    switch(tempCase)
    {
        case 1: // "In5patches"
        {
            break;
        }
        case 2: // "In11patches"
        {
            break;
        }
        case 3: // "Out4patches"
        {
            break;
        }
        case 4: // "Out6patches"
        {
            break;
        }
        case 5: // "Out7patches"
        {
            break;
        }
        case 6: // "Out8patches"
        {
            break;
        }
        case 10: // "UniKl22patches"
        {
            real_t coord = 0.0;
            if ( splitPoints(0, 18) - math::abs(splitPoints(0, 17)-splitPoints(0, 18)) <= splitPoints(0, 19))
            {
                coord = 0.5*(splitPoints(0, 18)+splitPoints(0, 19));
            }
            else
            {
                coord = splitPoints(0, 18) - math::abs(splitPoints(0, 17)-splitPoints(0, 18));
            }
            segmPts << 0.5*(splitPoints(0, 15)+splitPoints(0, 16)), splitPoints(0, 16), splitPoints(0, 17), splitPoints(0, 18),                         coord,             splitPoints(0, 18),
                       splitPoints(1, 0),                           splitPoints(1, 0),  splitPoints(1, 0),  0.5*(splitPoints(1, 0)+splitPoints(1, 18)), splitPoints(1, 0), 0.5*(splitPoints(1, 0)+splitPoints(1, 20));
            break;
        }
        default:
        {
            gsInfo << "Invalid template case" << std::endl;
        }
    };

    return segmPts;
}

/**************************************************/
// Returns connectivity matrix
gsMatrix<> connectivityMatrix(const int tempCase)
{
    gsMatrix<> result;

    switch(tempCase)
    {
        case 1: // "In5patches"
        {
            const index_t numSegmLines = 2;
            gsMatrix<> connectMatrix(numSegmLines, 2);
            // x, y = point No. x should be connected with point No. y
            connectMatrix <<  0, 0,
                              0, 0;
            result = connectMatrix.transpose();
            break;
        }
        case 2: // "In11patches"
        {
            const index_t numSegmLines = 2;
            gsMatrix<> connectMatrix(numSegmLines, 2);
            // x, y = point No. x should be connected with point No. y
            connectMatrix <<  0, 0,
                              0, 0;
            result = connectMatrix.transpose();
            break;
        }
        case 3: // "Out4patches"
        {
            const index_t numSegmLines = 2;
            gsMatrix<> connectMatrix(numSegmLines, 2);
            // x, y = point No. x should be connected with point No. y
            connectMatrix <<  0, 0,
                              0, 0;
            result = connectMatrix.transpose();
            break;
        }
        case 4: // "Out6patches"
        {
            const index_t numSegmLines = 2;
            gsMatrix<> connectMatrix(numSegmLines, 2);
            // x, y = point No. x should be connected with point No. y
            connectMatrix <<  0, 0,
                              0, 0;
            result = connectMatrix.transpose();
            break;
        }
        case 5: // "Out7patches"
        {
            const index_t numSegmLines = 2;
            gsMatrix<> connectMatrix(numSegmLines, 2);
            // x, y = point No. x should be connected with point No. y
            connectMatrix <<  0, 0,
                              0, 0;
            result = connectMatrix.transpose();
            break;
        }
        case 6: // "Out8patches"
        {
            const index_t numSegmLines = 2;
            gsMatrix<> connectMatrix(numSegmLines, 2);
            // x, y = point No. x should be connected with point No. y
            connectMatrix <<  0, 0,
                              0, 0;
            result = connectMatrix.transpose();
            break;
        }
        case 10: // "UniKl22patches"
        {
            const index_t numSegmLines = 27;
            gsMatrix<> connectMatrix(numSegmLines, 2);
            // x, y = point No. x should be connected with point No. y
            connectMatrix <<  14, 0,
                              14, 2,
                              15, 3,
                              16, 4,
                              17, 5,
                              18, 6,
                              19, 7,
                              20, 8,
                              21, 9,
                              22, 10,
                              23, 11,
                              14, 12,
                              24, 15,
                              24, 23,
                              25, 24,
                              25, 16,
                              25, 22,
                              26, 25,
                              26, 17,
                              26, 21,
                              27, 26,
                              29, 26,
                              27, 18,
                              29, 20,
                              28, 27,
                              28, 29,
                              19, 28;
            result = connectMatrix.transpose();
            break;
        }
        default:
        {
            gsInfo << "Invalid template case" << std::endl;
        }
    };

    return result;
}

/**************************************************/
// Returns multipatch of segmentation lines
gsMultiPatch<>* segmentationLines(const gsMatrix<> pts,
                                  const gsMatrix<> connect)
{
    gsMultiPatch<>* mpSegm = new gsMultiPatch<>();
    const index_t numLines = connect.cols();
    const gsKnotVector<> kv (0, 1, 0, 2); // {0, 0, 1, 1};
    const gsBSplineBasis<> basis( kv );
    for (index_t i = 0; i < numLines; i++)
    {
        gsMatrix<> ctrlPts(2, 2);
        ctrlPts << pts(0, connect(0, i)), pts(1, connect(0, i)),
                   pts(0, connect(1, i)), pts(1, connect(1, i));
        const gsBSpline<> segmLine( basis, ctrlPts );
        mpSegm->addPatch(gsBSpline<>(segmLine));
    }

    return mpSegm;
}

/**************************************************/
// Splits curves at points corresponding to given parameters
gsMultiPatch<> splittedCurveNew(const gsGeometry<>& geom,
                             const gsMatrix<>& pars)
{
    // B-spline curves
    gsBSpline<> inputGeometry = dynamic_cast< const gsBSpline<>& >(geom);
    gsMatrix<> ctrlPts = inputGeometry.coefs();

    int deg = inputGeometry.degree();
    gsKnotVector<> kv = inputGeometry.knots();

    gsMatrix<> points = inputGeometry.eval(pars);

    index_t numPoints = points.cols();

    // Adding knots
    for (int i = 0; i < numPoints; i++)
    {
        if (pars(i) > 0)
        {
            gsBoehm(kv, ctrlPts, pars(i), kv.degree());
        }
    }

    gsMatrix<index_t> indices(2, numPoints);
    indices.setZero();

    indices(0, 0) = 0;
    for (int i = 1; i < numPoints; i++)
    {
        index_t ind = 0;
        for(ind = 0; kv.at(ind) < pars(0, i); ind++); //
        ind--;
        indices(1, i-1)   = ind;
        if (i <= numPoints-1)
        {
            indices(0, i) = ind;
        }
    }
    index_t ind = 0;
    for(ind = 0; kv.at(ind) < 1.0; ind++); //
    ind--;
    indices(1, numPoints-1) = ind;
    // Adding to multipatch all curves
    gsMultiPatch<> mp;
    for (int i = 0; i < numPoints; i++)
    {
        gsMatrix<> ctrlPts3 = ctrlPts.block(indices(0, i), 0, indices(1, i)-indices(0, i)+1, 2);
        gsKnotVector<> kv3(0.0, 1.0, ctrlPts3.rows()-deg-1, deg+1);
        gsBSpline<> extractedGeom(kv3, ctrlPts3);
        mp.addPatch(extractedGeom);
    }
    return mp;
}

/**************************************************/
// Splits template
gsMultiPatch<>* splittedTemplate(const gsMultiPatch<>& temp,
                                 const int tempCase)
{
    gsMultiPatch<>* mpTempSplitted = new gsMultiPatch<>();

    switch(tempCase)
    {
        case 10: // "UniKl22patches"
        {
            // Splitting rectangle
            const gsMatrix<> splitPts = splittingPoints(temp, tempCase);
            const gsKnotVector<> kv (0, 1, 0, 2); // {0, 0, 1, 1};
            const gsBSplineBasis<> basis( kv );
            for (int i = 0; i < 14; i++)
            {
                gsMatrix<> ctrlPts(2, 2);
                if (i < 13)
                {
                    ctrlPts << splitPts(0, i),   splitPts(1, i),
                               splitPts(0, i+1), splitPts(1, i+1);
                }
                else
                {
                    ctrlPts << splitPts(0, i), splitPts(1, i),
                               splitPts(0, 0), splitPts(1, 0);
                }
                const gsBSpline<> tempLine( basis, ctrlPts );
                mpTempSplitted->addPatch(gsBSpline<>(tempLine));
            }

            // Splitting ellipsis
            const gsGeometry<>& tempEllipse = temp.patch(1);
            gsMatrix<> splitParamsEllipse = splittingParameters(tempCase);
            gsMultiPatch<> mpSplittedEllipse = splittedCurveNew(tempEllipse, splitParamsEllipse);
            for (int i = 0; i < 10; i++)
            {
                mpTempSplitted->addPatch(mpSplittedEllipse[i]);
            }
            break;
        }
        default:
        {
            gsInfo << "Invalid template case" << std::endl;
        }
    };

    return mpTempSplitted;
}

/**************************************************/
// Returns a vector of multipatches which are boundaries of the template patches
std::vector< gsMultiPatch<> > templateBoundaries(const gsMultiPatch<>& mpTemp,
                                                 const gsMultiPatch<>& mpSegm,
                                                 const int tempCase)
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

    switch(tempCase)
    {
        case 10: // "UniKl22patches"
        {
            boundary0.addPatch(mpTemp[0]);
            boundary0.addPatch(mpTemp[1]);
            boundary0.addPatch(mpSegm[1]);
            boundary0.addPatch(mpSegm[0]);
            boundaries.push_back(boundary0);

            boundary1.addPatch(mpTemp[2]);
            boundary1.addPatch(mpSegm[2]);
            boundary1.addPatch(mpTemp[14]);
            boundary1.addPatch(mpSegm[1]);
            boundaries.push_back(boundary1);

            boundary2.addPatch(mpTemp[3]);
            boundary2.addPatch(mpSegm[3]);
            boundary2.addPatch(mpTemp[15]);
            boundary2.addPatch(mpSegm[2]);
            boundaries.push_back(boundary2);

            boundary3.addPatch(mpTemp[4]);
            boundary3.addPatch(mpSegm[4]);
            boundary3.addPatch(mpTemp[16]);
            boundary3.addPatch(mpSegm[3]);
            boundaries.push_back(boundary3);

            boundary4.addPatch(mpTemp[5]);
            boundary4.addPatch(mpSegm[5]);
            boundary4.addPatch(mpTemp[17]);
            boundary4.addPatch(mpSegm[4]);
            boundaries.push_back(boundary4);

            boundary5.addPatch(mpTemp[6]);
            boundary5.addPatch(mpSegm[6]);
            boundary5.addPatch(mpTemp[18]);
            boundary5.addPatch(mpSegm[5]);
            boundaries.push_back(boundary5);

            boundary6.addPatch(mpTemp[7]);
            boundary6.addPatch(mpSegm[7]);
            boundary6.addPatch(mpTemp[19]);
            boundary6.addPatch(mpSegm[6]);
            boundaries.push_back(boundary6);

            boundary7.addPatch(mpTemp[8]);
            boundary7.addPatch(mpSegm[8]);
            boundary7.addPatch(mpTemp[20]);
            boundary7.addPatch(mpSegm[7]);
            boundaries.push_back(boundary7);

            boundary8.addPatch(mpTemp[9]);
            boundary8.addPatch(mpSegm[9]);
            boundary8.addPatch(mpTemp[21]);
            boundary8.addPatch(mpSegm[8]);
            boundaries.push_back(boundary8);

            boundary9.addPatch(mpTemp[10]);
            boundary9.addPatch(mpSegm[10]);
            boundary9.addPatch(mpTemp[22]);
            boundary9.addPatch(mpSegm[9]);
            boundaries.push_back(boundary9);

            boundary10.addPatch(mpTemp[11]);
            boundary10.addPatch(mpSegm[11]);
            boundary10.addPatch(mpTemp[23]);
            boundary10.addPatch(mpSegm[10]);
            boundaries.push_back(boundary10);

            boundary11.addPatch(mpTemp[12]);
            boundary11.addPatch(mpTemp[13]);
            boundary11.addPatch(mpSegm[0]);
            boundary11.addPatch(mpSegm[11]);
            boundaries.push_back(boundary11);

            boundary12.addPatch(mpTemp[14]);
            boundary12.addPatch(mpSegm[12]);
            boundary12.addPatch(mpSegm[13]);
            boundary12.addPatch(mpTemp[23]);
            boundaries.push_back(boundary12);

            boundary13.addPatch(mpTemp[15]);
            boundary13.addPatch(mpSegm[15]);
            boundary13.addPatch(mpSegm[14]);
            boundary13.addPatch(mpSegm[12]);
            boundaries.push_back(boundary13);

            boundary14.addPatch(mpTemp[22]);
            boundary14.addPatch(mpSegm[16]);
            boundary14.addPatch(mpSegm[14]);
            boundary14.addPatch(mpSegm[13]);
            boundaries.push_back(boundary14);

            boundary15.addPatch(mpTemp[16]);
            boundary15.addPatch(mpSegm[18]);
            boundary15.addPatch(mpSegm[17]);
            boundary15.addPatch(mpSegm[15]);
            boundaries.push_back(boundary15);

            boundary16.addPatch(mpTemp[21]);
            boundary16.addPatch(mpSegm[19]);
            boundary16.addPatch(mpSegm[17]);
            boundary16.addPatch(mpSegm[16]);
            boundaries.push_back(boundary16);

            boundary17.addPatch(mpTemp[17]);
            boundary17.addPatch(mpSegm[22]);
            boundary17.addPatch(mpSegm[20]);
            boundary17.addPatch(mpSegm[18]);
            boundaries.push_back(boundary17);

            boundary18.addPatch(mpTemp[20]);
            boundary18.addPatch(mpSegm[23]);
            boundary18.addPatch(mpSegm[21]);
            boundary18.addPatch(mpSegm[19]);
            boundaries.push_back(boundary18);

            boundary19.addPatch(mpTemp[18]);
            boundary19.addPatch(mpSegm[26]);
            boundary19.addPatch(mpSegm[24]);
            boundary19.addPatch(mpSegm[22]);
            boundaries.push_back(boundary19);

            boundary20.addPatch(mpTemp[19]);
            boundary20.addPatch(mpSegm[23]);
            boundary20.addPatch(mpSegm[25]);
            boundary20.addPatch(mpSegm[26]);
            boundaries.push_back(boundary20);

            boundary21.addPatch(mpSegm[24]);
            boundary21.addPatch(mpSegm[25]);
            boundary21.addPatch(mpSegm[21]);
            boundary21.addPatch(mpSegm[20]);
            boundaries.push_back(boundary21);

            break;
        }
        default:
        {
            gsInfo << "Invalid template case" << std::endl;
            break;
        }
    };

    return boundaries;
}

}
