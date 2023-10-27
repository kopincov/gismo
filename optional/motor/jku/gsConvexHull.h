/** @file gsConvexHull.h

    @brief Functions used for the convex hull construction.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/

#include <gismo.h>
#include "gsMotorUtils.h"

using namespace gismo;

gsGeometry<>* convexHull(const gsGeometry<>& geom,
                         const int numPoints,
                         const int numCtrlPts,
                         const int deg,
                         const int numKnots,
                         const bool debug);
void markPoints(const gsMatrix<>& sampledPts,
                const gsMatrix<>& derivs,
                gsMatrix<index_t> &markedPts);
void identifyReplacementRegions(const gsMatrix<index_t> &markedPts,
                                gsMatrix<index_t> &replRegions);
index_t getReplacementRegionsNumber(const gsMatrix<index_t> &replRegions);
int getReplacementRegionLength(const index_t ind1,
                               const index_t ind2,
                               const index_t numPoints);
gsGeometry<>* fitCurve(gsMatrix<> const &vals,
                      gsMatrix<> const &pts,
                      const int deg,
                      int numKnots);


gsGeometry<>* convexHull(const gsGeometry<> &geom,
                         const int numPoints,
                         const int numCtrlPts,
                         const int deg,
                         const int numKnots,
                         const bool debug)
{
    // Step 1: Sampling points
    const gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints);
    const gsMatrix<> sampledPts = geom.eval(parameters);
    const gsMatrix<> derivs = geom.deriv(parameters);

    // Step 2: Marking points
    gsMatrix<index_t> markedPts(1, numPoints);
    markedPts.setZero();
    markPoints(sampledPts, derivs, markedPts);

    // Step 3: Replacing sample points
    // Identifying regions which need replacement
    gsMatrix<index_t> replRegions(2, numPoints);
    replRegions.setConstant(-1);
    identifyReplacementRegions(markedPts, replRegions);
    index_t replRegionsNum = getReplacementRegionsNumber(replRegions);

    // Replacing points
    gsMatrix<> convexHullPts(2, numPoints);
    convexHullPts = sampledPts;
    index_t replRegionBegin = 0;
    index_t replRegionEnd   = 0;
    gsMatrix<> replRegionBeginPoint;
    gsMatrix<> replRegionEndPoint;
    gsMatrix<> parametersSegment;

    int replRegionLength = 0;
    for (index_t i = 0; i != replRegionsNum; i++)
    {
        replRegionBegin = replRegions(0, i);
        replRegionEnd   = replRegions(1, i);
        if (replRegionBegin == 0)
        {
            replRegionBeginPoint = sampledPts.col(numPoints-1); // This point is not marked for replacement
        }
        else
        {
            replRegionBeginPoint = sampledPts.col(replRegionBegin-1); // This point is not marked for replacement
        }
        if (replRegionEnd == numPoints-1)
        {
            replRegionEndPoint = sampledPts.col(0); // This point is not marked for replacement
        }
        else
        {
            replRegionEndPoint = sampledPts.col(replRegionEnd+1); // This point is not marked for replacement
        }
        replRegionLength = getReplacementRegionLength(replRegionBegin, replRegionEnd, numPoints);

        if (replRegionBegin < replRegionEnd)
        {
            parametersSegment = uniformParameters<2>(0.0, 1.0, replRegionLength+2);
        }
        else
        {
            parametersSegment = uniformParameters<2>(0.0, 1.0, replRegionLength+1);
        }

        index_t iii = 1;
        for (index_t ii = replRegionBegin; ii != replRegionEnd+1; ii++)
        {
            convexHullPts(0, ii) = (replRegionEndPoint(0, 0)-replRegionBeginPoint(0, 0))*parametersSegment(iii) + replRegionBeginPoint(0, 0);
            convexHullPts(1, ii) = (replRegionEndPoint(1, 0)-replRegionBeginPoint(1, 0))*parametersSegment(iii) + replRegionBeginPoint(1, 0);
            if ( (ii ==  numPoints-1) && (ii != replRegionEnd))
            {
                // We reach the last point of the curve, but we need to continue replacing (ii != replRegionEnd)
                // Now the first point must be treated.
                ii = -1;
                // We do not increase iii, since the last and first point corespond to the same parameter
            }
            else
            {
                // We use next parameter
                iii++;
            }
        }
    }

    if (debug)
    {
        gsWriteParaviewPoints(convexHullPts, "convexHullPts");
    }

    // Step 4: Fitting
    /*
    gsKnotVector<> kv ( 0, 1,  numPoints-deg-1, deg+1 );
    gsBSplineBasis<> basis ( kv );
    gsGeometry<>* convexHull = *(fitbasis.interpolateData(convexHullPts, parameters);
    convexHull->coefs().row(0) = convexHull->coefs().row(convexHull->coefs().rows()-1);
    */
    /*
    gsKnotVector<> kv ( 0, 1,  numCtrlPts-deg-1, deg+1 );
    gsBSplineBasis<> basis( kv );
    gsFitting<> fitting(parameters, convexHullPts, basis);
    fitting.compute(1e-9);
    gsBSpline<>* convexHull = static_cast<gsBSpline<>*>(fitting.result()->clone());
    convexHull->coefs().row(0) = convexHull->coefs().row(convexHull->coefs().rows()-1);
    */

    /*
    gsKnotVector<> kv = getPeriodicKnotVector(deg, numKnots);
    gsCurveFitting<> fitting(parameters.transpose(), convexHullPts.transpose(), kv, true);
    fitting.compute_periodic();
    // calculate error
    // if error > threshold
    // increase number of knots
    // gsCurveFitting<> fitting(parameters.transpose(), convexHullPts.transpose(), kv, true);
    // fitting.compute_periodic();
    gsGeometry<>* convexHull = fitting.curve().clone();
    */

    gsGeometry<>* convexHull = fitCurve(convexHullPts, parameters, deg, numKnots);

    return convexHull;
}

void markPoints(const gsMatrix<>& sampledPts,
                const gsMatrix<>& derivs,
                gsMatrix<index_t> &markedPts)
{
    // Tangent line at point (c_x, c_y):
    //  p_y*(x-c_x) - p_x*(y-c_y) = 0

    const int numPoints = sampledPts.cols();
    int sign_prev = 0;
    int sign      = 0;

    // Going through all points and marking them
    for (index_t i = 0; i != numPoints-1; i++)
    {
        // First sample point "sign" wrt i-th tangent line
        sign_prev = sgn( derivs(1, i)*(sampledPts(0, 0)-sampledPts(0, i)) - derivs(0, i)*(sampledPts(1, 0)-sampledPts(1, i)) );
        // Checking rest of the sample point (until sign change is detected)
        for (index_t ii = 1; ii != numPoints; ii++)
        {
            if (i != ii)
            {
                sign = sgn( derivs(1, i)*(sampledPts(0, ii)-sampledPts(0, i)) - derivs(0, i)*(sampledPts(1, ii)-sampledPts(1, i)) );
                if ( (sign_prev != sign ) && (sign_prev != 0) && (sign != 0) )
                {
                    markedPts(0, i) = 1;
                    break; // Q: How to break out of the for loop (by i)
                    break;
                }
                else
                {
                    markedPts(0, i) = 0;
                }
            }
            sign_prev = sign;
        }
    }
}

void identifyReplacementRegions(const gsMatrix<index_t> &markedPts,
                                gsMatrix<index_t>& replRegions)
{

    const index_t numPoints = markedPts.cols();

    bool replacement = false;
    index_t replRegionNum = -1;
    index_t replRegionBegin = 0;
    index_t replRegionEnd   = 0;

    for (index_t i = 0; i != numPoints; i++)
    {
        if ( (markedPts(0, i) == 1) && (replacement == false) )
        {
            replRegionBegin = i;
            replacement = true;
            replRegionNum++;
            replRegions(0, replRegionNum) = replRegionBegin;
        }
        if ( (markedPts (0, i) == 0) && (replacement == true) )
        {
            replRegionEnd   = i-1;
            replacement = false;
            replRegions(1, replRegionNum) = replRegionEnd;
        }
        else if ( (i == numPoints-1) && (replacement == true) )
        {
            replRegionEnd   = i;
            replacement = false;
            replRegions(1, replRegionNum) = replRegionEnd;
        }
    }
    // If the first and the last sampled points should be modified,
    // then do not add new replacement region but modify beginning
    // point of the first region
    if ( (replRegionNum > 0) && (replRegions(0, 0) == 0) && (replRegions(1, replRegionNum) == numPoints-1) )
    {
         replRegions(0, 0) = replRegions(0, replRegionNum);
         replRegions(0, replRegionNum) = -1;
         replRegions(1, replRegionNum) = -1;
    }

}

index_t getReplacementRegionsNumber(const gsMatrix<index_t>& replRegions)
{
    index_t replRegionsNum = 0;

    for (index_t i = 0; i != replRegions.cols(); i++)
    {
        if ( (replRegions(0, i) != -1) && (replRegions(1, i) != -1) )
        {
            replRegionsNum++;
        }
        else
        {
            break;
        }
    }

    return replRegionsNum;
}

int getReplacementRegionLength(const index_t ind1,
                               const index_t ind2,
                               const index_t numPoints)
{
    if (ind2 < ind1)
    {
        return ind2 - ind1 + numPoints + 1;
    }
    else if (ind2 > ind1)
    {
        return ind2 - ind1 + 1;
    }
    else
    {
        return 0;
    }
}

gsGeometry<>* fitCurve(gsMatrix<> const &vals,
                      gsMatrix<> const &pts,
                      const int deg,
                      int numKnots)
{

    // Fitting with interpolateData

    const int numPoints = pts.cols();
    gsKnotVector<> kv ( 0, 1,  numPoints-deg-1, deg+1 );
    gsBSplineBasis<> basis ( kv );
    gsGeometry<>* geom = basis.interpolateData(vals, pts).release();
    geom->coefs().row(0) = geom->coefs().row(geom->coefs().rows()-1);

    // Fitting with gsFitting
    /*
    int numPoints = pts.cols();
    const real_t lambda = 1.0e-6;
    gsKnotVector<> kv ( 0, 1,  numPoints-deg-1, deg+1 );
    gsBSplineBasis<> basis( kv );
    gsFitting<> fitting(pts, vals, basis);
    fitting.compute(lambda);
    gsGeometry<>* geom = fitting.result();
    geom->coefs().row(0) = geom->coefs().row(geom->coefs().rows()-1);
    */

    // Fitting with gsCurveFitting
    /*
    real_t approxError = 1.0e+3; // Initial value (large)
    real_t threshold = 1.0e-4;
    gsKnotVector<> kv = getPeriodicKnotVector(deg, numKnots);
    gsCurveFitting<> fitting(pts.transpose(), vals.transpose(), kv, true);
    fitting.compute_periodic();
    fitting.computeApproxError(approxError);
    gsGeometry<>* geom = static_cast<gsGeometry<>*>(fitting.curve().clone());
    while (threshold < approxError)
    {
        numKnots *= 2;
        kv = getPeriodicKnotVector(deg, numKnots);
        gsCurveFitting<> fitting2(pts.transpose(), vals.transpose(), kv, true);
        fitting2.compute_periodic();
        fitting2.computeApproxError(approxError);
        geom = static_cast<gsGeometry<>*>(fitting2.curve().clone());
    }
    */

    return geom;
}
