/** @file gsMinimizedLengthOffset.cpp

    @brief Construction of minimzed length ("rubber band") offset for non-convex smooth closed curve.
    Some modifications of the offset are possible.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius, J. Speh
*/

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsConvexHull.h"


using namespace gismo;

gsGeometry<>* offset(const gsGeometry<> &geom,
                     const int numPoints,
                     const int numKnots,
                     const real_t distance);
void removeSelfIntersections(const gsGeometry<> &geom,
                             const int numPoints);
void modifyOffset(const gsGeometry<>& geom,
                  gsGeometry<>& offset,
                  const int numPoints,
                  const int numKnots,
                  const real_t percentage);
gsMatrix<> pseudoNormals(const gsGeometry<>& geom1,
                         const gsGeometry<>& geom2,
                         const int numPoints,
                         const bool debug);
gsMatrix<> smoothedPseudoNormals(const gsGeometry<>& geom1,
                                 const gsGeometry<>& geom2,
                                 const int numPoints,
                                 int numKnots,
                                 const bool debug);
void changePseudoNormals(const gsGeometry<> &geom1,
                         gsGeometry<> &geom2,
                         const gsMatrix<>& pseudoNrmls,
                         int numKnots,
                         const bool debug);
//std::vector<real_t> getPeriodicKnots(const int deg,
//                                     const int numKnots);
//gsKnotVector<> getPeriodicKnotVector(const int deg,
//                                     const int numKnots);

int main(int argc, char *argv[])
{
    // Options with defaul t values
    std::string inputFile(MOTOR_DATA_DIR "jku/target_input.xml");
    index_t numPoints  = 1000;
    index_t numCtrlPts = 100;
    index_t numKnots   = 100;
    index_t degree     = 2;
    real_t distance = 0.1;
    real_t percentage = 1.0;
    real_t smoothing  = 1.0e-8;
    bool debug = false;
    
    std::string outputFile("result_offsetting");

    gsCmdLine cmd("Construction of minimized length (\"rubber band\") offset for non-convex smooth closed curve");
    cmd.addString("f", "input", "Name of input file", inputFile);
    cmd.addInt("n", "npts", "Number of sampling points (for sampling)", numPoints);
    cmd.addInt("N", "ncpts", "Number of control points", numCtrlPts);
    cmd.addInt("K", "nknots", "Number of knots used fo curve fitting", numKnots);
    cmd.addInt("D", "degree", "Degree", degree);
    cmd.addReal("d", "dist", "The value of offsetting distance", distance);
    cmd.addReal("p", "perc", "...", percentage);
    cmd.addReal("s", "smoothing", "...", smoothing);
    cmd.addString("o", "output", "Name of output file", outputFile);
    cmd.addSwitch("debug", "If debug, it outputs more data.", debug);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input file:  " << inputFile << "\n"
           << "numPoints:   " << numPoints << "\n"
           << "numCtrlPts:  " << numCtrlPts << "\n"
           << "numKnots:    " << numKnots << "\n"
           << "degree:      " << degree << "\n"
           << "distance:    " << distance << "\n"
           << "percentage:  " << percentage << "\n"
           << "Output file: " << outputFile << "\n"
           << "Debug:       " << debug << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsGeometry<>* p_geometry = readGeometry(inputFile);
    gsGeometry<>& geometry = *p_geometry;

    // (1) CONSTRUCTING THE CONVEX HULL
    gsGeometry<>* CH = convexHull(geometry, numPoints, numCtrlPts, degree, numKnots, debug);

    // (2) CONSTRUCTING (A GENERALIZED) OFFSET OF THE CONVEX HULL
    gsGeometry<>* offsetCH = offset(*CH, numPoints, numKnots, distance);
    if ( percentage < 1.0 )
    //
    {
        modifyOffset(geometry, *offsetCH, numPoints, numKnots, percentage);
    }

    if (debug)
    {
        writeGeometry(*offsetCH, outputFile + "_offsetConvexHull");
    }

    /*
    // (3) SMOOTHING NORMAL ANGLES
    gsMatrix<> smoothedPseudoNrmls = smoothedPseudoNormals(geometry, *offsetCH, numPoints,
                                                           numKnots, debug);

    // FINALIZING
    changePseudoNormals(geometry, *offsetCH, smoothedPseudoNrmls, numKnots, debug);
    */
    writeGeometry(*offsetCH, outputFile + "_offsetConvexHull_FINAL");


    if (debug)
    {
        writeGeometry(*CH, outputFile + "_convexHull");

    }
    
    return 0;

}

gsGeometry<>* offset(const gsGeometry<>& geom,
                     const int numPoints,
                     const int numKnots,
                     const real_t distance)
{
    gsBSpline<> inputGeom = dynamic_cast< const gsBSpline<>& >(geom);
    gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints);
    gsMatrix<> geomPts = inputGeom.eval(parameters);
    // Normal vectors
    gsMatrix<> derivs = inputGeom.deriv(parameters);
    gsMatrix<> nrmlVects = rotateVectors(derivs, -EIGEN_PI/2.0);
    normalize(nrmlVects);
    // Offset
    gsMatrix<> offsetPts(2, numPoints);
    offsetPts.setZero();
    for (index_t i = 0; i != numPoints; i++)
    {
        offsetPts.col(i) = geomPts.col(i) + distance*nrmlVects.col(i);
    }

    // Final fitting
    /*
    int deg = inputGeom.degree();
    gsKnotVector<> kv ( 0, 1,  numPoints-deg-1, deg+1 );
    gsBSplineBasis<> basis ( kv );
    gsGeometry<>* offset = basis.interpolateData(offsetPts, parameters);
    offset->coefs().row(0) = offset->coefs().row(offset->coefs().rows()-1);
    */
    /*
    int deg = inputGeom.degree();
    gsKnotVector<> kv = getPeriodicKnotVector(deg, numKnots);
    gsCurveFitting<> fitting(parameters.transpose(), offsetPts.transpose(), kv, true);
    fitting.compute_periodic();
    gsGeometry<>* offset = fitting.curve().clone();
    */

    int deg = inputGeom.degree();
    gsGeometry<>* offset = fitCurve(offsetPts, parameters, deg, numKnots);

    return offset;
}

void removeSelfIntersections(const gsGeometry<> &geom,
                             const int numPoints)
{
    gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints);
    //gsMatrix<> sampledPts = inputGeometry.eval(parameters);
    //gsWriteParaviewPoints(sampledPts, "sampledPts");
    gsMatrix<> derivs = geom.deriv(parameters);
    gsMatrix<> derivs2 = geom.deriv2(parameters);

    gsMatrix<> curvs = getCurvatures(derivs, derivs2);
    gsInfo << "curvs:\n" << curvs << std::endl;



}


void modifyOffset(const gsGeometry<> &geom,
                  gsGeometry<>& offset,
                  const int numPoints,
                  const int numKnots,
                  const real_t percentage)
{
    gsBSpline<> inputGeometry = dynamic_cast< const gsBSpline<>& >(geom);
    gsBSpline<> inputOffset = dynamic_cast< const gsBSpline<>& >(offset);

    gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints);
    gsMatrix<> geomPts   = inputGeometry.eval(parameters);
    gsMatrix<> offsetPts = inputOffset.eval(parameters);

    for (index_t i = 0; i < numPoints; i++)
    {
        offsetPts.col(i) =  geomPts.col(i) + percentage*(offsetPts.col(i) - geomPts.col(i));
    }
    /*
    int deg = inputGeometry.degree();
    gsKnotVector<> kv ( 0, 1,  numPoints-deg-1, deg+1 );
    gsBSplineBasis<> basis ( kv );
    offset = *(basis.interpolateData(offsetPts, parameters));
    offset.coefs().row(0) = offset.coefs().row(offset.coefs().rows()-1);
    */
    /*
    int deg = inputGeometry.degree();
    gsKnotVector<> kv = getPeriodicKnotVector(deg, numKnots);
    gsCurveFitting<> fitting(parameters.transpose(), offsetPts.transpose(), kv, true);
    fitting.compute_periodic();
    offset = fitting.curve();
    */

    int deg = inputGeometry.degree();
    offset = *(fitCurve(offsetPts, parameters, deg, numKnots));

}

gsMatrix<> pseudoNormals(const gsGeometry<>& geom1,
                         const gsGeometry<>& geom2,
                         const int numPoints,
                         const bool debug)
{
    gsBSpline<> inputGeom1 = dynamic_cast< const gsBSpline<>& >(geom1);
    gsBSpline<> inputGeom2 = dynamic_cast< const gsBSpline<>& >(geom2);

    gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints);
    gsMatrix<> geom1Pts   = inputGeom1.eval(parameters);
    gsMatrix<> geom2Pts   = inputGeom2.eval(parameters);

    // Smoothing normal angles of convex hull offset
    gsMatrix<> pseudoNrmls(2, numPoints);
    for (index_t i = 0; i != numPoints; i++)
    {
        pseudoNrmls.col(i) = geom2Pts.col(i) - geom1Pts.col(i);
    }

    if (debug)
    {
        gsWriteParaviewPoints(pseudoNrmls, "pseudoNrmls");
        plotNormals(geom1.eval(parameters), pseudoNrmls, "pseudoNrmls_lines", 1.0);
    }

    return pseudoNrmls;
}

gsMatrix<> smoothedPseudoNormals(const gsGeometry<>& geom1,
                                 const gsGeometry<>& geom2,
                                 const int numPoints,
                                 int numKnots,
                                 const bool debug)
{
    gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints);

    gsMatrix<> pseudoNrmls = pseudoNormals(geom1, geom2, numPoints, debug);

    int deg = 2;
    if (numPoints < numKnots)
    {
        gsInfo << "smoothedPseudoNormals: numKnots is decreased" << std::endl;
        numKnots = numPoints;
    }

    //int deg = inputGeometry.degree();
    /*
    gsKnotVector<> kv ( 0, 1,  numPoints-deg-1, deg+1 );
    gsBSplineBasis<> basis ( kv );
    gsGeometry<>* curve = (basis.interpolateData(pseudoNrmls, parameters));
    curve->coefs().row(0) = curve->coefs().row(curve->coefs().rows()-1);
    gsMatrix<> smoothedPseudoNrmls = curve->eval(parameters);
    */
    /*
    gsKnotVector<> kv ( 0, 1,  numPoints-deg-1, deg+1 );
    gsBSplineBasis<> basis( kv );
    gsFitting<> fitting(parameters, pseudoNrmls, basis);
    fitting.compute(1e-6);
    gsBSpline<>* curve = static_cast<gsBSpline<>*>(fitting.result()->clone());
    curve->coefs().row(0) = curve->coefs().row(curve->coefs().rows()-1);
    gsMatrix<> smoothedPseudoNrmls = curve->eval(parameters);
    */
    /*
    gsKnotVector<> kv = getPeriodicKnotVector(deg, numKnots);
    gsCurveFitting<> fitting(parameters.transpose(), pseudoNrmls.transpose(), kv, true);
    fitting.compute_periodic();
    gsBSpline<> curve = fitting.curve();
    gsMatrix<> smoothedPseudoNrmls = curve.eval(parameters);
    */

    gsGeometry<>* curve = fitCurve(pseudoNrmls, parameters, deg, numKnots);
    gsMatrix<> smoothedPseudoNrmls = curve->eval(parameters);

    if (debug)
    {
        gsWriteParaviewPoints(smoothedPseudoNrmls, "smoothedPseudoNrmls");
        plotNormals(geom1.eval(parameters), smoothedPseudoNrmls, "smoothedPseudoNrmls_lines", 1.0);
    }

    return smoothedPseudoNrmls;
}

void changePseudoNormals(const gsGeometry<>& geom1,
                         gsGeometry<>& geom2,
                         const gsMatrix<>& pseudoNrmls,
                         int numKnots,
                         const bool debug)
{
    gsBSpline<> inputGeom1 = dynamic_cast< const gsBSpline<>& >(geom1);
    gsBSpline<> inputGeom2 = dynamic_cast< const gsBSpline<>& >(geom2);

    int numPoints = pseudoNrmls.cols();
    if (numPoints < numKnots)
    {
        gsInfo << "changePseudoNormals: numKnots is decreased" << std::endl;
        numKnots = numPoints;
    }
    gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints);
    gsMatrix<> geom1Pts   = inputGeom1.eval(parameters);
    gsMatrix<> geom2Pts   = inputGeom2.eval(parameters);

    for (index_t i = 0; i < numPoints; i++)
    {
        geom2Pts.col(i) =  geom1Pts.col(i) + pseudoNrmls.col(i);
    }

    if (debug)
    {
        gsWriteParaviewPoints(geom2Pts, "smoothedPseudoNrmls_pts");
    }

    int deg = inputGeom2.degree();
    geom2 = *(fitCurve(geom2Pts, parameters, deg, numKnots));

}
