/** @file gsMotorUtils.h

    @brief This file contains various auxilary routines used by Geometry workpackage
    in the MOTOR project (http://motor-project.eu)

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Sajavicius, J. Speh
*/

#pragma once

#include <gsModeling/gsCoonsPatch.h>
#include "gsMotorIOUtils.h"
#include "gsUtils/gsExportMatrix.h"

namespace gismo
{

gsMatrix<> rotateVectors(const gsMatrix<>& vects,
                         const real_t rad);
gsMatrix<> rotateVectors(const gsMatrix<>& vects,
                         const real_t rad,
                         const gsVector<>& center);
gsMatrix<> scaleVectors(const gsMatrix<>& vects,
                        const real_t sFact1,
                        const real_t sFact2);
gsMatrix<> translateVectors(const gsMatrix<>& vects,
                            const real_t tConst1,
                            const real_t tConst2);
void normalize(gsMatrix<>& vects);
real_t getAngle(gsMatrix<> vector);
gsMatrix<> getAngles(const gsMatrix<>& vectors);
gsMatrix<> getCurvatures(const gsMatrix<>& der,
                         const gsMatrix<>& der2);
gsVector<> curvatures_ctrlpts(const gsMatrix<>& ctrlPts,
                              const gsMatrix<>& nrmlVects);
void make_c0(gsMatrix<>& ctrlPts);
void make_c0(gsGeometry<>& geom);
void make_g1(gsMatrix<>& ctrlPts);
void make_g1(gsGeometry<>& geom);
//real_t min(real_t a, real_t b);
//real_t max(real_t a, real_t b);
real_t determinant(const gsVector<>& vec1,
                   const gsVector<>& vec2);
index_t findClosestIndex(const gsMatrix<>& point,
                         const gsMatrix<>& pointCloud);
gsMatrix<> binarySearch(const gsMatrix<>& point,
                        const gsGeometry<>& offset,
                        const gsMatrix<>& lowStart,
                        const gsMatrix<>& uppStart,
                        const real_t tolerance,
                        const index_t maxSteps);
gsMatrix<> findParameter(const gsGeometry<>& geom,
                         const gsMatrix<>& point,
                         const real_t tolerance);
real_t percentPointBelowTreshold(const std::vector<real_t>& errors,
                                 const real_t treshold);
gsMatrix<> reverseMatrixRows(const gsMatrix<> mat);
gsMatrix<> reverseMatrixCols(const gsMatrix<> mat);
gsMatrix<> getNormal(const real_t angle);
gsMatrix<> getNormals(const gsMatrix<>& angles);
void plotNormals(const gsMatrix<>& points,
                 const gsMatrix<>& normals,
                 const std::string& file,
                 const real_t length);
gsGeometry<>::uPtr fitAngles(const gsMatrix<>& parameters,
                        const gsMatrix<>& angles,
                        const int numCtrlPts,
                        const int degree);
//gsGeometry<>* makeCurve(const gsMatrix<>& parameters,
//                        const gsMatrix<>& points,
//                        const int numCtrlPts,
//                        const int degree);
gsMatrix<index_t> generatePermutationMatrix(const index_t num);
gsMatrix<> wachspressCoordinates(const gsMatrix<> anchors,
                                 const gsMatrix<> points);
real_t triangleArea(const gsMatrix<> points);
gsMatrix<> mapPointsWachspress(const gsMatrix<> templateAnchors,
                               const gsMatrix<> targetAnchors,
                               const gsMatrix<> points);
std::size_t rectangleIndex(const std::string& inputFile);
std::size_t rectangleIndex(const gsMultiPatch<>& multipatch);
gsMatrix<> detectVertices(const std::string& inputFile);
gsMatrix<> detectVertices(const gsGeometry<>& geom);
gsMultiPatch<>* removePatch(gsMultiPatch<>& mp,
                            const std::size_t id);
gsMatrix<> arcLengthParameters(gsMatrix<>& points); // This function is commented
std::vector<real_t> getPeriodicKnots(const int deg,
                                     const int numKnots);
gsKnotVector<> getPeriodicKnotVector(const int deg,
                                     const int numKnots);
gsMultiPatch<> splittedCurve (const gsGeometry<>& geom,
                              const gsMatrix<>& points);
gsMultiPatch<> splittedCurveParameters (const gsGeometry<>& geom,
                                        const gsMatrix<>& pars);
// Splits B-Spline curve at points corresponding
// to given parameteric values
gsMultiPatch<>* splitOpenBSpline(const gsBSpline<>& geom,
                                 const gsMatrix<>& pars);

// Takes two multipatches and joins them
gsMultiPatch<>* joinMultipatches(const std::string& mpFile1,
                                 const std::string& mpFile2);
gsMultiPatch<>* joinMultipatches(const gsMultiPatch<>& mp1,
                                 const gsMultiPatch<>& mp2);

gsMatrix<> sumMatrixRows (const gsMatrix<>& mat);
gsMatrix<> sumMatrixCols (const gsMatrix<>& mat);

// Functions from gsMotorIOUtils.h
gsMultiPatch<>* readMultipatch(const std::string& inputFile);
//void writeMultipatch(const gsMultiPatch<>& multipatch,
//                     const std::string& outputFile);

gsBSpline<> constructEllipse(const real_t x_center,
                              const real_t y_center,
                              const real_t r1,
                              const real_t r2,
                              const int numCtrlPts,
                              const int degree);
// The mapping m maps the multipatch mp
gsMultiPatch<> mapMultipatch(const gsGeometry<>& map,
                             const gsMultiPatch<>& mp);
// Map points given in ASCII file
/*
void mapPoints(const std::string& inputFile,
               const gsGeometry<>& map,
               const std::string& outputPrefix);
*/
// Map points given in gsMatrix format
void mapPoints(const gsMatrix<> pts,
               const gsGeometry<>& map,
               const std::string& outputPrefix);

gsMatrix<> stdVectorToMatrix(const std::vector< gsMatrix<> >& vector);
//std::vector<> MatrixToStdVectorRowwise(const gsMatrix<> & mat);
//std::vector<> MatrixToStdVectorColwise(const gsMatrix<> & mat);

real_t detJacobian(const gsGeometry<>& geom,
                   const gsMatrix<> point);
void checkJacobianDeterminant(const gsGeometry<>& geom,
                              const int points,
                              const bool savePoints,
                              const std::string& output,
                              const int number);
gsMultiPatch<>* extrudeMultipatch(const gsMultiPatch<>& mp, const real_t length);

// Constructs an std::vector of Coon's patches defined by boundaries stored in std::vector
std::vector< gsTensorBSpline<3, real_t> > getCoonsPatchesVolumetric(const std::vector< gsMultiPatch<> > boundaries);
std::vector< gsTensorBSpline<2, real_t> > getCoonsPatchesPlanar( const std::vector< gsMultiPatch<> > boundaries);

// TO DO: change this d-1 -> d
template <unsigned d>
gsMatrix<> uniformParameters(const real_t lower,
                             const real_t upper,
                             const int numPoints = 100)
{
    gsVector<> start(d-1);
    start.setConstant(lower);
    gsVector<> end(d-1);
    end.setConstant(upper);
    return uniformPointGrid(start, end, numPoints);
}

template <unsigned d>
gsMatrix<> samplePoints(const gsGeometry<>& geom,
                        const int numPoints = 100)
{
    const gsMatrix<> parameters = uniformParameters<d>(0.0, 1.0, numPoints); // different
    gsMatrix<> points = geom.eval(parameters);
    return points;
}

template <unsigned d>
gsMatrix<> samplePoints(const gsMultiPatch<>& mp,
                        const int numPoints = 100)
{
    const gsMatrix<> parameters = uniformParameters<d>(0.0, 1.0, numPoints); // different
    gsMatrix<> points(d, parameters.cols() * mp.nPatches());
    for (std::size_t i = 0; i != mp.nPatches(); i++)
    {
        const gsGeometry<>& geom = mp[i];
        points.block(0, i * parameters.cols(), d, parameters.cols()) = samplePoints<d>(geom, numPoints);
    }

    return points;
}

template <unsigned d>
gsMatrix<> samplePoints(const std::string& inputFile,
                        const int numPoints = 100)
{
    gsMultiPatch<>* mp = readMultipatch(inputFile);
    gsMatrix<> parameters = uniformParameters<d>(0.0, 1.0, numPoints); // different
    gsMatrix<> points(d, parameters.cols() * mp->nPatches());
    for (std::size_t i = 0; i != mp->nPatches(); i++)
    {
        const gsGeometry<>& geom = (*mp)[i];
        gsMatrix<> geomPts = geom.eval(parameters);
        points.block(0, i * parameters.cols(), d, parameters.cols()) = geomPts;
    }

    delete mp;

    return points;
}

template <unsigned d>
void mapMultipatch(const std::string& inputFile,
                   const gsGeometry<>& map,
                   const int numPoints = 100,
                   const std::string& outputPrefix = "result_")
{
    gsMultiPatch<>* inputMultipatch = readMultipatch(inputFile);
    gsMultiPatch<> resultMultipatch;
    int n = numPoints;
    if(d==3)
    {
        n = numPoints*numPoints;
    }
    gsMatrix<> parameters = uniformParameters<d>(0.0, 1.0, n);
    for (std::size_t i = 0; i !=  inputMultipatch->nPatches(); i++)
    {
        const gsGeometry<>& inputPatch = (*inputMultipatch)[i];
        gsMatrix<> inputPatchPoints = inputPatch.eval(parameters);

        gsMatrix<> inputPatchPointsMapped = map.eval(inputPatchPoints);
        gsWriteParaviewPoints(inputPatchPointsMapped, outputPrefix + "_Points" + util::to_string(i));
        if(d == 2)
        {
            gsKnotVector<> KV(0.0, 1.0, numPoints-2, 2);
            gsBSpline<>::uPtr inputPatchMapped(new gsBSpline<>(KV, inputPatchPointsMapped.transpose()));
            gsWriteParaview(*inputPatchMapped, outputPrefix + "_Curve" + util::to_string(i));
            resultMultipatch.addPatch(give(inputPatchMapped));
        }
        else if(d == 3)
        {
            gsKnotVector<> KV(0.0, 1.0, numPoints-2, 2);
            gsTensorBSpline<2>::uPtr inputPatchMapped(new gsTensorBSpline<2>(KV, KV, inputPatchPointsMapped.transpose()));
            gsWriteParaview(*inputPatchMapped, outputPrefix + "_Surface" + util::to_string(i));
            resultMultipatch.addPatch(give(inputPatchMapped));
        }
        else
        {
            gsWarn << "Multipatch mapping is implemented in 2D and 3D cases only" << std::endl;
        }
    }

    writeMultipatch(resultMultipatch, outputPrefix);

    delete inputMultipatch;
}

template <unsigned d>
gsMatrix<> boundingBox(const gsMatrix<>& points,
                       const real_t percent)
{
    // first column lower corner, second column upper corner
    gsMatrix<> bbox(d, 2);
    bbox.col(0) = points.col(0);
    bbox.col(1) = points.col(0);

    for (index_t col = 1; col != points.cols(); col++)
    {
        for (index_t row = 0; row != points.rows(); row++)
        {
            if ( points(row, col) < bbox(row, 0) )
            {
                bbox(row, 0) = points(row, col);
            }
            if ( bbox(row, 1) < points(row, col) )
            {
                bbox(row, 1) = points(row, col);
            }
        }
    }

    for (index_t row = 0; row != bbox.rows(); row++)
    {
        real_t dx = bbox(row, 1) - bbox(row, 0);
        bbox(row, 0) = bbox(row, 0) - percent * dx;
        bbox(row, 1) = bbox(row, 1) + percent * dx;
    }

    return bbox;
}

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

gsMatrix<> rotateVectors(const gsMatrix<>& vects,
                         const real_t rad)
{
    //GISMO_ENSURE(vects.rows() == 2, "Only 2D vectors can be rotated");
    const size_t d = vects.rows();
    gsMatrix<> transformationMat(d, d);
    transformationMat.setZero();
    // 2D rotation matrix
    transformationMat << math::cos(rad), -math::sin(rad),
                         math::sin(rad),  math::cos(rad);

    return transformationMat*vects;
}

gsMatrix<> rotateVectors(const gsMatrix<>& vects,
                         const real_t rad,
                         const gsVector<>& center)
{
    const size_t d = vects.rows();
    gsMatrix<> transformationMat(d, d);
    transformationMat.setZero();
    // 2D rotation matrix
    transformationMat << math::cos(rad), -math::sin(rad),
                         math::sin(rad),  math::cos(rad);

    gsMatrix<> pts(d, vects.cols());
    pts = vects;
    //gsInfo << "pts:\n" << pts << std::endl;
    pts.colwise() -= center;
    //gsInfo << "pts:\n" << pts << std::endl;
    pts = transformationMat*pts;
    //gsInfo << "pts:\n" << pts << std::endl;
    pts.colwise() += center;
    //gsInfo << "pts:\n" << pts << std::endl;


    /*
    gsMatrix<> transformationMat2(d, d);
    transformationMat2.setZero();
    transformationMat2 << 1-math::cos(rad),  math::sin(rad),
                          1-math::sin(rad), -math::cos(rad);

    gsMatrix<> mat(d, vects.cols());
    mat.setZero();
    mat = transformationMat*vects;

    gsVector<> vect(2);
    vect = transformationMat2*center;
    mat.colwise() += vect;

    */

    return pts;
}

gsMatrix<> scaleVectors(const gsMatrix<>& vects,
                        const real_t sFact1,
                        const real_t sFact2)
{
    //GISMO_ENSURE(vects.rows() == 2, "Only 2D vectors can be scaled");
    const size_t d = vects.rows();
    gsMatrix<> transformationMat(d, d);
    transformationMat.setZero();
    // 2D scaling matrix
    transformationMat << sFact1, 0.0,
                         0.0,    sFact2;

    return transformationMat*vects;
}

gsMatrix<> translateVectors(const gsMatrix<>& vects,
                            const real_t tConst1,
                            const real_t tConst2)
{
    //GISMO_ENSURE(vects.rows() == 2, "Only 2D vectors can be translated");
    const size_t d = vects.rows();
    gsMatrix<> transformationVect(d, 1);
    transformationVect.setZero();
    // 2D translating vector
    transformationVect(0, 0) = tConst1,
    transformationVect(1, 0) = tConst2;

    gsMatrix<> translatedVect(d, vects.cols());
    translatedVect.setZero();
    for (index_t i = 0; i < vects.cols(); i++)
    {
         translatedVect.col(i) += vects.col(i) + transformationVect;
    }

    return translatedVect;
}


void normalize(gsMatrix<>& vects)
{
    for (index_t i = 0; i != vects.cols(); i++)
    {
        vects.col(i) /= vects.col(i).norm();
    }
}

real_t getAngle(gsMatrix<> vector)
{
    const real_t x = vector(0, 0);
    const real_t y = vector(1, 0);

    real_t phi = math::acos(x);
    if (y < 0)
    {
        return 2*EIGEN_PI - phi;
    }
    else
    {
        return phi;
    }
}

gsMatrix<> getAngles(const gsMatrix<>& vectors)
{
    gsMatrix<> angles(1, vectors.cols());
    gsMatrix<> normalizedVectors = vectors;
    normalize(normalizedVectors);
    real_t delta =  0.0;
    real_t angle_old = 0.0;
    real_t angle_new = 0.0;
    const real_t eps = 1.0e-12;

    for (index_t col = 0; col != normalizedVectors.cols(); col++)
    {
        angle_new = getAngle(normalizedVectors.col(col));

        if ( (col > 0) &&  (math::abs(angle_old - angle_new - 2*EIGEN_PI) < eps) ) // Jump down from 2PI to 0
        {
            delta += angle_old;
        }
        else if ( (col > 0) &&  (math::abs(angle_new - angle_old - 2*EIGEN_PI) < eps) ) // Jump up from 0 to 2PI
        {
            delta -= angle_new;
        }
        angles(0, col) =  angle_new + delta;
        angle_old = angle_new;
    }

    return angles;
}

gsMatrix<> getCurvatures(const gsMatrix<>& der,
                         const gsMatrix<>& der2)
{
    gsMatrix<> curvs(1, der.cols());

    for (index_t col = 0; col != der.cols(); col++)
    {
        const gsVector<>& dC = der.col(col);
        const gsVector<>& ddC = der2.col(col);
        const real_t norm = dC.norm();
        const real_t norm3 = norm * norm * norm;
        curvs(col) = determinant(dC, ddC) / norm3;
    }

    return curvs;
}

gsVector<> curvatures_ctrlpts(const gsMatrix<>& ctrlPts,
                              const gsMatrix<>& nrmlVects)
{
    GISMO_ENSURE(ctrlPts.dim() == nrmlVects.dim(), "Matrices of control points and normal vectors must have the same dimension");
    size_t d = ctrlPts.cols();
    size_t num = ctrlPts.rows();
    gsVector<> kappa(num);
    kappa.setZero();
    gsVector<> tmp(d);
    tmp.setZero();

    tmp = ctrlPts.row(1) - 2.0*ctrlPts.row(0) + ctrlPts.row(num-2);
    kappa(0) = tmp.dot(nrmlVects.row(0));
    for(std::size_t i = 1; i != num-1; i++)
    {
        tmp = ctrlPts.row(i+1) - 2.0*ctrlPts.row(i) + ctrlPts.row(i-1);
        kappa(i) = tmp.dot(nrmlVects.row(i));
    }
    tmp = ctrlPts.row(1) - 2.0*ctrlPts.row(num-1) + ctrlPts.row(num-2);
    //tmp /= tmp.norm();
    kappa(num-1) = tmp.dot(nrmlVects.row(num-1));

    return kappa;
}

void make_c0(gsMatrix<>& ctrlPts)
{

    const index_t end = ctrlPts.rows() - 1;

    ctrlPts.row(1) = ctrlPts.row(end);
}

void make_c0(gsGeometry<>& geom)
{
    gsMatrix<> ctrlPts = geom.coefs();
    make_c0(ctrlPts);
    geom.setCoefs(ctrlPts);
}

void make_g1(gsMatrix<>& ctrlPts)
{
    // Assumptions:
    // 1. We have more than 4 control points
    // 2. First and last control points coincide

    const index_t end = ctrlPts.rows() - 1;

    const gsVector<> tangent1 = ctrlPts.row(1) - ctrlPts.row(0);
    const gsVector<> tangent2 = ctrlPts.row(end) - ctrlPts.row(end - 1);
    const gsVector<> average = (tangent1 + tangent2) * 0.5;

    ctrlPts.row(1) = ctrlPts.row(0) + average.transpose();
    ctrlPts.row(end - 1) = ctrlPts.row(end) - average.transpose();
}

void make_g1(gsGeometry<>& geom)
{
    gsMatrix<> ctrlPts = geom.coefs();
    make_g1(ctrlPts);
    geom.setCoefs(ctrlPts);
}

/*
real_t min(real_t a, real_t b)
{
    return a < b ? a : b;
}

real_t max(real_t a, real_t b)
{
    return a < b ? b : a;
}
*/

real_t determinant(const gsVector<>& vec1,
                   const gsVector<>& vec2)
{
    return vec1(0) * vec2(1) - vec1(1) * vec2(0);
}


index_t findClosestIndex(const gsMatrix<>& point,
                         const gsMatrix<>& pointCloud)
{
    real_t distance = -1;
    index_t index = -1;

    for (index_t i = 0; i != pointCloud.cols(); i++)
    {
        const real_t localDist = (point - pointCloud.col(i)).norm();

        if (i == 0 || localDist < distance)
        {
            distance = localDist;
            index = i;
        }
    }

    return index;
}

gsMatrix<> binarySearch(const gsMatrix<>& point,
                        const gsGeometry<>& offset,
                        const gsMatrix<>& lowStart,
                        const gsMatrix<>& uppStart,
                        const real_t tolerance,
                        const index_t maxSteps)
{
    index_t step = 1;

    gsMatrix<> low = lowStart;

    gsMatrix<> upp = uppStart;
    gsMatrix<> uppPoint = offset.eval(upp);
    real_t uppDist = (point - uppPoint).norm();

    gsMatrix<> mid = (low + upp) / 2;
    gsMatrix<> midPoint = offset.eval(mid);
    real_t midDist = (point - midPoint).norm();

    while (true)
    {
        if (midDist < tolerance || maxSteps < step)
        {
            return mid;
        }

        if (midDist < uppDist)
        {
            upp = mid;
            uppDist = midDist;
        }
        else
        {
            low = mid;
            // lowDist = midDist;
        }

        mid = (low + upp) / 2;
        midPoint = offset.eval(mid);
        midDist = (point - midPoint).norm();

        step++;
    }

    return mid;
}

gsMatrix<> findParameter(const gsGeometry<>& geom,
                         const gsMatrix<>& point,
                         const real_t tolerance)
{

    gsMatrix<> uniParams = uniformParameters<2> (0.0, 1.0, 1000);
    gsMatrix<> uniPoints = geom.eval(uniParams);
    index_t index = findClosestIndex(point, uniPoints);

    if (index < 0)
    {
        std::cout << "Index is out of range." << std::endl;
    }

    gsMatrix<> param(1, 1);
    param = uniParams.col(index);

    if (tolerance < (uniPoints.col(index) - point).norm())
    {
        index_t lowIndex = (index != 0) ? index - 1: index;
        index_t uppIndex = (index != uniParams.cols() - 1) ? index + 1: index;

        gsMatrix<> lowPoint = uniParams.col(lowIndex);
        gsMatrix<> uppPoint = uniParams.col(uppIndex);
        param = binarySearch(point, geom, lowPoint, uppPoint, tolerance, 100);
    }

    return param;
}

real_t percentPointBelowTreshold(const std::vector<real_t>& errors,
                                 const real_t treshold)
{
    int numErrors = errors.size();
    int belowTreshold = 0;
    for (std::size_t i = 0; i != errors.size(); i++)
    {
        if (errors[i] <= treshold)
        {
            belowTreshold++;
        }
    }

    return belowTreshold * 1.0 / numErrors;
}


gsMatrix<> reverseMatrixRows(const gsMatrix<> mat)
{
    const int numRows = mat.rows();
    const int numCols = mat.cols();

    gsMatrix<> reversedMat(numRows, numCols);

    for (index_t i = 0; i != numRows; i++)
    {
        reversedMat(i, 0) = mat(numRows-1-i, 0);
        reversedMat(i, 1) = mat(numRows-1-i, 1);
    }

    return reversedMat;
}

gsMatrix<> reverseMatrixCols(const gsMatrix<> mat)
{
    const int numRows = mat.rows();
    const int numCols = mat.cols();

    gsMatrix<> reversedMat(numRows, numCols);

    for (index_t i = 0; i != numCols; i++)
    {
        reversedMat(0, i) = mat(0, numCols-1-i);
        reversedMat(1, i) = mat(1, numCols-1-i);
    }

    return reversedMat;
}


gsMatrix<> getNormal(const real_t angle)
{
    gsMatrix<> normal(2, 1);
    normal(0, 0) = math::cos(angle);
    normal(1, 0) = math::sin(angle);

    return normal;
}

gsMatrix<> getNormals(const gsMatrix<>& angles)
{
    gsMatrix<> normals(2, angles.cols());
    normals.setZero();

    for (index_t col = 0; col != angles.cols(); col++)
    {
        gsMatrix<> n = getNormal(angles(0, col));
        normals.col(col) = n.col(0);
    }

    return  normals;
}

void plotNormals(const gsMatrix<>& points,
                 const gsMatrix<>& normals,
                 const std::string& file,
                 const real_t length)
{
    gsMultiPatch<> mpNormals;

    index_t dim = points.rows();

    for (index_t col = 0; col != points.cols(); col++)
    {
        gsMatrix<> controlPoints(2, dim);
        controlPoints.row(0) = points.col(col).transpose();

        controlPoints.row(1) = (points.col(col) + length * normals.col(col)).transpose();

        gsKnotVector<> kv(0.0, 1.0, 0, 2);
        mpNormals.addPatch(gsBSpline<>(kv, controlPoints));
    }
    gsInfo << "Plotting normals" << std::endl;
    gsWriteParaview(mpNormals, file);
}


gsGeometry<>::uPtr fitAngles(const gsMatrix<>& parameters,
                        const gsMatrix<>& angles,
                        const int numCtrlPts,
                        const int degree)
{
    gsKnotVector<> kv ( 0, 1, numCtrlPts-degree-1, degree + 1 );
    gsBSplineBasis<> basis( kv );

    gsFitting<> fitting(parameters, angles, basis);

    // usual fitting
    const int size = basis.size();
    gsSparseMatrix<> A(size, size);
    gsMatrix<> B(size, 1);
    fitting.assembleSystem(A, B);

    // Lagrange constraints
    gsMatrix<> C(1, size);
    C.setZero();
    C(0, 0) = 1.0;
    C(0, size - 1) = -1;

    // making bigger matrix AA
    gsMatrix<> A_dense(A);
    gsMatrix<> AA(size + 1, size + 1);
    AA.block(0, 0, size, size) = A_dense;
    AA.block(size, 0, 1, size) = C;
    AA.block(0, size, size, 1) = C.transpose();
    AA(size, size) = 0;

    // making bigger right hand side
    gsMatrix<> BB(size + 1, 1);
    BB.block(0, 0, size, 1) = B;
    BB(size, 0) = -2 * EIGEN_PI;

    // solving a system
    gsSparseMatrix<> AA_sparse = AA.sparseView();
    gsSparseSolver<>::BiCGSTABILUT solver(AA_sparse);

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
    }

    gsMatrix<> x = solver.solve(BB); //toDense()
    gsMatrix<> coefs = x.block(0, 0, size, 1);

    return basis.makeGeometry(coefs);
}
/*
gsGeometry<>* makeCurve(const gsMatrix<>& parameters,
                        const gsMatrix<>& points,
                        const int numCtrlPts,
                        const int degree)
{
    gsKnotVector<> kv ( 0, 1,  numCtrlPts-degree-1, degree+1 );
    gsBSplineBasis<> basis( kv );

    gsFitting<> fitting(parameters, points, basis);

    // usual fitting
    const int size = basis.size();
    gsSparseMatrix<> A(size, size);
    gsMatrix<> B(size, 2);

    fitting.assembleSystem(A, B);
    fitting.applySmoothing(1e-6, A);
    // Lagrange constraints
    gsMatrix<> C(1, size);
    C.setZero();
    C(0, 0) = 1.0;
    C(0, size - 1) = -1;

    // making bigger matrix AA
    gsMatrix<> A_dense(A);
    gsMatrix<> AA(size + 1, size + 1);
    AA.block(0, 0, size, size) = A_dense;
    AA.block(size, 0, 1, size) = C;
    AA.block(0, size, size, 1) = C.transpose();
    AA(size, size) = 0;

    // making bigger right hand side
    gsMatrix<> BB(size + 1, 2);
    BB.block(0, 0, size, 2) = B;
    BB(size, 0) = 0;
    BB(size, 1) = 0;

    // solving a system
    gsSparseMatrix<> AA_sparse = AA.sparseView();
    gsSparseSolver<>::BiCGSTABILUT solver(AA_sparse);

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
    }

    gsMatrix<> x = solver.solve(BB); //toDense()
    gsMatrix<> coefs = x.block(0, 0, size, 2);

    return basis.makeGeometry(coefs);
}
*/
gsMatrix<index_t> generatePermutationMatrix(const index_t num)
{
    gsMatrix<index_t> mat(num, 3);
    mat(0, 0) = num-1;
    mat(0, 1) = 0;
    mat(0, 2) = 1;

    for (index_t i = 1; i != num-1; i++)
    {
        mat(i, 0) = i-1;
        mat(i, 1) = i;
        mat(i, 2) = i+1;
    }
    mat(num-1, 0) = num-2;
    mat(num-1, 1) = num-1;
    mat(num-1, 2) = 0;

    return mat;
}

gsMatrix<> wachspressCoordinates(const gsMatrix<> anchors, const gsMatrix<> points)
{
    const int numVertices = anchors.cols();
    const int dim         = anchors.rows();
    const int numPoints   = points.cols();
    const real_t eps = 1.0e-10;

    gsMatrix<> coords(numVertices, numPoints);
    coords.setZero();

    gsMatrix<> w(1, numVertices);
    w.setZero();

    gsMatrix<> mat(dim, 3); // Two-dimensional case only
    mat.setZero();
    gsMatrix<> mat1(dim, 3);
    mat1.setZero();
    gsMatrix<> mat2(dim, 3);
    mat2.setZero();

    gsMatrix<index_t> indices = generatePermutationMatrix(numVertices);

    for (index_t ii = 0; ii != numPoints; ii++)
    {
    // Calculate barycentric coordinates for each given point
        bool onEdge = false; // Is the point located on the edge of the polygon?
        for (index_t i = 0; i != numVertices; i++)
        {
            gsMatrix<index_t> idx = indices.row(i);
            for (index_t j = 0; j != 3; j++)
            {
                mat.block(0, 0, 2, 3).col(j) = anchors.col(idx(j));
            }
            mat1.block(0, 0, 2, 3).col(0) = points.col(ii);
            mat2.block(0, 0, 2, 3).col(0) = points.col(ii);
            for (index_t j = 1; j != 3; j++)
            {
                mat1.block(0, 0, 2, 3).col(j) = anchors.col(idx(j-1));
                mat2.block(0, 0, 2, 3).col(j) = anchors.col(idx(j));
            }

            real_t area  = triangleArea(mat);
            real_t area1 = triangleArea(mat1);
            real_t area2 = triangleArea(mat2);

            if ( (math::abs(area1) > eps) && (math::abs(area2) > eps) )
            // Point is in the middle of polygon
            {
                if (onEdge == false)
                {
                    w(0, i) = area/(area1*area2);
                }
            }
            else if ( (math::abs(area1) < eps) && (math::abs(area2) > eps) )
            // Point is on the edge of polygon
            {
                if (onEdge == false)
                {
                    w.setZero();
                    onEdge = true;
                }
                w(0, i) = area/(area2);
            }
            else if ( (math::abs(area1) > eps) && (math::abs(area2) < eps) )
            // Point is on the edge of polygon
            {
                if (onEdge == false)
                {
                    w.setZero();
                    onEdge = true;
                }
                w(0, i) = area/(area1);
            }
            else if ( (math::abs(area1) < eps) && (math::abs(area2) < eps) )
            // Point coincides with the vertex of polygon
            {
                // All coordinates are equal to zero,
                // except the one which corresponds to the vertex
                // (which is equal to one)
                w.setZero();
                w(0, i) = 1.0;
                break;
            }
        }

        for (index_t i = 0; i != numVertices; i++)
        {
            coords(i, ii) = w(0, i)/w.sum();
        }
    }

    return coords;
}

real_t triangleArea(const gsMatrix<> points)
{
    const int numPoints = points.cols();
    const int dim = points.rows();

    gsMatrix<> mat(dim+1, numPoints);
    mat.setOnes();  // Ones will stay in the first row
    mat.block(1, 0, dim, numPoints) = points;

    return 0.5*(mat.determinant());
}

gsMatrix<> mapPointsWachspress(const gsMatrix<> templateAnchors, const gsMatrix<> targetAnchors, const gsMatrix<> points)
{
    gsMatrix<> wc = wachspressCoordinates(templateAnchors, points);
    // Applying barycentric coordinates to sampled points and target anchors
    gsMatrix<> mappedPoints = targetAnchors*wc;

    return mappedPoints;
}

std::size_t rectangleIndex(const std::string& inputFile)
{

    const gsMultiPatch<>* multipatch = readMultipatch(inputFile);

    return rectangleIndex(*multipatch);
}

std::size_t rectangleIndex(const gsMultiPatch<>& multipatch)
{
    std::size_t idxRectangle = 0;

    for (std::size_t i = 0; i != multipatch.nPatches(); i++)
    {
        const gsGeometry<>& geom = multipatch[i];
        if (geom.basis().degree(0) == 1)
        {
            idxRectangle = i;
            break;
        }
    }

    return idxRectangle;
}

gsMatrix<> detectVertices(const std::string& inputFile)
{
    const gsMultiPatch<>* mp = readMultipatch(inputFile);
    const gsGeometry<>& geom = (*mp)[rectangleIndex(*mp)];

    return detectVertices(geom);
}

gsMatrix<> detectVertices(const gsGeometry<>& geom)
{
    GISMO_ASSERT(geom.basis().degree(0)==1 && geom.parDim()==1, "Geometry should be polygon");
    gsMatrix<> coefs = geom.coefs();
    int num = coefs.rows();
    int dim = coefs.cols();

    gsMatrix<> vertices;
    vertices.setZero();

    real_t eps = 1.0e-6;

    index_t idx = 0;
    for (index_t i = 2; i != num; i++)
    {
        gsMatrix<> triangle(2, 3);
        triangle.col(0) = coefs.row(i-2);
        triangle.col(1) = coefs.row(i-1);
        triangle.col(2) = coefs.row(i);
        real_t plot = triangleArea(triangle);
        if ( eps < math::abs(plot) )
        {
            vertices.conservativeResize(dim, idx+1);
            vertices.col(idx) = coefs.row(i-1).transpose();
            idx++;
        }
    }
    return vertices;
}

// This functions removes patch (id) from multipatch
gsMultiPatch<>* removePatch(gsMultiPatch<>& mp,
                            const std::size_t id)
{
    gsMultiPatch<>* mp_new = new gsMultiPatch<>();

    for (std::size_t i = 0; i != mp.nPatches(); i++)
    {
        if (i != id)
        {
            const gsGeometry<>& geom = mp[i];
            mp_new->addPatch(geom);
        }
    }

    return mp_new;
}
/*
gsMatrix<> arcLengthParameters(gsMatrix<>& points)
{
     int numPoints = points.cols();
     gsMatrix<> parameters(1, numPoints);
     parameters.setZero();

     parameters(0, 0) = 0.0;
     for (int i = 1; i != numPoints; i++)
     {
         parameters(0, i) = parameters(0, i) + norm(points.col(i) - points.col(i-1));
     }
     parameters = parameters.colwise()/parameters(0, numPoints-1);

     return parameters;
}
*/

std::vector<real_t> getPeriodicKnots(const int deg,
                                     const int numKnots)
{
    real_t h = 1.0 / numKnots;

    std::vector<real_t> knots;

    for (int d = deg; 0 < d; d--)
    {
        knots.push_back(-1 * d * h);
    }

    for (int i = 0; i != numKnots + 1; i++)
    {
        knots.push_back(i * h);
    }

    for (int d = 1; d != deg + 1; d++)
    {
        knots.push_back(1 + d * h);
    }

    return knots;
}

gsKnotVector<> getPeriodicKnotVector(const int deg,
                                     const int numKnots)
{
    std::vector<real_t> knots = getPeriodicKnots(deg, numKnots);

    return gsKnotVector<>(knots, deg);
}

gsMultiPatch<> splittedCurve (const gsGeometry<>& geom,
                              const gsMatrix<>& points)
{
    // B-spline curves
    gsBSpline<> inputGeometry = dynamic_cast< const gsBSpline<>& >(geom);
    gsMatrix<> ctrlPts = inputGeometry.coefs();
    int deg = inputGeometry.degree();
    gsKnotVector<> kv = inputGeometry.knots();

    index_t numPoints = points.cols();
    gsMatrix<> pars(1, numPoints);

    for (int i = 0; i < numPoints; i++)
    {
        pars.col(i) = findParameter(inputGeometry, points.col(i), 1.0e-12);
    }

    // Adding knots
    for (int i = 0; i < numPoints; i++)
    {
        gsBoehm(kv, ctrlPts, pars(i), kv.degree());
    }

    // Extracting curves
    int i_start = 0;
    int i_end = 0;

    gsMatrix<index_t> indices(2, numPoints);
    indices.setZero();

    gsMultiPatch<> mp;
    for (int i = 0; i < numPoints-1; i++)
    {
        // Finding index of knot with value less the parameter
        // corresponding to the CURRENT splitting point
        index_t ind1 = 0;
        for(ind1 = 0; kv.at(ind1) < pars(0, i); ind1++); // Empty loop
        ind1--;
        indices(0, i+1) = ind1;
        //gsInfo << "ind1: " << ind1 << std::endl;
        // If the CURRENT splitting point is THE FIRST in the list, we mark it
        if (i == 0)
        {
            i_start = ind1; // i_start and i_end are unnecessary
            indices(1, 0) = ind1;
            //gsInfo << "indices(1, 0): " << indices(1, 0) << std::endl;
        }

        // Finding index of knot with value less the parameter
        // corresponding to the NEXT splitting point
        index_t ind2 = 0;
        for(ind2 = 0; kv.at(ind2) < pars(0, i+1); ind2++); // Empty loop
        ind2--;
        indices(1, i+1) = ind2; // i_start and i_end are unnecessary
        // If the NEXT splitting point is THE LAST in the list, we mark it
        if (i == numPoints-2)
        {
            i_end = ind2;
            indices(0, 0) = ind2;
            //gsInfo << "indices(1, 0): " << indices(1, 0) << std::endl;
        }

    }

    // Fixing the order of the curves in multipatch
    // Extracting the curve to which the start point belongs
    // and adding it to the multipatch
    gsMatrix<> ctrlPts2(ctrlPts.rows()-i_end-2 + i_start+2, 2);
    ctrlPts2.setZero();
    ctrlPts2.block(0, 0, ctrlPts.rows()-i_end-1, 2) = ctrlPts.block(i_end, 0, ctrlPts.rows()-i_end-1, 2);
    ctrlPts2.block(ctrlPts.rows()-i_end-1, 0, i_start+1, 2) = ctrlPts.block(0, 0, i_start+1, 2);
    gsKnotVector<> kv2(0.0, 1.0, ctrlPts2.rows()-deg-1, deg+1);
    mp.addPatch(gsBSpline<>(kv2, ctrlPts2));

    // Adding to multipatch all other curves
    for (int i = 1; i < numPoints; i++)
    {
        gsMatrix<> ctrlPts3 = ctrlPts.block(indices(0, i), 0, indices(1, i)-indices(0, i)+1, 2);
        gsKnotVector<> kv3(0.0, 1.0, ctrlPts3.rows()-deg-1, deg+1);
        mp.addPatch(gsBSpline<>(kv3, ctrlPts3));
    }

    return mp;

}

gsMultiPatch<> splittedCurveParameters (const gsGeometry<>& geom,
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

    // Extracting curves
    int i_start = 0;
    int i_end = 0;

    gsMatrix<index_t> indices(2, numPoints);
    indices.setZero();

    gsMultiPatch<> mp;
    for (int i = 0; i < numPoints-1; i++)
    {
        // Finding index of knot with value less the parameter
        // corresponding to the CURRENT splitting point
        index_t ind1 = 0;
        for(ind1 = 0; kv.at(ind1) < pars(0, i); ind1++); // Empty loop
        ind1--;
        indices(0, i+1) = ind1;
        // If the CURRENT splitting point is THE FIRST in the list, we mark it
        if (i == 0)
        {
            i_start = ind1; // i_start and i_end are unnecessary
            indices(1, 0) = ind1;
        }

        // Finding index of knot with value less the parameter
        // corresponding to the NEXT splitting point
        index_t ind2 = 0;
        for(ind2 = 0; kv.at(ind2) < pars(0, i+1); ind2++); // Empty loop
        ind2--;
        indices(1, i+1) = ind2; // i_start and i_end are unnecessary
        // If the NEXT splitting point is THE LAST in the list, we mark it
        if (i == numPoints-2)
        {
            i_end = ind2;
            indices(0, 0) = ind2;
        }

    }

    // Fixing the order of the curves in multipatch
    // Extracting the curve to which the start point belongs
    // and adding it to the multipatch
    gsMatrix<> ctrlPts2(ctrlPts.rows()-i_end-2 + i_start+2, 2);
    ctrlPts2.setZero();
    ctrlPts2.block(0, 0, ctrlPts.rows()-i_end-1, 2) = ctrlPts.block(i_end, 0, ctrlPts.rows()-i_end-1, 2);
    ctrlPts2.block(ctrlPts.rows()-i_end-1, 0, i_start+1, 2) = ctrlPts.block(0, 0, i_start+1, 2);
    gsKnotVector<> kv2(0.0, 1.0, ctrlPts2.rows()-deg-1, deg+1);
    mp.addPatch(gsBSpline<>(kv2, ctrlPts2));

    // Adding to multipatch all other curves
    for (int i = 1; i < numPoints; i++)
    {
        gsMatrix<> ctrlPts3 = ctrlPts.block(indices(0, i), 0, indices(1, i)-indices(0, i)+1, 2);
        gsKnotVector<> kv3(0.0, 1.0, ctrlPts3.rows()-deg-1, deg+1);
        mp.addPatch(gsBSpline<>(kv3, ctrlPts3));
    }

    return mp;
}

gsMultiPatch<> splittedRectangle (const gsGeometry<>& geom)
{
    // Identifying corners of the rectangle
    gsMatrix<> points = detectVertices(geom);
    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(points.col(3), points.col(0)));
    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(points.col(0), points.col(1)));
    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(points.col(1), points.col(2)));
    mp.addPatch(gsNurbsCreator<>::BSplineLineSegment(points.col(2), points.col(3)));

    return mp;
}

std::vector< gsMultiPatch<> > splittedTemplate1 (const gsMultiPatch<>& multipatch,
                                                const gsMatrix<>& points)
{
    // This function splits the template of the type 1:
    // Assume the template multipatch consists of rectangle (curve 0 in multipatch)
    // and some other smooth curve, e. g. circle (curve 1 in multipatch)
    const gsGeometry<>& rectangle = multipatch[0];
    const gsGeometry<>& curve = multipatch[1];

    std::vector< gsMultiPatch<> > templateMultipatches;

    // Rectangle is splitted at corners
    gsMultiPatch<> mp_rectangle = splittedRectangle (rectangle);
    templateMultipatches.push_back(mp_rectangle);

    // Smooth curve is splitting according to given points
    gsMultiPatch<> mp_curve = splittedCurve (curve, points);
    templateMultipatches.push_back(mp_curve);

    return templateMultipatches;
}

std::vector< gsMultiPatch<> > splittedTemplate2 (const gsMultiPatch<>& multipatch,
                                                 const gsMatrix<>& points)
{
    // This function splits the template of the type 2:
    // Assume the template multipatch consists of a smooth curve, e. g. circle (curve 0 in multipatch)
    const gsGeometry<>& curve = multipatch[1];

    std::vector< gsMultiPatch<> > templateMultipatches;

    // Smooth curve is splitting according to given points
    gsMultiPatch<> mp_curve = splittedCurve (curve, points);
    templateMultipatches.push_back(mp_curve);

    return templateMultipatches;
}

gsMultiPatch<> fixedBoundary(gsMultiPatch<>& boundary)
{
    gsMultiPatch<>* fBoundary = new gsMultiPatch<>();
    const gsGeometry<>& tmp_face0 = boundary[0];
    const gsGeometry<>& tmp_face1 = boundary[1];
    const gsGeometry<>& tmp_face2 = boundary[2];
    const gsGeometry<>& tmp_face3 = boundary[3];

    gsBSpline<> face0 = dynamic_cast< const gsBSpline<>& >(tmp_face0);
    gsBSpline<> face1 = dynamic_cast< const gsBSpline<>& >(tmp_face1);
    gsBSpline<> face2 = dynamic_cast< const gsBSpline<>& >(tmp_face2);
    gsBSpline<> face3 = dynamic_cast< const gsBSpline<>& >(tmp_face3);

    face0.coef(0)                   = face1.coef(0); //
    face0.coef(face0.coefsSize()-1) = face3.coef(0);
    face2.coef(0)                   = face1.coef(face1.coefsSize()-1); //
    face2.coef(face2.coefsSize()-1) = face3.coef(face3.coefsSize()-1);

    // Making all degrees equal
    //int deg0 = face0.basis().degree();
    //int deg1 = face1.basis().degree();
    //int deg2 = face2.basis().degree();
    //int deg3 = face3.basis().degree();

    //int maxDegree = math::max( math::max( math::max(deg0, deg1), deg2), deg3);
    //face0.degreeElevate(maxDegree);
    //face1.degreeElevate(maxDegree);
    //face2.degreeElevate(maxDegree);
    //face3.degreeElevate(maxDegree);

    fBoundary->addPatch(face0);
    fBoundary->addPatch(face1);
    fBoundary->addPatch(face2);
    fBoundary->addPatch(face3);

    return *fBoundary;
}

gsCoonsPatch<real_t> coonsPatch(gsMultiPatch<>& boundary)
{
    boundary.computeTopology(1e-3);

    if (boundary.nBoundary() != 0)
    {
        // Making the boundary closed
        boundary = fixedBoundary(boundary);
    }

    gsCoonsPatch<real_t> coons(boundary);
    return coons;
}

// Returns a skeleton of the template (template type 1)
std::vector< gsMultiPatch<> > templateSkeletonBoundaries1(const gsMultiPatch<>& mpRectangle,
                                                         const gsMultiPatch<>& mpCurve,
                                                         const gsMultiPatch<>& mpSegmentation)
{
    if ((mpRectangle.nPatches() == mpCurve.nPatches()) && (mpRectangle.nPatches() == mpSegmentation.nPatches()) && (mpCurve.nPatches() == mpSegmentation.nPatches()))
    { }
    else
    {
        gsWarn << "templateSkeletonBoundaries1: Numbers of patches dismatch" << std::endl;
    }

    std::vector< gsMultiPatch<> > boundaries;

    int numPatches = mpRectangle.nPatches();
    for (int i = 0; i < numPatches; i++)
    {
        gsMultiPatch<> boundary;

        if (i == 0)
        {
            boundary.addPatch(mpSegmentation[numPatches-1]);
        }
        else
        {
            boundary.addPatch(mpSegmentation[i-1]);
        }

        boundary.addPatch(mpRectangle[i]);
        boundary.addPatch(mpSegmentation[i]);
        boundary.addPatch(mpCurve[i]);

        boundaries.push_back(boundary);
    }

    return boundaries;
}

// Returns a skeleton of the template (template type 2)
std::vector< gsMultiPatch<> > templateSkeletonBoundaries2(const gsMultiPatch<>& mpCurve,
                                                          const gsMultiPatch<>& mpSegmentation)
{
    std::vector< gsMultiPatch<> > boundaries;

    int numPatches = mpCurve.nPatches() + 1; // = 5

    // Patches around rectangle
    for (int i = 0; i < numPatches-1; i++)
    {
        gsMultiPatch<> boundary;
        // Line from rectangle vertices to curve
        if (i == 0)
        {
            boundary.addPatch(mpSegmentation[numPatches-2]);
        }
        else
        {
            boundary.addPatch(mpSegmentation[i-1]);
        }

        boundary.addPatch(mpCurve[i]);          // Curve segment
        boundary.addPatch(mpSegmentation[i]);   // Line from rectangle to curve
        boundary.addPatch(mpSegmentation[i+4]); // Inner rectange side
        boundaries.push_back(boundary);
    }

    // Rectangle patch
    gsMultiPatch<> boundary;
    for ( int i = 4; i < 8; i++ )
    {
        boundary.addPatch(mpSegmentation[i]);   // Rectangle edge
    }
    boundaries.push_back(boundary);

    return boundaries;
}


template<typename T>
void printStdVector(T vector)
{
    std::cout << vector << " ";
}

void unifyKnots (gsGeometry<>& geom1,
                 gsGeometry<>& geom2)
{
    gsBSpline<> inputGeometry1 = dynamic_cast< const gsBSpline<>& >(geom1);
    gsBSpline<> inputGeometry2 = dynamic_cast< const gsBSpline<>& >(geom2);

    // Degree unification (if necessary)
    int deg1 = inputGeometry1.knots().degree();
    int deg2 = inputGeometry2.knots().degree();
    if (deg1 < deg2)
    {
        inputGeometry1.degreeElevate(deg2-deg1);
    }
    else if (deg2 < deg1)
    {
        inputGeometry2.degreeElevate(deg1-deg2);
    }

    gsKnotVector<> kv1 = inputGeometry1.knots();
    gsKnotVector<> kv2 = inputGeometry2.knots();

    if (kv1.degree() != kv2.degree())
    {
        gsWarn << "Degrees mismatch" << std::endl;
    }

    std::vector<real_t> kv_temp;
    std::set_union(kv1.begin(), kv1.end(),
                   kv2.begin(), kv2.end(),
                   std::back_inserter(kv_temp));
    gsKnotVector<> kv(kv_temp, kv1.degree());

    std::vector<real_t> diff1;
    std::set_difference(kv.begin(), kv.end(), kv1.begin(), kv1.end(),
                        std::inserter(diff1, diff1.begin()));

    std::vector<real_t> diff2;
    std::set_difference(kv.begin(), kv.end(), kv2.begin(), kv2.end(),
                        std::inserter(diff2, diff2.begin()));

    if (diff1.size() != 0)
    {
        inputGeometry1.insertKnots(diff1.begin(), diff1.end());
    }
    if (diff2.size() != 0)
    {
        inputGeometry2.insertKnots(diff2.begin(), diff2.end());
    }

    geom1 = inputGeometry1;
    geom2 = inputGeometry2;

}

gsMultiPatch<> unifyKnotsBoundary (gsMultiPatch<>& boundary)
{
    gsMultiPatch<> uMultipatch;

    gsBSpline<> geom0 = dynamic_cast< const gsBSpline<>& >(boundary[0]);
    gsBSpline<> geom1 = dynamic_cast< const gsBSpline<>& >(boundary[1]);
    gsBSpline<> geom2 = dynamic_cast< const gsBSpline<>& >(boundary[2]);
    gsBSpline<> geom3 = dynamic_cast< const gsBSpline<>& >(boundary[3]);

    unifyKnots(geom0, geom2);
    unifyKnots(geom1, geom3);

    uMultipatch.addPatch(geom0);
    uMultipatch.addPatch(geom1);
    uMultipatch.addPatch(geom2);
    uMultipatch.addPatch(geom3);

    return uMultipatch;
}


std::vector< gsMultiPatch<> > unifyKnotsBoundaries (std::vector< gsMultiPatch<> > boundaries)
{
    //std::vector< gsMultiPatch<> > uMultipatch;
    for (std::size_t i = 0; i != boundaries.size(); i++)
    {
        gsMultiPatch<> mpTmp = boundaries.at(i);
        mpTmp = unifyKnotsBoundary (mpTmp);
        boundaries.at(i) = mpTmp;
    }

    return boundaries;
}


std::vector< gsMultiPatch<> > splittedMultipatch(const gsMultiPatch<>& multipatch,
                                                 const gsFileData<>& points)
// This function splits B-spline curves given in a multipatch
// Input:
// - Multipatch containing curves
// - Filedata containing matrices with coordinates of splitting points
// Output:
// - Splitted curves outputed in separate files
{

    std::size_t numPatches = multipatch.nPatches();
    std::vector< gsMultiPatch<> > multipatches;
    for (std::size_t i = 0; i != numPatches; i++)
    {
        const gsGeometry<>& geom = multipatch[i];
        gsMatrix<>* pts = points.getId< gsMatrix<> >(i).release();
        gsMultiPatch<> mp = splittedCurve (geom, *pts);

        multipatches.push_back(mp);
    }

    return multipatches;
}

/**************************************************/

gsMultiPatch<>* joinMultipatches(const std::string& mpFile1, const std::string& mpFile2)
{
    gsMultiPatch<>* mp1 = readMultipatch(mpFile1);
    gsMultiPatch<>* mp2 = readMultipatch(mpFile2);

    gsMultiPatch<>* mp = joinMultipatches(*mp1, *mp2);

    return mp;
}

gsMultiPatch<>* joinMultipatches(const gsMultiPatch<>& mp1, const gsMultiPatch<>& mp2)
{
    gsMultiPatch<>* mp = new gsMultiPatch<>();
    // Adding mp1 patches to mp
    for (std::size_t i = 0; i != mp1.nPatches(); i++)
    {
        const gsGeometry<>& patch = mp1[i];
        mp->addPatch(patch);
    }
    // Adding mp2 patches to mp
    for (std::size_t i = 0; i != mp2.nPatches(); i++)
    {
        const gsGeometry<>& patch = mp2[i];
        mp->addPatch(patch);
    }
    return mp;
}

template<unsigned d>
std::vector< gsTensorBSpline<d> > coonsPatches(std::vector< gsMultiPatch<> > boundaries)
{
    std::vector< gsTensorBSpline<d, real_t> > cpatches;
    for (std::size_t i = 0; i != boundaries.size(); i++)
    {
        gsMultiPatch<> boundary = boundaries.at(i);
        gsCoonsPatch<real_t> cpatch = coonsPatch(boundary);
        cpatch.compute();

        const gsTensorBSpline<d>& cp = dynamic_cast< const gsTensorBSpline<d>& >(cpatch.result());
        //writeMap(cpatch.result(), i, "CP_map");
        //writeKnots(cpatch.result(), i, "CP_knots");

        cpatches.push_back(cp);
    }

    return cpatches;
}


/**************************************************/
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

gsMatrix<> sumMatrixRows (const gsMatrix<>& mat)
// Sums the rows of the matrix
// For example, if
// A =
// a_11 a_12
// a_21 a_22
// then sumMatrixRows(A) is a matrix with elements
// a_11 + a_21 a12 + a22
// in one row
{
    const index_t numRows = mat.rows();
    const index_t numCols = mat.cols();
    gsMatrix<> result(1, numCols);
    result.setZero();

    for (index_t j = 0; j < numCols; j++)
    {
        for (index_t i = 0; i < numRows; i++)
        {
            result(0, j) += mat(i, j);
        }
    }

    return result;
}

gsMatrix<> sumMatrixCols (const gsMatrix<>& mat)
// Sums the columns of the matrix
// For example, if
// A =
// a_11 a_12
// a_21 a_22
// then sumMatrixCols(A) is a matrix with elements
// a_11 + a_11
// a_21 + a22
// in one column
{
    const index_t numRows = mat.rows();
    const index_t numCols = mat.cols();
    gsMatrix<> result(numRows, 1);
    result.setZero();

    for (index_t i = 0; i < numRows; i++)
    {
        for (index_t j = 0; j < numCols; j++)
        {
            result(0, i) += mat(i, j);
        }
    }

    return result;
}

gsBSpline<> constructEllipse(const real_t x_center,
                             const real_t y_center,
                             const real_t r1,
                             const real_t r2,
                             const int numCtrlPts,
                             const int degree)
{

    gsMatrix<> ctrlPts(numCtrlPts, 2);
    ctrlPts.setZero();

    for (int i = 0; i != numCtrlPts; i++)
    {
        const real_t param = (2.0 * EIGEN_PI * i) / (numCtrlPts - 1);
        ctrlPts(i, 0) = r1 * math::cos(param) + x_center;
        ctrlPts(i, 1) = r2 * math::sin(param) + y_center;
    }
    make_c0(ctrlPts);

    gsKnotVector<> KV(0.0, 1.0, numCtrlPts - degree - 1, degree + 1);

    // Constructing ellipse as a B-Spline
    gsBSpline<> ellipsis(KV, ctrlPts);

    return ellipsis;
}

gsMultiPatch<> mapMultipatch(const gsGeometry<>& map,
                             const gsMultiPatch<>& mp)
{
    gsMultiPatch<> newMp(mp);

    for (std::size_t i = 0; i != newMp.nPatches(); i++)
    {
        gsGeometry<>& geom = newMp.patch(i);
        gsMatrix<> coefs = geom.coefs();
        coefs.transposeInPlace();

        gsMatrix<> newCoefs = map.eval(coefs);
        newCoefs.transposeInPlace();
        geom.setCoefs(newCoefs);
    }

    return newMp;
}

/*
void mapPoints(const std::string& inputFile,
               const gsGeometry<>& map,
               const std::string& outputPrefix = "result_")
{
    gsMatrix<> pts;
    importMatrixFromASCII (inputFile, pts);
    mapPoints(pts, map, outputPrefix);
}
*/

void mapPoints(const gsMatrix<> pts,
               const gsGeometry<>& map,
               const std::string& outputPrefix = "result_")
{
    gsMatrix<> mappedPts = map.eval(pts);

    gsFileData<> fd;
    fd << pts;
    fd.dump(outputPrefix);

    gsWriteParaviewPoints(mappedPts, outputPrefix);
}


gsMatrix<> stdVectorToMatrix(const std::vector< gsMatrix<> >& vector)
{
    index_t numCols = 0;
    index_t numRows = vector.at(0).rows();
    for (std::size_t i = 0; i != vector.size(); i++)
    {
        numCols += vector.at(i).cols();
        if (numRows != vector.at(i).rows())
        {
            gsInfo <<  "All matrices must have the same number of rows";
        }
    }
    gsMatrix<> mat (numRows, numCols);

    index_t ind = 0;
    for (std::size_t i = 0; i != vector.size(); i++)
    {
        const gsMatrix<>& vector_element = vector.at(i);
        const std::size_t numColsi = vector_element.cols();
        mat.block(0, ind, numRows, numColsi) = vector_element.block(0, 0, numRows, numColsi);
        ind += numColsi;
    }

    return mat;
}

/*
std::vector<> MatrixToStdVectorRowwise(const gsMatrix<> & mat)
{
    index_t size = mat.rows();
    std::vector<> vec;
    for (std::size_t i = 0; i != size; i++)
    {
        vec.push_back(mat.row(i));
    }

    return vec;
}

std::vector<> MatrixToStdVectorColwise(const gsMatrix<> & mat)
{
    index_t size = mat.cols();
    std::vector<> vec ();
    for (std::size_t i = 0; i != size; i++)
    {
        vec.push_back(mat.col(i));
    }

    return vec;
}
*/

real_t detJacobian(const gsGeometry<>& geom,
                   const gsMatrix<> point)
{
    return geom.jacobian(point).determinant();
}

// Returns minimal absolute value of Jacobian
real_t minDetJacobian(const gsGeometry<>& geom,
                     const int points = 100)
{
    gsMatrix<> para  = geom.support();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> pts = uniformPointGrid(c0, c1, points);

    // calculating the determinant of the Jacobian
    real_t minDetJacobian = detJacobian(geom, pts.col(1));
    for (int col = 1; col < pts.cols(); col++)
    {
        const real_t det = math::abs(detJacobian(geom, pts.col(col)));
        if (det < minDetJacobian)
        {
            minDetJacobian = det;
        }
    }
    return minDetJacobian;
}

void checkJacobianDeterminant(const gsGeometry<>& geom,
                              const int points = 100,
                              const bool savePoints = false,
                              const std::string& output = "",
                              const int number = 0)
{
    gsInfo << "Checking Jacobian determinant ..." << "\n";

    gsMatrix<> para  = geom.support();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> pts = uniformPointGrid(c0, c1, points);

    gsMatrix<> plus(pts.rows(), pts.cols());
    plus.setZero();
    gsMatrix<> minus(pts.rows(), pts.cols());
    minus.setZero();
    gsMatrix<> zero(pts.rows(), pts.cols());
    zero.setZero();

    int plusCounter = 0;
    int minusCounter = 0;
    int zeroCounter = 0;

    // calculating the determinant of the Jacobian
    for (int col = 0; col < pts.cols(); col++)
    {
        const real_t determinant = geom.jacobian(pts.col(col)).determinant();
        gsMatrix<> tmp = geom.eval(pts.col(col));

        if (determinant < 0)
        {
            minus.col(minusCounter) = tmp.col(0);
            minusCounter++;
        }
        else if (determinant > 0)
        {
            plus.col(plusCounter) = tmp.col(0);
            plusCounter++;
        }
        else
        {
            zero.col(zeroCounter) = tmp.col(0);
            zeroCounter++;
        }
    }

    gsInfo << "Number of points with \n"
              << "  - positive sign: " << plusCounter << "\n"
              << "  - negative sign: " << minusCounter << "\n"
              << "  - zero sign:     " << zeroCounter << "\n"
              << "  - perc. of neg.: " << 100.0*minusCounter/points << "\n" << "\n";

    if (savePoints)
    {
        if (0 < plusCounter)
        {
            std::string out = output + "_PositivePoints_" + util::to_string(number);
            gsWriteParaviewPoints(plus, out);
        }

        if (0 < minusCounter)
        {
            std::string out = output + "_NegativePoints_" + util::to_string(number);
            gsWriteParaviewPoints(minus, out);
        }

        if (0 < zeroCounter)
        {
            std::string out = output + "_ZeroPoints_" + util::to_string(number);
            gsWriteParaviewPoints(zero, out);
        }
    }
}

gsMultiPatch<>* extrudeMultipatch(const gsMultiPatch<>& mp, const real_t length)
{
    gsMultiPatch<>* mp_new = NULL;
    mp_new = new gsMultiPatch<>();

    gsKnotVector<> knots_v(0.0, 1.0, 0, 2);

    for (std::size_t i = 0; i != mp.nPatches(); i++)
    {
        gsGeometry<>& patch = mp.patch(i);
        gsBSpline<>& curve = dynamic_cast< gsBSpline<>& >(patch);

        int N = curve.coefsSize();
        gsMatrix<> C_0(N, 3);
        C_0.setZero();
        gsMatrix<> C_1(N, 3);
        C_1.setConstant(length);
        gsMatrix<> C(2*N, 3);
        C.setZero();

        C_0.block(0, 0, N, 2) = curve.coefs();
        C_1.block(0, 0, N, 2) = curve.coefs();
        C.block(0, 0, N, 3) = C_0;
        C.block(N, 0, N, 3) = C_1;

        gsKnotVector<> knots_u = curve.knots();

        gsTensorBSplineBasis<2> BSplineBasis( knots_u, knots_v );
        gsTensorBSpline<2> BSpline( BSplineBasis, C );

        mp_new->addPatch(gsTensorBSpline<2>(BSpline));
    }

    return mp_new;
}

/**************************************************/
template<unsigned d>
void mapCoonsPatches(const gsGeometry<>* map,
                     const std::vector< gsTensorBSpline<d, real_t> > cpatches,
                     const int numPts)
{
    gsMultiPatch<> mpCoonsPatches;
    gsMatrix<> parameters = uniformParameters<d+1>(0.0, 1.0, numPts);
    gsFileData<> fdPointsParameters;

        gsInfo << "mapCoonsPatches / cpatches.size(): " << cpatches.size() << std::endl;
    for (std::size_t i = 0; i != cpatches.size(); i++)
    {
        gsInfo << "mapCoonsPatches / patch: " << i << std::endl;
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

/**************************************************/
std::vector< gsTensorBSpline<2, real_t> > getCoonsPatchesPlanar(const std::vector< gsMultiPatch<> > boundaries)
{
    std::vector< gsMultiPatch<> > bound = unifyKnotsBoundaries(boundaries);
    std::vector< gsTensorBSpline<2, real_t> > cpatches = coonsPatches<2>(bound);
    return cpatches;
}

/**************************************************/
std::vector< gsTensorBSpline<3, real_t> > getCoonsPatchesVolumetric(const std::vector< gsMultiPatch<> > boundaries)
{
    std::vector< gsMultiPatch<> > bound = boundaries;
    std::vector< gsTensorBSpline<3, real_t> > cpatches = coonsPatches<3>(bound);
    return cpatches;
}

/*
void adaptiveFitting()
// Input:
// - vector of gsCoonsPatch objects
// - the map from template geometry to taget geometry (constructed with gsTemplateMappingBarycentric)
// -
{

    real_t lambda = 1e-6;
    int numIterations = 7;
    int extension = 1;
    int degree = 3;
    int interiorKnots = 2;


    std::vector< gsMatrix<> > xy(cpatches.size());
    gsMatrix<> uv = uniformParameters<3>(0.0, 1.0, 100);
    for (std::size_t i = 0; i != cpatches.size(); i++)
    {
        xy.at(i) = (cpatches.at(i)).eval(uv);
    }

    gsMatrix<> bbox = boundingBox<d>(uv, 0.1);

    std::vector< gsKnotVector<> > KV;
    for(unsigned i = 0; i != d; i++)
    {
        KV.push_back(gsKnotVector<> ( bbox(i, 0), bbox(i, 1), interiorKnots, degree + 1 ));
    }

    gsTensorBSplineBasis<d> tbasis( KV );
    gsTHBSplineBasis<d> thb(tbasis);
    std::vector<unsigned> ext(d, extension);
    gsHFitting<d, real_t> fitting( paruv, xy, thb, 1.0, ext, lambda);

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

}
*/


} // namespace gismo

