/** @file gsSmoothCurveOffset.cpp

    @brief Construction of inward and outward offsets for smooth curve
    (WORK IN PROGRESS...)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius, J. Speh
*/

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

gsMatrix<> tangentVectors(const gsMatrix<>& vects);
real_t func_f(const real_t arg,
              const real_t slope,
              const real_t kappa_max,
              const real_t f_lower_threshold,
              const real_t f_upper_threshold,
              const bool inwardOffset);
real_t func_f(const real_t arg,
              const real_t slope,
              const real_t kappa_max,
              const real_t f_min,
              const bool inwardOffset);
real_t replacement_f(const real_t arg,
                     const real_t k1,
                     const real_t k2,
                     const real_t a0,
                     const real_t b0,
                     const real_t a1,
                     const real_t b1,
                     const real_t a2,
                     const real_t b2,
                     const bool inwardOffset);
gsMatrix<> replaceControlPoints(const gsMatrix<>& ctrlPts,
                                const gsMatrix<>& nrmlVects,
                                const gsVector<>& kappa,
                                const real_t slope,
                                const real_t kappa_max,
                                const real_t f_lower_threshold,
                                const real_t f_upper_threshold,
                                const bool inwardOffset);
gsMatrix<> replaceControlPoints(const gsMatrix<>& ctrlPts,
                                const gsMatrix<>& nrmlVects,
                                const gsVector<>& kappa,
                                const real_t slope,
                                const real_t kappa_max,
                                const real_t f_min,
                                const bool inwardOffset);
gsMatrix<> replaceControlPoints(const gsMatrix<>& ctrlPts,
                                const gsMatrix<>& nrmlVects,
                                const gsVector<>& kappa,
                                const real_t k1,
                                const real_t k2,
                                const real_t a0,
                                const real_t b0,
                                const real_t a1,
                                const real_t b1,
                                const real_t a2,
                                const real_t b2,
                                const bool inwardOffset);

int main(int argc, char *argv[])
{
    // This code should be revised (problems with cdash)
    return 0;

//    // Options with default values
//    std::string inputFile(MOTOR_DATA_DIR "jku/target_input.xml");
//    real_t k1 = 0.0;
//    real_t k2 = 0.0;
//    real_t a0 = 0.0;
//    real_t b0 = 0.0;
//    real_t a1 = 0.0;
//    real_t b1 = 0.0;
//    real_t a2 = 0.0;
//    real_t b2 = 0.0;
//    int refine = 0;
//    int numIterations = 1;
//    bool inwardOffset = false;

//    gsCmdLine cmd("Construction of offset for smooth curve");
//    cmd.addString("f", "file", "file", inputFile);
//    cmd.addReal("A", "k1", "...", k1);
//    cmd.addReal("B", "k2", "...", k2);
//    cmd.addReal("C", "a0", "...", a0);
//    cmd.addReal("D", "b0", "...", b0);
//    cmd.addReal("E", "a1", "...", a1);
//    cmd.addReal("F", "b1", "...", b1);
//    cmd.addReal("G", "a2", "...", a2);
//    cmd.addReal("H", "b2", "...", b2);

//    cmd.addInt("r", "ref", "Uniformly refines the input curve.", refine);
//    cmd.addInt("i", "numIterations", "Number of iterations (offset constructing)", numIterations);
//    cmd.addSwitch("inward", "Constructs iward offset, if this variable is true", inwardOffset);

//    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

//    gsInfo << "Input arguments:\n"
//           << "Input file:        " << inputFile << "\n"
//           << "k1:                " << k1 << "\n"
//           << "k2:                " << k2 << "\n"
//           << "a0:                " << a0 << "\n"
//           << "b0:                " << b0 << "\n"
//           << "a1:                " << a1 << "\n"
//           << "b1:                " << b1 << "\n"
//           << "a2:                " << a2 << "\n"
//           << "b2:                " << b2 << "\n"
//           << "refine:            " << refine << "\n"
//           << "numIterations:     " << numIterations << "\n"
//           << "inwardOffset:       " << inwardOffset << "\n"
//           << "--------------------------------------------------\n" << std::endl;
// /*
//    if (inputFile == "")
//    {
//        gsInfo << "Didn't specify the input file." << std::endl;
//        return 0;
//    }
// */

//    size_t d = 2;

//    gsFileData<> fd(inputFile);
//    gsGeometry<>* p_geometry = fd.getAnyFirst< gsGeometry<> >().release();
//    gsGeometry<>& geometry = *p_geometry;
//    if (geometry.parDim() == 1)
//        d = 2;
//    else if (geometry.parDim() == 2)
//        d = 3;

//    gsBSpline<> inputGeometry = dynamic_cast< const gsBSpline<>& >(geometry);

//    if (0 < refine)
//    {
//        for (int r = 0; r != refine; r++)
//        {
//            gsKnotVector<> KV = inputGeometry.knots();

//            gsKnotVector<>::knotContainer newKnots;
//            KV.getUniformRefinementKnots(1, newKnots);

//            inputGeometry.insertKnots(newKnots.begin(), newKnots.end());
//        }
//    }

//    gsMatrix<> inCtrlPts = inputGeometry.coefs();
//    gsKnotVector<> KV = inputGeometry.knots();
//    size_t num = inCtrlPts.rows();

//    gsMatrix<> replacedCtrlPts(num, d);
//    replacedCtrlPts.setZero();

//    // Constructing curve as a B-Spline
//    gsBSpline<> bsc(KV, inCtrlPts);

//    gsInfo << "Writing to ParaView..." << std::endl;
//    gsWrite(bsc, "bsplinecurve");
//    gsWriteParaview(bsc, "bsplinecurve", 100);
//    writeControlPoints(inCtrlPts, "bsplinecurveCtrlPts");

//    gsMatrix<> tangentVects(num, d);  // Normalized tangent vectors
//    tangentVects.setZero();

//    gsMatrix<> nrmlVects(num, d);  // Normal vectors
//    nrmlVects.setZero();

//    gsVector<> kappa(num);            // Curvatures
//    kappa.setZero();

//    for (int iter = 0; iter != numIterations; iter++)
//    {
//        gsMatrix<> CtrlPts = inputGeometry.coefs();
//        gsKnotVector<> KV = inputGeometry.knots();

//        tangentVects = tangentVectors(inCtrlPts);

//        if(inwardOffset)
//            nrmlVects  = rotateVectors(tangentVects,  EIGEN_PI/2.0);
//        else // Outer offset
//        {
//            gsInfo << "dim" << tangentVects.rows() << std::endl;
//            nrmlVects = rotateVectors(tangentVects, -EIGEN_PI/2.0);
//        }
//        kappa = curvatures_ctrlpts(CtrlPts, nrmlVects);

//        //replacedCtrlPts = replaceControlPoints(ctrlPts, nrmlVects, kappa, slope, kappa_max, f_lower_threshold, f_upper_threshold, inwardOffset);
//        //replacedCtrlPts = replaceControlPoints(CtrlPts, nrmlVects, kappa, slope, kappa_max, f_min, inwardOffset);
//        replacedCtrlPts = replaceControlPoints(CtrlPts, nrmlVects, kappa, k1, k2, a0, b0, a1, b1, a2, b2, inwardOffset);

//        gsBSpline<> offset(KV, replacedCtrlPts);

//        gsInfo << "Writing to ParaView ... ..." << std::endl;
//        std::string prefix("");
//        if(inwardOffset)
//            prefix = "in";
//        else
//            prefix = "out";
//        gsWrite(offset, prefix + "Offsetcurve");
//        gsWriteParaview(offset, prefix + "Offsetcurve", 100);
//        writeControlPoints(replacedCtrlPts, prefix + "OffsetcurveCtrlPts");
//        writeDisplacements(inCtrlPts, replacedCtrlPts, "displacements");

//        gsInfo << "Min. curv." << kappa.minCoeff() << std::endl;
//        gsInfo << "Max. curv." << kappa.maxCoeff() << std::endl;

//        inputGeometry = offset;
//    }

//    return system("paraview bsplinecurve.pvd &");
}

gsMatrix<> tangentVectors(const gsMatrix<>& vects)
{
    size_t d = vects.cols();
    size_t num = vects.rows();
    gsMatrix<> tangentVects(num, d);
    for(std::size_t i = 0; i != num-1; i++)
    {
        tangentVects.row(i) = vects.row(i+1) - vects.row(i);
    }
    tangentVects.row(num-1) = tangentVects.row(0);
    normalize(tangentVects);

    return tangentVects;
}

real_t func_f(const real_t arg,
              const real_t slope,
              const real_t kappa_max,
              const real_t f_lower_threshold,
              const real_t f_upper_threshold,
              const bool inwardOffset)
{
    if (inwardOffset) // Inward offset (therefore, we have '-' before min)
        return -math::min(math::abs(arg) <= kappa_max ?  slope*sgn(arg)*arg + (f_lower_threshold - slope*kappa_max) : f_lower_threshold, f_upper_threshold);
    else // Outward offset
        return math::max(math::abs(arg) <= kappa_max ? -slope*sgn(arg)*arg + (f_upper_threshold + slope*kappa_max) : f_upper_threshold, f_lower_threshold);
}

real_t func_f(const real_t arg,
              const real_t slope,
              const real_t kappa_max,
              const real_t f_min,
              const bool inwardOffset)
{
    if (inwardOffset) // Inward offset (therefore, we have '-' before min)
        return 0.0;
    else // Outer offset
        return arg <= kappa_max ? 0.005  : 0.01;
        //return arg <= kappa_max ? -slope*arg  + f_min + slope*kappa_max : f_min;
}

real_t replacement_f(const real_t arg,
                     const real_t k1,
                     const real_t k2,
                     const real_t a0,
                     const real_t b0,
                     const real_t a1,
                     const real_t b1,
                     const real_t a2,
                     const real_t b2,
                     const bool inwardOffset)
{
    if(arg < k1)
    {
        return a0*arg + b0;
    }
    else if(arg < k2)
    {
        return a1*arg + b1;
    }
    else
    {
        return a2*arg + b2;
    }
}

gsMatrix<> replaceControlPoints(const gsMatrix<>& ctrlPts,
                                const gsMatrix<>& nrmlVects,
                                const gsVector<>& kappa,
                                const real_t slope,
                                const real_t kappa_max,
                                const real_t f_lower_threshold,
                                const real_t f_upper_threshold,
                                const bool inwardOffset)
{
    //GISMO_ENSURE(ctrlPts.dim() == nrmlVects.dim(), "Matrices of control points and normal vectors must have the same dimension");
    size_t d = ctrlPts.cols();
    size_t num = ctrlPts.rows();
    gsMatrix<> replacedPts(num, d); // Matrix for replaced control points
    replacedPts.setZero();

    for(std::size_t i = 0; i != num-1; i++)
    {
        replacedPts.row(i) = ctrlPts.row(i) + func_f(kappa(i), slope, kappa_max, f_lower_threshold, f_upper_threshold, inwardOffset)*nrmlVects.row(i);
    }
    replacedPts.row(num-1) = replacedPts.row(0);

    return replacedPts;
}

gsMatrix<> replaceControlPoints(const gsMatrix<>& ctrlPts,
                                const gsMatrix<>& nrmlVects,
                                const gsVector<>& kappa,
                                const real_t slope,
                                const real_t kappa_max,
                                const real_t f_min,
                                const bool inwardOffset)
{
    //GISMO_ENSURE(ctrlPts.dim() == nrmlVects.dim(), "Matrices of control points and normal vectors must have the same dimension");
    size_t d = ctrlPts.cols();
    size_t num = ctrlPts.rows();
    gsMatrix<> replacedPts(num, d); // Matrix for replaced control points
    replacedPts.setZero();

    for(std::size_t i = 0; i != num-1; i++)
    {
        replacedPts.row(i) = ctrlPts.row(i) + func_f(kappa(i), slope, kappa_max, f_min, inwardOffset)*nrmlVects.row(i);
    }
    replacedPts.row(num-1) = replacedPts.row(0);

    return replacedPts;

}

gsMatrix<> replaceControlPoints(const gsMatrix<>& ctrlPts,
                                const gsMatrix<>& nrmlVects,
                                const gsVector<>& kappa,
                                const real_t k1,
                                const real_t k2,
                                const real_t a0,
                                const real_t b0,
                                const real_t a1,
                                const real_t b1,
                                const real_t a2,
                                const real_t b2,
                                const bool inwardOffset)
{
    //GISMO_ENSURE(ctrlPts.dim() == nrmlVects.dim(), "Matrices of control points and normal vectors must have the same dimension");
    size_t d = ctrlPts.cols();
    size_t num = ctrlPts.rows();
    gsMatrix<> replacedPts(num, d); // Matrix for replaced control points
    replacedPts.setZero();

    for(std::size_t i = 0; i != num-1; i++)
    {

        const real_t displ = replacement_f(kappa(i), k1, k2, a0, b0, a1, b1, a2, b2, inwardOffset);
        replacedPts.row(i) = ctrlPts.row(i) + displ * nrmlVects.row(i);
    }
    
    replacedPts.row(num-1) = replacedPts.row(0);

    return replacedPts;
}

