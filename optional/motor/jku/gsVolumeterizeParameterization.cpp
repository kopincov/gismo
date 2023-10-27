/** @file gsVolumeterizeParameterization.cpp

    @brief Volumeterizes planar multi-patch parameterization

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Sajavicius
*/

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsParameterLines.h"

using namespace gismo;

gsMultiPatch<>* extrudeParameterization(const gsMultiPatch<>& mp,
                                        const real_t length);
gsMultiPatch<>* rotateParameterization(const gsMultiPatch<>& mp,
                                       const real_t dist,
                                       const real_t angle,
                                       const gsVector<>& center,
                                       const int numKnots,
                                       const int degree);

int main(int argc, char *argv[])
{
    // File containing planar parameterization
    std::string paramFile = MOTOR_DATA_DIR "geo/unitSquare.xml";
    // Distance
    real_t dist = 1.0;
    // Total wrap angle which planar parameterization makes in dist
    real_t angle = 0;
    // Coordinates of the rotation center
    real_t c1 = 0.5;
    real_t c2 = 0.5;
    // Number of knots in z-direction
    int numKnots = 5;
    // Degree in z-direction
    int degree = 2;

    int npts = 1000;
    // Prefix for all output files
    std::string outputPrefix("results");

    gsCmdLine cmd("Volumeterizes planar multi-patch parameterization");

    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);
    cmd.addInt("n", "npts", "Number of sample points (for ParaView visualisation)", npts);
    cmd.addInt("d", "degree", "", degree);
    cmd.addInt("k", "numKnots", "", numKnots);
    cmd.addReal("", "c2", "Coordinate of the rotating center (y)", c2);
    cmd.addReal("", "c1", "Coordinate of the rotating center (x)", c1);
    cmd.addReal("", "angle", "Total wrap angle", angle);
    cmd.addReal("", "dist", "Distance", dist);
    cmd.addString("P", "paramFile", "File containing planar parameterization", paramFile);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "paramFile:     " << paramFile << "\n"
           << "dist:          " << dist << "\n"
           << "angle:         " << angle << "\n"
           << "c1:            " << c1 << "\n"
           << "c2:            " << c2 << "\n"
           << "numKnots:      " << numKnots << "\n"
           << "degree:        " << degree << "\n"
           << "npts:          " << npts << "\n"
           << "outputPrefix:  " << outputPrefix << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsMatrix<> center(2, 1);
    center << c1, c2;

    gsFileData<> fd;
    fd.read( paramFile );
    if( !fd.hasAny< gsMultiPatch<> >() )
    {
        gsInfo << "Convert " << paramFile << " to gsMultipatch format (see gsGeometryToMultipatch)" << std::endl;
        return -1;
    }

    const gsMultiPatch<>* paramM = readMultipatch(paramFile);
    const gsMultiPatch<>& paramMultipatch = *paramM;

    gsMultiPatch<>* volParamMultipatch = NULL;
    volParamMultipatch = new gsMultiPatch<>();

    if (angle !=0)
    // Rotation
    {
        volParamMultipatch = rotateParameterization(paramMultipatch, dist, angle,
                                                    center, numKnots, degree);
    }
    else
    // Extrusion
    {
        volParamMultipatch = extrudeParameterization(paramMultipatch, dist);
    }

    // saving results
    writeMultipatch(*volParamMultipatch, outputPrefix + "_mp", npts);
    writeParameterLines(*volParamMultipatch, outputPrefix + "_param_lines", 500, 7);

    return 0;
}


gsMultiPatch<>* extrudeParameterization(const gsMultiPatch<>& mp, const real_t length)
{
    gsMultiPatch<>* mp_new = NULL;
    mp_new = new gsMultiPatch<>();

    const gsKnotVector<> knots_2(0.0, 1.0, 0, 2);

    for (std::size_t i = 0; i != mp.nPatches(); i++)
    {
        gsGeometry<>& patchGeom = mp.patch(i);
        gsTensorBSpline<2>& patch = dynamic_cast< gsTensorBSpline<2>& >(patchGeom);

        const gsMatrix<> coefs = patch.coefs();
        const int N = patch.coefsSize();
        // The bottom layer of control points
        gsMatrix<> C_0(N, 3);
        C_0.setZero();
        // The top layer of control points
        gsMatrix<> C_1(N, 3);
        C_1.setConstant(length);

        // Control points of volumetric patch
        gsMatrix<> C(2*N, 3);
        C.setZero();

        C_0.block(0, 0, N, 2) = coefs; // The values of the third coordinate stays equal to 0
        C_1.block(0, 0, N, 2) = coefs; // The values of the third coordinate stays equal to "length"

        C.block(0, 0, N, 3) = C_0;
        C.block(N, 0, N, 3) = C_1;

        // Creating volumetric patch
        const gsKnotVector<> knots_0 = patch.basis().knots(0);
        const gsKnotVector<> knots_1 = patch.basis().knots(1);

        const gsTensorBSplineBasis<3> BSplineBasis( knots_0, knots_1, knots_2 );
        const gsTensorBSpline<3> BSpline( BSplineBasis, C );

        mp_new->addPatch(gsTensorBSpline<3>(BSpline));
    }

    return mp_new;
}


gsMultiPatch<>* rotateParameterization(const gsMultiPatch<>& mp,
                                       const real_t dist,
                                       const real_t angle,
                                       const gsVector<>& center,
                                       const int numKnots,
                                       const int degree)
{
    gsMultiPatch<>* mp_new = NULL;
    mp_new = new gsMultiPatch<>();

    const gsKnotVector<> knots_2(0.0, 1.0, numKnots, degree+1);
    const int numLayers = knots_2.size()-degree-1; // Number of control point layers in z-direction
    const real_t h = 1.0*dist / (numLayers-1); // Distance between two successive control point layers
    const real_t h_angle = 1.0*angle / (numLayers-1); // Rotation angle between two layers

    for (std::size_t i = 0; i != mp.nPatches(); i++)
    {
        gsGeometry<>& patchGeom = mp.patch(i);
        gsTensorBSpline<2>& patch = dynamic_cast< gsTensorBSpline<2>& >(patchGeom);

        const gsMatrix<> coefs = patch.coefs();
        const int N = patch.coefsSize();
        // The bottom layer of control points
        gsMatrix<> C_0(N, 3);
        C_0.setZero();
        // The intermediate layers of control points
        gsMatrix<> C_1(N, 3);
        gsMatrix<> C(numLayers*N, 3);
        // Control points of volumetric patch
        C.setZero();

        C_0.block(0, 0, N, 2) = coefs; // The values of the third coordinate stays equal to 0
        C.block(0, 0, N, 3) = C_0;

        // Calculating intermediate layers of control points
        for (int l = 1; l < numLayers; l++)
        {
            // Rotating control point by angle "l*h_angle"
            gsMatrix<> rotatedCoefs = rotateVectors(coefs.transpose(), l*h_angle, center).transpose();
            C_1.setConstant(l*h);
            C_1.block(0, 0, N, 2) = rotatedCoefs;  // The values of the third coordinate stays equal to l*h
            C.block(l*N, 0, N, 3) = C_1;
        }

        // Creating volumetric patch
        const gsKnotVector<> knots_0 = patch.basis().knots(0);
        const gsKnotVector<> knots_1 = patch.basis().knots(1);

        const gsTensorBSplineBasis<3> BSplineBasis( knots_0, knots_1, knots_2 );
        const gsTensorBSpline<3> BSpline( BSplineBasis, C );

        mp_new->addPatch(gsTensorBSpline<3>(BSpline));

    }

    return mp_new;
}
