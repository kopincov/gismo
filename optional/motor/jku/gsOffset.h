/** @file gsOffset.h

    @brief ...

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Sajavicius, J. Speh
*/

#pragma once

namespace gismo
{

gsMultiPatch<> *readMultipatch(const std::string& inputFile);
void writeMultipatch(const gsMultiPatch<>& multipatch,
                     const std::string& outputFile);
gsGeometry<>* readGeometry(const std::string& inputFile);
gsMatrix<> rotateVectors(const gsMatrix<>& vects,
                         const real_t rad);
void writeControlPoints(const gsMatrix<>& ctrlPts,
                        const std::string& filename);
void make_g1(gsMatrix<>& ctrlPts);
void make_g1(gsGeometry<>& geom);

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
    gsMatrix<> parameters = uniformParameters<d>(0.0, 1.0, numPoints); // different
    gsMatrix<> points = geom.eval(parameters);
    return points;
}

template <unsigned d>
gsMatrix<> samplePoints(const gsMultiPatch<>& mp,
                        const int numPoints = 100)
{
    gsMatrix<> parameters = uniformParameters<d>(0.0, 1.0, numPoints); // different
    gsMatrix<> points(d, parameters.cols() * mp.nPatches());
    for (std::size_t i = 0; i != mp.nPatches(); i++)
    {
        const gsGeometry<>& geom = mp[i];
        points.block(0, i * parameters.cols(), d, parameters.cols()) = samplePoints<d>(geom, numPoints);
    }

    return points;
}

gsMultiPatch<>* readMultipatch(const std::string& inputFile)
{
    gsFileData<> fd(inputFile);
    gsMultiPatch<>* mp = NULL;

    if ( fd.has< gsMultiPatch<> >() )
    {
        mp = fd.getAnyFirst< gsMultiPatch<> >().release();
    }
    else
    {
        gsInfo << "There is no multipatch in the file " << inputFile << std::endl;
    }

    return mp;
}

void writeMultipatch(const gsMultiPatch<>& multipatch,
                     const std::string& outputFile)
{
    gsFileData<> fd;
    gsInfo << "Writing " << outputFile << ".xml" << std::endl;
    fd << multipatch;
    fd.dump(outputFile);
    gsInfo << "Writing " << outputFile << ".pvd" << std::endl;
    gsWriteParaview(multipatch, outputFile);
}

gsGeometry<>* readGeometry(const std::string& inputFile)
{
    gsFileData<> fd(inputFile);
    gsGeometry<>* geom = NULL;

    if ( fd.has< gsGeometry<> >() )
    {
        geom = fd.getAnyFirst< gsGeometry<> >().release();
    }
    else
    {
        gsInfo << "There is no geometry in the file " << inputFile << std::endl;
    }

    return geom;
}

gsMatrix<> rotateVectors(const gsMatrix<>& vects,
                         const real_t rad)
{
    GISMO_ENSURE(vects.rows() == 2, "Only 2D vectors can be rotated");
    size_t d = vects.rows();
    gsMatrix<> transformM(d, d);
    transformM.setZero();
    // 2D rotation matrix
    transformM << math::cos(rad), -math::sin(rad),
                  math::sin(rad),  math::cos(rad);

    return transformM*vects;
}

void writeControlPoints(const gsMatrix<>& ctrlPts,
                        const std::string& filename)
{
    const gsMatrix<>& points = ctrlPts.transpose();
    gsWriteParaviewPoints(points, filename);
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


} // namespace gismo

