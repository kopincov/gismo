/** @file gsParameterLines.h

    @brief Outputs parameter lines to paraview.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <fstream>
#include <sstream>

#include <gsIO/gsIOUtils.h>
#include <gsIO/gsWriteParaview.h>

#include <gsCore/gsMultiPatch.h>
#include <gsUtils/gsPointGrid.h>


#pragma once

using namespace gismo;

void addParameterLinesToMesh(gsMesh<>& mesh,
			     const gsGeometry<>& geom,
			     const int resolution,
			     int dir,
			     const int num_parameter_lines)
{
    const gsMatrix<>& support  = geom.support();
    const gsVector<>& supStart = support.col(0);
    const gsVector<>& supEnd = support.col(1);

    int dir2 = (dir == 0) ? 1 : 0;
    const gsMatrix<>& samples = gsPointGrid(supStart(dir2), supEnd(dir2), resolution);


    gsMatrix<real_t> ones(1, samples.cols());
    ones.setOnes();

    gsMatrix<> xParams(1, num_parameter_lines);

    for (int i = 0; i != num_parameter_lines; ++i)
    {
	const real_t t = i * 1.0 / (num_parameter_lines - 1);
	xParams(0, i) = supStart(dir) * (1 - t) + supEnd(dir) * t;
    }

    typedef gsMesh<real_t>::VertexHandle Vertex;
    for (index_t index = 0; index != xParams.size(); index++)
    {
        const real_t constant = xParams(0, index);
        gsMatrix<> params(2, samples.cols());

        params.row(dir) = constant * ones;
        params.row(dir2) = samples;

        gsMatrix<> result;
        geom.eval_into(params, result);

        gsVector<> vertex(result.rows());
        vertex = result.col(0);
        mesh.addVertex(vertex);
        for (int col = 1; col != result.cols(); ++col)
        {
            Vertex first = mesh.vertices().back();
            vertex = result.col(col);
            mesh.addVertex(vertex);
            Vertex second = mesh.vertices().back();
            mesh.addEdge(first, second);
        }
    }
}

void addParameterLinesToMesh3D(gsMesh<>& mesh,
			       const gsGeometry<>& geom,
			       const int resolution,
			       int dir,
			       const int num_parameter_lines)
{
    const gsMatrix<>& support  = geom.support();
    const gsVector<>& supStart = support.col(0);
    const gsVector<>& supEnd = support.col(1);

    const gsMatrix<>& samples = gsPointGrid(supStart(dir), supEnd(dir), resolution);

    int dir1 = (dir == 0) ? 1 : 0;
    int dir2 = (dir == 2) ? 1 : 2;

    gsVector<> a(2);
    a << supStart(dir1), supStart(dir2);

    gsVector<> b(2);
    b << supEnd(dir1), supEnd(dir2);

    gsVector<unsigned> np(2);
    np << static_cast<unsigned>(num_parameter_lines), static_cast<unsigned>(num_parameter_lines);
    gsMatrix<> params = gsPointGrid(a, b, np);

    gsMatrix<real_t> ones(1, samples.cols());
    ones.setOnes();

    typedef gsMesh<real_t>::VertexHandle Vertex;
    for (index_t index = 0; index != params.cols(); ++index)
    {
        const real_t constant1 = params(0, index);
        const real_t constant2 = params(1, index);
        gsMatrix<> params(3, samples.cols());

        params.row(dir1) = constant1 * ones;
	params.row(dir2) = constant2 * ones;
        params.row(dir) = samples;

        gsMatrix<> result;
        geom.eval_into(params, result);

        gsVector<> vertex(result.rows());
        vertex = result.col(0);
        mesh.addVertex(vertex);
        for (int col = 1; col != result.cols(); ++col)
        {
            Vertex first = mesh.vertices().back();
            vertex = result.col(col);
            mesh.addVertex(vertex);
            Vertex second = mesh.vertices().back();
            mesh.addEdge(first, second);
        }
    }
}

void writeParameterLines_(gsMesh<>& mesh,
			  const gsGeometry<>& geom,
			  const int resolution,
			  const int num_parameter_lines)
{
    if (geom.parDim() == 2)
    {
	addParameterLinesToMesh(mesh, geom, resolution, 0, num_parameter_lines);
	addParameterLinesToMesh(mesh, geom, resolution, 1, num_parameter_lines);
    }
    else if (geom.parDim() == 3)
    {
	addParameterLinesToMesh3D(mesh, geom, resolution, 0, num_parameter_lines);
	addParameterLinesToMesh3D(mesh, geom, resolution, 1, num_parameter_lines);
	addParameterLinesToMesh3D(mesh, geom, resolution, 2, num_parameter_lines);
    }
}

void writeParameterLines(const gsGeometry<>& geom,
			 const std::string& output,
			 const int resolution,
			 const int num_parameter_lines)
{

    gsMesh<> mesh;
    writeParameterLines_(mesh, geom, resolution, num_parameter_lines);
    gsWriteParaview(mesh, output);
}


void writeParameterLines(const gsMultiPatch<>& multipatch,
			 const std::string& output,
			 const int resolution,
			 const int num_parameter_lines)
{
    gsMesh<> mesh;
    for (size_t i = 0; i != multipatch.nPatches(); i++)
    {
        gsGeometry<>& geom = multipatch.patch(i);
	writeParameterLines_(mesh, geom, resolution, num_parameter_lines);
    }

    gsWriteParaview(mesh, output);
}
