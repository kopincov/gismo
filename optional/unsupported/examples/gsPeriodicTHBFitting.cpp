/** @file gsPeriodicTHBFitting

    @brief Demonstrates periodic THB-spline fitting via gsMappedBasis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
n
    Author(s): J. Speh
*/

#include <iostream>
#include <string>

#include <gismo.h>
#include <gismo_dev.h>
#include <gsIO/gsIOUtils.h>

#include "gsMSplines/gsMappedBasis.h"
#include "gsSmoothPatches/gsCompositeBSplineBasis.h"
#include "gsSmoothPatches/gsCompositeHBasis.h"

using namespace gismo;

void computeErrors(std::vector<real_t>& errors,
		   const gsMatrix<>& uv,
		   const gsMatrix<>& xyz,
		   const gsGeometry<>& geom)
{
    errors.clear();
    gsMatrix<> eval;
    geom.eval_into(uv, eval);

    gsVector<real_t, 3> dist;
    for (index_t col = 0; col != eval.cols(); col++)
    {
	dist = xyz.col(col) - eval.col(col);

	errors.push_back(dist.norm());
    }
}


int main(int argc, char* argv[])
{
    std::string input("fitting/discontinuous_cylinder.xml");
    std::string output("out");
    index_t iterations = 5;
    index_t degree = 1;


    gsCmdLine cmd("BSpline Periodic fitting");
    cmd.addString("i", "input", "Input file", input);
    cmd.addString("o", "output", "Output file", output);
    cmd.addInt("I", "iter", "Number of iterations", iterations);
    cmd.addInt("d", "degree", "Degree of surface", degree);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << " \n\nInput arguments: \n\n"
	      << "input: " << input << "\n\n"
	      << "output: " << output << "\n\n"
	      << "iterations: " << iterations << "\n\n"
	      << "degree: " << degree << "\n\n"
	      << "--------------------------------------------------\n" << "\n";

    gsFileData<> fd(input);

    if (!fd.has< gsMatrix<> >())
    {
        gsWarn << "Invalid input.\n";
        return -1;
    }

    gsMatrix<> uv, xyz;
    fd.getId<gsMatrix<> >(0, uv );
    fd.getId<gsMatrix<> >(1, xyz);

    gsKnotVector<> kv1(0, 1, 2, degree + 1, degree);
    gsKnotVector<> kv2(0, 1, 2, degree + 1, degree);

    gsTensorBSplineBasis<2> tensorBasis(kv1, kv2);
    gsTHBSplineBasis<2>  thbBasis(tensorBasis);

    std::vector<unsigned> extension;
    extension.push_back(1);
    extension.push_back(1);

    gsBoxTopology topology(2, 1);
    topology.addInterface(0, boundary::left, 0, boundary::right);
    std::vector< gsHTensorBasis<2>* > vectorOfBasis;
    vectorOfBasis.push_back(&thbBasis);

    gsCompositeHBasis<2, real_t> periodicBasis(vectorOfBasis, topology);

    for (int it = 0; it != iterations; it++)
    {
	gsInfo << "Iteration: " << it + 1 << " / " << iterations << "\n"
		  << "DOF: " << periodicBasis.size() << "\n";

	gsTHBSplineBasis<2> basis = dynamic_cast< const gsTHBSplineBasis<2>& >
	    (periodicBasis.getBase(0));

    gsSparseMatrix<> global2local = periodicBasis.getMapper().asMatrix();

	const int numBasis = basis.size();

	gsSparseMatrix<> A_local(numBasis, numBasis);
	gsMatrix<> B_local(numBasis, 3);

	A_local.setZero();
	B_local.setZero();

	gsHFitting<2, real_t> fitting(uv, xyz, basis, 1, extension, 1e-9);
	fitting.assembleSystem(A_local, B_local);
	fitting.applySmoothing(1e-9, A_local);

	gsSparseMatrix<> A = global2local.transpose() * A_local * global2local;
	gsMatrix<> B = global2local.transpose() * B_local;

	A.makeCompressed();
	gsSparseSolver<>::BiCGSTABILUT solver(A);
	if ( solver.preconditioner().info() != Eigen::Success )
	{
	    gsWarn<<  "The preconditioner failed. Aborting.\n";
	}

	gsMatrix<> coefs = global2local * solver.solve(B);


	gsTHBSpline<2> surface(basis, coefs);


	std::vector<real_t> errors;
	computeErrors(errors, uv, xyz, surface);
	const real_t maxError = *std::max_element(errors.begin(), errors.end());
	gsInfo << "Maximum error: " << maxError << "\n";

	std::string out = output + "_geom" + util::to_string(it);
	gsInfo << "  Writing geometry: " << out << " ..." << "\n";
	gsWriteParaview(surface, out, 10000);

	gsMesh<> mesh;
	makeMesh<>(surface.basis(), mesh);
	surface.evaluateMesh(mesh);
	out = output + "_mesh" + util::to_string(it);
	gsInfo << "  Writing mesh: " << out << " ..." << "\n";
	gsWriteParaview(mesh, out);

	if (maxError < 0.01)
	{
	    gsInfo << "Maximum error is below trashold." << "\n";
	    return 0;
	}


	if (it != iterations - 1)
	{
	    std::vector<index_t> boxes = fitting.getBoxes(errors, 0.01);
	    periodicBasis.refineElements(0, boxes);

	    gsInfo << "Inserted " << boxes.size() / 5 << " boxes." << "\n";
	}

    }


    return 0;

}

