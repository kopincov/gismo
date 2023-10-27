/** @file gsFittingUtils.h

    @brief Utility functions for fitting

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#pragma once

#include "gsMultipatchFitting.h"

namespace gismo
{

gsMultiBasis<> makeMultiBasisBSplines(const gsMultiPatch<>& mp,
				      const int degree = 2,
				      const int num_internal_knots = 0)
{
    gsMultiBasis<>::BasisContainer bases;
    
    for (unsigned i = 0; i != mp.nPatches(); i++)
    {
        gsKnotVector<> kv1(0.0, 1.0, num_internal_knots, degree + 1);
        gsKnotVector<> kv2(0.0, 1.0, num_internal_knots, degree + 1);
        gsTensorBSplineBasis<2> tenBasis(kv1, kv2);
        bases.push_back(new gsTensorBSplineBasis<2>(tenBasis));
    }
    std::vector< patchSide > boundaries = mp.boundaries();
    std::vector< boundaryInterface > interfaces = mp.interfaces();
        
    gsBoxTopology topology(2, mp.nPatches(), boundaries, interfaces);

    gsMultiBasis<> m(bases, topology);

    return m;    
}

int number_of_patches(const std::string& filename)
{
    gsFileData<> fd(filename);
    
    int i = 0;
    while (true)
    {
	std::cerr.setstate(std::ios_base::failbit);
	gsMatrix<>* matrix = fd.getId< gsMatrix<> >(i).release();
	std::cerr.clear();
	if (matrix == NULL)
	{
	    break;
	}
	else
	{
	    ++i;
	    delete matrix;
	}
    }
    return i / 2;
}


void getPoints(std::vector< gsMatrix<> >& params,
               std::vector< gsMatrix<> >& points,
               const std::string& filename)
{
    const int num_of_patches = number_of_patches(filename);
    
    gsFileData<> fd(filename);

    for (int i = 0; i != num_of_patches; i++)
    {
        const int id_par = 2 * i;
        const int id_pts = 2 * i + 1;

        gsMatrix<>* par = fd.getId< gsMatrix<> >(id_par).release();
        gsMatrix<>* pts = fd.getId< gsMatrix<> >(id_pts).release();

        params.push_back(*par);
        points.push_back(*pts);
    }
}


gsMatrix<> solve(const gsSparseMatrix<>& A,
                 const gsMatrix<>& B)
{
    gsSparseSolver<>::BiCGSTABILUT solver(A);
    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
        return gsMatrix<>();
    }
    return solver.solve(B);
}


gsMultiPatch<> fitting(const gsMultiBasis<>& basis,
                       const std::vector< gsMatrix<> >& params,
                       const std::vector< gsMatrix<> >& points,
                       const real_t lambda)
{
    gsSparseMatrix<> A;
    gsMatrix<> B;

    fittingMatrix(basis, params, points, lambda, A, B);

    const gsMatrix<> X = solve(A, B);

    return constructSolution(basis, X);
}


}
