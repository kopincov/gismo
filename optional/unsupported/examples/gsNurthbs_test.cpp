/** @file gsNurthbs_test.cpp

    @brief A basic test for checking the gsNurthbs geometry

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

#include <gismo.h>

#include <gsHSplines/gsNurthbs.h>
#include <gsHSplines/gsNurthbsBasis.h>

using namespace gismo;


/*
This is just a very rough check of a NUR-THB-Spline geometry.

A NURTHBS-geometry is created and sampling points are mapped to
the physical domain.

The same geometry map is "reconstructed" using tensor-product
B-splines, which are multiplied with the corresponding
weights and geometry coefficients.

The NURTHBS-geometry is also refined locally.

The differences between these three versions of the same geometry map
are then computed and displayed.
*/

int main (int argc, char** argv)
{
    int result = 0;

    bool plot = false;
    index_t grid_n(10);
    // Command line arguments
    gsCmdLine cmd("Example testing a Nur-THB-spline geometry");
    gsInfo << "Type \"-h\" to see the description for all arguments.\n\n";
    cmd.addInt   ("n", "gridpoints", "Number of grid-points per coordinate direction ", grid_n);
    cmd.addSwitch("p", "plot", "Plot (otherwise no plot)", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // compute sampling points
    gsVector<real_t> vec_a(2);
    vec_a << 0, 0;
    gsVector<real_t> vec_b(2);
    vec_b << 1, 1;
    gsVector<unsigned> vec_n(2);
    vec_n << unsigned(grid_n), unsigned(grid_n);
    gsMatrix<real_t> u = gsPointGrid( vec_a, vec_b, vec_n );

    // create knot vector
    // (start, end, #interior, multip end, mult interior, degree)
    gsKnotVector<real_t> kv(0,1,0,3);

    // create 2D Tensor-B-Spline basis
    gsTensorBSplineBasis<2,real_t> tens( kv, kv);

    // create THB-Spline basis
    gsTHBSplineBasis<2,real_t> thb( tens );

    // coefficients of the geometry mapping
    // quarter annulus, inner radius = 1, outer radius = 3
    gsMatrix<real_t> coeffs(9,2);
    coeffs << 1, 0,
            2, 0,
            3, 0,
            1, 1,
            2, 2,
            3, 3,
            0, 1,
            0, 2,
            0, 3;

    // weights for the rational basis
    // for exact representation of the quarter annulus
    gsMatrix<real_t> weights( thb.size(), 1);
    weights.setOnes();
    weights(3,0) = (real_t)(1.0)/(real_t)(std::sqrt(2.0));
    weights(4,0) = (real_t)(1.0)/(real_t)(std::sqrt(2.0));
    weights(5,0) = (real_t)(1.0)/(real_t)(std::sqrt(2.0));


    // set up rational basis from THB-spline basis and weights
    gsNurthbsBasis<2,real_t> nurthb(thb.clone().release(), weights);

    // get NurThbs-Geometry from NurThbsBasis and coefficients
    gsNurthbs<2,real_t> nurthbGeo(nurthb, coeffs );

    // ---------- NURTHBS
    // map the sampling points to the physical domain with the
    // NurTHBs-mapping
    gsMatrix<real_t> U_nurthbs;
    nurthbGeo.eval_into( u, U_nurthbs );

    // ---------- NURTHBS, Refined
    // adding some random local refinement
    std::vector<index_t> boxes;
    boxes.push_back(5);
    boxes.push_back(4);
    boxes.push_back(0);
    boxes.push_back(12);
    boxes.push_back(6);

    boxes.push_back(4);
    boxes.push_back(8);
    boxes.push_back(8);
    boxes.push_back(16);
    boxes.push_back(12);

    boxes.push_back(3);
    boxes.push_back(4);
    boxes.push_back(2);
    boxes.push_back(6);
    boxes.push_back(8);
    nurthbGeo.refineElements( boxes );

    gsMatrix<real_t> U_nurthbsRef;
    nurthbGeo.eval_into( u, U_nurthbsRef );


    // ---------- own mapping with tensor-product splines

    //gsMatrix<real_t> U_tens( u->rows(), u->cols() );
    gsMatrix<real_t> U_tens_wWeights( u.rows(), u.cols() );

    // matrix for denominator
    gsMatrix<real_t> denom( 1, u.cols() );
    //U_tens.setZero();
    U_tens_wWeights.setZero();
    denom.setZero();

    gsMatrix<index_t> actTens;
    gsMatrix<real_t> bTens;

    // get indices and basis function values of
    // tensor-product basis functions
    tens.active_into( u, actTens );
    tens.eval_into(   u, bTens );

    for( size_t j = 0; j < size_t( actTens.cols() ); j++)
        for( size_t i = 0; i < size_t( actTens.rows() ); i++)
        {
            unsigned a = actTens(i,j);
            U_tens_wWeights(0,j) += weights(a,0) * coeffs( a,0 ) * bTens(i,j);
            U_tens_wWeights(1,j) += weights(a,0) * coeffs( a,1 ) * bTens(i,j);
            denom(0,j) += weights(a,0) * bTens(i,j);
        }

    for( size_t j = 0; j < size_t( U_tens_wWeights.cols() ); j++ )
    {
        U_tens_wWeights(0,j) /= denom(0,j);
        U_tens_wWeights(1,j) /= denom(0,j);
    }

    real_t diffSq1(0.0);
    real_t diffSq2(0.0);

    for( size_t j = 0; j < size_t( u.cols() ); j++)
    {
        diffSq1 += (U_tens_wWeights(0,j)-U_nurthbs(0,j))*(U_tens_wWeights(0,j)-U_nurthbs(0,j));
        diffSq1 += (U_tens_wWeights(1,j)-U_nurthbs(1,j))*(U_tens_wWeights(1,j)-U_nurthbs(1,j));

        diffSq2 += (U_nurthbs(0,j)-U_nurthbsRef(0,j))*(U_nurthbs(0,j)-U_nurthbsRef(0,j));
        diffSq2 += (U_nurthbs(1,j)-U_nurthbsRef(1,j))*(U_nurthbs(1,j)-U_nurthbsRef(1,j));
    }

    std::cout << "number of sampling points: " << u.cols() << std::endl;
    std::cout << "difference NURTHBS-Geometry and 'manual' computation       : " << diffSq1 << std::endl;
    std::cout << "difference non-refined and locally refined NURTHBS-Geometry: " << diffSq2 << std::endl;
    std::cout << std::endl;

    if( plot )
    {
        gsMultiPatch<real_t> patches;
        patches.addPatch( nurthbGeo );

        gsWriteParaview<>( patches, "Nurthbs_physMesh", 1000, true);
        result = system("paraview Nurthbs_physMesh.pvd &");
    }

    return result;
}
