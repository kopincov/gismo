/** @file gsCompatibleStokesSolver.cpp

    @brief Test the convergence rate for Stokes using divergence perserving elements

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/



// Author:Jarle Sogn
//
// Verification test for the 2D Stokes solver with the gsStokesCompatibleSolver
// structure. The test uses the Method of Manufactured Solutions (MMS)
// and measure the convergence rate of the L2 error. If the convergence
// rate is not sufficient, the test fails.
// The test used the divergence perserving transformation with generalized
// Ravier-Thomas elements.
//
// Test inclueds:
// Inhomogeneous source term (right hand side)
// Inhomogeneous Dirichlet boundary conditions (Nitsch)
// Inhomogeneous Neumann boundary conditions
//
//     x,y \in (0,1)    \nu = 0.25
// Source function:
//
//      f(x,y)_x =  0.5*pi^3*sin(pi*y)*(2*(sin(pi*x*0.5)^2)- cos(pi*x)) + pi*cos(pi*x)
//      f(x,y)_y = -0.5*pi^3*sin(pi*x)*(2*(sin(pi*y*0.5)^2)- cos(pi*y))"

//
// Solution:
//
//     u(x,y)_x =  pi*sin(pi*y)*(sin(pi*x*0.5)^2)
//     u(x,y)_y = -pi*sin(pi*x)*(sin(pi*y*0.5)^2)
//       p(x,y) = -sin(pi*x)
//
// Note: In gsStokesAssembler.hpp the pressure sign is flipped from the classical
// Stokes equation to obtain a linear system. That is, we are now solving:
//
// -\nu \Delta u - nabla p = f
//          /nabla /cdot u = 0

#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsAssembler/gsStokesAssembler.h>
#include <gsSolver/gsSolverUtils.h>


;
using namespace gismo;

#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n");\

bool passed=true;

int main()
{
    int numRefine = 1;     // Lowest number of refinement: numRefine + 1
    int maxIterations = 3; // Highes number of refinement: numRefine + maxIterations
    real_t convratetol_u = 2.90; // Convergence rate should be around 2. For P2 elements
    real_t convratetol_p = 1.90; // Convergence rate should be around 1. For P1 elements

    // List of element sizes
    std::vector <real_t> h_list;

    // List of L2-errors
    std::vector <real_t> error_list_u;
    std::vector <real_t> error_list_p;

    // Source function

    //Mardal's MMS
    gsFunctionExpr<> f(" 0.25*2*pi^3*sin(2*pi*y*0.5)*(2*(sin(pi*x*0.5)^2)- cos(2*pi*x*0.5)) + 0.5*2*pi*cos(2*pi*x*0.5)",
                        "-0.25*2*pi^3*sin(2*pi*x*0.5)*(2*(sin(pi*y*0.5)^2)- cos(2*pi*y*0.5))",2) ;
    gsFunctionExpr<> g("pi*sin(2*pi*y*0.5)*(sin(pi*x*0.5)^2)", "-pi*sin(2*pi*x*0.5)*(sin(pi*y*0.5)^2)",2);

    gsFunctionExpr<> h_neu("sin(pi*x)-0.5*pi*pi*sin(pi*x)*sin(pi*y)",
                             "pi*pi*cos(pi*x)*sin(pi*0.5*y)*sin(pi*0.5*y)",2);

    //Sogn's MMS (Homogenius BC on Annulus) (NON-DIVERGENCE FREE)
    //gsFunctionExpr<> g("((0.5*(y-x+1))^2+x-1)*(0.5*(0.5*(y-x)+1)^2+x-2)x*y","((0.5*(y-x+1))^2+x-1)*(0.5*(0.5*(y-x)+1)^2+x-2)x*y",2);
    //gsFunctionExpr<> f("-0.125*(-2*x^4+x^3*(14*y-3)+x*x*(-24y*y+9*y+23)+x*(14*y^3+9*y*y-21*y-18)-y*(2*y^3+3*y*y-23*y+18))",
    //                          " -0.125*(-2*x^4+x^3*(14*y-3)+x*x*(-24y*y+9*y+23)+x*(14*y^3+9*y*y-21*y-18)-y*(2*y^3+3*y*y-23*y+18))",2);


    //Standard
    //gsFunctionExpr<> f(" 2*4*pi*pi*sin(2*pi*x)*cos(2*pi*y)",
    //                     "-2*4*pi*pi*cos(2*pi*x)*sin(2*pi*y)",2) ;
    //gsFunctionExpr<> g("sin(2*pi*x)*cos(2*pi*y)", "-cos(2*pi*x)*sin(2*pi*y)",2);


    //Quadratic flow
    //gsFunctionExpr<> g("y*0.5*(1-y*0.5)","0.0",2);
    //gsFunctionExpr<> f("0.0", "0.5",2) ;
    //gsFunctionExpr<> g("0.0","x*0.5*(1-x*0.5)",2);

    gsFunctionExpr<> p0("-sin(pi*x)",2);

    gsInfo<<"Source function: "<< f << "\n";
    gsInfo<<"Exact solution: "<< g <<".\n" << "\n";


    // Define Geometry
    gsMultiPatch<> * patches;
    patches = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    //patches = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineSquare(2.0,0,0) );
    //patches = gsReadFile<>( "../planar/linear_map1.xml");

    // Define Geometry (Unit square with 4 patches)
    //patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);


    // ---------------Define Boundary conditions---------------
    gsBoundaryConditions<> bcInfo;

    //Dirichlet BCs
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &g);
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &g);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &g);
    bcInfo.addCondition( boundary::north, condition_type::neumann, &h_neu);

    // Start Loop
    // For each iteration we h-refine the basis, then set up and solve the
    // Poisson problem. Then find the L2 error and element size.
    for (int it = 0; it < maxIterations; ++it)
    {
        std::vector< gsMultiBasis<>* >  refine_bases;
        refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for velocity_x
        refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for velocity_y
        refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for pressure
        gsMultiBasis<> * basisNorm =  new gsMultiBasis<>( *patches );//Basis used for integrating when finding the L2 Norm

        //The Annulus geometry has different order on the x and y directions
        //p-refine to get equal polynomial degree s,t directions
        refine_bases[0]->degreeElevate(1,0);
        refine_bases[1]->degreeElevate(1,0);
        refine_bases[2]->degreeElevate(1,0);

        // Define discretization space by h-refining the basis of the geometry
        for (int i = 0; i < numRefine + it + 1; ++i)
        {
            refine_bases[0]->uniformRefine();
            refine_bases[1]->uniformRefine();
            refine_bases[2]->uniformRefine();
            basisNorm->uniformRefine();
        }

        //Taylor-Hood element
        //refine_bases[0]->degreeElevate();
        //refine_bases[1]->degreeElevate();
        //refine_bases[2]->degreeElevate();
        //basisNorm->degreeElevate();

        //Raviart-Thomas
        refine_bases[0]->degreeElevate(1,0);
        refine_bases[1]->degreeElevate(1,1);
        basisNorm->degreeElevate();

        gsInfo << "Degrees of basis for velocity x-component: "<< refine_bases[0]->degree(0,0)<< " and "<< refine_bases[0]->degree(0,1)<<"\n";
        gsInfo << "Degrees of basis for velocity y-component: "<< refine_bases[1]->degree(0,0)<< " and "<< refine_bases[1]->degree(0,1)<<"\n";
        gsInfo << "Degrees of basis for pressure: "<< refine_bases[2]->degree(0,0) << " and "<< refine_bases[2]->degree(0,1)<< "\n";
        gsInfo << "Degrees of basis for intergration: "<< basisNorm->degree(0,0)<< " and "<< basisNorm->degree(0,1)<<"\n";

        // Initilize Solver
        gsStokesAssembler<> StokesSolver(*patches, bcInfo, refine_bases, f);

        StokesSolver.setDirichletStrategy(dirichlet::nitsche);
        StokesSolver.setInterfaceStrategy(iFace::glue);
        StokesSolver.setGeometryTransformation(DIV_CONFORMING);
        //StokesSolver4.setGeometryTransformation(INVERSE_COMPOSITION);

        // Assemble and solve
        StokesSolver.solve();

        // Access the solutions
        const gsField<> & sol_u = StokesSolver.solution(0);
        const gsField<> & sol_p = StokesSolver.solution(1);

        // Find the element size
        real_t h = math::pow( (real_t) refine_bases[0]->size(0), -1.0 / refine_bases[0]->dim() );
        h_list.push_back(h);

        // Find the l2 error
        real_t l2error_u = sol_u.distanceL2(g, *basisNorm);
        real_t l2error_p = sol_p.distanceL2(p0);
        error_list_u.push_back(l2error_u);
        error_list_p.push_back(l2error_p);

        gsInfo << "Iteration number "<< it <<"  L2 error velocity : " << l2error_u << "  L2 error pressure : " << l2error_p << "   mesh size: " << h <<"\n";
        freeAll( refine_bases );
        delete (basisNorm);
    }

    // Print out the L2 errors and element sizes
    gsInfo << "L2: error_u     error_p   Element sizes" << "\n";
    for (unsigned k= 0; k<error_list_u.size(); ++k)
         gsInfo << "    "<<error_list_u[k] << " "<< error_list_p[k] << "  "<< h_list[k] <<"\n";

    // Finding convergence rate in to differnt ways

    // Convergence rate found by last to errors are elemet sizes
    real_t convratelast_u = math::log(error_list_u[error_list_u.size()-2]/error_list_u[error_list_u.size()-1]) /
            math::log(h_list[h_list.size()-2]/h_list[h_list.size()-1]);
    real_t convratelast_p = math::log(error_list_p[error_list_p.size()-2]/error_list_p[error_list_p.size()-1]) /
            math::log(h_list[h_list.size()-2]/h_list[h_list.size()-1]);

    // Convergence rate found by a least square method
    // PS: This method of finding convergence rate mights require some initialt
    // refinement such that convergence have startet for the largest element
    // size value.
    real_t convrateavg_u =  gsSolverUtils<>::convergenceRateLS(error_list_u,h_list);
    real_t convrateavg_p =  gsSolverUtils<>::convergenceRateLS(error_list_p,h_list);

    gsInfo << "The convergence rate for velocity is (fitted line)        : "<< convrateavg_u<< "   ";
    TEST(convrateavg_u > convratetol_u);
    gsInfo << "The convergence rate for velocity is (last two iterations): "<< convratelast_u<< "   ";
    TEST(convratelast_u > convratetol_u);

    gsInfo << "The convergence rate for pressure is (fitted line)        : "<< convrateavg_p  << "   ";
    TEST(convrateavg_p  > convratetol_p);
    gsInfo << "The convergence rate for pressure is (last two iterations): "<< convratelast_p << "   ";
    TEST(convratelast_p > convratetol_p);

    gsInfo << "\n";

    //Clean up before next test
    error_list_u.clear();
    error_list_p.clear();
    h_list.clear();


    //real_t conditionnumber =  gsSolverUtils<>::conditionNumber(StokesSolver4.systemMatrix());
    //gsInfo << "the condotion number is : " << conditionnumber << "\n";


    delete patches;
    gsInfo << "Test is done: ";
    if (passed)
        gsInfo << "Success!\n";
    else
        gsInfo << "Failure!\n";
    return passed==1?0:1;
}
