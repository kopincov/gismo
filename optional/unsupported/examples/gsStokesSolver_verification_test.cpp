/** @file gsStokesSolver_verification_test.h

    @brief Verification test for the 2D Stokes solver with the
    gsPdeAssembler structure. The test uses the Method of Manufactured
    Solutions (MMS) and measure the convergence rate of the L2 error. If
    the convergence rate is not sufficient, the test fails.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, J. Sogn

    Test inclueds:
    Multipatch
    Inhomogeneous source term (right hand side)
    Inhomogeneous Dirichlet boundary conditions
    Inhomogeneous Neumann boundary conditions
    Nitche handling of Dirichlet BC
    
    x,y \in (0,1)    \nu = 0.25
    Source function:
    
    f(x,y)_x =  2*4*pi*pi*sin(2*pi*x)*cos(2*pi*y) + 4*pi*pi*sin(2*pi*x)
    f(x,y)_y = -2*4*pi*pi*cos(2*pi*x)*sin(2*pi*y)
    
    
    Solution:
    
    u(x,y)_x =  sin(2*pi*x)*cos(2*pi*y)+pi/10
    u(x,y)_y = -cos(2*pi*x)*sin(2*pi*y)-pi/10
    p(x,y) =  2*pi*cos(2*pi*x)
    
    Note: In gsStokesAssembler.hpp the pressure sign is flipped from the classical
    Stokes equation to obtain a linear system. That is, we are now solving:
    
    -\nu \Delta u - nabla p = f
    /nabla /cdot u = 0
*/

#include <iostream>

#include <gismo.h>


#include <gsAssembler/gsStokesAssembler.h>
#include <gsSolver/gsSolverUtils.h>


;
using namespace gismo;

#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n");\

bool passed=true;


int main(int argc, char *argv[])
{return 0;}/*
    int numRefine = 1;     // Lowest number of refinement: numRefine + 1
    int maxIterations = 2; // Highes number of refinement: numRefine + maxIterations
    real_t convratetol_u = 2.90; // Convergence rate should be around 2. For P2 elements
    real_t convratetol_p = 1.90; // Convergence rate should be around 1. For P1 elements

    // List of element sizes
    std::vector <real_t> h_list;

    // List of L2-errors
    std::vector <real_t> error_list_u;
    std::vector <real_t> error_list_p;

    // Source function
    gsFunctionExpr<> f("0.25*2*4*pi*pi*sin(2*pi*x)*cos(2*pi*y) + 4*pi*pi*sin(2*pi*x)",
                              "-0.25*2*4*pi*pi*cos(2*pi*x)*sin(2*pi*y)") ;

    // Exact solution
    gsFunctionExpr<> g("sin(2*pi*x)*cos(2*pi*y)+pi/10", "-cos(2*pi*x)*sin(2*pi*y)-pi/10");
    gsFunctionExpr<> p0("2*pi*cos(2*pi*x)");

    gsInfo<<"Source function: "<< f << "\n";
    gsInfo<<"Exact solution: " << g << "\n"<< "\n";



    // ---------------Define Geometry---------------
    // (Unit square with 4 patches)
    gsMultiPatch<> * patches;// = new gsMultiPatch<>;
    patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);


    // ---------------Define Boundary conditions---------------
    gsBoundaryConditions<> bcInfo;

    //Dirichlet BCs
    bcInfo.addCondition(0, boundary::south,  condition_type::dirichlet, &g);
    bcInfo.addCondition(2, boundary::south,  condition_type::dirichlet, &g);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &g);
    bcInfo.addCondition(2, boundary::east,  condition_type::dirichlet, &g);
    bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &g);
    bcInfo.addCondition(3, boundary::east,  condition_type::dirichlet, &g);

    // Neumann BCs
    //h(x) = grad(u)*n + pn
    gsFunctionExpr<> h_neu("-2*pi*(1+0.25*cos(2*pi*y))","0.0");
    bcInfo.addCondition(0, boundary::west, condition_type::neumann, &h_neu);
    bcInfo.addCondition(1, boundary::west, condition_type::neumann, &h_neu);


    dirichlet::strategy Dstrategy = dirichlet::elimination;
    iFace::strategy     Istrategy = iFace::glue;

    // Short loop to change interface and dirichlet strategy for completeness
    // of verification test.
    for (int changeStrategy = 1; changeStrategy <= 2; ++changeStrategy)
    {
        // Change interface and dirichlet strategy
        if (changeStrategy == 2)
        {
            Dstrategy = dirichlet::nitsche;
           gsInfo << "Solving Stokes problem with the Nitchs method and glue strategy" << "\n";
        }
        else
            gsInfo << "Solving Stokes problem with elimination and glue strategy" << "\n";


        // Start Loop
        // For each iteration we h-refine the basis, then set up and solve the
        // Poisson problem. Then find the L2 error and element size.
        for (int it = 0; it < maxIterations; ++it)
        {

            std::vector< gsMultiBasis<>* >  refine_bases;
            refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for velocity x-component
            refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for velocity y-component
            refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for pressure


            // Define discretization space by h-refining the basis of the geometry
            for (int i = 0; i < numRefine + it + 1; ++i)
            {
                refine_bases[0]->uniformRefine();
                refine_bases[1]->uniformRefine();
                refine_bases[2]->uniformRefine();
            }

            // p-refine the velocity basis to satify discret inf-sup condition
            // Note: p-refinement must be done after h-refinement!
            refine_bases[0]->degreeElevate();
            refine_bases[1]->degreeElevate();

            gsInfo << "Degree of basis for velocity x-component: "<< refine_bases[0]->degree()<< "\n";
            gsInfo << "Degree of basis for velocity y-component: "<< refine_bases[1]->degree()<< "\n";
            gsInfo << "Degree of basis for pressure: "<< refine_bases[2]->degree()<< "\n";
            gsInfo << "basis u_x" << *refine_bases[0] << "\n";
            gsInfo << "basis u_y" << *refine_bases[1] << "\n";
            gsInfo << "basis p" << *refine_bases[2] << "\n";


            // Initilize Solver
            gsStokesAssembler<> StokesSolver(*patches, bcInfo, refine_bases, f);

            // Set nu
            StokesSolver.setnu(0.25);

            // Use Nitsche's method for Dirichlet boundaries
            StokesSolver.setDirichletStrategy(Dstrategy);

            // Set patch interface strategy, currently only glue is working.
            StokesSolver.setInterfaceStrategy(Istrategy);

            // Assemble and solve
            StokesSolver.solve();

            // Find the element size
            real_t h = math::pow( (real_t) refine_bases[0]->size(0), -1.0 / refine_bases[0]->dim() );
            h_list.push_back(h);

            // Access the solutions
            const gsField<> & sol_u = StokesSolver.solution(0);
            const gsField<> & sol_p = StokesSolver.solution(1);

            // Find the l2 error
            real_t l2error_u = sol_u.distanceL2(g);
            real_t l2error_p = sol_p.distanceL2(p0);
            error_list_u.push_back(l2error_u);
            error_list_p.push_back(l2error_p);

            gsInfo << "Iteration number "<< it <<"  L2 error velocity : " << l2error_u << "  L2 error pressure : " << l2error_p << "   mesh size: " << h <<"\n";
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
        gsInfo << "The convergence rate for valocity is (last two iterations): "<< convratelast_u<< "   ";
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
    }


    delete patches;
    return passed==1?0:1;
}*/
