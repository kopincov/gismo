/** @file stokesMultiPatch.cpp

    @brief Example on how the use one of the stokes assemblers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <iostream>

#include <gismo.h>
#include <gsAssembler/gsStokesAssembler.h>
#include <gsAssembler/gsPdeAssembler.h>
#include <gsSolver/gsSolverUtils.h>


;
using namespace gismo;

int main(int argc, char *argv[])
{ //return 0;}/*

    int numRefine = 1;
    if (argc == 2)
        numRefine = atoi(argv[1]);//3

    bool plot = false; // If set to true, paraview files is generated and stored

    // Source function
    gsFunctionExpr<> f("2*4*pi*pi*sin(2*pi*x)*cos(2*pi*y) + 4*pi*pi*sin(2*pi*x)", "-2*4*pi*pi*cos(2*pi*x)*sin(2*pi*y)",2) ;

    // Exact velocity solution
    gsFunctionExpr<> g("sin(2*pi*x)*cos(2*pi*y)", "-cos(2*pi*x)*sin(2*pi*y)",2);

    // Exact pressure solution
    gsFunctionExpr<> p0("2*pi*cos(2*pi*x)",2);

    gsInfo<<"Source function: "<< f << "\n";
    gsInfo<<"Exact solution: " << g << "\n";


    // Define Geometry
    gsMultiPatch<> patches;
    // Define Geometry (Unit square with 4 patches)
    patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);

    std::vector< gsMultiBasis<>* >  refine_bases;
    refine_bases.push_back(new gsMultiBasis<>( patches ));//Basis for the first component of the velocity
    refine_bases.push_back(new gsMultiBasis<>( patches ));//Basis for the second component of the velocity
    refine_bases.push_back(new gsMultiBasis<>( patches ));//Basis for pressure

    // Define discretization space by refining the basis of the geometry

    for (int i = 0; i < numRefine; ++i)
    {
        refine_bases[0]->uniformRefine();
        refine_bases[1]->uniformRefine();
        refine_bases[2]->uniformRefine();
    }

    // Degree elavateing to get Taylor-Hood elements
    refine_bases[0]->degreeElevate();
    refine_bases[1]->degreeElevate();

    gsInfo << "Degree of basis for velocity x-component: "<< refine_bases[0]->degree(0,0)<< "\n";
    gsInfo << "Degree of basis for velocity y-component: "<< refine_bases[1]->degree(0,0)<< "\n";
    gsInfo << "Degree of basis for pressure: "<< refine_bases[2]->degree(0,0)<< "\n";

    gsInfo << "The number of patches is " << patches.nPatches() << "\n";

    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    //Boundary conditions for each patch
    gsFunctionExpr<> g_dir("sin(2*pi*x)*cos(2*pi*y)", "-cos(2*pi*x)*sin(2*pi*y)",2);
    gsFunctionExpr<> h_neu("-2*pi*(1+cos(2*pi*y))","0.0",2);

    bcInfo.addCondition(0, boundary::south,  condition_type::dirichlet, &g_dir);
    bcInfo.addCondition(2, boundary::south,  condition_type::dirichlet, &g_dir);
    //bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g_dir);
    //bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, &g_dir);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &g_dir);
    bcInfo.addCondition(2, boundary::east,  condition_type::dirichlet, &g_dir);
    bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &g_dir);
    bcInfo.addCondition(3, boundary::east,  condition_type::dirichlet, &g_dir);

    bcInfo.addCondition(0, boundary::west, condition_type::neumann, &h_neu);
    bcInfo.addCondition(1, boundary::west, condition_type::neumann, &h_neu);

    // Initilize assembler
    gsStokesAssembler<> StokesSolver(patches, bcInfo, refine_bases, f);

    // Assmbler options
    StokesSolver.setDirichletStrategy(dirichlet::nitsche);
    StokesSolver.setInterfaceStrategy(iFace::glue);
    StokesSolver.setGeometryTransformation(INVERSE_COMPOSITION);

    StokesSolver.initialize();
    StokesSolver.assemble();
    StokesSolver.solveSystem();
    StokesSolver.reconstructSolution();

    // Access the solutions
    const gsField<> & sol_u = StokesSolver.solution(0);
    const gsField<> & sol_p = StokesSolver.solution(1);

    // Write approximate and exact solution to paraview files
    if (plot)
    {
        gsInfo<<"Writing to Paraview...\n";
        gsWriteParaview<>( sol_u, "Stokes2dvelocity", 10000);
        gsWriteParaview<>( sol_p, "Stokes2dpressure", 10000);
        //gsField<> exact( StokesSolver.patches(), g, false );
        gsField<> exact( patches , g_dir , false ) ;
        gsWriteParaview<>( exact, "Stokes2d_exact", 10000);

    }


    //Clean up!
    freeAll(refine_bases);
    gsInfo << "Test is done: Exiting" << "\n";
    return  0;

}//*/
