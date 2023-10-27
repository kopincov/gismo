/** @file stokes_Piola_2D.cpp

    @brief Test solving the stokes problem with div preserving transformation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <iostream>

#include <gismo.h>

#include <gsAssembler/gsStokesAssembler.h>
//#include <gsAssembler/gsStokesAssemblerAssembler.h>
#include <gsCore/gsDivConSolution.h>
#include <gsAssembler/gsPdeAssembler.h>
#include <gsSolver/gsSolverUtils.h>
#include <gsSolver/gsStokesIterativeSolver.h>


;
using namespace gismo;


int main(int argc, char *argv[])
{

    int numRefine = 2;
    int numDegree = 0;
    if (argc >= 2)
        numRefine = atoi(argv[1]);
    if (argc >= 3)
        numDegree = atoi(argv[2]);

    bool plot = true; // If set to true, paraview file is generated and launched on exit
    bool onlyAssemble = false;

    // Source function
    //gsFunctionExpr<> f("2*pi^2*sin(pi*x)*sin(pi*y)",2) ;

    //Mardal's MMS
    gsFunctionExpr<> f(" 0.25*2*pi^3*sin(2*pi*y*0.5)*(2*(sin(pi*x*0.5)^2)- cos(2*pi*x*0.5)) + 0.5*2*pi*cos(2*pi*x*0.5)",
                              "-0.25*2*pi^3*sin(2*pi*x*0.5)*(2*(sin(pi*y*0.5)^2)- cos(2*pi*y*0.5))",2) ;
    gsFunctionExpr<> g("pi*sin(2*pi*y*0.5)*(sin(pi*x*0.5)^2)", "-pi*sin(2*pi*x*0.5)*(sin(pi*y*0.5)^2)",2);

    //Sogn's MMS (Homogenius BC on Annulus) (NON-DIVERGENCE FREE)
    //gsFunctionExpr<> g("((0.5*(y-x+1))^2+x-1)*(0.5*(0.5*(y-x)+1)^2+x-2)x*y","((0.5*(y-x+1))^2+x-1)*(0.5*(0.5*(y-x)+1)^2+x-2)x*y",2);
    //gsFunctionExpr<> f("-0.125*(-2*x^4+x^3*(14*y-3)+x*x*(-24y*y+9*y+23)+x*(14*y^3+9*y*y-21*y-18)-y*(2*y^3+3*y*y-23*y+18))",
    //                          " -0.125*(-2*x^4+x^3*(14*y-3)+x*x*(-24y*y+9*y+23)+x*(14*y^3+9*y*y-21*y-18)-y*(2*y^3+3*y*y-23*y+18))",2);


    //Standard
    //gsFunctionExpr<> f(" 2*4*pi*pi*sin(2*pi*x)*cos(2*pi*y)",
    //                          "-2*4*pi*pi*cos(2*pi*x)*sin(2*pi*y)",2) ;
    //gsFunctionExpr<> g("sin(2*pi*x)*cos(2*pi*y)", "-cos(2*pi*x)*sin(2*pi*y)",2);

    // Exact solution

    //Quadratic flow
    //gsFunctionExpr<> g("y*0.5*(1-y*0.5)","0.0",2);
    //gsFunctionExpr<> f("0.0", "0.5",2) ;
    //gsFunctionExpr<> g("0.0","x*0.5*(1-x*0.5)",2);

    gsFunctionExpr<> p0("-sin(pi*x)",2);

    //gsInfo<<"Source function "<< f << "\n";
    //gsInfo<<"Exact solution "<< g <<".\n" << "\n";


    // Define Geometry
    gsMultiPatch<> * patches;
    patches = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    //patches = new gsMultiPatch<>(gsNurbsCreator<>::BSplineSquare(2.0,0,0));
    // Unit cube
    //patches = gsReadFile<>( "../planar/linear_map1.xml");

    // Define Geometry (Unit square with 4 patches)
    //patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);

    std::vector< gsMultiBasis<>* >  refine_bases;
    refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for velocity_x
    refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for velocity_y
    refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for pressure
    gsMultiBasis<> * basisNorm =  new gsMultiBasis<>( *patches );//Basis used for integrating when finding the L2 Norm
    // Define discretization space by refining the basis of the geometry

    //p-refine to get equal polynomial degree s,t directions
    refine_bases[0]->degreeElevate(1,0);
    refine_bases[1]->degreeElevate(1,0);
    refine_bases[2]->degreeElevate(1,0);

    gsInfo << "The number of uniform h-refinements is: " <<numRefine <<"\n";
    for (int i = 0; i < numRefine; ++i)
    {
        refine_bases[0]->uniformRefine();
        refine_bases[1]->uniformRefine();
        refine_bases[2]->uniformRefine();
        basisNorm->uniformRefine();
    }
    gsDebug <<  " refine_bases.size() " << refine_bases.size() << "\n";

    //Taylor-Hood element
    //refine_bases[0]->degreeElevate();
    //refine_bases[1]->degreeElevate();
    //basisNorm->degreeElevate();


    for (int i = 0; i < numDegree; ++i)
    {
        refine_bases[0]->degreeElevate();
        refine_bases[1]->degreeElevate();
        refine_bases[2]->degreeElevate();
        basisNorm->degreeElevate();
    }

    //Raviart-Thomas
    refine_bases[0]->degreeElevate(1,0);
    refine_bases[1]->degreeElevate(1,1);
    basisNorm->degreeElevate();

    gsInfo << "Degrees of basis for velocity x-component: "<< refine_bases[0]->degree(0,0)<< " and "<< refine_bases[0]->degree(0,1)<<"\n";
    gsInfo << "Degrees of basis for velocity y-component: "<< refine_bases[1]->degree(0,0)<< " and "<< refine_bases[1]->degree(0,1)<<"\n";
    gsInfo << "Degrees of basis for pressure: "<< refine_bases[2]->degree(0,0) << " and "<< refine_bases[2]->degree(0,1)<< "\n";

    gsInfo << "Number of patches is " << patches->nPatches() << "\n";

    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    //Dirichlet BCs
    gsFunctionExpr<> h_neu("sin(pi*x)-0.5*pi*pi*sin(pi*x)*sin(pi*y)",
                             "pi*pi*cos(pi*x)*sin(pi*0.5*y)*sin(pi*0.5*y)",2);

    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &g);//Annulus: small arch lenght
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &g);//Annulus: Large arch lenght
    //bcInfo.addCondition( boundary::north, condition_type::dirichlet, &g);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &g);
    //bcInfo.addCondition( boundary::front, condition_type::dirichlet, &g);
    //bcInfo.addCondition( boundary::back, condition_type::dirichlet, &g);

    bcInfo.addCondition( boundary::north,  condition_type::neumann, &h_neu);
    //bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &g);


    /////////////////// Setup solver ///////////////////
    //Initilize Solver
    gsStokesAssembler<> StokesAssembler(*patches, bcInfo, refine_bases, f);

    gsDebug << "-------------Testing << operator--------------\n";
    gsDebug << "\n" << "Print Stokes solver 42:\n" << StokesAssembler<< "\n";

    dirichlet::strategy m_dirStrategy = dirichlet::nitsche;
    if (argc == 4)
    {
        if (atoi(argv[3]) == 0)
        {
            m_dirStrategy = dirichlet::nitsche;
            gsInfo << "Using Full nitsche" << "\n";
        }
        if (atoi(argv[3]) == 1)
        {
            m_dirStrategy = dirichlet::eliminatNormal;
            gsInfo << "Using strong normal component" << "\n";
        }
    }
    StokesAssembler.setDirichletStrategy(m_dirStrategy);

    iFace::strategy m_interfaceStrategy = iFace::glue;
    StokesAssembler.setInterfaceStrategy(m_interfaceStrategy);

    ValueTransformationType m_geoTrans = DIV_CONFORMING; //INVERSE_COMPOSITION
    StokesAssembler.setGeometryTransformation(m_geoTrans);

    gsDebug << "-------------Testing initialize()--------------\n";
    StokesAssembler.initialize();

    gsDebug << "-------------Testing assemble()--------------\n";
    gsStopwatch time;
    StokesAssembler.assemble();
    const double assmTime = time.stop();
    gsInfo << "Assembling time: " <<  assmTime << " s" << "\n";

    if (onlyAssemble)
    {
        freeAll( refine_bases );
        delete patches;
        gsInfo << "Test is done (only assemble): Exiting" << "\n";
        return  0;
    }

    gsDebug << "-------------Testing solveSystem()--------------\n";
    gsStokesIterativeSolver<real_t> StokesSolver(&StokesAssembler);
    if (refine_bases[2]->basis(0).size() < 50) //Using Canonical Preconditiors with MINRES
    {
        gsInfo << "Solving with canonical preconditioner..." << "\n";
        StokesSolver.solveWithExactInverse();
        gsInfo << "\nCondition number: \n" << StokesSolver.conditionNumber(true) << "\n";
    }
    else
    {
        gsInfo << "Solving with multigrid and Gauss-Seidel preconditioner..." << "\n";
        StokesSolver.solveWithMGandGS();
    }
    StokesAssembler.setSolution(StokesSolver.systemSolution());

    gsDebug << "-------------Testing recontructionSystem()--------------\n";
    StokesAssembler.reconstructSolution();


    gsDebug << "------------- Access the solutions --------------\n";
    // Access the solutions
    const gsField<> & sol_u = StokesAssembler.solution(0);
    const gsField<> & sol_p = StokesAssembler.solution(1);



    // Find the l2 error
    gsInfo << "Degrees of basis for intergration: "<< basisNorm->degree(0,0)<< " and "<< basisNorm->degree(0,1)<<"\n";
    real_t l2error_u = sol_u.distanceL2(g,*basisNorm);
    real_t l2error_p = sol_p.distanceL2(p0);
    gsInfo << "The L2 Error for velocity error: " << l2error_u << "\n";
    gsInfo << "The L2 Error for pressure error: " << l2error_p << "\n";

    //sol_u.function(0)
    //real_t l2error_p = sol_p.distanceL2(p0);
    //error_list_u.push_back(l2error_u);
    //error_list_p.push_back(l2error_p);
    //real_t h = math::pow( (real_t) refine_bases[0]->size(), -1.0 / refine_bases[0]->dim() );

    // Plot solution in paraview
    //gsFunctionExpr<> gtest("((1-x)+y-2*sqrt(abs(1-x)))*((2-x)+y-2*sqrt(2)*sqrt(2-x))x*y","0.0",2);



    int result = 0;
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( sol_u, "velocityStokes", 10000);
        gsWriteParaview<>( sol_p, "pressureStokes", 10000);
        //gsField<> exact( StokesAssembler.patches(), g, false );
        gsField<> exact( *patches , g , false ) ;
        gsWriteParaview<>( exact, "Stokes2d_exact", 10000);

        // Run paraview
        //result = system("paraview Stokes2dpressure.pvd &");
    }


    freeAll( refine_bases );
    delete basisNorm;
    delete patches;

    gsInfo << "Test is done: Exiting" << "\n";
    return  result;

}
