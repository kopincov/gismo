/** @file gsConvergence_test.cpp

    @brief Produce Paraview file output from XML input, fo Visualizing  G+Smo objects

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, S. Moore
*/

#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsTransform.h>
#include <gsCore/gsFuncData.h>



using namespace gismo;


int main(int argc, char *argv[])
{   
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 0;
    index_t minDegree = 0;
    bool plot       = false;

    // Multipatch object
    gsMultiPatch<> mp;

    // Pde
    memory::unique_ptr<gsPoissonPde<> > pde;
    gsFunctionExpr<>::Ptr exactSol;
    
    int result = 0;
    std::string fn("planar/two_squares.xml");
    std::string fn_pde("");
    
    gsCmdLine cmd("Testing a multipatch problem.");
    cmd.addInt("r","hRefine", 
                   "Number of dyadic h-refinement (bisection) steps to perform before solving",
                   numHref);
    cmd.addInt("p","degree",
                   "Degree of the basis functions to use for solving (will elevate or reduce the input)",
                   minDegree);
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)",
                      fn);
    cmd.addString("q","pde","File containing a poisson PDE (.xml)", fn_pde);
    cmd.addInt("e","degreeElevation", 
               "Number of degree elevation steps to perform on the Geometry's basis before solving", 
               numElevate);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    gsReadFile<>(fn, mp);
    // Read PDE (and known exact solution)
    if ( fn_pde.empty() )
        switch ( mp.geoDim() )
        {
        case 1:
            fn_pde = "pde/poisson1d_sin.xml" ;
            break;
        case 2:
            fn_pde = "pde/poisson2d_sin.xml";
            break;
        case 3:
            fn_pde = "pde/poisson3d_sin.xml";
            break;
        default:
            gsInfo <<"Got "<< mp;
        }
    pde = gsReadFile<>(fn_pde);
    exactSol = gsReadFile<>(fn_pde);
    // Boundary conditions
    gsBoundaryConditions<> BCs;
    // Create Dirichlet boundary conditions for all boundaries
    for (gsMultiPatch<>::const_biterator 
             bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        BCs.addCondition( *bit, condition_type::dirichlet, exactSol );
    }
    
    // Set up and elevate discretization bases
    gsMultiBasis<> bases(mp);

    if (minDegree)
        for (short_t k=0; k!=bases.dim(); ++k) // Elevate if requested
        {
            const short_t deg = bases[0].degree(k);
            if ( deg < minDegree)
            {
                const index_t numEl = minDegree - deg;
                bases.degreeIncrease(numEl, k);
            }
        }
    else if (numElevate)
        bases.degreeElevate(numElevate);
//               bases.degreeReduce(1);	

    // Run tests
    gsMatrix<> testsL2(numHref+1,7);
    testsL2.setZero();

    gsMatrix<> testsH1(numHref+1,7);
    testsH1.setZero();
    
    //gsSparseSolver<>::CGDiagonal solver; // identity preconditioner
    gsSparseSolver<>::CGDiagonal solver;
    solver.setTolerance(1e-10);
    //gsSparseSolver<>::CGDiagonal solver;
    //gsSparseSolver<>::BiCGSTABILUT solver( galerkin.matrix() );
    //gsSparseSolver<>::BiCGSTABILUT solver( galerkin.fullMatrix() );

    int i = 0;
    do
    {
        // Setup the assemblers
        gsPoissonAssembler<real_t> galerkin(mp,bases,BCs,*pde->rhs(),
                                            dirichlet::elimination,iFace::glue);
        galerkin.options().setInt("DirichletValues", dirichlet::l2Projection);
        gsPoissonAssembler<real_t> galerkin_wBC(mp,bases,BCs,*pde->rhs(),
                                                dirichlet::nitsche,iFace::glue);
        galerkin_wBC.options().setInt("DirichletValues", dirichlet::l2Projection);
        gsPoissonAssembler<real_t> galerkin_dg(mp,bases,BCs,*pde->rhs(),
                                               dirichlet::nitsche,iFace::dg);
        galerkin_dg.options().setInt("DirichletValues", dirichlet::l2Projection);
        //gsPoissonAssemler<> gp( bvp, bases, elimination, dg  );

        gsInfo<<"Discretization Space for patch 0: \n"<< bases[0] << "\n";
        gsInfo<<"Initial DoFs             : "<< galerkin_dg.numDofs() << "\n";  
        gsInfo<<"Penalty constant         : "<< galerkin_dg.penalty(0) << "\n";  
        gsInfo<<"Gauss nodes per direction: "<< "\n";  
        gsInfo << "---------------------------------------\n";
        gsInfo<<"System  size (elim. BCs)  : "<< galerkin.numDofs() << "\n";
        gsInfo<<"Coupled size (elim. BCs)  : "<< galerkin.system().colMapper(0).coupledSize() << "\n";  
        gsInfo<<"System size (Nitsche BCs) : "<< galerkin_wBC.numDofs() << "\n";  
        gsInfo<<"Coupled size(Nitsche BCs) : "<< galerkin_wBC.system().colMapper(0).coupledSize() << "\n";  
        gsInfo<<"System size (dg)          : "<< galerkin_dg.numDofs() << "\n";
        gsInfo << "---------------------------------------\n";
        
        gsInfo<< "Computing conforming C^0 solution..\n";
        gsStopwatch time;
        galerkin.assemble();
        const double assTime = time.stop();
        gsInfo << "Assembly time (elim. BCs): " << assTime << " s" << "\n";
        time.restart();
        
        gsMatrix<> solVector;

        double solvTime(0.0);
        if ( galerkin.numDofs() )
        {
            solver.compute( galerkin.matrix() );
            solVector = solver.solve( galerkin.rhs() );
            solvTime = time.stop();
            gsInfo << "Solving time (elim. BCs): " << solvTime << " s" << "\n";
            gsInfo << "residual error: " << solver.error() << "\n";
            gsInfo << "    iterations: " << solver.iterations() << "\n";
            gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
        }
        gsField<> sol = galerkin.constructSolution(solVector);
            
        gsInfo << "Computing solution with weakly imposed BCs..\n";
        time.restart();
        galerkin_wBC.assemble();
        const double assTime_wBC = time.stop();
        gsInfo << "Assembly time (weak  BCs): " << assTime_wBC << " s" << "\n";
        time.restart();
        solver.compute( galerkin_wBC.matrix() );
        solVector = solver.solve( galerkin_wBC.rhs() );
        const double solvTime_wBC = time.stop();
        gsInfo << "Solving time (weak  BCs): " << solvTime_wBC << " s" << "\n";
        gsInfo << "residual error: " << solver.error() << "\n";
        gsInfo << "    iterations: " << solver.iterations() << "\n";
        gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
        gsField<> sol_wBC = galerkin_wBC.constructSolution(solVector);
        
        gsField<> sol_dg;
        double solvTime_dg(0), assTime_dg(0);
        if ( mp.nPatches() > 1 )
        {
            gsInfo << "Computing solution with patch-wise Disc. Galerkin method..\n";
            time.restart();
            galerkin_dg.assemble();
            assTime_dg = time.stop();
            gsInfo << "Assembly time (full DG): " << assTime_dg << " s" << "\n";
            time.restart();
            solver.compute( galerkin_dg.matrix() );
            solVector = solver.solve( galerkin_dg.rhs() );
            solvTime_dg = time.stop();
            gsInfo << "Solving time (full DG): " << solvTime_dg << " s" << "\n";
            gsInfo << "residual error: " << solver.error() << "\n";
            gsInfo << "    iterations: " << solver.iterations() << "\n";
            gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
            sol_dg = galerkin_dg.constructSolution(solVector);
        }

        // Collect data
        testsL2(i,0)  =
            testsH1(i,0)= bases.totalSize();
        if ( exactSol )
        {
            testsL2(i,1)= sol.distanceL2    ( *exactSol ) ;
            testsL2(i,3)= sol_wBC.distanceL2( *exactSol ) ;
            
            testsH1(i,1) = sol.distanceH1    ( *exactSol ) ;
            testsH1(i,3) = sol_wBC.distanceH1( *exactSol ) ;
        }
                                          
        if ( mp.nPatches() > 1 )
        {
            testsL2(i,5)= sol_dg.distanceL2( *exactSol ) ;
            testsH1(i,5)= igaFieldDGDistance(sol_dg, *exactSol, false);
        }

        if (i > 0)
        {
            testsL2(i,2)= testsL2(i-1,1) / testsL2(i,1);
            testsL2(i,4)= testsL2(i-1,3) / testsL2(i,3);
            testsL2(i,6)= testsL2(i-1,5) / testsL2(i,5);

            testsH1(i,2)= testsH1(i-1,1) / testsH1(i,1);
            testsH1(i,4)= testsH1(i-1,3) / testsH1(i,3);
            testsH1(i,6)= testsH1(i-1,5) / testsH1(i,5);
        }
        
        gsInfo << "Total time (elim. BCs): " << assTime+solvTime   << " s" << "\n";
        gsInfo << "Total time (weak  BCs): " << assTime_wBC+solvTime_wBC << " s" << "\n";
        if ( mp.nPatches() > 1 )
            gsInfo << "Total time (full DG)  : " << assTime_dg+solvTime_dg  << " s" << "\n";    
        gsInfo << "---------------------------------------\n";    
        gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING      |        WEAK BCs        |         FULL DG            \n";
        gsInfo << "    Dofs   |  L2 error  | err. ratio|  L2 error  | err. ratio|  L2 error  | err. ratio    \n" << testsL2.row(i)  << "\n";  
        gsInfo << "           |  H1 error  | err. ratio|  H1 error  | err. ratio|  H1 error  | err. ratio    \n" << testsH1.row(i)  << "\n";  

        if (plot && i== numHref)
        {
            // Write approximate and exact solution to paraview files
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>( sol, "poisson_problem", 1000);

            const gsField<> exact( mp, *exactSol, false );
            gsWriteParaview<>(exact, "poisson_exact", 1000);
            
            // Run paraview
                result = system("paraview poisson_problem.pvd &");
        }
        
        bases.uniformRefine();
    } 
    while ( i++ < numHref );
    
    for(i = 1; i<= numHref; ++i)
    {   // Compute convergence rates
        testsL2(i,2)= math::log(testsL2(i,2))/std::log(2.0);
        testsL2(i,4)= math::log(testsL2(i,4))/std::log(2.0);
        testsL2(i,6)= math::log(testsL2(i,6))/std::log(2.0);

        testsH1(i,2)= math::log(testsH1(i,2))/std::log(2.0);
        testsH1(i,4)= math::log(testsH1(i,4))/std::log(2.0);
        testsH1(i,6)= math::log(testsH1(i,6))/std::log(2.0);
    }
    
    if ( mp.nPatches() == 1 )
    {
        testsL2.col(6).setZero();
        testsH1.col(6).setZero();
    }
    
    gsInfo << "Summary:\n\n";
    gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING      |        WEAK BCs        |         FULL DG            \n";
    gsInfo << "    Dofs   |  L2 error  | conv. rate|  L2 error  | conv. rate|  L2 error  | conv. rate    \n" << testsL2  << "\n";


    gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING      |        WEAK BCs        |         FULL DG            \n";
    gsInfo << "    Dofs   |  H1 error  | conv. rate|  H1 error  | conv. rate|  H1 error  | conv. rate    \n" << testsH1  << "\n";


/*
    gsDebugVar( exactSol );
    gsMatrix<> x(2,1);
    x << 0.2, 0.3;
    gsDebugVar( pde->solution()->eval(x).transpose() );
    gsDebugVar( pde->solution()->deriv(x).transpose() );
*/

    return result;
}
