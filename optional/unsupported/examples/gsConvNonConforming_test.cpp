/** @file gsConvNonConforming.cpp

    @brief Testing convergence for non conforming nested meshes

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris. S. E. Moore
*/

#include <iostream>

#include <gismo.h>
#include <gsUtils/gsNorms.h>

using namespace gismo;

int main(int argc, char *argv[])
{   
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 0;
    index_t numNHref    = 0;
    index_t basisDegree = 0;
    bool plot       = false;

    // Multipatch object
    gsMultiPatch<> mp;

    // Pde
    memory::unique_ptr< gsPoissonPde<> > pde;

    int result = 0;
    std::string fn( "planar/two_squares.xml");
    std::string fn_pde("");
    
    gsCmdLine cmd("Testing a multipatch problem.");
    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("R","NonMatch", "Number of non-matching elements steps to perform before solving",
               numNHref);
    cmd.addInt("p","degree", 
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)", fn);
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
    pde = gsReadFile<>(fn_pde) ;
    gsFunctionExpr<>::uPtr solution = gsReadFile<>(fn_pde);

    // Create Dirichlet boundary conditions for all boundaries
    gsBoundaryConditions<> BCs;
    for (gsMultiPatch<>::const_biterator
             bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        BCs.addCondition( *bit, condition_type::dirichlet, *solution );
    }
    
    // Set up and elevate discretization bases
    gsMultiBasis<> bases(mp);
    if (basisDegree)
        bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate);

    gsInfo<<"Solving "<< *pde << "\n";  

    // NON- matching
    // uniformRefine(NHref) creates N+1 elements
    bases[0].uniformRefine(numNHref-1); //
    bases[1].setDegree(basisDegree+1);

    // Run tests 
    gsMatrix<> tests(numHref+1,6);
    tests.setZero();

    gsSparseSolver<>::CGDiagonal solver;
    
    int i = 0;
    do
    {
        // Setup the Galerkin solvers
        gsPoissonAssembler<real_t> galerkin_dg(mp,bases, BCs,*pde->rhs(),
                                                dirichlet::nitsche, iFace::dg);

        gsInfo<<"Discretization Space for patch 0: \n"<< bases[0] << "\n";
        gsInfo<<"Discretization Space for patch 1: \n"<< bases[1] << "\n";
        gsInfo<<"Penalty constant         : "<< galerkin_dg.penalty(0) << "\n";  
        gsInfo<<"Gauss nodes per direction: "<< "\n";
        gsInfo << "---------------------------------------\n";
        gsInfo<<"System size (dg)         : "<< galerkin_dg.numDofs() << "\n";
        gsInfo << "---------------------------------------\n";

        gsStopwatch time;        
        gsInfo << "Computing solution with patch-wise Disc. Galerkin method..\n";
        time.restart();
        galerkin_dg.assemble();
        solver.compute( galerkin_dg.matrix() );
        gsMatrix<> solVector = solver.solve( galerkin_dg.rhs() );
        //double totalTime_dg = time.stop();
        tests(i,1) = solver.iterations();

        gsField<> sol_dg = galerkin_dg.constructSolution(solVector);
            
        // Collect data
        tests(i,0) = galerkin_dg.numDofs();
	
        tests(i,2)= sol_dg.distanceL2( *solution ) ;
        tests(i,4)= igaFieldDGDistance( sol_dg, *solution,false);

        if (i > 0)
        {
            tests(i,3)= tests(i-1,2) / tests(i,2);
            tests(i,5)= tests(i-1,4) / tests(i,4);
        }
        if ( mp.nPatches() > 1 )

            if (plot && i== numHref)
            {
                // Write approximate and exact solution to paraview files
                gsInfo<<"Plotting in Paraview...\n";
                gsWriteParaview<>( sol_dg, "poisson_problem", 1000, true);
            
                // Run paraview
                if ( mp.nPatches() == 1 )
                    result = system("paraview poisson_problem.vts &");
                else
                    result = system("paraview poisson_problem.pvd &");
            }
        
        bases.uniformRefine();
    } 
    while ( i++ < numHref );
    
    for(i = 1; i<= numHref; ++i)
    {   // Compute convergence rates
        tests(i,3)= math::log(tests(i,3))/std::log(2.0);
        tests(i,5)= math::log(tests(i,5))/std::log(2.0);
    }
        
    gsInfo << "Summary:\n\n";  
    gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        Non-Matching DG            \n";
    gsInfo  << "    Dofs   | Iter     |   L2 error     |    conv. rate    |    DG error  | conv. rate  \n" << tests  << "\n";

    return result;
}
