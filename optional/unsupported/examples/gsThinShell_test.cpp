/** @file gsThinShell_test.cpp

    @brief Example testing thin shell solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Goyal, A. Mantzaflaris
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsThinShell/gsShellAssembler.h>

#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;


// Choose among various shell examples, default = Thin Plate
int main (int argc, char** argv)
{
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 0;
    bool plot       = false;
    index_t testCase    = 1;
    bool nonlinear  = false;
 
    real_t thickness     = 0.05;
    real_t E_modulus     = 1.2E6;
    real_t poisson_ratio = 0.0;
    real_t rho = 1.0;
    gsMultiPatch<> mp;

    int result = 0;
//     std::string fn("thin_plate2.xml");
    
    gsCmdLine cmd("Thin shell plate example.");
    cmd.addInt("t",  "testcase", "Choose a test case " 
               "(1: Thin plate, 2: Scordelis lo roof, 3: Quarter hemisphere, 4: Pinched cylinder)", 
               testCase);
    cmd.addInt("r","hRefine", 
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
//     cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)",fn);
    cmd.addInt("e","degreeElevation", 
               "Number of degree elevation steps to perform on the Geometry's basis before solving", 
               numElevate);
    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
        
    if (testCase == 1)
        gsReadFile<>("thin_plate2.xml", mp);
    else if (testCase == 2)
    {
        thickness = 0.125;
        E_modulus = 4.32E8;
        gsReadFile<>("scordelis_lo_roof.xml", mp);
    }
    else if (testCase == 3)
    {
        thickness = 0.02;
        E_modulus = 6.825E7;
        poisson_ratio = 0.3;
        gsReadFile<>("quarter_hemisphere.xml", mp);
    }
    else if (testCase == 4)
    {
        thickness = 1.5;
        E_modulus = 3E6;
        poisson_ratio = 0.3;
        gsReadFile<>("pinched_cylinder.xml", mp);
    }
    
    if (testCase == 1)
        mp.patch(0).degreeElevate();    // Elevate the degree
    else if (testCase == 2)
    {
        mp.patch(0).degreeElevate(2);   // Elevate the degree
        numHref = 3;
    }
    else if (testCase == 3)
    {
        mp.patch(0).degreeElevate(1);   // Elevate the degree
        numHref = 4;
    }
    else if (testCase == 4)
    {
        mp.patch(0).degreeElevate(1);   // Elevate the degree
        numHref = 4;
    }
    
    // h-refine
    for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();

    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    std::vector< std::pair<patchSide,index_t> > clamped;
    gsBoundaryConditions<> BCs;
    gsVector<> tmp(3);
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    
    // NOTE because neumann data is only present in case 1, else would not work
    tmp << 0, 0, -0.2;
    gsConstantFunction<> neuData(tmp,3);
    
    if (testCase == 1)
    {
        BCs.addCondition(boundary::east, condition_type::neumann,&neuData );
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 2 ); // unknown 2 - z

        // Surface forces
        tmp.setZero();

        // z-clamped side 
        clamped.push_back( std::make_pair(patchSide(0,boundary::west),2) );

        // Point loads
        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (3); load << -0.30, -0.20, -0.2 ;
        pLoads.addLoad(point, load, 0 ); 
    }
    else if (testCase == 2)
    {
        // Diaphragm conditions
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 2 ); // unknown 2 - z
        BCs.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)
        
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 2 ); // unknown 2 - z 
        
        // Surface forces
        tmp << 0, 0, -90;
    }
    else if (testCase == 3)
    {
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 );
        clamped.push_back( std::make_pair(patchSide(0,boundary::east),1) ); //clamp at y
        clamped.push_back( std::make_pair(patchSide(0,boundary::east),2) ); //clamp at z
        
        // Symmetry in y-direction:
        clamped.push_back( std::make_pair(patchSide(0,boundary::west),0) ); //clamp at x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 1 );
        clamped.push_back( std::make_pair(patchSide(0,boundary::west),2) ); //clamp at z
        
        // Surface forces
        tmp.setZero();
        
        // Point loads
        gsVector<> point(2); 
        gsVector<> load (3); 
        point<< 0.0, 0.0 ; load << 1.0, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 ); 
        point<< 1.0, 0.0 ; load << 0.0, -1.0, 0.0 ;
        pLoads.addLoad(point, load, 0 ); 
    }
    else if (testCase == 4)
    {
        // Symmetry in y-direction for back side
        clamped.push_back( std::make_pair(patchSide(0,boundary::north),0) ); //clamp at x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 1 );
        clamped.push_back( std::make_pair(patchSide(0,boundary::north),2) ); //clamp at z
        
        // Diaphragm conditions for left side
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 2 ); // unknown 2 - z
        
        // Symmetry in x-direction: for right side
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 );
        clamped.push_back( std::make_pair(patchSide(0,boundary::east),1) ); //clamp at y
        clamped.push_back( std::make_pair(patchSide(0,boundary::east),2) ); //clamp at z
        
        // Symmetry in z-direction:for the front side
        clamped.push_back( std::make_pair(patchSide(0,boundary::south),0) ); //clamp at x
        clamped.push_back( std::make_pair(patchSide(0,boundary::south),1) ); //clamp at y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 2 );
        
        // Surface forces
        tmp.setZero();
        
        // Point loads
        gsVector<> point(2); point<< 1.0, 1.0 ;
        gsVector<> load (3); load << 0.0, 0.0, -0.25 ;
        pLoads.addLoad(point, load, 0 ); 
    }
    
    gsConstantFunction<> surfForce(tmp,3);
    
    // Construct assembler object
    gsShellAssembler<real_t> assembler(mp, thickness, rho, E_modulus, 
                                       poisson_ratio, BCs, surfForce, clamped, pLoads);
//     gsDebug << "Mass matrix: \n" << assembler.massMatrix() << "\n";

    gsMultiPatch<> solution;
    
    if (!nonlinear)
    {
        // The linear case
          assembler.assemble();
          gsSparseSolver<>::LU solver;
          solver.compute( assembler.matrix() );
          gsMatrix<> solVector = solver.solve( assembler.rhs() );
          assembler.constructSolution(solVector, solution);
    }
    else
    {
        // The nonlinear case
        // Initial solution for the Newton iteration
        gsMultiPatch<> deformed = assembler.patches();
        gsNewtonIterator<real_t> newtonSolver(assembler, deformed);
        //newtonSolver.setMaxIterations(20);
        newtonSolver.solve();

        if ( newtonSolver.converged() )
        {
            gsInfo <<"Converged after "<<newtonSolver.numIterations()
                    <<" iterations with tolerance "<<newtonSolver.tolerance() <<".\n";
        }
        else
            gsInfo <<"Newton iteration did not converge.\n";

        solution = newtonSolver.solution();
    }
    
    // compute the deformation spline
    gsMultiPatch<> deformation = solution;
    deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
    
    //gsInfo <<"Deformation norm       : "<< deformation.patch(0).coefs().norm() <<".\n";
    gsInfo <<"Maximum deformation coef: "
           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";

    if (plot)
    {            
        // Write solution to paraview file
        gsField<> solField(solution, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "thin_shell_def", 1000);
        
        // Run paraview on exit
        result = system("paraview thin_shell_def.pvd &");
    }
    
    return result;
}
