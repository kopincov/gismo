/** @file gsCDR_TrilinosML.cpp

    @brief Convection-diffusion-reaction equation in 2D 

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moeller, R. Tielen
*/

#include <gismo.h>
#include <gismo_dev.h>
#include <gsTrilinos/gsTrilinos.h>

using namespace gismo;

/** @brief
 *  Configure the solver based on the configuration string
 *
 *  The config string must have the form "NAME:TYPE=VALUE", where NAME
 *  is the name of the option, TYPE denotes its type, and VALUE is the
 *  value. Admissible values for TYPE are BOOL, INT, DOUBLE, and
 *  STRING. Multiple config triples are separated by ";", e.g.
 *  "NAME1:TYPE1=VALUE1;NAME2:TYPE2=VALUE2".
 */
template<typename Solver>
void configureSolver(Solver & solver, const std::string & config);

/** @brief
 *  Main routine
 *
 *  To get an overview of the command line arguments run
 *
 *  ./gsCDR_TrilinosML -s <SOLVER> -p <PRECONDITIONER> --verbose
 *
 *  where <SOLVER> is one of "Amesos", "Aztec", and "Belos" and
 *  <PRECONDITIONER> is either "" or "ML". This will run the selected
 *  solver with the default options and parameters. The --verbose flag
 *  tells the solver to output detailed information about the options
 *  and parameters used. This gives a hint of which options and
 *  parameters can be passed to the solver.
 *
 *  Many solvers can be further controlled by passing a specific typer
 *  type, again using the particular default option, e.g.
 *
 *  ./gsCDR_TrilinosML -s "Belos:BiCGStab"
 *
 *  which selects the BiCGStab solver from the Belos package.
 *
 *  More fine-grained control is possible by passing solver-specific
 *  parameters, e.g.
 *
 *  ./gsCDR_TrilinosML -s "Belos:BiCGStab" -B "Convergence Tolerance:DOUBLE=1e-12"
 *
 *  This will set the convergence tolerance to 1e-12 of type
 *  double. The expression follows the CMake pattern
 *  NAME:TYPE=VALUE. Multiple parameters can be separated by ";", e.g.
 *
 *  ./gsCDR_TrilinosML -s "Belos:BiCGStab" -B "Convergence Tolerance:DOUBLE=1e-12;Maximum Iterations:INT=100"
 */
int main(int argc, char* argv[])
{    
    // Parameters
    bool plot = true;
    bool verbose = false;
    bool status = false;
    bool timing = false;
    index_t numRefine = 5;
    index_t numDegree = 1;
    std::string solver = "";
    std::string preconditioner = "";
    std::string amesos_options = "";
    std::string aztec_options = "";
    std::string belos_options = "";
    std::string ml_options = "";

    // Exact solution
    std::string exact = "sin(pi*(x+1/2))*sin(pi*(y+1/2))";

    // Load vector
    std::string load  = "1.2*pi*pi*sin(pi*(x+1/2))*sin(pi*(y+1/2))+0.9*pi*pi*sin(pi*(x+1/2))*sin(pi*(y+1/2))+0.7*pi*pi*cos(pi*(x+1/2))*cos(pi*(y+1/2)) + 0.4*pi*pi*cos(pi*(x+1/2))*cos(pi*(y+1/2)) +0.4*pi*cos(pi*(x+1/2))*sin(pi*(y+1/2))-0.2*pi*sin(pi*(x+1/2))*cos(pi*(y+1/2))+0.3*sin(pi*(x+1/2))*sin(pi*(y+1/2))";

    // Diffusion coefficients
    std::string diffusion = "1.2;-0.7;-0.4;0.9";

    // Convection coefficient
    std::string convection = "0.4;-0.2";

    // Reaction coefficient
    std::string reaction = "0.3";
    
    // Command line argument parser
    gsCmdLine cmd("This file solves the convection-diffusion-reaction equation with different solvers to test the Trilinos solver backend");

    // Add command line arguments
    cmd.addInt("r", "refine", "Number of global refinements", numRefine);
    cmd.addInt("d", "degree", "Number of order elevation steps", numDegree);

    cmd.addSwitch("verbose", "Verbose output", verbose);
    cmd.addSwitch("status", "Status output", status);
    cmd.addSwitch("timing", "Timing output", timing);
    
    cmd.addString("C", "convection", "Expression of convection coefficient", convection);
    cmd.addString("D", "diffusion", "Expression of diffusion coefficient", diffusion);
    cmd.addString("E", "exact", "Expression of exact solution", exact);
    cmd.addString("L", "load", "Expression of load vector", load);
    cmd.addString("R", "reaction", "Expression of reaction coefficient", reaction);
    cmd.addString("p", "preconditioner", "Type of external preconditioner", preconditioner);
    cmd.addString("s", "solver", "Type of solver", solver);

    cmd.addString("A", "amesos", "Additional options for Amesos solver", amesos_options);
    cmd.addString("B", "belos",  "Additional options for Belos solver", belos_options);
    cmd.addString("M", "ml",     "Additional options for ML solver", ml_options);
    cmd.addString("Z", "aztec",  "Additional options for Aztec solver", aztec_options);
    
    // Read parameters from command line
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // Exact solution
    gsFunctionExpr<> sol_exact(exact,2);
    
    // Diffusion, convection and reaction coefficient: 
    gsFunctionExpr<> coeff_diff(util::tokenize(diffusion, ";", 0),
                                util::tokenize(diffusion, ";", 1),
                                util::tokenize(diffusion, ";", 2),
                                util::tokenize(diffusion, ";", 3), 2);
    gsFunctionExpr<> coeff_conv(util::tokenize(convection, ";", 0),
                                util::tokenize(convection, ";", 1), 2);
    gsFunctionExpr<> coeff_reac(reaction, 2);
    gsFunctionExpr<> f(load, 2);

    // Create a square patch of size 1 x 1 and a multibasis
    gsMultiPatch<> patch(*gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0));
    gsMultiBasis<> basis(patch);

    // Apply refinement and order elevation
    for (int i = 0; i < numDegree; ++i)
    {
        basis.degreeElevate();
    }

    for (int i = 0; i < numRefine; ++i)
    {
        basis.uniformRefine();
    }

    // Print some information about the basis
    gsInfo << "Number of patches: \n" << patch.nPatches() << "\n \n";
    gsInfo << "Degree of basis functions: \n" << basis.degree() << "\n \n";
    gsInfo << "Number of basis functions: \n" << basis.totalSize() << "\n \n";  

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;

    // Dirichlet conditions 
    bcInfo.addCondition(boundary::west,  condition_type::dirichlet, &sol_exact);
    bcInfo.addCondition(boundary::east,  condition_type::dirichlet, &sol_exact);
    bcInfo.addCondition(boundary::south, condition_type::dirichlet, &sol_exact);
    bcInfo.addCondition(boundary::north, condition_type::dirichlet, &sol_exact);

    // Define the PDE
    gsConvDiffRePde<real_t> cdrPde(patch, bcInfo, & coeff_diff, & coeff_conv, & coeff_reac, & f);
  
    // Set up the assembler (and a patch to store solution)
    gsCDRAssembler<real_t> pa( cdrPde, basis);
    gsMultiPatch<> mpsol;

    // Assemble
    gsInfo << "Assembling... \n";
    pa.assemble();
    gsInfo << "Done. \n \n";
    
    // Solve the system with solver from library
    gsInfo << "Solving... \n";

    // Define solution vector
    gsVector <> solVector;

#ifdef GISMO_WITH_TRILINOS
    // Initialize the MPI environment and obtain the world communicator
    gsMpiComm comm = gsMpi::init(argc,argv).worldComm();

    trilinos::SparseMatrix t_A(pa.matrix());
    trilinos::Vector       t_b(pa.rhs(), t_A);

    //Get the rank of the processor
    int _rank = comm.rank();
    gsInfo << "Processor "<<_rank <<" has "<< t_b.mySize()
           <<" out of "<<t_b.size()<<" equations.\n";
    comm.barrier();

    if (solver.substr(0,6) == "Amesos")
    {
        /// Sparse direct solver: Amesos
        ///
        /// Valid solver options are:
        /// 
        /// -s "Amesos:Lapack"
        /// -s "Amesos:KLU"
        /// -s "Amesos:Umfpack"
        /// -s "Amesos:Pardiso"
        /// -s "Amesos:Taucs"
        /// -s "Amesos:SuperLU"
        /// -s "Amesos:SuperLUDist"
        /// -s "Amesos:Mumps"
        /// -s "Amesos:Dscpack"
        ///
        /// Note that not all options are available yet and,
        /// moreover, some of them need to be enabled explicitly
        /// in the G+Smo CMake configuration.
        trilinos::solver::AmesosSolver t_solver(t_A,
                                                solver.substr(7) == "Lapack" ?
                                                trilinos::solver::AmesosSolvers::Lapack :
                                                solver.substr(7) == "KLU" ?
                                                trilinos::solver::AmesosSolvers::KLU :
                                                solver.substr(7) == "Umfpack" ?
                                                trilinos::solver::AmesosSolvers::Umfpack :
                                                solver.substr(7) == "Pardiso" ?
                                                trilinos::solver::AmesosSolvers::Pardiso :
                                                solver.substr(7) == "Taucs" ?
                                                trilinos::solver::AmesosSolvers::Taucs :
                                                solver.substr(7) == "SuperLU" ?
                                                trilinos::solver::AmesosSolvers::SuperLU :
                                                solver.substr(7) == "SuperLUDist" ?
                                                trilinos::solver::AmesosSolvers::SuperLUDist :
                                                solver.substr(7) == "Mumps" ?
                                                trilinos::solver::AmesosSolvers::Mumps :
                                                solver.substr(7) == "Dscpack" ?
                                                trilinos::solver::AmesosSolvers::Dscpack :
                                                0);
        configureSolver(t_solver, amesos_options);

        if (verbose) gsInfo << t_solver.currentParams();
      
        /// Compute solution
        const trilinos::Vector & t_vx = t_solver.solve(t_b);
        GISMO_UNUSED(t_vx);
            
        t_solver.getSolution(solVector); // collect solution at Proc 0

        /// Get solver status and timing
        if (status) gsInfo << t_solver.status();
        if (timing) gsInfo << t_solver.timing();
    }
    else if (solver.substr(0,5) == "Aztec" &&
             (preconditioner == "" ||
              (preconditioner.length() >= 5 && preconditioner.substr(0,5) == "Aztec")))
    {
        /// Iterative solver: Aztec with or without built-in preconditioner
        ///
        /// Valid solver options are:
        ///
        /// -s "Aztec" -> Use option "AZ_solver" 
        ///               to specify solver and internal preconditioner
        /// -s "Aztec:Analyze"
        /// -s "Aztec:BiCGStab"
        /// -s "Aztec:CG"
        /// -s "Aztec:CGS"
        /// -s "Aztec:CG_Condnum"
        /// -s "Aztec:Fixed_Pt"
        /// -s "Aztec:Gmres"
        /// -s "Aztec:GmresR"
        /// -s "Aztec:Gmres_Condnum"
        /// -s "Aztec:LU"
        /// -s "Aztec:SLU"
        /// -s "Aztec:SymmLQ"
        /// -s "Aztec:TFQMR"
        ///
        /// Valid preconditioner options are:
        ///
        /// -s "Aztec" -> Use option "AZ_precond" 
        ///               to specify solver and internal preconditioner
        /// -p "Aztec:Dom_Decomp"
        /// -p "Aztec:Jacobi"
        /// -p "Aztec:LS"
        /// -p "Aztec:Neumann"
        /// -p "Aztec:Sym_GS"
        trilinos::solver::AztecSolver t_solver(t_A);
        configureSolver(t_solver, aztec_options);

        /// Set Aztec solver
        if (solver.length() > 5) {
            t_solver.set("AZ_solver", solver.substr(6));
        }

        /// Set Aztec preconditioner
        if (preconditioner.length() > 5) {
            t_solver.set("AZ_precond", preconditioner.substr(6));
        }
        
        if (verbose) gsInfo << t_solver.currentParams();

        /// Compute solution
        const trilinos::Vector & t_vx = t_solver.solve(t_b);
        GISMO_UNUSED(t_vx);
            
        t_solver.getSolution(solVector); // collect solution at Proc 0

        /// Get solver status and timing
        if (status) gsInfo << t_solver.status();
        if (timing) gsInfo << t_solver.timing();
    }
    else if (solver.substr(0,5) == "Aztec" &&
             preconditioner.substr(0,2) == "ML")
    {
        /// Iterative solver: Aztec with ML preconditioner
        ///
        /// Valid solver options are:
        ///
        /// -s "Aztec" -> Use option "AZ_solver" and "AZ_precond"
        ///               to specify solver and internal preconditioner
        /// -s "Aztec:Analyze"
        /// -s "Aztec:BiCGStab"
        /// -s "Aztec:CG"
        /// -s "Aztec:CGS"
        /// -s "Aztec:CG_Condnum"
        /// -s "Aztec:Fixed_Pt"
        /// -s "Aztec:Gmres"
        /// -s "Aztec:GmresR"
        /// -s "Aztec:Gmres_Condnum"
        /// -s "Aztec:LU"
        /// -s "Aztec:SLU"
        /// -s "Aztec:SymmLQ"
        /// -s "Aztec:TFQMR"
        trilinos::solver::AztecSolver t_solver(t_A);
        configureSolver(t_solver, aztec_options);

        /// Set Aztec solver
        if (solver.length() > 5) {
            t_solver.set("AZ_solver", solver.substr(6));
        }
        
        if (verbose) gsInfo << t_solver.currentParams();
      
        /// Multilevel preconditioner
        ///
        /// Valid preconditioner options are:
        ///
        /// -p "ML:SA"
        /// -p "ML:NSSA"
        /// -p "ML:DD"
        /// -p "ML:DDLU"
        /// -p "ML:DDML"
        /// -p "ML:DDMLLU"
        trilinos::solver::MLSolver t_precond(t_A,
                                             preconditioner.substr(3) == "SA" ?
                                             trilinos::solver::MLSolvers::SA :
                                             preconditioner.substr(3) == "NSSA" ?
                                             trilinos::solver::MLSolvers::NSSA :
                                             preconditioner.substr(3) == "DD" ?
                                             trilinos::solver::MLSolvers::DD :
                                             preconditioner.substr(3) == "DDLU" ?
                                             trilinos::solver::MLSolvers::DDLU :
                                             preconditioner.substr(3) == "DDML" ?
                                             trilinos::solver::MLSolvers::DDML :
                                             preconditioner.substr(3) == "DDMLLU" ?
                                             trilinos::solver::MLSolvers::DDMLLU :
                                             0);
        configureSolver(t_precond, ml_options);

        if (verbose) gsInfo << t_precond.currentParams();
      
        /// Set ML-preconditioner
        t_solver.setPreconditioner(t_precond);

        /// Compute solution
        const trilinos::Vector & t_vx = t_solver.solve(t_b);
        GISMO_UNUSED(t_vx);
            
        t_solver.getSolution(solVector); // collect solution at Proc 0

        /// Get solver status and timing
        if (status) gsInfo << t_solver.status();
        if (timing) gsInfo << t_solver.timing();

        /// Get preconditioner status and timing
        if (status) gsInfo << t_precond.status();
        if (timing) gsInfo << t_precond.timing();
    }
    else if (solver.substr(0,5) == "Belos")
    {
        /// Iterative solver: Belos
        ///
        /// Valid solver options are:
        ///
        /// -s "Belos:BiCGStab"
        /// -s "Belos:BlockCG"
        /// -s "Belos:BlockGCRODR"
        /// -s "Belos:BlockGmres"
        /// -s "Belos:FixedPoint"
        /// -s "Belos:GCRODR"
        /// -s "Belos:GmresPoly"
        /// -s "Belos:LSQR"
        /// -s "Belos:Minres"
        /// -s "Belos:PCPG"
        /// -s "Belos:PseudoBlockCG"
        /// -s "Belos:PseudoBlockGmres"
        /// -s "Belos:PseudoBlockStochasticCG"
        /// -s "Belos:PseudoBlockTFQMR"
        /// -s "Belos:RCG"
        /// -s "Belos:TFQMR"
        trilinos::solver::BelosSolver t_solver(t_A,
                                               solver.substr(6) == "BiCGStab" ?
                                               trilinos::solver::BelosSolvers::BiCGStab :
                                               solver.substr(6) == "BlockCG" ?
                                               trilinos::solver::BelosSolvers::BlockCG :
#                                              ifdef Belos_ENABLE_Experimental
                                               solver.substr(6) == "BlockGCRODR" ?
                                               trilinos::solver::BelosSolvers::BlockGCRODR :
#                                              endif
                                               solver.substr(6) == "BlockGmres" ?
                                               trilinos::solver::BelosSolvers::BlockGmres :
                                               solver.substr(6) == "FixedPoint" ?
                                               trilinos::solver::BelosSolvers::FixedPoint :
                                               solver.substr(6) == "GCRODR" ?
                                               trilinos::solver::BelosSolvers::GCRODR :
                                               solver.substr(6) == "GmresPoly" ?
                                               trilinos::solver::BelosSolvers::GmresPoly :
                                               solver.substr(6) == "LSQR" ?
                                               trilinos::solver::BelosSolvers::LSQR :
                                               solver.substr(6) == "Minres" ?
                                               trilinos::solver::BelosSolvers::Minres :
                                               solver.substr(6) == "PCPG" ?
                                               trilinos::solver::BelosSolvers::PCPG :
                                               solver.substr(6) == "PseudoBlockCG" ?
                                               trilinos::solver::BelosSolvers::PseudoBlockCG :
                                               solver.substr(6) == "PseudoBlockGmres" ?
                                               trilinos::solver::BelosSolvers::PseudoBlockGmres :
                                               solver.substr(6) == "PseudoBlockStochasticCG" ?
                                               trilinos::solver::BelosSolvers::PseudoBlockStochasticCG :
                                               solver.substr(6) == "PseudoBlockTFQMR" ?
                                               trilinos::solver::BelosSolvers::PseudoBlockTFQMR :
                                               solver.substr(6) == "RCG" ?
                                               trilinos::solver::BelosSolvers::RCG :
                                               solver.substr(6) == "TFQMR" ?
                                               trilinos::solver::BelosSolvers::TFQMR :
                                               0);
                                             
        configureSolver(t_solver, belos_options);

        if (verbose) gsInfo << t_solver.currentParams();

        /// Compute solution
        const trilinos::Vector & t_vx = t_solver.solve(t_b);
        GISMO_UNUSED(t_vx);
            
        t_solver.getSolution(solVector); // collect solution at Proc 0

        /// Get solver status and timing
        if (status) gsInfo << t_solver.status();
        if (timing) gsInfo << t_solver.timing();
    }
    
    /// ML as a solver is not yet working !!!
    
    // else if (solver == "ML")
    //     {
    //         /// Iterative solver: ML
    //         trilinos::solver::MLSolver t_solver(t_A);
    //         configureSolver(t_solver, ml_options);

    //         if (verbose) gsInfo << t_solver.currentParams();

    //         /// Compute solution
    //         const trilinos::Vector & t_vx = t_solver.solve(t_b);
    //         GISMO_UNUSED(t_vx);
            
    //         t_solver.getSolution(solVector); // collect solution at Proc 0

    //         /// Get solver status and timing
    //         if (status) gsInfo << t_solver.status();
    //         if (timing) gsInfo << t_solver.timing();
    //     }
    else
    {
        gsWarn << "Invalid choice of solver. Exiting.\n";
        gsWarn << "\nValid choices for the solver are:\n";
        gsWarn << "-s \"Amesos:Dscpack\"\n";
        gsWarn << "-s \"Amesos:KLU\"\n";
        gsWarn << "-s \"Amesos:Lapack\"\n";
        gsWarn << "-s \"Amesos:Mumps\"\n";
        gsWarn << "-s \"Amesos:Pardiso\"\n";
        gsWarn << "-s \"Amesos:SuperLUDist\"\n";
        gsWarn << "-s \"Amesos:SuperLU\"\n";
        gsWarn << "-s \"Amesos:Taucs\"\n";
        gsWarn << "-s \"Aztec\"\n";
        gsWarn << "-s \"Aztec:Analyze\"\n";
        gsWarn << "-s \"Aztec:BiCGStab\"\n";
        gsWarn << "-s \"Aztec:CGS\"\n";
        gsWarn << "-s \"Aztec:CG\"\n";
        gsWarn << "-s \"Aztec:CG_Condnum\"\n";
        gsWarn << "-s \"Aztec:Fixed_Pt\"\n";
        gsWarn << "-s \"Aztec:GmresR\"\n";
        gsWarn << "-s \"Aztec:Gmres\"\n";
        gsWarn << "-s \"Aztec:Gmres_Condnum\"\n";
        gsWarn << "-s \"Aztec:LU\"\n";
        gsWarn << "-s \"Aztec:SLU\"\n";
        gsWarn << "-s \"Aztec:SymmLQ\"\n";
        gsWarn << "-s \"Aztec:TFQMR\"\n";
        gsWarn << "-s \"Belos:BiCGStab\"\n";
        gsWarn << "-s \"Belos:BlockCG\"\n";
#       ifdef Belos_ENABLE_Experimental
        gsWarn << "-s \"Belos:BlockGCRODR\"\n";
#       endif
        gsWarn << "-s \"Belos:BlockGmres\"\n";
        gsWarn << "-s \"Belos:FixedPoint\"\n";
        gsWarn << "-s \"Belos:GCRODR\"\n";
        gsWarn << "-s \"Belos:GmresPoly\"\n";
        gsWarn << "-s \"Belos:LSQR\"\n";
        gsWarn << "-s \"Belos:Minres\"\n";
        gsWarn << "-s \"Belos:PCPG\"\n";
        gsWarn << "-s \"Belos:PseudoBlockCG\"\n";
        gsWarn << "-s \"Belos:PseudoBlockStochasticCG\"\n";
        gsWarn << "-s \"Belos:PseudoBlockTFQMR\"\n";
        gsWarn << "-s \"Belos:RCG\"\n";
        gsWarn << "-s \"Belos:TFQMR\"\n";
        gsWarn << "\nValid choices for the preconditioner are:\n";
        gsWarn << "-p \"Aztec\"\n";
        gsWarn << "-p \"Aztec:Dom_Decomp\"\n";
        gsWarn << "-p \"Aztec:Jacobi\"\n";
        gsWarn << "-p \"Aztec:LS\"\n";
        gsWarn << "-p \"Aztec:Neumann\"\n";
        gsWarn << "-p \"Aztec:Sym_GS\"\n";
        gsWarn << "-p \"ML:SA\"\n";
        gsWarn << "-p \"ML:NSSA\"\n";
        gsWarn << "-p \"ML:DD\"\n";
        gsWarn << "-p \"ML:DDLU\"\n";
        gsWarn << "-p \"ML:DDML\"\n";
        gsWarn << "-p \"ML:DDMLLU\"\n";
        
        return 1;
    }
#else
    gsSparseSolver<>::CGDiagonal slv(pa.matrix());
    solVector = slv.solve(pa.rhs());
    gsInfo << "Done in " << slv.iterations() << " iterations. \n \n";
#endif

    // Construct the solution from solVector
    pa.constructSolution(solVector, mpsol);
    gsField<> sol(pa.patches(), mpsol);
  
    // Determine L2 norm of the error
    gsNormL2<real_t> norm(sol, sol_exact);
    norm.compute();
    real_t error_L2 = norm.value();
  
    // Print L2 norm of the error
    gsInfo << "The L2 error is equal to: \n" << error_L2 << "\n \n";
  
    //Plot solution in Paraview
    if (plot)
    {
        // Write approximate and exact solution in Paraview
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(sol, "gsCDR_TrilinosML",5000);
        const gsField<> exact_sol(patch, sol_exact, false);
        gsWriteParaview<>(exact_sol, "gsCDR_TrilinosML_exact", 5000);
          
    }
    return 0;
}


template<typename Solver>
void configureSolver(Solver & solver, const std::string & str)
{
    size_t npos=0, nlen=0;
    while (nlen < str.length())
    {
        // Get token NAME:TYPE=VALUE
        nlen = str.substr(npos).find(";");
        std::string token = str.substr(npos,nlen);
        npos+=nlen+1;
            
        // Split token into NAME
        size_t mlen = token.find_last_of(":");
        std::string name = token.substr(0,mlen);
        size_t mpos = mlen+1;
        
        /// ... TYPE
        mlen = token.substr(mpos).find_last_of("=");
        std::string type = token.substr(mpos,mlen);
        mpos+= mlen+1;
        
        /// ... and VALUE
        std::string value = token.substr(mpos);
        
        if (type == "BOOL")
            solver.set(name, value != "0");
        else if (type == "INT")
            solver.set(name, util::stoi(value));
        else if (type == "DOUBLE")
            solver.set(name, util::stod(value));
        else if (type == "STRING")
            solver.set(name, value);
        else
            GISMO_ERROR("Error : Invalid parameter in solver configuration.");    
    }
}
