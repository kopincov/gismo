/** @file gsMagnetostaticAssembler_test.cpp

    The aim of this file is to test the magnetostatic assembler for IETI 
    on simpler domains than the electric motor
*/

#include <gismo.h>
#include <gsIETI/gsIETIUtils.h>
#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIScaledDirichlet.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsAssembler/gsExpressions.h>
#include <gsAssembler/gsMagnetostaticAssembler.h>
#include <iostream>

using namespace gismo;

gsMultiPatch<> solvePDE(gsMagnetostaticPde<real_t> magPDE, gsMultiPatch<> mp, gsOptionList opt, gsPiecewiseFunction<> f);
std::pair<gsMultiPatch<>, gsMatrix<> > solvePDE(gsMultiPatch<> mp, gsOptionList opt, gsBoundaryConditions<> bcInfo, gsFunctionSet<>* uh, gsFunctionExpr<> Bd, gsPiecewiseFunction<real_t> stepfunc);
void createBoundaryMP(gsMultiPatch<> & boundaryMP);

int main()
{
    gsMultiPatch<> testPatch( gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5) );
  
    gsWriteParaview(testPatch, "initialSquare");

    // Create the options for the problem
    dirichlet::values dVals = dirichlet::homogeneous;
    dirichlet::strategy dStrat = dirichlet::elimination;
    iFace::strategy iStrat = iFace::glue;
    
    gsOptionList opt = gsAssembler<>::defaultOptions();
    opt.setInt("DirichletValues"  , dVals);
    opt.setInt("DirichletStrategy", dStrat );
    opt.setInt("InterfaceStrategy", iStrat  );

    
    typedef gsExprAssembler<real_t>::space space;
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::variable variable;
    typedef gsExprAssembler<real_t>::solution        solution;
    
    // Material coefficients
    // nu_air corresponds to reluctivity in air
    // nu_ferr corresponds to reluctivity in ferromagnetic material
    gsFunctionExpr<> nu_air("10^7/(4 * pi)", 2);
    gsFunctionExpr<> nu_ferr("10^5/(204*pi)", 2);
    
    //gsConstantFunction<> constOne(1.0, 1.0, 2);
    gsConstantFunction<> constOne(1.0, 1.0, 2); //for the tests against PoissonHeterogeneous
    
    // Create piecewise function for the material coefficient 
    gsPiecewiseFunction<real_t> nu;

    for(size_t i = 0; i < testPatch.nPatches(); i++)
    {
      if( i == 0 || i == 3 )
        nu.addPiece(nu_air);
      else
      {
        nu.addPiece(nu_ferr);
	  }
    }

    // Create piecewise function for the load vector, corresponds to M perpenticular
    gsPiecewiseFunction<real_t> f;
    for(size_t i = 0; i < testPatch.nPatches(); i++)
    {
      f.addPiece(constOne);
    }
    

    // Set the boundary conditions
    // Boundary values for the constraints
    // Homogeneous Dirichlet conditions on the inside and the outside of the motor.
    // The remaining boundary conditions are periodic -> Already implemented in the .xml file
    gsConstantFunction<> bc(0.0, 2);
    gsConstantFunction<> nc(0.0, 2);
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &bc);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &bc);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, &bc);
    bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &bc);
    
    bcInfo.addCondition(0, boundary::south, condition_type::neumann, &nc);
    bcInfo.addCondition(2, boundary::east, condition_type::neumann, &nc);
    bcInfo.addCondition(2, boundary::south, condition_type::neumann, &nc);
    bcInfo.addCondition(3, boundary::east, condition_type::neumann, &nc);

    ///
    // Define the number of max. iterations and the error bound
    unsigned numRef = 4;
    
    gsMultiBasis<> multiBasis(testPatch);
    
    
    //==============================Start main loop===========================================
    for(unsigned k = 0; k < numRef; k++)
    {
      // Extracting the patches of the domain
      //std::vector<gsGeometry<>* > patches;
      //patches = testPatch->patches();
      multiBasis.uniformRefine();
    }
    gsMultiPatch<> solution_uh, solution_rh;
    gsMatrix<> u_h, p_h;

    // Solve the state equation
    gsMagnetostaticPde<real_t> magPDE(testPatch, bcInfo, f, nu);
    solution_uh = solvePDE(magPDE, testPatch, opt, f);

    gsField<> vec(testPatch, solution_uh);
    // Plot the solution of the state equation
    gsWriteParaview(vec, "gsMagnetostaticAssembler_test", 1000);
    gsInfo << "Finished process! \n";

    return 0;
}

/*
Function to solve the state equation for a scalar valued function
*/
gsMultiPatch<> solvePDE(gsMagnetostaticPde<real_t> magPDE, gsMultiPatch<> mp, gsOptionList opt, gsPiecewiseFunction<> f)
{
    unsigned numRef = 4;
    gsMultiBasis<> multiBasis(mp);
    for(unsigned k = 0; k < numRef; k++)
    {
        // Extracting the patches of the domain
        //std::vector<gsGeometry<>* > patches;
        //patches = testPatch->patches();
        multiBasis.uniformRefine();
    }
    gsMagnetostaticAssembler<real_t> magAss(magPDE, multiBasis, opt, f);
    
    magAss.assembleFull();
    gsSparseMatrix<> K = magAss.matrix();
    gsSparseSolver<>::LU solver;
    solver.compute(K);
    // !!! magAss.rhs() gives only the 0th piece of the function !!!!!
    gsMatrix<> sol = solver.solve( magAss.getRhsFull() );
    
    gsIETIAssembler<real_t> EM_IETIAss(magAss);

    gsOptionList defaultOpt = gsIETIAssembler<real_t>::defaultOptions();

    EM_IETIAss.setOptions(defaultOpt);
    EM_IETIAss.init();

    EM_IETIAss.assemble();

    //Setup the linearoperator for the iterative solver
    gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(EM_IETIAss);
    solv->init();

    //Setup the preconditioner
    gsScaledDirichletPrecond<real_t>::Ptr precDir = gsScaledDirichletPrecond<real_t>::make(EM_IETIAss);
    //Setup the CG
    gsConjugateGradient<> PCG(solv,precDir);
    PCG.setMaxIterations(100);
    PCG.setCalcEigenvalues(true);

    //Intial guess
    gsMatrix<> solVector(EM_IETIAss.systemSize(),EM_IETIAss.numberRhs());
    solVector.setZero();

        //Solve the IETI system, we obtain the lagrange multipliers
    PCG.solve(solv->getRhs(),solVector);

    //Information about the eigenvalues of the preconditioned system
    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo<<"eigenvalues: "<<"\n"<<eigs.minCoeff()<<" -- "<<eigs.maxCoeff()<<"\n";
    gsInfo<<"Number of iterations: "<<PCG.iterations()<<"\n";

    //Backsubstituting the lagrange multipliers in order to obtain the solution from Ku=f
    gsMatrix<> solution;

    solv->calculateSolution(solVector,solution);
    gsMultiPatch<> mpSol;
    
    gsInfo << "The error in the solutions is: " << (solution - sol).norm() << "\n";
    
    magAss.constructSolution(solution, mpSol);
    
    //gsField<> solIETI = (magAss.constructSolution(solution)); 

    return mpSol;
  
}



