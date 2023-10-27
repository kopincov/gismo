/** @file ElectricMotorOpti_IETI.cpp

    File to optimize the magnetic flux in an electric motor with the IETI method
*/

# include <gismo.h>
# include <gsIETI/gsIETIUtils.h>
# include <gsIETI/gsIETIAssembler.h>
# include <gsIETI/gsIETISolver.h>
# include <gsIETI/gsIETIScaledDirichlet.h>
# include <gsAssembler/gsMagnetostaticAssembler.h>
# include <gsAssembler/gsMagnetostaticAdjointAssembler.h>
# include <gsAssembler/gsMagnetostaticShapeDerivAssembler.h>
# include <iostream>
#include <gsAssembler/gsMagnetostaticShapeDerivAssembler_decoupled.h>

#include <gsOptimizer/gsOptElectricMotor.h>

//static const real_t pi = (3.141592653589793238462643383279502);

using namespace gismo;

gsMultiPatch<> solvePDE(gsMagnetostaticPde<real_t> magPDE, gsMultiPatch<> mp, gsOptionList opt, gsPiecewiseFunction<real_t> f);
gsMultiPatch<> solveAdjoint(gsMagnetostaticAdjointPde<real_t> magAss, gsOptionList opt, gsBoundaryConditions<> bcInfo);
gsMultiPatch<> calcShapeDerivative(gsMagnetostaticShapeDerivPde<real_t> shapeD, gsOptionList opt, gsBoundaryConditions<> bcInfo, real_t & oldNorm);
std::pair<std::vector<bool>, std::vector<bool> > calculateCriticalAngle(std::vector<gsGeometry<>* > & patches);
void createInterfaceTopology(std::vector<std::pair<unsigned, unsigned> > & interfaces);
void createInterfacePatch(gsMultiPatch<> & InterfaceMP);

int main()
{
  std::string in("motor_conforming.xml");
  gsFileData<> input(in);
  
  if(input.has< gsMultiPatch<> >() )
  {
    gsMultiPatch<> mp;
    input.getFirst(mp);

    //gsWriteParaview(mp, "motorInit");

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
    // nu_mag corresponds to reluctivity in magnet
    gsFunctionExpr<> nu_air("10^7/(4 * pi)", 2);
    gsFunctionExpr<> nu_ferr("10^5/(204*pi)", 2);
    gsFunctionExpr<> nu_mag("10^7/(4*pi*1.086)", 2);

    //gsFunctionExpr<> nu_air("1.", 2);
    //gsFunctionExpr<> nu_ferr("1.", 2);
    //gsFunctionExpr<> nu_mag("1.", 2);

    // Source functions, also piecewise defined, already orthogonal
    gsFunctionExpr<> f_p13("0.9239*1.28*"+nu_mag.expression(0), "-0.3827*1.28*"+nu_mag.expression(0), 2);
    gsFunctionExpr<> f_p14("-0.3827*1.28*"+nu_mag.expression(0), "0.9239*1.28*"+nu_mag.expression(0), 2);
    gsConstantFunction<> constZero(0.0, 0.0, 2);

    // Constant function for the penalty function in the lhs of the auxiliary problem
    gsFunctionExpr<> penalizeDesignDom("1.0", 2);
    gsFunctionExpr<> penalizeRotorDom("1000.0", 2);
    
    // Create piecewise function for the material coefficient 
    gsPiecewiseFunction<real_t> nu;

    for(size_t i = 0; i < mp.nPatches(); i++)
    {
      if(i == 63 || i == 67 || i == 74 || i == 75)
	nu.addPiece(nu_mag);
      else
      {
	//if(i == 0 || i == 2 || i == 5 || i == 7 || i == 8 || i == 9 || i == 10 || i == 11 || i == 12 || i == 15 || i == 18 || i == 19 || i == 22 || i == 26) //old motor model
	if(i == 86 || i == 87 || i == 88 || i == 89 || i == 90 || i == 91 || i == 92 ||
	   i == 79 || i == 80 || i == 81 || i == 82 || i == 83 || i == 84 || i == 85 ||
	   i == 71 || i == 72 || i == 73 ||
	   i == 78 ||
	   i == 66 ||
	   i == 54 || i == 65 ||
	   i == 45 || i == 46 || i == 47 || i == 48 || i == 49 || i == 50 || i == 51 || i == 52 || i == 53 || //corresponds to ferromagnetic ring
	   i == 8 || i == 9 || i == 10 || i == 11 || i == 12 || i == 13 || i == 14 || i == 15 || i == 25 ||
	   i == 62 || i == 77 ||
	   i == 58 || i == 69 || 
	   i == 60 || //corresponds to upper design domain
	   i == 56) //corresponds to lower design domain
	  nu.addPiece(nu_ferr);
        else
	  nu.addPiece(nu_air);
      }
    }

    // Create piecewise function for the load vector, corresponds to M perpenticular
    gsPiecewiseFunction<real_t> f;
    for(size_t i = 0; i < mp.nPatches(); i++)
    {
      if(i == 74 || i == 75)
	f.addPiece(f_p13);
      else
	if(i == 63 || i == 67)
	  f.addPiece(f_p14);
	else
	  f.addPiece(constZero);
    }

    // Create the penalty functon 
    gsPiecewiseFunction<real_t> a;
    for(size_t i = 0; i < mp.nPatches(); i++)
    {
      if(i == 55 || i == 56 || i == 57 || i == 59 || i == 60 || i == 61)
	a.addPiece(penalizeDesignDom);
      else
	a.addPiece(penalizeRotorDom);
    }

    // Set the boundary conditions
    // Boundary values for the constraints
    gsConstantFunction<> bc(0.0, 2);
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(24, boundary::north, condition_type::dirichlet, &bc);
    bcInfo.addCondition(7, boundary::east, condition_type::dirichlet, &bc);
    bcInfo.addCondition(6, boundary::south, condition_type::dirichlet, &bc);
    bcInfo.addCondition(5, boundary::south, condition_type::dirichlet, &bc);
    bcInfo.addCondition(4, boundary::north, condition_type::dirichlet, &bc);
    bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &bc);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, &bc);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &bc);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &bc);
    // end outer boundary
    bcInfo.addCondition(86, boundary::south, condition_type::dirichlet, &bc);
    bcInfo.addCondition(87, boundary::north, condition_type::dirichlet, &bc);
    bcInfo.addCondition(88, boundary::north, condition_type::dirichlet, &bc);
    bcInfo.addCondition(89, boundary::south, condition_type::dirichlet, &bc);
    bcInfo.addCondition(90, boundary::north, condition_type::dirichlet, &bc);
    bcInfo.addCondition(91, boundary::east, condition_type::dirichlet, &bc);
    bcInfo.addCondition(92, boundary::south, condition_type::dirichlet, &bc);
    // end inner boundary

    // Set the boundary conditions
    // Boundary values for the final problem
    gsConstantFunction<> bc_2D(0.0, 0.0, 2);
    gsBoundaryConditions<> bcInfo_2D;
    bcInfo_2D.addCondition(24, boundary::north, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(7, boundary::east, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(6, boundary::south, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(5, boundary::south, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(4, boundary::north, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(3, boundary::north, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(2, boundary::north, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(1, boundary::north, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(0, boundary::north, condition_type::dirichlet, &bc_2D);
    // end outer boundary
    bcInfo_2D.addCondition(86, boundary::south, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(87, boundary::north, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(88, boundary::north, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(89, boundary::south, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(90, boundary::north, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(91, boundary::east, condition_type::dirichlet, &bc_2D);
    bcInfo_2D.addCondition(92, boundary::south, condition_type::dirichlet, &bc_2D);
    // end inner boundary

    // Define the optimal function
    gsFunctionExpr<> Bd("0.005*sin( 4. *  atan2(y,x) )", 2);

    //Define the unit functions
    gsFunctionExpr<> unit1("1.0", "0.0", 2);
    gsFunctionExpr<> unit2("0.0", "1.0", 2);

    ///
    // Define the number of max. iterations and the error bound
    unsigned maxiter = 10;

    // Initialize oldGrad, newGrad to check the norm of the gradient and the initial length for line search
    real_t newGrad = 0, oldGrad = 0;

    real_t h = 1.0;
    
    //==============================Start main loop===========================================
    for(unsigned k = 0; k < maxiter; k++)
    {
      
      // Extracting the patches of the domain
      std::vector<gsGeometry<>* > patches;
      patches = mp.patches();

      // Create the Multipatch domain with just the required interface
      gsMultiPatch<> InterfaceMP(mp);
      createInterfacePatch(InterfaceMP);
      
      // Get the basis of the current domain
      gsMultiBasis<> multiBasis(mp);
    
      gsMultiPatch<> solution_uh, solution_ph, shapeGradient;
      gsMatrix<> u_h, p_h;
      
      // Compute the coefficients of the solution Vector
      gsMagnetostaticPde<real_t> magPDE(mp, bcInfo, f, nu);
      solution_uh = solvePDE(magPDE, mp, opt, f);
      gsField<> vec(mp, solution_uh);
      //gsWriteParaview(vec, "state");
      
      // Compute the solution of the adjoint equation
      gsMagnetostaticAdjointPde<real_t> magAdjPDE(mp, bcInfo, solution_uh, Bd, nu);
      solution_ph = solveAdjoint(magAdjPDE, opt, bcInfo);
   
      // Initialize the solver for the shape gradient
      gsMagnetostaticShapeDerivPde<real_t> shapeDerivPDE(mp, bcInfo, a, solution_uh, solution_ph, nu, f, 2);
      shapeGradient = calcShapeDerivative(shapeDerivPDE, opt, bcInfo, newGrad);


      // Algorithm from the paper
      gsExprEvaluator<real_t> J_old;
      J_old.setIntegrationElements(multiBasis);
      geometryMap G_ref_old = J_old.getMap(mp);
      variable u_df = J_old.getVariable(Bd, G_ref_old );
      variable u_f = J_old.getVariable(solution_uh );
      J_old.integralInterface( (((fjac(u_f).tr() * jac(G_ref_old).ginv()) * (tv(G_ref_old)/tv(G_ref_old).norm())).tr() - u_df).sqr() * nv(G_ref_old).norm(), InterfaceMP.interfaces());
      real_t currVal = J_old.value();
      real_t newVal = 0;

      // Get the maximum value of the solution vector. Avoid division by 0
      // M = sqrt( \nabla J(row, 0)^2 + \nabla J(row, 1)^2 )
      real_t M = 0;

      for (size_t i = 0; i < patches.size(); i++)
      {
	    for (index_t row = 0; row < patches[i]->coefs().rows(); row++)
	    {
	        real_t n = sqrt(pow(shapeGradient[i].coefs()(row, 0), 2) + pow(shapeGradient[i].coefs()(row, 1), 2));

	        if (n > M)
	        M = n;
	    }
      }

      if(M == 0)
	    M = 1;

      unsigned t = 0;

      if (newGrad > oldGrad)
	    h *= 1.0;

      //h = 0.05;

      while(currVal <= newVal || t == 0)
      {
	    gsMultiPatch<> mpCopy(mp);

        std::vector<gsGeometry<>* > patchesCopy;
        patchesCopy = mpCopy.patches();

        for (size_t i = 0; i < patches.size(); i++)
        {
          gsMatrix<> cf = (*(patchesCopy[i])).coefs();

          //if(i != 86 && i != 92 && i != 0 && i != 24 && i != 16 && i != 26 && i != 27 && i != 35 && i != 36 && i != 44 && i != 45 && i != 53 && i != 54 && i != 62 && i != 65 && i != 77 && i != 79 && i != 85 && i != 8 && i != 25) {
              for (index_t row = 0; row < patches[i]->coefs().rows(); row++) {
                  cf(row, 0) = cf(row, 0) - (h / (pow(2, t))) *
                                            (shapeGradient[i].coefs()(row, 0) / M); // usually M, but 1 for test reasons
                  cf(row, 1) = cf(row, 1) - (h / (pow(2, t))) * (shapeGradient[i].coefs()(row, 1) / M); // M
              }
          //}
	        (patchesCopy)[i]->setCoefs(cf);

	    }

        //mpCopy.closeGaps(); // TODO: try to set the periodic interfaces explicitely again
        /*mp.clearTopology();
        for(int i = 0; i < mp.nInterfaces(); i++)
            mpCopy.addInterface(mp.interfaces()[i]);

        for(int i = 0; i < mp.nBoundary(); i++)
            mpCopy.addBoundary(mp.boundaries()[i]);
            */

        gsMultiBasis<> multiBasisCopy(mpCopy);

        gsMagnetostaticPde<real_t> newMagPDE(mpCopy, bcInfo, f, nu); //mpCopy
        gsMultiPatch<> solution_new = solvePDE(newMagPDE, mpCopy, opt, f); //mpCopy

          // Create a topology with just the specific interfaces
          gsMultiPatch<> InterfaceMPCopy(mpCopy);
          createInterfacePatch(InterfaceMPCopy);

        gsExprEvaluator<real_t> J_new;
        geometryMap G_ref_new = J_new.getMap(mpCopy); //mpCopy
        variable u_dnew = J_new.getVariable(Bd, G_ref_new);
        J_new.setIntegrationElements(multiBasisCopy); //multiBasisCopy
        variable u_new = J_new.getVariable(solution_new );
        J_new.integralInterface( (((fjac(u_new).tr() * jac(G_ref_new).ginv()) * (tv(G_ref_new)/tv(G_ref_new).norm())).tr() - u_dnew).sqr() * nv(G_ref_new).norm(), InterfaceMPCopy.interfaces() );
        newVal = J_new.value();
	
	    gsInfo << "New functional value: " << J_new.value() << " and old one: " <<  J_old.value() <<  "\n";
	    //gsInfo << "Gradient: " << solVector.norm() << "\n";
      
        if( J_new.value() < J_old.value())
        {
            mp = mpCopy;
            break;
        }
        gsInfo << "T: " << t << "\n";
        t++;

        if (t == 50)
        {
          k = maxiter - 1;
          break;
        }

      }

      oldGrad = newGrad;
      
      if (k == maxiter-1)
      {
        gsInfo << "Old: " << currVal << "    New: " << newVal << "\n";
        gsInfo << "Finished: " << currVal - newVal << "\n";
        gsInfo << "Iterations: " << k << "\n";
	
        //gsField<> vec(mp, shapeGradient);
        //gsWriteParaview(vec, "Gradient");
	
	    break;
      }

      if ( k == 0 || k == 1 || k == 2 || k == 4 || k == 9 || k==14 || k == 19)
      {	
	    std::stringstream num;
	    num << k;

	    gsField<> vec(mp, shapeGradient);
	    gsWriteParaview(vec, "Gradient"+num.str() );
      }
      //gsField<> vec(mp, shapeGradient);
      //gsWriteParaview(vec, "Gradient");

    }

  }

  return 0;
}

/*
Function to solve the state equation for a scalar valued function
*/
gsMultiPatch<> solvePDE(gsMagnetostaticPde<real_t> magPDE, gsMultiPatch<> mp, gsOptionList opt, gsPiecewiseFunction<real_t> f)
{
    gsMultiBasis<> multiBasis(mp);

    //for(size_t i = 0; i < 3; i++)
    //    multiBasis.uniformRefine();

    gsMagnetostaticAssembler<real_t> magAss(magPDE, multiBasis, opt, f);

    gsStopwatch time;
    magAss.assembleFull();
    gsSparseMatrix<> K = magAss.matrix();
    //gsInfo << "Eigenvalues of the full matrix: \n";
    //gsInfo << K.toDense().eigenvalues().real().minCoeff() << "\n";
    //gsInfo<<"Start calculating eigenvalues, #dofs: "<<K.cols()<<"\n";
    //Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<real_t> > eigenSolver(K.cols());
    //eigenSolver.compute(K);
    //gsInfo<<"Condition number of K:"<<eigenSolver.eigenvalues().maxCoeff()/eigenSolver.eigenvalues().minCoeff()<<"\n\n";

    gsSparseSolver<>::LU solver;
    solver.compute(K);
    // !!! magAss.rhs() gives only the 0th piece of the function !!!!!
    gsMatrix<> sol = solver.solve( magAss.getRhsFull() );
    gsInfo << "Time needed to solve the equation on this domain with a direct solver: " << time << "\n";

    time.restart();
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
    time.stop();
    gsInfo << "Time needed to solve the equation on this domain with the IETI solver: " << time << "\n";

    //Information about the eigenvalues of the preconditioned system
    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo<<"eigenvalues: "<<"\n"<<eigs.minCoeff()<<" -- "<<eigs.maxCoeff()<<"\n";
    gsInfo<<"Number of iterations: "<<PCG.iterations()<<"\n";

    //Backsubstituting the lagrange multipliers in order to obtain the solution from Ku=f
    gsMatrix<> solution;
    gsMultiPatch<> result;

    solv->calculateSolution(solVector,solution);

    magAss.constructSolution(solution, result);

    return result;


}

/*
Function to solve the adjoint equation for a scalar valued function
 */
gsMultiPatch<> solveAdjoint(gsMagnetostaticAdjointPde<real_t> magPDE,
                            gsOptionList                      opt,
                            gsBoundaryConditions<> ) // bcInfo)
{
    gsMultiBasis<> multiBasis(magPDE.domain());

    //for(size_t i = 0; i < 3; i++)
    //    multiBasis.uniformRefine();

    std::vector<std::pair<unsigned, unsigned> > ints;
    createInterfaceTopology(ints);
    gsMagnetostaticAdjointAssembler<real_t> magAss(magPDE, multiBasis, opt, ints);

    gsStopwatch time;

    magAss.assemble();
    gsSparseMatrix<> K = magAss.matrix();
    gsSparseSolver<>::LU solver;
    solver.compute(K);
    // !!! magAss.rhs() gives only the 0th piece of the function !!!!!
    gsMatrix<> sol = solver.solve( magAss.getRhsFull() );
    gsInfo << "Time needed to solve the equation on this domain with a direct solver: " << time << "\n";

    time.restart();
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
    time.stop();
    gsInfo << "Time needed to solve the equation on this domain with the IETI solver: " << time << "\n";

    //Information about the eigenvalues of the preconditioned system
    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo<<"eigenvalues: "<<"\n"<<eigs.minCoeff()<<" -- "<<eigs.maxCoeff()<<"\n";
    gsInfo<<"Number of iterations: "<<PCG.iterations()<<"\n";

    //Backsubstituting the lagrange multipliers in order to obtain the solution from Ku=f
    gsMatrix<> solution;
    gsMultiPatch<> result;

    solv->calculateSolution(solVector,solution);
    //gsInfo << "The difference between the solution techniques is: " << (solution - sol).norm() << "\n";

    magAss.constructSolution(solution, result);
    return result;


}

gsMultiPatch<> calcShapeDerivative(gsMagnetostaticShapeDerivPde<real_t> shapeD,
                                   gsOptionList                         opt,
                                   gsBoundaryConditions<>               /*bcInfo*/,
                                   real_t                             & oldNorm)
{
    gsMultiBasis<> multiBasis(shapeD.domain());

    //for(size_t i = 0; i < 1; i++)
    //    multiBasis.uniformRefine();

    gsMagnetostaticShapeDerivAssembler_decoupled<real_t> magAss(shapeD, multiBasis, opt);

    gsStopwatch time;
    magAss.assemble();
    gsSparseMatrix<> K = magAss.matrix();
    gsSparseSolver<>::LU solver;
    solver.compute(K);
    // !!! magAss.rhs() gives only the 0th piece of the function !!!!!
    gsMatrix<> sol = solver.solve(magAss.getRhsFull());
    gsInfo << "Time needed to solve the equation on this domain with a direct solver: " << time << "\n";

    time.restart();
    gsIETIAssembler<real_t> EM_IETIAss(magAss);

    gsOptionList defaultOpt = gsIETIAssembler<real_t>::defaultOptions();
    defaultOpt.setInt("nRhs", 2);
    defaultOpt.setString("Scaling", "coeff");
    defaultOpt.setString("Strategy", "C");

    EM_IETIAss.setOptions(defaultOpt);
    EM_IETIAss.init();

    EM_IETIAss.assemble();

    //Setup the linearoperator for the iterative solver
    gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(EM_IETIAss);
    solv->init();

    //Setup the preconditioner
    gsScaledDirichletPrecond<real_t>::Ptr precDir = gsScaledDirichletPrecond<real_t>::make(EM_IETIAss);
    //Setup the CG
    gsConjugateGradient<> PCG(solv, precDir);
    PCG.setMaxIterations(100);
    PCG.setCalcEigenvalues(true);

    //Intial guess
    gsMatrix<> solVector(EM_IETIAss.systemSize(), EM_IETIAss.numberRhs());
    solVector.setZero();

    //Solve the IETI system, we obtain the lagrange multipliers
    for (int i = 0; i < 2; i++) {
        gsMatrix<> tmpSolution;
        PCG.solve(solv->getRhs(i), tmpSolution);
        solVector.col(i) = tmpSolution;
    }
    time.stop();
    gsInfo << "Time needed to solve the equation on this domain with the IETI solver: " << time << "\n";

    //Information about the eigenvalues of the preconditioned system
    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo << "eigenvalues: " << "\n" << eigs.minCoeff() << " -- " << eigs.maxCoeff() << "\n";
    gsInfo << "Number of iterations: " << PCG.iterations() << "\n";

    //Backsubstituting the lagrange multipliers in order to obtain the solution from Ku=f
    gsMatrix<> solution;
    gsMultiPatch<> result;

    solv->calculateSolution(solVector, solution);
    gsInfo << "The difference between the solution techniques is: " << (solution - sol).norm() << "\n";
    oldNorm = solution.norm();

    magAss.constructSolution(solution, result, 0);

    //for(size_t n = 0; n < result.nPatches(); n++)
    //    gsInfo << "The coefficients: " << result[n].coefs() << "\n";

/*
    gsVector<index_t > unk(2);
    unk(0) = 0;
    unk(1) = 1;
    magAss.constructSolution(solution, result, unk);
*/
    return result;
}

// For the vector valued problem
/*
gsMultiPatch<> calcShapeDerivative(gsMagnetostaticShapeDerivPde<real_t> shapeD, gsOptionList opt, gsBoundaryConditions<> bcInfo, real_t & oldNorm)
{
    gsMultiBasis<> multiBasis(shapeD.domain());

    for(size_t i = 0; i < 1; i++)
        multiBasis.uniformRefine();

    gsMagnetostaticShapeDerivAssembler<real_t> magAss(shapeD, multiBasis, opt);

    gsStopwatch time;
    magAss.assemble();
    gsSparseMatrix<> K = magAss.matrix();
    gsSparseSolver<>::LU solver;
    solver.compute(K);
    // !!! magAss.rhs() gives only the 0th piece of the function !!!!!
    gsMatrix<> sol = solver.solve( magAss.getRhsFull() );
    gsInfo << "Time needed to solve the equation on this domain with a direct solver: " << time << "\n";

    time.restart();
    gsIETIAssembler<real_t> EM_IETIAss(magAss);

    gsOptionList defaultOpt = gsIETIAssembler<real_t>::defaultOptions();
    defaultOpt.setString("Scaling", "none"); //stiff
    defaultOpt.setString("Strategy", "A");

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
    time.stop();
    gsInfo << "Time needed to solve the equation on this domain with the IETI solver: " << time << "\n";

    //Information about the eigenvalues of the preconditioned system
    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo<<"eigenvalues: "<<"\n"<<eigs.minCoeff()<<" -- "<<eigs.maxCoeff()<<"\n";
    gsInfo<<"Number of iterations: "<<PCG.iterations()<<"\n";

    //Backsubstituting the lagrange multipliers in order to obtain the solution from Ku=f
    gsMatrix<> solution;
    gsMultiPatch<> result;

    solv->calculateSolution(solVector,solution);
    gsInfo << "The difference between the solution techniques is: " << (solution - sol).norm() << "\n";
    oldNorm = solution.norm();

    magAss.constructSolution(solution, result);

    //for(size_t n = 0; n < result.nPatches(); n++)
    //    gsInfo << "The coefficients on patch: " << n << "\n" << result[n].coefs() << "\n";


    //gsVector<index_t > unk(2);
    //unk(0) = 0;
    //unk(1) = 1;
    //magAss.constructSolution(solution, result, unk);

    return result;
}
*/

void createInterfaceTopology(std::vector<std::pair<unsigned, unsigned > > & interfaces)
{
    interfaces.clear();

    std::pair<unsigned, unsigned> interface_pair(44, 27);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(43, 28);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(42, 29);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(41, 30);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(40, 31);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(32, 39);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(38, 33);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(37, 34);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(36, 35);
    interfaces.push_back(interface_pair);

}

void createInterfacePatch(gsMultiPatch<> & InterfaceMP)
{
    InterfaceMP.clearTopology();
    InterfaceMP.addInterface(&InterfaceMP.patch(44), 2, &InterfaceMP.patch(27), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(43), 4, &InterfaceMP.patch(28), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(42), 4, &InterfaceMP.patch(29), 4);
    InterfaceMP.addInterface(&InterfaceMP.patch(41), 4, &InterfaceMP.patch(30), 4);
    InterfaceMP.addInterface(&InterfaceMP.patch(40), 4, &InterfaceMP.patch(31), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(32), 3, &InterfaceMP.patch(39), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(38), 4, &InterfaceMP.patch(33), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(37), 4, &InterfaceMP.patch(34), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(36), 4, &InterfaceMP.patch(35), 1);
}

