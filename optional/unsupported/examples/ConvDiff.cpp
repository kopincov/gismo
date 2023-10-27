// example solving the convection-diffusion equation with gismo

// requires cleaning up and organizing the examples

#include <iostream>

#include <gismo.h>


using namespace gismo;


int main(int argc, char *argv[])
{   
	// Input options
	index_t numElevate  = 0;
	index_t numHref     = 0;
	index_t basisDegree = 0;
	// Flag for SUPG-Stabilization
	stabilizerCDR::method Flag_Stabilization = stabilizerCDR::none;
	//Flag_Stabilization = stabilizerCDR::SUPG;
	bool plot       = false;

	int result = 0;

	gsCmdLine cmd("Testing a multipatch problem.");
	cmd.addInt("r","hRefine", 
		"Number of dyadic h-refinement (bisection) steps to perform before solving",
		numHref);
	cmd.addInt("p","degree", 
		"Degree of the basis functions to use for solving (will elevate or reduce the input)",
		basisDegree);
	cmd.addInt("e","degreeElevation", 
		"Number of degree elevation steps to perform on the Geometry's basis before solving", 
		numElevate);
	cmd.addSwitch("plot", "Plot result in ParaView format", plot);

	try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


	// --- Example 1: with complex geometry and anisotropic diffusion square_with_sq_hole ------------------------------
	/*
    gsFunctionExpr<real_t>  f("0", 2);
    gsFunctionExpr<real_t>  h("1", 2);
    gsFunctionExpr<real_t>  g("if( (y<0.6) and (y>0.4) and (x>0.4) and (x<0.6),1,-1)", 2);

    gsFunctionExpr<real_t>  coeff_A("75.2500","42.8683","42.8683","25.7500", 2);
    gsFunctionExpr<real_t>  coeff_b("0","0", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
	gsBoundaryConditions<> BCs;
	// -- first benchmark

	gsMultiPatch<real_t> mp;
	gsReadFile<>("planar/square_with_sq_hole.xml", mp);
	for (gsMultiPatch<>::const_biterator 
	bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
	{
	BCs.addCondition( *bit, condition_type::dirichlet, &g);
	}
	//*/

	// --- Example 2: with complex geometry and anisotropic diffusion rectangle_with_2_sq_holes ------------------------------
	/*
    gsFunctionExpr<real_t>  f("0", 2);
    gsFunctionExpr<real_t>  h("1", 2);
    gsFunctionExpr<real_t>  g("if( (y<0.6) and (y>0.4) and (x>0.4) and (x<0.6),1,-1)", 2);

    gsFunctionExpr<real_t>  coeff_A("75.2500","42.8683","42.8683","25.7500", 2);
    gsFunctionExpr<real_t>  coeff_b("0","0", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
	gsBoundaryConditions<> BCs;
	// -- second benchmark
	gsMultiPatch<real_t> mp;
    gsReadFile<>("planar/rectangle_with_2_sq_holes.xml", mp);
	BCs.addCondition(0,3,condition_type::neumann, &f);
	BCs.addCondition(4,3,condition_type::neumann, &f);
	BCs.addCondition(5,2,condition_type::neumann, &f);
	BCs.addCondition(6,4,condition_type::neumann, &f);
	BCs.addCondition(2,4,condition_type::neumann, &f);
	BCs.addCondition(3,1,condition_type::neumann, &f);

	BCs.addCondition(0,4,condition_type::dirichlet, &f);
	BCs.addCondition(1,1,condition_type::dirichlet, &f);
	BCs.addCondition(2,3,condition_type::dirichlet, &f);
	BCs.addCondition(3,2,condition_type::dirichlet, &f);

	BCs.addCondition(4,4,condition_type::dirichlet, &h);
	BCs.addCondition(5,1,condition_type::dirichlet, &h);
	BCs.addCondition(6,3,condition_type::dirichlet, &h);
	BCs.addCondition(7,2,condition_type::dirichlet, &h);
	//*/

	// --- Example 3: Unit square, advection-diffusion, internal layer skew to mesh --------------------------------------------
	/*
    gsFunctionExpr<real_t>  f("0", 2);
    gsFunctionExpr<real_t>  g("if( y<=0.2-0.2*x,1,0)", 2);

    gsFunctionExpr<real_t>  coeff_A("0.00001","0","0","0.00001", 2);
    gsFunctionExpr<real_t>  coeff_b("sqrt(2)","sqrt(2)", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
	gsMultiPatch<real_t> mp( *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0) );
	gsBoundaryConditions<> BCs;
	for (gsMultiPatch<>::const_biterator 
	bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
	{
	BCs.addCondition( *bit, condition_type::dirichlet, &g);
	}
	//*/

	// --- Example 4:Unit square, diffusion only, test of convergence ----------------------------------------------------
	/*
    gsFunctionExpr<real_t>  f("2*pi^2*sin(pi*x)*sin(pi*y)", 2);
    gsFunctionExpr<real_t>  g("sin(pi*x) * sin(pi*y)", 2);
	gsMultiPatch<real_t> mp( *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0) );
    gsFunctionExpr<real_t>  coeff_A("1.0","0","0","1.0", 2);
    gsFunctionExpr<real_t>  coeff_b("0","0", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
	gsBoundaryConditions<> BCs;
	for (gsMultiPatch<>::const_biterator 
	bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
	{
	BCs.addCondition( *bit, condition_type::dirichlet, &g);
	}
	//*/

	// --- Example 5:Unit square, advection-diffusion, test of convergence ----------------------------------------------------
	///*
    gsFunctionExpr<real_t>  f("-1*(60*x^2+24*x+6)+0.2*(20*x^3+12*x^2+6*x+4)", 2);
    gsFunctionExpr<real_t>  g("5*x^4+4*x^3+3*x^2+2*x+y", 2);
    gsMultiPatch<real_t> mp( *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0) );
    gsFunctionExpr<real_t>  coeff_A("1.0","0","0","1.0", 2);
    gsFunctionExpr<real_t>  coeff_b("0.2","0.4", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
	gsBoundaryConditions<> BCs;
	for (gsMultiPatch<>::const_biterator 
		bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
	{
		BCs.addCondition( *bit, condition_type::dirichlet, &g);
	}
	//*/

	// Set up and elevate discretization bases
	gsMultiBasis<> bases(mp);

	if (basisDegree)
		bases.setDegree(basisDegree);
	else if (numElevate)
		bases.degreeElevate(numElevate);
	//               bases.degreeReduce(1);	

	// Run tests
	gsMatrix<> testsL2(numHref+1,3);
	testsL2.setZero();

	gsMatrix<> testsH1(numHref+1,3);
	testsH1.setZero();

	//gsSparseSolver<>::CGDiagonal solver; // identity preconditioner
	//gsSparseSolver<>::CGDiagonal solver;
	gsSparseSolver<>::LU solver;	
	//gsSparseSolver<>::CGDiagonal solver;
	//gsSparseSolver<>::BiCGSTABILUT solver( galerkin.matrix() );
	//gsSparseSolver<>::BiCGSTABILUT solver( galerkin.fullMatrix() );

	int i = 0;
	do
	{
		// Setup the assemblers
		gsCDRAssembler<real_t> galerkin(mp,bases,BCs,f,coeff_A,coeff_b,coeff_c,
			dirichlet::elimination, iFace::glue, Flag_Stabilization);

		gsInfo<<"Discretization Space for patch 0: \n"<< bases[0] << "\n";
		gsInfo<<"System  size (elim. BCs)  : "<< galerkin.numDofs() << "\n";
//		gsInfo<<"Coupled size (elim. BCs)  : "<< galerkin.dofMapper().coupledSize() << "\n";   
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
			/*
			gsInfo << "residual error: " << solver.error() << "\n";
			gsInfo << "    iterations: " << solver.iterations() << "\n";
			gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
			//*/
		}
		gsField<> sol = galerkin.constructSolution(solVector);
        
		// Collect data
		testsL2(i,0)  =
			testsH1(i,0)= bases.totalSize();

		testsL2(i,1)= sol.distanceL2(g) ;

		testsH1(i,1) = sol.distanceH1(g) ;


		if (i > 0)
		{
			testsL2(i,2)= testsL2(i-1,1) / testsL2(i,1);

			testsH1(i,2)= testsH1(i-1,1) / testsH1(i,1);
		}

		gsInfo << "Total time (elim. BCs): " << assTime+solvTime   << " s" << "\n";
		gsInfo << "---------------------------------------\n";    
		gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING      \n";
		gsInfo << "    Dofs   |  L2 error  | err. ratio   \n" << testsL2.row(i)  << "\n";  
		gsInfo << "           |  H1 error  | err. ratio   \n" << testsH1.row(i)  << "\n";  

		if (plot && i== numHref)
		{
			// Write approximate and exact solution to paraview files
			gsInfo<<"Plotting in Paraview...\n";
			gsWriteParaview<>(sol, "aaa__CD", 1000);

			// Run paraview
			result = system("paraview poisson_problem.vts &");
		}

		bases.uniformRefine();
	} 
	while ( i++ < numHref );

	for(i = 1; i<= numHref; ++i)
	{   // Compute convergence rates
		testsL2(i,2)= math::log(testsL2(i,2))/math::log(2.0);

		testsH1(i,2)= math::log(testsH1(i,2))/math::log(2.0);
	}

	gsInfo << "Summary:\n\n";
	gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING           \n";
	gsInfo << "    Dofs   |  L2 error  | conv. rate   \n" << testsL2  << "\n";


	gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING           \n";
	gsInfo << "    Dofs   |  H1 error  | conv. rate   \n" << testsH1  << "\n";

	return result;
}

