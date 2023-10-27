/** @file gsElasticityMixedTH_test.cpp

    @brief Example testing 2D/3D (near) incompressible linear and nonlinear elasticity solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger
*/

#ifndef GISMO_EXT_ELASTICITY

int main(int argc, char** argv)
{
    return 1;
}

#else
#include <gismo.h>
#include <gismo_dev.h>
#include <gsIO/gsMatrixIO.h>

#include <gsElasticity/gsElasticityMixedTHAssembler.h>
#include <gsElasticity/gsElasticityMixedTHNewton.h>

using namespace gismo;

;


/*
	Comments:
	- Validated for Cook's membrane with nu=0.5 and nu=0.4999 (2D linear and nonlinear)
    - Validated for 2-patch plate with hole (2D nonlinear)
	- Valdiated for Alenia pipe (3D nonlinear)
*/

int main (int argc, char** argv)
{
	int result = 0;
	
	// Input options
    int numElevate    = 1;
    int numHref       = 2;
	bool DO_NONLINEAR = false;
    bool plot         = false;
	int TEST_CASE = 0;	// 0: Cook's membrane

	// Material paramaters
    real_t E_modulus     = 1.2E6;
    real_t poisson_ratio = 0.5;
	real_t rho           = 1e3;

	// Command line arguments
        gsCmdLine cmd("Example testing 2D/3D elasticity solver");
        gsInfo << "Type \"-h\" to see the description for all arguments.\n\n";
        cmd.addInt   ("t",  "testcase", "Choose a test case " 
                        "(0: Cook\'s membrane, 1: Plate with hole as 2 patches, 2: Alenia pipe)", TEST_CASE);
	cmd.addInt   ("r", "refine", "Number of h-refinement steps", numHref);
	cmd.addInt   ("e", "elevate", "Number of p-refinement steps", numElevate);
	cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", DO_NONLINEAR);
	cmd.addSwitch("p", "plot", "Plot (otherwise no plot)", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

	// Read geometry
	gsMultiPatch<>::uPtr mp;
	std::string problStr;

	if (TEST_CASE == 0)		// Cook's membrane
	{
                mp = gsReadFile<>("surfaces/cooks_membrane.xml" );
		problStr.assign("cooks_membrane"); 
		gsWriteParaview<>( *mp, problStr + "_geo", 1000);
	}
	else if (TEST_CASE == 1)		// Plate with hole as 2 patches
	{
                mp = gsReadFile<>("planar/plate_with_hole_2p2.xml" );
		problStr.assign("plate_with_whole_2p"); 
		gsWriteParaview<>( *mp, problStr + "_geo", 1000);
	}
	else if (TEST_CASE == 2)		// Alenia pipe
	{
                mp = gsReadFile<>("volumes/alenia_pipe_mod.g2" );
		problStr.assign("alenia_pipe"); 
		gsWriteParaview<>( *mp, problStr + "_geo", 1000);
	}
	else
	{	
		gsInfo << "TEST_CASE " << TEST_CASE << " not available - abort!" << "\n";
		return result;
	}
		
	gsInfo << mp->parDim() << " " << mp->geoDim() << " " << mp->coDim() << "\n";

	mp->computeTopology(1e-4);
	mp->checkConsistency();

	gsInfo << "\n" << "Topology:" << "\n";
	gsInfo << "Interfaces (" << mp->nInterfaces() << "):" << "\n";
	for (int ii = 0; ii < mp->nInterfaces(); ii++)
		gsInfo << mp->interfaces().at(ii) << "\n";

	// Basis refinement
	/*
	gsMultiBasis<> bases( *mp );

	while (numElevate-- > 0)
		bases.degreeElevate();
    
    while (numHref-- > 0)
		bases.uniformRefine();
	*/
	/*
	std::vector<unsigned> numHrefP, numElevateP ;
	numElevateP.resize(mp->nPatches(), numElevate);
	numHrefP.resize(mp->nPatches(), numHref);
	*/
	gsMultiBasis<> bases_p(*mp);
		
	if (TEST_CASE == 2)			
		bases_p.patchBases().at(0)->degreeElevate(1,0);// (elevationAmount,direction)

	while (numElevate-- > 0)
		for (index_t i = 0; i < mp->nPatches(); i++)
			bases_p.patchBases().at(i)->degreeElevate();
    
    while (numHref-- > 0)
		for (index_t i = 0; i < mp->nPatches(); i++)
			bases_p.patchBases().at(i)->uniformRefine();

	if (TEST_CASE == 2)			
		bases_p.patchBases().at(0)->component(2).uniformRefine();

	// Taylor-Hood basis for deformation and pressure
	gsMultiBasis<> bases_u(bases_p);
	bases_u.degreeElevate();

	/*
	std::vector< std::vector<unsigned> > numElevateP; //, numHrefP;
	numElevateP.resize(mp->nPatches());
	//numHrefP.resize(mp->nPatches());

	for (index_t i = 0; i < mp->nPatches(); i++)
	{
		numElevateP[i].resize(mp->dim(), numElevate);
		//numHrefP[i].resize(mp->dim(), numHref);
		// Plate:
		numElevateP[i][1] = numElevate+1;
	}

	for (index_t i = 0; i < mp->nPatches(); i++)
		for (size_t di = 0; di < mp->dim(); di++)
			mp->patch(i).degreeElevate(di, numElevateP[i][di]);
    
    for (index_t i = 0; i < mp->nPatches(); i++)
		while (numHrefP[i]-- > 0)
			mp->patch(i).uniformRefine();
	*/

	gsInfo << "\n" << "Bases:" << "\n";
	for (index_t ip = 0; ip < mp->nPatches(); ip++)
	{
		gsInfo << "Basis_u (patch " << ip << "): "<< bases_u.basis(ip) << "\n";
		gsInfo << "Basis_p (patch " << ip << "): "<< bases_p.basis(ip) << "\n";
	}

	// Boundary conditions

	gsBoundaryConditions<> BCs;
	gsVector<> tmpForce;
	//gsConstantFunction<> * neuData = NULL;
	gsFunction<> * neuData = NULL;

	if (TEST_CASE == 0)					// Cook's membrane
	{
		E_modulus     = 240.565e6;
		poisson_ratio = 0.5; //0.4999; // 
		rho           = 1.0e3;
		
		gsVector<> tmp(2);
		tmp << 0, 62.5e5; 
		neuData = new gsConstantFunction<>(tmp,2);
		BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0);
		BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 1);
		BCs.addCondition(boundary::east, condition_type::neumann, neuData);
    
		tmpForce.resize(2);
		tmpForce << 0, 0;
	}
	else if (TEST_CASE == 1)			// Plate with hole as 2 patches
	{
		E_modulus     = 71.72e+9;
		poisson_ratio = 0.5;
		rho           = 2800.0;
		
		BCs.addCondition(0, boundary::west , condition_type::dirichlet, 0, 1);	
		BCs.addCondition(1, boundary::east , condition_type::dirichlet, 0, 0);	
	
		neuData = new gsConstantFunction<>(1.e10,0.,2);		
		BCs.addCondition(0, boundary::north, condition_type::neumann, neuData);
    
		tmpForce.resize(2);
		tmpForce << 0, 0;
	}
	else if (TEST_CASE == 2)			// Alenia pipe
	{
		E_modulus     = 240.565e6;
		poisson_ratio = 0.5; //0.4999; // 
		rho           = 1.0e3;
		
		gsVector<> tmp(3);
		tmp << 1e3, 0, 0; 
		neuData = new gsConstantFunction<>(tmp,2);
		BCs.addCondition(5, condition_type::dirichlet, 0, 0);
		BCs.addCondition(5, condition_type::dirichlet, 0, 1);
		BCs.addCondition(5, condition_type::dirichlet, 0, 2);
		BCs.addCondition(3, condition_type::dirichlet, 0, 0);
		BCs.addCondition(4, condition_type::dirichlet, 0, 0);
		BCs.addCondition(6, condition_type::dirichlet, 0, 0);
		BCs.addCondition(6, condition_type::dirichlet, 0, 1);
		BCs.addCondition(6, condition_type::dirichlet, 0, 2);
		BCs.addCondition(1, condition_type::neumann, neuData);
    
		tmpForce.resize(3);
		tmpForce << 0, 0, 0;
	}

	// Body forces
	gsConstantFunction<> bodyForce(tmpForce,mp->dim());

	gsInfo << "\n" << "Initialize Assembler" << "\n";

	// Construct assembler object
	gsElasticityMixedTHAssembler<real_t> elasticityAssembler( *mp, bases_u, bases_p, E_modulus, poisson_ratio, rho, BCs, bodyForce );
	
	gsMultiPatch<> deformation;
	gsMultiPatch<> pressure;

	if (!DO_NONLINEAR)
	{
		// * * * Solve linear elasticity problem * * *

		gsInfo << "\n" << "Linear mixed solve" << "\n";

		gsInfo << "\n" << "Assemble" << "\n";
		elasticityAssembler.assemble();
	
		saveMarket(elasticityAssembler.matrix(),"matK.mat");
		saveMarket(elasticityAssembler.rhs(),"matb.mat");

		gsInfo << "\n" << "Solve system" << "\n";
	
		//gsSparseSolver<>::CGDiagonal solver( elasticityAssembler.matrix() );
		gsSparseSolver<>::LU  solver( elasticityAssembler.matrix() );
		gsMatrix<> solVector = solver.solve( elasticityAssembler.rhs() );

		gsInfo << "Internal energy: " << solVector.transpose() * elasticityAssembler.matrix() * solVector << "\n";
		gsInfo << "External energy: " << elasticityAssembler.rhs().transpose() * solVector << "\n";
	
		saveMarket(solVector,"matx.mat");

		elasticityAssembler.constructSolution(solVector, deformation, 0);
	    elasticityAssembler.constructSolution(solVector, pressure, 1);

		/* // Solver experiments
		gsElasticityMassAssembler<real_t> elasticityMassAssembler( *mp, bases, rho, BCs );
		elasticityMassAssembler.assemble();
		saveMarket(elasticityMassAssembler.matrix(),"matM.mat");
		gsGaussSeidelOp<gsSparseMatrix<>,gsGaussSeidel::symmetric> gspc(elasticityMassAssembler.matrix());
		gsConjugateGradient<> solver( elasticityAssembler.matrix() );
                solver.setMaxIterations(500);
                solver.setTolerance(1e-6);

		gsMatrix<> solVector(elasticityAssembler.rhs());
		solVector.setZero();
		solver.solve( elasticityAssembler.rhs(), solVector, gspc );
		*/

	}
	else
	{
		// * * * Solve nonlinear elasticity problem * * * 

		gsInfo << "\n" << "Nonlinear solve" << "\n";

		gsElasticityMixedTHNewton<real_t> newtonSolver(elasticityAssembler, deformation, pressure);
		newtonSolver.setMaxIterations(20);
		newtonSolver.setTolerance(1e-6);
		
		gsMatrix<real_t> solVector; 
		//newtonSolver.solVector(solVector);
		newtonSolver.solve();
		newtonSolver.solVector(solVector);

		if ( newtonSolver.converged() )
		{
			gsInfo <<"Converged after "<<newtonSolver.numIterations()
					<<" iterations with tolerance "<<newtonSolver.tolerance() <<".\n";
			deformation = newtonSolver.solution(0);
			pressure = newtonSolver.solution(1);
		}
		else
			gsInfo <<"Newton iteration did not converge.\n";
	}

	// Post-processing

	for (index_t ip = 0; ip < mp->nPatches(); ip++)
	{
		gsInfo <<"Maximum deformation coef: "
			   << deformation.patch(ip).coefs().colwise().maxCoeff() <<".\n";
		gsInfo <<"Minimum deformation coef: "
			   << deformation.patch(ip).coefs().colwise().minCoeff() <<".\n";
		
		gsVector<> eval_x(deformation.parDim());
		gsMatrix<> eval_u(deformation.dim(),1);
		//eval_x << 0.5, 0.5, 1.0;
		eval_x << 1.0, 1.0;			// Cook's membrane
		eval_x << 1.0, 0.0; 		// Plate with hole 2p
		deformation.patch(ip).eval_into(eval_x, eval_u);
		gsInfo << "Eval: " << eval_u <<".\n";
	}

    if (plot)
    {            
        // Write solution to paraview file
        gsInfo << "Plotting in Paraview...\n";
		
		gsField<> solField(*mp, deformation);
        gsWriteParaview<>( solField, problStr + "_def", 1000 );

		gsField<> prexField(*mp, pressure);
        gsWriteParaview<>( prexField, problStr + "_prex", 1000 );
        
        // Run paraview on exit
        //result = system("paraview thin_shell_def.pvd &");
    }


    // Clean up and exit
	delete neuData;

	return result;
}

#endif //GISMO_EXT_ELASTICITY
