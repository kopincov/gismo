/* 	@file gsExprNewtonIteratorTest.cpp
	
   	@brief solves the a nonlinear elliptic equation showing the usage of gsExprNewtonIterator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Beiser

*/

////////////////////////////////////////
// Including gismo header 
////////////////////////////////////////

#include <gismo.h>
#include <gismo_dev.h>
#include <gsPde/gsExprNewtonIterator.h>
#include <gsAssembler/gsL2Projection.h>

using namespace gismo;


int main (int argc, char *argv[])
{

    ////////////////////////////////////////
    // Use Command Line Arguments
    ////////////////////////////////////////

    index_t numRefine = 4;
    index_t numDegreeIncrease = 2;

    // Command line argument parser
    gsCmdLine cmd("This file solves the Poisson equation as basic example of G+Smo");

    // Add command line arguments
    cmd.addInt("r", "Refinement", "Number of global refinements", numRefine);
    cmd.addInt("p", "Degree", "Number of order elevation steps", numDegreeIncrease);

    // Read parameters from command line
    try { cmd.getValues(argc,argv);  } catch (int rv) { return rv; }


	////////////////////////////////////////
	// Create gsMultiPatch-geometry
	////////////////////////////////////////

	// Initialize MultiPatch
	gsMultiPatch<> MultiPatch;

	// Define MulitPatch via NURBS-Creator
	MultiPatch = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);


	////////////////////////////////////////
	// Define gsMultiBasis
	////////////////////////////////////////

	// Construct MultiBasis on MulitPatch
	gsMultiBasis<> MultiBasis(MultiPatch);

    // p-refinement
    MultiBasis.degreeIncrease( numDegreeIncrease-1);

    // h-refinement
    for (int i = 0; i < numRefine; i++){
        MultiBasis.uniformRefine( );
    }


    ////////////////////////////////////////
    // Define BoundaryConditions
    ////////////////////////////////////////

    gsFunctionExpr<> BoundaryCondition_g("0.0", 2); 

    gsBoundaryConditions<> BoundaryConditions;

    for (gsMultiPatch<>::const_biterator bit = MultiPatch.bBegin(); bit != MultiPatch.bEnd(); ++bit)
    {
        BoundaryConditions.addCondition( *bit, condition_type::dirichlet, BoundaryCondition_g);
    }


	////////////////////////////////////////
	// Assembling
	////////////////////////////////////////

	// Initialize and construct assembler
	gsExprAssembler<> A(1,1); // 1-dim test functions, 1-dim trial functions

	// Set options for assembler
    gsOptionList OptionList = gsAssembler<>::defaultOptions(); // for more options see poisson2_example 

    // Set elements for quadrature
    A.setIntegrationElements( MultiBasis );

    // Set geometryMap
    gsExprAssembler<>::geometryMap G = A.getMap( MultiPatch );

    // Set discretization space 
    gsExprAssembler<>::space H = A.getSpace( MultiBasis ); // Define space by MultiBasis
    H.setInterfaceCont(0); // Set continuity over patch boundaries to C^0
    H.addBc( BoundaryConditions.get("Dirichlet") ); // Invoke Dirichlet boundary part into the space

    // Define RightHandSide
    gsFunctionExpr<> RightHandSide_f( "- 2*pi^2*sin(pi*x)*cos(pi*x)^2*sin(pi*y)^3 - 2*pi^2*sin(pi*x)^3*sin(pi*y)*cos(pi*y)^2 + 2*pi^2*(sin(pi*x)^2*sin(pi*x)^2 + 0.1^2)*sin(pi*x)*sin(pi*y)", 2 ); 
    RightHandSide_f = gsFunctionExpr<>( "0.0", 2 );

    gsMultiPatch<> RightHandSide_MultiPatch;
    gsL2Projection<real_t> L2Projection( RightHandSide_f, MultiPatch, MultiBasis, RightHandSide_MultiPatch );
    gsExprAssembler<>::variable RightHandSide_variable = A.getCoeff( RightHandSide_MultiPatch );

    // Define Initial Solution 
    gsFunctionExpr<> initialSolution( "1.0", 2 );

    gsMultiPatch<> initialSolution_MultiPatch;
    L2Projection = gsL2Projection<real_t>( initialSolution, MultiPatch, MultiBasis, initialSolution_MultiPatch );
    gsExprAssembler<>::variable u_n = A.getCoeff( initialSolution_MultiPatch );

    // Define eps constant function
    gsFunctionExpr<> epsFunction( "0.1", 2 );

    gsMultiPatch<> eps_MultiPatch;
    L2Projection = gsL2Projection<real_t>( epsFunction, MultiPatch, MultiBasis, eps_MultiPatch );
    gsExprAssembler<>::variable eps_variable = A.getCoeff( eps_MultiPatch );


    // Define variational formulation for NewtonIterator
    auto exprJacobian = 
            ( 2*u_n.val()*igrad(u_n, G)[0].val() * igrad(H, G)[0] * H.tr() + 2*u_n.val()*igrad(u_n, G)[1].val() * igrad(H, G)[1] * H.tr() ) * meas(G) 
             + (u_n.val()*u_n.val() + eps_variable.val()*eps_variable.val()) * igrad(H, G) * igrad(H, G).tr() * meas(G) //;
             + eps_variable.val() * igrad(H, G) * igrad(H, G).tr() * meas(G);

    auto exprResidual = 
            H * RightHandSide_variable.val() *  meas(G) 
            - (u_n.val()*u_n.val() + eps_variable.val()*eps_variable.val()) * (igrad(u_n, G)[0].val() * igrad(H, G)[0] + igrad(u_n, G)[1].val() * igrad(H, G)[1]) * meas(G);

    typedef decltype(exprJacobian) E1;
    typedef decltype(exprResidual) E2;

    gsExprNewtonIterator<real_t, E1, E2> ExprNewtonIterator(
                                                            MultiPatch,
                                                            A, 
                                                            G,
                                                            H,
                                                            exprJacobian,
                                                            exprResidual,
                                                            initialSolution_MultiPatch
                                                            );

    ExprNewtonIterator.solve();

    gsMultiPatch<> solutionMultiPatch = ExprNewtonIterator.solution();

    ////////////////////////////////////////
    // Post processing
    ////////////////////////////////////////

    gsField<> solutionField ( MultiPatch, solutionMultiPatch );
    gsWriteParaview( solutionField, "solution" );

    // Exact solution
    gsFunctionExpr<> ExactSolution("sin(pi*x)*sin(pi*y)", 2);
    ExactSolution = gsFunctionExpr<>( "1.0", 2 );
    gsExprAssembler<>::variable H_ex = A.getCoeff( ExactSolution, G );

    // Determine L2 error 
    gsExprEvaluator<> ev( A );
    gsExprAssembler<>::variable H_sol_class = A.getCoeff(solutionMultiPatch);
    real_t errorL2 = math::sqrt( ev.integral( (H_sol_class - H_ex).sqNorm() * meas(G) ) );
    

    gsInfo << "DoF :" << MultiBasis.totalSize() << "\n";
    gsInfo << "L2 error: " << errorL2 << "\n";


}
