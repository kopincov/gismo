/** @file gsRecipeAssemblerPoissonWithMultiGrid.cpp

    @brief Poisson example using the recipe assembler and multigrid.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/


#include <gismo.h>
#include <gismo_dev.h>
#include <gsMultiGrid/gsGridHierarchy.h>
#include <gsSolver/gsSolverUtils.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>
#include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>
#include <gsPde/gsPoissonPde.h>


using namespace gismo;

void gsIterativeSolverInfo(const gsIterativeSolver<real_t> &method, std::string methodName, real_t time)
{
    gsInfo << methodName +": System size         : " << method.size() << "\n";
    gsInfo << methodName +": Tolerance           : " << method.tolerance() << "\n";
    gsInfo << methodName +": Residual error      : " << method.error() << "\n";
    gsInfo << methodName +": Number of iterations: " << method.iterations() << "\n";
    gsInfo << methodName +": Time to solve:      : " << time << "\n";
}

int main(int argc, char *argv[])
{
    index_t numRefine = 3;
    index_t numDegree = 0;
    if (argc >= 2)
        numRefine = atoi(argv[1]);
    if (argc >= 3)
        numDegree = atoi(argv[2]);


    gsInfo<< std::setprecision(6);

    gsFunctionExpr<> solSin ("sin(2*pi*x) * sin(4*pi*y)",2);
    gsFunctionExpr<> fSin("20*pi**2*sin(2*pi*x)*sin(4*pi*y)",2);

    gsMultiPatch<> * geo = NULL;
    geo = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );


    gsMultiBasis<> geoBases(*geo);
    geoBases.degreeElevate(1,0);

    for (int i = 0; i < numDegree; ++i)
        geoBases.degreeElevate();
    for (int i = 0; i < numRefine; ++i)
        geoBases.uniformRefine();


    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &solSin);//Annulus: small arch lenght
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &solSin);//Annulus: Large arch lenght
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &solSin);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &solSin);


    //PDE info container
    gsPoissonPde<real_t> pdeInfo(*geo,bcInfo,fSin, &solSin);

    gsRecipeAssemblerPoisson assembler(pdeInfo);
    assembler.setDirichletStrategy(dirichlet::elimination);
    std::vector<gsPhysicalSpace*> phySpace;
    phySpace.push_back(new gsPhysicalSpaceScalar(geoBases,*geo,INVERSE_COMPOSITION));
    assembler.setSpace(phySpace);

    assembler.assemble();
    gsMatrix<>  eli;
    gsMatrix<>  rhs = assembler.getSystemRhs();

    if (true)
    {
        gsSparseSolver<>::CGDiagonal solver;
        eli  = solver.compute(assembler.getEliminatedMatrix()).solve(assembler.getEliminatedRhs());
        rhs -= assembler.getRhsModMatrix()*eli;
    }

    gsSparseMatrix<>sysSparse = assembler.getSystemMatrix();
    gsInfo << "Matrix size: " << sysSparse.rows() << " X " << sysSparse.cols()<<"\n";


    //Setting up Multigrid preconditioner
    gsInfo << "Setting up multigrid" << "\n";
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;
    gsGridHierarchy<>::buildByCoarsening(gsMultiBasis<>(geoBases.basis(0)), bcInfo, gsAssembler<>::defaultOptions(), 100, 60)
        .moveTransferMatricesTo(transferMatrices)
        .clear();


    //Set the smoother
    gsMultiGridOp<> mg(sysSparse,transferMatrices);
   
    for (int i = 1; i < mg.numLevels(); ++i)
        mg.setSmoother(i, makeGaussSeidelOp(mg.matrix(i)));

    mg.setNumPreSmooth(2);
    mg.setNumPostSmooth(2);
    //End inite mg

    //Initialize the CG solver
    gsMultiGridOp<>::Ptr mg_ptr = memory::make_shared_not_owned(&mg);
    gsConjugateGradient<> CGSolver(sysSparse,mg_ptr);
    CGSolver.setMaxIterations(2000);
    CGSolver.setTolerance(1e-12);

    //Set the initial guess to zero
    gsMatrix<> x0;
    x0.setZero(sysSparse.cols(),1);

    //Solve system with given preconditioner (solution is stored in x0)
    gsInfo << "\nCG: Started solving..."  << "\n";
    gsStopwatch clock;
    CGSolver.solve(rhs, x0);
    gsIterativeSolverInfo(CGSolver, "CG", clock.stop());
    gsMatrix<> solVector= x0;


    // Compute the error
    std::vector<gsMatrix<real_t> > coefs;
    for (size_t s=0; s< phySpace.size(); ++s)
        coefs.push_back(assembler.reconstructSolution(s,solVector,eli));

    std::vector<gsFunction<real_t>*> svec(1,
    const_cast<gsFunctionExpr<real_t>*>(&solSin));
    gsRecipeDistance dist(*geo,phySpace,coefs,svec);
    gsMatrix<> norms(2,2);
    norms<<0,2,1,2;
    dist.setSeminorms(0,norms);
    dist.assemble();

    gsMatrix<> error = dist.getDistance(0).rowwise().sum();
    real_t errorL2 = math::sqrt(error(0,0));
    real_t errorH1 = math::sqrt(error(1,0));

    gsInfo << "The L2 error is: " << errorL2 << "\n";
    gsInfo << "The H1 error is: " << errorH1 << "\n";



    gsInfo << "Test is done: Cleaning up..." << "\n"; //freeAll(m_bconditions);

    delete geo;
    freeAll(phySpace);

    gsInfo << "Test is done: Exiting" << "\n";

    return  0;

}


