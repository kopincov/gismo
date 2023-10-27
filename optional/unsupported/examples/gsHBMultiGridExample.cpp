/** @file gsHBMultiGridExample.cpp

    @brief Provides test examples for multigrid algorithms for THB-splines

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither

    Date: 2014
*/

#include <iostream>
#include <ctime>
#include <cstdlib>

#include <gismo.h>
#include <gismo_dev.h>
#include <gsUtils/gsNorms.h>
#include <gsSolver/gsTimedOp.h>


using namespace gismo;

namespace Smoother {
    enum type {
        Richardson,
        Jacobi,
        GaussSeidel
    };
}

void testTransferMatrix(const gsBasis<real_t>& basis0, const gsBasis<real_t>& basis1, const gsSparseMatrix<real_t>& transfer)
{
    // create random coefficients for coarse basis
    gsMatrix<> coefs0;
    coefs0.setRandom( basis0.size(), 1 );

    // apply transfer matrix
    gsMatrix<> coefs1 = transfer * coefs0;

    // construct two geometries
    gsGeometry<real_t>::uPtr geo0 = basis0.makeGeometry( coefs0 );
    gsGeometry<real_t>::uPtr geo1 = basis1.makeGeometry( coefs1 );

    gsMatrix<> paramRange = geo0->parameterRange();

    const real_t dist = computeMaximumDistance<real_t>(*geo0, *geo1, paramRange.col(0), paramRange.col(1));

    if (dist > 1e-14)
    {
        gsInfo << "\n";
        gsInfo << "   ***** WARNING! Transfer matrix test has failed! *****" << "\n";
        gsInfo << "   ***** Maximum distance: " << dist << "\n" << "\n";
    }
}

/// Provides the name of the specified smoother
std::string smootherName( Smoother::type smoother )
{
    switch (smoother)
    {
        case Smoother::Richardson: return "a Richardson smoother";
        case Smoother::Jacobi: return "a Jacobi smoother";
        case Smoother::GaussSeidel: return "a Gauss-Seidel smoother";
    }
    return "an unknwon smoother";
}

int main(int argc, char *argv[])
{
    index_t degree = 2;
    index_t numLevels = 1;
    Smoother::type smoother;
    index_t numPreSmooth = 1;
    index_t numPostSmooth = 1;
    index_t cycles = 1;
    bool useCG = false;
    bool useFMG = false;
    bool useCascadic = false;
    bool monitorL2 = false;
    bool writeLog = false;
    real_t tol = std::pow(10.0, - REAL_DIG * 0.5);
    real_t damping = -1.0;
    bool plot = false;
    std::string smooth("gs");
    //bool spectral = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addInt("p", "degree", "Degree of the B-spline discretization space", degree);
    cmd.addInt("l", "levels", "Number of levels to use for multigrid iteration", numLevels);
    cmd.addInt("", "presmooth", "Number of pre-smoothing steps", numPreSmooth);
    cmd.addInt("", "postsmooth", "Number of post-smoothing steps", numPostSmooth);
    cmd.addString("s", "smoother", "Smoothing method", smooth);
    cmd.addReal("", "damping", "Damping factor for the smoother", damping);
    cmd.addInt("c", "cycles", "Number of multi-grid cycles", cycles);
    cmd.addSwitch("log", "Write results to log file", writeLog);
    cmd.addSwitch("cg", "Use CG iteration", useCG);
    cmd.addSwitch("fmg", "Use full multigrid cycle", useFMG);
    cmd.addSwitch("cascadic", "Use cascadic multigrid", useCascadic);
    cmd.addReal("", "tol", "Tolerance for multigrid solver stopping criterion", tol);
    //cmd.addSwitch( "spectral-test", "Perform a numerical spectral test", spectral);
    cmd.addSwitch("monitor-L2", "Monitor the L2 errors over the iteration", monitorL2);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (numLevels < 1)
    {
        std::cerr << "Number of levels must be positive.\n"; return -1;
    }

    if (smooth == "r" || smooth == "richardson")
        smoother = Smoother::Richardson;
    else if (smooth == "j" || smooth == "jacobi")
        smoother = Smoother::Jacobi;
    else if (smooth == "gs" || smooth == "gauss-seidel")
        smoother = Smoother::GaussSeidel;
    else
    {
        gsWarn << "Unknown smoother \"" << smooth << "\".\n"; 
        gsWarn << "Allowed are: richardson (r), jacobi (j), and gauss-seidel (gs).\n";
        return -1;
    }
    if (cycles < 1)
    {
        std::cerr << "Number of cycles must be positive.\n"; return -1;
    }
    if (useCG && useFMG)
    {
        std::cerr << "Cannot use CG for full multigrid.\n"; return -1;
    }
    if (damping < 0)
    {
        // Here, one could add appropriate parameter-choice rules (depending on the degree)
        switch (smoother)
        {
            case Smoother::Richardson:                                     damping = 0.20; break;
            case Smoother::Jacobi:                                         damping = 0.80; break;
            case Smoother::GaussSeidel:                                    break;
        }
    }

    /******************** Define Geometry ********************/
    gsStopwatch time;

    gsGeometry<>* geo = ((gsGeometry<>::uPtr)gsReadFile<>("planar/adapt_geo.xml")).release();
    gsMultiPatch<> mp( *geo );

    //gsFunctionExpr<> f = gsFunctionExpr<>("(2*pi^2 ) * sin(pi*x) * sin(pi*y)",2);
    //gsFunctionExpr<> g = gsFunctionExpr<>("sin(pi*x) * sin(pi*y)",2);

    // example: source field with disk-shaped singularity
    gsFunctionExpr<>  f("if( (x-0.25)^2 + (y-0.6)^2 < 0.2^2, 1, 0 )", 2);
    gsFunctionExpr<>  g("0", 2);

    // example: L-shape with singularity at reentrant corner
    //gsFunctionExpr<>  f("0");
    //gsFunctionExpr<>  g("if( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) +3*pi)/3.0 ) )");

    //gsFunctionExpr<> f("-((100*exp(-sqrt(13+20*x*(-2+5*x)+20*y*(-3+5*y))) * (-1+sqrt(13+20*x*(-2+5*x)+20*y*(-3+5*y)))) / sqrt(13+20*x*(-2+5*x)+20*y*(-3+5*y)))");
    //gsFunctionExpr<> g("1/exp(sqrt((10*x-2)^2+(10*y-3)^2))");

    gsInfo << "Geometry: " << *geo << "\n";
    gsInfo << "Source function: " << f << ".\n" << "\n";
    gsInfo << "Exact solution:  " << g << ".\n" << "\n";

    // set up the boundary value problem
    //gsBVProblem<> bvp( geo, new gsPoissonPde<>(f, geo->geoDim()) ); // -Delta u = f

    const condition_type::type bc_type         = condition_type::dirichlet;
    gsFunction<> * const bc_func   = &g;
    //const boundary::type bc_type         = condition_type::neumann;
    //gsFunction<> * const bc_func   = &zero;

    gsBoundaryConditions<> bc;
    bc.addCondition( boundary::west,  bc_type, bc_func );
    bc.addCondition( boundary::east,  bc_type, bc_func );
    bc.addCondition( boundary::south, bc_type, bc_func );
    bc.addCondition( boundary::north, bc_type, bc_func );


    std::vector<gsBasis<>*> bases;
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;

    for (int i = 0; i < numLevels; ++i)
    {
        std::string fname = "basis2d/adapt_basis_" + util::to_string(i) + ".xml";
        bases.push_back(
            ((gsBasis<>::uPtr)gsReadFile<>(fname)).release()
        );
    }

    numLevels = (int)bases.size();

    transferMatrices.resize( numLevels - 1 );

    // compute HB transfer matrices
    for (int i = 0; i < numLevels - 1; ++i)
    {
        gsHTensorBasis<2,real_t> * bc2 = dynamic_cast<gsHTensorBasis<2,real_t>*>(bases[i]);
        gsHTensorBasis<2,real_t> * bf = dynamic_cast<gsHTensorBasis<2,real_t>*>(bases[i+1]);

        gsSparseMatrix<> t;
        bf->transfer(bc2->getXmatrix(), t);
        transferMatrices[i] = t;

        testTransferMatrix(*bases[i], *bases[i+1], transferMatrices[i]);
    }


    gsInfo << (useFMG ? "Full multigrid" : (useCascadic ? "Cascadic multigrid" : (useCG ? "CG preconditioned by multigrid" : "Multigrid")));
    gsInfo << " with " << numLevels << " levels using " << smootherName(smoother);
    gsInfo << "(" << numPreSmooth << "," << numPostSmooth << ") smoother" << " and " << cycles << "-cycle.\n";
    gsInfo << "Discretization space: dim=" << bases[0]->dim() << " deg=" << bases[0]->degree(0) << "\n";

    
    // construct a vector of dof-mappers and a vector of pointers pointing onto them
    std::vector< gsDofMapper > dofMappers( bases.size() );

    gsOptionList assemblerOptions = gsAssembler<>::defaultOptions();

    for (int i = 0; i < numLevels; ++i)
    {
        gsInfo << "Level " << i << ": " << bases[i]->size() << " dofs\n";

        gsMultiBasis<>(*bases[i]).getMapper((dirichlet::strategy)assemblerOptions.getInt("DirichletStrategy"),
                                       (iFace::strategy)assemblerOptions.getInt("InterfaceStrategy"),bc,dofMappers[i],0); 
    }

    real_t timeStartAsm = time.stop();

    // assemble fine-grid problem
    gsPoissonAssembler<> assembler(mp, gsMultiBasis<>(*bases[numLevels-1]), bc, f);
    assembler.assemble();
    const gsMatrix<>& rhs = assembler.rhs();

    const real_t timeStartMgSetup = time.stop();
    
    // set up the multigrid solver
    //gsMultiGridOp<> mg(bvp, bases, transferMatrices, useFMG || useCascadic);
    
    // incorporate boundary conditions
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatricesWithBC(transferMatrices.size());
    for ( unsigned i=0; i<transferMatrices.size(); ++i )
        gsMultiBasis<>::combineTransferMatrices( std::vector< gsSparseMatrix<real_t, RowMajor> >(1,give(transferMatrices[i])), dofMappers[i], dofMappers[i+1], transferMatricesWithBC[i] );
    transferMatrices.clear();
    
    gsMultiGridOp<> mg( assembler.fullMatrix(), transferMatricesWithBC );
    
    // Determine coarse solve time:
    gsTimedOp<>::Ptr coarseSolver = gsTimedOp<>::make("CoarseSolver",mg.coarseSolver(),false);
    mg.setCoarseSolver(coarseSolver);

    mg.setNumPreSmooth( numPreSmooth );
    mg.setNumPostSmooth( numPostSmooth );
    mg.setNumCycles( cycles );

    for (int i = 0; i < mg.numLevels(); ++i)
    {
         switch( smoother )
         {
            case Smoother::Richardson:                            mg.setSmoother(i, makeRichardsonOp(mg.matrix(i),damping)); break;
            case Smoother::Jacobi:                                mg.setSmoother(i, makeJacobiOp(mg.matrix(i),damping)); break;
            case Smoother::GaussSeidel:                           mg.setSmoother(i, makeGaussSeidelOp(mg.matrix(i))); break;
         }
    }

    gsMatrix<> x;
    x.setRandom( mg.nDofs(), 1 );

    const real_t timeStartSolve = time.stop();

    const real_t resNorm0 = (rhs - mg.matrix() * x).norm();
    gsInfo << "Residual norm:     " << resNorm0 << "\n";
    real_t resNorm(0), oldResNorm = resNorm0;

    gsField<> sol;

    real_t l2Err, oldL2Err(0);
    if (monitorL2)
    {
        sol = assembler.constructSolution(x);
        oldL2Err = computeL2Distance(sol, g, false, 3*mg.nDofs());
        gsInfo << "L2 error: " << oldL2Err << "\n";
    }

    int numIter = 0;
    real_t minReduction = 1e6;

    gsConjugateGradient<> cg( mg.underlyingOp(), memory::make_shared_not_owned(&mg) );
    cg.setTolerance(tol);

    if (useCG)
        cg.initIteration( rhs, x );

    do
    {
        if (useCG)
            cg.step(x);
        else if (numLevels == 1)
            mg.smoothingStep(rhs, x);
        else
            mg.step(rhs, x);

        resNorm = (rhs - mg.matrix() * x).norm();
        gsInfo << "Residual norm:     " << std::left << std::setw(15) << resNorm << "          reduction:  1 / " << std::setprecision(3) << (oldResNorm/resNorm) << std::setprecision(6) << "\n";
        minReduction = math::min(minReduction, oldResNorm/resNorm);
        oldResNorm = resNorm;

        if (monitorL2)
        {
            sol = assembler.constructSolution(x);
            l2Err = computeL2Distance(sol, g, false, 3*mg.nDofs());
            gsInfo << "                                                                   |  L2 error:     "
                << std::left << std::setw(15) << l2Err << "          reduction:  1 / " << std::setprecision(3) << (oldL2Err/l2Err) << std::setprecision(6) << "\n";
            oldL2Err = l2Err;
        }

        ++numIter;
    } while (/*resNorm > tol && */ resNorm / resNorm0 > tol);
    sol = assembler.constructSolution(x);

    const real_t timeTotal = time.stop();
    const real_t timeAssembling = timeStartMgSetup - timeStartAsm;
    const real_t timeSetup = timeStartSolve - timeStartMgSetup;
    const real_t timeSolve = timeTotal - timeStartSolve;

    if (!useFMG && !useCascadic) {
        gsInfo << "Converged in " << numIter << " iterations.\n";
        gsInfo << "Average convergence factor:  1 / " << std::setprecision(3) << math::pow(resNorm0 / resNorm, 1.0 / numIter) << std::setprecision(6) << "\n";
        gsInfo << "Worst   convergence factor:  1 / " << std::setprecision(3) << minReduction << std::setprecision(6) << "\n";
        gsInfo << "\n";
    }
    gsInfo << "Assembling time: "; formatTime(gsInfo, timeAssembling);        gsInfo << "\n";
    gsInfo << "MG time:         "; formatTime(gsInfo, timeSetup+timeSolve);   gsInfo << "\n";
    gsInfo << " setup:          "; formatTime(gsInfo, timeSetup);             gsInfo << "\n";
    gsInfo << " solving:        "; formatTime(gsInfo, timeSolve);
    if (!useFMG && !useCascadic) {
        gsInfo << "         (avg. "; formatTime(gsInfo, timeSolve/numIter);   gsInfo << " per iteration)" << "\n";
        gsInfo << "  coarse solver: "; formatTime(gsInfo, coarseSolver->getTime());
    }
    gsInfo << "\n";
    gsInfo << "Total time:      "; formatTime(gsInfo, timeTotal);             gsInfo << "\n";
    gsInfo << "\n";

    // Compute L2 error
    const real_t error = computeL2Distance(sol, g, false, 3*mg.nDofs());
    gsInfo << "L2 error: " << error << "\n";

    if (writeLog)
    {
        std::fstream log("out.txt", std::fstream::out | std::fstream::app);
        log << degree << "\t" << damping << "\t" << numPreSmooth << "\t" << numPostSmooth << "\t" << numIter << "\n";
    }

    if (plot)
    {
        // Plot solution in Paraview
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(sol, "multigrid", 1000);

        // Plot exact solution in Paraview
        gsField<> exact( mp, g, false );
        gsWriteParaview<>( exact, "multigrid_exact", 1000);
    }

    delete geo;
    freeAll( bases );

    return 0;
};
