/** @file gsMultiGridPreconditionerBiharmonicExample.cpp

    @brief Provides test examples for multigrid algorithms, where the multigrid method
    on the parameter domain is used as a preconditioner for the problem on the physical
    domain.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, J. Sogn, S. Takacs
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsRankOneAssembler.h>
#include <gsSolver/gsKronecker.h>
#include <gsSolver/gsTimedOp.h>
#include <gsSolver/gsSolverUtils.h>
#include <gsMultiGrid/gsMassSmoother.h>

#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>
#include <gsRecipeAssembler/gsRecipeAssemblerBiharmonicSogn.h>
#include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>

#include <gsPde/gsBiharmonicPde.h>

#include <gsIO/gsCmdLineWithEnumSupport.h>

namespace Smoother {
    enum type {
        Richardson,
        Jacobi,
        GaussSeidel,
        MassRichardson,
        MassRichardsonBoundaryCorrection,
        MassRichardsonSubspaceCorrection,
        MassRichardsonSubspaceCorrectionBiharmonic,
        MassRichardsonSubspaceCorrectionBiharmonicRank1,
        MassRichardsonSubspaceCorrectionFullBiharmonic
    };
}

using namespace std;
using namespace gismo;


gsMatrix<> solve(const gsSparseMatrix<> &sys, const gsMatrix<> &rhs)
{
    gsSparseSolver<>::QR  solver;
    solver.analyzePattern( sys );
    solver.factorize     ( sys );
    return solver.solve( rhs );
}




index_t maxMRSLevels( index_t p, index_t r )
{
    if(p<3)
        return r+1;
    else if(p<5)
        return r;
    else if(p<9)
        return r-1;
    else
        return r-2;
}

int main(int argc, char *argv[])
{
    index_t dim = 2;
    index_t numRefine = 3;
    index_t mult = 1;
    index_t degree = 3;
    index_t numLevels = 0;
    Smoother::type smoother = Smoother::GaussSeidel;
    string smoother_name;
    index_t numPreSmooth = 1;
    index_t numPostSmooth = 1;
    index_t cycles = 1;
    real_t alpha = 0;
    bool writeLog = false;
    real_t tol = 1e-8;
    real_t damping = -1.0;
    bool compEigs = false;

    gsCmdLineWithEnumSupport cmd("Solves a PDE with an isogeometric discretization using a conjugate gradient method, preconditioned with multigrid.");
    cmd.addInt("d", "dimension",
               "Geometric dimension of the problem (2 or 3)", dim);
    cmd.addInt("r", "uniformRefine",
                "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("m", "multiplicity",
               "Multiplicity of knots to insert when refining", mult);
    cmd.addInt("p", "degree",
               "Degree of the B-spline discretization space", degree);
    cmd.addInt("l", "levels",
               "Number of levels to use for multigrid iteration", numLevels);
    cmd.addInt("", "presmooth",
               "Number of pre-smoothing steps", numPreSmooth);
    cmd.addInt("", "postsmooth",
               "Number of post-smoothing steps", numPostSmooth);
    cmd.addEnum("s", "smoother","Smoothing method", smoother)
        .add(Smoother::Richardson,                                     "r",      "richardson"                                     )
        .add(Smoother::Jacobi,                                         "j",      "jacobi"                                         )
        .add(Smoother::GaussSeidel,                                    "gs",     "gauss-seidel"                                   )
        .add(Smoother::MassRichardson,                                 "mr",     "mass-richardson"                                )
        .add(Smoother::MassRichardsonBoundaryCorrection,               "mrb",    "mass-richardson-boundary-correction"            )
        .add(Smoother::MassRichardsonSubspaceCorrection,               "mrs",    "mass-richardson-subspace-correction"            )
        .add(Smoother::MassRichardsonSubspaceCorrectionBiharmonic,     "mrsb",   "mass-richardson-subspace-correction-biharmonic" )
        .add(Smoother::MassRichardsonSubspaceCorrectionBiharmonicRank1,     "mrsbR1",   "mass-richardson-subspace-correction-biharmonic-rank-1-approximation" )
        .add(Smoother::MassRichardsonSubspaceCorrectionFullBiharmonic, "mrsbf",  "mass-richardson-subspace-correction-biharmonic" )
        .writeDescOfChosenOptTo(smoother_name);
    cmd.addReal("", "damping",
                "Damping factor for the smoother", damping);
    cmd.addInt("c", "cycles",
               "Number of multi-grid cycles", cycles);
    cmd.addReal("a", "alpha", "alpha in \"-div(A grad(u)) + alpha u = f\"", alpha);
    cmd.addSwitch("log", "Write results to log file", writeLog);
    cmd.addSwitch("eigs", "Compute eigenvalues of the preconditioned system.", compEigs);
    cmd.addReal("", "tol", "Tolerance for multigrid solver stopping criterion", tol);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (!(dim == 2 || dim == 3))
    {
        cerr << "Only 2- and 3-dimensional problems implemented." << endl;
        return -1;
    }

    /******************** Setup the reference problem ********************/
    
    gsStopwatch time;

    gsGeometry<>::uPtr geo;
    switch (dim) {
        case 2: geo = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(degree)); break;
        case 3: geo = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(degree)); break;
    }
    cout << *geo << endl;

    gsFunctionExpr<> f, g;

    switch (dim)
    {
        case 1:
            f = gsFunctionExpr<>("(pi^4 ) * sin(pi*x)",1);
            g = gsFunctionExpr<>("sin(pi*x)",1);
            break;
        case 2:
            f = gsFunctionExpr<>("(4*pi^4 ) * sin(pi*x) * sin(pi*y)",2);
            g = gsFunctionExpr<>("sin(pi*x) * sin(pi*y)",2);
            break;
        case 3:
            f = gsFunctionExpr<>("(9*pi^4 ) * sin(pi*x) * sin(pi*y) * sin(pi*z)",3);
            g = gsFunctionExpr<>("sin(pi*x) * sin(pi*y) * sin(pi*z)",3);
            break;
        default:
            cerr << "Invalid geometry dimension.\n";
            return -1;
    }

    if (numRefine < 1)
    {
        cerr << "Number of refinements must be positive.\n"; return -1;
    }
    if (mult < 1)
    {
        cerr << "Multiplicity must be positive.\n"; return -1;
    }
    if (numLevels == 0)
    {
        numLevels = numRefine + 1;
        if( smoother == Smoother::MassRichardsonSubspaceCorrection
            || smoother == Smoother::MassRichardsonSubspaceCorrectionBiharmonic
            || smoother == Smoother::MassRichardsonSubspaceCorrectionBiharmonicRank1
            || smoother == Smoother::MassRichardsonSubspaceCorrectionFullBiharmonic

        )
            numLevels = maxMRSLevels(degree,numRefine);
        gsInfo << "The number of levels was chosen to be " << numLevels << ".\n";
    }
    if (numLevels < 1)
    {
        cerr << "Number of levels must be positive.\n"; return -1;
    }
    if (numRefine - numLevels + 1 < 0)
    {
        cerr << "Not enough refinements for the desired number of levels.\n"; return -1;
    }
    
    
    if (cycles < 1)
    {
        cerr << "Number of cycles must be positive.\n"; return -1;
    }
    if (damping < 0)
    {
        // Here, one could add appropriate parameter-choice rules (depending on the degree)
        switch (smoother)
        {
            case Smoother::Richardson:                                      damping = 0.80*0.80; break; //JS2: Don't know if this is a good daping
            case Smoother::Jacobi:                                          damping = 0.80*0.80; break;
            case Smoother::GaussSeidel:                                     break;
            case Smoother::MassRichardson:                                  damping = 0.25*0.25 / ( (degree+1.)*(degree+1.) ); break;
            case Smoother::MassRichardsonBoundaryCorrection:                damping = 0.09*0.09; break;
            case Smoother::MassRichardsonSubspaceCorrection:                damping = 0.09*0.09; break;
            case Smoother::MassRichardsonSubspaceCorrectionBiharmonic:      damping = 0.09*0.09; break;
            case Smoother::MassRichardsonSubspaceCorrectionBiharmonicRank1: damping = 0.09*0.09; break;
            case Smoother::MassRichardsonSubspaceCorrectionFullBiharmonic:  damping = 0.09*0.09; break;
        }
    }

    gsInfo << "Source function: " << f << ".\n" << "\n";
    gsInfo << "Exact solution:  " << g << ".\n" << "\n";

    // set up boundary conditions

    gsConstantFunction<> zero(0.0, dim);
    gsConstantFunction<> one (1.0, dim);

    gsBoundaryConditions<> bc;
    gsBoundaryConditions<> bc2;

    
    bc.addCondition( boundary::west,  condition_type::dirichlet, &zero );
    bc.addCondition( boundary::east,  condition_type::dirichlet, &zero ); 
    bc.addCondition( boundary::south, condition_type::dirichlet, &zero );
    bc.addCondition( boundary::north, condition_type::dirichlet, &zero );
    bc2.addCondition( boundary::west,  condition_type::neumann, &zero );
    bc2.addCondition( boundary::east,  condition_type::neumann, &zero );
    bc2.addCondition( boundary::south, condition_type::neumann, &zero );
    bc2.addCondition( boundary::north, condition_type::neumann, &zero );

    if (dim >= 3)
    {
        bc.addCondition( boundary::front,  condition_type::dirichlet, &zero );
        bc.addCondition( boundary::back,   condition_type::dirichlet, &zero );
        bc2.addCondition( boundary::front, condition_type::neumann, &zero );
        bc2.addCondition( boundary::back,  condition_type::neumann, &zero );
    }

    gsBasis<> * tbasis = geo->basis().clone().release();
    vector< gsMultiBasis<> > bases;
    vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;
    
    // Define coarse discretization space by refining the basis of the geometry
    for (int i = 0; i < numRefine - numLevels + 1; ++i)
        tbasis->uniformRefine();       // refine until coarsest level

    if (mult>1)
    {
        gsInfo << "Reduce continuity by " << mult-1 << std::endl;
        tbasis->reduceContinuity(mult-1);
    }
    
    gsInfo << "Coarsest basis: " << *tbasis << std::endl;    

    // set up the hierarchy of spaces and transfer matrices between them
    const index_t refineKnots = 1;
    gsGridHierarchy<>::buildByRefinement(gsMultiBasis<>(*tbasis), bc, gsAssembler<>::defaultOptions(), numLevels, refineKnots, mult)
        .moveMultiBasesTo(bases)
        .moveTransferMatricesTo(transferMatrices)
        .clear();

    gsInfo << "CG preconditioned by multigrid";
    gsInfo << " with " << numLevels << " levels using " << smoother_name;
    gsInfo << "(" << numPreSmooth << "," << numPostSmooth << ") smoother" << " and " << cycles << "-cycle.\n";
    gsInfo << "Coarse discretization space: dim=" << tbasis->dim() << " deg=" << tbasis->degree(0) << " dofs=" << tbasis->size() << "\n";



    /******************* Setup the "true" problem ******************/
    
    const double timeAssembling1 = time.stop(); time.restart();
    
    real_t a = 4;
    gsGeometry<>::Ptr geoTrue;
    gsMultiPatch<>::Ptr geoTrueMp;
    switch (dim)
    {
        //case 2: geoTrue = gsNurbsCreator<>::BSplineSquareDeg(1); break;
        //case 2: geoTrue = gsNurbsCreator<>::NurbsQuarterAnnulus(0.5, 1.0); break;
        case 2: geoTrue = gsNurbsCreator<>::BSplineRectangle(0,0,1,a) ; break;
        //case 2: geoTrue = gsNurbsCreator<>::BSplineFatQuarterAnnulus(0.5, 1.0); break;
        //case 3: geoTrue = gsNurbsCreator<>::NurbsCube(); break;
        case 3: {
             //Read 3D geometry
            std::string input("volumes/twistedFlatQuarterAnnulus.xml");
            gsFileData<> fileData(input);
            if (!fileData.has< gsMultiPatch<> >()) { cerr << "No multipatch object found.\n"; return 1; }
            geoTrueMp = fileData.getFirst< gsMultiPatch<> >();
            GISMO_ASSERT( geoTrueMp->nPatches() == 1, "Multigrid only works for single-patch domains so far." );
            geoTrue = memory::make_shared_not_owned(const_cast<gsGeometry<>*>(&((*geoTrueMp)[0])));
            break;
        }
    }

    // use the finest basis from the reference domain for the true problem

    gsInfo << "Setup gsRecipeAssembler for true problem... " << flush;
    gsMultiPatch<> mpTrue(*geoTrue);
    
    gsBiharmonicPde<real_t> pdeBiharmonic(mpTrue,bc,bc2,f);

    gsRecipeAssemblerBiharmonicSogn assemblerTrue(pdeBiharmonic);
    assemblerTrue.setDirichletStrategy(dirichlet::elimination);
    assemblerTrue.setZeroAverage(false);
    vector<gsPhysicalSpace*> phySpace;
    phySpace.push_back(new gsPhysicalSpaceScalar(bases.back(),mpTrue,INVERSE_COMPOSITION));

    assemblerTrue.setSpace(phySpace);
    
    gsInfo << "done." << endl;

    gsInfo << "Assembling true system matrix... " << flush;
    assemblerTrue.assemble();
    gsMatrix<> rhs = assemblerTrue.getSystemRhs();

    gsMatrix<> eli  = solve(assemblerTrue.getEliminatedMatrix(),assemblerTrue.getEliminatedRhs());
    rhs -= assemblerTrue.getRhsModMatrix()*eli;


    gsSparseMatrix<real_t> BfineTrue = assemblerTrue.getSystemMatrix();
    gsInfo << "done, " << BfineTrue.rows() << " dofs." << "\n";

    gsInfo << "Setup assembler... " << flush;

    gsSparseMatrix<> Bfine;
    const gsBasis<>& bs = bases.back()[0];

    if (smoother == Smoother::MassRichardsonSubspaceCorrectionBiharmonicRank1)
    {
        gsInfo << "\n Using Rank 1 B in multigrid  \n";
        assembleRankOneSimpleBiharmonic(* geoTrue, bs, bc, Bfine);
    }
    else if (smoother == Smoother::GaussSeidel)
    {
        gsInfo << "\n Using true B in multigrid \n";
        Bfine = BfineTrue;
    }
    else
    {
        gsInfo << "\n Using parametric B in multigrid \n";
        gsInfo << "Assembling mass and stiffness matrix in 1D... " << flush;

        gsSparseMatrix<> Bx, By, Bz, Mx, My, Mz;


        //assemble1DMass(*geoTrue, qp, 0, bs.component(0), Mx);
        //assemble1D2ndDer(*geoTrue, qp, 0, bs.component(0), Bx);
        assembleParameter2ndDer1D(bs.component(0), Bx);
        assembleParameterMass(bs.component(0), Mx);
        handleDirichletConditions(Mx, bc, boundary::west, boundary::east);
        handleDirichletConditions(Bx, bc, boundary::west, boundary::east);

        //assemble1DMass(*geoTrue, qp, 1, bs.component(1), My);
        //assemble1D2ndDer(*geoTrue, qp, 1, bs.component(1), By);
        assembleParameter2ndDer1D(bs.component(1), By);
        assembleParameterMass(bs.component(1), My);
        handleDirichletConditions(My, bc, boundary::south, boundary::north);
        handleDirichletConditions(By, bc, boundary::south, boundary::north);

        if (dim == 3)
        {
            //assembleParameterStiffness(bs.component(2), Kz);
            assembleParameter2ndDer1D(bs.component(2), Bz);
            assembleParameterMass(bs.component(2), Mz);
            handleDirichletConditions(Mz, bc, boundary::front, boundary::back);
            handleDirichletConditions(Bz, bc, boundary::front, boundary::back);
        }


        gsInfo << "done." << endl;

        gsInfo << "Determining Kronecker product... " << flush;

        gsSparseMatrix<> B1, B2, B3;
        
        if (dim == 2)
        {
            B1 = By.kron(Mx);
            B2 = My.kron(Bx);

            Bfine = B1 + B2; // \bar{B}
        }
        else if (dim == 3)
        {
            B1 = Bz.kron(My).kron(Mx);
            B2 = Mz.kron(By).kron(Mx);
            B3 = Mz.kron(My).kron(Bx);

            Bfine = B1 + B2 + B3; // \bar{B}
        }
        if  (smoother == Smoother::MassRichardsonSubspaceCorrectionFullBiharmonic)
        {
            gsInfo << "Assembling full parametric B\n";
            gsSparseMatrix<> Kx, Ky, Kz;
            assembleParameterStiffness(bs.component(0), Kx);
            handleDirichletConditions(Kx, bc, boundary::west, boundary::east);
            assembleParameterStiffness(bs.component(1), Ky);
            handleDirichletConditions(Ky, bc, boundary::south, boundary::north);
            if (dim == 3)
            {
                assembleParameterStiffness(bs.component(2), Kz);
                handleDirichletConditions(Kz, bc, boundary::front, boundary::back);
            }
            gsSparseMatrix<> K1, K2, K3;

            if (dim == 2)
            {
                K1 = Ky.kron(Kx);

                Bfine += 2*K1;
            }
            else if (dim == 3)
            {
                K1 = Kz.kron(Ky).kron(Mx);
                K2 = Kz.kron(My).kron(Kx);
                K3 = Mz.kron(Ky).kron(Kx);

                Bfine += 2*K1 + 2*K2 + 2*K3;
            }

        }

    }

    gsInfo << "done." << endl;


    
    /**************************** Setup mg ****************************/

    gsMatrix<> K3 = Bfine;
    gsMatrix<> KfineInv= K3.inverse();
    gsMatrix<> K4 = BfineTrue;
    //gsInfo << "\n" << K3 << "\n\n\n";
    //gsInfo << "\n" << K4 << "\n\n\n";
    gsInfo << "The condition number is :              :" << gsSolverUtils<real_t>::conditionNumber(K4*KfineInv, false) << "\n";

    const double timeAssembling2 = time.stop(); time.restart();
    
    // set up the multigrid solver
    gsInfo << "Constructing smoothers... \n";

    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make(Bfine, transferMatrices);
    // Determine coarse solve time:
    gsTimedOp<>::Ptr coarseSolver = gsTimedOp<>::make("CoarseSolver",mg->coarseSolver(),false);
    mg->setCoarseSolver(coarseSolver);

    mg->setNumPreSmooth( numPreSmooth );
    mg->setNumPostSmooth( numPostSmooth );
    mg->setNumCycles( cycles );
    gsInfo << "Constructing smoothers... \n";

    gsInfo << "Constructing smoothers... " << flush;
    
    for (int i = mg->numLevels() == 1 ? 0 : 1; i < mg->numLevels(); ++i)
    {
         switch( smoother ) {
            case Smoother::Richardson:                            mg->setSmoother(i, makeRichardsonOp(mg->matrix(i),damping)); break;
            case Smoother::Jacobi:                                mg->setSmoother(i, makeJacobiOp(mg->matrix(i),damping)); break;
            case Smoother::GaussSeidel:                           mg->setSmoother(i, makeGaussSeidelOp(mg->matrix(i))); break;
            case Smoother::MassRichardson:                        mg->setSmoother(i, gsPreconditionerFromOp<>::make(mg->underlyingOp(i),makeMassSmootherOperator(bases[i][0], damping, bc) )); break;
            case Smoother::MassRichardsonBoundaryCorrection:      mg->setSmoother(i, gsPreconditionerFromOp<>::make(mg->underlyingOp(i),makeBoundaryCorrectedMassSmootherOperator(bases[i][0], damping, bc) )); break;
            case Smoother::MassRichardsonSubspaceCorrection:      mg->setSmoother(i, gsPreconditionerFromOp<>::make(mg->underlyingOp(i),makeSubspaceCorrectedMassSmootherOperator(bases[i][0], damping, bc) )); break;
            case Smoother::MassRichardsonSubspaceCorrectionBiharmonic:          mg->setSmoother(i, gsPreconditionerFromOp<>::make(mg->underlyingOp(i),makeSubspaceCorrectedMassSmootherOperatorBiharmonic(bases[i][0], damping, bc) )); break;
            case Smoother::MassRichardsonSubspaceCorrectionBiharmonicRank1:     mg->setSmoother(i, gsPreconditionerFromOp<>::make(mg->underlyingOp(i),makeSubspaceCorrectedMassSmootherOperatorBiharmonic(*geoTrue, bases[i][0], damping, bc) )); break;
            case Smoother::MassRichardsonSubspaceCorrectionFullBiharmonic:      mg->setSmoother(i, gsPreconditionerFromOp<>::make(mg->underlyingOp(i),makeSubspaceCorrectedMassSmootherOperatorFullBiharmonic(bases[i][0], damping, bc) )); break;

        }
    }
    gsInfo << "done." << "\n";

    //gsMatrix<> mgMatrix;
    //mg->toMatrix(mgMatrix);
    //gsInfo << "The mg condition number is :              :" << gsSolverUtils<real_t>::conditionNumber(mgMatrix*K3, false) << "\n";


    gsMatrix<> x;
    x.setRandom( mg->nDofs(), 1 );

    const double timeSetup = time.stop(); time.restart();

    const real_t resNorm0 = (rhs - BfineTrue * x).norm();
    gsInfo << "Residual norm:     " << resNorm0 << "\n";
    real_t resNorm(0), oldResNorm = resNorm0;

    int numIter = 0;
    real_t minReduction = 1e6;

    gsConjugateGradient<> cg( BfineTrue, mg );
    cg.setTolerance(tol);

    if (compEigs)
        cg.setCalcEigenvalues(true);

    cg.initIteration( rhs, x );
    
    //real_t errorH2;
    //real_t oldErrorH2=1;
    do
    {
        cg.step(x);
 
        resNorm = (rhs - BfineTrue * x).norm();
        //
        // compute the error
        std::vector<gsMatrix<real_t> > coefs;
        for (size_t s=0; s< phySpace.size(); ++s)
            coefs.push_back(assemblerTrue.reconstructSolution(s,x,eli));
        std::vector<gsFunction<real_t>*> svec(1, &g); //const_cast<gsFunctionExpr<real_t>*>(&data.solution));
        /*gsRecipeDistance dist(pdeBiharmonic.domain(),phySpace,coefs,svec);
        gsMatrix<> norms(3,2);
        norms<<0,2,1,2,2,2;
        dist.setSeminorms(0,norms);
        dist.assemble();

        gsMatrix<> error=dist.getDistance(0).rowwise().sum();
        errorH2 =math::sqrt(error(2,0));
        */
        //
        //gsInfo << "Residual norm:   " << left << setw(15) << resNorm
        //       << "     reduction:  1 / " << setprecision(3) << left << setw(8)
        //       << (oldResNorm/resNorm) << setprecision(6)<< "\n";
        //       /*<< " H2 semi norm error:   " << left << setw(15) << errorH2
        //       << "     reduction:  1 / " << setprecision(3)
        //       << (oldErrorH2/errorH2) << setprecision(6) << "\n";*/
        minReduction = math::min(minReduction, oldResNorm/resNorm);
        oldResNorm = resNorm;
        //oldErrorH2 = errorH2;
        ++numIter;
    } while (/*resNorm > tol && */ resNorm / resNorm0 > tol && gsIsfinite(resNorm));

    const double timeSolve = time.stop();
    gsInfo << "########################\n";
    if (resNorm / resNorm0 > tol || !gsIsfinite(resNorm))
        gsInfo << "Did not converge.\n";
    else    
        gsInfo << "Converged in " << numIter << " iterations.\n";
    gsInfo << "########################\n";
    gsInfo << "Average convergence factor:  1 / " << setprecision(3) << math::pow(resNorm0 / resNorm, 1.0 / numIter) << setprecision(6) << "\n";
    gsInfo << "Worst   convergence factor:  1 / " << setprecision(3) << minReduction << setprecision(6) << "\n";
    gsInfo << "\n";

    gsInfo << "Assembling time: "; formatTime(gsInfo, timeAssembling1+timeAssembling2); gsInfo << "\n";
    gsInfo << "  aux problem:   "; formatTime(gsInfo, timeAssembling1);                 gsInfo << "\n";
    gsInfo << "  the problem:   "; formatTime(gsInfo, timeAssembling2);                 gsInfo << "\n";
    gsInfo << "MG time:         "; formatTime(gsInfo, timeSetup+timeSolve);             gsInfo << "\n";
    gsInfo << " setup:          "; formatTime(gsInfo, timeSetup);                       gsInfo << "\n";
    gsInfo << " solving:        "; formatTime(gsInfo, timeSolve);
    gsInfo << "         (avg. ";   formatTime(gsInfo, timeSolve/numIter);               gsInfo << " per iteration)" << "\n";
    gsInfo << "  coarse solver: "; formatTime(gsInfo, coarseSolver->getTime());
    gsInfo << "\n";
    gsInfo << "Total time:      "; formatTime(gsInfo, timeAssembling1+timeAssembling2+timeSetup+timeSolve); gsInfo << "\n";
    gsInfo << "\n";



    //gsSparseMatrix<> sys = assemblerTrue.getSystemMatrix();
    gsMatrix<> solVec = x;
    //gsInfo << x <<"\n";

    // compute the error
    /*
    std::vector<gsMatrix<real_t> > coefs;
    for (size_t s=0; s< phySpace.size(); ++s)
        coefs.push_back(assemblerTrue.reconstructSolution(s,solVec,eli));
    std::vector<gsFunction<real_t>*> svec(1, &g); //const_cast<gsFunctionExpr<real_t>*>(&data.solution));
    gsRecipeDistance dist(pdeBiharmonic.domain(),phySpace,coefs,svec);
    gsMatrix<> norms(3,2);
    norms<<0,2,1,2,2,2;
    dist.setSeminorms(0,norms);
    dist.assemble();

    gsMatrix<> error=dist.getDistance(0).rowwise().sum();
    real_t errorL2 =math::sqrt(error(0,0));
    real_t errorH1 =math::sqrt(error(1,0));
    real_t errorH22 =math::sqrt(error(2,0));

    gsInfo << " The L2 error is:      " << errorL2 << "\n";
    gsInfo << " The H1 semi error is: " << errorH1 << "\n";
    gsInfo << " The H2 semi error is: " << errorH22 << "\n";
    */
     //errorH1=math::sqrt(error(1,0));
    //freeAll(phySpace);
    
    real_t max = 0., min = 1.e12, essMin = 1.e12;
    if (compEigs)
    {
        //gsInfo << "Condition number: " << cg.getConditionNumber() << std::endl;
        gsMatrix<real_t> eigs;
        cg.getEigenvalues(eigs);
        //gsInfo << "Eigenvalues: " << eigs.transpose() << std::endl;
        for ( index_t i = 0; i < eigs.rows(); ++i )
        {
            if (eigs(i) > max)
                max = eigs(i);
            if (eigs(i) < min)
                min = eigs(i);
            if (eigs(i) < essMin && eigs(i) > 1.e-5)
                essMin = eigs(i);
        }
        gsInfo << "Eigenvalues: Minimum, essentialMinimum, maximum, essentialConditionNumber: " << min << ", " << essMin << ", " << max << ", " << (max/essMin) << std::endl;
    }
    
    if (writeLog)
    {
        fstream log("out.txt", fstream::out | fstream::app);

        log << "multigrid_preconditioner_example\t"
            << alpha << "\t"
            << cycles << "\t"
            << degree << "\t"
            << numRefine << "\t"
            << numLevels << "\t"
            << damping << "\t"
            << smoother_name << "\t"
            << numPreSmooth << "\t"
            << numPostSmooth << "\t"
            << numIter << "\t"
            << math::pow(resNorm0 / resNorm, 1.0 / numIter) << "\t"
            << timeAssembling2 << "\t"
            << (timeAssembling1+timeSetup+timeSolve) /*MG time*/;
            
            if (compEigs)
                log << "\t" << min << "\t" << essMin << "\t" << max << "\t" << (max/essMin) ;
            
            log << "\n";
    }

    delete tbasis;
    
    freeAll( phySpace );

    return (resNorm / resNorm0 > tol) ? 1 : 0;
};
