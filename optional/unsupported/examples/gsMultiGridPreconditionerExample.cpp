/** @file gsMultiGridExample.cpp

    @brief Provides test examples for multigrid algorithms, where the multigrid method
    on the parameter domain is used as a preconditioner for the problem on the physical
    domain.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, S. Takacs
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsSolver/gsKronecker.h>
#include <gsSolver/gsTimedOp.h>
#include <gsMultiGrid/gsMassSmoother.h>
#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>
#include <gsIO/gsCmdLineWithEnumSupport.h>

namespace Smoother {
    enum type {
        Richardson,
        Jacobi,
        GaussSeidel,
        MassRichardson,
        MassRichardsonBoundaryCorrection,
        MassRichardsonSubspaceCorrection
    };
}

using namespace std;
using namespace gismo;


class gsRecipeAssemblerGenericSecondOrderOp : public gsRecipeAssemblerPoisson
{
protected:
    gsFunction<real_t> *m_A; 
    gsFunction<real_t> *m_b; 
    gsFunction<real_t> *m_c;
public:
     gsRecipeAssemblerGenericSecondOrderOp( const gsPoissonPde<real_t> &pde, gsFunction<real_t> *A=NULL,gsFunction<real_t> *b=NULL,gsFunction<real_t> *c=NULL)
        : gsRecipeAssemblerPoisson(pde), m_A(A), m_b(b), m_c(c)
    {}

protected:
    virtual gsRecipe<real_t> getPatchRecipe(index_t ) // patch)
    {
        gsRecipe<real_t>            result;
        gsRecipeIngredient<real_t>  ingr;
    
        ingr.setOperator(new gsGenericSecondOrderOp<real_t>(m_A,m_b,m_c));
        //ingr.setOperator(new gsGradGradOp<real_t>());
        ingr.setTestSpace(0);
        ingr.setUnknownSpace(0);
        ingr.setRule(getSysWriter());
        result.add(ingr);
    
        ingr.setOperator(new gsL2TestOp<real_t>(*m_pde.rhs()));
        ingr.setTestSpace(0);
        ingr.setUnknownSpace(0);
        ingr.setRule(getRhsWriter());
        result.add(ingr);

        return result;
    }

};

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
    index_t degree = 2;
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
    bool matrixFree = false;

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
        .add(Smoother::Richardson,                               "r",      "richardson"                                                       )
        .add(Smoother::Jacobi,                                   "j",      "jacobi"                                                           )
        .add(Smoother::GaussSeidel,                              "gs",     "gauss-seidel"                                                     )
        .add(Smoother::MassRichardson,                           "mr",     "mass-richardson"                                                  )
        .add(Smoother::MassRichardsonBoundaryCorrection,         "mrb",    "mass-richardson-boundary-correction"                              )
        .add(Smoother::MassRichardsonSubspaceCorrection,         "mrs",    "mass-richardson-subspace-correction"                              )
        .writeDescOfChosenOptTo(smoother_name);
    cmd.addReal("", "damping",
                "Damping factor for the smoother", damping);
    cmd.addInt("c", "cycles",
               "Number of multi-grid cycles", cycles);
    cmd.addReal("a", "alpha", "alpha in \"-div(A grad(u)) + alpha u = f\"", alpha);
    cmd.addSwitch("log", "Write results to log file", writeLog);
    cmd.addSwitch("eigs", "Compute eigenvalues of the preconditioned system.", compEigs);
    cmd.addSwitch("matrixfree", "Matrix free implementation of the solver", matrixFree);
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
            f = gsFunctionExpr<>("(pi^2 ) * sin(pi*x)",1);
            g = gsFunctionExpr<>("sin(pi*x)",1);
            break;
        case 2:
            f = gsFunctionExpr<>("(2*pi^2 ) * sin(pi*x) * sin(pi*y)",2);
            g = gsFunctionExpr<>("sin(pi*x) * sin(pi*y)",2);
            break;
        case 3:
            f = gsFunctionExpr<>("(3*pi^2 ) * sin(pi*x) * sin(pi*y) * sin(pi*z)",3);
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
            case Smoother::Richardson:                                     damping = 0.80; break;
            case Smoother::Jacobi:                                         damping = 0.80; break;
            case Smoother::GaussSeidel:                                    break;
            case Smoother::MassRichardson:                                 damping = 0.25 / ( (degree+1.)*(degree+1.) ); break;
            case Smoother::MassRichardsonBoundaryCorrection:               damping = 0.09; break;
            case Smoother::MassRichardsonSubspaceCorrection:               damping = 0.09; break;
        }
    }

    gsInfo << "Source function: " << f << ".\n" << "\n";
    gsInfo << "Exact solution:  " << g << ".\n" << "\n";

    // set up boundary conditions

    gsConstantFunction<> zero(0.0, dim);
    gsConstantFunction<> one (1.0, dim);

    gsBoundaryConditions<> bc;
    
    bc.addCondition( boundary::west,  condition_type::dirichlet, &zero );
    bc.addCondition( boundary::east,  condition_type::dirichlet, &one );
    
    bc.addCondition( boundary::south, condition_type::neumann, &zero );
    bc.addCondition( boundary::north, condition_type::neumann, &zero );

    if (dim >= 3)
    {
        bc.addCondition( boundary::front, condition_type::neumann, &zero );
        bc.addCondition( boundary::back,  condition_type::neumann, &zero );
    }

    gsBasis<>::uPtr tbasis = geo->basis().clone();
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

    gsInfo << "Setup assembler... " << flush;
    gsMultiPatch<> mp(*geo);
    gsGenericAssembler<> assm(mp, bases.back(), gsAssembler<>::defaultOptions(), &bc);
    assm.computeDirichletDofs();
    gsInfo << "done." << endl;
   
    gsSparseMatrix<> Kfine;
    
    vector< gsLinearOperator<>::Ptr > allKops(0), allProlongs(0), allRestricts(0);
    vector< gsSparseMatrix<> > allKx, allKy, allKz,
                               allMx, allMy, allMz;
    
    if ( !matrixFree )
    {
        gsInfo << "Assembling mass and stiffness matrix in 1D... " << flush;
        const gsBasis<>& bs = bases.back()[0];
              
        gsSparseMatrix<> Kx, Ky, Kz, Mx, My, Mz;

        assembleParameterStiffness(bs.component(0), Kx);
        assembleParameterMass(bs.component(0), Mx);
        handleDirichletConditions(Mx, bc, boundary::west, boundary::east);
        handleDirichletConditions(Kx, bc, boundary::west, boundary::east);

        assembleParameterStiffness(bs.component(1), Ky);
        assembleParameterMass(bs.component(1), My);
        handleDirichletConditions(My, bc, boundary::south, boundary::north);
        handleDirichletConditions(Ky, bc, boundary::south, boundary::north);

        if (dim == 3)
        {
            assembleParameterStiffness(bs.component(2), Kz);
            assembleParameterMass(bs.component(2), Mz);
            handleDirichletConditions(Mz, bc, boundary::front, boundary::back);
            handleDirichletConditions(Kz, bc, boundary::front, boundary::back);
        }
        
        gsInfo << "done." << endl;
        
        gsInfo << "Determining Kronecker product... " << flush;
        
        gsSparseMatrix<> M, K1, K2, K3;

        if (dim == 2)
        {
            M  = My.kron(Mx);
            K1 = Ky.kron(Mx);
            K2 = My.kron(Kx);

            Kfine = M + K1 + K2; // -Laplace u + u
        }
        else if (dim == 3)
        {
            M  = Mz.kron(My).kron(Mx);
            K1 = Kz.kron(My).kron(Mx);
            K2 = Mz.kron(Ky).kron(Mx);
            K3 = Mz.kron(My).kron(Kx);

            Kfine = M + K1 + K2 + K3; // -Laplace u + u
        }

        gsInfo << "done." << endl;
    }
    else    // matrix-free multigrid implementation
    {
        const index_t bsz = bases.size();

        allKx.resize(bsz);
        allKy.resize(bsz);
        allMx.resize(bsz);
        allMy.resize(bsz);
        if (dim == 3)
        {
            allKz.resize(bsz);
            allMz.resize(bsz);
        }
        allKops.resize(bsz);
        allProlongs.resize(bsz-1);
        allRestricts.resize(bsz-1);
        
        for ( index_t lv=0; lv<bsz; ++lv )
        {

            gsBasis<>& bs = bases[lv][0];

            gsInfo << "Assembling mass and stiffness matrix in 1D for level" << lv << "... " << flush;
                
            assembleParameterStiffness(bs.component(0), allKx[lv]);
            assembleParameterMass     (bs.component(0), allMx[lv]);

            assembleParameterStiffness(bs.component(1), allKy[lv]);
            assembleParameterMass     (bs.component(1), allMy[lv]);

            handleDirichletConditions(allMx[lv], bc, boundary::west, boundary::east);
            handleDirichletConditions(allKx[lv], bc, boundary::west, boundary::east);

            handleDirichletConditions(allMy[lv], bc, boundary::south, boundary::north);
            handleDirichletConditions(allKy[lv], bc, boundary::south, boundary::north);

            if (dim == 3)
            {
                assembleParameterStiffness(bs.component(2), allKz[lv]);
                assembleParameterMass     (bs.component(2), allMz[lv]);

                handleDirichletConditions(allMz[lv], bc, boundary::front, boundary::back);
                handleDirichletConditions(allKz[lv], bc, boundary::front, boundary::back);
            }

            gsInfo << "done." << endl;
            
            gsInfo << "Determining Kronecker product... " << flush;
            
            gsSumOp<>::Ptr op = gsSumOp<>::make();
            if (dim == 2)
            {
                op->addOperator( gsKroneckerOp<>::make( makeMatrixOp(allMy[lv]), makeMatrixOp(allMx[lv]) ) );     // MM
                op->addOperator( gsKroneckerOp<>::make( makeMatrixOp(allKy[lv]), makeMatrixOp(allMx[lv]) ) );     // KM
                op->addOperator( gsKroneckerOp<>::make( makeMatrixOp(allMy[lv]), makeMatrixOp(allKx[lv]) ) );     // MK
            }
            else if (dim == 3)
            {
                op->addOperator( gsKroneckerOp<>::make( makeMatrixOp(allMz[lv]), makeMatrixOp(allMy[lv]), makeMatrixOp(allMx[lv]) ) );    // MMM
                op->addOperator( gsKroneckerOp<>::make( makeMatrixOp(allKz[lv]), makeMatrixOp(allMy[lv]), makeMatrixOp(allMx[lv]) ) );    // KMM
                op->addOperator( gsKroneckerOp<>::make( makeMatrixOp(allMz[lv]), makeMatrixOp(allKy[lv]), makeMatrixOp(allMx[lv]) ) );    // MKM
                op->addOperator( gsKroneckerOp<>::make( makeMatrixOp(allMz[lv]), makeMatrixOp(allMy[lv]), makeMatrixOp(allKx[lv]) ) );    // MMK
            }
            allKops[lv] = op;

            gsInfo << "done." << endl;
        }
        
        for ( index_t lv=0; lv<bsz-1; ++lv )
        {
            allProlongs[lv] = makeMatrixOp(transferMatrices[lv]);
            allRestricts[lv] = makeMatrixOp(transferMatrices[lv].transpose());
        }
    }

    /******************* Setup the "true" problem ******************/
    
    const double timeAssembling1 = time.stop(); time.restart();
    
    gsGeometry<>::Ptr geoTrue;
    gsMultiPatch<>::Ptr geoTrueMp;
        switch (dim)
    {
        case 2: geoTrue = gsNurbsCreator<>::NurbsQuarterAnnulus(0.5, 1.0); break;
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

    gsInfo << "Setup gsRecipeAssembler for true problem... " << flush;
    gsMultiPatch<> mpTrue(*geoTrue);
    
    gsPoissonPde<real_t> pde(mpTrue,bc,f,&g);
    
    // set diffusion coefficient matrix
    gsFunction<>::Ptr A;
    switch (dim)
    {
        case 2:
            A = memory::make_shared(new gsFunctionExpr<>(
                // matrix with
                //  - eigenvectors: (x,y)/r, (-y,x)/r
                //  - eigenvalues :    2   ,     1
                    //"(2*x^2 + y^2)/(x^2 + y^2)", "(x*y)/(x^2 + y^2)",
                    //"(x*y)/(x^2 + y^2)", "(x^2 + 2*y^2)/(x^2 + y^2)",

                // matrix with
                //  - eigenvectors (y,x),   (-x,y)
                //  - eigenvalues    1    1+x^2+y^2
                    "1 + x^2", "-x*y",
                    "-x*y", "1 + y^2",
                2));
            break;

        case 3:
            A = memory::make_shared(new gsFunctionExpr<>(
                    "1 + x^2", "-x*y/3", "-x*z/3",
                    "-x*y/3", "1 + y^2", "-y*z/3",
                    "-x*z/3", "-y*z/3", "1 + z^2",
                3));
            break;
    }
    gsFunction<>::Ptr c = memory::make_shared(new gsConstantFunction<>( alpha, dim ));
    gsRecipeAssemblerGenericSecondOrderOp assemblerTrue(pde, &*A, NULL, &*c);
    assemblerTrue.setDirichletStrategy(dirichlet::elimination);
    assemblerTrue.setZeroAverage(false);
    vector<gsPhysicalSpace*> phySpace;
    // use the finest basis from the reference domain for the true problem
    phySpace.push_back(new gsPhysicalSpaceScalar(bases.back(),mpTrue,INVERSE_COMPOSITION));
    assemblerTrue.setSpace(phySpace);
    
    gsInfo << "done." << endl;
            
    gsInfo << "Assembling true system matrix... " << flush;
    assemblerTrue.assemble();
    gsSparseMatrix<real_t> KfineTrue = assemblerTrue.getSystemMatrix();
    gsInfo << "done, " << KfineTrue.rows() << " dofs." << "\n";  
    
    /**************************** Setup mg ****************************/
     
    const double timeAssembling2 = time.stop(); time.restart();
    
    // set up the multigrid solver
    gsMultiGridOp<>::Ptr mg;
    if ( !matrixFree )
    {
        mg = gsMultiGridOp<>::make(Kfine, transferMatrices);
    }
    else    // matrix-free multigrid implementation
    {
        mg = gsMultiGridOp<>::make(allKops, allProlongs, allRestricts);
        if ( smoother == Smoother::Richardson || smoother == Smoother::Jacobi || smoother == Smoother::GaussSeidel )
        {
            gsInfo << "The chosen smoother is not matrix free." << endl;
            return -1;
        }
    }
    
    // Determine coarse solve time:
    gsTimedOp<>::Ptr coarseSolver = gsTimedOp<>::make("CoarseSolver",mg->coarseSolver(),false);
    mg->setCoarseSolver(coarseSolver);

    mg->setNumPreSmooth( numPreSmooth );
    mg->setNumPostSmooth( numPostSmooth );
    mg->setNumCycles( cycles );

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
        }
    }
    gsInfo << "done." << "\n";

    gsMatrix<> x;
    x.setRandom( mg->nDofs(), 1 );

    const double timeSetup = time.stop(); time.restart();

    const real_t resNorm0 = (assemblerTrue.getSystemRhs() - KfineTrue * x).norm();
    gsInfo << "Residual norm:     " << resNorm0 << "\n";
    real_t resNorm(0), oldResNorm = resNorm0;

    int numIter = 0;
    real_t minReduction = 1e6;

    gsConjugateGradient<> cg( KfineTrue, mg );
    cg.setTolerance(tol);

    if (compEigs)
        cg.setCalcEigenvalues(true);

    cg.initIteration( assemblerTrue.getSystemRhs(), x );
    
    do
    {
        cg.step(x);
 
        resNorm = (assemblerTrue.getSystemRhs() - KfineTrue * x).norm();
        gsInfo << "Residual norm:     " << left << setw(15) << resNorm << "          reduction:  1 / " << setprecision(3) << (oldResNorm/resNorm) << setprecision(6) << "\n";
        minReduction = math::min(minReduction, oldResNorm/resNorm);
        oldResNorm = resNorm;
        ++numIter;
    } while (/*resNorm > tol && */ resNorm / resNorm0 > tol && gsIsfinite(resNorm));

    const double timeSolve = time.stop();
    if (resNorm / resNorm0 > tol || !gsIsfinite(resNorm))
        gsInfo << "Did not converge.\n";
    else    
        gsInfo << "Converged in " << numIter << " iterations.\n";
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
    gsInfo << "Total time:      "; formatTime(gsInfo, timeAssembling1+timeAssembling2+timeSetup+timeSolve);                       gsInfo << "\n";
    gsInfo << "\n";

    // Compute L2 error
    if( !true )
    {
        // TODO: how to do this with the gsRecipeAssembler ?
        // TODO: that might work wrong
        gsField<> sol;
        gsGenericAssembler<real_t> assmTrue(mpTrue, bases.back(), gsAssembler<>::defaultOptions(), &bc);
        assmTrue.computeDirichletDofs();
        sol = assmTrue.constructSolution(x);
        const real_t error = computeL2Distance(sol, g, false, 3*mg->nDofs());
        gsInfo << "L2 error: " << error << "\n";
    }
    
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
    
    freeAll( phySpace );

    return (resNorm / resNorm0 > tol) ? 1 : 0;
};
