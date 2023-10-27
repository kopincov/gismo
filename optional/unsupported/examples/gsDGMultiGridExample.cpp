/** @file gsDGMultiGridExample.cpp

    @brief Provides test examples for multigrid algorithms

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <iostream>
#include <iomanip>
#include <ctime>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsCore/gsLinearCombinationOfFunctionsFunction.h>
#include <gsMultiGrid/gsMassSmoother.h>
#include <gsSolver/gsTimedOp.h>
#include <gsTensor/gsTensorTools.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsGeneralizedPoissonAssembler.h>
#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>

#include <gsIO/gsCmdLineWithEnumSupport.h>


using namespace std;
using namespace gismo;

namespace Smoother {
    enum type {
        Richardson,
        Jacobi,
        GaussSeidel,
        MassRichardsonSubspaceCorrectionPiecewise,
        MassRichardsonSubspaceCorrectionPatchwise,
        MassRichardsonSubspaceCorrectionPiecewiseQuasiExact,
        MassRichardsonSubspaceCorrectionPatchwiseQuasiExact
    };
}

std::vector<gsLinearOperator<>::Ptr> makePiecewiseExactSolvers(const gsSparseMatrix<>& A,
                            const gsMultiBasis<>& mb,
                            // const gsBoundaryConditions<>& bc,
                            const gsOptionList& opt)
{
    const index_t n = mb.nBases();
    std::vector<gsLinearOperator<>::Ptr> ops(n);

    gsBoundaryConditions<> localbc; //TODO: THAT IS WRONG FOR NEUMANN BC!
    for (index_t i=0; i<n; ++i)
    {
        const index_t d = mb[i].dim();
        for ( index_t ps=0; ps < (1<<d); ++ps )
            localbc.addCondition( i, 1+ps, condition_type::dirichlet, NULL );
    }

    gsDofMapper dm;
    mb.getMapper(
       (dirichlet::strategy)opt.askInt("DirichletStrategy",11),
       (iFace    ::strategy)opt.askInt("InterfaceStrategy", 1),
       localbc,
       dm,
       0
    );

    for (index_t i=0; i<n; ++i)
    {
        const index_t nDofs = mb[i].size();

        gsSparseEntries<> tmp;
        tmp.reserve(nDofs);
        index_t counter = 0;
        for (index_t j=0; j<nDofs; ++j)
        {
            const index_t dof_idx = dm.index(j,i);
            if (dm.is_free_index(dof_idx)) //TODO: check this!
            {
                tmp.add(counter,dof_idx,1);
                counter++;
            }
        }
        gsSparseMatrix<> transfer( counter, A.rows() );
        transfer.setFromTriplets(tmp.begin(), tmp.end());
        transfer.makeCompressed();

        gsSparseMatrix<> local_mat = transfer * A * transfer.transpose();

        ops[i] = makeSparseCholeskySolver( local_mat );
    }
    return ops;
}

std::vector< gsLinearOperator<>::Ptr > printSizes( const std::vector< gsLinearOperator<>::Ptr >& in )
{
    for (size_t i=0; i<in.size(); ++i)
        gsInfo << "local solver: " << in[i]->rows() << "x" << in[i]->cols() << "\n";
    return in;
}

/// Allows to setup boundary conditions from command line
void bcChoose( char bc, condition_type::type& bc_type, gsFunction<> *& bc_func, gsFunction<> * bc_func_dirichlet, gsFunction<> * bc_func_neumann )
{
    if( bc == 'd' )
    {
        bc_type = condition_type::dirichlet;
        bc_func = bc_func_dirichlet;
    }
    else if( bc == 'n' )
    {
        bc_type = condition_type::neumann;
        bc_func = bc_func_neumann;
    }
    else
    {
        cerr << "Invalid boundary condition. Allowed are: dirichlet (d), neumann (n).\n";
        exit(-1);
    }
}

void computeMatrixLevels(const std::vector< gsSparseMatrix<real_t,RowMajor> >& transfer, const gsSparseMatrix<real_t>& M, std::vector< gsSparseMatrix<real_t> >& Mlevels)
{
    const index_t sz = transfer.size() + 1;
    Mlevels.resize( sz );
    Mlevels[ sz-1 ] = M;

    for (index_t i = sz-2; i >= 0; --i)
    {
        Mlevels[i] = transfer[i].transpose() * Mlevels[i+1] * transfer[i];
    }
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


gsPiecewiseFunction<real_t> makeUniformPiecewiseFunction( const gsFunction<real_t> & func, index_t numberOfPieces )
{
    gsPiecewiseFunction<real_t>::FunctionContainer funcs( numberOfPieces );
    for( index_t i=0; i<numberOfPieces; ++i)
        funcs[i] = func.clone().release();
    return gsPiecewiseFunction<real_t>( funcs );
}

int main(int argc, char *argv[])
{
    string geometry("yeti_mp2.xml"); //yeti_mp2.xml //two_squares.xml
    string boundaryCondition("d");
    index_t numRefine = 4;
    index_t degree = 2;
    index_t numLevels = -1;
    index_t numSplit = 0;
    Smoother::type smoother = Smoother::MassRichardsonSubspaceCorrectionPiecewise;
    std::string smoother_name;
    index_t numPreSmooth = 1;
    index_t numPostSmooth = 1;
    index_t cycles = 1;
    real_t gamma = 1.;
    bool useCG = false;
    bool monitorL2 = false;
    bool monitorDG = false;
    bool writeLog = false;
    real_t tol = 1e-8;
    real_t damping = -1;
    real_t outerDamping = -1;
    real_t coarseDamping = 1.;
    bool fancy = false;
    index_t maxIter = 1000;
    bool notGalerkin = false;

    gsCmdLineWithEnumSupport cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addString("g", "geometry",           "Specification of the geometry (overrides dimension)",             geometry         );
    cmd.addString("b", "boundary-condition", "Boundary condition",                                              boundaryCondition);
    cmd.addInt   ("r", "uniformRefine",      "Number of uniform h-refinement steps to perform before solving",  numRefine        );
    cmd.addInt   ("p", "degree",             "Degree of the B-spline discretization space",                     degree           );
    cmd.addInt   ("l", "levels",             "Number of levels to use for multigrid iteration",                 numLevels        );
    cmd.addInt   ("",  "split",              "Split every patch uniformly into 2^d patches (default: 0)",       numSplit         );
    cmd.addInt   ("",  "presmooth",          "Number of pre-smoothing steps",                                   numPreSmooth     );
    cmd.addInt   ("",  "postsmooth",         "Number of post-smoothing steps",                                  numPostSmooth    );
    cmd.addEnum  ("s", "smoother",           "Smoothing method",                                                smoother         )
        .add(Smoother::Richardson,                                          "r",      "Richardson smoother"                                              )
        .add(Smoother::Jacobi,                                              "j",      "Jacobi smoother"                                                  )
        .add(Smoother::GaussSeidel,                                         "gs",     "GaussSeidel smoother"                                             )
        .add(Smoother::MassRichardsonSubspaceCorrectionPiecewise,           "mrs-ad", "mass smoother with piecewise subspace correction"                 )
        .add(Smoother::MassRichardsonSubspaceCorrectionPatchwise,           "mrs-a",  "mass smoother with patchwise subspace correction"                 )
        .add(Smoother::MassRichardsonSubspaceCorrectionPiecewiseQuasiExact, "mrs-adx","mass smoother with piecewise subspace correction (quasi exact)"   )
        .add(Smoother::MassRichardsonSubspaceCorrectionPatchwiseQuasiExact, "mrs-ax", "mass smoother with patchwise subspace correction (quasi exact)"   )
        .writeDescOfChosenOptTo(smoother_name);
    cmd.addReal  ("",  "damping",            "Damping factor for the smoother (handed over to smoother)",       damping          );
    cmd.addReal  ("",  "outerdamping",       "Damping factor for the smoother (globally)",                      outerDamping     );
    cmd.addReal  ("",  "coarsedamping",      "Damping factor for the coarse grid correction",                   coarseDamping    );
    cmd.addInt   ("c", "cycles",             "Number of multi-grid cycles",                                     cycles           );
    cmd.addInt   ("",  "maxiter",            "Maximum number of iterations",                                    maxIter          );
    cmd.addReal  ("",  "gamma",              "gamma in \"- gamma LAPLACE u  = f\"",                             gamma            );
    cmd.addSwitch(     "log",                "Write results to log file",                                       writeLog         );
    cmd.addSwitch(     "cg",                 "Use CG iteration",                                                useCG            );
    cmd.addSwitch(     "monitor-L2",         "Monitor the L2 error",                                            monitorL2        );
    cmd.addSwitch(     "monitor-dG",         "Monitor the errors in the DG norm",                               monitorDG        );
    cmd.addReal  ("",  "tol",                "Tolerance for multigrid solver stopping criterion",               tol              );
    cmd.addSwitch(     "fancy",              "Make the discretization non-matching (variate degree&size)",      fancy            );
    cmd.addSwitch(     "not-galerkin",       "Do not construct coarser level matrices using Galrkin projection",notGalerkin      );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    /******************** Define Geometry ********************/

    gsStopwatch time;
    gsMultiPatch<> mp;

    if ( gsFileManager::fileExists(geometry) )
    {
        geometry = gsFileManager::find(geometry);
        gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
        if (!mpPtr) { cerr << "No multipatch object found in file " << geometry << ".\n"; return 1; }
        mp = *mpPtr;
    }
    else
    {
        cerr << "Could not find file.\n";
        return -1;
    }

    for (index_t i=0; i<numSplit; ++i)
    {
        gsInfo << "Split multipatch object uniformly... " << flush;
        mp = mp.uniformSplit();
        gsInfo << "done." << endl;
    }

    gsFunction<>::Ptr f0;
    gsFunction<>::Ptr g;

    switch (mp.geoDim())
    {
        case 2:
            f0 = memory::make_shared( new gsFunctionExpr<>("(2*pi^2 ) * sin(pi*x) * sin(pi*y)",2) );
            g  = memory::make_shared( new gsFunctionExpr<>("sin(pi*x) * sin(pi*y)",2) );
            //f0 = memory::make_shared( new gsFunctionExpr<>("pi^2 * sin(pi*x) * sqrt((y-3)*(y-3))",2) );
            //g  = memory::make_shared( new gsFunctionExpr<>("sin(pi*x) * sqrt((y-3)^2)",2) );
            break;
        case 3:
            f0 = memory::make_shared( new gsFunctionExpr<>("(3*pi^2 ) * sin(pi*x) * sin(pi*y) * sin(pi*z)",3) );
            g  = memory::make_shared( new gsFunctionExpr<>("sin(pi*x) * sin(pi*y) * sin(pi*z)",3) );
            break;
        default:
            cerr << "Invalid geometry dimension:"<<mp.geoDim()<<"\n";
            return -1;
    }

    gsFunction<>::Ptr f = gsLinearCombinationOfFunctionsFunction<>::make(1,f0,0,g);

    if (numRefine < 0)
    {
        cerr << "Number of refinements must be positive.\n"; return -1;
    }
    if (numLevels < 1)
    {
        numLevels = numRefine + 1;
        if( smoother == Smoother::MassRichardsonSubspaceCorrectionPiecewise
            || smoother == Smoother::MassRichardsonSubspaceCorrectionPatchwise
            || smoother == Smoother::MassRichardsonSubspaceCorrectionPiecewiseQuasiExact
            || smoother == Smoother::MassRichardsonSubspaceCorrectionPatchwiseQuasiExact
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
            case Smoother::Richardson:                                           damping = 0.80; break;
            case Smoother::Jacobi:                                               damping = 0.80; break;
            case Smoother::GaussSeidel:                                          damping = 1.00; break;
            case Smoother::MassRichardsonSubspaceCorrectionPiecewise:            damping = 0.25; break;
            case Smoother::MassRichardsonSubspaceCorrectionPatchwise:            damping = 0.25; break;
            case Smoother::MassRichardsonSubspaceCorrectionPiecewiseQuasiExact:  damping = 1.00; break;
            case Smoother::MassRichardsonSubspaceCorrectionPatchwiseQuasiExact:  damping = 1.00; break;
        }
    }
    if (outerDamping < 0)
    {
        // Here, one could add appropriate parameter-choice rules (depending on the degree)
        switch (smoother)
        {
            case Smoother::Richardson:                                           outerDamping = 1.00; break;
            case Smoother::Jacobi:                                               outerDamping = 1.00; break;
            case Smoother::GaussSeidel:                                          outerDamping = 1.00; break;
            case Smoother::MassRichardsonSubspaceCorrectionPiecewise:            outerDamping = 0.75; break;
            case Smoother::MassRichardsonSubspaceCorrectionPatchwise:            outerDamping = 0.75; break;
            case Smoother::MassRichardsonSubspaceCorrectionPiecewiseQuasiExact:  outerDamping = 0.70; break;
            case Smoother::MassRichardsonSubspaceCorrectionPatchwiseQuasiExact:  outerDamping = 0.75; break;
        }
    }

    // ---------------------------------------------------------------------------------------------------------------------------

    gsInfo << "Source function: " << *f << ".\n" << "\n";
    gsInfo << "Exact solution:  " << *g << ".\n" << "\n";


    // set up boundary conditions

    gsConstantFunction<> zero(0.0, mp.geoDim());
    //gsConstantFunction<> one (1.0, mp.geoDim());
    //

    gsBoundaryConditions<> bc;
    condition_type::type bc_type;
    gsFunction<> * bc_func;

    if (boundaryCondition.length() != 1)
    {
        gsWarn << "Only one boundary condition is acceptable for multipatch domains.\n";
        return -1;
    }

    {
        index_t i = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
        {
            ++i;
            bcChoose( boundaryCondition[0], bc_type, bc_func, &*g, &zero );
            bc.addCondition( *it, bc_type, bc_func );
        }
        gsInfo << "Added " << i << " boundary condtions.\n";
    }

    gsMultiBasis<> mb(mp);

    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    vector< gsMultiBasis<> > bases;
    vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;

    gsOptionList assemblerOptions = gsAssembler<>::defaultOptions();

    assemblerOptions.setInt( "InterfaceStrategy", iFace::dg );

    for (int i = 0; i < numRefine - numLevels + 1; ++i)
        mb.uniformRefine();       // refine until coarsest level


    if (fancy)
    {
        gsInfo << "Make uniform refinement for every third patch...\n";
        for (size_t i = 0; i < mb.nBases(); ++i)
            if ( i%3 == 0 )
                mb[i].uniformRefine();
        gsInfo << "done.\n";

        gsInfo << "Increase spline degree for every third patch...\n";
        for (size_t i = 0; i < mb.nBases(); ++i)
            if ( i%3 == 1 )
                mb[i].setDegreePreservingMultiplicity(degree+1);
        gsInfo << "done.\n";

    }


    // set up the hierarchy of spaces and transfer matrices between them
    const index_t refineKnots = 1, mult = 1;
    gsGridHierarchy<>::buildByRefinement(give(mb), bc, assemblerOptions, numLevels, refineKnots, mult)
        .moveMultiBasesTo(bases)
        .moveTransferMatricesTo(transferMatrices);

    gsInfo << (useCG ? "CG preconditioned by multigrid" : "Multigrid");
    gsInfo << " with " << numLevels << " levels using " << smoother_name;
    gsInfo << "(" << numPreSmooth << "," << numPostSmooth << ") smoother" << " and " << cycles << "-cycle.\n";
    for (size_t i=0; i<mb.nBases(); ++i)
        gsInfo << "Coarse discretization space: dim=" << mb[i].dim() << " deg=" << mb[i].degree(0) << " dofs=" << mb[i].size() << "\n";

    gsInfo << "Setup assembler... " << flush;
    gsConstantFunction<> gammafunction (gamma, mp.geoDim());
    gsPiecewiseFunction<real_t> f_pw = makeUniformPiecewiseFunction(*f,mp.nPatches());
    gsPiecewiseFunction<real_t> gammafunction_pw = makeUniformPiecewiseFunction(gammafunction,mp.nPatches());
    gsPoissonHeterogeneousPde<real_t> pde(mp, bc, f_pw, gammafunction_pw );
    gsPoissonHeterogeneousAssembler<real_t> assm(
        pde,
        bases.back(),
        (dirichlet::strategy)assemblerOptions.getInt("DirichletStrategy"),
        (iFace::strategy)assemblerOptions.getInt("InterfaceStrategy")
    );
    //gsGeneralizedPoissonAssembler<real_t> assm(mp, bases.back(), bc, *f, 0, assemblerOptions);
    gsInfo << "done." << endl;

    gsInfo << "Assembling stiffness matrix and rhs... " << flush;
    assm.assemble();
    gsSparseMatrix<real_t> Kfine = assm.matrix();
    gsMatrix<real_t> rhs = assm.rhs();
    gsInfo << "done, " << Kfine.rows() << "x" << Kfine.cols() << endl;


    std::vector< gsSparseMatrix<> > Mlevels;
    if (smoother == Smoother::MassRichardsonSubspaceCorrectionPiecewiseQuasiExact
        || smoother == Smoother::MassRichardsonSubspaceCorrectionPatchwiseQuasiExact
    )
    {
        gsInfo << "Assembling mass matrix... " << flush;
        const index_t sz = transferMatrices.size() + 1;
        Mlevels.resize( sz );

        // Assemble mass
        gsGenericAssembler<real_t> gassm(mp, bases.back(), assemblerOptions, &bc);
        Mlevels[ sz-1 ] = gassm.assembleMass();
        gsInfo << "done, " << Mlevels[ sz-1 ].rows() << "x" << Mlevels[ sz-1 ].cols() << endl;

        gsInfo << "Restrict mass matrix to levels... " << flush;
        for (index_t i = sz-2; i >= 0; --i)
            Mlevels[i] = transferMatrices[i].transpose() * Mlevels[i+1] * transferMatrices[i];
        gsInfo << "done." << endl;
    }

    /*************************************************************************/

    const double timeAssembling = time.stop(); time.restart();

    // set up the multigrid solver
    gsMultiGridOp<> mg(give(Kfine), transferMatrices);

    // Determine coarse solve time:
    gsTimedOp<>::Ptr coarseSolver = gsTimedOp<>::make("CoarseSolver",mg.coarseSolver(),false);
    mg.setCoarseSolver(coarseSolver);

    mg.setNumPreSmooth( numPreSmooth );
    mg.setNumPostSmooth( numPostSmooth );
    mg.setNumCycles( cycles );
    mg.setCoarseGridCorrectionDamping( coarseDamping );

    if ( notGalerkin )
    {
        gsInfo << "Reassemble on the coarser levels... " << flush;
        const index_t sz = transferMatrices.size(); // not +1 as the finest grid is already perefct!
        for (index_t i=0; i<sz; ++i)
        {
            gsInfo << "level " << i << "... " << flush;
            gsPoissonHeterogeneousPde<real_t> pde0(mp, bc, f_pw, gammafunction_pw );
            gsPoissonHeterogeneousAssembler<real_t> assm0(
                pde0,
                bases[i],
                (dirichlet::strategy)assemblerOptions.getInt("DirichletStrategy"),
                (iFace::strategy)assemblerOptions.getInt("InterfaceStrategy")
            );
            assm0.assemble();
            gsSparseMatrix<real_t> K = assm0.matrix();
            gsInfo << "done, " << K.rows() << "x" << K.cols() << ". ";
            mg.setUnderlyingOp( i, makeMatrixOp(K.moveToPtr()) );
        }

        gsInfo << "done." << endl;
    }

    GISMO_ASSERT( outerDamping == 1 || smoother != Smoother::GaussSeidel, "GaussSeidel does not support --outerdamping" );

    gsInfo << "Constructing smoothers... " << flush;
    for (int i = mg.numLevels() == 1 ? 0 : 1; i < mg.numLevels(); ++i)
    {
        switch( smoother ) {
            case Smoother::Richardson:                            mg.setSmoother(i, makeRichardsonOp(mg.matrix(i),damping*outerDamping)); break;
            case Smoother::Jacobi:                                mg.setSmoother(i, makeJacobiOp(mg.matrix(i),damping*outerDamping)); break;
            case Smoother::GaussSeidel:                           mg.setSmoother(i, makeGaussSeidelOp(mg.matrix(i))); break;
            case Smoother::MassRichardsonSubspaceCorrectionPiecewise:
            {
                mg.setSmoother(
                    i,
                    gsAdditiveSmoother::make(
                        mg.underlyingOp(i),
                        setupPiecewisePreconditioner(
                            mg.matrix(i),
                            makeSubspaceCorrectedMassSmootherOperatorsDirichlet(bases[i],damping),
                            bases[i],
                            bc,
                            assemblerOptions
                        ),
                        outerDamping
                    )
                );
                break;
            }
            case Smoother::MassRichardsonSubspaceCorrectionPatchwise:
            {
                GISMO_ENSURE ( false, "Not implemented" );
                mg.setSmoother(
                    i,
                    gsAdditiveSmoother::make(
                        mg.underlyingOp(i),
                        getPatchwiseTransfers(bases[i], bc, assemblerOptions),
                        makeSubspaceCorrectedMassSmootherOperators(bases[i],damping,bc),
                        outerDamping
                    )
                );
                break;
            }
            case Smoother::MassRichardsonSubspaceCorrectionPiecewiseQuasiExact:
            {
                GISMO_ENSURE ( boundaryCondition == std::string("d"), "Not implemented" );
                const real_t h = bases[i][0].getMinCellLength();
                mg.setSmoother(
                    i,
                    gsAdditiveSmoother::make(
                        mg.underlyingOp(i),
                        setupPiecewisePreconditioner(
                            mg.matrix(i),
                            makePiecewiseExactSolvers(
                                (1/damping) * mg.matrix(i) + (12/(h*h*damping)) * Mlevels[i],
                                bases[i],
                                // bc,
                                assemblerOptions
                            ),
                            bases[i],
                            bc,
                            assemblerOptions
                        ),
                        outerDamping
                    )
                );
                break;
            }
            case Smoother::MassRichardsonSubspaceCorrectionPatchwiseQuasiExact:
            {
                const real_t h = bases[i][0].getMinCellLength();
                std::vector< gsSparseMatrix<> > patchwise_transfers = getPatchwiseTransfers(bases[i], bc, assemblerOptions);
                std::vector< gsLinearOperator<>::Ptr > patchwise_solvers = getLocalExactSolvers<real_t>(
                            (1/damping) * mg.matrix(i) + (12/(h*h*damping)) * Mlevels[i],
                            patchwise_transfers
                        );
                mg.setSmoother(
                    i,
                    gsAdditiveSmoother::make(
                        mg.underlyingOp(i),
                        give(patchwise_transfers),
                        patchwise_solvers,
                        outerDamping
                    )
                );

                break;
            }

        }
    }
    gsInfo << "done." << "\n";

    gsMatrix<> x;
    //x.setRandom( mg.nDofs(), 1 );
    x.setZero( mg.nDofs(), 1 );

    const double timeSetup = time.stop(); time.restart();

    const real_t resNorm0 = (rhs - mg.matrix() * x).norm();
    gsInfo << "Residual norm:     " << resNorm0 << "\n";
    real_t resNorm(0), oldResNorm = resNorm0;

    int numIter = 0;
    real_t minReduction = 1e6;

    gsConjugateGradient<> cg( mg.underlyingOp(), memory::make_shared_not_owned(&mg) );
    cg.setTolerance( tol );

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
        gsInfo << "Residual norm:     " << left << setw(15) << resNorm << "          reduction:  1 / " << setprecision(3) << (oldResNorm/resNorm) << setprecision(6) << "\n";
        minReduction = math::min(minReduction, oldResNorm/resNorm);
        oldResNorm = resNorm;
        ++numIter;
    } while (resNorm / resNorm0 > tol && numIter < maxIter && gsIsfinite(resNorm));

    const double timeSolve = time.stop();

    if (resNorm / resNorm0 > tol || ! gsIsfinite(resNorm))
        gsInfo << "Did not converge.\n";
    else
        gsInfo << "Converged in " << numIter << " iterations.\n";
    gsInfo << "Average convergence factor:  1 / " << setprecision(3) << math::pow(resNorm0 / resNorm, 1.0 / numIter) << setprecision(6) << "\n";
    gsInfo << "Worst   convergence factor:  1 / " << setprecision(3) << minReduction << setprecision(6) << "\n";
    gsInfo << "\n";

    gsInfo << "Assembling time: "; formatTime(gsInfo, timeAssembling);               gsInfo << "\n";
    gsInfo << "MG time:         "; formatTime(gsInfo, timeSetup+timeSolve);          gsInfo << "\n";
    gsInfo << " setup:          "; formatTime(gsInfo, timeSetup);                    gsInfo << "\n";
    gsInfo << " solving:        "; formatTime(gsInfo, timeSolve);
    gsInfo << "         (avg. "; formatTime(gsInfo, timeSolve/numIter);          gsInfo << " per iteration)" << "\n";
    gsInfo << "  coarse solver: "; formatTime(gsInfo, coarseSolver->getTime());
    gsInfo << "\n";
    gsInfo << "Total time:      "; formatTime(gsInfo, timeAssembling+timeSetup+timeSolve);                    gsInfo << "\n";
    gsInfo << "\n";

    real_t l2_error = 0;
    if (monitorL2)
    {
        gsField<> sol = assm.constructSolution(x);
        l2_error = computeL2Distance(sol, *g, false, 3*mg.nDofs());
        gsInfo << "L2 error: " << l2_error << "\n";
    }

    real_t dg_error = 0;
    if (monitorDG)
    {
        gsField<> sol = assm.constructSolution(x);
        dg_error = igaFieldDGDistance(sol, *g, false);
        gsInfo << "DG error: " << dg_error << "\n";
    }

    if (writeLog)
    {
        fstream log("out.txt", fstream::out | fstream::app);

        log << "dg_multigrid_example" << (useCG ? "_cg" : "" ) << "\t"
            << geometry << "\t"
            << 0 << "\t"
            << cycles << "\t"
            << degree << (fancy ? "_fancy" : "") << "\t"
            << numRefine << "\t"
            << numLevels << "\t"
            << numSplit << "\t"
            << boundaryCondition << "\t"
            << damping << "\t"
            << outerDamping << "\t"
            << (int)smoother << "\t"
            << numPreSmooth << "\t"
            << numPostSmooth << "\t"
            << numIter << "\t"
            << math::pow(resNorm0 / resNorm, 1.0 / numIter) << "\t"
            << std::fixed << std::setprecision(6)
            << timeAssembling << "\t"
            << timeSetup << "\t"
            << timeSolve << "\t"
            << l2_error << "\t"
            << dg_error << "\n";
    }

    bool plot = !true;
    if (plot)
    {
        gsField<> sol = assm.constructSolution(x);
        // Plot solution in Paraview
        gsInfo << "Plotting in Paraview: multigrid.pvd.\n";
        gsWriteParaview<>(sol, "multigrid", 1000);
        gsFileManager::open("multigrid.pvd");
    }

    return (resNorm / resNorm0 > tol) ? 1 : 0;
};
