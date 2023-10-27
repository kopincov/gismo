/** @file gsMultiGridCourantExample.cpp

    @brief Provides test examples for multigrid algorithms with non-nested grids.

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

#include <gsSolver/gsTwoLevel.h>
#include <gsSolver/gsConjugateGradient.h>
#include <gsCore/gsLinearCombinationOfFunctionsFunction.h>
#include <gsSolver/gsKronecker.h>
#include <gsSolver/gsTimedOp.h>
#include <gsMultiGrid/gsMassSmoother.h>

//#include <gsIO/gsMatrixIO.h>
#include <gsTensor/gsTensorTools.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsGeneralizedPoissonAssembler.h>

namespace Smoother {
    enum type {
        Richardson,
        Jacobi,
        GaussSeidel
    };
}



using namespace std;
using namespace gismo;

/// A locally used geometry
gsTensorBSpline<2>::uPtr BSplineMySquare(short_t deg)
{
    gsTensorBSpline<2>::uPtr res(gsNurbsCreator<>::BSplineSquareDeg(deg));
    res->insertKnot( 0.5, 0 );
    return res;
}

/// A locally used geometry
gsGeometry<>::uPtr approximateQuarterAnnulus(int deg)
{
    gsGeometry<>::uPtr quann = gsNurbsCreator<>::NurbsQuarterAnnulus();

    gsKnotVector<> KV1(0,1, 0, deg+1);        // no interior knots in x direction
    gsKnotVector<> KV2(0,1, 1, deg+1);        // 1 interior knot in y direction

    gsTensorBSplineBasis<2> tbsp (new gsBSplineBasis<>(KV1), new gsBSplineBasis<>(KV2));
    return tbsp.interpolateAtAnchors( quann->eval(tbsp.anchors()) );
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
        cerr << "Invalid boundary condition. Allowed are: dirichlet (d), neumann (n) and mixed (dd, dn, nd, nn, ddd, ddn, dnd, dnn, ndd, ndn, nnd, nnn).\n";
        exit(-1);
    }
}

int main(int argc, char *argv[])
{
    index_t geoIndex = 2;
    string boundaryCondition("n");
    index_t numRefine = 3;
    index_t degree = 2;
    index_t numLevels = 0;
    Smoother::type smoother;
    index_t numPreSmooth = 1;
    index_t numPostSmooth = 1;
    index_t cycles = 1;
    real_t alpha = 1.;
    bool monitorl2 = false;
    bool writeLog = false;
    real_t tol = 1e-8;
    real_t damping = -1.0;
    real_t fineDamping = 0.09;
    real_t coarseDamping = 1.;
    index_t maxIter = 1000;
    string smooth("gs");
    bool tildeEmbed = false;
    bool useCG = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");

    cmd.addInt("g", "geometry",
               "Specification of the geometry", geoIndex);
    cmd.addString("b", "boundary-condition",
               "Boundary condition", boundaryCondition);
    cmd.addInt("r", "uniformRefine",
                "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("p", "degree",
               "Degree of the B-spline discretization space", degree);
    cmd.addInt("l", "levels",
               "Number of levels to use for multigrid iteration", numLevels);
    cmd.addInt("", "presmooth",
               "Number of pre-smoothing steps", numPreSmooth);
    cmd.addInt("", "postsmooth",
               "Number of post-smoothing steps", numPostSmooth);
    cmd.addString("s", "smoother",
                  "Smoothing method", smooth);
    cmd.addReal("", "fine-damping",
                "Damping factor for the smoother on the finest grid", fineDamping);
    cmd.addReal("", "coarse-damping",
                "Damping factor for the coarse-grid-correction (on the finest grid)", coarseDamping);
    cmd.addReal("", "damping",
                "Damping factor for the smoother", damping);
    cmd.addInt("c", "cycles",
               "Number of multi-grid cycles", cycles);
    cmd.addInt("", "maxiter", "Maximum number of iterations", maxIter);
    cmd.addReal("a", "alpha", "alpha in \"- LAPLACE u + alpha u = f\"", alpha);
    cmd.addSwitch("tildeEmbed", "Embed the coarse-grid correction in the tilde space", tildeEmbed);
    cmd.addSwitch("cg", "Use method as preconditioner for cg", useCG);
    cmd.addSwitch("log", "Write results to log file", writeLog);
    cmd.addReal("", "tol",
                       "Tolerance for multigrid solver stopping criterion", tol);
    cmd.addSwitch("monitor-l2", "Monitor the discrete l2 errors over the iteration", monitorl2);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    /******************** Define Geometry ********************/

    gsStopwatch time;    
    
    gsGeometry<>::Ptr geo;
    gsMultiPatch<>::Ptr geomp;
       
    switch (geoIndex)
    {
        case 1: geo = gsNurbsCreator<>::BSplineUnitInterval(static_cast<short_t>(degree)); break;
        case 2: geo = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(degree));    break;
        case 3: geo = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(degree));         break;
        case 4: geo = approximateQuarterAnnulus(static_cast<short_t>(degree));             break;
        case 5: geo = BSplineMySquare(static_cast<short_t>(degree));                       break;
        case 6: {
            //Read 3D geometry
            gsFileData<> fileData("volumes/twistedFlatQuarterAnnulus.xml");
            geomp = fileData.getFirst< gsMultiPatch<> >();
            GISMO_ASSERT( geomp->nPatches() == 1, "Multigrid only works for single-patch domains so far." );
            geo = memory::make_shared_not_owned( const_cast<gsGeometry<>*>(&((*geomp)[0])) );
            break;
        }
        case 7: {
            gsFileData<> fileData("yeti_mp2.xml");
            geomp = fileData.getFirst< gsMultiPatch<> >();
            geo =  memory::make_shared_not_owned( const_cast<gsGeometry<>*>(&((*geomp)[0])) );
            break;
        }
        default: cerr << "Invalid geometry. Allowed are:\n"
            << "1: unit interval\n"
            << "2: unit square\n"
            << "3: unit cube\n"
            << "4: approximate quarter annulus\n"
            << "5: unit square with one additional refinement in x-direction\n"
            << "6: twistedFlatQuarterAnnulus.xml\n"
            << "7: yeti_mp2.xml\n";
        return -1;
    }
    
    gsFunction<>::Ptr f0, g;

    switch (geo->geoDim())
    {
        case 1:
            f0 = memory::make_shared(new gsFunctionExpr<>("(pi^2 ) * sin(pi*x)",1));
            g = memory::make_shared(new gsFunctionExpr<>("sin(pi*x)",1));
            break;
        case 2:
            f0 = memory::make_shared(new gsFunctionExpr<>("(2*pi^2 ) * sin(pi*x) * sin(pi*y)",2));
            g = memory::make_shared(new gsFunctionExpr<>("sin(pi*x) * sin(pi*y)",2));
            break;
        case 3:
            f0 = memory::make_shared(new gsFunctionExpr<>("(3*pi^2 ) * sin(pi*x) * sin(pi*y) * sin(pi*z)",3));
            g = memory::make_shared(new gsFunctionExpr<>("sin(pi*x) * sin(pi*y) * sin(pi*z)",3));
            break;
        default:
            cerr << "Invalid geometry dimension.\n";
            return -1;
    }
    
    gsFunction<>::Ptr f = gsLinearCombinationOfFunctionsFunction<>::make(1,f0,alpha,g);

    if (numRefine < 1)
    {
        cerr << "Number of refinements must be positive.\n"; return -1;
    }
    if (numLevels == 0)
        numLevels = numRefine;
    if (numLevels < 1)
    {
        cerr << "Number of levels must be positive.\n"; return -1;
    }
    if (numRefine - numLevels + 1 < 0)
    {
        cerr << "Not enough refinements for the desired number of levels.\n"; return -1;
    }
    
    if (smooth == "r" || smooth == "richardson")
        smoother = Smoother::Richardson;
    else if (smooth == "j" || smooth == "jacobi")
        smoother = Smoother::Jacobi;
    else if (smooth == "gs" || smooth == "gauss-seidel")
        smoother = Smoother::GaussSeidel;
    else
    {
        cerr << "Unknown smoother \"" << smooth << "\".\n"; 
        cerr << "Allowed are: richardson (r), jacobi (j), and gauss-seidel (gs).\n";
        return -1;
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
        }
    }

    gsInfo << "Source function: " << *f << ".\n" << "\n";
    gsInfo << "Exact solution:  " << *g << ".\n" << "\n";


    // set up boundary conditions

    gsConstantFunction<> zero(0.0, geo->geoDim());
    //gsConstantFunction<> one (1.0, geo->geoDim());

    gsBoundaryConditions<> bc;
    condition_type::type bc_type;
    gsFunction<> * bc_func;
    
    // if only single BC given, use it in all coordinate directions
    if (boundaryCondition.length() == 1)
        boundaryCondition = string(geo->geoDim(), boundaryCondition[0]);

    if( (index_t)boundaryCondition.length() != geo->geoDim() )
        boundaryCondition = "x"; // Let the bcChoose do the work

    if ( boundaryCondition != string(geo->geoDim(), 'n') )
    {
        gsInfo << "**********\n* WARNING: The tensor assembler does not treat the inhomogenous Dirichlet boundary conditions properly.\n**********\n";
    }
    
    bcChoose( boundaryCondition[0], bc_type, bc_func, &*g, &zero );
    bc.addCondition( boundary::west,  bc_type, bc_func );
    bc.addCondition( boundary::east,  bc_type, bc_func );
    if (geo->geoDim() >= 2)
    {
        bcChoose( boundaryCondition[1], bc_type, bc_func, &*g, &zero );
        bc.addCondition( boundary::south, bc_type, bc_func );
        bc.addCondition( boundary::north, bc_type, bc_func );
    }
    if (geo->geoDim() >= 3)
    {
        bcChoose( boundaryCondition[2], bc_type, bc_func, &*g, &zero );
        bc.addCondition( boundary::front, bc_type, bc_func );
        bc.addCondition( boundary::back,  bc_type, bc_func );
    }

    gsBasis<>::uPtr tbasis = geo->basis().clone();
    gsBasis<>::uPtr tbasisCourant = tbasis->clone();
    tbasisCourant->degreeDecrease(tbasisCourant->totalDegree()-1);
    if(geoIndex == 6 || geoIndex == 7)
        tbasis->degreeIncrease(degree - tbasis->totalDegree() );

    
    gsInfo << "Original degree " << tbasis->totalDegree() << " reduced to degree " << tbasisCourant->totalDegree() << "\n";
    
    vector< gsMultiBasis<> > basesCourant;
    vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;

    gsOptionList assemblerOptions = gsAssembler<>::defaultOptions();

    // Define coarse discretization space by refining the basis of the geometry
    for (int i = 0; i < numRefine - numLevels + 1; ++i)
    {
        tbasisCourant->uniformRefine();       // refine until coarsest level
    }

    // set up the hierarchy of spaces and transfer matrices between them
    gsGridHierarchy<>::buildByRefinement(gsMultiBasis<>(*tbasisCourant), bc, assemblerOptions, numLevels)
        .moveMultiBasesTo(basesCourant)
        .moveTransferMatricesTo(transferMatrices)
        .clear();

    // Define coarse discretization space by refining the basis of the geometry
    for (int i = 0; i < numRefine; ++i)
    {
        tbasis->uniformRefine(); // refine until FINEST level
    }
    
    gsInfo << "Multigrid";
    gsInfo << " with " << numLevels << " levels using " << smootherName(smoother);
    gsInfo << "(" << numPreSmooth << "," << numPostSmooth << ") smoother" << " and " << cycles << "-cycle.\n";
    gsInfo << "Coarse discretization space: dim=" << tbasis->dim() << " deg=" << tbasis->degree(0) << " dofs=" << tbasis->size() << "\n";
    gsInfo << "Coarse discretization space: dim=" << tbasisCourant->dim() << " deg=" << tbasisCourant->degree(0) << " dofs=" << tbasisCourant->size() << "\n";


    /*************************************************************************/
    
    gsSparseMatrix<real_t> Korig;
    gsMatrix<real_t> rhs;
    if (false)
    {
        gsInfo << "Assemble original stiffness... " << flush;
        assembleGeneralizedParameterStiffnessForTensorProductSpace<>(*tbasis, bc, (real_t)1., alpha, Korig);
        assembleParameterMomentsForTensorProduct(*tbasis, bc, *f, rhs);
        gsInfo << "done, " << Korig.rows() << " dofs." << endl;
    }
    else
    {
        gsInfo << "Setup gsGeneralizedPoissonAssembler... " << flush;
        gsMultiPatch<> mp(*geo);
        gsMultiBasis<> mb(*tbasis);
        gsGeneralizedPoissonAssembler<real_t> assm(mp, mb, bc, *f, alpha, gsAssembler<>::defaultOptions());
        gsInfo << "done." << endl;

        gsInfo << "Assembling original stiffness... " << flush;
        assm.assemble();
        Korig = assm.matrix();
        rhs = assm.rhs();
        gsInfo << "done, " << Korig.rows() << " dofs." << endl;
        
    }
    
    gsSparseMatrix<real_t> Kfine;
    
    if (false)
    {
        gsInfo << "Assemble Courant stiffness... " << flush;
        assembleGeneralizedParameterStiffnessForTensorProductSpace<>(basesCourant.back()[0], bc, (real_t)1., alpha, Kfine);
        gsInfo << "done, " << Kfine.rows() << " dofs." << endl;
    }
    else
    {
        gsInfo << "Setup gsGeneralizedPoissonAssembler... " << flush;
        gsMultiPatch<> mp(*geo);
        gsGeneralizedPoissonAssembler<real_t> assm(mp, basesCourant.back(), bc, *f, alpha, gsAssembler<>::defaultOptions());
        gsInfo << "done." << endl;

        gsInfo << "Assembling Courant stiffness... " << flush;
        assm.assemble();
        Kfine = assm.matrix();
        gsInfo << "done, " << Kfine.rows() << " dofs." << endl;
        
    }
    
    gsInfo << "Assemble original mass... " << flush;
    std::vector< gsSparseMatrix<real_t> > Morigs(geo->geoDim());
    for ( index_t i=0; i<geo->geoDim(); ++i )
    {
        assembleParameterMass(tbasis->component(i), Morigs[i]);
        handleDirichletConditions(Morigs[i],bc,1+2*i,2+2*i);
    }
    gsInfo << "done." << endl;
    
    gsInfo << "Assemble rectangular mass..." << flush;
    gsSparseMatrix<real_t> Mrect;
    assembleParameterMassForTensorProductSpace<>( *tbasis, basesCourant[basesCourant.size()-1][0], bc, Mrect);
    gsInfo << "done, " << Mrect.rows() << "x" << Mrect.cols() << " dofs." << endl;
    
    gsSparseMatrix<> B_tilde_full;
    if (tildeEmbed)
    {
        gsInfo << "Setup tilde space basis..." << flush;
        std::vector< gsSparseMatrix<> > B_tilde;
        std::vector< gsSparseMatrix<> > B_tilde2;
        std::vector< gsSparseMatrix<> > B_l2compl;
        constructTildeSpaceBasis( *tbasis, bc, B_tilde, B_l2compl ); //TODO: should this return a different ordering?
        gsInfo << "done." << endl;

        gsInfo << "Multiply the mass matrices with the tilde space basis..." << flush;    
        for ( index_t i=0; i<geo->geoDim(); ++i )
            Morigs[i] = B_tilde[i].transpose() * Morigs[i] * B_tilde[i];
        gsInfo << "done." << endl;    

        gsInfo << "Multiply the rectangular mass matrices with the tilde space basis..." << flush;

        B_tilde2.resize(geo->geoDim());
        for ( index_t i=0; i<geo->geoDim(); ++i )
            B_tilde2[i] = B_tilde[geo->geoDim()-1-i];
        
        B_tilde_full = kroneckerProduct( B_tilde2 );
        Mrect = B_tilde_full.transpose() * Mrect;
        //for ( index_t i=0; i<geo->geoDim(); ++i )
        //    Mrect[i] = B_tilde[i].transpose() * Mrect[i];
        gsInfo << "done." << endl;
    }
    
    gsInfo << "Setup mass solver... " << flush;
    std::vector< gsLinearOperator<>::Ptr> massInvs(geo->geoDim());
    for ( index_t i=0; i<geo->geoDim(); ++i )
        massInvs[geo->geoDim()-1-i] = makeSparseCholeskySolver(Morigs[i]);
    
    gsLinearOperator<>::Ptr massInv = gsKroneckerOp<>::make( massInvs ); 
    gsInfo << "done." << endl;
    
    /*************************************************************************/
    
    const double timeStartMgSetup = time.stop();
    
    // set up the multigrid solver
    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make(Kfine, transferMatrices);
    
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
        }
    }
    gsInfo << "done." << "\n";
  

    // TODO: we might be interested in telling gsMultigrid not to smooth on the finest grid level
    //       we can just set "damping=0" there (with Richardson/Jacobi smoother), however it would be unefficient


    gsMatrix<> x;
    //x.setRandom( Korig.rows(), 1 );
    x.setZero( Korig.rows(), 1 );
    
    const double timeStartSolve = time.stop();
    /*****************************************************************************/
    
    gsMatrix<> res = rhs - Korig * x;
    gsMatrix<> update, tmp, tmp2;
    
    const real_t resNorm0 = res.norm();
    gsInfo << "Residual norm:     " << resNorm0 << "\n";
    real_t resNorm(0), oldResNorm = resNorm0;
    real_t eucl_error, old_eucl_error(0);

    int numIter = 0;
    real_t minReduction = 1e6;

    gsMatrix<> exactDiscreteSol;
    if (monitorl2)
    {
        Eigen::SparseLU< gsSparseMatrix<real_t> > directsolver( Korig );
        exactDiscreteSol = directsolver.solve( rhs );
        old_eucl_error = (exactDiscreteSol - x).norm();
        gsInfo << "Euclidean error: " << old_eucl_error << "\n";
    }
    
    gsLinearOperator<>::Ptr fineSm = makeSubspaceCorrectedMassSmootherOperator(*tbasis, fineDamping, bc, alpha );
    //fineSm = makeGaussSeidelOp(Korig);
    //gsLinearOperator<>::Ptr fineSm = makeBoundaryCorrectedMassSmootherOperator(*tbasis, fineDamping, bc );
     
    gsTwoLevel::Ptr tl;
    
    if (tildeEmbed)
    {
        tl = gsTwoLevel::make(Korig,fineSm,Mrect,massInv,B_tilde_full,mg);
    }
    else
    {
        tl = gsTwoLevel::make(Korig,fineSm,Mrect,massInv,mg);
    }
    
    tl->setPreSmooth( numPreSmooth );
    tl->setPostSmooth( numPostSmooth );
    tl->setCycles( cycles );
    tl->setCoarseDamping( coarseDamping );
    
    gsConjugateGradient<> cg( Korig, tl );
    if (useCG)
    {
        gsInfo << "Setup cg..." << flush;
        cg.setCalcEigenvalues(true);
        if (cg.initIteration( rhs, x )) { gsInfo << "Reached goal." << endl; return -1; }
        gsInfo << "done." << endl;
    }
    
    do
    {
        if (useCG)
            cg.step(x);
        else
            tl->step(rhs,x);        
        
       
        // SOME CHECKS
        
        res = rhs - Korig * x;
        resNorm = res.norm();
        gsInfo << "Residual norm:     " << left << setw(15) << resNorm << "          reduction:  1 / " << setprecision(3) << (oldResNorm/resNorm) << setprecision(6) << "\n";
        minReduction = math::min(minReduction, oldResNorm/resNorm);
        oldResNorm = resNorm;

        if (monitorl2)
        {
            eucl_error = (exactDiscreteSol - x).norm();
            gsInfo << "                                                                   |  l2 error:     "
                    << left << setw(15) << eucl_error << "          reduction:  1 / " << setprecision(3) << (old_eucl_error/eucl_error) << setprecision(6) << "\n";
            old_eucl_error = eucl_error;
        }

        ++numIter;
    } while (resNorm / resNorm0 > tol && numIter < maxIter && gsIsfinite(resNorm));

    const double timeTotal = time.stop();
    const double timeAssembling = timeStartMgSetup;
    const double timeSetup = timeStartSolve - timeStartMgSetup;
    const double timeSolve = timeTotal - timeStartSolve;

    
    if (resNorm / resNorm0 > tol || !gsIsfinite(resNorm))
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
    gsInfo << "Total time:      "; formatTime(gsInfo, timeTotal);                    gsInfo << "\n";
    gsInfo << "\n";
    
    if (useCG)
    {
        gsMatrix<real_t> eigs;
        cg.getEigenvalues(eigs);
        gsInfo << "Eigenvalues: " << eigs.transpose() << std::endl;
        real_t max = 0, min = 1.e+10, essMin = 1.e+10;
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
    
    if (monitorl2)
    {
        const gsVector<> f_err = exactDiscreteSol - x;
        const real_t eucl_error_f = f_err.norm();
        gsInfo << "l2 error: " << eucl_error_f << "\n";
        gsInfo << "Discrete energy error: " << f_err.dot(Korig*f_err)  << "\n";
    }
    
    if (writeLog)
    {
        fstream log("out.txt", fstream::out | fstream::app);

        log << "gsMultiGridCourantExample_" << geoIndex << "\t"
            << alpha << "\t"
            << cycles << "\t"
            << degree << "\t"
            << numRefine << "\t"
            << numLevels << "\t"
            << boundaryCondition << "\t"
            << damping << "\t"
            << smooth << "\t"
            << numPreSmooth << "\t"
            << numPostSmooth << "\t"
            << (tildeEmbed ? "Stilde" : "S") << "\t"
            << numIter << "\t"
            << math::pow(resNorm0 / resNorm, 1.0 / numIter) << "\t"
            << timeAssembling << "\t"
            << (timeSetup+timeSolve) /*MG time*/;
            
            
            log << "\n";
    }
    
    return (resNorm / resNorm0 > tol) ? 1 : 0;
};
