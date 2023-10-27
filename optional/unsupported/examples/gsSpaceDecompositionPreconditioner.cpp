/** @file gsSpaceDecompositionPreconditioner.cpp

    @brief Provides test examples for multigrid algorithms with non-nested grids.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, S. Takacs
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsSolver/gsAdditiveSchwarzOp.h>
#include <gsSolver/gsKronecker.h>
#include <gsCore/gsLinearCombinationOfFunctionsFunction.h>
#include <gsSolver/gsTimedOp.h>
#include <gsMultiGrid/gsMassSmoother.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsVisitorPoisson.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorNitsche.h>

namespace Smoother {
    enum type {
        Richardson,
        Jacobi,
        GaussSeidel
    };
}

using namespace std;
using namespace gismo;



/// Constructs a matrix for swapping a tensor product
/// from A (x) B (x) C (x) D (x) E  to   A (x) D (x) C (x) B (x) E,
/// where only the dimensions of those matrices have to be given.
/// Note that also those A, B, etc. could be tensor products (here as dimension just provide the products)
/// Note that also those A, B, etc. could vanish. Then just provide a 1 as dimension, i.e., a scalar.
/// So, literally every thinkable swap is possible.
gsSparseMatrix<> kroneckerSwap( index_t e, index_t d, index_t c, index_t b, index_t a )
{
    const index_t sz = a*b*c*d*e;
    gsSparseMatrix<real_t> result(sz,sz);
    gsSparseEntries<real_t> entries;
    entries.reserve(sz);
    for ( index_t i=0; i<a; ++i )
        for ( index_t j=0; j<b; ++j )
            for ( index_t k=0; k<c; ++k )
               for ( index_t l=0; l<d; ++l )
                   for ( index_t m=0; m<e; ++m )
                       entries.add( i+a*(j+b*(k+c*(l+d*m))), i+a*(l+d*(k+c*(j+b*m))), 1. );

    result.setFrom(entries);
    result.makeCompressed();

    return result;
}

template< index_t d >
gsLinearOperator<>::Ptr makeSpaceDecompositionPreconditionerTensor(gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc, real_t alpha, gsLinearOperator<>::Ptr interior)
{

    gsTensorBSplineBasis<d,real_t> * tb = dynamic_cast<gsTensorBSplineBasis<d,real_t>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This smoother only works for gsTensorBSplineBasis.");

    const real_t h = tb->getMinCellLength();

    
    
    // setup the basis
    std::vector< gsSparseMatrix<> > B_tilde(d), B_l2compl(d), B_compl(d);
    constructTildeSpaceBasis(basis, bc, B_tilde, B_l2compl);
        
    std::vector< gsSparseMatrix<> > M_compl(d), K_compl(d);
    std::vector< gsLinearOperator<>::Ptr > M_tilde_inv(d);
    for ( index_t i=0; i<d; ++i )
    {
        // assemble
        gsSparseMatrix<> M, K;
        assembleParameterMass(tb->component(i), M);
        assembleParameterStiffness(tb->component(i), K);
        // boundary::west = 1, east = 2, south = 3, north = 4, front = 5, back = 6, stime = 7, etime = 8 (cf. gsCore/gsBoundary.h)
        handleDirichletConditions(M,bc,1+2*i,2+2*i);
        handleDirichletConditions(K,bc,1+2*i,2+2*i);

        // transform the complement

        gsLinearOperator<>::Ptr M_inv = makeSparseCholeskySolver(M);
        gsMatrix<> B_compl_dense;
        M_inv->apply( B_l2compl[i], B_compl_dense );
        B_compl[i] = B_compl_dense.sparseView();
        
        // setup the matrices and corresponding solvers

        gsSparseMatrix<> M_tilde = B_tilde[i].transpose() * M * B_tilde[i];
        M_tilde_inv[i] = makeSparseCholeskySolver( M_tilde );
        M_compl[i] = B_compl[i].transpose() * M * B_compl[i];
        K_compl[i] = B_compl[i].transpose() * K * B_compl[i];
    }
        
    // setup the final operator
    gsAdditiveSchwarzOp<>::Ptr result = gsAdditiveSchwarzOp<>::make();
    
    for ( index_t type = 0; type < (1<<d); ++ type )
    {
        std::vector< gsLinearOperator<>::Ptr > correction(0);
        gsSparseMatrix<> transfer, tmp;
        
        std::vector< gsSparseMatrix<>* > transfers(d);
 
        index_t numberInteriors = 0;
        
        // setup the transfer       
        for ( index_t j = 0; j<d; ++ j )
        {
            if ( type & ( 1 << j ) )
                transfers[j] = &(B_compl[j]);
            else
            {
                transfers[j] = &(B_tilde[j]);
                ++numberInteriors;
            }
            
            if ( j == 0 )
                transfer = *(transfers[j]);
            else
                transfer = transfers[j]->kron(transfer);
            
        }
        
        // If the subspace is not present, ignore it.
        if ( transfer.cols() == 0 )
            continue;
          
        // setup the swap, where the boundary part is shifted to the begin
        
        index_t left = 1, current = transfers[d-1]->cols(), right = 1;
        for ( index_t j = 0; j < d-1; ++j )
            left *= transfers[j]->cols();

        for ( index_t j = d-1; j >= 0; --j )
        {
            if ( type & ( 1 << j ) )
            {
                transfer = transfer * kroneckerSwap( right, current, left, 1, 1 );
                if ( j > 0 )
                {
                    left *= current;
                    current = transfers[j-1]->cols();
                    left /= current;
                }
            }
            else
            { 
                if ( j > 0 )
                {
                    right *= current;
                    current = transfers[j-1]->cols();
                    left /= current;
                }
            }

        }

        // setup the interior correction
        for ( index_t j = d-1; j>=0; --j )
        {
            if ( ! ( type & ( 1 << j ) ) )
                correction.push_back( M_tilde_inv[j] );
        }
        
        // If we are in the interior, we have to do the scaling here as there is no boundary correction
        if ( numberInteriors == d )
            correction[0] = gsScaledOp<>::make( correction[0], 1./( alpha + numberInteriors/(damping*h*h) ) );
        
        // Setup the bondary correction, like  K(x)M(x)M+M(x)K(x)M+M(x)M(x)K+(alpha + (d-3)/(damping*h*h)) M(x)M(x)M
        // \sigma from the paper equals 1/(daming*h*h) here.
        if ( numberInteriors < d )
        {
            gsSparseMatrix<> bc_matrix;

            {
                //~gsInfo << "Setup the BC contriber pure-M" << std::endl;
                gsSparseMatrix<> s(0,0);
                for ( index_t k = d-1; k>=0; --k )
                {
                    if ( type & ( 1 << k ) )
                    {
                        if ( s.rows() == 0 )
                            s = M_compl[k];
                        else
                            s = M_compl[k].kron(s);
                    }
                }
                bc_matrix = ( alpha + numberInteriors/(damping*h*h) ) * s;
            }
            
            
            for ( index_t j = d-1; j>=0; --j )
            {
                if ( type & ( 1 << j ) )
                {
                    gsSparseMatrix<> s(0,0);
                    for ( index_t k = d-1; k>=0; --k )
                    {
                        if ( type & ( 1 << k ) )
                        {
                            gsSparseMatrix<>* chosenMatrix;
                            if ( j == k )
                                chosenMatrix = &(K_compl[k]);
                            else
                                chosenMatrix = &(M_compl[k]);
                            
                            if ( s.rows() == 0 )
                                s = *chosenMatrix;
                            else
                                s = s.kron(*chosenMatrix);
                        }
                    }
                    bc_matrix += s;
                }
            }
        
            correction.push_back(makeSparseCholeskySolver(bc_matrix));
        }
               
        // setup the whole operator
        // the correction is the Kronecker-product of the operators in the vector correction
        if( numberInteriors < d )
            result->addSubspace( transfer, gsKroneckerOp<>::make( correction ) );
        else 
            result->addSubspace( transfer, interior ); //TODO: in this case, we just throw away what we have constructed so far...
    }

    return result;
   
}

gsLinearOperator<>::Ptr makeSpaceDecompositionPreconditioner(gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc, real_t alpha, gsLinearOperator<>::Ptr interior)
{
    if (basis.dim() == 1)
        return makeSpaceDecompositionPreconditionerTensor<1>(basis, damping, bc, alpha, interior);
    else if (basis.dim() == 2)
        return makeSpaceDecompositionPreconditionerTensor<2>(basis, damping, bc, alpha, interior);
    else if (basis.dim() == 3)
        return makeSpaceDecompositionPreconditionerTensor<3>(basis, damping, bc, alpha, interior);
    else if (basis.dim() == 4)
        return makeSpaceDecompositionPreconditionerTensor<4>(basis, damping, bc, alpha, interior);
    else
    {
        GISMO_ERROR ("This smoother is instanciated only for up to 4 dimensions.");
    }
}





template <typename T>
int estimateNonzerosPerRow(const gsBasis<T>& basis1, const gsBasis<T>& basis2)
{
    int nnz = 1;
    for (int i = 0; i < basis1.dim(); ++i) // to do: improve
        nnz *= basis1.degree(i) + basis2.degree(i) + 1;
    return nnz;
}


template <typename T>
void localToGlobal(const gsMatrix<T>       & localMat,
                   const gsMatrix<index_t> & localDofs1,
                   const gsMatrix<index_t> & localDofs2,
                   gsSparseMatrix<T>       & globalMat)
{
    const int numActive1 = localDofs1.rows();
    const int numActive2 = localDofs2.rows();

    for (index_t i = 0; i < numActive1; ++i)
    {
        const int ii = localDofs1(i,0);
        for (index_t j = 0; j < numActive2; ++j)
        {
            const int jj = localDofs2(j,0);
            globalMat.coeffRef(ii, jj) += localMat(i, j);
        }
    }
}


/// A locally used geometry
gsTensorBSpline<2>::uPtr BSplineMySquare(short_t deg)
{
    gsTensorBSpline<2>::uPtr res(gsNurbsCreator<>::BSplineSquareDeg(deg));
    res->insertKnot( 0.5, 0 );
    return res;
}

/// A locally used geometry
gsTensorBSpline<2>::uPtr BSplineQuarterAnnulus()
{
    gsKnotVector<> KV(0,1, 0,3);

    gsMatrix<> C(9,2);
    C  << 1,   0,
          1.5, 0,
          2,   0,
          1,   1,
          1.5, 1.5,
          2,   2,
          0,   1,
          0,   1.5,
          0,   2;

    return memory::make_unique(new gsTensorBSpline<2>(KV,KV, give(C)));
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
    real_t fineDamping = 1.;
    real_t coarseDamping = 1.2;
    index_t maxIter = 1000;
    bool useCG = false;
    string smooth("gs");

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a space decomposition preconditioner.");

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
    cmd.addSwitch("cg", "Use conjugate gradient", useCG);
    cmd.addSwitch("log", "Write results to log file", writeLog);
    cmd.addReal("", "tol",
                       "Tolerance for multigrid solver stopping criterion", tol);
    cmd.addSwitch("monitor-l2", "Monitor the discrete l2 errors over the iteration", monitorl2);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    /******************** Define Geometry ********************/

    gsStopwatch time;    
    
    gsGeometry<>::uPtr geo;
    gsGeometry<>::uPtr geoCourant;
       
    switch (geoIndex)
    {
        case 1:
            geo = gsNurbsCreator<>::BSplineUnitInterval(static_cast<short_t>(degree));
            geoCourant = gsNurbsCreator<>::BSplineUnitInterval(static_cast<short_t>(1));
            break;
        case 2:
            geo = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(degree));
            geoCourant = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(1));
            break;
        case 3:
            geo = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(degree));
            geoCourant = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(1));
            break;
        case 4:
            geo = approximateQuarterAnnulus(static_cast<short_t>(degree));
            geoCourant = approximateQuarterAnnulus(static_cast<short_t>(1));
            break;
        default: cerr << "Invalid geometry. Allowed are:\n"
            << "1: unit interval\n"
            << "2: unit square\n"
            << "3: unit cube\n"
            << "4: approximateQuarterAnnulus\n";
        return -1;
    }
    
    gsFunctionExpr<>::Ptr f0, g;

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
    gsBasis<>::uPtr tbasisCourant = geoCourant->basis().clone();
    
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
    
    gsInfo << "Space decomposition preconditioner. Use multigrid";
    gsInfo << " with " << numLevels << " levels using " << smootherName(smoother);
    gsInfo << "(" << numPreSmooth << "," << numPostSmooth << ") smoother" << " and " << cycles << "-cycle.\n";
    gsInfo << "Coarse discretization space: dim=" << tbasis->dim() << " deg=" << tbasis->degree(0) << " dofs=" << tbasis->size() << "\n";
    gsInfo << "Coarse discretization space: dim=" << tbasisCourant->dim() << " deg=" << tbasisCourant->degree(0) << " dofs=" << tbasisCourant->size() << "\n";


    /*************************************************************************/
    
    gsInfo << "Assemble Courant stiffness... " << flush;
    gsSparseMatrix<real_t> Kfine;
    assembleGeneralizedParameterStiffnessForTensorProductSpace<>(basesCourant.back()[0], bc, (real_t)1., alpha, Kfine);
    gsInfo << "done, " << Kfine.rows() << " dofs." << endl;
    
    gsInfo << "Assemble original stiffness... " << flush;
    gsSparseMatrix<real_t> Korig;
    gsMatrix<real_t> rhs;
    assembleGeneralizedParameterStiffnessForTensorProductSpace<>(*tbasis, bc, (real_t)1., alpha, Korig);
    assembleParameterMomentsForTensorProduct(*tbasis, bc, *f, rhs);
    gsInfo << "done, " << Korig.rows() << " dofs." << endl;
    
    gsInfo << "Assemble original mass... " << flush;
    std::vector< gsSparseMatrix<real_t> > Morigs(geo->geoDim());
    for ( index_t i=0; i<geo->geoDim(); ++i )
        assembleParameterMass(tbasis->component(i), Morigs[i]);
    gsInfo << "done." << endl;
    
    gsInfo << "Assemble rectangular mass..." << flush;
    gsSparseMatrix<real_t> Mrect;
    assembleParameterMassForTensorProductSpace<>( *tbasis, basesCourant.back()[0], bc, Mrect);
    gsInfo << "done, " << Mrect.rows() << "x" << Mrect.cols() << " dofs." << endl;

    gsSparseMatrix<> B_tilde_full;
    {
        gsInfo << "Setup tilde space basis..." << flush;
        std::vector< gsSparseMatrix<> > B_tilde;
        std::vector< gsSparseMatrix<> > B_l2compl;
        constructTildeSpaceBasis( *tbasis, bc, B_tilde, B_l2compl );
        gsInfo << "done." << endl;

        gsInfo << "Multiply the mass matrices with the tilde space basis..." << flush;    
        for ( index_t i=0; i<geo->geoDim(); ++i )
            Morigs[i] = B_tilde[i].transpose() * Morigs[i] * B_tilde[i];
        gsInfo << "done." << endl;    

        gsInfo << "Multiply the rectangular mass matrices with the tilde space basis..." << flush;
        B_tilde_full = kroneckerProduct( B_tilde );
        Mrect = B_tilde_full.transpose() * Mrect;
        //for ( index_t i=0; i<geo->geoDim(); ++i )
        //    Mrect[i] = B_tilde[i].transpose() * Mrect[i];
        gsInfo << "done." << endl;
    }
    
    gsInfo << "Setup mass solver... " << flush;
    std::vector< gsLinearOperator<>::Ptr> massInvs(geo->geoDim());
    for ( index_t i=0; i<geo->geoDim(); ++i )
        massInvs[i] = makeSparseCholeskySolver(Morigs[i]);
    
    gsLinearOperator<>::Ptr massInv = gsKroneckerOp<>::make( massInvs ); 
    gsInfo << "done." << endl;
    
    /*************************************************************************/
    
    const real_t timeStartMgSetup = time.stop();
    
    // set up the multigrid solver
    gsMultiGridOp<> mg(Kfine, transferMatrices);
    
    // Determine coarse solve time:
    gsTimedOp<>::Ptr coarseSolver = gsTimedOp<>::make("CoarseSolver",mg.coarseSolver(),false);
    mg.setCoarseSolver(coarseSolver);

    mg.setNumPreSmooth( numPreSmooth );
    mg.setNumPostSmooth( numPostSmooth );
    mg.setNumCycles( cycles );

    gsInfo << "Constructing smoothers... " << flush;
    for (int i = mg.numLevels() == 1 ? 0 : 1; i < mg.numLevels(); ++i)
    {
         switch( smoother ) {
            case Smoother::Richardson:                            mg.setSmoother(i, makeRichardsonOp(mg.matrix(i),damping)); break;
            case Smoother::Jacobi:                                mg.setSmoother(i, makeJacobiOp(mg.matrix(i),damping)); break;
            case Smoother::GaussSeidel:                           mg.setSmoother(i, makeGaussSeidelOp(mg.matrix(i))); break;
        }
    }
    gsInfo << "done." << "\n";
  

    // TODO: we might be interested in telling gsMultigrid not to smooth on the finest grid level
    //       we can just set "damping=0" there (with Richardson/Jacobi smoother), however it would be unefficient


    gsMatrix<> x;
    x.setRandom( Korig.rows(), 1 );
    //x.setZero( Korig.rows(), 1 );
    
    
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
    
    std::vector< gsLinearOperator<>::Ptr > interiorProducts;
    
    interiorProducts.push_back( massInv );
    interiorProducts.push_back( makeMatrixOp(Mrect.transpose()) );
    interiorProducts.push_back( memory::make_shared_not_owned(&mg) );
    interiorProducts.push_back( makeMatrixOp(Mrect) );
    interiorProducts.push_back( massInv );
    
    gsLinearOperator<>::Ptr interior = gsScaledOp<>::make( gsProductOp<>::make(interiorProducts), coarseDamping );
    
    gsLinearOperator<>::Ptr precond = makeSpaceDecompositionPreconditioner(*tbasis, fineDamping, bc, alpha, interior );
    
    gsConjugateGradient<> cg( Korig, precond );

    if (useCG)
        cg.initIteration( rhs, x );
    
    const real_t timeStartSolve = time.stop();
    do
    {
        if (useCG)
            cg.step(x);
        else
        {
            precond->apply(res, update);
            x+= update;
        }
        res = rhs - Korig * x;
                
        
        // SOME CHECKS
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

    const real_t timeTotal = time.stop();
    const real_t timeAssembling = timeStartMgSetup;
    const real_t timeSetup = timeStartSolve - timeStartMgSetup;
    const real_t timeSolve = timeTotal - timeStartSolve;

    
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

        log << "gsSpaceDecompositionPreconditioner_" << geoIndex << "\t"
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
            << fineDamping << "\t"
            << coarseDamping << "\t"
            << "Stilde\t"
            << numIter << "\t"
            << (useCG ? "cg" : "grad" ) << "\t"
            << math::pow(resNorm0 / resNorm, 1.0 / numIter) << "\t"
            << timeAssembling << "\t"
            << (timeSetup+timeSolve) /*MG time*/;
            
            
            log << "\n";
    }

    return (resNorm / resNorm0 > tol) ? 1 : 0;
};
