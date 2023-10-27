/** @file gsMassSmoother.cpp

    @brief Provides Multigrid smoothers based on mass matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, C. Hofreither, J. Sogn
*/

#include <gsMultiGrid/gsMassSmoother.h>
#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsRankOneAssembler.h>
#include <gsSolver/gsLowRankCorrectedOp.h>
#include <gsSolver/gsKronecker.h>
#include <gsSolver/gsAdditiveSchwarzOp.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsTensor/gsTensorTools.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsSolver/gsPowerIteration.h>

namespace gismo {

template< index_t d >
gsLinearOperator<>::Ptr makeMassSmootherTensor(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    const gsTensorBSplineBasis<d,real_t> * tb = dynamic_cast<const gsTensorBSplineBasis<d,real_t>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This smoother only works for gsTensorBSplineBasis.");

    const real_t h = tb->getMinCellLength();

    gsSparseMatrix<> temp;
    std::vector< gsSparseMatrix<> > M(d);
    for( index_t i=0; i<d; ++i )
        assembleParameterMass(tb->component(i), M[i]);

    for( index_t i=0; i<d; ++i )
        handleDirichletConditions(M[i],bc,1+2*i,2+2*i);
        // boundary::west = 1, east = 2, south = 3, north = 4, front = 5, back = 6, stime = 7, etime = 8 (cf. gsCore/gsBoundary.h)
    M[0] *= (1.0 / (h*h*damping));

    std::vector<gsLinearOperator<>::Ptr> solvers(d);
    for( index_t i=0; i<d; ++i )
        solvers[d-1-i] = makeSparseCholeskySolver( M[i] );

    return gsKroneckerOp<>::make(solvers);
}

gsLinearOperator<>::Ptr makeMassSmootherOperator(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    if (basis.dim() == 1)
        return makeMassSmootherTensor<1>(basis, damping, bc);
    else if (basis.dim() == 2)
        return makeMassSmootherTensor<2>(basis, damping, bc);
    else if (basis.dim() == 3)
        return makeMassSmootherTensor<3>(basis, damping, bc);
    else if (basis.dim() == 4)
        return makeMassSmootherTensor<4>(basis, damping, bc);
    else
    {
        GISMO_ERROR ("This smoother is instanciated only for up to 4 dimensions.");
    }
}

void makeKtilde(const gsSparseMatrix<>& K, int p, gsSparseMatrix<>& Ktilde)
{
    const int N = K.rows();
    const int NI = N - 2*p;
    GISMO_ENSURE( NI > 0,"Space too small, no interior dofs!");
    //gsInfo << "Making 1D Schur complement, N = " << N << ", p = " << p << ", NI = " << NI << "\n";

    gsMatrix<> K_gg(2*p, 2*p);
    K_gg.block(0, 0, p, p) = K.block(0, 0, p, p);
    K_gg.block(0, p, p, p) = K.block(0, N-p, p, p);
    K_gg.block(p, 0, p, p) = K.block(N-p, 0, p, p);
    K_gg.block(p, p, p, p) = K.block(N-p, N-p, p, p);

    gsMatrix<> K_ii = K.block(p, p, NI, NI);

    gsMatrix<> K_gi(2*p, NI);
    K_gi.block(0, 0, p, NI) = K.block(0, p, p, NI);
    K_gi.block(p, 0, p, NI) = K.block(N-p, p, p, NI);

    Eigen::PartialPivLU< Eigen::Matrix<real_t, Dynamic, Dynamic> > solver(K_ii);

    gsMatrix<> schur = K_gg - K_gi * solver.solve(K_gi.transpose()).eval();

    gsSparseEntries<real_t> entries;
    entries.reserve(4*p*p);

    // put the four blocks of the Schur complement into the four corners
    // of the correction matrix
    for (int i = 0; i < p; ++i)
        for (int j = 0; j < p; ++j)
        {
            entries.add(i    ,     j, schur(i  ,   j));
            entries.add(N-p+i,     j, schur(p+i,   j));
            entries.add(    i, N-p+j, schur(i  , p+j));
            entries.add(N-p+i, N-p+j, schur(p+i, p+j));
        }

    Ktilde.resize(N, N);
    Ktilde.setFrom(entries);
    Ktilde.makeCompressed();
}

// Supposed to implement arXiv:1512.07091
gsLinearOperator<>::Ptr makeBoundaryCorrectedMassSmoother1D(const gsBSplineBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    const int p = basis.maxDegree();
    const real_t h = basis.getMinCellLength();

    gsSparseMatrix<> M;
    assembleParameterMass(basis, M);

    gsSparseMatrix<> K;
    assembleParameterStiffness(basis, K);

    gsSparseMatrix<> sm;

    makeKtilde(K, p, sm);

    sm += (1.0 / (h*h * damping)) * M;

    handleDirichletConditions(sm,bc,boundary::west,boundary::east);

    return makeSparseCholeskySolver(sm);
}

// Supposed to implement arXiv:1512.07091
gsLinearOperator<>::Ptr makeBoundaryCorrectedMassSmoother2D(const gsTensorBSplineBasis<2,real_t>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{

    GISMO_ASSERT(basis.dim() == 2, "This only works for two dimensions.");

    const int p = basis.maxDegree();
    const real_t h = basis.getMinCellLength();

    gsSparseMatrix<> Mx, My;
    assembleParameterMass(basis.component(0), Mx);
    assembleParameterMass(basis.component(1), My);

    gsSparseMatrix<> Kx, Ky, Ktx, Kty;

    assembleParameterStiffness(basis.component(0), Kx);
    assembleParameterStiffness(basis.component(1), Ky);

    // add mass matrix for -Laplace u + u
    Kx += Mx;
    Ky += My;

    makeKtilde(Kx, p, Ktx);
    makeKtilde(Ky, p, Kty);

    // make 1D smoothers

    //                     Smx = (1/(h*h)) * Mx + Ktx;
    //               SmxDamped = ((h*h)/damping) * Smx;
    gsSparseMatrix<> SmxDamped = (1/damping) * Mx + ((h*h) / damping) * Ktx;
    gsSparseMatrix<>       Smy = (1/(h*h)) * My + Kty;

    // setup transfer-matrices Px and Py
    const index_t Nx = Ktx.rows();
    const index_t Ny = Kty.rows();

    gsSparseMatrix<> Px(Nx,2*p), Py(Ny,2*p);

    gsSparseEntries<real_t> entriesPx, entriesPy;
    entriesPx.reserve(2*p);
    entriesPy.reserve(2*p);

    for (int j = 0; j < p; ++j)
    {
        entriesPx.add(     j,  j,1. );
        entriesPx.add(Nx-p+j,p+j,1. );
        entriesPy.add(     j,  j,1. );
        entriesPy.add(Ny-p+j,p+j,1. );
    }

    Px.resize(Nx, 2*p);
    Px.setFrom(entriesPx);
    Px.makeCompressed();

    Py.resize(Ny, 2*p);
    Py.setFrom(entriesPy);
    Py.makeCompressed();

    // incorporate Dirichlet boundary conditions

    handleDirichletConditions(SmxDamped,bc,boundary::west,boundary::east);
    handleDirichletConditions(Px,bc,boundary::west,boundary::east);
    handleDirichletConditions(Ktx,bc,boundary::west,boundary::east);

    handleDirichletConditions(Smy,bc,boundary::south,boundary::north);
    handleDirichletConditions(Py,bc,boundary::south,boundary::north);
    handleDirichletConditions(Kty,bc,boundary::south,boundary::north);

    // setup low-rank formulation of Ktx(x)Kty = U Q^(-1) V
    gsSparseMatrix<> temp = Ktx * Px;
    gsSparseMatrix<> U = Py.kron(temp);

    temp = Kty * Py;
    gsSparseMatrix<> V = temp.kron(Px);

    const index_t rank = U.cols();
    gsSparseMatrix<> Q(rank,rank);
    Q.setIdentity();
    Q *= (damping/(h*h));

    // this represents the inverse of
    //   SmxDamped(x)Smy - U Q^(-1) V^T
    //   =  ((h*h)/damping) * ( Smx(x)Smy - Ktx(x)Kty )

    return gsLowRankCorrectedOp<>::make(
            gsKroneckerOp<>::make( makeSparseCholeskySolver(Smy), makeSparseCholeskySolver(SmxDamped) ),
            Q,
            U,
            V
    );

}

// Supposed to implement arXiv:1512.07091
gsLinearOperator<>::Ptr makeBoundaryCorrectedMassSmootherOperator(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    // just combine the above two
    if (basis.dim() == 1)
    {
        // only works for B-spline bases!
        const gsBSplineBasis<real_t> * tb = dynamic_cast<const gsBSplineBasis<real_t>*>(&basis);
        if( tb )
            return makeBoundaryCorrectedMassSmoother1D(*tb, damping, bc);
        else
            GISMO_ERROR("This smoother only works for gsBSplineBasis.");
    }
    else if (basis.dim() == 2)
    {
        // only works for tensor product B-spline bases!
        const gsTensorBSplineBasis<2, real_t> * tb = dynamic_cast<const gsTensorBSplineBasis<2,real_t>*>(&basis);
        if( tb )
            return makeBoundaryCorrectedMassSmoother2D(*tb, damping, bc);
        else
            GISMO_ERROR("This smoother only works for gsTensorBSplineBasis.");
    }
    else
    {
        GISMO_ERROR("This smoother only works in 1 or 2 dimensions.");
    }
}



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
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherTensor(const gsGeometry<>& domain, const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc, real_t alpha, bool truncate)
{
    //TODO: use domain!!
    const gsTensorBSplineBasis<d,real_t> * tb = dynamic_cast<const gsTensorBSplineBasis<d,real_t>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This smoother only works for gsTensorBSplineBasis.");

    //const real_t h = tb->getMinCellLength();
    gsMatrix<> grd_sz(d,1);

    gsMatrix<> qp(d,1);
    {
        gsMatrix<> supp = tb->support();
        for (index_t i=0; i<d; ++i)
            qp(i,0) = (supp(i,0)+supp(i,1))/2;
    }

    // setup the basis
    std::vector< gsSparseMatrix<> > B_tilde(d), B_l2compl(d), B_compl(d);
    constructTildeSpaceBasis(basis, bc, B_tilde, B_l2compl); //TODO: want/need domain?

    std::vector< gsSparseMatrix<> > M_compl(d), K_compl(d);
    std::vector< gsLinearOperator<>::Ptr > M_tilde_inv(d);
    for ( index_t i=0; i<d; ++i )
    {
        // assemble
        gsSparseMatrix<> M, K;
        assemble1DMass(domain,qp,i,tb->component(i), M);
        assemble1DStiffness(domain,qp,i,tb->component(i), K);
        // boundary::west = 1, east = 2, south = 3, north = 4, front = 5, back = 6, stime = 7, etime = 8 (cf. gsCore/gsBoundary.h)
        handleDirichletConditions(M,bc,1+2*i,2+2*i);
        handleDirichletConditions(K,bc,1+2*i,2+2*i);

        // transform the complement

        gsLinearOperator<>::Ptr M_inv = makeSparseCholeskySolver(M);
        gsMatrix<> B_compl_dense;
        M_inv->apply( B_l2compl[i], B_compl_dense );
        if( truncate )
        {
            // Truncate to only use the first cc basis functions on the particular side

            B_compl[i].resize(B_compl_dense.rows(),B_compl_dense.cols());
            const index_t cc = B_compl_dense.cols()/2; // half of the cols belong to each side, which are treated seperately below
            index_t rr = 2*tb->component(i).degree();
            const index_t rows = B_compl_dense.rows();
            if ( rr > rows ) rr = rows;
            for ( index_t r=0; r<rr; ++r )
                for ( index_t c=0; c<cc; ++c )
                {
                    B_compl[i](r,c) = B_compl_dense(r,c); // left side
                    B_compl[i](rows-r-1,cc+c) = B_compl_dense(rows-r-1,cc+c); // right side
                }
        }
        else
            B_compl[i] = B_compl_dense.sparseView();

        // setup the matrices and corresponding solvers

        gsSparseMatrix<> M_tilde = B_tilde[i].transpose() * M * B_tilde[i];
        M_tilde_inv[i] = makeSparseCholeskySolver( M_tilde );
        gsSparseMatrix<> K_tilde = B_tilde[i].transpose() * K * B_tilde[i];

        // estimate grid size on physical domain using some kind of estimate for the inverse inequality on the tilde spaces
        /*{
            real_t h = 1;
            const index_t len = M_tilde.rows();
            for ( index_t j=0; j<len; ++j )
            {
                const real_t h_local = M_tilde(i,i) / K_tilde(i,i);
                if (h_local < h) h = h_local;
            }
            grd_sz(i,0) = math::sqrt(6*h);
        }*/
        grd_sz(i,0) = math::sqrt(1/powerIteration<real_t>( makeMatrixOp(K_tilde), M_tilde_inv[i], 10)); //TODO: might be expensive
        gsDebug << "grd_sz(i,0)=" << grd_sz(i,0) << "\n";

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

        index_t numberInteriors = 0; //type = alpha (0,1)

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
            {
                tmp = transfers[j]->kron(transfer);
                transfer.swap( tmp );
            }

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
        {
            real_t scaling = alpha;
            for ( index_t k = d-1; k>=0; --k )
                scaling += 1/(damping*grd_sz(k,0)*grd_sz(k,0));
            correction[0] = gsScaledOp<>::make( correction[0], 1/scaling ); //( alpha + numberInteriors/(damping*h*h) )
        }

        // Setup the bondary correction, like  K(x)M(x)M+M(x)K(x)M+M(x)M(x)K+(alpha + (d-3)/(damping*h*h)) M(x)M(x)M
        // \sigma from the paper equals 1/(damping*h*h) here.
        if ( numberInteriors < d )
        {
            gsSparseMatrix<> bc_matrix;

            {
                //~gsInfo << "Setup the BC contriber pure-M" << std::endl;
                gsSparseMatrix<> s(0,0);
                real_t scaling = alpha;
                for ( index_t k = d-1; k>=0; --k )
                {
                    if ( type & ( 1 << k ) )
                    {
                        if ( s.rows() == 0 )
                            s = M_compl[k];
                        else
                        {
                            tmp = M_compl[k].kron(s);
                            s.swap(tmp);
                        }
                    }
                    else
                    {
                        scaling += 1/(damping*grd_sz(k,0)*grd_sz(k,0));
                    }
                }
                bc_matrix = scaling * s; //( alpha + numberInteriors/(damping*h*h) )
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
                            {
                                tmp = s.kron(*chosenMatrix);
                                s.swap(tmp);
                            }
                        }
                    }
                    bc_matrix += s;
                }
            }

            correction.push_back(makeSparseCholeskySolver(bc_matrix));
        }

        // setup the whole operator
        // the correction is the Kronecker-product of the operators in the vector correction
        result->addSubspace( transfer, gsKroneckerOp<>::make( correction ) );
    }

    return result;

}

// A simplified biharmonic version of the smoother
template< index_t d >
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherTensorBiharmonic(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    const gsTensorBSplineBasis<d,real_t> * tb = dynamic_cast<const gsTensorBSplineBasis<d,real_t>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This smoother only works for gsTensorBSplineBasis.");

    const real_t h = tb->getMinCellLength();



    // setup the basis
    std::vector< gsSparseMatrix<> > B_tilde(d), B_l2compl(d), B_compl(d);
    constructTildeSpaceBasis(basis, bc, B_tilde, B_l2compl, false);

    std::vector< gsSparseMatrix<> > M_compl(d), K_compl(d), Bi_compl(d);
    std::vector< gsLinearOperator<>::Ptr > M_tilde_inv(d);
    for ( index_t i=0; i<d; ++i )
    {
        // assemble
        gsSparseMatrix<> M, K, B;
        assembleParameterMass(tb->component(i), M);
        assembleParameterStiffness(tb->component(i), K);
        assembleParameter2ndDer1D(tb->component(i), B);
        // boundary::west = 1, east = 2, south = 3, north = 4, front = 5, back = 6, stime = 7, etime = 8 (cf. gsCore/gsBoundary.h)
        handleDirichletConditions(M,bc,1+2*i,2+2*i);
        handleDirichletConditions(K,bc,1+2*i,2+2*i);
        handleDirichletConditions(B,bc,1+2*i,2+2*i); //JS2: Only works for 1st kind of DBC.


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
        Bi_compl[i] = B_compl[i].transpose() * B * B_compl[i];
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
            {
                tmp = transfers[j]->kron(transfer);
                transfer.swap( tmp );
            }

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
            correction[0] = gsScaledOp<>::make( correction[0], (damping*h*h*h*h)/numberInteriors  );

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
                        {
                            tmp = M_compl[k].kron(s);
                            s.swap(tmp);
                        }
                    }
                }
                bc_matrix = ( numberInteriors/(damping*h*h*h*h) ) * s;
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
                                chosenMatrix = &(Bi_compl[k]);
                            else
                                chosenMatrix = &(M_compl[k]);

                            if ( s.rows() == 0 )
                                s = *chosenMatrix;
                            else
                            {
                                tmp = s.kron(*chosenMatrix);
                                s.swap(tmp);
                            }
                        }
                    }
                    bc_matrix += s;

                }
            }

            correction.push_back(makeSparseCholeskySolver(bc_matrix));
        }

        // setup the whole operator
        // the correction is the Kronecker-product of the operators in the vector correction
        result->addSubspace( transfer, gsKroneckerOp<>::make( correction ) );
    }

    return result;
}

// A simplified biharmonic version of the smoother with 1 rank geometry approimation
template< index_t d >
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherTensorBiharmonic(const gsGeometry<>& domain, const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    const gsTensorBSplineBasis<d,real_t> * tb = dynamic_cast<const gsTensorBSplineBasis<d,real_t>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This smoother only works for gsTensorBSplineBasis.");

    //const real_t h = tb->getMinCellLength();
    gsMatrix<> grd_sz(d,1);

    gsMatrix<> qp(d,1);
    {
        gsMatrix<> supp = tb->support();
        for (index_t i=0; i<d; ++i)
            qp(i,0) = (supp(i,0)+supp(i,1))/2;
    }


    // setup the basis
    std::vector< gsSparseMatrix<> > B_tilde(d), B_l2compl(d), B_compl(d);
    constructTildeSpaceBasis(basis, bc, B_tilde, B_l2compl, false);

    std::vector< gsSparseMatrix<> > M_compl(d), K_compl(d), Bi_compl(d);
    std::vector< gsLinearOperator<>::Ptr > M_tilde_inv(d);

    for ( index_t i=0; i<d; ++i )
    {
        // assemble
        gsSparseMatrix<> M, K, B;
        assemble1DMass(domain, qp, i, tb->component(i), M);
        assemble1DStiffness(domain, qp, i, tb->component(i), K);
        assemble1D2ndDer(domain, qp, i, tb->component(i), B);

        // boundary::west = 1, east = 2, south = 3, north = 4, front = 5, back = 6, stime = 7, etime = 8 (cf. gsCore/gsBoundary.h)
        handleDirichletConditions(M,bc,1+2*i,2+2*i);
        handleDirichletConditions(K,bc,1+2*i,2+2*i);
        handleDirichletConditions(B,bc,1+2*i,2+2*i); //JS2: Only works for 1st kind of DBC.

        // transform the complement

        gsLinearOperator<>::Ptr M_inv = makeSparseCholeskySolver(M);
        gsMatrix<> B_compl_dense;
        M_inv->apply( B_l2compl[i], B_compl_dense );

        B_compl[i] = B_compl_dense.sparseView();

        // setup the matrices and corresponding solvers

        gsSparseMatrix<> M_tilde = B_tilde[i].transpose() * M * B_tilde[i];
        M_tilde_inv[i] = makeSparseCholeskySolver( M_tilde );
        M_compl[i]  = B_compl[i].transpose() * M * B_compl[i];
        K_compl[i]  = B_compl[i].transpose() * K * B_compl[i];
        Bi_compl[i] = B_compl[i].transpose() * B * B_compl[i];

        //gsSparseMatrix<> K_tilde = B_tilde[i].transpose() * K * B_tilde[i];
        gsSparseMatrix<> Bi_tilde = B_tilde[i].transpose() * B * B_tilde[i];

        //grd_sz(i,0) = math::sqrt(1/powerIteration( makeMatrixOp(K_tilde), M_tilde_inv[i], 10)); //TODO: might be expensive
        grd_sz(i,0) = math::sqrt(math::sqrt(1/powerIteration<real_t>( makeMatrixOp(Bi_tilde), M_tilde_inv[i], 10)));
        //grd_sz(i,0) = tb->getMinCellLength();//

        //gsInfo << "\n h: " << tb->getMinCellLength() << " " << grd_sz(i,0) << "  " << math::sqrt(math::sqrt(1/powerIteration( makeMatrixOp(Bi_tilde), M_tilde_inv[i], 10))) <<" " << i <<"\n";

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
            {
                tmp = transfers[j]->kron(transfer);
                transfer.swap( tmp );
            }

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
        {
            real_t scaling = 0;
            for ( index_t k = d-1; k>=0; --k )
                scaling += 1/(damping*grd_sz(k,0)*grd_sz(k,0)*grd_sz(k,0)*grd_sz(k,0));
            correction[0] = gsScaledOp<>::make( correction[0], 1/scaling );
        }

        if ( numberInteriors < d )
        {

            gsSparseMatrix<> bc_matrix;

            {
                //gsInfo << "Setup the BC contriber pure-M" << std::endl;
                gsSparseMatrix<> s(0,0);
                real_t scaling = 0;

                for ( index_t k = d-1; k>=0; --k )
                {
                    if ( type & ( 1 << k ) )
                    {
                        if ( s.rows() == 0 )
                            s = M_compl[k];

                        else
                        {
                            tmp = M_compl[k].kron(s);
                            s.swap(tmp);
                        }
                    }
                    else
                    {
                        scaling += 1/(damping*grd_sz(k,0)*grd_sz(k,0)*grd_sz(k,0)*grd_sz(k,0));
                    }
                }
                bc_matrix = scaling * s;
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
                                chosenMatrix = &(Bi_compl[k]);
                            else
                                chosenMatrix = &(M_compl[k]);

                            if ( s.rows() == 0 )
                                s = *chosenMatrix;
                            else
                            {
                                tmp = s.kron(*chosenMatrix);
                                s.swap(tmp);
                            }
                        }
                    }
                    bc_matrix += s;

                }
            }
            correction.push_back(makeSparseCholeskySolver(bc_matrix));
        }
        // setup the whole operator
        // the correction is the Kronecker-product of the operators in the vector correction
        result->addSubspace( transfer, gsKroneckerOp<>::make( correction ) );
    }

    return result;
}


// A full biharmonic version of the smoother
template< index_t d >
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherTensorFullBiharmonic(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    GISMO_ASSERT(d==2 || d==3, "Only implemented for dimension 2 and 3");

    real_t damping2 = math::sqrt(damping);
    const gsTensorBSplineBasis<d,real_t> * tb = dynamic_cast<const gsTensorBSplineBasis<d,real_t>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This smoother only works for gsTensorBSplineBasis.");

    const real_t h = tb->getMinCellLength();



    // setup the basis
    std::vector< gsSparseMatrix<> > B_tilde(d), B_l2compl(d), B_compl(d);
    constructTildeSpaceBasis(basis, bc, B_tilde, B_l2compl, false);

    std::vector< gsSparseMatrix<> > M_compl(d), K_compl(d), Bi_compl(d);
    std::vector< gsLinearOperator<>::Ptr > M_tilde_inv(d);
    for ( index_t i=0; i<d; ++i )
    {
        // assemble
        gsSparseMatrix<> M, K, B;
        assembleParameterMass(tb->component(i), M);
        assembleParameterStiffness(tb->component(i), K);
        assembleParameter2ndDer1D(tb->component(i), B);
        // boundary::west = 1, east = 2, south = 3, north = 4, front = 5, back = 6, stime = 7, etime = 8 (cf. gsCore/gsBoundary.h)
        handleDirichletConditions(M,bc,1+2*i,2+2*i);
        handleDirichletConditions(K,bc,1+2*i,2+2*i);
        handleDirichletConditions(B,bc,1+2*i,2+2*i); //JS2: Only works for 1st kind of DBC.


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
        Bi_compl[i] = B_compl[i].transpose() * B * B_compl[i];
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
            {
                tmp = transfers[j]->kron(transfer);
                transfer.swap( tmp );
            }

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
            correction[0] = gsScaledOp<>::make( correction[0], (damping*h*h*h*h)/numberInteriors*numberInteriors  );


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
                        {
                            M_compl[k] = s.kron(tmp);
                            s.swap(tmp);
                        }
                    }
                }
                bc_matrix = ( numberInteriors*numberInteriors/(damping*h*h*h*h) ) * s;
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
                                chosenMatrix = &(Bi_compl[k]);
                            else
                                chosenMatrix = &(M_compl[k]);

                            if ( s.rows() == 0 )
                                s = *chosenMatrix;
                            else
                            {
                                tmp = s.kron(*chosenMatrix);
                                s.swap(tmp);
                            }
                        }
                    }
                    bc_matrix += s;
                }
            }
            if (numberInteriors != d)
            {
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
                                {
                                    tmp = s.kron(*chosenMatrix);
                                    s.swap(tmp);
                                }
                            }
                        }
                        bc_matrix += ( 2*numberInteriors/(damping2*h*h) ) * s;
                    }
                }
            }
            else
            {
                for ( index_t j = (d*(d-1))/2-1; j<=0; --j )
                {
                    gsSparseMatrix<> s(0,0);
                    for ( index_t k = d-1; k>=0; --k )
                    {
                        if ( type & ( 1 << k ) )
                        {
                            gsSparseMatrix<>* chosenMatrix;
                            if ( d == 2 )
                                chosenMatrix = &(K_compl[k]);
                            else if ( j == k )
                                chosenMatrix = &(M_compl[k]);
                            else
                                chosenMatrix = &(K_compl[k]);

                            if ( s.rows() == 0 )
                                s = *chosenMatrix;
                            else
                            {
                                tmp = s.kron(*chosenMatrix);
                                s.swap(tmp);
                            }
                        }
                    }
                    bc_matrix +=   2*s;
                }
            }

            correction.push_back(makeSparseCholeskySolver(bc_matrix));
        }

        // setup the whole operator
        // the correction is the Kronecker-product of the operators in the vector correction
        result->addSubspace( transfer, gsKroneckerOp<>::make( correction ) );
    }

    return result;
}

gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherOperator(const gsGeometry<>& domain, const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc, real_t alpha, bool truncate)
{
    if (basis.dim() == 1)
        return makeSubspaceCorrectedMassSmootherTensor<1>(domain, basis, damping, bc, alpha, truncate);
    else if (basis.dim() == 2)
        return makeSubspaceCorrectedMassSmootherTensor<2>(domain, basis, damping, bc, alpha, truncate);
    else if (basis.dim() == 3)
        return makeSubspaceCorrectedMassSmootherTensor<3>(domain, basis, damping, bc, alpha, truncate);
    else if (basis.dim() == 4)
        return makeSubspaceCorrectedMassSmootherTensor<4>(domain, basis, damping, bc, alpha, truncate);
    else
    {
        GISMO_ERROR ("This smoother is instanciated only for up to 4 dimensions.");
    }
}

/// Supposed to implement an extension of G+S-Report 45/2015
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherOperatorBiharmonic(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    if (basis.dim() == 1)
        return makeSubspaceCorrectedMassSmootherTensorBiharmonic<1>(basis, damping, bc);
    else if (basis.dim() == 2)
        return makeSubspaceCorrectedMassSmootherTensorBiharmonic<2>(basis, damping, bc);
    else if (basis.dim() == 3)
        return makeSubspaceCorrectedMassSmootherTensorBiharmonic<3>(basis, damping, bc);
    else if (basis.dim() == 4)
        return makeSubspaceCorrectedMassSmootherTensorBiharmonic<4>(basis, damping, bc);
    else
    {
        GISMO_ERROR ("This smoother is instanciated only for up to 4 dimensions.");
    }
}

/// Supposed to implement an extension of G+S-Report 45/2015
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherOperatorBiharmonic(const gsGeometry<>& domain, const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    if (basis.dim() == 1)
        return makeSubspaceCorrectedMassSmootherTensorBiharmonic<1>(domain, basis, damping, bc);
    else if (basis.dim() == 2)
        return makeSubspaceCorrectedMassSmootherTensorBiharmonic<2>(domain, basis, damping, bc);
    else if (basis.dim() == 3)
        return makeSubspaceCorrectedMassSmootherTensorBiharmonic<3>(domain, basis, damping, bc);
    else if (basis.dim() == 4)
        return makeSubspaceCorrectedMassSmootherTensorBiharmonic<4>(domain, basis, damping, bc);
    else
    {
        GISMO_ERROR ("This smoother is instanciated only for up to 4 dimensions.");
    }
}

/// Supposed to implement an extension of G+S-Report 45/2015
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherOperatorFullBiharmonic(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc)
{
    if (basis.dim() == 2)
        return makeSubspaceCorrectedMassSmootherTensorFullBiharmonic<2>(basis, damping, bc);
    else if (basis.dim() == 3)
        return makeSubspaceCorrectedMassSmootherTensorFullBiharmonic<3>(basis, damping, bc);
    else
    {
        GISMO_ERROR ("This smoother is instanciated only for 2 and 3 dimensions.");
    }
}

} // end of namespace gismo
