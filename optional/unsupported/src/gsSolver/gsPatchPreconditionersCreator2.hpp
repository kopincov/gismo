/** @file gsPatchPreconditionersCreator2.h

    @brief Provides robust preconditioners for single patch geometries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, C. Hofreither
*/
#pragma once

#include <gsSolver/gsPatchPreconditionersCreator.h>
#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsRankOneAssembler.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsKroneckerOp.h>

namespace gismo
{

namespace
{

template<typename T, index_t d>
typename gsLinearOperator<T>::uPtr fastDiagonalizationWithGeoOp_impl(gsGeometry<T>& geo, gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, T alpha)
{
    gsTensorBSplineBasis<d,T> * tb = dynamic_cast<gsTensorBSplineBasis<d,T>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This solver only works for gsTensorBSplineBasis.");

    typedef typename gsLinearOperator<T>::Ptr OpPtr;

    // Determine overall size first
    index_t sz = 1;
    for ( index_t i=0; i<d; ++i )
    {
        index_t mySize = tb->component(i).size() ;
        // boundary::west = 1, east = 2, south = 3, north = 4, front = 5, back = 6, stime = 7, etime = 8 (cf. gsCore/gsBoundary.h)
        patchSide mywest(0,1+2*i), myeast(0,2+2*i);
        if (bc.getConditionFromSide( mywest )!=NULL && bc.getConditionFromSide( mywest )->type() == condition_type::dirichlet ) mySize -= 1;
        if (bc.getConditionFromSide( myeast )!=NULL && bc.getConditionFromSide( myeast )->type() == condition_type::dirichlet ) mySize -= 1;
        sz *= mySize;
    }

    // Initialize the diagonal with 1
    gsMatrix<T> diag;
    diag.setConstant(sz,1,alpha);// This is the pure-mass part!

    index_t glob = 1; // Indexing value for setting up the Kronecker product

    gsSparseMatrix<T> M, K;
    typedef typename gsMatrix<T>::GenSelfAdjEigenSolver EVSolver;
    typedef typename EVSolver::EigenvectorsType evMatrix;
    typedef typename EVSolver::RealVectorType evVector;
    EVSolver ges;

    gsMatrix<> qp(d,1);
    {
        gsMatrix<> supp = tb->support();
        for (index_t i=0; i<d; ++i)
            qp(i,0) = (supp(i,0)+supp(i,1))/2;
    }

    std::vector<OpPtr> Qop(d);
    std::vector<OpPtr> QTop(d);
    gsMatrix<T> ev;

    // Now, setup the Q's and update the D's
    for ( index_t i=0; i<d; ++i )
    {
        // Assemble
        assemble1DMass(geo,qp,i,tb->component(i), M);
        assemble1DStiffness(geo,qp,i,tb->component(i), K);
        // boundary::west = 1, east = 2, south = 3, north = 4, front = 5, back = 6, stime = 7, etime = 8 (cf. gsCore/gsBoundary.h)
        handleDirichletConditions(M,bc,1+2*i,2+2*i);
        handleDirichletConditions(K,bc,1+2*i,2+2*i);

        // Solve generalized eigenvalue problem
        ges.compute(K, M, Eigen::ComputeEigenvectors);
        // Q^T M Q = I, or M = Q^{-T} Q^{-1}
        // Q^T K Q = D, or K = Q^{-T} D Q^{-1}

        // From the eigenvalues, we setup the matrix D already in an Kroneckerized way.
        const evVector & D = ges.eigenvalues();

        const index_t loc = D.rows();
        const index_t glob2 = sz / loc / glob;

        for ( index_t l=0; l<loc; ++l )
            for ( index_t m=0; m<glob; ++m )
                for ( index_t n=0; n<glob2; ++n )
                    diag( m + l*glob + n*loc*glob, 0 ) += D(l,0);

        glob *= loc;

        // Finally, we store the eigenvectors
        ev.swap(const_cast<evMatrix&>(ges.eigenvectors()));

        // These are the operators representing the eigenvectors
        typename gsMatrixOp< gsMatrix<T> >::Ptr matrOp = makeMatrixOp( ev.moveToPtr() );
        Qop [d-1-i] = matrOp;
        // Here we are safe as long as we do not want to apply QTop after Qop got destroyed.
        QTop[d-1-i] = makeMatrixOp( matrOp->matrix().transpose() );
    }

    for ( index_t l=0; l<sz; ++l )
        diag( l, 0 ) = 1/diag( l, 0 );

    memory::unique_ptr< Eigen::DiagonalMatrix<T,Dynamic> > diag_mat( new Eigen::DiagonalMatrix<T,Dynamic>( give(diag) ) );

    return gsProductOp<T>::make(
        gsKroneckerOp<T>::make(QTop),
        makeMatrixOp(give(diag_mat)),
        gsKroneckerOp<T>::make(Qop)
    );

}

}

template<typename T>
typename gsPatchPreconditionersCreator2<T>::OpUPtr gsPatchPreconditionersCreator2<T>::fastDiagonalizationWithGeoOp(gsGeometry<T>& geo, gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, T alpha)
{
    if (basis.dim() == 1)
        return fastDiagonalizationWithGeoOp_impl<T,1>(geo, basis, bc, alpha);
    else if (basis.dim() == 2)
        return fastDiagonalizationWithGeoOp_impl<T,2>(geo, basis, bc, alpha);
    else if (basis.dim() == 3)
        return fastDiagonalizationWithGeoOp_impl<T,3>(geo, basis, bc, alpha);
    else if (basis.dim() == 4)
        return fastDiagonalizationWithGeoOp_impl<T,4>(geo, basis, bc, alpha);
    else
        GISMO_ERROR ("This solver is instanciated only for up to 4 dimensions.");
}


} // namespace gismo
