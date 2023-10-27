/** @file gsRankOneAssembler.hpp

    @brief Provides very easy assemblers for mass and stiffness matrix, living
    on the parameter domain. Used in multigrid testing

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <gsAssembler/gsRankOneAssembler.h>
#include <gsSolver/gsKronecker.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsGeometryEvaluator.h>
#include <gsCore/gsFuncData.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

template <typename T>
void localToGlobal(const gsMatrix<T>& localMat, const gsMatrix<index_t>& localDofs, gsSparseMatrix<T>& globalMat)
{
    const int numActive = localDofs.rows();

    for (index_t i = 0; i < numActive; ++i)
    {
        const int ii = localDofs(i,0);
        for (index_t j = 0; j < numActive; ++j)
        {
            const int jj = localDofs(j,0);
            globalMat.coeffRef(ii, jj) += localMat(i, j);
        }
    }
}

template <typename T>
void assemble1DMass(const gsGeometry<T>& domain, const gsMatrix<T>& qp, index_t direction, const gsBasis<T>& basis, gsSparseMatrix<T>& M)
{
    GISMO_ASSERT ( basis.dim() == 1, "This is for 1D only." );

    const int n = basis.size();
    const int dim = qp.size();

    M.resize(n, n);
    M.reserve( gsVector<int>::Constant(n, 2 * basis.degree(0) + 1) );

    gsMatrix<T> localMass;

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Set the number of integration points for each element
    gsVector<index_t> numQuadNodes( 1 );
    numQuadNodes[0] = basis.degree(0) + 1;

    gsGaussRule<T> quRule( numQuadNodes );
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsFuncData<T> bdata(NEED_VALUE|NEED_ACTIVE);

    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_MEASURE, domain));

    gsMatrix<T> quadratureAt;

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
        basis.compute(quNodes, bdata);
        const index_t numActive = bdata.actives.rows();

        quadratureAt.resize(dim,numQuadNodes[0]);
        for (int i = 0; i<dim; ++i)
        {
            if (i == direction)
                for (int j = 0; j<numQuadNodes[0]; ++j)
                    quadratureAt(i,j) = quNodes(j);
            else
                for (int j = 0; j<numQuadNodes[0]; ++j)
                    quadratureAt(i,j) = qp(i,0);
        }
        geoEval->evaluateAt(quadratureAt);

        localMass.setZero(numActive, numActive);

        for (index_t k = 0; k < quNodes.cols(); ++k)      // loop over quadrature nodes
        {
            const T weight = quWeights[k] * geoEval->measure(k);
            localMass.noalias() += weight * (bdata.values[0].col(k) * bdata.values[0].col(k).transpose());
        }

        localToGlobal(localMass, bdata.actives, M);
    }

    M.makeCompressed();

    if (direction==0) // should be done exactly once
    {
        geoEval->evaluateAt(qp);
        M /= math::pow(geoEval->measure(0),dim-1);
    }
}

template <typename T>
void assemble1DStiffness(const gsGeometry<T>& domain, const gsMatrix<T>& qp, index_t direction, const gsBasis<T>& basis, gsSparseMatrix<T>& K)
{
    GISMO_ASSERT ( basis.dim() == 1, "This is for 1D only." );
    const int n = basis.size();
    const int dim = qp.size();
    K.resize(n, n);
    K.reserve( gsVector<int>::Constant(n, 2 * basis.degree(0) + 1) );
    gsMatrix<T> trf_grads_k;
    gsMatrix<T> tmp;
    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element
    gsVector<index_t> numQuadNodes( 1 );
    numQuadNodes[0] = basis.degree(0) + 1;

    gsGaussRule<T> quRule( numQuadNodes );
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsFuncData<T> bdata(NEED_VALUE|NEED_DERIV|NEED_ACTIVE);

    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM, domain));
    gsMatrix<T> quadratureAt;
    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
        basis.compute(quNodes, bdata);
        const index_t numActive = bdata.actives.rows();

        quadratureAt.resize(dim,numQuadNodes[0]);
        for (int i = 0; i<dim; ++i)
        {
            if (i == direction)
                for (int j = 0; j<numQuadNodes[0]; ++j)
                    quadratureAt(i,j) = quNodes(j);
            else
                for (int j = 0; j<numQuadNodes[0]; ++j)
                    quadratureAt(i,j) = qp(i,0);
        }
        geoEval->evaluateAt(quadratureAt);
        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);

        for (index_t k = 0; k < quNodes.cols(); ++k)      // loop over quadrature nodes
        {
            const T weight = quWeights[k] * geoEval->measure(k);

            const index_t numGrads = bdata.values[1].rows() / 1;
            const gsAsConstMatrix<T> grads_k(bdata.values[1].col(k).data(), 1, numGrads);

            // Extend the 1D-grads (=derivatives) to real gradients
            gsMatrix<T> ext_grads_k;
            ext_grads_k.setZero(numGrads*dim,1);
            for (index_t i=0; i<numGrads; ++i)
                ext_grads_k(direction+dim*i,0) = grads_k(0,i);
            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval->transformGradients(0, ext_grads_k, trf_grads_k);
            localStiffness.noalias() += weight * (trf_grads_k.transpose() * trf_grads_k);
        }  // end loop Gauss nodes
        localToGlobal(localStiffness, bdata.actives, K);
    } //end loop over all domain elements
    K.makeCompressed();
    if (direction==0) // should be done exactly once
    {
        geoEval->evaluateAt(qp);
        K /= math::pow(geoEval->measure(0),dim-1);
    }

}

template <typename T>
void assemble1D2ndDer(const gsGeometry<T>& domain, const gsMatrix<T>& qp, index_t direction, const gsBasis<T>& basis, gsSparseMatrix<T>& B)
{
    GISMO_ASSERT ( basis.dim() == 1, "This is for 1D only." );
    const int n = basis.size();
    const int dim = qp.size();
    const int num2der = dim + (dim*(dim-1))/2;
    B.resize(n, n);
    B.reserve( gsVector<int>::Constant(n, 2 * basis.degree(0) + 1) );
    gsMatrix<T> trf_secondDerivs_k;
    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element
    gsVector<index_t> numQuadNodes( 1 );
    numQuadNodes[0] = basis.degree(0) + 1;

    gsGaussRule<T> quRule(numQuadNodes);
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsFuncData<T> bdata(NEED_VALUE|NEED_ACTIVE|NEED_DERIV|NEED_2ND_DER);

    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM |NEED_2ND_DER, domain));
    gsMatrix<T> quadratureAt;
    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
        basis.compute(quNodes, bdata);

        const index_t numActive = bdata.actives.rows();

        quadratureAt.resize(dim,numQuadNodes[0]);
        for (int i = 0; i<dim; ++i)
        {
            if (i == direction)
                for (int j = 0; j<numQuadNodes[0]; ++j)
                    quadratureAt(i,j) = quNodes(j);
            else
                for (int j = 0; j<numQuadNodes[0]; ++j)
                    quadratureAt(i,j) = qp(i,0);
        }
        geoEval->evaluateAt(quadratureAt);
        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);

        for (index_t k = 0; k < quNodes.cols(); ++k)      // loop over quadrature nodes
        {
            const T weight = quWeights[k] * geoEval->measure(k);

            const index_t numPoints = bdata.values[1].rows() / 1;
            const gsAsConstMatrix<T> double_derv_k(bdata.values[2].col(k).data(), 1, numPoints);
            const gsAsConstMatrix<T> grads_k(bdata.values[1].col(k).data(), 1, numPoints);

            // Extend the 1D-grads (=derivatives) to real gradients
            gsMatrix<T> ext_double_derv_k;
            ext_double_derv_k.setZero(numPoints*num2der,1);
            gsMatrix<T> ext_grads_k;
            ext_grads_k.setZero(numPoints*dim,1);

            for (index_t i=0; i<numPoints; ++i)
                ext_double_derv_k(direction+num2der*i,0) = double_derv_k(0,i);
            for (index_t i=0; i<numPoints; ++i)
                ext_grads_k(direction+dim*i,0) = grads_k(0,i);
            //
            geoEval->transformLaplaceHgrad(0, ext_grads_k, ext_double_derv_k, trf_secondDerivs_k);
            localStiffness.noalias() += weight * (trf_secondDerivs_k.transpose() * trf_secondDerivs_k);
        }  // end loop Gauss nodes
        localToGlobal(localStiffness, bdata.actives, B);
    } //end loop over all domain elements
    B.makeCompressed();
    if (direction==0) // should be done exactly once
    {
        geoEval->evaluateAt(qp);
        B /= math::pow(geoEval->measure(0),dim-1);
    }

}


template <index_t d, typename T>
void assembleRankOneGeneralizedStiffness(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, T alpha, T beta, gsSparseMatrix<T>& K)
{
    const gsTensorBSplineBasis<d,T> * tb = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This assembler only works for gsTensorBSplineBasis.");

    gsMatrix<T> qp(d,1);
    {
        gsMatrix<T> supp = tb->support();
        for (index_t i=0; i<d; ++i)
            qp(i,0) = (supp(i,0)+supp(i,1))/2;
    }

    gsSparseMatrix<T> M;

    assemble1DMass(domain, qp, 0, tb->component(0), M);
    handleDirichletConditions(M,bc,1+2*0,2+2*0);

    assemble1DStiffness(domain, qp, 0, tb->component(0), K);
    handleDirichletConditions(K,bc,1+2*0,2+2*0);

    for ( index_t i=1; i<d; ++i )
    {
        gsSparseMatrix<T> K0, M0;

        assemble1DMass(domain, qp, i, tb->component(i), M0);
        handleDirichletConditions(M0,bc,1+2*i,2+2*i);

        assemble1DStiffness(domain, qp, i, tb->component(i), K0);
        handleDirichletConditions(K0,bc,1+2*i,2+2*i);

        K  = M0.kron(K);
        K += K0.kron(M);
        M  = M0.kron(M);
    }

    if ( alpha != 1. )
        K *= alpha;

    if ( beta != 0. )
        K += beta * M;

}

template <index_t d, typename T>
void assembleRankOneSimpleBiharmonic(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& B)
{
    const gsTensorBSplineBasis<d,T> * tb = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This assembler only works for gsTensorBSplineBasis.");

    gsMatrix<T> qp(d,1);
    {
        gsMatrix<T> supp = tb->support();
        for (index_t i=0; i<d; ++i)
            qp(i,0) = (supp(i,0)+supp(i,1))/2;
    }

    gsSparseMatrix<T> M;

    assemble1DMass(domain, qp, 0, tb->component(0), M);
    handleDirichletConditions(M,bc,1+2*0,2+2*0);

    assemble1D2ndDer(domain, qp, 0, tb->component(0), B);
    handleDirichletConditions(B,bc,1+2*0,2+2*0);

    for ( index_t i=1; i<d; ++i )
    {
        gsSparseMatrix<T> B0, M0;

        assemble1DMass(domain, qp, i, tb->component(i), M0);
        handleDirichletConditions(M0,bc,1+2*i,2+2*i);

        assemble1D2ndDer(domain, qp, i, tb->component(i), B0);
        handleDirichletConditions(B0,bc,1+2*i,2+2*i);

        B  = M0.kron(B);
        B += B0.kron(M);
        M  = M0.kron(M);

    }
}

template <index_t d, typename T>
void assembleRankOneMass(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& M)
{
    const gsTensorBSplineBasis<d,T> * tb = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This assembler only works for gsTensorBSplineBasis.");

    gsMatrix<T> qp(d,1);
    {
        gsMatrix<T> supp = tb->support();
        for (index_t i=0; i<d; ++i)
            qp(i,0) = (supp(i,0)+supp(i,1))/2;
    //gsDebug << "Support is: " << supp << "\n\n";
    }


    assemble1DMass(domain, qp, 0, tb->component(0), M);
    handleDirichletConditions(M,bc,1+2*0,2+2*0);

    for ( index_t i=1; i<d; ++i )
    {
        gsSparseMatrix<T> M0, tmp;

        assemble1DMass(domain, qp, i, tb->component(i), M0);
        handleDirichletConditions(M0,bc,1+2*i,2+2*i);
        M = M0.kron(M);
    }

}

template <index_t d, typename T>
void assembleRankOneMassInverse(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, typename gsKroneckerOp<T>::Ptr& M)
{
    const gsTensorBSplineBasis<d,T> * tb = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This assembler only works for gsTensorBSplineBasis.");

    gsMatrix<T> qp(d,1);
    {
        gsMatrix<T> supp = tb->support();
        for (index_t i=0; i<d; ++i)
            qp(i,0) = (supp(i,0)+supp(i,1))/2;
    }

    std::vector< typename gsLinearOperator<T>::Ptr > vec(d);

    for ( index_t i=0; i<d; ++i )
    {
        gsSparseMatrix<T> M0;
        assemble1DMass(domain, qp, i, tb->component(i), M0);
        handleDirichletConditions(M0,bc,1+2*i,2+2*i);
        vec[d-i-1] = makeSparseCholeskySolver(M0);
    }

    M = gsKroneckerOp<T>::make( vec );
}

template <typename T>
void assembleRankOneGeneralizedStiffness(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, T alpha, T beta, gsSparseMatrix<T>& K)
{
    if (basis.dim() == 1)
        assembleRankOneGeneralizedStiffness<1,T>(domain, basis, bc, alpha, beta, K);
    else if (basis.dim() == 2)
        assembleRankOneGeneralizedStiffness<2,T>(domain, basis, bc, alpha, beta, K);
    else if (basis.dim() == 3)
        assembleRankOneGeneralizedStiffness<3,T>(domain, basis, bc, alpha, beta, K);
    else if (basis.dim() == 4)
        assembleRankOneGeneralizedStiffness<4,T>(domain, basis, bc, alpha, beta, K);
    else
    {
        GISMO_ERROR ("This assembler is instanciated only for up to 4 dimensions.");
    }
}

template <typename T>
void assembleRankOneSimpleBiharmonic(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& B)
{
    if (basis.dim() == 1)
        assembleRankOneSimpleBiharmonic<1,T>(domain, basis, bc, B);
    else if (basis.dim() == 2)
        assembleRankOneSimpleBiharmonic<2,T>(domain, basis, bc, B);
    else if (basis.dim() == 3)
        assembleRankOneSimpleBiharmonic<3,T>(domain, basis, bc, B);
    else if (basis.dim() == 4)
        assembleRankOneSimpleBiharmonic<4,T>(domain, basis, bc, B);
    else
    {
        GISMO_ERROR ("This assembler is instanciated only for up to 4 dimensions.");
    }
}

template <typename T>
void assembleRankOneMass(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& M)
{
    if (basis.dim() == 1)
        assembleRankOneMass<1,T>(domain, basis, bc, M);
    else if (basis.dim() == 2)
        assembleRankOneMass<2,T>(domain, basis, bc, M);
    else if (basis.dim() == 3)
        assembleRankOneMass<3,T>(domain, basis, bc, M);
    else if (basis.dim() == 4)
        assembleRankOneMass<4,T>(domain, basis, bc, M);
    else
    {
        GISMO_ERROR ("This assembler is instanciated only for up to 4 dimensions.");
    }
}

template <typename T>
void assembleRankOneMassInverse(const gsGeometry<T>& domain, const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, typename gsKroneckerOp<T>::Ptr& M)
{
    if (basis.dim() == 1)
        assembleRankOneMassInverse<1,T>(domain, basis, bc, M);
    else if (basis.dim() == 2)
        assembleRankOneMassInverse<2,T>(domain, basis, bc, M);
    else if (basis.dim() == 3)
        assembleRankOneMassInverse<3,T>(domain, basis, bc, M);
    else if (basis.dim() == 4)
        assembleRankOneMassInverse<4,T>(domain, basis, bc, M);
    else
    {
        GISMO_ERROR ("This assembler is instanciated only for up to 4 dimensions.");
    }
}


}
