/** @file gsParameterDomainAssembler.hpp

    @brief Provides very easy assemblers for mass and stiffness matrix, living
    on the parameter domain. Used in multigrid testing

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsSolver/gsKronecker.h>
#include <gsCore/gsBasis.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <typename T>
int estimateNonzerosPerRow(const gsBasis<T>& basis)
{
    int nnz = 1;
    for (int i = 0; i < basis.dim(); ++i) // to do: improve
        nnz *= 2 * basis.degree(i) + 1;
    return nnz;
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
void localToGlobal(const gsMatrix<T>& localMat,
                   const gsMatrix<index_t>& localDofs1,
                   const gsMatrix<index_t>& localDofs2,
                   gsSparseMatrix<T>& globalMat)
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


template <typename T>
void localToGlobal(const gsMatrix<T>& localMat, const gsMatrix<index_t>& localDofs, gsSparseMatrix<T>& globalMat)
{
    localToGlobal(localMat, localDofs, localDofs, globalMat);
}

template <typename T>
void assembleParameterMass(const gsBasis<T>& basis1, const gsBasis<T>& basis2, gsSparseMatrix<T>& M)
{
    const int n1 = basis1.size();
    const int n2 = basis2.size();

    M.resize(n1, n2);

    M.reserve( gsVector<int>::Constant(n1, estimateNonzerosPerRow(basis1,basis2)) );

    gsMatrix<T> localMass;

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt1 = basis1.makeDomainIterator();
    typename gsBasis<T>::domainIter domIt2 = basis2.makeDomainIterator();

    // Set the number of integration points for each element
    gsVector<index_t> numQuadNodes( basis1.dim() );
    for (int i = 0; i < basis1.dim(); ++i)
        numQuadNodes[i] = basis1.degree(i) + 1;

    gsGaussRule<T> quRule1( numQuadNodes );
    gsGaussRule<T> quRule2( numQuadNodes );
    gsMatrix<T> quNodes1, quNodes2;
    gsVector<T> quWeights1, quWeights2;

    gsFuncData<T> bdata1(NEED_VALUE|NEED_ACTIVE);
    gsFuncData<T> bdata2(NEED_VALUE|NEED_ACTIVE);

    // Start iteration over elements
    for (; domIt1->good() && domIt2->good(); domIt1->next() && domIt2->next())
    {
        GISMO_ENSURE ( domIt1->centerPoint() == domIt2->centerPoint(), "Yes, this somehow does not work. Sorry." );

        quRule1.mapTo(domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights1);
        basis1.compute(quNodes1, bdata1);
        const index_t numActive1 = bdata1.actives.rows();

        quRule2.mapTo(domIt2->lowerCorner(), domIt1->upperCorner(), quNodes2, quWeights2);
        basis2.compute(quNodes2, bdata2);
        const index_t numActive2 = bdata2.actives.rows();

        localMass.setZero(numActive1, numActive2);

        for (index_t k = 0; k < quNodes1.cols(); ++k)      // loop over quadrature nodes
        {
            const T weight = quWeights1[k];
            localMass.noalias() += weight * (bdata1.values[0].col(k) * bdata2.values[0].col(k).transpose());
        }

        localToGlobal(localMass, bdata1.actives, bdata2.actives, M);
    }

    M.makeCompressed();
}


template <typename T>
void assembleParameterStiffness(const gsBasis<T>& basis1, const gsBasis<T>& basis2, gsSparseMatrix<T>& K)
{
    const int n1 = basis1.size();
    const int n2 = basis2.size();
    const int dim = basis1.dim();
    GISMO_ASSERT( basis1.dim() == basis2.dim(), "Bases do not fit." );

    K.resize(n1, n2);
    K.reserve( gsVector<int>::Constant(n1, estimateNonzerosPerRow(basis1,basis2)) );

    gsMatrix<T> trf_grads_k1;    // transformed (physical) gradients of all basis functions at one quadrature node
    gsMatrix<T> trf_grads_k2;    // transformed (physical) gradients of all basis functions at one quadrature node
    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt1 = basis1.makeDomainIterator();
    typename gsBasis<T>::domainIter domIt2 = basis2.makeDomainIterator();

    // Set the number of integration points for each element
    gsVector<index_t> numQuadNodes( dim );
    for (int i = 0; i < basis1.dim(); ++i)
        numQuadNodes[i] = basis1.degree(i) + 1;

    gsGaussRule<T> quRule1( numQuadNodes );
    gsGaussRule<T> quRule2( numQuadNodes );
    gsMatrix<T> quNodes1, quNodes2;
    gsVector<T> quWeights1, quWeights2;

    gsFuncData<T> bdata1(NEED_VALUE|NEED_ACTIVE|NEED_DERIV);
    gsFuncData<T> bdata2(NEED_VALUE|NEED_ACTIVE|NEED_DERIV);

    // Start iteration over elements
    for (; domIt1->good() && domIt2->good(); domIt1->next() && domIt2->next())
    {

        GISMO_ENSURE ( domIt1->centerPoint() == domIt2->centerPoint(), "Yes, this somehow does not work. Sorry." );

        quRule1.mapTo(domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights1);
        basis1.compute(quNodes1, bdata1);
        const index_t numActive1 = bdata1.actives.rows();

        quRule2.mapTo(domIt2->lowerCorner(), domIt1->upperCorner(), quNodes2, quWeights2);
        basis2.compute(quNodes2, bdata2);
        const index_t numActive2 = bdata2.actives.rows();

        // initialize local linear system to 0
        localStiffness.setZero(numActive1, numActive2);

        for (index_t k = 0; k < quNodes1.cols(); ++k)      // loop over quadrature nodes
        {
            const T weight = quWeights1[k];

            // compute physical gradients at k as a Dim x NumActive matrix
            // geoEval->transformGradients(k, bdata1.values[1], trf_grads_k);

            const index_t numGrads = bdata1.values[1].rows() / dim;
            GISMO_ASSERT( bdata1.values[1].rows() == bdata2.values[1].rows(), "Bases do not fit." );
            const gsAsConstMatrix<T> grads_k1(bdata1.values[1].col(k).data(), dim, numGrads);
            trf_grads_k1.noalias() = grads_k1;
            const gsAsConstMatrix<T> grads_k2(bdata2.values[1].col(k).data(), dim, numGrads);
            trf_grads_k2.noalias() = grads_k2;

            localStiffness.noalias() += weight * (trf_grads_k1.transpose() * trf_grads_k2);
        }  // end loop Gauss nodes

        localToGlobal(localStiffness, bdata1.actives, bdata2.actives, K);

    } //end loop over all domain elements

    K.makeCompressed();

}

template <typename T>
void assembleParameter2ndDer1D(const gsBasis<T>& basis, gsSparseMatrix<T>& B)
{
    //JS2: See comments to know where to change if dim != 1.
    GISMO_ASSERT(basis.dim() == 1, "This only works for one dimension.");
    const int n = basis.size();
    const int dim = basis.dim();

    B.resize(n, n);
    B.reserve( gsVector<int>::Constant(n, 2*estimateNonzerosPerRow(basis)) );

    gsMatrix<T> noTrf_secondDerivs_k;    // transformed (physical) gradients of all basis functions at one quadrature node
    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Set the number of integration points for each element
    gsVector<index_t> numQuadNodes( dim );
    for (int i = 0; i < basis.dim(); ++i)
        numQuadNodes[i] = basis.degree(i) + 1;

    gsGaussRule<T> quRule(numQuadNodes);
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsFuncData<T> bdata(NEED_VALUE|NEED_ACTIVE|NEED_DERIV|NEED_2ND_DER);

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
        basis.compute(quNodes, bdata);

        const index_t numActive = bdata.actives.rows();

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);

        for (index_t k = 0; k < quNodes.cols(); ++k)      // loop over quadrature nodes
        {
            const T weight = quWeights[k];

            const index_t numGrads = bdata.values[2].rows() / dim;

            const gsAsConstMatrix<T> secondDerivs_k(bdata.values[2].col(k).data(), dim, numGrads); //JS2: dimension = 1 dependent
            noTrf_secondDerivs_k.noalias() = secondDerivs_k; //JS2: dimension = 1 dependent

            localStiffness.noalias() += weight * (noTrf_secondDerivs_k.transpose() * noTrf_secondDerivs_k); //JS2: dimension = 1 dependent
        }  // end loop Gauss nodes

        localToGlobal(localStiffness, bdata.actives, B);

    } //end loop over all domain elements

    B.makeCompressed();
}

template <typename T>
void handleDirichletConditions(gsSparseMatrix<T>& matrix, const gsBoundaryConditions<T>& bc, const boxSide& west, const boxSide& east)
{
    patchSide mywest(0,west), myeast(0,east);

    int i = 0;

    if (bc.getConditionFromSide( mywest )!=NULL && bc.getConditionFromSide( mywest )->type() == condition_type::dirichlet ) i += 1;
    if (bc.getConditionFromSide( myeast )!=NULL && bc.getConditionFromSide( myeast )->type() == condition_type::dirichlet ) i += 2;

    if (i%2 + i/2 >= matrix.rows() || i%2 + i/2 >= matrix.cols())
    {
        gsWarn << "Boundary conditions force matrix to be empty." << std::endl;
        matrix.resize(0,0); return;
    }
    
    switch ( i )
    {
        case 0: return;
        case 1: matrix = matrix.block( 1, 1, matrix.rows()-1, matrix.cols()-1 ); return;
        case 2: matrix = matrix.block( 0, 0, matrix.rows()-1, matrix.cols()-1 ); return;
        case 3: matrix = matrix.block( 1, 1, matrix.rows()-2, matrix.cols()-2 ); return;
    }
}

template <index_t d, typename T>
void assembleGeneralizedParameterStiffnessForTensorProductSpace(const gsBasis<T>& basis1, const gsBasis<T>& basis2, const gsBoundaryConditions<T>& bc, T alpha, T beta, gsSparseMatrix<T>& K)
{
    const gsTensorBSplineBasis<d,T> * tb1 = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis1);
    if( !tb1 )
        GISMO_ERROR ("This assembler only works for gsTensorBSplineBasis.");

    const gsTensorBSplineBasis<d,T> * tb2 = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis2);
    if( !tb2 )
        GISMO_ERROR ("This assembler only works for gsTensorBSplineBasis.");


    gsSparseMatrix<T> M;

    assembleParameterMass(tb1->component(0), tb2->component(0), M);
    handleDirichletConditions(M,bc,1+2*0,2+2*0);

    assembleParameterStiffness(tb1->component(0), tb2->component(0), K);
    handleDirichletConditions(K,bc,1+2*0,2+2*0);

    for ( index_t i=1; i<d; ++i )
    {
        gsSparseMatrix<T> K0, M0;

        assembleParameterMass(tb1->component(i), tb2->component(i), M0);
        handleDirichletConditions(M0,bc,1+2*i,2+2*i);

        assembleParameterStiffness(tb1->component(i), tb2->component(i), K0);
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

template <typename T>
void assembleGeneralizedParameterStiffnessForTensorProductSpace(const gsBasis<T>& basis1, const gsBasis<T>& basis2, const gsBoundaryConditions<T>& bc, T alpha, T beta, gsSparseMatrix<T>& K)
{
    if (basis1.dim() == 1)
        assembleGeneralizedParameterStiffnessForTensorProductSpace<1,T>(basis1, basis2, bc, alpha, beta, K);
    else if (basis1.dim() == 2)
        assembleGeneralizedParameterStiffnessForTensorProductSpace<2,T>(basis1, basis2, bc, alpha, beta, K);
    else if (basis1.dim() == 3)
        assembleGeneralizedParameterStiffnessForTensorProductSpace<3,T>(basis1, basis2, bc, alpha, beta, K);
    else if (basis1.dim() == 4)
        assembleGeneralizedParameterStiffnessForTensorProductSpace<4,T>(basis1, basis2, bc, alpha, beta, K);
    else
    {
        GISMO_ERROR ("This assembler is instanciated only for up to 4 dimensions.");
    }
}

template <index_t d, typename T>
void assembleParameterMassForTensorProductSpace(const gsBasis<T>& basis1, const gsBasis<T>& basis2, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& M)
{
    const gsTensorBSplineBasis<d,T> * tb1 = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis1);
    if( !tb1 )
        GISMO_ERROR ("This assembler only works for gsTensorBSplineBasis.");
    const gsTensorBSplineBasis<d,T> * tb2 = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis2);
    if( !tb2 )
        GISMO_ERROR ("This assembler only works for gsTensorBSplineBasis.");

    assembleParameterMass(tb1->component(0), tb2->component(0), M);
    handleDirichletConditions(M,bc,1+2*0,2+2*0);

    for ( index_t i=1; i<d; ++i )
    {
        gsSparseMatrix<T> M0;        

        assembleParameterMass(tb1->component(i), tb2->component(i), M0);
        handleDirichletConditions(M0,bc,1+2*i,2+2*i);

        M = M0.kron( M );
    }
}

template <typename T>
void assembleParameterMassForTensorProductSpace(const gsBasis<T>& basis1, const gsBasis<T>& basis2, const gsBoundaryConditions<T>& bc, gsSparseMatrix<T>& M)
{
    if (basis1.dim() == 1)
        assembleParameterMassForTensorProductSpace<1,T>(basis1, basis2, bc, M);
    else if (basis1.dim() == 2)
        assembleParameterMassForTensorProductSpace<2,T>(basis1, basis2, bc, M);
    else if (basis1.dim() == 3)
        assembleParameterMassForTensorProductSpace<3,T>(basis1, basis2, bc, M);
    else if (basis1.dim() == 4)
        assembleParameterMassForTensorProductSpace<4,T>(basis1, basis2, bc, M);
    else
    {
        GISMO_ERROR ("This assembler is instanciated only for up to 4 dimensions.");
    }
}


template <typename T>
void localToGlobalMoments(const gsMatrix<T>& localVec, const gsMatrix<index_t>& localDofs, gsMatrix<T>& globalVec)
{
    const int numActive = localDofs.rows();

    for (index_t i = 0; i < numActive; ++i)
    {
        const int ii = localDofs(i,0);
        globalVec.coeffRef(ii, 0) += localVec(i, 0);
    }
}


template <typename T>
void assembleParameterMoments(const gsBasis<T>& basis, gsFunction<T>& f, gsMatrix<T>& fh, index_t a, index_t b)
{
    //GISMO_ERROR("The function assembleParameterMoments is not working. Ask Juergen or Angelos.");
    //return;

    const index_t n = basis.size();
    const index_t d = basis.dim();

    fh.resize(n, 1);
    fh.fill(0.);
    gsMatrix<T> localMoments;
    gsMatrix<T> fvals;
    gsMatrix<T> pos;
    pos.resize(a+d+b,1);
    pos.fill(0.5); // TODO: is this appropriate?

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Set the number of integration points for each element
    gsVector<index_t> numQuadNodes(d);
    for (int i = 0; i < d; ++i)
        numQuadNodes[i] = basis.degree(i) + 1;

    gsGaussRule<T> quRule( numQuadNodes );
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsFuncData<T> bdata(NEED_VALUE|NEED_ACTIVE);

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
        basis.compute(quNodes, bdata);

        const index_t numActive = bdata.actives.rows();

        localMoments.setZero(numActive, 1);

        for (index_t k = 0; k < quNodes.cols(); ++k)      // loop over quadrature nodes
        {
            const T weight = quWeights[k];
            pos.block(a,0,d,1) = quNodes.col(k);
            f.eval_into(pos, fvals);
            localMoments.noalias() += weight * fvals(0,0) * bdata.values[0].col(k);
        }

        localToGlobalMoments(localMoments, bdata.actives, fh);
    }
}

template <typename T>
void handleDirichletConditionsForMoments(gsMatrix<T>& matrix, const gsBoundaryConditions<T>& bc, const boxSide& west, const boxSide& east)
{
    patchSide mywest(0,west), myeast(0,east);

    int i = 0;

    if ( bc.getConditionFromSide( mywest )->type() == condition_type::dirichlet ) i += 1;
    if ( bc.getConditionFromSide( myeast )->type() == condition_type::dirichlet ) i += 2;

    switch ( i )
    {
        case 0: return;
        case 1: matrix = matrix.block( 1, 0, matrix.rows()-1, matrix.cols() ); return;
        case 2: matrix = matrix.block( 0, 0, matrix.rows()-1, matrix.cols() ); return;
        case 3: matrix = matrix.block( 1, 0, matrix.rows()-2, matrix.cols() ); return;
    }
}


template <index_t d, typename T>
void assembleParameterMomentsForTensorProduct(const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsFunction<T>& f, gsMatrix<T>& fh)
{
    const gsTensorBSplineBasis<d,T> * tb = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&basis);
    if( !tb )
        GISMO_ERROR ("This assembler only works for gsTensorBSplineBasis.");

    assembleParameterMoments(tb->component(0), f, fh, 0, d-1);
    handleDirichletConditionsForMoments(fh,bc,1+2*0,2+2*0);

    for ( index_t i=1; i<d; ++i )
    {
        gsMatrix<T> fhloc;
        assembleParameterMoments(tb->component(i), f, fhloc, i, d-1-i);
        handleDirichletConditionsForMoments(fhloc,bc,1+2*i,2+2*i);

        fh = fhloc.kron(fh);
    }

    if ( d>1 )
    {
        gsMatrix<T> pos, fvals;
        pos.resize(d,1);
        pos.fill(0.5); // TODO: is this appropriate?
        f.eval_into(pos, fvals);
        fh /= math::pow(fvals(0,0),d-1);
    }
}

template <typename T>
void assembleParameterMomentsForTensorProduct(const gsBasis<T>& basis, const gsBoundaryConditions<T>& bc, gsFunction<T>& f, gsMatrix<T>& fh)
{
    if (basis.dim() == 1)
        assembleParameterMomentsForTensorProduct<1,T>(basis, bc, f, fh);
    else if (basis.dim() == 2)
        assembleParameterMomentsForTensorProduct<2,T>(basis, bc, f, fh);
    else if (basis.dim() == 3)
        assembleParameterMomentsForTensorProduct<3,T>(basis, bc, f, fh);
    else if (basis.dim() == 4)
        assembleParameterMomentsForTensorProduct<4,T>(basis, bc, f, fh);
    else
    {
        GISMO_ERROR ("This assembler is instanciated only for up to 4 dimensions.");
    }
}


}

