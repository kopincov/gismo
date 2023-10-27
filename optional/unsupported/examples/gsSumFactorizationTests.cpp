/** @file gsSumFactorizationTests.cpp

    @brief Provides test examples for sum factorization ideas

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, A. Bressan
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsSolver/gsKronecker.h>

using namespace gismo;

template <index_t dim, bool useSumFactorizationApproach, bool useSumFactorizationForGeo>
void assembleMass(const gsGeometry<real_t>& domain, const gsTensorBSplineBasis<dim,real_t>& tbasis, gsSparseMatrix<real_t>& M);
template<index_t dim, bool useSumFactorizationForGeo>
gsSparseMatrix<real_t> evaluateTensorGeometry( const gsGeometry<real_t>& domain, const std::vector< gsMatrix<real_t> >& quNodes );
// Evaluation functions
gsSparseMatrix<real_t> evaluateBasis( const gsBasis<real_t>& basis );
template<index_t dim>
std::vector< gsSparseMatrix<real_t> > evaluateTensorBasis( const gsTensorBSplineBasis<dim,real_t>& tbasis );
gsSparseMatrix<real_t> evaluateBasisDerivatives( const gsBasis<real_t>& basis );
template<index_t dim>
std::vector< gsSparseMatrix<real_t> > evaluateTensorBasisDerivatives( const gsTensorBSplineBasis<dim,real_t>& tbasis );
gsSparseMatrix<real_t> evaluateBasis( const gsBasis<real_t>& basis, const gsMatrix<real_t>& quNodes );
template<index_t dim>
std::vector< gsSparseMatrix<real_t> > evaluateTensorBasis( const gsTensorBSplineBasis<dim,real_t>& tbasis, const std::vector< gsMatrix<real_t> >& quNodes );
gsSparseMatrix<real_t> evaluateBasisDerivatives( const gsBasis<real_t>& basis, const gsMatrix<real_t>& quNodes );
template<index_t dim>
std::vector< gsSparseMatrix<real_t> > evaluateTensorBasisDerivatives( const gsTensorBSplineBasis<dim,real_t>& tbasis, const std::vector< gsMatrix<real_t> >& quNodes );
template <index_t dim>
gsSparseMatrix<real_t> getQuWeights( const gsTensorBSplineBasis<dim,real_t>&  basis );
template<index_t dim>
std::vector< gsSparseMatrix<real_t> > getTensorQuWeights( const gsTensorBSplineBasis<dim,real_t>& tbasis );
template <index_t dim>
gsMatrix<real_t> getQuNodes( const gsTensorBSplineBasis<dim,real_t>&  basis );
template<index_t dim>
std::vector< gsMatrix<real_t> > getTensorQuNodes( const gsTensorBSplineBasis<dim,real_t>& tbasis );
// Linear algebra functions
gsSparseMatrix<real_t> getSelectiveKroneckerProduct(index_t idx, const std::vector< gsSparseMatrix<real_t> >& matricesA, const std::vector< gsSparseMatrix<real_t> >& matricesB);
gsSparseMatrix<real_t> getKroneckerProductWithIdentities(index_t pre, const gsSparseMatrix<real_t>& A, index_t post);
gsSparseMatrix<real_t> smartTensorProduct( const std::vector<gsSparseMatrix<real_t> > & a, const gsSparseMatrix<real_t> & b, const std::vector<gsSparseMatrix<real_t> > & c );
gsSparseMatrix<real_t> getDiagonalMatrix( const gsVector<real_t>& data );
gsMatrix<real_t> getComponentwiseDeterminants( const std::vector< std::vector< gsMatrix<real_t> > >& data );
std::vector< gsLinearOperator<real_t>::Ptr > getGradients( std::vector< gsSparseMatrix<real_t> > basisValues, std::vector< gsSparseMatrix<real_t> > basisDerivs );

    void assemble1 (const gsSparseMatrix<> &M_reference, gsGeometry<>::uPtr &geo, gsBasis<real_t>::uPtr &tbasis);
    void assemble2 (const gsSparseMatrix<> &M_reference, gsGeometry<>::uPtr &geo, gsBasis<real_t>::uPtr &tbasis);
    void assemble3 (const gsSparseMatrix<> &M_reference, gsGeometry<>::uPtr &geo, gsBasis<real_t>::uPtr &tbasis);

int main(int argc, char *argv[])
{
    std::string geometry("4");
    index_t numIntervals    = 50 ;
    index_t degree          = 5  ;

    gsCmdLine cmd("Assembles the mass matrix with a few approaches and compares the result.");
    cmd.addString("g", "geometry",  "Geometry index",               geometry                );
    cmd.addInt   ("n", "intervals", "Number of intervals",          numIntervals            );
    cmd.addInt   ("p", "degree",    "Spline degree to be used",     degree                  );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (numIntervals < 1) { gsInfo << "Error: Number of intervals must be positive.\n"; return 1; }
    if (degree < 1)       { gsInfo << "Error: Degree must be positive.\n";              return 1; }

    gsGeometry<>::uPtr geo;

    if (geometry=="1")
    {
        gsInfo << "Geometry:              BSplineUnitInterval\n";
        geo = gsNurbsCreator<>::BSplineUnitInterval(static_cast<short_t>(degree));
    }
    else if (geometry=="2")
    {
        gsInfo << "Geometry:              BSplineSquare\n";
        geo = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(degree));
    }
    else if (geometry=="3")
    {
        gsInfo << "Geometry:              BSplineCube\n";
        geo = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(degree));
    }
    else if (geometry=="4")
    {
        gsInfo << "Geometry:              BSplineQuarterAnnulus\n";
        geo = gsNurbsCreator<>::BSplineQuarterAnnulus(static_cast<short_t>(2));
    }
    else if ( gsFileManager::fileExists(geometry) )
    {
        geometry = gsFileManager::fileExists(geometry);
        gsInfo << "Geometry:              " << geometry << "\n";
        gsMultiPatch<> mp;
        gsFileData<> fileData(geometry);
        if (!fileData.has< gsMultiPatch<> >())
        { gsInfo << "Error: No multipatch object found.\n"; return 1; }
        fileData.getFirst< gsMultiPatch<> >(mp);
        if ( mp.nPatches() != 1 )
        { gsInfo << "Error: There is not exactly one patch.\n"; return 1; }
        geo = mp[0].clone();
    }
    else
    { gsInfo << "Error: Invalid geometry.\n"; return 1; }


    gsKnotVector<real_t> KV(0.0, 1.0, numIntervals-1, static_cast<short_t>(degree+1));
    gsInfo << "Knots:                 " << KV << std::endl;
    gsBasis<real_t>::uPtr tbasis;
    switch (geo->geoDim())
    {
    case 1: tbasis = memory::make_unique<gsBasis<real_t> >(new gsBSplineBasis<real_t>( KV )); break;
    case 2: tbasis = memory::make_unique<gsBasis<real_t> >(new gsTensorBSplineBasis<2,real_t>( KV, KV )); break;
    case 3: tbasis = memory::make_unique<gsBasis<real_t> >(new gsTensorBSplineBasis<3,real_t>( KV, KV, KV )); break;
    default: return 1;
    }

    gsInfo << "Discretization space:  "
           << "dim=" << tbasis->dim()
           << ", deg=" << tbasis->degree(0)
           << ", numElements=" << tbasis->numElements()
           << ", dofs=" << tbasis->size() << "\n\n";

    gsSparseMatrix<real_t> M_reference;
    {
        gsInfo << "Assembling mass matrix with gsGenericAssembler... " << std::endl;
        gsStopwatch time;
        gsBoundaryConditions<real_t> bc;
        gsGenericAssembler<real_t> genassm(gsMultiPatch<>(*geo), gsMultiBasis<>(*tbasis), gsGenericAssembler<real_t>::defaultOptions(), &bc);
        genassm.assembleMass();
        M_reference = genassm.fullMatrix();
        double duration = time.stop();
        gsInfo << "Done in " << duration << " secs.\n" << std::flush;
        gsInfo << "Obtain a " << M_reference.rows() << " x " << M_reference.cols() << " matrix.\n\n";
    }


    assemble1 (M_reference, geo, tbasis);
    assemble2 (M_reference, geo, tbasis);
    assemble3 (M_reference, geo, tbasis);

    return 0;
}


void assemble1 (const gsSparseMatrix<> &M_reference, gsGeometry<>::uPtr &geo, gsBasis<real_t>::uPtr &tbasis)
{
    gsInfo << "Assembling mass matrix with Approach 1...         " << std::endl;
    gsStopwatch time;
    gsSparseMatrix<real_t> M;
    switch (geo->geoDim())
    {
    case 1: assembleMass<1,false,false>(*geo, dynamic_cast<gsTensorBSplineBasis<1,real_t>&>(*tbasis), M); break;
    case 2: assembleMass<2,false,false>(*geo, dynamic_cast<gsTensorBSplineBasis<2,real_t>&>(*tbasis), M); break;
    case 3: assembleMass<3,false,false>(*geo, dynamic_cast<gsTensorBSplineBasis<3,real_t>&>(*tbasis), M); break;
    }
    double duration = time.stop();
    gsInfo << "Done in " << duration << " secs.\n" << std::flush;
    gsInfo << "Obtain a " << M.rows() << " x " << M.cols() << " matrix.\n";
    gsInfo << "Relative error: " << (M-M_reference).norm() / M_reference.norm() << "\n\n";
}

void assemble2 (const gsSparseMatrix<> &M_reference, gsGeometry<>::uPtr &geo, gsBasis<real_t>::uPtr &tbasis)
{
    gsInfo << "Assembling mass matrix with Approach 2...         " << std::endl;
    gsStopwatch time;
    gsSparseMatrix<real_t> M;
    switch (geo->geoDim())
    {
    case 1: assembleMass<1,true,false>(*geo, dynamic_cast<gsTensorBSplineBasis<1,real_t>&>(*tbasis), M); break;
    case 2: assembleMass<2,true,false>(*geo, dynamic_cast<gsTensorBSplineBasis<2,real_t>&>(*tbasis), M); break;
    case 3: assembleMass<3,true,false>(*geo, dynamic_cast<gsTensorBSplineBasis<3,real_t>&>(*tbasis), M); break;
    }
    double duration = time.stop();
    gsInfo << "Done in " << duration << " secs.\n" << std::flush;
    gsInfo << "Obtain a " << M.rows() << " x " << M.cols() << " matrix.\n";
    gsInfo << "Relative error: " << (M-M_reference).norm() / M_reference.norm() << "\n\n";
}

void assemble3 (const gsSparseMatrix<> &M_reference, gsGeometry<>::uPtr &geo, gsBasis<real_t>::uPtr &tbasis)
{
    gsInfo << "Assembling mass matrix with Approach 3...         " << std::endl;
    gsStopwatch time;
    gsSparseMatrix<real_t> M;
    switch (geo->geoDim())
    {
    case 1: assembleMass<1,true,true>(*geo, dynamic_cast<gsTensorBSplineBasis<1,real_t>&>(*tbasis), M); break;
    case 2: assembleMass<2,true,true>(*geo, dynamic_cast<gsTensorBSplineBasis<2,real_t>&>(*tbasis), M); break;
    case 3: assembleMass<3,true,true>(*geo, dynamic_cast<gsTensorBSplineBasis<3,real_t>&>(*tbasis), M); break;
    }
    double duration = time.stop();
    gsInfo << "Done in " << duration << " secs.\n" << std::flush;
    gsInfo << "Obtain a " << M.rows() << " x " << M.cols() << " matrix.\n";
    gsInfo << "Relative error: " << (M-M_reference).norm() / M_reference.norm() << "\n\n";
}




template <index_t dim, bool useSumFactorizationApproach, bool useSumFactorizationForGeo>
void assembleMass(const gsGeometry<real_t>& domain, const gsTensorBSplineBasis<dim,real_t>& tbasis, gsSparseMatrix<real_t>& M)
{

    gsStopwatch time1;
    std::vector< gsSparseMatrix<real_t> > quWeights = getTensorQuWeights<dim>( tbasis );
    gsInfo << "  Computed weights in " << time1.stop() << " secs.\n" << std::flush;

    gsStopwatch time2;
    std::vector< gsSparseMatrix<real_t> > basisValues = evaluateTensorBasis<dim>( tbasis );
    gsInfo << "  Computed basis in " << time2.stop() << " secs.\n" << std::flush;

    gsStopwatch time3;
    std::vector< gsMatrix<real_t> > quNodes = getTensorQuNodes<dim>(tbasis);
    gsInfo << "  Computed nodes in " << time3.stop() << " secs.\n" << std::flush;

    gsStopwatch time4;
    gsInfo << "  Compute geometry...\n";
    gsSparseMatrix<real_t> geo = evaluateTensorGeometry<dim,useSumFactorizationForGeo>(domain, quNodes);
    gsInfo << "  Computed geometry in " << time4.stop() << " secs.\n" << std::flush;

    gsStopwatch time5;
    if (useSumFactorizationApproach)
        M = smartTensorProduct( basisValues, geo * kroneckerProduct( quWeights ), basisValues );
    else
    {
        gsSparseMatrix<real_t> tmp = kroneckerProduct( basisValues );
        gsSparseMatrix<real_t> tmp2 = geo * kroneckerProduct( quWeights );
        M = tmp.transpose() * tmp2 * tmp;
    }
    gsInfo << "  Computed products in " << time5.stop() << " secs.\n" << std::flush;

}


template<index_t dim, bool useSumFactorizationForGeo>
gsSparseMatrix<real_t> evaluateTensorGeometry( const gsGeometry<real_t>& domain, const std::vector< gsMatrix<real_t> >& quNodes )
{
    gsMatrix<real_t> coefs = domain.coefs();
    typedef gsTensorBSplineBasis<dim, real_t> TensorBasisType;
    const TensorBasisType* tbasis = dynamic_cast<const TensorBasisType*>(& domain.basis());
    GISMO_ENSURE( tbasis, "The underlying basis is not appropriate." );

    gsStopwatch time1;
    std::vector< gsSparseMatrix<real_t> > basisValues = evaluateTensorBasis<dim>( *tbasis, quNodes );
    gsInfo << "    Computed basis in " << time1.stop() << " secs.\n" << std::flush;

    gsStopwatch time2;
    std::vector< gsSparseMatrix<real_t> > basisDerivs = evaluateTensorBasisDerivatives<dim>( *tbasis, quNodes );
    gsInfo << "    Computed basis derivatives in " << time2.stop() << " secs.\n" << std::flush;

    gsStopwatch time3;
    std::vector< std::vector< gsMatrix<real_t> > > jacobi( dim );
    if (!useSumFactorizationForGeo)
    {
        for (index_t i=0; i<dim; ++i)
        {
            jacobi[i].resize(dim);
            for (index_t j=0; j<dim; ++j)
                jacobi[i][j] = getSelectiveKroneckerProduct( i, basisValues, basisDerivs ) * coefs.col(dim-1-j);
        }
    }
    else
    {
        std::vector< gsLinearOperator<real_t>::Ptr > grads = getGradients( basisValues, basisDerivs );
        for (index_t i=0; i<dim; ++i)
        {
            jacobi[i].resize(dim);
            for (index_t j=0; j<dim; ++j)
                grads[i]->apply( coefs.col(dim-1-j), jacobi[i][j] );
        }
    }
    gsInfo << "    Computed Jacobi matrix in " << time3.stop() << " secs.\n" << std::flush;

    gsStopwatch time4;
    gsSparseMatrix<real_t> result = getDiagonalMatrix( getComponentwiseDeterminants( jacobi ) );
    gsInfo << "    Computed Jacobi determinant in " << time4.stop() << " secs.\n" << std::flush;

    return result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Evaluation functions

gsSparseMatrix<real_t> evaluateBasis( const gsBasis<real_t>& basis )
{
    GISMO_ASSERT( basis.dim() == 1, "This is only implemented for the univariate case." );
    gsBasis<real_t>::domainIter domIt = basis.makeDomainIterator();

    const index_t nrQuNodes   = basis.degree(0) + 1;
    const index_t nrActiveFct = basis.degree(0) + 1;
    const index_t nrElements  = domIt->numElements();
    const index_t nrBasisFct  = basis.size();

    gsSparseEntries<real_t> se;
    se.reserve( nrQuNodes * nrElements * nrActiveFct );

    gsVector<index_t> numQuadNodes(1);
    numQuadNodes[0] = nrQuNodes;

    gsGaussRule<> quRule(numQuadNodes);
    gsMatrix<> quNodes;
    gsVector<> quWeights;
    gsFuncData<> bdata(NEED_VALUE|NEED_ACTIVE|SAME_ELEMENT);

    index_t n = 0;
    for (; domIt->good(); domIt->next())
    {
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
        basis.compute(quNodes, bdata);

        for (index_t i=0; i<nrQuNodes; ++i)
            for (index_t j=0; j<nrActiveFct; ++j)
                se.add( n*nrQuNodes + i, n+j, bdata.values[0](j,i) );

        n++;
    }
    gsSparseMatrix<real_t> result( nrQuNodes * nrElements, nrBasisFct );
    result.setFromTriplets( se.begin(), se.end() );
    result.makeCompressed();
    return result;
}

template<index_t dim>
std::vector< gsSparseMatrix<real_t> > evaluateTensorBasis( const gsTensorBSplineBasis<dim,real_t>& tbasis )
{
    std::vector< gsSparseMatrix<real_t> > result(dim);
    for (index_t d=0; d<dim; ++d)
        result[dim-1-d] = evaluateBasis(tbasis.component(d));
    return result;
}

gsSparseMatrix<real_t> evaluateBasisDerivatives( const gsBasis<real_t>& basis )
{
    GISMO_ASSERT( basis.dim() == 1, "This is only implemented for the univariate case." );
    gsBasis<real_t>::domainIter domIt = basis.makeDomainIterator();

    const index_t nrQuNodes   = basis.degree(0) + 1;
    const index_t nrActiveFct = basis.degree(0) + 1;
    const index_t nrElements  = domIt->numElements();
    const index_t nrBasisFct  = basis.size();

    gsSparseEntries<real_t> se;
    se.reserve( nrQuNodes * nrElements * nrActiveFct );

    gsVector<index_t> numQuadNodes(1);
    numQuadNodes[0] = nrQuNodes;

    gsGaussRule<> quRule(numQuadNodes);
    gsMatrix<> quNodes;
    gsVector<> quWeights;
    gsFuncData<> bdata(NEED_DERIV|NEED_ACTIVE|SAME_ELEMENT);

    index_t n = 0;
    for (; domIt->good(); domIt->next())
    {
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
        basis.compute(quNodes, bdata);

        for (index_t i=0; i<nrQuNodes; ++i)
            for (index_t j=0; j<nrActiveFct; ++j)
                se.add( n*nrQuNodes + i, n+j, bdata.values[1](j,i) );

        n++;
    }
    gsSparseMatrix<real_t> result( nrQuNodes * nrElements, nrBasisFct );
    result.setFromTriplets( se.begin(), se.end() );
    result.makeCompressed();
    return result;
}

template<index_t dim>
std::vector< gsSparseMatrix<real_t> > evaluateTensorBasisDerivatives( const gsTensorBSplineBasis<dim,real_t>& tbasis )
{
    std::vector< gsSparseMatrix<real_t> > result(dim);
    for (index_t d=0; d<dim; ++d)
        result[dim-1-d] = evaluateBasisDerivatives(tbasis.component(d));
    return result;
}

gsSparseMatrix<real_t> evaluateBasis( const gsBasis<real_t>& basis, const gsMatrix<real_t>& quNodes )
{
    GISMO_ASSERT( basis.dim() == 1, "This is only implemented for the univariate case." );
    gsBasis<real_t>::domainIter domIt = basis.makeDomainIterator();

    const index_t totNrQuNodes = quNodes.rows();
    const index_t nrActiveFct  = basis.degree(0) + 1;
    const index_t nrBasisFct   = basis.size();

    gsSparseEntries<real_t> se;
    se.reserve( totNrQuNodes * nrActiveFct );

    gsMatrix<real_t> vals;
    gsMatrix<index_t> actives;
    for (index_t n=0; n<totNrQuNodes; ++n)
    {
        basis.eval_into( quNodes.row(n), vals );
        basis.active_into( quNodes.row(n), actives );
        GISMO_ASSERT( vals.cols() == actives.cols() && vals.rows() == actives.rows() && vals.cols() == 1, "Inconsistent basis response." );
        for (index_t j=0; j<vals.rows(); ++j)
            se.add( n, actives(j,0), vals(j,0) );
    }

    gsSparseMatrix<real_t> result( totNrQuNodes, nrBasisFct );
    result.setFromTriplets( se.begin(), se.end() );
    result.makeCompressed();
    return result;
}

template<index_t dim>
std::vector< gsSparseMatrix<real_t> > evaluateTensorBasis( const gsTensorBSplineBasis<dim,real_t>& tbasis, const std::vector< gsMatrix<real_t> >& quNodes )
{
    GISMO_ASSERT( quNodes.size() == dim, "Sizes do not agree." );
    std::vector< gsSparseMatrix<real_t> > result(dim);
    for (index_t d=0; d<dim; ++d)
        result[dim-1-d] = evaluateBasis(tbasis.component(d),quNodes[d]);
    return result;
}

gsSparseMatrix<real_t> evaluateBasisDerivatives( const gsBasis<real_t>& basis, const gsMatrix<real_t>& quNodes )
{
    GISMO_ASSERT( basis.dim() == 1, "This is only implemented for the univariate case." );
    gsBasis<real_t>::domainIter domIt = basis.makeDomainIterator();

    const index_t totNrQuNodes = quNodes.rows();
    const index_t nrActiveFct  = basis.degree(0) + 1;
    const index_t nrBasisFct   = basis.size();

    gsSparseEntries<real_t> se;
    se.reserve( totNrQuNodes * nrActiveFct );

    gsMatrix<real_t> vals;
    gsMatrix<index_t> actives;
    for (index_t n=0; n<totNrQuNodes; ++n)
    {
        basis.deriv_into( quNodes.row(n), vals );
        basis.active_into( quNodes.row(n), actives );
        GISMO_ASSERT( vals.cols() == actives.cols() && vals.rows() == actives.rows() && vals.cols() == 1, "Inconsistent basis response." );
        for (index_t j=0; j<vals.rows(); ++j)
            se.add( n, actives(j,0), vals(j,0) );
    }

    gsSparseMatrix<real_t> result( totNrQuNodes, nrBasisFct );
    result.setFromTriplets( se.begin(), se.end() );
    result.makeCompressed();
    return result;
}

template<index_t dim>
std::vector< gsSparseMatrix<real_t> > evaluateTensorBasisDerivatives( const gsTensorBSplineBasis<dim,real_t>& tbasis, const std::vector< gsMatrix<real_t> >& quNodes )
{
    GISMO_ASSERT( quNodes.size() == dim, "Sizes do not agree." );
    std::vector< gsSparseMatrix<real_t> > result(dim);
    for (index_t d=0; d<dim; ++d)
        result[dim-1-d] = evaluateBasisDerivatives(tbasis.component(d),quNodes[d]);
    return result;
}

template <index_t dim>
gsSparseMatrix<real_t> getQuWeights( const gsTensorBSplineBasis<dim,real_t>&  basis )
{
    gsBasis<real_t>::domainIter domIt = basis.makeDomainIterator();
    index_t sz = domIt->numElements();

    gsVector<index_t> numQuadNodes( basis.dim() );
    for (int i = 0; i < basis.dim(); ++i)
    {
        numQuadNodes[i] = basis.degree(i) + 1;
        sz *= numQuadNodes[i];  // TODO: probably not correct.
    }

    gsSparseEntries<real_t> se;
    se.reserve( sz );

    gsGaussRule<> quRule(numQuadNodes);
    gsMatrix<> quNodes;
    gsVector<> quWeights;

    index_t idx = 0;
    for (; domIt->good(); domIt->next())
    {
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

        for (index_t i=0; i < quWeights.size(); ++i)
        {
            se.add( idx, idx, quWeights[i] );
            idx++;
        }
    }
    gsDebugIf(idx != sz, sz);// << "se.reserve(sz) needs to be optimized, only " << idx << " is needed.";

    gsSparseMatrix<real_t> result( sz, sz );
    result.setFromTriplets( se.begin(), se.end() );
    result.makeCompressed();
    return result;
}

template<index_t dim>
std::vector< gsSparseMatrix<real_t> > getTensorQuWeights( const gsTensorBSplineBasis<dim,real_t>& tbasis )
{
    std::vector< gsSparseMatrix<real_t> > result(dim);
    for (index_t d=0; d<dim; ++d)
        result[dim-1-d] = getQuWeights<1>(tbasis.component(d));
    return result;
}

template <index_t dim>
gsMatrix<real_t> getQuNodes( const gsTensorBSplineBasis<dim,real_t>&  basis )
{
    gsBasis<real_t>::domainIter domIt = basis.makeDomainIterator();
    index_t sz = domIt->numElements();

    gsVector<index_t> numQuadNodes( basis.dim() );
    for (int i = 0; i < basis.dim(); ++i)
    {
        numQuadNodes[i] = basis.degree(i) + 1;
        sz *= numQuadNodes[i];
    }

    gsMatrix<real_t> result( sz, dim );

    //domIt->computeQuadratureRule( numQuadNodes );
    gsGaussRule<> quRule(numQuadNodes);
    gsMatrix<> quNodes;
    gsVector<> quWeights;

    index_t idx = 0;
    for (; domIt->good(); domIt->next())
    {
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

        for (index_t i=0; i < quNodes.cols(); ++i)
        {
            result.row(idx) = quNodes.col(i).transpose();
            idx++;
        }
    }

    return result;
}

template<index_t dim>
std::vector< gsMatrix<real_t> > getTensorQuNodes( const gsTensorBSplineBasis<dim,real_t>& tbasis )
{
    std::vector< gsMatrix<real_t> > result(dim);
    for (index_t d=0; d<dim; ++d)
        result[dim-1-d] = getQuNodes<1>(tbasis.component(d));
    return result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Linear algebra functions


gsSparseMatrix<real_t> getSelectiveKroneckerProduct(index_t idx, const std::vector< gsSparseMatrix<real_t> >& matricesA, const std::vector< gsSparseMatrix<real_t> >& matricesB)
{
    GISMO_ASSERT( matricesA.size() == matricesB.size(), "Sizes do not agree" );

    if ( matricesA.size() == 0 )
        return gsSparseMatrix<real_t>(0,0);

    gsSparseMatrix<real_t> result, temp;

    if ( idx == 0 )
        result = matricesB[0];
    else
        result = matricesA[0];

    index_t sz = matricesA.size();

    for ( index_t i=1; i<sz; ++i )
    {
        result.swap(temp);
        if ( idx == i )
            result=temp.kron(matricesB[i]);
        else
            result=temp.kron(matricesA[i]);

    }
    return result;

}

gsSparseMatrix<real_t> getKroneckerProductWithIdentities(index_t pre, const gsSparseMatrix<real_t>& A, index_t post)
{
    const index_t Ar = A.rows(), Ac = A.cols();

    gsSparseMatrix<real_t> result(pre*Ar*post, pre*Ac*post);
    result.resizeNonZeros(0);

    if( pre*Ar*post == 0 || pre*Ac*post == 0 )
        return result;

    //typedef gsSparseMatrix<real_t> Dest;
    typedef typename gsSparseMatrix<real_t>::InnerIterator InnerIterator;


    //TODO faster without gsSparseEntries; we know the number of entries per row/col?
    gsSparseEntries<real_t> se;
    se.reserve(A.nonZeros()*pre*post);

    for (index_t kA=0; kA < A.outerSize(); ++kA)
    {
        for (InnerIterator itA(A,kA); itA; ++itA)
        {
            const index_t itAr = itA.row(),
                    itAc = itA.col();
            const real_t  itAv = itA.value();

            for (index_t a=0; a<pre; ++a)
            {
                for (index_t b=0; b<post; ++b)
                {
                    const index_t i = (a*Ar + itAr) * post + b,
                            j = (a*Ac + itAc) * post + b;
                    se.add(i,j,itAv);
                }
            }
        }
    }
    result.setFromTriplets( se.begin(), se.end() );
    result.makeCompressed();
    return result;
}

gsSparseMatrix<real_t> smartTensorProduct( const std::vector<gsSparseMatrix<real_t> > & a, const gsSparseMatrix<real_t> & b, const std::vector<gsSparseMatrix<real_t> > & c )
{
    GISMO_ASSERT( a.size() == c.size(), "Assume that the sizes match." );
    const index_t dim = a.size();

    gsSparseMatrix<real_t> result = b;

    index_t corr1 = 1, corr2 = 1, corr3 = 1, corr4 = 1;
    for (index_t d = 0; d<dim; ++d) { corr2 *= a[d].rows(); corr4 *= c[d].rows(); }

    for (index_t d = 0; d<dim; ++d)
    {
        corr2 /= a[d].rows();
        corr4 /= c[d].rows();
        if ( false )
        {
            gsSparseMatrix<real_t> id1( corr1, corr1 ); id1.setIdentity();
            gsSparseMatrix<real_t> id2( corr2, corr2 ); id2.setIdentity();
            gsSparseMatrix<real_t> id3( corr3, corr3 ); id3.setIdentity();
            gsSparseMatrix<real_t> id4( corr4, corr4 ); id4.setIdentity();
            result = id2.kron(a[d]).kron(id1).transpose() * result * id4.kron(c[d]).kron(id3);
        }
        else
            result = getKroneckerProductWithIdentities( corr2, a[d], corr1 ).transpose() * result * getKroneckerProductWithIdentities( corr4, c[d], corr3 );
        corr1 *= a[d].cols();
        corr3 *= c[d].cols();
    }

    return result;
}

gsSparseMatrix<real_t> getDiagonalMatrix( const gsVector<real_t>& data )
{
    const index_t sz = data.size();
    gsSparseMatrix<real_t> result( sz, sz );
    result.reserve(1);
    for (index_t i=0; i<sz; ++i)
        result(i,i) = data[i];
    result.makeCompressed();
    return result;
}

gsMatrix<real_t> _getComponentwiseDeterminants( const std::vector< std::vector< const gsMatrix<real_t>* > >& data )
{
    GISMO_ASSERT( data.size() > 0, "Zero blocks are not allowed." );

    if (data.size() == 1)
    {
        GISMO_ASSERT( data[0].size() == 1, "Sizes inconsistent." );
        return *(data[0][0]);
    }
    else if (data.size() == 2)
    {
        GISMO_ASSERT( data[0].size() == 2 && data[1].size() == 2, "Sizes inconsistent." );
        return (data[0][0]->array() * data[1][1]->array() - data[0][1]->array() * data[1][0]->array()).matrix();
    }
    else
    {
        GISMO_ASSERT( data[0].size() == data.size(), "Sizes inconsistent." );
        const index_t sz = data.size();

        gsMatrix<real_t> result( data[0][0]->rows(), data[0][0]->cols() );
        result.setZero();

        for (index_t i=0; i<sz; ++i)
        {
            std::vector< std::vector< const gsMatrix<real_t>* > > submat( sz-1 );
            for (index_t a=0; a<sz-1; ++a)
            {
                submat[a].resize(sz-1);
                for (index_t b=0; b<sz-1; ++b)
                    submat[a][b] = data[a+1][ (b<i) ? b : b+1 ];
            }
            if (i%2)
                result -= (data[0][i]->array() * _getComponentwiseDeterminants(submat).array()).matrix();
            else
                result += (data[0][i]->array() * _getComponentwiseDeterminants(submat).array()).matrix();
        }

        return result;
    }
}

gsMatrix<real_t> getComponentwiseDeterminants( const std::vector< std::vector< gsMatrix<real_t> > >& data )
{
    const index_t sz = data.size();
    std::vector< std::vector< const gsMatrix<real_t>* > > data_ptrs(sz);
    for (index_t i=0; i<sz; ++i)
    {
        data_ptrs[i].resize(data[i].size());
        for (index_t j=0; j<(index_t)data[i].size(); ++j)
            data_ptrs[i][j] = &(data[i][j]);
    }

    return _getComponentwiseDeterminants( data_ptrs );
}


std::vector< gsLinearOperator<real_t>::Ptr > getGradients( std::vector< gsSparseMatrix<real_t> > basisValues, std::vector< gsSparseMatrix<real_t> > basisDerivs )
{
    GISMO_ASSERT( basisValues.size() == basisDerivs.size(), "Sizes inconsistent." );

    const index_t dim = basisValues.size();
    std::vector< gsLinearOperator<real_t>::Ptr > result( dim ), basisValuesPtr( dim ), basisDerivsPtr( dim );

    for (index_t i=0; i<dim; ++i)
        basisValuesPtr[i] = makeMatrixOp( basisValues[i].moveToPtr() );

    for (index_t i=0; i<dim; ++i)
        basisDerivsPtr[i] = makeMatrixOp( basisDerivs[i].moveToPtr() );

    for (index_t i=0; i<dim; ++i)
    {
        std::vector< gsLinearOperator<real_t>::Ptr > tmp;
        for (index_t j=0; j<dim; ++j)
            if (i==j)
                tmp.push_back( basisDerivsPtr[j] );
            else
                tmp.push_back( basisValuesPtr[j] );
        result[i] = gsKroneckerOp<real_t>::make(tmp);
    }
    return result;
}
