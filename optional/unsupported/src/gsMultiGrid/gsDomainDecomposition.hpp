/** @file gsDomainDecomposition.hpp

    @brief Provides multi-patch preconditioners, particularly smoothers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsUtils/gsSortedVector.h>
//#include <gsCore/gsMemory.h>


namespace gismo
{


/// gsBasis

template<class T>
gsMatrix<index_t> gsBasisHelper<T>::boundary(const gsBasis<T>& basis, boxComponent2 b, bool strict)
{

    //gsInfo << "gsBasisHelper<T>::boundary(" << &basis << "," << b.m_index << "," << strict << ")\n" << std::flush;

    gsMatrix<index_t> boundary_indices;
    BasisUPtr result;
    index_t dim = basis.dim();

    if ( b.m_index == 0 )
    {
        const index_t sz = basis.size();
        boundary_indices.resize(sz,1);
        for (index_t i=0;i<sz;++i)
            boundary_indices(i,0) = i;
    }
    else
    {

        index_t d=0;
        while (b.m_index > 0)
        {
            if (b.m_index%3)
            {
                if (result)
                {
                    gsMatrix<index_t> tmp = result->boundary( boxSide( (b.m_index%3)+2*d ) );
                    for (index_t i=0; i<tmp.size(); ++i)
                        tmp(i,0) = boundary_indices(tmp(i,0),0);
                    tmp.swap(boundary_indices);
                    result = result->boundaryBasis( boxSide( (b.m_index%3)+2*d ) );
                }
                else
                {
                    boundary_indices = basis.boundary( boxSide( (b.m_index%3)+2*d ) );
                    result = basis.boundaryBasis( boxSide( (b.m_index%3)+2*d ) );
                }
                --dim;
            }
            else
            {
                ++d;
            }
            b.m_index /= 3;
        }
    }


    if (strict && dim > 0)
    {

        gsMatrix<index_t> fin_bdy;
        if (result)
            fin_bdy = result->allBoundary();
        else
            fin_bdy = basis.allBoundary();

        const index_t boundary_indices_sz = boundary_indices.rows();
        const index_t fin_bdy_sz = fin_bdy.rows();

        // Copy all entries from boundary_indices to boundary_indices_cleaned except
        // those with indices in fin_bdy

        gsMatrix<index_t> boundary_indices_cleaned(boundary_indices_sz - fin_bdy_sz,1);
        index_t j=0, t=0;
        for (index_t i=0; i<boundary_indices_sz; ++i)
        {
            if (i>fin_bdy(j,0) && j < fin_bdy_sz)
                ++j;
            if (i<fin_bdy(j,0) || j == fin_bdy_sz)
            {
                boundary_indices_cleaned(t,0) = boundary_indices(i,0);
                ++t;
            }
        }
        GISMO_ASSERT( t == boundary_indices_cleaned.rows(), "Internal error." );
        return boundary_indices_cleaned;
    }
    else
        return boundary_indices;
}

template<class T>
typename gsBasisHelper<T>::BasisUPtr gsBasisHelper<T>::boundaryBasis(const gsBasis<T>& basis, boxComponent2 b)
{
    if ( b.m_index == 0 )
        return basis.clone();

    BasisUPtr result;
    index_t d=1;
    while ( b.m_index > 0 )
    {
        if ( b.m_index%3 )
        {
            if (result)
                result = result->boundaryBasis( boxSide(b.m_index+2*d) );
            else
                result = basis.boundaryBasis( boxSide(b.m_index+2*d) );
        }
        b.m_index /= 3;
        ++d;
    }
    return result;
}

template<class T>
gsSparseMatrix<T,RowMajor> gsBasisHelper<T>::setupTransferMatrix( const gsMatrix<index_t>& coefs, index_t totalSize )
{
    GISMO_ASSERT(coefs.cols() == 1, "Assume a column vector." );
    index_t sz = coefs.rows();
    gsSparseEntries<T> se;
    se.reserve(sz);
    for (index_t i=0; i<sz; ++i)
    {
        GISMO_ASSERT( coefs.at(i) >= 0 && coefs.at(i) < static_cast<index_t>(totalSize), "Out of range." );
        se.add(coefs.at(i),i,T(1));
    }
    gsSparseMatrix<T,RowMajor> result(totalSize,sz);
    result.setFrom(se);
    return result;
}

template<class T>
gsSparseMatrix<T,RowMajor> gsBasisHelper<T>::setupTransferMatrix( const std::vector<index_t>& coefs, index_t totalSize )
{
    index_t sz = coefs.size();
    gsSparseEntries<T> se;
    se.reserve(sz);
    for (index_t i=0; i<sz; ++i)
    {
        GISMO_ASSERT( coefs[i] >= 0 && coefs[i] < static_cast<index_t>(totalSize), "Out of range." );
        se.add(coefs.at(i),i,T(1));
    }
    gsSparseMatrix<T,RowMajor> result(totalSize,sz);
    result.setFrom(se);
    return result;
}


struct patchCornerSort {

    bool operator() (const patchCorner& c1, const patchCorner& c2) const
    { return c1.patch<c2.patch || (c1.patch==c2.patch && c1.m_index<c2.m_index); }

};

/// gsMultiBasis
template<class T>
patchCorner getCanonicCorner( const patchCorner& c, const gsMultiBasis<T>& mb )
{
    std::vector< patchCorner > corners;
    mb.topology().getCornerList(c,corners);
    std::sort(corners.begin(), corners.end(), patchCornerSort());
    return corners[0];
}

template<class T>
std::vector<patchCorner> getCanonicCorners( const std::vector<patchCorner>& c, const gsMultiBasis<T>& mb )
{
    const index_t sz = c.size();
    std::vector< patchCorner > corners;
    corners.reserve(sz);
    for (index_t i=0; i<sz; ++i)
        corners.push_back( getCanonicCorner(c[i],mb) );
    std::sort(corners.begin(),corners.end());
    return corners;
}




std::vector<index_t> getCornerIndices( const std::vector<patchCorner>& corner, index_t dim )
{
    const index_t sz = corner.size();
    std::vector<index_t> result(sz);
    for (index_t i=0; i<sz; ++i)
        result[i] = corner[i].patch*(1u<<(dim)) + corner[i].m_index;
    return result;
}

template<class T>
void print(const typename gsMultiBasisHelper<T>::component& g)
{
    gsInfo << "[glob, contributes to patches ";
    for (size_t i=0; i<g.components.size(); ++i)
    {
        gsInfo << g.components[i].patch << ":";
        index_t factor = 1;
        for (index_t j=0; j<g.components[i].m_total_dim; ++j)
            factor *= 3;
        for (index_t j=0; j<g.components[i].m_total_dim; ++j)
        {
            gsInfo << (g.components[i].m_index % factor) / ( factor / 3 );
            factor /= 3;
        }
        gsInfo << " ";
    }
    gsInfo << "and consisting of corners ";
    for (size_t i=0; i<g.corners.size(); ++i)
        gsInfo << g.corners[i] << " ";
    gsInfo << "]\n";

}

template<typename T>
std::vector<typename gsMultiBasisHelper<T>::component> gsMultiBasisHelper<T>::getComponents(
    const gsMultiBasis<T>& mb,
    bool combineCorners
)
{
    const index_t nPatches = mb.nBases();
    const index_t dim = mb.dim();

    index_t cnr = 1;
    for (index_t i=0; i<dim; ++i) cnr *= 3;

    typedef std::map< std::vector<index_t>, component>   MapT;

    std::vector<MapT> globs(dim+1);

    for (index_t i = 0; i<nPatches; ++i)
    {
        for (index_t j = 0; j<cnr; ++j)
        {
            patchComponent2 pc(i, j, dim);
            const index_t d = pc.dim();
            std::vector< patchCorner > crns = getCanonicCorners(pc.containedCorners(),mb);
            component& g = globs[d][getCornerIndices(crns, mb.dim())];
            g.corners = give(crns);
            g.components.push_back(pc);
        }
    }
    index_t sz = 0;
    for (index_t i=0; i<dim+1; ++i)
        sz += globs[i].size();

    std::vector<component> result;
    result.reserve(sz);
    for (index_t i=dim; i>=0; --i)
    {
        if (!combineCorners || i>0)
            for( typename MapT::iterator it = globs[i].begin(); it != globs[i].end(); ++it )
                result.push_back( it->second );
        else
        {
            component last;
            for( typename MapT::iterator it = globs[i].begin(); it != globs[i].end(); ++it )
            {
                const index_t nrcr = it->second.corners.size();
                for (index_t j=0; j<nrcr; ++j)
                    last.corners.push_back(it->second.corners[i]);
                const index_t nrcp = it->second.components.size();
                for (index_t j=0; j<nrcp; ++j)
                    last.components.push_back(it->second.components[i]);
            }
        }
    }
/*
    {
        gsInfo << "DIMENSIONS: ";
        for (index_t i=dim; i>=0; --i)
            gsInfo << globs[i].size() << " // ";
        gsInfo << std::endl;
    }

    for (size_t i=0; i<result.size();++i)
        { gsInfo << i << ":"; print(result[i]); }
*/

    return result;
}

template<typename T>
gsSparseMatrix<T, RowMajor> gsMultiBasisHelper<T>::getTransfer(
    const std::vector<patchComponent2>& pc,
    const gsMultiBasis<T>& mb,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    bool strict
)
{
    gsDofMapper dm;
    mb.getMapper(
       (dirichlet::strategy)opt.askInt("DirichletStrategy",11),
       (iFace    ::strategy)opt.askInt("InterfaceStrategy", 1),
       bc,
       dm,
       0
    );
    const index_t nTotalDofs = dm.freeSize();

    std::vector<index_t> global_coefs;

    const index_t nrpc = pc.size();
    for (index_t n=0; n<nrpc; ++n)
    {
        gsMatrix<index_t> coefs = gsBasisHelper<T>().boundary( mb[pc[n].patch], pc[n], strict);

        const index_t sz = coefs.size();

        global_coefs.reserve(sz*pc.size()); //TODO

        for (index_t i=0; i<sz; ++i)
            if (dm.is_free(coefs(i,0), pc[n].patch))
            {
                const index_t new_idx = dm.index(coefs(i,0), pc[n].patch);
                if ( std::find(global_coefs.begin(), global_coefs.end(), new_idx) == global_coefs.end() )
                    global_coefs.push_back(new_idx);
            }
    }

    return gsBasisHelper<T>::setupTransferMatrix(global_coefs, nTotalDofs);
}

} // namespace gismo
