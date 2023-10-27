/** @file gsMultiPatchPreconditioners.hpp

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

template<typename T>
void gsAdditivePreconditionerOp<T>::step(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_A->rows() == x.rows() && x.rows() == f.rows() && m_A->cols() == m_A->rows() && x.cols() == f.cols(),
        "The dimensions do not fit." );

    const index_t n = m_ops.size();
    m_A->apply( x, m_res );
    m_res -= f;

    for (index_t i=0; i<n; ++i)
    {
        m_res_local.noalias() = m_transfers[i]*m_res;
        m_ops[i]->apply(m_res_local, m_corr_local);
        m_corr_local *= m_damping;
        x.noalias() -= m_transfers[i].transpose()*m_corr_local;
    }
}

template<typename T>
void gsMultiplicativePreconditionerOp<T>::step(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_A->rows() == x.rows() && x.rows() == f.rows() && m_A->cols() == m_A->rows() && x.cols() == f.cols(),
        "The dimensions do not fit." );

    const index_t n = m_ops.size();
    for (index_t i=0; i<n; ++i)
    {
        m_A->apply( x, m_res ); //TODO: recomputung the residual here is inefficient; we could do better for matrices
        m_res -= f;
        m_ops[i]->apply(m_transfers[i]*m_res, m_corr_local);
        m_corr_local *= m_damping;
        x.noalias() -= m_transfers[i].transpose()*m_corr_local;
    }
}

template<typename T>
std::vector< gsSparseMatrix<T> > getPatchwiseTransfers(
                                    const gsMultiBasis<T>& mb,
                                    const gsBoundaryConditions<T>& bc,
                                    const gsOptionList& opt
                                )
{
    std::vector< gsSparseMatrix<T> > transfers;

    const index_t n = mb.nBases();
    gsDofMapper dm;
    mb.getMapper(
       (dirichlet::strategy)opt.askInt("DirichletStrategy",11),
       (iFace    ::strategy)opt.askInt("InterfaceStrategy", 1),
       bc,
       dm,
       0
    );
    const index_t nTotalDofs = dm.freeSize();
    for (index_t i=0; i<n; ++i)
    {
        const index_t nDofs = mb[i].size();
        gsSparseEntries<T> transfer_se;
        transfer_se.reserve(nDofs);
        for (index_t j=0; j<nDofs; ++j)
        {
            const index_t dof_idx = dm.index(j,i);
            if (dm.is_free_index(dof_idx)) //TODO: check this!
                transfer_se.add(j,dof_idx,1);
        }
        gsSparseMatrix<T> transfer(nDofs, nTotalDofs);
        transfer.setFrom(transfer_se);
        transfer.makeCompressed();
        transfers.push_back(give(transfer));
    }
    return transfers;
}


namespace {
template<typename T>
void _setupDofMapHelper(index_t i, index_t index_for_corners, patchSide patch_side, const gsMultiBasis<T>& mb, const gsDofMapper& dm, std::vector<index_t>& dof_map)
{
    index_t patch_index = patch_side.patch;

    gsMatrix<index_t> bdy = mb.basis(patch_index).boundary(patch_side);
    const index_t nDof = bdy.rows();

    for (index_t j=0; j<nDof; ++j)
    {
        const index_t local_idx = bdy(j,0);
        const index_t idx = dm.index(local_idx,patch_index);

        if ( dm.is_free_index( idx ) )
        {
            if ( dof_map[idx] == -2 )
                dof_map[idx] = i;
            else if ( dof_map[idx] != i )
                dof_map[idx] = index_for_corners;
        }
    }

}
} // anonymous namespace

template<class T>
std::vector< std::vector< gsSparseMatrix<T> > > getPiecewiseTransfers(
                                    const gsMultiBasis<T>& mb,
                                    const gsBoundaryConditions<T>& bc,
                                    const gsOptionList& opt
                                )
{
    std::vector< std::vector< gsSparseMatrix<T> > > result( mb.dim()+1 );
    std::vector<typename gsMultiBasisHelper<T>::component> data = gsMultiBasisHelper<T>::getComponents(mb,true);
    const index_t sz = data.size();
    for (index_t i = 0; i<sz; ++i)
    {
        gsSparseMatrix<T> sm = gsMultiBasisHelper<T>::getTransfer(data[i].components,mb,bc,opt,true).transpose();
        // TODO: reserve for result[...]
        result[ mb.dim() - data[i].components[0].dim() ].push_back(give(sm));
    }
    return result;
}


/*template<typename T>
std::vector< std::vector< gsSparseMatrix<T> > > getPiecewiseTransfers(
                                    const gsMultiBasis<T>& mb,
                                    const gsBoundaryConditions<T>& bc,
                                    const gsOptionList& opt
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

    index_t nPatches = mb.nBases();

    std::vector<index_t> dof_map(nTotalDofs,-2);

    // interfaces
    std::vector<boundaryInterface> intf = mb.topology().interfaces();
    std::vector<patchSide> bdy = mb.topology().boundaries();
    const index_t nIntf = intf.size();
    const index_t nBdy = bdy.size();
    for (index_t i = 0; i<nIntf; ++i)
    {
        _setupDofMapHelper(nPatches+i, nPatches+nIntf+nBdy, intf[i].first(), mb, dm, dof_map);
        _setupDofMapHelper(nPatches+i, nPatches+nIntf+nBdy, intf[i].second(), mb, dm, dof_map);
    }
    // boundary
    for (index_t i = 0; i<nBdy; ++i)
    {
        _setupDofMapHelper(nPatches+nIntf+i, nPatches+nIntf+nBdy, bdy[i], mb, dm, dof_map);
    }
    // interiors
    index_t max_nDofs = 0;
    for (index_t i=0; i<nPatches; ++i)
    {
        const index_t nDofs = mb[i].size();
        if (nDofs>max_nDofs) max_nDofs=nDofs;

        for (index_t j=0; j<nDofs; ++j)
        {
            const index_t idx = dm.index(j,i);
            if (dm.is_free_index(idx) && dof_map[idx] == -2 )
                dof_map[idx] = i;
        }
    }

    // Setup transfers
    const index_t nTransfers = nPatches+nIntf+nBdy+1;

    std::vector< gsSparseEntries<T> > transfers_se(nTransfers);
    for (index_t i=0; i<nTransfers; ++i)
        transfers_se[i].reserve(3*math::sqrt(max_nDofs));
    std::vector<index_t> counter(nTransfers,0);
    for (index_t j=0; j<nTotalDofs;++j)
    {
        const index_t i = dof_map[j];
        transfers_se[i].add( counter[i], j, 1 );
        counter[i]++;
    }

    std::vector< std::vector< gsSparseMatrix<T> > > transfers(3);
    transfers[0].reserve(nPatches);
    transfers[1].reserve(nIntf+nBdy);
    transfers[2].reserve(1);
    for (index_t i=0; i<nTransfers; ++i)
    {
        index_t d;
        if (i<nPatches)
            d = 0;
        else if (i<nPatches+nIntf+nBdy)
            d = 1;
        else
            d = 2;

        if (counter[i]>0)
        {
            gsSparseMatrix<T> transfer( counter[i], nTotalDofs );
            transfer.setFrom(transfers_se[i]);
            transfer.makeCompressed();
            transfers[d].push_back(give(transfer));
        }
        else
        {
            GISMO_ENSURE(d>0, "Empty patches are forbidden");
        }
    }

    return transfers;
}*/


template<typename T>
std::pair< std::vector< gsSparseMatrix<T> >, std::vector< typename gsLinearOperator<T>::Ptr > > setupPiecewisePreconditioner(
                                    gsSparseMatrix<T> A,
                                    const std::vector< typename gsLinearOperator<T>::Ptr >& ops_in,
                                    const gsMultiBasis<T>& mb,
                                    const gsBoundaryConditions<T>& bc,
                                    const gsOptionList& opt
                                )
{
    std::vector< std::vector< gsSparseMatrix<T> > > transfers_in = getPiecewiseTransfers(mb, bc, opt);
    GISMO_ASSERT( transfers_in.size() == static_cast<size_t>(mb.dim()+1), "getPiecewiseTransfers did not do what expected." );
    GISMO_ASSERT( transfers_in[0].size() == ops_in.size(), "Number of preconditioners does not agree with number of interiors." );

    std::vector< typename gsLinearOperator<T>::Ptr > ops;
    std::vector< gsSparseMatrix<T> > transfers;

    ops.reserve( transfers_in[0].size() + transfers_in[1].size() + transfers_in[2].size() );
    transfers.reserve( transfers_in[0].size() + transfers_in[1].size() + transfers_in[2].size() );

    for (size_t i=0; i<transfers_in[0].size(); ++i)
    {
        ops.push_back(ops_in[i]);
        transfers.push_back(give(transfers_in[0][i]));
    }
    for (size_t d=1; d<transfers_in.size(); ++d)
        for (size_t i=0; i<transfers_in[d].size(); ++i)
        {
            gsSparseMatrix<T> local_mat = transfers_in[d][i] * A * transfers_in[d][i].transpose();
            ops.push_back(makeSparseCholeskySolver(local_mat));
            transfers.push_back(give(transfers_in[d][i]));
        }
    return std::pair< std::vector< gsSparseMatrix<T> >, std::vector< typename gsLinearOperator<T>::Ptr > >(transfers, ops);
}


template<typename T>
std::vector<typename gsLinearOperator<T>::Ptr> getLocalExactSolvers(
                            const gsSparseMatrix<T>& A,
                            const std::vector< gsSparseMatrix<T> >& transfers
                            )
{
    std::vector<typename gsLinearOperator<T>::Ptr> ops;
    const size_t n=transfers.size();
    for (size_t i=0; i<n; ++i)
    {
        gsSparseMatrix<T> local_mat = transfers[i] * A * transfers[i].transpose();
        ops.push_back(makeSparseCholeskySolver(local_mat));
    }
    return ops;
}

namespace {
bool removeInvalids( gsVector<index_t>& vec )
{
    GISMO_ASSERT( vec.cols() == 1, "Only works for colum-vectors." );
    const index_t sz = vec.rows();
    index_t i = 0;
    index_t j = 0;
    while (j<sz)
    {
        if ( vec[j] != (index_t)(-1) )
        {
            vec[i] = vec[j];
            i++;
        }
        j++;
    }
    bool found = i>0;
    while (i<sz) { vec[i] = -1; i++; }
    return found;
}

template<typename T>
struct Piece {
    gsVector<index_t> indices;
    typename gsBasis<T>::Ptr basis;
// We want to have move semantics only
#if EIGEN_HAS_RVALUE_REFERENCES
    Piece(const Piece&) = delete;
    Piece(Piece&&) = default;
    Piece& operator=(const Piece&) = delete;
    Piece& operator=(Piece&&) = default;
#endif
    void swap(Piece& other) { indices.swap(other.indices); basis.swap(other.basis);  }
    Piece() {} // needed for give in C++98
    Piece( gsVector<index_t> i, typename gsBasis<T>::Ptr b ) { i.swap(indices); b.swap(basis); }
};

template<typename T>
std::vector< Piece<T> > decomposeIntoPieces( typename gsBasis<T>::Ptr basis, gsVector<index_t> all )
{
    std::vector< Piece<T> > result;
    const index_t d = basis ? basis->dim() : 0;

    result.reserve(2*d);

    for (index_t i=0;i<2*d;++i)
    {
        boxSide bs(i/2,i%2);
        std::vector< Piece<T> > tmp;
        if(d>1)
            tmp = decomposeIntoPieces<T>( basis->boundaryBasis(bs), basis->boundary(bs) );
        else
            tmp = decomposeIntoPieces<T>(typename gsBasis<T>::Ptr(), basis->boundary(bs) );
        const index_t sz = tmp.size();

        for (index_t j=0;j<sz;++j)
        {
            for (index_t l=0; l<tmp[j].indices.rows(); ++l)
            {
                const index_t loc = tmp[j].indices[l];
                if (loc != -1)
                {
                    tmp[j].indices[l] = all[loc] != -1 ? all[loc] : -1 ;
                    all[loc] = -1;
                }
            }
            if (removeInvalids(tmp[j].indices))
                result.push_back(give(tmp[j]));
        }

    }
    if (removeInvalids(all))
        result.push_back( Piece<T>(give(all), give(basis) ) );
    return result;
}

template<typename T>
gsSparseMatrix<T> setupTransfer( index_t cols, const gsVector<index_t>& indices )
{
    const index_t sz0 = indices.rows();
    index_t sz=0;
    for (; sz<sz0 && indices[sz] != -1; ++sz) {}
    gsSparseEntries<T> se;
    se.reserve(sz);
    for (index_t i=0; i<sz; ++i)
        se.add(i,indices[i],1);
    gsSparseMatrix<T> result(sz,cols);
    result.setFrom(se);
    return result;
}

template<typename T>
bool first_is_lower( const Piece<T>& a, const Piece<T>& b )
{ return a.indices[0]<b.indices[0]; }

template<typename T>
bool first_is_equal( const Piece<T>& a, const Piece<T>& b )
{ return a.indices[0]==b.indices[0]; }

} // anonymous namespace

template<typename T>
std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<T> > > > constructPieces(
                                    const gsMultiBasis<T>& mb,
                                    const std::vector< gsVector<index_t> >& locals,
                                    index_t totalNumberDof,
                                    bool combineVertices
                                )
{
    typedef typename gsBasis<T>::Ptr BasisPtr;

    const index_t dim = mb.dim();

    std::vector< Piece<T> > pd;

    for (size_t i=0;i<mb.nBases();++i)
    {
        typename gsBasis<T>::Ptr basis = mb[i].clone();
        std::vector< Piece<T> > tmp = decomposeIntoPieces<T>(basis, locals[i]);
        //pd.insert( pd.end(), tmp.begin(), tmp.end() ); // would make copy, so we do...
        for( typename std::vector< Piece<T> >::iterator it=tmp.begin(); it<tmp.end(); ++it)
            pd.push_back(give(*it));
    }
    std::sort(pd.begin(), pd.end(), first_is_lower<T>);
    typename std::vector< Piece<T> >::iterator it = std::unique(pd.begin(), pd.end(), first_is_equal<T>);
    pd.erase(it, pd.end());

    std::vector< std::vector< std::pair< BasisPtr, gsSparseMatrix<T> > > > result(dim+1);

    const index_t sz = pd.size();

    if (combineVertices)
    {
        gsVector<index_t> vertices(sz);
        index_t counter = 0;
        for (index_t i=0; i<sz; ++i)
        {
            const index_t d = pd[i].basis ? pd[i].basis->dim() : 0;
            GISMO_ASSERT( d<=dim, "Internal error." );
            if (d>0)
            {
                result[d].push_back(
                    std::pair< BasisPtr, gsSparseMatrix<T> >(
                        pd[i].basis,
                        setupTransfer<T>(totalNumberDof, pd[i].indices)
                    )
                );
            }
            else
            {
                GISMO_ASSERT( pd[i].indices.rows()==1 || (pd[i].indices.rows()>1&&pd[i].indices[2]==-1), "Internal error." );
                vertices[counter] = pd[i].indices[0];
                counter++;
            }
        }
        if (counter)
        {
            vertices[counter] = -1;
            result[0].push_back(
                std::pair< BasisPtr, gsSparseMatrix<T> >(
                    BasisPtr(),
                    setupTransfer<T>(totalNumberDof, vertices)
                )
            );
        }

    }
    else
    {
        for (index_t i=0; i<sz; ++i)
        {
            const index_t d = pd[i].basis ? pd[i].basis->dim() : 0;
            GISMO_ASSERT( d<=dim, "Internal error." );
            result[d].push_back(
                std::pair< BasisPtr, gsSparseMatrix<T> >(
                    pd[i].basis,
                    setupTransfer<T>(totalNumberDof, pd[i].indices)
                )
            );
        }
    }
    return result;

}


template<typename T>
std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<T> > > > constructPieces(
                                    const gsMultiBasis<T>& mb,
                                    const gsDofMapper& dm,
                                    bool combineVertices
                                )
{
    std::vector< gsVector<index_t> > locals(mb.nBases());
    for (size_t i=0;i<mb.nBases();++i)
    {
        const index_t sz = mb[i].size();
        gsVector<index_t> local(sz);
        for (index_t j=0; j<sz; ++j)
            local[j] = dm.is_free(j,i) ? dm.index(j, i) : -1;
        locals[i].swap(local);
    }

    return constructPieces( mb, locals, dm.freeSize(), combineVertices );

}

template<typename T>
std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<T> > > > constructPieces(
                                    const gsMultiBasis<T>& mb,
                                    const gsBoundaryConditions<T>& bc,
                                    const gsOptionList& opt,
                                    bool combineVertices
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

    return constructPieces( mb, dm, combineVertices );

}

} // namespace gismo
