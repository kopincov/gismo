/** @file gsParallelMultiPatchSmoothers.cpp

    @brief Provides Multigrid smoothers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#include <gsMultiGrid/gsParallelMultiPatchPreconditioners.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsCore/gsMultiBasis.h>
#include <gsUtils/gsSortedVector.h>
#include <gsIETI/gsParallelOperator.h>

namespace gismo
{

namespace {
template <typename T>
memory::shared_ptr< gsSparseMatrix<T> > getSparseMatrixPtr( index_t rows, index_t cols, const gsSparseEntries<T>& se )
{
    memory::shared_ptr< gsSparseMatrix<T> > result( new gsSparseMatrix<T>( rows, cols ) );
    result->setFrom(se);
    result->makeCompressed();
    return result;
}
}

template <typename T>
void gsParallelAdditivePreconditionerOp<T>::step(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_A->rows() == x.rows() && x.rows() == f.rows() && m_A->cols() == m_A->rows() && x.cols() == f.cols(),
                  "The dimensions do not fit." );

    const index_t n = m_ops.size();
    m_A->apply( x, m_res );
    m_res -= f;

    m_corr_global.setZero(x.rows(), x.cols());
    m_corr_global_buffer.setZero(x.rows(), x.cols());
    for (index_t i=0; i<n; ++i)
    {
        m_restriction[i]->apply(m_res, m_local_res);
        m_ops[i]->apply(m_local_res, m_corr);
        m_corr *= m_damping;
        m_prolongation[i]->apply(m_corr, m_corr_global);
        m_corr_global_buffer+=m_corr_global;
    }
    m_A->accumulate(m_corr_global_buffer,m_corr_global);
    x -= m_corr_global;
}

/*
template <typename T>
std::vector< gsSparseMatrix<T> > getPatchwiseTransfers(const gsMultiBasis<T>& mb,const gsBoundaryConditions<T>& bc,const gsOptionList& opt,
                                                       const gsSortedVector<size_t>& myPatches,const gsMpiComm comm)
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
*/


namespace {
template <typename T>
void _setupDofMapHelper(index_t i, patchSide patch_side, const gsMultiBasis<T>& mb, const gsDofMapper& dm, std::vector<index_t>& dof_map)
{
    index_t patch_index = patch_side.patch;

    gsMatrix<unsigned> bdy = mb.basis(patch_index).boundary(patch_side);
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
                dof_map[idx] = -1;
        }
    }

}
}

template <typename T>
std::vector< std::vector< gsSparseMatrix<T> > > getPiecewiseTransfers(const gsMultiPatch<T>& mp,const gsMultiBasis<T>& mb,
                                                                      const gsBoundaryConditions<T>& bc,const gsOptionList& opt,
                                                                      const gsSortedVector<size_t>& myPatches,
                                                                      const gsParallelGlobalLocalHandler& A_handler,
                                                                      const gsMpiComm comm, gsSortedVector<index_t>& interfaceDofs)
{
   // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  start.\n"<<std::flush;
    std::vector< std::vector< gsSparseMatrix<T> > > transfers(3);

    gsDofMapper dm;
    mb.getMapper(
                (dirichlet::strategy)opt.askInt("DirichletStrategy",11),
                (iFace    ::strategy)opt.askInt("InterfaceStrategy", 1),
                bc,
                dm,
                0
                );
    const index_t nTotalDofs = dm.freeSize();
    const size_t nProcLocalDofs = A_handler.localSize();
    std::vector<index_t> dof_map(nTotalDofs,-2);
  //  gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  got mapper.\n"<<std::flush;
    //get ProcLocal indices
    gsMatrix<index_t> inp(A_handler.localSize(),1);
    inp.col(0) = gsVector<index_t>::LinSpaced(A_handler.localSize(),0,A_handler.localSize()-1);
    gsMatrix<index_t> inp2(A_handler.globalSize(),1);
    inp2.col(0) = gsVector<index_t>::LinSpaced(A_handler.globalSize(),0,A_handler.globalSize()-1);
    gsMatrix<index_t> globToLoc,locToGlob;
    globToLoc.setZero(A_handler.globalSize(),1);
    locToGlob.setZero(A_handler.localSize(),1);
    A_handler.addLocalVectorToGlobal(inp,globToLoc);
    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  got globToLoc."<<globToLoc.transpose()<<"\n"<<std::flush;
    A_handler.extractLocalVector(inp2,locToGlob);
    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  got locToGlob."<<locToGlob.transpose()<<"\n"<<std::flush;


    std::vector<boundaryInterface> intf = mp.interfaces();
    const index_t nIntf = intf.size();
    index_t nProcLocIntf =0;
    for (index_t i = 0; i<nIntf; ++i)
    {
        if(myPatches.bContains(intf[i].first().patch)||myPatches.bContains(intf[i].second().patch))
        {
            _setupDofMapHelper(nProcLocIntf, intf[i].first(), mb, dm, dof_map);
            _setupDofMapHelper(nProcLocIntf, intf[i].second(), mb, dm, dof_map);
            nProcLocIntf++;
        }
    }
    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  setup interfaces."+ util::to_string(nIntf)+"|"+util::to_string(nProcLocIntf)+"\n"<<std::flush;

    std::vector<patchSide> bdy = mp.boundaries();
    const index_t nBdy = bdy.size();
    index_t nProcLocBdy =0;
    for (index_t i = 0; i<nBdy; ++i)
    {
        if(myPatches.bContains(bdy[i].patch))
        {
            _setupDofMapHelper(nProcLocIntf+nProcLocBdy, bdy[i], mb, dm, dof_map);
            nProcLocBdy++;
        }
    }
    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  setup boundary. "+ util::to_string(nBdy)+"|"+util::to_string(nProcLocBdy)+"\n"<<std::flush;
    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  got transfer for interfaces "<<gsAsVector<index_t>(dof_map).transpose()<<std::flush;


/*
    std::vector<index_t> corners_;
    corners_.reserve(myPatches.size()*math::exp2(mb.dim()));
    for(index_t i = 0; i< dof_map.size();++i)
        if(dof_map[i]==-1) corners_.push_back(i);
    gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  setup corners.\n"<<std::flush;

    std::vector<index_t> c_size(comm.size());
    c_size[comm.rank()]=(index_t)corners_.size();
    comm.allgather(MPI_IN_PLACE,1,c_size.data());

    std::vector<index_t> disp(comm.size());
    size_t sum = c_size[0];
    disp[0]=0;
    for(index_t p =1; p<comm.size();++p)
    {
        disp[p] = disp[p-1]+c_size[p-1];
        sum+=c_size[p];
    }

    std::vector<index_t> corners(sum);
    comm.allgatherv(corners_.data(),corners_.size(),corners.data(),c_size.data(),disp.data());
    for(size_t i=0; i<corners.size();++i)
        dof_map[corners[i]]=-1; //maybe this is not needed!!! TODO
    */

    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  finished corners.\n"<<std::flush;
    index_t max_nDofs = 0;
    for (size_t ii=0; ii<myPatches.size(); ++ii)
    {
        index_t i = myPatches[ii];
        const index_t nDofs = mb[i].size();
        if (nDofs>max_nDofs) max_nDofs=nDofs;
        gsSparseEntries<T> transfer_se;
        transfer_se.reserve(nDofs);
        index_t counter = 0;
        for (index_t j=0; j<nDofs; ++j)
        {
            const index_t idx = dm.index(j,i);
            if (dm.is_free_index(idx) && dof_map[idx] == -2 )
            {
                transfer_se.add(counter,globToLoc(idx,0),1);
                counter++;
            }
        }
        gsSparseMatrix<T> transfer(counter, nProcLocalDofs);
        transfer.setFrom(transfer_se);
        transfer.makeCompressed();

        transfers[0].push_back(give(transfer));
    }

    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  got transfer for interior.\n"<<std::flush;


    const index_t nInterfaces = nProcLocIntf+nProcLocBdy;

    std::vector< gsSparseEntries<T> > transfers_se(nInterfaces+1);
    for (index_t i=-1; i<nInterfaces; ++i)
        transfers_se[i+1].reserve(3*math::sqrt(max_nDofs));
    std::vector<index_t> counter(nInterfaces+1,0);

    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  setup SE.\n"<<std::flush;

    size_t nIntfDofs = 0;
    for (size_t j=0; j<nProcLocalDofs;++j)
    {
        const index_t i = dof_map[locToGlob(j,0)];
        if ( i>-2 )
        {
            transfers_se[i+1].add( counter[i+1], j, 1 );
            counter[i+1]++;
            nIntfDofs++;
        }
    }
    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  got transfer for interfaces.\n"<<std::flush;
    interfaceDofs.reserve(nIntfDofs);
    for (size_t j=0; j<nProcLocalDofs;++j)
        if ( dof_map[locToGlob(j,0)]>-2 )
            interfaceDofs.push_sorted_unique(j);

    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  got interface dofs.\n"<<std::flush;

    for (index_t i=-1; i<nInterfaces; ++i)
    {
        if (counter[i+1]>0)
        {
            gsSparseMatrix<T> transfer( counter[i+1], nProcLocalDofs );
            transfer.setFrom(transfers_se[i+1]);
            transfer.makeCompressed();

            if (i==-1)
                transfers[2].push_back(give(transfer));
            else
                transfers[1].push_back(give(transfer));
        }
    }
    return transfers;
}


template <typename T>
std::pair< std::vector< gsSparseMatrix<T> >, std::vector< typename gsLinearOperator<T>::Ptr > > setupPiecewisePreconditioner(
        const gsParallelOperator<T>& A,
        const gsParallelGlobalLocalHandler& A_handler,
        const std::vector< typename gsLinearOperator<T>::Ptr >& ops_in,
        const gsMultiPatch<T>& mp,
        const gsMultiBasis<T>& mb,
        const gsBoundaryConditions<T>& bc,
        const gsOptionList& opt,
        const gsSortedVector<size_t>& myPatches,
        const gsMpiComm comm
        )
{
    gsSortedVector<index_t> interfaceDofs;
    std::vector< std::vector< gsSparseMatrix<T> > > transfers_in = getPiecewiseTransfers(mp, mb, bc, opt, myPatches,A_handler, comm,interfaceDofs);

 //   gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  got transfer.\n"<<std::flush;

    std::vector<typename gsLinearOperator<T>::Ptr > opsInput(myPatches.size());
    if(ops_in.size() == mb.nBases())
    {
        for(size_t i=0; i<myPatches.size();++i)
            opsInput[i] = ops_in[myPatches[i]];
    }
    else
        opsInput = ops_in;

    GISMO_ASSERT( transfers_in.size() == 3, "getPiecewiseTransfers did not do what expected." );
    GISMO_ASSERT( transfers_in[0].size() == opsInput.size(), "Number of preconditioners does not agree with number of interiors." );

    std::vector< typename gsLinearOperator<T>::Ptr > ops;
    std::vector< gsSparseMatrix<T> > transfers;

    ops.reserve( transfers_in[0].size() + transfers_in[1].size() + transfers_in[2].size() );
    transfers.reserve( transfers_in[0].size() + transfers_in[1].size() + transfers_in[2].size() );

    gsSparseMatrix<T> intfA;
    const gsPatchSubassambledLocalOperator<T>* A_ptr = dynamic_cast<const gsPatchSubassambledLocalOperator<T>*>(A.getLinearOperator().get());
    const gsMatrixOp<gsSparseMatrix<T> >* A_ptr2 = dynamic_cast<const gsMatrixOp<gsSparseMatrix<T> >*>(A.getLinearOperator().get());
    GISMO_ASSERT(A_ptr!=NULL || A_ptr2 !=NULL, "You need a PatchSubassembledLocalOperator or a gsMatrixOp.");

    if(A_ptr!=NULL)
        intfA = A_ptr->buildSparseMatrix(interfaceDofs);
    else if(A_ptr2!=NULL)
        intfA = A_ptr2->matrix(); //TODO: optimize by not copying, e.g., accumulate A and distribute at the end again.


 //   gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  got restricted matrix.\n"<<std::flush;
/*
    for(int i=0; i<comm.size();++i)
    {
        comm.barrier();
        if(comm.rank()==i)
            gsInfo<<"restricted A:\n"<<intfA.toDense()<<"\n\n";
    }
*/
    A.getConnectionHandlerTestSpace()->accumulateSparseMatrix(intfA);
    /*
    gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  accumulated matrix.\n"<<std::flush;
    for(int i=0; i<comm.size();++i)
    {
        comm.barrier();
        if(comm.rank()==i)
            gsInfo<<"restricted A:\n"<<intfA.toDense()<<"\n\n";
    }
*/


    for (size_t i=0; i<transfers_in[0].size(); ++i)
    {
        ops.push_back(opsInput[i]);
        transfers.push_back(give(transfers_in[0][i]));
    }
    /*
    for(int i=0; i<comm.size();++i)
    {
        comm.barrier();
        if(comm.rank()==i)
        {
            for (size_t i=0; i<transfers_in[1].size(); ++i)
            {
                gsSparseMatrix<T> local_mat = transfers_in[1][i]* intfA * transfers_in[1][i].transpose();
                gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  accumulated matrix #"+util::to_string(i)+"  \n"<<local_mat.toDense()<<"\n"<<std::flush;
            }
            for (size_t i=0; i<transfers_in[2].size(); ++i)
            {
                gsSparseMatrix<T> local_mat = transfers_in[2][i]* intfA * transfers_in[2][i].transpose();
                gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  accumulated matrix 2#"+util::to_string(i)+"  \n"<<local_mat.toDense()<<"\n"<<std::flush;
            }

        }
    }
    */

    for (size_t i=0; i<transfers_in[1].size(); ++i)
    {
        gsSparseMatrix<T> local_mat = transfers_in[1][i]* intfA * transfers_in[1][i].transpose();

        ops.push_back(makeSparseCholeskySolver(local_mat));
        transfers.push_back(give(transfers_in[1][i]));
    }
    for (size_t i=0; i<transfers_in[2].size(); ++i)
    {
        gsSparseMatrix<T> local_mat = transfers_in[2][i] * intfA * transfers_in[2][i].transpose();
        // For 2D, local_mat is diagonal, so we could do something smarter...
        // But for 3D, this contains all vertices _and_ edges, so here it's not diagonal...
        ops.push_back(makeSparseCholeskySolver(local_mat));
        transfers.push_back(give(transfers_in[2][i]));
    }
    return std::pair< std::vector< gsSparseMatrix<T> >, std::vector<typename gsLinearOperator<T>::Ptr > >(transfers, ops);
}

/*
template <typename T>
std::vector<gsLinearOperator<T>::Ptr> getLocalExactSolvers(const gsSparseMatrix<T>& A, const std::vector< gsSparseMatrix<T> >& transfers)
{
    std::vector<gsLinearOperator<T>::Ptr> ops;
    const size_t n=transfers.size();
    for (size_t i=0; i<n; ++i)
    {
        gsSparseMatrix<T> local_mat = transfers[i] * A * transfers[i].transpose();
        ops.push_back(makeSparseCholeskySolver(local_mat));
    }
    return ops;
}
*/
} // namespace gismo
