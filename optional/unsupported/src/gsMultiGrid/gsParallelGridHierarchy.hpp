/**  gsParallelGridHierarchy.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on:  2018-02-23
*/

#pragma once

#include <gsSolver/gsMatrixOp.h>
#include <gsIETI/gsParallelOperator.h>
#include <gsMultiGrid/gsParallelGridHierarchy.h>
#include <gsAssembler/gsAssembler.h>


namespace gismo {

template<typename T>
gsParallelGridHierarchy<T> gsParallelGridHierarchy<T>::buildByRefinement(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& options,
        gsSortedVector<size_t> myPatches,
        gsMpiComm comm,
        index_t levels,
        index_t numberOfKnotsToBeInserted,
        index_t multiplicityOfKnotsToBeInserted
        )
{
    bool noSubassembledOps = options.askSwitch("NoSubassembledOperators",false);
    gsParallelGridHierarchy<T> result;
    result.m_boundaryConditions = boundaryConditions;
    result.m_options = options;
    result.m_mBases.resize(levels);
    result.m_transferMatrices.resize(levels-1);
    result.m_mBases[0] = give(mBasis);
    result.m_myPatches = give(myPatches);
    result.m_comm = comm;
    result.m_subassTopology.resize(levels);
    result.m_connectionHandlers.resize(levels);
    result.m_globLocHandlers.resize(levels);



    size_t nPatches = result.m_mBases.front().nBases();

    std::vector<gsDofMapper> localMappers(nPatches);
    gsDofMapper globalMapper;

    std::vector<gsBoundaryConditions<real_t> > patchBC(nPatches);
    for(size_t j=0; j<nPatches;++j)
        boundaryConditions.getConditionsForPatch(j,patchBC[j]);

    //Get initial coarse mappers
    result.m_mBases[0].getMapper(
                (dirichlet::strategy)options.askInt("DirichletStrategy",11),
                (iFace    ::strategy)options.askInt("InterfaceStrategy", 1),
                boundaryConditions,
                globalMapper,
                0
                );

    //all local mappers are required
    for(size_t j=0; j<nPatches;++j)
        gsMultiBasis<T>(result.m_mBases[0][j]).getMapper(
                    (dirichlet::strategy)options.askInt("DirichletStrategy",11),
                    iFace::none,
                    patchBC[j],
                    localMappers[j],
                    0
                    );

    result.m_coarseGlobMapper = globalMapper;
    result.m_coarseLocMappers = localMappers;

    result.m_subassTopology[0] = gsPatchSubassembledTopology<T>::make(result.m_myPatches,result.m_mBases.front(),localMappers,globalMapper);
    result.m_subassTopology[0]->reorderLike(globalMapper);

    gsPatchInterfaceConnections<real_t> patchConnection(result.m_myPatches,result.m_mBases[0],localMappers,result.m_subassTopology[0]->getLocalPatchMapper(),globalMapper,comm);
    patchConnection.init();

    result.m_connectionHandlers[0] = gsConnectionHandler<T>::make(patchConnection.getConnectionPairs(),comm);
    std::vector<gsDofMapper> procLocalMapper(comm.size());
    procLocalMapper[comm.rank()] = result.m_subassTopology[0]->getLocalMapper();
    result.m_globLocHandlers[0] = gsParallelGlobalLocalHandler::make(patchConnection.generateProcGlobalMapper(),procLocalMapper,comm);

    for ( index_t i=1; i<levels; ++i )
    {
        result.m_mBases[i] = result.m_mBases[i-1];
        result.m_transferMatrices[i-1].resize(nPatches);

        // Refine
        if(noSubassembledOps)
        {
            for (size_t np = 0; np < nPatches; ++np)
            {
                if(result.m_myPatches.bContains(np))
                    result.m_mBases[i][np].uniformRefine_withTransfer(result.m_transferMatrices[i-1][np],
                            numberOfKnotsToBeInserted,multiplicityOfKnotsToBeInserted);
                else
                    result.m_mBases[i][np].uniformRefine(numberOfKnotsToBeInserted,multiplicityOfKnotsToBeInserted);
            }
        }
        else
        {
            for (size_t np = 0; np < nPatches; ++np)
            {
                if(result.m_myPatches.bContains(np))
                    gsMultiBasis<T>(result.m_mBases[i][np]).uniformRefine_withTransfer(result.m_transferMatrices[i-1][np],
                            patchBC[np],options,numberOfKnotsToBeInserted,multiplicityOfKnotsToBeInserted);
                // else
                result.m_mBases[i][np].uniformRefine(numberOfKnotsToBeInserted,multiplicityOfKnotsToBeInserted);
            }
        }

        //Get initial coarse mappers
        result.m_mBases[i].getMapper(
                    (dirichlet::strategy)options.askInt("DirichletStrategy",11),
                    (iFace    ::strategy)options.askInt("InterfaceStrategy", 1),
                    boundaryConditions,
                    globalMapper,
                    0
                    );

        //all local mappers are required
        for(size_t j=0; j<nPatches;++j)
            gsMultiBasis<T>(result.m_mBases[i][j]).getMapper(
                        (dirichlet::strategy)options.askInt("DirichletStrategy",11),
                        iFace::none,
                        patchBC[j],
                        localMappers[j],
                        0
                        );



        result.m_subassTopology[i] = gsPatchSubassembledTopology<T>::make(result.m_myPatches,result.m_mBases[i],localMappers,globalMapper);
        result.m_subassTopology[i]->reorderLike(globalMapper);

        gsPatchInterfaceConnections<real_t> patchConnection(result.m_myPatches,result.m_mBases[i],localMappers,result.m_subassTopology[i]->getLocalPatchMapper(),globalMapper,comm);
        patchConnection.init();

        result.m_connectionHandlers[i] = gsConnectionHandler<T>::make(patchConnection.getConnectionPairs(),comm);

        procLocalMapper[comm.rank()] = result.m_subassTopology[i]->getLocalMapper();
        result.m_globLocHandlers[i] = gsParallelGlobalLocalHandler::make(patchConnection.generateProcGlobalMapper(),procLocalMapper,comm);

    }
    return result;
}


template <typename T>
gsParallelGridHierarchy<T> gsParallelGridHierarchy<T>::buildByCoarsening(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& options,
        gsSortedVector<size_t> myPatches,
        gsMpiComm comm,
        index_t levels,
        index_t degreesOfFreedom
        )
{

    bool noSubassembledOps = options.askSwitch("NoSubassembledOperators",false);
    gsParallelGridHierarchy<T> result;
    result.m_boundaryConditions = boundaryConditions;
    result.m_options = options;
    result.m_myPatches = give(myPatches);
    result.m_comm = comm;
    result.m_transferMatrices.reserve(levels-1);
    result.m_mBases.reserve(levels);
    result.m_subassTopology.reserve(levels);
    result.m_connectionHandlers.reserve(levels);
    result.m_globLocHandlers.reserve(levels);

    result.m_mBases.push_back(give(mBasis));


    index_t lastSize = result.m_mBases[0].totalSize();
    size_t nPatches = result.m_mBases.front().nBases();

    std::vector<gsDofMapper> localMappers(nPatches);
    gsDofMapper globalMapper;

    std::vector<gsBoundaryConditions<real_t> > patchBC(nPatches);
    for(size_t j=0; j<nPatches;++j)
        boundaryConditions.getConditionsForPatch(j,patchBC[j]);

    //Get initial coarse mappers
    result.m_mBases[0].getMapper(
                (dirichlet::strategy)options.askInt("DirichletStrategy",11),
                (iFace    ::strategy)options.askInt("InterfaceStrategy", 1),
                boundaryConditions,
                globalMapper,
                0
                );

    //all local mappers are required
    for(size_t j=0; j<nPatches;++j)
        gsMultiBasis<T>(result.m_mBases[0][j]).getMapper(
                    (dirichlet::strategy)options.askInt("DirichletStrategy",11),
                    iFace::none,
                    patchBC[j],
                    localMappers[j],
                    0
                    );


    result.m_subassTopology.push_back(gsPatchSubassembledTopology<T>::make(result.m_myPatches,result.m_mBases.front(),localMappers,globalMapper));
    result.m_subassTopology.front()->reorderLike(globalMapper);


    gsPatchInterfaceConnections<real_t> patchConnection(result.m_myPatches,result.m_mBases.front(),localMappers,result.m_subassTopology[0]->getLocalPatchMapper(),globalMapper,comm);
    patchConnection.init();

    result.m_connectionHandlers.push_back(gsConnectionHandler<T>::make(patchConnection.getConnectionPairs(),comm));
    std::vector<gsDofMapper> procLocalMapper(comm.size());
    procLocalMapper[comm.rank()] = result.m_subassTopology[0]->getLocalMapper();
    result.m_globLocHandlers.push_back(gsParallelGlobalLocalHandler::make(patchConnection.generateProcGlobalMapper(),procLocalMapper,comm));


    for (int i = 0; i < levels -1 && lastSize > degreesOfFreedom; ++i)
    {
        result.m_transferMatrices.push_back(std::vector<gsSparseMatrix<T, RowMajor> >(nPatches));
        gsMultiBasis<T> coarseMBasis = result.m_mBases[i];

        // Coarsen
        if(noSubassembledOps)
        {
            for (size_t np = 0; np < nPatches; ++np)
            {
                if(result.m_myPatches.bContains(np))
                    coarseMBasis.basis(np).uniformCoarsen_withTransfer(result.m_transferMatrices[i][np]);
                else
                    coarseMBasis.basis(np).uniformCoarsen();
            }
        }
        else
        {
            for (size_t np = 0; np < nPatches; ++np)
            {
                if(result.m_myPatches.bContains(np))
                    gsMultiBasis<T>(coarseMBasis.basis(np)).uniformCoarsen_withTransfer(result.m_transferMatrices[i][np],patchBC[np],options);
                //       else
                coarseMBasis.basis(np).uniformCoarsen();
            }
        }


        //Get initial coarse mappers
        coarseMBasis.getMapper(
                    (dirichlet::strategy)options.askInt("DirichletStrategy",11),
                    (iFace    ::strategy)options.askInt("InterfaceStrategy", 1),
                    boundaryConditions,
                    globalMapper,
                    0
                    );

        //all local mappers are required
        for(size_t j=0; j<nPatches;++j)
            gsMultiBasis<T>(coarseMBasis[j]).getMapper(
                        (dirichlet::strategy)options.askInt("DirichletStrategy",11),
                        iFace::none,
                        patchBC[j],
                        localMappers[j],
                        0
                        );

        result.m_subassTopology.push_back(gsPatchSubassembledTopology<T>::make(result.m_myPatches,coarseMBasis,localMappers,globalMapper));
        result.m_subassTopology.back()->reorderLike(globalMapper);

        gsPatchInterfaceConnections<real_t> patchConnection(result.m_myPatches,coarseMBasis,localMappers,result.m_subassTopology.back()->getLocalPatchMapper(),globalMapper,comm);
        patchConnection.init();

        result.m_connectionHandlers.push_back(gsConnectionHandler<T>::make(patchConnection.getConnectionPairs(),comm));

        procLocalMapper[comm.rank()] = result.m_subassTopology.back()->getLocalMapper();
        result.m_globLocHandlers.push_back(gsParallelGlobalLocalHandler::make(patchConnection.generateProcGlobalMapper(),procLocalMapper,comm));

        index_t newSize = coarseMBasis.totalSize();
        // If the number of dofs could not be decreased, then cancel. However, if only the number
        // of levels was specified, then this should be ignored (the caller might need to have a
        // fixed number of levels).
        if (lastSize <= newSize && degreesOfFreedom > 0)
            break;
        lastSize = newSize;

        result.m_mBases.push_back(give(coarseMBasis));
    }

    std::reverse( result.m_mBases.begin(), result.m_mBases.end() );
    std::reverse( result.m_transferMatrices.begin(), result.m_transferMatrices.end() );
    std::reverse( result.m_subassTopology.begin(), result.m_subassTopology.end() );
    std::reverse( result.m_connectionHandlers.begin(), result.m_connectionHandlers.end() );
    std::reverse( result.m_globLocHandlers.begin(), result.m_globLocHandlers.end() );

    result.m_coarseGlobMapper = globalMapper;
    result.m_coarseLocMappers = localMappers;

    return result;
}


template<typename T>
std::pair<std::pair<std::vector< typename gsParallelOperator<T>::Ptr >, std::vector< typename gsParallelOperator<T>::Ptr > > ,
std::vector<gsParallelGlobalLocalHandler::Ptr> >
gsParallelGridHierarchy<T>::getRestrictionAndProlongationOperators()
{
    std::pair<std::pair<std::vector< typename gsParallelOperator<T>::Ptr >, std::vector< typename gsParallelOperator<T>::Ptr > > ,
            std::vector<gsParallelGlobalLocalHandler::Ptr> > result;
    result.first.first.resize(m_mBases.size()-1);
    result.first.second.resize(m_mBases.size()-1);
    result.second.resize(m_mBases.size());

    std::vector<typename gsLinearOperator<T>::Ptr> localOpsT(m_mBases.front().nBases());
    std::vector<typename gsLinearOperator<T>::Ptr> localOps(m_mBases.front().nBases());

    bool noSubassembledOps = m_options.askSwitch("NoSubassembledOperators",false);

    if(noSubassembledOps)
        m_transferMatricesAss.resize(m_mBases.size()-1);

    for (size_t i=1; i<m_mBases.size(); ++i ) //levels
    {
        gsLinearOperator<real_t>::Ptr prolongationOp,restrictionOp;

        if(noSubassembledOps)
        {
            const gsMultiBasis<real_t>& mb = m_mBases[i];
            gsSparseMatrix<real_t,RowMajor> transfer;
            mb.combineTransferMatrices(m_transferMatrices[i-1],m_subassTopology[i-1]->getLocalPatchMapper(),m_subassTopology[i]->getLocalPatchMapper(),transfer);
            m_transferMatricesAss[i-1] = transfer.moveToPtr();
            prolongationOp = makeMatrixOp(m_transferMatricesAss[i-1]);
            restrictionOp = makeMatrixOp(m_transferMatricesAss[i-1]->transpose());
        }
        else
        {
            for(size_t j=0; j<m_myPatches.size();++j)
            {
                localOps[m_myPatches[j]] = makeMatrixOp(m_transferMatrices[i-1][m_myPatches[j]]);
                localOpsT[m_myPatches[j]] = makeMatrixOp(m_transferMatrices[i-1][m_myPatches[j]].transpose());
            }
            prolongationOp = gsPatchSubassambledLocalOperator<real_t,true,false>::make(give(localOps),m_subassTopology[i],m_subassTopology[i-1]);
            restrictionOp = gsPatchSubassambledLocalOperator<real_t,false,true>::make(give(localOpsT),m_subassTopology[i-1],m_subassTopology[i]);
        }

        result.first.second[i-1] = gsParallelOperator<real_t>::make(m_connectionHandlers[i],m_connectionHandlers[i-1],prolongationOp);
        result.first.first[i-1] = gsParallelOperator<real_t>::make(m_connectionHandlers[i-1],m_connectionHandlers[i],restrictionOp);


        result.second[i-1] = m_globLocHandlers[i-1];
        if(i == m_mBases.size()-1 )
            result.second[i] = m_globLocHandlers[i];

        localOps.clear(); localOpsT.clear();
        localOps.resize(m_mBases.front().nBases()); localOpsT.resize(m_mBases.front().nBases());
    }

    return result;
}

template<typename T>
std::pair<std::vector<typename gsParallelOperator<T>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > gsParallelGridHierarchy<T>::generateGalerkinProjection(const std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr >& A, std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr >& A_coarse, std::vector<gsDofMapper>& localMappers_c, gsDofMapper& globalMapper_c) const
{
    A_coarse.resize(m_mBases.front().nBases());

    std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr> localOps(m_mBases.front().nBases());
    std::vector<typename gsLinearOperator<T>::Ptr > localOps_(m_mBases.front().nBases());

    std::pair<std::vector<typename gsParallelOperator<T>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > result;
    result.first.resize(m_mBases.size());
    result.second.resize(m_mBases.size());


    for ( index_t i=m_mBases.size()-1; i>=0; --i ) //levels
    {
        if((size_t)i==m_mBases.size()-1)
        {
            localOps=A;
        }
        else
        {
            for( size_t k =0; k<m_myPatches.size();++k)
            {
                size_t np = m_myPatches[k];
                typename gsSparseMatrix<T>::Ptr X = gsSparseMatrix<T>(m_transferMatrices[i][np].transpose()*localOps[np]->matrix()*m_transferMatrices[i][np]).moveToPtr();

                // gsDebug<<"local tr\n"<<m_transferMatrices[i][np].toDense()<<"\n\n";
                typename gsMatrixOp<gsSparseMatrix<T> >::Ptr op = makeMatrixOp(X);
                localOps[np] = op;
                if(i==0) //coarsest grid
                    A_coarse[np] = op;
            }
        }

        //setup subassembled operators
        for( size_t k =0; k<m_myPatches.size();++k)
        {
            size_t np = m_myPatches[k];
            localOps_[np] = localOps[np]; //cast to gsLinearOperator<T>::Ptr
        }

        typename gsPatchSubassambledLocalOperator<T>::Ptr subassLocOp = gsPatchSubassambledLocalOperator<T>::make(localOps_,m_subassTopology[i]);
        result.first[i] = gsParallelOperator<T>::make(m_connectionHandlers[i],m_connectionHandlers[i],subassLocOp);
        result.second[i] = m_globLocHandlers[i];
    }

    localMappers_c = m_coarseLocMappers;
    globalMapper_c = m_coarseGlobMapper;


    return result;
}


template<typename T>
std::pair<std::vector<typename gsParallelOperator<T>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > gsParallelGridHierarchy<T>::generateGalerkinProjection(const typename gsMatrixOp<gsSparseMatrix<T> >::Ptr & A, typename gsMatrixOp<gsSparseMatrix<T> >::Ptr  & A_coarse, std::vector<gsDofMapper>& localMappers_c, gsDofMapper& globalMapper_c) const
{

    typename gsMatrixOp<gsSparseMatrix<T> >::Ptr localOps;

    std::pair<std::vector<typename gsParallelOperator<T>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > result;
    result.first.resize(m_mBases.size());
    result.second.resize(m_mBases.size());


    for ( index_t i=m_mBases.size()-1; i>=0; --i ) //levels
    {
        if((size_t)i==m_mBases.size()-1)
        {
            localOps=A;
        }
        else
        {
                typename gsSparseMatrix<T>::Ptr X =  gsSparseMatrix<T>(m_transferMatricesAss[i]->transpose()*localOps->matrix()**m_transferMatricesAss[i]).moveToPtr();

                // gsDebug<<"local tr\n"<<m_transferMatrices[i][np].toDense()<<"\n\n";
                localOps = makeMatrixOp(X);
                if(i==0) //coarsest grid
                    A_coarse = localOps;
        }

        result.first[i] = gsParallelOperator<T>::make(m_connectionHandlers[i],m_connectionHandlers[i],localOps);
        result.second[i] = m_globLocHandlers[i];
    }

    localMappers_c = m_coarseLocMappers;
    globalMapper_c = m_coarseGlobMapper;


    return result;
}

template<typename T>
typename gsCoarseSolverAdapter<T>::Ptr constructCoarseSolver(const std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr >& mat_coarse, const gsParallelGlobalLocalHandler::Ptr& handler_coarse, const gsSortedVector<size_t>& myPatches, const std::vector<gsDofMapper>& patchLocalMappers, const gsDofMapper& globalMapper )
{
    gsMpiComm comm = handler_coarse->getComm();
    //serialize the data
    //  std::vector<std::vector<T> > buffer(comm.size());
    //   std::vector<std::vector<index_t> >sizes(comm.size());
    //   sizes[comm.rank()].resize(myPatches.size());
    //   for(size_t k=0; k<myPatches.size();++k)
    //       sizes[comm.rank()][k] = mat_coarse[np].rows()*mat_coarse[np].cols();

    std::vector<size_t> sizes(mat_coarse.size());
    for(size_t k=0; k<myPatches.size();++k)
        sizes[myPatches[k]] = mat_coarse[myPatches[k]]->rows()*mat_coarse[myPatches[k]]->cols();

    comm.sum(sizes.data(),sizes.size());
    size_t sum=0;
    for(size_t i=0; i<sizes.size();++i)
        sum+=sizes[i];

    gsMatrix<T> buffer;
    buffer.setZero(sum,1);
    size_t shift =0;
    for(size_t np=0; np< mat_coarse.size();++np)
    {
        if(myPatches.bContains(np))
        {
            /*
           gsMatrix<T> mat = mat_coarse[np].matrix().toDense().resize(mat.cols()*mat.rows(),1);
               for(index_t r =0; r<mat.rows();++r)
                   buffer[shift+r]=mat(r,0);
                   */
            gsMatrix<T> temp = mat_coarse[np]->matrix().toDense();
            temp.resize(mat_coarse[np]->cols()*mat_coarse[np]->rows(),1);
            buffer.block(shift,0,sizes[np],1) = give(temp);
        }
        shift+=sizes[np];

    }
    comm.sum(buffer.data(), buffer.rows()); //TODO: could be optimized using allgatherv

    std::vector<gsMatrix<T> > mats(mat_coarse.size());
    shift = 0;
    for(size_t np=0; np< mat_coarse.size();++np)
    {
        gsMatrix<T> temp = buffer.block(shift,0,sizes[np],1);
        temp.resize(math::sqrt(sizes[np]),math::sqrt(sizes[np]));
        mats[np] =give(temp);
        shift+=sizes[np];
    }

    gsSparseMatrix<T> coarseMat(handler_coarse->globalSize(),handler_coarse->globalSize());
    gsSparseEntries<T> entries;
    entries.reserve(sum);

    //coarseMat.reserve(sum);
    //Assemble Matrix
    for(size_t np=0; np<mat_coarse.size();++np)
    {
        gsVector<index_t> invMap = patchLocalMappers[np].inverseAsVector();
        for(index_t r=0; r<mats[np].rows();++r)
            for(index_t c=0; c<mats[np].cols();++c)
                entries.add(globalMapper.index(invMap[r],np),globalMapper.index(invMap[c],np),mats[np](r,c));
        //coarseMat.coeffRef(globalMapper.index(invMap[r],np),globalMapper.index(invMap[c],np)) += mats[np](r,c);
    }
    coarseMat.setFrom(entries);

    return gsCoarseSolverAdapter<T>::make(makeSparseCholeskySolver(coarseMat),handler_coarse);
}

template<typename T>
typename gsCoarseSolverAdapter<T>::Ptr constructCoarseSolver(const typename gsMatrixOp<gsSparseMatrix<T> >::Ptr & mat_coarse,
                                                             const gsParallelGlobalLocalHandler::Ptr& handler_coarse)
{
    gsMpiComm comm = handler_coarse->getComm();
    //serialize the data
    //  std::vector<std::vector<T> > buffer(comm.size());
    //   std::vector<std::vector<index_t> >sizes(comm.size());
    //   sizes[comm.rank()].resize(myPatches.size());
    //   for(size_t k=0; k<myPatches.size();++k)
    //       sizes[comm.rank()][k] = mat_coarse[np].rows()*mat_coarse[np].cols();

    std::vector<index_t> sizes(comm.size());
    sizes[comm.rank()] = mat_coarse->matrix().nonZeros();

    comm.sum(sizes.data(),sizes.size());
    size_t sum=0;
    for(size_t i=0; i<sizes.size();++i)
        sum+=sizes[i];

    std::vector<index_t> displ(comm.size(),0);
    for(index_t p=1; p<comm.size();++p)
    {
        displ[p] = displ[p-1]+sizes[p-1];
    }

    std::vector<typename gsConnectionHandler<T>::Triplet> buffer, recvBuffer;
    buffer.reserve(sizes[comm.rank()]);;
    recvBuffer.resize(sum);
    for(index_t p=0; p< comm.size();++p)
    {
        if(p==comm.rank())
        {
            for (int j=0; j<mat_coarse->matrix().outerSize(); ++j)
              for (typename gsSparseMatrix<T>::InnerIterator it(mat_coarse->matrix(),j); it; ++it)
              {
                  typename gsConnectionHandler<T>::Triplet t;
                  t.index[0] = it.row();
                  t.index[1] = it.col();
                  t.value = it.value();
                  buffer.push_back(give(t));
              }
        }
    }
    MPI_Aint offsets[2];
    offsets[0] = offsetof(typename gsConnectionHandler<T>::Triplet, index) ; //most likely going to be 0
    offsets[1] = offsetof(typename gsConnectionHandler<T>::Triplet, value) ;
    MPI_Datatype types[2];
    types[0] = MPITraits<index_t>::getType();
    types[1] = MPITraits<T>::getType();
    int lengths[2];
    lengths[0] = 2;
    lengths[1] = 1;


    MPI_Datatype tempType, tripletType;
    MPI_Aint lb, extent;
    MPI_Type_create_struct(2, lengths, offsets, types,&tempType);

    MPI_Type_get_extent( tempType, &lb, &extent );
    MPI_Type_create_resized( tempType, lb, extent, &tripletType );
    MPI_Type_commit(&tripletType);

    MPI_Allgatherv(buffer.data(),buffer.size(),tripletType,recvBuffer.data(),sizes.data(),displ.data(),tripletType,comm);

    gsSparseMatrix<T> coarseMat(handler_coarse->globalSize(),handler_coarse->globalSize());
    gsSparseEntries<T> entries;
    entries.reserve(sum);

    gsMatrix<index_t> inp2(handler_coarse->globalSize(),1);
    inp2.col(0) = gsVector<index_t>::LinSpaced(handler_coarse->globalSize(),0,handler_coarse->globalSize()-1);
    gsMatrix<index_t>locToGlob;
    locToGlob.setZero(handler_coarse->localSize(),1);
    handler_coarse->extractLocalVector(inp2,locToGlob);
    gsMatrix<index_t> locToGlobTotal;
    locToGlobTotal.setZero(sum,1);


    //coarseMat.reserve(sum);
    //Assemble Matrix


    comm.allgatherv(locToGlob.data(),locToGlob.rows(),locToGlobTotal.data(),sizes.data(),displ.data());


    for(index_t p=0; p<comm.size();++p)
    {
        for(index_t i=0; i<sizes[p];++i)
                entries.add(locToGlobTotal(displ[p]+recvBuffer[displ[p]+i].index[0],0),
                        locToGlobTotal(displ[p]+recvBuffer[displ[p]+i].index[1],0),
                        recvBuffer[displ[p]+i].value);
    }
    coarseMat.setFrom(entries);

    MPI_Type_free(&tripletType);
    MPI_Type_free(&tempType);
    return gsCoarseSolverAdapter<T>::make(makeSparseCholeskySolver(coarseMat),handler_coarse);
}




} // namespace gismo
