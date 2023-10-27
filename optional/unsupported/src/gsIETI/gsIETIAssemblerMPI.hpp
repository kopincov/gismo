/**  gsIETIAssemblerMPI.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on:  2014-12-03
*/


#include <gsCore/gsConfig.h>


#ifdef GISMO_WITH_MPI
#include <gsIETI/gsIETIAssemblerMPI.h>
#include <gsCore/gsField.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsAssembler/gsGaussRule.h>

namespace gismo {


template<class T>
void gsIETI_MPI_Base<T>::init_MPI(const gsIETIInfo& info, const std::vector<std::pair<patchDof,patchDof> >& lagrangeTable, std::vector<index_t>& freeElimLagr)
{
    int rank = m_comm.rank();
    m_patch2proc.resize(info.numberPatches,0);
    m_proc2patch.resize(m_comm.size());


    if(m_patchIdx.size()==0)
    {
    // Assign based on the # of dofs and try to have similar numbered patches
    // together (assuming them to be close)
    const index_t nPatches   = info.numberPatches;
    index_t totalDofNr = 0;
    for (index_t np=0; np<nPatches; ++np)
        totalDofNr += info.dofsB[np]+info.dofsI[np];
    const index_t nProc      = m_comm.size();
    index_t dofs             = 0;
    index_t asignee          = 0;
    for (index_t np=0; np<nPatches; ++np)
    {
        if (dofs*nProc > (asignee+1)*totalDofNr)
            asignee++;

        dofs += info.dofsB[np]+info.dofsI[np];
        // assign to processor "asignee"
        if (asignee == rank)
        {
            m_patchIdx.push_sorted_unique(np);
            m_patch2proc[np] = rank;
        }
    }
    }
    else
    {
        for (size_t np=0; np<m_patchIdx.size(); ++np)
        {
            m_patch2proc[m_patchIdx[np]] = rank;
        }
    }
/*
    //This is a very simple distribution of patches.
    unsigned ratio = cast<T,unsigned>(info.numberPatches/m_comm.size());
    for(size_t np=0; np<ratio*m_comm.size();np++)
        if(np /ratio == rank)
        {
            m_patchIdx.push_sorted(np);
            m_patch2proc[np] = rank;
        }
    for(size_t np=ratio*m_comm.size(); np<info.numberPatches;np++)
        if(np % m_comm.size() == rank)
        {
            m_patchIdx.push_sorted(np);
            m_patch2proc[np] = rank;
        }
        */
    m_comm.sum(m_patch2proc.data(),m_patch2proc.size());


    for(size_t np=0; np<info.numberPatches;np++)
        m_proc2patch[m_patch2proc[np]].push_sorted(np);

    if(rank==0)
    {
        gsInfo<<"Patches to proc (IETI): ";
    for(size_t i=0; i< m_patch2proc.size();i++)
        gsInfo<<m_patch2proc[i]<<" "<<std::flush;
      gsInfo<<"\n";
    }


    //Do the inverse
    {
        size_t iter=0;
        m_req2patch.reserve(info.numberPatches-m_patchIdx.size());
        m_patch2req.resize(info.numberPatches,info.numberPatches);
        for(size_t np=0; np<info.numberPatches;np++)
            if(m_patch2proc[np] != m_comm.rank())
            {
                m_req2patch.push_back(np);
                m_patch2req[np]=iter;
                iter++;
            }
    }
    //Fill the reduced Lagrange table (only containing adjacent lagrage Dofs)
    int iter=0;
    std::vector<index_t> reducedFreeElim;
    m_reducedLagrangeTable.reserve(lagrangeTable.size());
    m_lagLocToGlob.reserve(lagrangeTable.size());
    m_originalLagrangeSize = lagrangeTable.size();
    for(size_t i=0; i<lagrangeTable.size();++i)
        if(m_patch2proc[lagrangeTable[i].first.first]==rank || m_patch2proc[lagrangeTable[i].second.first]==rank)
        {
            m_reducedLagrangeTable.push_back(lagrangeTable[i]);
            if(std::find(freeElimLagr.begin(),freeElimLagr.end(),i)!=freeElimLagr.end())
                reducedFreeElim.push_back(iter);
            iter++;
            m_lagLocToGlob.push_back(i);
        }



    freeElimLagr=reducedFreeElim;

    infoMPI.procNeigbour.clear();
    m_mult.resize(m_reducedLagrangeTable.size(),1);
    for(size_t i =0; i<m_reducedLagrangeTable.size();++i)
    {
        const patchDof& a =  m_reducedLagrangeTable[i].first;
        const patchDof& b = m_reducedLagrangeTable[i].second;
        //  gsInfo<<"rank: "<<rank<<" ( "<<a.first<<" , "<<a.second<<" ) <--> ( "<<b.first<<" , "<<b.second<<" )\n";
        if(m_patch2proc[a.first] != rank)
        {
            infoMPI.procNeigbour.push_sorted_unique(m_patch2proc[a.first]);
            infoMPI.patchNeigbour.push_sorted_unique(a.first);
        }
        else if(m_patch2proc[b.first] != rank)
        {
            infoMPI.procNeigbour.push_sorted_unique(m_patch2proc[b.first]);
            infoMPI.patchNeigbour.push_sorted_unique(b.first);
        }

        if((m_patch2proc[a.first] !=rank) ^ (m_patch2proc[b.first] !=rank))
            m_mult[i]++;
    }

    infoMPI.nNeigbours = infoMPI.procNeigbour.size();
    infoMPI.lagrangeMultReduce = m_reducedLagrangeTable.size();


    /*
    MPI_Comm_group(m_comm, &m_world_group);

    std::vector<int> neig(infoMPI.nNeigbours);
    for(size_t i =0; i<infoMPI.procNeigbour.size(); ++i)
    {
        neig[i] = (int)infoMPI.procNeigbour[i];
    }
    neig.push_back(m_comm.rank());

    MPI_Group_incl(m_world_group, neig.size(),neig.data() , &m_neighbourGoup);
    MPI_Comm_create(m_comm,m_neighbourGoup,&m_neigbours);

    int r=0;
    MPI_Comm_rank(m_neigbours,&r);
    gsInfo<<"My rank in Group: "<<r<<"\n";

    infoMPI.procNeigbourHood.resize(infoMPI.nNeigbours);
    neig.pop_back();
    MPI_Group_translate_ranks(m_world_group,infoMPI.nNeigbours,neig.data(),m_neighbourGoup,infoMPI.procNeigbourHood.data());


    for(int i=0; i<infoMPI.procNeigbourHood.size();++i)
        gsInfo<<infoMPI.procNeigbourHood[i]<<" , ";
*/
    std::vector<std::vector<int> > blocklengths(m_comm.size());
    std::vector<std::vector<int> > displacements(m_comm.size());
    int proc=0;

    //  gsInfo<<"LagSize - rank: "<<m_comm.rank()<<" "<<m_reducedLagrangeTable.size()<<" \n"<<std::flush;
    for(size_t i =0; i<m_reducedLagrangeTable.size();++i)
    {
        const patchDof& a =  m_reducedLagrangeTable[i].first;
        const patchDof& b = m_reducedLagrangeTable[i].second;

        const patchDof& other = m_patchIdx.bContains(a.first) ? b : a;
        //init!
        if(i==0)
        {
            proc = m_patch2proc[other.first];
            for(gsSortedVector<size_t>::iterator pIt=infoMPI.procNeigbour.begin();pIt!=infoMPI.procNeigbour.end();++pIt)
            {
                int p = (*pIt);
                if(p==proc)
                {
                    blocklengths[p].push_back(0);
                    displacements[p].push_back(0);
                }
                else
                    displacements[p].push_back(0);
            }

        }


        if(proc == m_patch2proc[other.first])
        {
            for(gsSortedVector<size_t>::iterator pIt=infoMPI.procNeigbour.begin();pIt!=infoMPI.procNeigbour.end();++pIt)
            {
                if((*pIt)==m_patch2proc[other.first])
                    blocklengths[(*pIt)].back()++;
                else
                    displacements[(*pIt)].back()++;
            }
        }
        else
        {
            for(gsSortedVector<size_t>::iterator pIt=infoMPI.procNeigbour.begin();pIt!=infoMPI.procNeigbour.end();++pIt)
            {
                if((*pIt)==m_patch2proc[other.first]) //new proc
                    blocklengths[(*pIt)].push_back(1);
                else if((*pIt)==proc) //old proc
                    displacements[(*pIt)].push_back(i+1);
                else
                    displacements[(*pIt)].back()++;
            }
        }
        proc = m_patch2proc[other.first];
    }
    for(gsSortedVector<size_t>::iterator pIt=infoMPI.procNeigbour.begin();pIt!=infoMPI.procNeigbour.end();++pIt)
        if(proc!=*pIt)
            displacements[*pIt].pop_back();

    iter=0;
    for(gsSortedVector<size_t>::iterator pIt=infoMPI.procNeigbour.begin();pIt!=infoMPI.procNeigbour.end();++pIt)
    {
        int p = (*pIt);
        //   gsInfo<<"Blocklengths from "<<m_comm.rank()<<" with size "<< blocklengths[p].size()<<" for "<<p<<": ";
        //    for(size_t i=0; i< blocklengths[p].size();i++)
        //       std::cout<<blocklengths[p][i]<<" "<<std::flush;
        //   gsInfo<<"\n";
        //    gsInfo<<"Displacements from "<<m_comm.rank()<<" with size "<< displacements[p].size()<<" for "<<p<<": ";
        //     for(size_t i=0; i< displacements[p].size();i++)
        //         std::cout<<displacements[p][i]<<" "<<std::flush;
        //   gsInfo<<"\n";
        GISMO_ASSERT(blocklengths[p].size() == displacements[p].size(), "number of blocks do not match");

        m_lagDataType.resize(infoMPI.procNeigbour.size());
        MPI_Type_indexed(blocklengths[p].size(),blocklengths[p].data(),displacements[p].data(),MPITraits<T>::getType(),&(m_lagDataType[iter]));
        MPI_Type_commit(&(m_lagDataType[iter]));
        iter++;
    }

    m_isSppHolder=false;
    if(m_options.nSppHolder==1)
    {
        if(rank==0)
            m_isSppHolder=true;

        m_SppMaster.push_sorted_unique(0);
    }
    else if(m_options.nSppHolder==m_comm.size())
    {
        m_isSppHolder=true;
        for(int i=0; i<m_comm.size();++i)
            m_SppMaster.push_sorted(i);
    }
    else
    {
        if(rank<m_options.nSppHolder)
            m_isSppHolder=true;
        for(int i=0; i<m_options.nSppHolder;++i)
            m_SppMaster.push_sorted(i);
    }

    if(m_options.nSppHolder!=m_comm.size() && m_options.nSppHolder!=1)
    {
        MPI_Comm_split(m_comm,rank%m_options.nSppHolder,rank,&m_embComm);
        if(m_isSppHolder)
            MPI_Comm_split(m_comm,0,rank,&m_embMasterComm);
        else
            MPI_Comm_split(m_comm,MPI_UNDEFINED,rank,&m_embMasterComm);
    }
    else
    {
        m_embComm = m_comm;
        m_embMasterComm = m_comm;
    }
    for(int i=0; i<m_comm.size();++i)
        if(i%m_options.nSppHolder == rank%m_options.nSppHolder)
            m_embGroup.push_sorted(i);

    m_commEmbT = gsMpiComm(m_embComm);
    m_commEmbTMaster = gsMpiComm(m_embMasterComm);

    m_idxHolder = -1;
    m_oldPos=0;
    m_nPDofs_group.setZero(m_commEmbT.size());
    m_pDispl_group.setZero(m_commEmbT.size());
    m_nSppDofs.setZero(m_comm.size());
    m_sppDispl.setZero(m_comm.size());
    for(size_t i=0; i<m_embGroup.size();++i)
    {
        int p=m_embGroup[i];
        for(size_t j = 0; j<m_proc2patch[p].size();++j)
            m_nPDofs_group[i]+=info.dofsP[m_proc2patch[p][j]];
    }
    for(int p =1; p<m_commEmbT.size();++p)
        m_pDispl_group[p]+=m_pDispl_group[p-1]+m_nPDofs_group[p-1];

    for(size_t np=0; np<info.numberPatches;np++)
        m_nSppDofs[m_patch2proc[np]]+=info.dofsP[np]*info.dofsP[np];
    for(int p =1; p<m_comm.size();++p)
        m_sppDispl[p]= m_sppDispl[p-1]+m_nSppDofs[p-1];

    if(m_isSppHolder)
        m_idxHolder = *m_SppMaster.find_ptr_or_fail(rank);


    m_lagrangeRequest = new MPI_Request[infoMPI.nNeigbours];
    m_Ibcast_request = new MPI_Request;

    if(m_options.nSppHolder==m_comm.size())
        m_Spp_request = new MPI_Request[1];
    else
        m_Spp_request = new MPI_Request[m_options.nSppHolder];

    m_embeddingT_request = new MPI_Request;

    m_embT_buffer.setZero(info.dofTotalP, info.nRhs);
    m_lagrangeBuff.resize(infoMPI.nNeigbours, gsMatrix<T>(infoMPI.lagrangeMultReduce, info.nRhs));
    m_emb_buffer.setZero(info.dofTotalP,info.nRhs);

}

template<class T>
void gsIETI_MPI_Base<T>::init_embT(const gsIETIInfo& info) const
{
    m_finishedEmbT = false;
    m_embT_buffer.setZero();
    /*
    if(m_isSppHolder)
    {
      //  gsInfo<<"rank "<<m_comm.rank()<<" is waiting for embT\n";

        for(size_t np =0; np<info.numberPatches;++np)
        {
            if(m_patch2proc[np]!=m_comm.rank())
            {
                m_embT_buffer[np].setZero(info.dofsP[np],info.nRhs);
                MPI_Irecv(m_embT_buffer[np].data(),info.dofsP[np]*info.nRhs,MPITraits<T>::getType(),m_patch2proc[np],0,m_comm,&m_embeddingT_request[m_patch2req[np]]);

            }
        }


    }
    */

}

template<class T>
void gsIETI_MPI_Base<T>::send_embT(const gsIETIInfo& info, size_t np, size_t nRhs) const
{
    if(np==m_patchIdx.back())
        send_embT(info, nRhs);
}

template<class T>
void gsIETI_MPI_Base<T>::send_embT(const gsIETIInfo& info, size_t nRhs) const
{
    int n = m_embT_buffer.rows()*m_embT_buffer.cols();
    if(m_options.nSppHolder==1)
    {
        m_comm.isum(m_embT_buffer.data(),n,0,m_embeddingT_request);
        //  m_comm.sum(m_embT_buffer.data(),n,0);
    }
    else if(m_options.nSppHolder==m_comm.size())
    {
        m_comm.isum(m_embT_buffer.data(),n,m_embeddingT_request);
        //    m_comm.sum(m_embT_buffer.data(),n);
    }
    else
    {
        m_commEmbT.sum(m_embT_buffer.data(),n);
        if(m_isSppHolder)
            m_commEmbTMaster.isum(m_embT_buffer.data(),n,m_embeddingT_request);
    }
}

template<class T>
void gsIETI_MPI_Base<T>::finalize_embT(gsMatrix<T>& result) const
{
    /*
    int idx;
    for(size_t i=0; i<m_patchIdx.size();++i)
    {
        size_t np=m_patchIdx[i];
        assPrim_fun(m_embT_buffer[np],np,result);
        //Base::assemblePrimal(m_embT_buffer[np],np,result);
    }
 //   gsInfo<<"first part2 Done!\n";
    for(size_t i=0; i<m_req2patch.size();++i)
    {
        MPI_Waitany(m_req2patch.size(),m_embeddingT_request,&idx,MPI_STATUS_IGNORE);
        size_t np=m_req2patch[idx];
        assPrim_fun(m_embT_buffer[np],np,result);
        //Base::assemblePrimal(m_embT_buffer[np],np,result);
    }
   //  gsInfo<<"second part2 Done!\n";
   */

    //so dont call MPI_WAIT if there are several SppHolder (but not all) and you are not a holder.
    if(m_comm.size()==m_options.nSppHolder || isSppHolder()) //TODO: check if m_comm.size()==1... why not all? see send_embT
        MPI_Wait(m_embeddingT_request,MPI_STATUS_IGNORE);

    if(isSppHolder())
        result = m_embT_buffer;
    m_finishedEmbT=true;
}

template<class T>
void gsIETI_MPI_Base<T>::init_Spp(const gsIETIInfo& info)
{
    /*
    if(m_comm.rank()==0)
    {
        for(size_t np =0; np<info.numberPatches;++np)
        {
            if(m_patch2proc[np]!=0)
            {
                m_Spp_buffer[np].setZero(info.dofsP[np],info.dofsP[np]);
                MPI_Irecv(m_Spp_buffer[np].data(),info.dofsP[np]*info.dofsP[np],MPITraits<T>::getType(),m_patch2proc[np],1,m_comm,&m_Spp_request[m_patch2req[np]]);
            }
        }
    }
    */
    size_t nP,tnP;
    nP=tnP =0;
    for(size_t i=0; i<m_patchIdx.size();i++)
        nP+=info.dofsP[m_patchIdx[i]]*info.dofsP[m_patchIdx[i]];

    for(size_t np =0; np<info.numberPatches;np++)
        tnP+=info.dofsP[np]*info.dofsP[np];
    m_Spp_sendBuffer.setZero(nP,1);
    if(m_isSppHolder)
        m_Spp_recvBuffer.setZero(tnP,1);
}



template<class T>
void gsIETI_MPI_Base<T>::sendSPP()
{
    if(m_options.nSppHolder==1)
    {
        /*
        if(m_comm.rank()!=0)
        {
            MPI_Request req;
            MPI_Isend(m_Spp_sendBuffer.data(),m_Spp_sendBuffer.rows()* m_Spp_sendBuffer.cols(),MPITraits<T>::getType(),0,1,m_comm,&req);
            MPI_Request_free(&req);
        }
        */
        MPI_Igatherv(m_Spp_sendBuffer.data(),m_Spp_sendBuffer.rows()*m_Spp_sendBuffer.cols(),MPITraits<T>::getType(),
                     m_Spp_recvBuffer.data(), m_nSppDofs.data(), m_sppDispl.data(), MPITraits<T>::getType(),0,
                     m_comm, &m_Spp_request[0]);
    }
    else if(m_options.nSppHolder==m_comm.size())
    {
        MPI_Iallgatherv(m_Spp_sendBuffer.data(), m_Spp_sendBuffer.rows(), MPITraits<T>::getType(),
                        m_Spp_recvBuffer.data(), m_nSppDofs.data(), m_sppDispl.data(), MPITraits<T>::getType(),
                        m_comm, &m_Spp_request[0]);
    }
    else
    {
        for(size_t i=0; i<m_SppMaster.size();++i)
            MPI_Igatherv(m_Spp_sendBuffer.data(),m_Spp_sendBuffer.rows()*m_Spp_sendBuffer.cols(),MPITraits<T>::getType(),
                         m_Spp_recvBuffer.data(), m_nSppDofs.data(), m_sppDispl.data(), MPITraits<T>::getType(),m_SppMaster[i],
                         m_comm, &m_Spp_request[i]);

    }
}

template<class T>
void gsIETI_MPI_Base<T>::finalizeSpp(const gsIETIInfo& info, const std::vector<std::vector<index_t> >& pDofsLoc2Glob,typename gsIETIAssembler<T>::sparseSPDfact*& LU_Spp)
{
    if(info.dofTotalP==0)
        return;

    // int idx=0;

    gsSparseMatrix<T> Spp(info.dofTotalP,info.dofTotalP);

    unsigned max = *max_element(info.dofsP.begin(),info.dofsP.end());
    //unsigned min = *min_element(info.dofsP.begin(),info.dofsP.end());
    // if(m_comm.rank()==0) gsInfo<<"pDofs: min: "<<min<<", max: "<<max<<"\n";

    Spp.reserve(gsVector<int>::Constant(info.dofTotalP, math::pow(3,info.dim)*max)); //TODO: improve

    // if(m_req2patch.size()!=0)
    //     MPI_Waitall(m_req2patch.size(),m_Spp_request,MPI_STATUSES_IGNORE);
    // MPI_Wait(&m_Spp_request[0],MPI_STATUSES_IGNORE);
    size_t np=0;
    size_t idx = 0;
    gsMatrix<T>  copy;

    for(size_t i=0; i<m_patchIdx.size();i++)
    {
        np=m_patchIdx[i];
        //  m_Spp_sendBuffer.resize(info.dofsP[np]*info.dofsP[np]);
        //  copy = m_Spp_recvBuffer.block(m_sppDispl[m_comm.rank()]+idx,0,info.dofsP[np]*info.dofsP[np],1);
        typename gsMatrix<T>::Block block = m_Spp_sendBuffer.block(idx,0,info.dofsP[np]*info.dofsP[np],1);
        // copy.resize(info.dofsP[np],info.dofsP[np]);

        for(index_t c=0; c<info.dofsP[np];++c)
            for(index_t r=0; r<info.dofsP[np];++r)
                //Spp.coeffRef(pDofsLoc2Glob[np][r],pDofsLoc2Glob[np][c])+= m_Spp_buffer[np](r,c);
                Spp.coeffRef(pDofsLoc2Glob[np][r],pDofsLoc2Glob[np][c])+= block(c*info.dofsP[np]+r,0);
        idx+=info.dofsP[np]*info.dofsP[np];
    }
    //  gsInfo<<"first part Done!\n";
    if(m_SppMaster.size()==(size_t)m_comm.size())
    {
        MPI_Wait(&m_Spp_request[0],MPI_STATUS_IGNORE);
        for(int p=0; p<m_comm.size();++p)
        {
            if(p==m_comm.rank())
                continue;
            idx=0;
            for(size_t i =0; i<m_proc2patch[p].size();i++)
            {
                np=m_proc2patch[p][i];
                typename gsMatrix<T>::Block block = m_Spp_recvBuffer.block(m_sppDispl[p]+idx,0,info.dofsP[np]*info.dofsP[np],1);
                // copy.resize(info.dofsP[np],info.dofsP[np]);
                //  MPI_Waitany(m_req2patch.size(), m_Spp_request,&idx,MPI_STATUS_IGNORE);

                for(index_t c=0; c<info.dofsP[np];++c)
                    for(index_t r=0; r<info.dofsP[np];++r)
                        //  Spp.coeffRef(pDofsLoc2Glob[np][r],pDofsLoc2Glob[np][c])+= m_Spp_buffer[np](r,c);
                        Spp.coeffRef(pDofsLoc2Glob[np][r],pDofsLoc2Glob[np][c])+= block(c*info.dofsP[np]+r,0);

                idx+=info.dofsP[np]*info.dofsP[np];
            }
        }
    }
    else
    {
        int proc;
        MPI_Wait(&m_Spp_request[m_idxHolder],MPI_STATUS_IGNORE);
        for(int p=0; p<m_comm.size();++p)
        {
            if(p==m_comm.rank())
                continue;
            idx=0;
            proc=p;
            //   MPI_Waitany(m_comm.size()-1, m_Spp_request,&proc,MPI_STATUS_IGNORE);
            for(size_t i =0; i<m_proc2patch[proc].size();i++)
            {
                np=m_proc2patch[proc][i];
                typename gsMatrix<T>::Block block = m_Spp_recvBuffer.block(m_sppDispl[proc]+idx,0,info.dofsP[np]*info.dofsP[np],1);
                // copy.resize(info.dofsP[np],info.dofsP[np]);
                //  MPI_Waitany(m_req2patch.size(), m_Spp_request,&idx,MPI_STATUS_IGNORE);

                for(index_t c=0; c<info.dofsP[np];++c)
                    for(index_t r=0; r<info.dofsP[np];++r)
                        //  Spp.coeffRef(pDofsLoc2Glob[np][r],pDofsLoc2Glob[np][c])+= m_Spp_buffer[np](r,c);
                        Spp.coeffRef(pDofsLoc2Glob[np][r],pDofsLoc2Glob[np][c])+= block(c*info.dofsP[np]+r,0);

                idx+=info.dofsP[np]*info.dofsP[np];
            }
        }
    }

    /*
    for(size_t i=0; i<m_req2patch.size();i++)
    {
      //  MPI_Waitany(m_req2patch.size(), m_Spp_request,&idx,MPI_STATUS_IGNORE);
        np=m_req2patch[idx];
        //   gsInfo<<"Waiting, got "<<idx<<"\n";
        for(index_t c=0; c<info.dofsP[np];++c)
            for(index_t r=0; r<info.dofsP[np];++r)
              //  Spp.coeffRef(pDofsLoc2Glob[np][r],pDofsLoc2Glob[np][c])+= m_Spp_buffer[np](r,c);
                Spp.coeffRef(pDofsLoc2Glob[np][r],pDofsLoc2Glob[np][c])+= m_Spp_recvBuffer.block(idx,0,info.dofsP[np]*info.dofsP[np],1).resize(info.dofsP[np]*info.dofsP[np])(r,c);
      idx+=info.dofsP[np]*info.dofsP[np];
    }
    */
    //  gsInfo<<"second part Done!\n";
    // gsInfo<<"Spp: \n"<<Spp.toDense()<<std::endl<<std::endl;
    Spp.makeCompressed();

    // Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real_t,Dynamic,Dynamic>  > eig(Spp.toDense());
    // gsInfo<<"Eigenvalues of Spp:\n"<<eig.eigenvalues()<<std::endl<<std::endl;

    //  gsInfo<<"Do third part!\n";
    LU_Spp=new typename gsIETIAssembler<T>::sparseSPDfact(Spp);
}

template<class T>
void gsIETI_MPI_Base<T>::initEmbedding(const gsIETIInfo& info, gsMatrix<T> & xP)const
{
    if(!m_isSppHolder)
        xP.setZero(info.dofTotalP,info.nRhs);

    if(m_comm.size()!=1 && m_options.nSppHolder!=m_comm.size())
    {
        // gsInfo<<"rank: "<<m_comm.rank()<<" Doing BCAST\n";
        MPI_Ibcast(xP.data(),xP.cols()*xP.rows(),MPITraits<T>::getType(),0,m_embComm, m_Ibcast_request);
    }
}

template<class T>
void gsIETI_MPI_Base<T>::initEmbedding_(const gsIETIInfo& info,gsMatrix<T>& xP, std::vector<gsMatrix<T> > & xPis)const
{
    //    int dummy;
    //GISMO_ASSERT(info.nRhs==1,"This method does not work for nRhs!=1");

    //gsInfo<<"Rank: "<<m_comm.rank()<<":  -- "<<xP.rows()<< "  "<<xP.cols()<< " -- " << std::flush;

    int offset = 0;
    if(!m_isSppHolder)
        xP.setZero(m_nPDofs_group[m_commEmbT.rank()], xP.cols() == 0 ? info.nRhs:xP.cols());


    for(int i = 0; i < xP.cols(); i++)
    {
        if (m_isSppHolder)
            MPI_Scatterv(xP.data() + offset, const_cast<int*>(m_nPDofs_group.begin()),const_cast<int*>(m_pDispl_group.begin()), MPITraits<T>::getType(),
                         MPI_IN_PLACE, m_nPDofs_group[m_commEmbT.rank()] , MPITraits<T>::getType(), 0, m_embComm);
        else
        {
            MPI_Scatterv(0, const_cast<int*>(m_nPDofs_group.begin()),const_cast<int*>(m_pDispl_group.begin()), MPITraits<T>::getType(),
                         xP.data() + offset, m_nPDofs_group[m_commEmbT.rank()], MPITraits<T>::getType(), 0, m_embComm);
        }
        offset += xP.rows();
    }

    int iter=0;
    //   gsInfo<<"rank: "<<m_comm.rank()<<"  xP_reorder is given by "<<xP.transpose()<<"\n";
    for(size_t i=0; i<m_patchIdx.size();++i)
    {
         xPis[i]= xP.block(iter,0,info.dofsP[m_patchIdx[i]],xP.cols());
         iter+=info.dofsP[m_patchIdx[i]];
    }
}


//**********************************************************************************
//**********************************************************************************
//**********************************************************************************

template<class T>
gsIETIAssemblerMPI<T>::gsIETIAssemblerMPI(gsAssembler<T> &assembler,gsSortedVector<size_t> myPatches, MPI_Comm comm)
    :Base(assembler), MPI_Base(myPatches,comm)
{
    m_IETIoptions.opt.update(defaultOptions(),gsOptionList::addIfUnknown);
}

template<class T>
gsOptionList gsIETIAssemblerMPI<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();

    opt.addInt   ("NumberSppHolder","number of Spp holder",1);

    return opt;
}

template<class T>
void gsIETIAssemblerMPI<T>::setOptions(const gsOptionList &opt)
{
    gsIETIAssembler<T>::setOptions(opt);
    m_IETIoptions.nSppHolder = m_IETIoptions.opt.getInt("NumberSppHolder");
    gsIETI_MPI_Base<T>::setOptions(m_IETIoptions);
}


template<class T>
void gsIETIAssemblerMPI<T>::init()
{
    Base::init();

    MPI_Base::init_MPI(info, m_lagrangeTable, m_freeElimLagr);
}

template<class T>
void gsIETIAssemblerMPI<T>::assemble(const gsMultiPatch<T> &curSol)
{
    Base::assembleInit();
    bool assembleRHS = !m_IETIoptions.opt.getSwitch("NoAssemblyOfRHS");
    MPI_Base::init_embT(info);
    MPI_Base::init_Spp(info);

    //   std::vector<gsMatrix<T> > spp_loc(info.numberPatches);
    std::vector<gsMatrix<T> > rhs_loc(info.numberPatches);

    std::vector<gsSparseMatrix<T> > matrices(info.numberPatches);
    std::vector<gsSparseMatrix<T> > matrices2(info.numberPatches);

    Base::reserveSpace(m_dirDofs,rhs_loc);

    gsStopwatch time, time2;
    Eigen::setNbThreads(1);
    size_t np;
    if(!m_IETIoptions.opt.getSwitch("ExtraTiming")) //The Fast One
    {
#pragma omp parallel
        {
            // gsSparseMatrix<T> matrix, matrix2;
            gsAssembler<T>* A = m_assembler->create();

#pragma omp for  schedule(static,1) nowait
            //for (std::vector<size_t>::const_iterator iNp=m_patchIdx.begin(); iNp!=m_patchIdx.end(); ++iNp )
            for (size_t i=0; i<m_patchIdx.size();++i)
            {
                np=m_patchIdx[i];
                //   gsInfo<<"I am on Patch "<<np<<"\n"<<std::flush;
                Base::assembleLocal(A, matrices[np], rhs_loc[np],m_dirDofs[np],np,curSol);
                // -- Reordering
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                {
                    matrices2[np] = matrices[np];
                    Base::makeReorderingPrimalRem(matrices2[np],np);
                }
                Base::makeReorderingBoundInt(matrices[np],np);

                // The main matrices for applying the system matrix
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                    Base::assembleKrrKrpKpp(matrices2[np],np);
                else
                {
                    Base::assembleC(np);
                    Base::assembleLUofKC(matrices[np],np);
                }

                //assemble local Spp
                //    Base::assembleSppLoc(MPI_Base::m_Spp_buffer[np], np);
               // m_comm.barrier(); gsInfo<<"rank: "<<m_comm.rank()<<" , doing patch "<<np<<"\n";
                Base::assembleSppLoc(MPI_Base::m_Spp_sendBuffer, np);
            }
            MPI_Base::sendSPP();


            for(size_t i=0; i<m_patchIdx.size();i++)
            {
                np=m_patchIdx[i];

                //For Preconditioning and processing the rhs and solution ICO minimum energy
                Base::assembleKiiKibKbb(matrices[np], np);

                //assemble local rhs
                if(assembleRHS) this->assembleRhs(rhs_loc, np);
                //m_comm.barrier();  gsInfo<<"rank: "<<m_comm.rank()<<" , doing patch "<<np<<"\n";

            } //end-for

            //free memory
            delete A;

        }//end-parallel

    }//end if
    else // The slow one
    {
#pragma omp parallel
        {
            gsAssembler<T>* A = m_assembler->create();

#pragma omp for  schedule(static, 1)
            for (size_t i=0; i<m_patchIdx.size();++i)
                np=m_patchIdx[i], Base::assembleLocal(A, matrices[np], rhs_loc[np],m_dirDofs[np], np,curSol);

            delete A;

            printTime(time,"Time for assemling all patch local matrices: ");

            // -- Reordering
#pragma omp for  schedule(static, 1)
            for (size_t i=0; i<m_patchIdx.size();++i)
            {
                np=m_patchIdx[i];
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                {
                    matrices2 = matrices;
                    Base::makeReorderingPrimalRem(matrices2[np],np);
                }
                Base::makeReorderingBoundInt(matrices[np],np);
            }

            printTime(time,"Time for reordering all patch local matrices: ");

            //------------------------------------------------------------------------------/
            // The main matrices for applying the system matrix
            if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
            {
#pragma omp for schedule(static, 1)
                for (size_t i=0; i<m_patchIdx.size();++i)
                    np=m_patchIdx[i],Base::assembleKrrKrpKpp(matrices2[np],np);

                printTime(time,"Time for computing LU factorization of Krr : ");
            }
            else
            {
#pragma omp for  schedule(static, 1)
                for (size_t i=0; i<m_patchIdx.size();++i)
                    np=m_patchIdx[i],Base::assembleC(np);

                printTime(time,"Time for computing C : ");

#pragma omp for  schedule(static, 1)
                for (size_t i=0; i<m_patchIdx.size();++i)
                    np=m_patchIdx[i],Base::assembleLUofKC(matrices[np],np);

                printTime(time,"Time for computing LU factorization of KC : ");
            }

#pragma omp for  schedule(static, 1)
            for (size_t i=0; i<m_patchIdx.size();++i)
                np=m_patchIdx[i],// Base::assembleSppLoc(MPI_Base::m_Spp_buffer[np], np);
                        Base::assembleSppLoc(MPI_Base::m_Spp_sendBuffer, np);

            MPI_Base::sendSPP();
            printTime(time,"Time for assemling local Spp: ");

            //For Preconditioning and processing the rhs and solution ICO minimum energy
#pragma omp for  schedule(static, 1)
            for (size_t i=0; i<m_patchIdx.size();++i)
                np=m_patchIdx[i],Base::assembleKiiKibKbb(matrices[np], np);

            printTime(time,"Time for calculating LU of Kii : ");


#pragma omp for  schedule(static, 1)
            for (size_t i=0; i<m_patchIdx.size();++i)
               { np=m_patchIdx[i]; if(assembleRHS) this->assembleRhs(rhs_loc, np);}

            printTime(time,"Time for assemling rhs: ");

        }
    }
    printTime(time2,"Time for doing parallel assembling part: ");

    //------------Do Serial stuff--------------
    Eigen::setNbThreads(0);

    //1. Rhs
    if(assembleRHS) assembleRhsFreeElim();
    // m_comm.barrier(); gsInfo<<"rank: "<<m_comm.rank()<<" finished elimRhs"<<"\n";
    time.restart();
    //2. Spp
    if(MPI_Base::isSppHolder()) //m_comm.rank()==0)
    {
        MPI_Base::finalizeSpp(info, m_pDofsLoc2Glob,m_LU_Spp);
        printTime(time,"Time for assemling Spp: ");
        // gsInfo<<"rank: "<<m_comm.rank()<<" finished Spp"<<"\n";
        /*
        MPI_Base::finalizeEmbeddingTrans(
        [=] (const gsMatrix<T>& a, int b, gsMatrix<T>& c) {
            Base::assemblePrimal(a,b,c);
        }
        ,m_rhs_p);
        */
        // finalizeEmbeddingTrans(MPI_Base::m_embT_buffer,m_rhs_p);
        // printTime(time,"Time for assemling final-rhs: ");
    }

    printTime(time2,"Time for doing serial assembling part: ");
}


template<class T>
void gsIETIAssemblerMPI<T>::storeSpp(const gsMatrix<T>& in, gsMatrix<T>& out)
{
    gsMatrix<T> temp = in;
    temp.resize(in.cols()*in.rows(),1);
    // gsInfo<<"rank: "<<m_comm.rank()<<" idx: "<<MPI_Base::m_oldPos<<"\n"<<in<<"\n";

    out.block(MPI_Base::m_oldPos,0,in.cols()*in.rows(),1)=temp;//resize(in.cols()*in.rows(),1);
    MPI_Base::m_oldPos+=  in.cols()*in.rows();
}

template<class T>
void gsIETIAssemblerMPI<T>::assembleRhs(const std::vector<gsMatrix<T> >& rhs_loc, size_t np)
{

    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy") )
    {
        //gsInfo<<"nRhs: "<<info.nRhs<<"\n";
        gsMatrix<T> rhs_b(info.dofsB[np], info.nRhs);
        m_rhs_d[np].setZero(info.dofsB[np],info.nRhs);
        m_rhs_i[np].setZero(info.dofsI[np],info.nRhs);

        unsigned idx;

        for(int k=0;k<rhs_loc[np].rows();k++)
        {
            idx = m_glob2BoundInteriorIndex[np][k];
            if(m_globIsBoundIndex[np][k] ==false)
                m_rhs_i[np].row(idx)=rhs_loc[np].row(k);
            else
                rhs_b.row(idx)=rhs_loc[np].row(k);
        }

        gsMatrix<T> temp;
        Base::template solveKii<true>(np,m_rhs_i[np], temp);
       // temp.transposeInPlace();
       // gsMatrix<T> temp2 = temp*m_Kib[np];
       // temp2.transposeInPlace();

        rhs_b=rhs_b-*m_Kbi[np]*temp; //fb - Kbi*Kii^-1*fi


        embeddingTransPartial(rhs_b,np,m_rhs_d[np]);
        MPI_Base::send_embT(info,np,rhs_b.cols());
    }
    else
    {
        GISMO_ERROR("This is not working any more!");
        //this part could be covered by embeddingTrans
        m_rhs_d[np].setZero(info.dofsR[np],info.nRhs);
        //   MPI_Base::m_embT_buffer[np].setZero(info.dofsP[np],info.nRhs);

        for(size_t k=0;k<m_primalVdofs[np].size();k++)
        {
            //     int c = Base::getComp(np,m_primalVdofs[np][k]);
            //     MPI_Base::m_embT_buffer[np].row(k)= rhs_loc[np].row( m_locDofsMapper[np][c].index(Base::compCalcBack(np,m_primalVdofs[np][k])));
        }

        // MPI_Base::send_embT(info,np,info.nRhs);

        for(size_t i=0;i<m_remDofs[np].size();i++)
        {
            int c = Base::getComp(np,m_remDofs[np][i]);
            m_rhs_d[np].row(i) = rhs_loc[np].row(m_locDofsMapper[np][c].index(Base::compCalcBack(np,m_remDofs[np][i])));
        }

    }
}

template <typename T>
void gsIETIAssemblerMPI<T>::giveAssembledMatrices(std::vector<gsSparseMatrix<T> >& matrices,const gsMatrix<T>& rhs)
{
    Base::m_useGivenMats = true;
    Base::m_tempMatrix.resize(matrices.size());
    Base::m_tempRhs.resize(matrices.size());

    Base::m_tempMatrix.swap(matrices);
    for (size_t i=0; i<m_patchIdx.size();++i)
    {
        size_t np=m_patchIdx[i];
        extractPatch(np,rhs,Base::m_tempRhs[np]);
    }
}


template <typename T>
void gsIETIAssemblerMPI<T>::extractPatch(size_t np, const gsMatrix<T>& rhs, gsMatrix<T>& rhsLocal) const
{
    const gsDofMapper& locMapper = m_locDofsMapper[np][0];
    const gsDofMapper& globMapper = m_assembler->system().colMapper(0);
    rhsLocal.setZero(locMapper.freeSize(),rhs.cols());

    for(size_t i=0; i<(size_t)locMapper.size();++i)
    {
        unsigned glIdxAss = globMapper.index(i,np);
        unsigned glIdx = m_stdMapper[0].index(i,np);
        unsigned idx = locMapper.index(i);
        if(!globMapper.is_free_index(glIdxAss)) continue;

        if(m_globIsBoundIndex[np][idx])
        {
            const std::set<patchDof>& set = (m_globalConnectionTable.find(glIdx))->second;
            int count =0;
            for(typename std::set<patchDof>::const_iterator it = set.begin();it!=set.end();++it)
                if(m_patchIdx.bContains(it->first))
                    count++;

            rhsLocal.row(idx) = rhs.row(glIdxAss)/count;
        }
        else
            rhsLocal.row(idx) = rhs.row(glIdxAss);
    }
}


template<typename T>
void gsIETIAssemblerMPI<T>::setNewRhs(const gsMatrix<T>& rhs)
{
    m_rhs_p.setZero(info.dofTotalP,info.nRhs);
    MPI_Base::init_embT(info);
    std::vector<gsMatrix<T> > rhs_loc(m_patches.nPatches());
  //  gsInfo<<"size of m_patchIdx: "<<m_patchIdx.size()<<"\n";
    for (size_t i=0; i<m_patchIdx.size();++i)
    {
        size_t np=m_patchIdx[i];
        extractPatch(np,rhs,rhs_loc[np]);
     //   gsInfo<<"Rhs on Patch: "<<np<<"\n "<<rhs_loc[np].transpose()<<"\n"<<std::flush;
        assembleRhs(rhs_loc, np);
    }

    /*
    for(size_t npi = 0; npi<info.numberPatches;++npi)
    {
        m_comm.barrier();
        if(gsMpi::worldRank() == m_patch2proc[npi])
            //gsInfo<<"patch: "<<np<<"  -- rhs_i: \t"<<m_rhs_i[np].transpose()<<"\n \n temp: \t"<<temp.transpose()<<"\n\n xI:\t"<<xI.transpose()<<"\n\n w:\t"<<w.transpose()<<"\n\n solVec:\t"<<solVec.transpose()<<"\n\n";
            gsInfo<<"rhs_i "<<npi<<":\t"<<rhs_loc[npi].transpose()<<"\n";
    }
    */
    assembleRhsFreeElim();
}

template<class T>
void gsIETIAssemblerMPI<T>::assembleRhsFreeElim()
{
    //assemble the rhs for the coupled eliminated and free dofs (as lagrange mult.)
    if(!Base::nCoupledElimDofs())
    {
        int sign;
        // m_rhs_dir.setZero(m_lagrangeTable.size(),info.nRhs);
        m_rhs_dir.setZero(MPI_Base::m_reducedLagrangeTable.size(),info.nRhs);
        for(std::vector<index_t>::iterator it = m_freeElimLagr.begin(); it!=m_freeElimLagr.end();it++)
        {
            patchDof p1= MPI_Base::m_reducedLagrangeTable[*it].first;
            patchDof p2 =MPI_Base:: m_reducedLagrangeTable[*it].second;
            if(m_patch2proc[p1.first]!= m_comm.rank() && m_patch2proc[p2.first]!= m_comm.rank())
                continue;

            p1.first > p2.first ? sign = 1 : sign = -1;

            int c = Base::getComp(p2.first,p2.second); // p1 and p2 must have the same component
            if(m_locDofsMapper[p1.first][c].is_free(Base::compCalcBack(p1.first,p1.second)) && m_patch2proc[p2.first]== m_comm.rank())
                m_rhs_dir.row(*it)=sign*m_dirDofs[p2.first].row(m_locDofsMapper[p2.first][c].bindex(Base::compCalcBack(p2.first,p2.second)));
            else if(m_locDofsMapper[p2.first][c].is_free(Base::compCalcBack(p2.first,p2.second))&& m_patch2proc[p1.first]== m_comm.rank())
                m_rhs_dir.row(*it)=-sign*m_dirDofs[p1.first].row(m_locDofsMapper[p1.first][c].bindex(Base::compCalcBack(p1.first,p1.second)));
            else
            {
                continue;
                //this error still holds, however, for MPI we have an additional case which can be handled easier if we do not care about
                //the case when both are ddofs. (we hope that this does not happen)
                //GISMO_ERROR("BOTH DOFS ARE DIRICHLET DOFS; THIS IS NOT ALLOWED")
            }
        }
    }
}


template<class T>
void gsIETIAssemblerMPI<T>::processSolution(const gsMatrix<T>& uP_, const std::vector<gsMatrix<T> >& u2, gsMatrix<T>& solVec) const
{


    unsigned totSize = 0;
    for(size_t c = 0;c<info.cDim;c++)
        totSize+= m_stdMapper[c].freeSize();

    solVec.setZero(totSize,info.nRhs);

    std::vector<gsDofMapper > stdMapper(info.cDim);
    for(size_t c=0; c<info.cDim;c++)
        stdMapper[c]= m_assembler->system().colMapper(c);

    gsMatrix<T> uP = uP_;
    std::vector<gsMatrix<T> > uPis(m_patchIdx.size());
    initEmbedding(info,uP,uPis);

    if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy") )
    {
#pragma omp parallel
        {
            gsMatrix<T> xPi;
#pragma omp for schedule(static, 1) nowait
            for (size_t i=0; i<m_patchIdx.size();++i)
            {
                size_t np=m_patchIdx[i];
                Base::distributePrimal(uP,np,xPi);

                for(size_t c=0;c<info.cDim;c++)
                {
                    int sz  = m_basis[c][np].size();
                    for (index_t i = 0; i < sz; ++i)
                    {
                        int idx = m_locDofsMapper[np][c].index(i);
                        int globI = stdMapper[c].index(i,np);

                        if(stdMapper[c].is_free(i,np))
                        {
                            if(m_globIsPrimalIndex[np][idx])
                                solVec.row(globI) = xPi.row(m_glob2PrimalRemainingIndex[np][idx]);// /m_globalConnectionTable.find(globI)->second.size();
                            else
                                solVec.row(globI) = u2[np].row(m_glob2PrimalRemainingIndex[np][idx]);
                        }
                    }
                }
            }

        }
    }
    else
#pragma omp parallel
    {
        gsMatrix<T> xI, temp, w;

        std::vector<int> count(solVec.rows(),0);
#pragma omp for schedule(static, 1) nowait
        for (size_t i=0; i<m_patchIdx.size();++i)
        {
            size_t np=m_patchIdx[i];

            embedding(uPis[i],u2[np],np,w);
            temp.noalias() = m_rhs_i[np]-*m_Kib[np]*w;

            Base::template solveKii<true>(np,temp, xI) ;



            for(size_t c=0;c<info.cDim;c++)
            {
                int sz  =  m_basis[c][np].size();
                for (index_t i = 0;i < sz ; ++i)
                {
                    int idx = m_locDofsMapper[np][c].index(i);
                    int globI = stdMapper[c].index(i,np);
                    if(stdMapper[c].is_free(i,np))
                    {
                        if(m_globIsBoundIndex[np][idx])
                        {
                            solVec.row(globI) += w.row(m_glob2BoundInteriorIndex[np][idx]);// /(T)m_globalConnectionTable.find(globI)->second.size();
                            count[globI]++;
                        }
                        else
                            solVec.row(globI)  = xI.row(m_glob2BoundInteriorIndex[np][idx]);
                    }
                }
            }

            /*
            for(size_t npi = 0; npi<info.numberPatches;++npi)
            {
                m_comm.barrier();
                if(np == npi)
                    //gsInfo<<"patch: "<<np<<"  -- rhs_i: \t"<<m_rhs_i[np].transpose()<<"\n \n temp: \t"<<temp.transpose()<<"\n\n xI:\t"<<xI.transpose()<<"\n\n w:\t"<<w.transpose()<<"\n\n solVec:\t"<<solVec.transpose()<<"\n\n";
                    gsInfo<<"solVec:\t"<<solVec.transpose()<<"\n";
            }
            */
        }

        for (size_t i=0; i<count.size();++i)
        if(count[i]>=2)
            solVec.row(i)/=T(count[i]);


    }

}

template<class T>
void gsIETIAssemblerMPI<T>::combineToCommonSolution(gsMatrix<T>& solVec) const
{
    bool minEnergy = !m_IETIoptions.opt.getSwitch("NoMinimumEnergy");
    std::vector<gsDofMapper > stdMapper(info.cDim);
    std::vector<bool> checked(solVec.rows(),false);
    for(size_t c=0; c<info.cDim;c++)
        stdMapper[c]= m_assembler->system().colMapper(c);


    std::set<int> procs;
    for (size_t pI=0; pI<m_patchIdx.size();++pI)
    {
        size_t np=m_patchIdx[pI];
        for(size_t c=0;c<info.cDim;c++)
        {
            int sz  = m_basis[c][np].size();
            for (index_t i = 0; i < sz; ++i)
            {
                int idx = m_locDofsMapper[np][c].index(i);
                int globI = m_stdMapper[c].index(i,np);
                int globIAss = stdMapper[c].index(i,np);

                if(stdMapper[c].is_free(i,np) && ( (  minEnergy && m_globIsBoundIndex[np][idx]) || ( !minEnergy &&  m_globIsPrimalIndex[np][idx]))&& !checked[globIAss])
                {
                    const std::set<patchDof>& set = m_globalConnectionTable.find(globI)->second;
                    procs.clear();
                    for(typename std::set<patchDof>::const_iterator pd = set.begin() ; pd!=set.end();  pd++)
                        procs.insert(MPI_Base::m_patch2proc[pd->first]);
                    solVec.row(globIAss) /= cast<size_t, T>(procs.size());
                    checked[globIAss]=true;
                }

            }

        }

    }
    m_comm.sum(solVec.data(), solVec.size());

}

/*
template<class T>
void gsIETIAssemblerMPI<T>::embedding(gsMatrix<T> & xP, gsMatrix<T> const & x2, size_t np, gsMatrix<T>& w) const
{
    int nRhs = x2.cols();
    gsMatrix<T> xPi;

    w.setZero(info.dofsB[np],nRhs);

    Base::distributePrimal(xP,np,xPi);
    if(m_options.isMinimumEnergy)
    {
        w=x2;

        if(info.dofsP[np]!=0)
            w.noalias()+=m_Phi[np]*xPi;
    }
    else
    {
        unsigned globI;

        for (unsigned i=0;i<info.dofsB[np]; ++i)
        {
            globI = m_locDofsMapper[np][Base::getComp(np,m_boundDofs[np][i])].index(Base::compCalcBack(np,m_boundDofs[np][i]));

            if(m_globIsPrimalIndex[np][globI])
                w.row(i) = xPi.row(m_glob2PrimalRemainingIndex[np][globI]);
            else
                w.row(i) = x2.row(m_glob2PrimalRemainingIndex[np][globI]);


        }

    }

}
*/
template<class T>
void gsIETIAssemblerMPI<T>::embedding(gsMatrix<T> & xPi, gsMatrix<T> const & x2, size_t np, gsMatrix<T>& w) const
{
    int nRhs = x2.cols();
    w.setZero(info.dofsB[np],nRhs);
    //  gsInfo<<"rank: "<<m_comm.rank()<<"  xPi on "<<np<<" is "<<xPi.transpose()<<"\n";
    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
    {
        w=x2;
        if(info.dofsP[np]!=0)
            w.noalias()+=m_Phi[np]*xPi;
    }
    else
    {
        unsigned globI;

        for (index_t i=0;i<info.dofsB[np]; ++i)
        {
            globI = m_locDofsMapper[np][Base::getComp(np,m_boundDofs[np][i])].index(Base::compCalcBack(np,m_boundDofs[np][i]));

            if(m_globIsPrimalIndex[np][globI])
                w.row(i) = xPi.row(m_glob2PrimalRemainingIndex[np][globI]);
            else
                w.row(i) = x2.row(m_glob2PrimalRemainingIndex[np][globI]);

        }
    }
}

//w has size info.dofsB[np]
template<class T>
void gsIETIAssemblerMPI<T>::embedding(std::vector<gsMatrix<T> > & xPis,std::vector<gsMatrix<T> >const & x2,std::vector<gsMatrix<T> >& w) const
{
    size_t np;

    for (size_t i=0; i<m_patchIdx.size();++i)
    {
        //gsInfo<<"Rank: "<<m_comm.rank()<<": "<<w[m_patchIdx[i]].rows()<<"  "<<w[m_patchIdx[i]].cols()<<" -- "<<xPis[i].rows()<< "  "<<xPis[i].cols()<< " -- "<<m_Phi[m_patchIdx[i]].rows()<< "  "<<m_Phi[m_patchIdx[i]].cols()<<"\n\n"<<std::flush;

        np=m_patchIdx[i], embedding(xPis[i],x2[np],np,w[np]);
     }
}


template<class T>
void gsIETIAssemblerMPI<T>::initEmbedding(const gsIETIInfo& info, gsMatrix<T>& xP, std::vector<gsMatrix<T> > & xPis)const
{
    xPis.resize(m_patchIdx.size());

    if(m_comm.size()==1 || m_IETIoptions.nSppHolder==m_comm.size())
    {
        for(size_t j=0; j<MPI_Base::m_proc2patch[m_comm.rank()].size();++j)
        {
            size_t np = MPI_Base::m_proc2patch[m_comm.rank()][j];
            Base::distributePrimal(xP,np,xPis[j]);
        }
        return;
    }
    gsMatrix<T> xP_reord;
    size_t nP =0;
    size_t iter =0;
    for(size_t i=0; i<MPI_Base::m_embGroup.size();++i)
        nP+=MPI_Base::m_nPDofs_group[i];


    if(MPI_Base::isSppHolder())
    {
        //  gsInfo<<"rank: "<<m_comm.rank()<<"  xP is given by "<<xP.transpose()<<"\n";
        gsMatrix<T> xPi;
        xP_reord.setZero(nP,xP.cols());
        for(size_t i=0; i<MPI_Base::m_embGroup.size();++i)
        {
            int p=MPI_Base::m_embGroup[i];
            for(size_t j=0; j<MPI_Base::m_proc2patch[p].size();++j)
            {
                size_t np = MPI_Base::m_proc2patch[p][j];
                Base::distributePrimal(xP,np,xPi);
                xP_reord.block(iter,0,xPi.rows(),xPi.cols())=xPi;
                iter+=xPi.rows();
            }
        }
    }
    MPI_Base::initEmbedding_(info,xP_reord,xPis);
}

template<class T>
void gsIETIAssemblerMPI<T>::embeddingTransPartial(const gsMatrix<T> & w,size_t np, gsMatrix<T> &  u2) const
{

    int nRhs = w.cols();
    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
    {
        /*
        u2.setZero(info.dofsB[np],nRhs);

        gsMatrix<T> w_trans= w.transpose(); //to circumwent const qualifier

        //MPI_Base::m_embT_buffer[np] = w_trans * m_Phi[np]; //phi.trans()*w == (w.trans *Phi).trans()
        // MPI_Base::m_embT_buffer[np].transposeInPlace(); //back to normal
        w_trans= w_trans*m_Phi[np];
        w_trans.transposeInPlace();
        Base::assemblePrimal(w_trans,np,MPI_Base::m_embT_buffer);

        // MPI_Base::send_embT(info,np,nRhs);

        //  gsMatrix<T> result = (MPI_Base::m_embT_buffer[np].transpose() *m_C[np]);
        gsMatrix<T> result = (w_trans.transpose() *m_C[np]);
        result.transposeInPlace();
        u2 = w - result;
*/
        u2=w;
        gsMatrix<T> temp_trans = (m_Phi[np].topRows(u2.rows())).transpose()*u2;

        Base::assemblePrimal(temp_trans,np,MPI_Base::m_embT_buffer);

        u2.topRows(info.dofsB[np]) -= (m_C[np].transpose())*temp_trans;

    }
    else
    {
        //   MPI_Base::m_embT_buffer[np].setZero(info.dofsP[np],nRhs);

        u2.setZero(info.dofsR[np],nRhs);

        //   for(size_t k=0;k<m_primalVdofs[np].size();k++)
        //     MPI_Base::m_embT_buffer[np].row(k) += w.row( m_glob2BoundInteriorIndex[np][m_locDofsMapper[np][Base::getComp(np,m_primalVdofs[np][k])].index(Base::compCalcBack(np,m_primalVdofs[np][k]))]);

        //   MPI_Base::send_embT(info,np,nRhs);

        for(size_t i=0;i<m_remDofs[np].size();i++)
            if(m_remDofIsAlsoBoundDof[np][i])
                u2.row(i) = w.row(m_glob2BoundInteriorIndex[np][m_locDofsMapper[np][Base::getComp(np,m_remDofs[np][i])].index(Base::compCalcBack(np,m_remDofs[np][i]))]);

    }

}




//w has size info.dofsB[np]
template<class T>
void gsIETIAssemblerMPI<T>::embeddingTrans(const std::vector<gsMatrix<T> >& w,gsMatrix<T>&  uP, std::vector<gsMatrix<T> > &  u2) const
{
    int nRhs = w.front().cols();
    uP.setZero(info.dofTotalP,nRhs);
    //  gsInfo<<"Rank: "<<m_comm.rank()<<"All have waited!!!!!"<<std::endl;
    MPI_Base::init_embT(info);

    for (size_t i=0; i<m_patchIdx.size();++i)
    {
        size_t np=m_patchIdx[i];
        embeddingTransPartial(w[np],np,u2[np]);
        MPI_Base::send_embT(info,np,nRhs);
    }
}

} // namespace gismo

#endif
