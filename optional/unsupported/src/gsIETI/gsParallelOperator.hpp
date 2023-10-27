/**  gsParallelOperator.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on:  2017-10-23
*/

#include <gsCore/gsConfig.h>
#ifdef GISMO_WITH_MPI

#include  <gsIETI/gsParallelOperator.h>


namespace gismo {

template<typename T>
void gsConnectionHandler<T>::init(const std::vector<connectionPair> &pairs)
{

    size_t size = pairs.size();
    int rank = m_comm.rank();

    for(size_t i=0; i<size; ++i)
    {
        const procDof& PD1 = pairs[i].first;
        const procDof& PD2 = pairs[i].second;
        if(PD1.first == rank)
            m_neighbor.push_sorted_unique(PD2.first);
        else if(PD2.first == rank)
            m_neighbor.push_sorted_unique(PD1.first);
    }
    m_nNeighbors = m_neighbor.size();

    /*
    for(int i=0; i< size;++i)
    {
        gsInfo<<"rank: "<<gsMpi::worldRank()<<"("<<pairs[i].first.first<<","<<pairs[i].first.second<<")-("<<pairs[i].second.first<<","<<pairs[i].second.second<<");\n";
    }
*/


    std::vector<std::vector<int> > blocklengths(m_comm.size());
    std::vector<std::vector<int> > displacements(m_comm.size());

    
    m_buffToIdx.resize(m_comm.size());
    for(int i=0; i<m_comm.size();++i)
    {
        m_buffToIdx[i].reserve(pairs.size());
        // TODO: optimize bound
        blocklengths[i].reserve(pairs.size()/(m_nNeighbors+1)*3);
        displacements[i].reserve(pairs.size()/(m_nNeighbors+1)*3);
    }
    
    int iter = 0;
    // int proc = 0;

    std::vector<int> lastDof(m_comm.size());
    for(size_t i=0; i< size;++i)
    {
        const procDof& PD1 = pairs[i].first;
        const procDof& PD2 = pairs[i].second;
        if((PD1.first == rank) ^ (PD2.first == rank))
        {
            bool isPD1Neighb = m_neighbor.bContains(PD1.first);
            const procDof& PDother = isPD1Neighb ? PD1 : PD2 ;
            const procDof& PDme = isPD1Neighb ? PD2 : PD1 ;

            m_buffToIdx[PDother.first].push_back(PDme.second);
            if(i==0)
            {
                //proc = PDother.first;
                for(gsSortedVector<int>::iterator pIt=m_neighbor.begin();pIt!=m_neighbor.end();++pIt)
                {
                    int p = (*pIt);
                    lastDof[p] = -2;//dummy value, dof >= 0
                }
            }

            if(lastDof[PDother.first]+1== PDme.second)
            {
                blocklengths[PDother.first].back()++;
            }
            else
            {
                blocklengths[PDother.first].push_back(1);
                displacements[PDother.first].push_back(PDme.second);
            }
            lastDof[PDother.first] = PDme.second;
            iter ++;
        }
    }

    // gsInfo<<"rank: "<<gsMpi::worldRank()<<": finished data structures\n"<<std::flush;
    iter=0;
    m_sendTypes.resize(m_comm.size());


    /*
    // debug variant
    for(int pp = 0; pp< m_comm.size();++pp)
    {
        m_comm.barrier();
        if(pp==m_comm.rank())
        {

            for(gsSortedVector<int>::iterator pIt=m_neighbor.begin();pIt!=m_neighbor.end();++pIt)
            {
                int p = (*pIt);

                gsInfo<<"rank: "<<gsMpi::worldRank()<<", Blocklengths to "<<p<<" with size "<< blocklengths[p].size()<<" for "<<p<<": ";
                for(size_t i=0; i< blocklengths[p].size();i++)
                    std::cout<<blocklengths[p][i]<<" "<<std::flush;
                gsInfo<<"\n";
                gsInfo<<"rank: "<<gsMpi::worldRank()<<", Displacements to "<<p<<" with size "<< displacements[p].size()<<" for "<<p<<": ";
                for(size_t i=0; i< displacements[p].size();i++)
                    std::cout<<displacements[p][i]<<" "<<std::flush;
                gsInfo<<"\n";

                GISMO_ASSERT(blocklengths[p].size() == displacements[p].size(), "number of blocks do not match");

                MPI_Type_indexed(blocklengths[p].size(),blocklengths[p].data(),displacements[p].data(),MPITraits<T>::getType(),&(m_sendTypes[p]));
                MPI_Type_commit(&(m_sendTypes[p]));
                iter++;
            }
            gsInfo<<"rank: "<<gsMpi::worldRank()<<": finished MPI_types\n"<<std::flush;
        }
    }
*/

    for(gsSortedVector<int>::iterator pIt=m_neighbor.begin();pIt!=m_neighbor.end();++pIt)
    {
        int p = (*pIt);
        GISMO_ASSERT(blocklengths[p].size() == displacements[p].size(), "number of blocks do not match");
        MPI_Type_indexed(blocklengths[p].size(),blocklengths[p].data(),displacements[p].data(),MPITraits<T>::getType(),&(m_sendTypes[p]));
        MPI_Type_commit(&(m_sendTypes[p]));
        iter++;
    }
    //gsInfo<<"rank: "<<gsMpi::worldRank()<<": finished MPI types\n"<<std::flush;

    for(size_t i=0; i< size;++i)
    {
        const procDof& PD1 = pairs[i].first;
        const procDof& PD2 = pairs[i].second;
        if(PD1.first == rank || PD2.first == rank)
        {
            bool isPD1Neighb = m_neighbor.bContains(PD1.first);
            const procDof& PDme = isPD1Neighb ? PD2 : PD1 ;
            const procDof& PDother = isPD1Neighb ? PD1 : PD2 ;

            if(m_multiplicity.count(PDme.second) > 0)
                m_multiplicity[PDme.second]++;
            else
                m_multiplicity[PDme.second] = 2;

            m_connPairs[PDme.second].insert(PDother);
        }
    }
    //gsInfo<<"rank: "<<gsMpi::worldRank()<<": finished init\n"<<std::flush;
}

template<typename T>
void gsConnectionHandler<T>::updateInterfaces(const gsMatrix<T> & buff, int p, gsMatrix<T>& inout) const
{

    // gsInfo<<"rank "<<gsMpi::worldRank()<<": buff.size()="<<buff.size()<<", inout.size()="<<inout.size()<<"; values of buffToIdx:";
    for(size_t i=0; i<m_buffToIdx[p].size();++i)
        inout.row(m_buffToIdx[p][i]) += buff.row(i);
}

template<typename T>
void gsConnectionHandler<T>::accumulateDistributedVector(gsMatrix<T> & inout) const
{
    GISMO_ASSERT(inout.cols()==1, "gsConnectionHandler allows only send/recv of vectors with one column");
    /*
    gsInfo<<"rank "<<gsMpi::worldRank()<<": nNeighb= "<<m_nNeighbors<<" with entries ";
    for(size_t n = 0; n<m_neighbor.size();++n)
        gsInfo<<m_neighbor[n]<<", ";
    gsInfo<<"\n";
    gsInfo<<"nBuffers: "<<m_buffer.size()<<"\n";

    */
    for(size_t n = 0; n<m_nNeighbors;++n)
    {
        int p= m_neighbor[n];
        MPI_Irecv(m_buffer[p].data(),m_buffer[p].size(),MPITraits<T>::getType(),p,11,m_comm,&m_requ[n]);
    }
    for(size_t n = 0; n<m_nNeighbors;++n)
    {
        int p= m_neighbor[n];
        MPI_Isend(inout.data(),1,m_sendTypes[p],p,11,m_comm,&m_sRequ[n]);
    }
    //   gsInfo<<"rank "<<gsMpi::worldRank()<<": send and recv posted, waiting\n"<<std::flush;
    //First wait that the data is sended
    MPI_Waitall(m_nNeighbors,m_sRequ,MPI_STATUS_IGNORE);
    
    //Then wait for the arrival of the data
    int idx = 0;
    for(size_t n = 0; n<m_nNeighbors;++n)
    {
        MPI_Waitany(m_nNeighbors,m_requ,&idx,MPI_STATUS_IGNORE);
        //     gsInfo<<"rank "<<gsMpi::worldRank()<<": got Index: "<<idx<<"\n"<<std::flush;
        int p = m_neighbor[idx];
        updateInterfaces(m_buffer[p],p,inout);
    }

}

template<typename T>
void gsConnectionHandler<T>::postAccumulate() const
{
    for(size_t n = 0; n<m_nNeighbors;++n)
    {
        int p= m_neighbor[n];
        MPI_Irecv(m_buffer[p].data(),m_buffer[p].size(),MPITraits<T>::getType(),p,111,m_comm,&m_requ[n]);
    }
}

template<typename T>
void gsConnectionHandler<T>::startAccumulate(const  gsMatrix<T> & input) const
{
    GISMO_ASSERT(input.cols()==1, "gsConnectionHandler allows only send/recv of vectors with one column");
    for(size_t n = 0; n<m_nNeighbors;++n)
    {
        int p= m_neighbor[n];
        MPI_Isend(const_cast<gsMatrix<T>&>(input).data(),1,m_sendTypes[p],p,111,m_comm,&m_sRequ[n]);
    }
}
template<typename T>
void gsConnectionHandler<T>::finishAccumulate(gsMatrix<T> & result) const
{
    GISMO_ASSERT(result.cols()==1, "gsConnectionHandler allows only send/recv of vectors with one column");
    int idx = 0;
    MPI_Waitall(m_nNeighbors,m_sRequ,MPI_STATUS_IGNORE);
    for(size_t n = 0; n<m_nNeighbors;++n)
    {
        MPI_Waitany(m_nNeighbors,m_requ,&idx,MPI_STATUS_IGNORE);
        // gsDebug<<"Rank: "+util::to_string(m_comm.rank())+" Got neighbor: "+util::to_string(idx)+" out of "+util::to_string(m_nNeighbors)+"\n";
        int p = m_neighbor[idx];
        updateInterfaces(m_buffer[p],p,result);
    }
}


template<typename T>
void gsConnectionHandler<T>::distributeAccumulatedVector(gsMatrix<T> & inout) const
{
    for(std::map<unsigned,unsigned>::const_iterator it = m_multiplicity.begin(); it!=m_multiplicity.end();++it)
        inout.row(it->first) /= it->second;
}

template<typename T>
void gsConnectionHandler<T>::accumulateSparseMatrix(gsSparseMatrix<T>& inout) const
{
    //first we use a not very suffisticated format for representing a SparseMatrix, namely (i,j,val(i,j))
    std::vector<size_t> sizes(inout.rows(),0);
    for(int k=0; k<inout.outerSize();++k)
        for(typename gsSparseMatrix<T>::InnerIterator it(inout,k); it;++it)
        {
            sizes[it.row()]++;
        }

    std::vector<size_t> sendSize(m_nNeighbors,0);
    for(size_t i=0; i<m_nNeighbors;++i)
        for(size_t j =0; j<m_buffToIdx[m_neighbor[i]].size();++j)
            sendSize[i] += sizes[m_buffToIdx[m_neighbor[i]][j]];

    /*
    std::vector<size_t> sendSizeProc(m_comm.size(),0);
    std::vector<size_t> recvSizeProc(m_comm.size(),0);
    for(size_t i=0; i<m_nNeighbors;++i)
        sendSizeProc[m_neighbor[i]] = sendSize[i];
    MPI_Alltoall(sendSizeProc.data(),1,MPITraits<size_t>::getType(),
                 recvSizeProc.data(),1,MPITraits<size_t>::getType(),
                 m_comm);
*/

    //    std::vector<gsMatrix<T> > buffer(m_comm.size());
    std::vector<std::vector<Triplet> > buffer(m_comm.size());
    MPI_Request* requests = new MPI_Request[m_nNeighbors];
    //send the data

    
    // gsDebug<<"rank "+util::to_string(m_com m.rank())+" sendSize: "<<gsAsVector<size_t>(sendSize).transpose()<<"\n"<<std::flush;
    //m_comm.barrier();

    //std::vector<gsMatrix<T> > sendData(m_comm.size());
    std::vector<std::vector<Triplet> > sendData(m_comm.size());
    for(size_t i=0; i<m_nNeighbors;++i)
        sendData[m_neighbor[i]].reserve(sendSize[i]);

    //  gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  start collecting data.\n"<<std::flush;m_comm.barrier();

   // std::vector<index_t> indices(m_comm.size(),0);
    for(int k=0; k<inout.outerSize();++k)
        for(typename gsSparseMatrix<T>::InnerIterator it(inout,k); it;++it)
        {
            std::map<unsigned, std::set<procDof> >::const_iterator rIt = m_connPairs.find(it.row());
            if(rIt!=m_connPairs.end())
            {
                const std::set<procDof>& r= rIt->second;
                std::map<unsigned, std::set<procDof> >::const_iterator cIt = m_connPairs.find(it.col());
                if(cIt!=m_connPairs.end())
                {
                    const std::set<procDof>& c= cIt->second;
                    for(std::set<procDof>::const_iterator rI = r.begin();rI!=r.end();++rI)
                        for(std::set<procDof>::const_iterator cI = c.begin();cI!=c.end();++cI)
                        {
                            if(rI->first == m_comm.rank() || cI->first == m_comm.rank())
                                continue;
                            if(rI->first == cI->first)
                            {
                                //   sendData[rI->first](0,indices[rI->first])=rI->second;
                                //   sendData[rI->first](1,indices[rI->first])=cI->second;
                                //   sendData[rI->first](2,indices[rI->first])=it.value();
                                //   indices[rI->first]++;
                                Triplet t;
                                t.index[0] = rI->second;
                                t.index[1] = cI->second;
                                t.value = it.value();
                                sendData[rI->first].push_back(give(t));
                            }
                        }
                }
            }

        }
    /*
    for(int p =0; p<m_comm.size();++p)
    {
        m_comm.barrier();
        if(p==m_comm.rank())
        {
            for(int pp =0; pp<m_comm.size();++pp)
                for(size_t i=0; i<sendData[pp].size();++i)
                    gsInfo<<"Rank "+util::to_string(gsMpi::worldRank())+": sendData["<<pp<<"]: "<<sendData[pp][i].index[0]<<" - "<<sendData[pp][i].index[1]<<": "<<sendData[pp][i].value<<"\n"<<std::flush;
        }
    }
*/
    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  collected data for Sp mat sending.\n"<<std::flush;m_comm.barrier();
    std::vector<size_t> sendSizeProc(m_comm.size(),0);
    std::vector<size_t> recvSizeProc(m_comm.size(),0);
    for(size_t i=0; i<m_nNeighbors;++i)
        sendSizeProc[m_neighbor[i]] = sendData[m_neighbor[i]].size();
    MPI_Alltoall(sendSizeProc.data(),1,MPITraits<size_t>::getType(),
                 recvSizeProc.data(),1,MPITraits<size_t>::getType(),
                 m_comm);
    for(size_t n = 0; n<m_nNeighbors;++n)
    {
        int p= m_neighbor[n];

        buffer[p].resize(recvSizeProc[p]);
        MPI_Irecv(buffer[p].data(),buffer[p].size(),m_tripletType,p,21,m_comm,&requests[n]);
    }
    
    MPI_Request requ[m_nNeighbors];
    for(size_t n = 0; n<m_nNeighbors;++n)
    {
        int p= m_neighbor[n];
        MPI_Isend(sendData[p].data(),sendData[p].size(),m_tripletType,p,21,m_comm,&requ[n]);
    }
    //   gsInfo<<"rank "<<gsMpi::worldRank()<<": send and recv posted, waiting\n"<<std::flush;
    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  sent the data.\n"<<std::flush;m_comm.barrier();
    
    /*
    for(int p =0; p<m_comm.size();++p)
    {
        m_comm.barrier();
        if(p==m_comm.rank())
        {
            gsInfo<<"Rank "+util::to_string(gsMpi::worldRank())+": sparse Matrix before: \n"<<inout.toDense()<<"\n";
        }
    }
    */

    //Wait for the arrival
    int idx = 0;
    for(size_t n = 0; n<m_nNeighbors;++n)
    {
        MPI_Waitany(m_nNeighbors,requests,&idx,MPI_STATUS_IGNORE);
        //     gsInfo<<"rank "<<gsMpi::worldRank()<<": got Index: "<<idx<<"\n"<<std::flush;
        int p = m_neighbor[idx];

        for(size_t i=0; i<recvSizeProc[p];++i)
            inout.addTo(buffer[p][i].index[0],buffer[p][i].index[1],buffer[p][i].value);
    }

    /*
    for(int p =0; p<m_comm.size();++p)
    {
        m_comm.barrier();
        if(p==m_comm.rank())
        {
            for(int pp =0; pp<m_comm.size();++pp)
                for(size_t i=0; i<recvSizeProc[pp];++i)
                    gsInfo<<"Rank "+util::to_string(gsMpi::worldRank())+": buffer["<<pp<<"]: "<<buffer[pp][i].index[0]<<" - "<<buffer[pp][i].index[1]<<": "<<buffer[pp][i].value<<"\n";
           // gsInfo<<"Rank "+util::to_string(gsMpi::worldRank())+": sparse Matrix: \n"<<inout.toDense()<<"\n";
            gsInfo<<"\n\n";
        }
    }
    */
    // gsDebug<<"Rank "+util::to_string(gsMpi::worldRank())+":  finished the data.\n"<<std::flush;m_comm.barrier();
    MPI_Waitall(m_nNeighbors,requ,MPI_STATUS_IGNORE); //leave the function only if all the data is sent.
    delete[] requests;
}

template<typename T>
gsMatrix<T> gsParallelOperator<T>::toDenseMatrix(gsParallelGlobalLocalHandler & handler) const
{
    gsMatrix<T> id = gsMatrix<T>::Identity(handler.globalSize(),handler.globalSize());
    gsMatrix<T> localMat, result,col;
    handler.extractLocalVector(id,localMat);
    for(int c=0; c<localMat.cols();++c)
    {
        m_locOp->apply(localMat.col(c),col);
        localMat.col(c) = col;
    }

    handler.buildGlobalVector(localMat,result);
    return result;

}

} // namespace gismo
#endif
