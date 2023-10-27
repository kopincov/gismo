/** @file gsParallelOperator.h

    @brief Implementation of a gsDistributedOperator, where each processor holds an
    gsLinearOperator.


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2017-10-23
*/


#pragma once

#include <gsCore/gsConfig.h>


#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>

#include <gsUtils/gsSortedVector.h>
#include <gsSolver/gsLinearOperator.h>
#include <gsCore/gsDofMapper.h>
#include <set>
#include <gsIETI/gsDistributedOperator.h>
#include <gsCore/gsBoxTopology.h>
#include <gsSolver/gsMatrixOp.h>

#include <stddef.h> //offsetof

namespace gismo {

template<typename T>
class gsConnectionHandler
{
public:
    struct Triplet
    {
        index_t index[2];
        T value;
    };
public:
    typedef memory::shared_ptr<gsConnectionHandler> Ptr;
    typedef memory::unique_ptr<gsConnectionHandler> uPtr;

    typedef std::pair<int,index_t> procDof;
    typedef std::pair<procDof,procDof> connectionPair;

    
    
    gsConnectionHandler(const std::vector<connectionPair> & pairs, MPI_Comm comm) : m_comm(comm)
    {
        init(pairs);
        m_buffer.resize(m_comm.size());
        for(size_t i=0; i<m_buffer.size();++i)
            m_buffer[i].resize(m_buffToIdx[i].size(),1);
        m_requ = new MPI_Request[m_nNeighbors];
        m_sRequ = new MPI_Request[m_nNeighbors];


        MPI_Aint offsets[2];
        offsets[0] = offsetof(Triplet, index) ; //most likely going to be 0
        offsets[1] = offsetof(Triplet, value) ;
        MPI_Datatype types[2];
        types[0] = MPITraits<index_t>::getType();
        types[1] = MPITraits<T>::getType();
        int lengths[2];
        lengths[0] = 2;
        lengths[1] = 1;


        MPI_Datatype tempType;
        MPI_Aint lb, extent;
        MPI_Type_create_struct(2, lengths, offsets, types,&tempType);

        MPI_Type_get_extent( tempType, &lb, &extent );
        MPI_Type_create_resized( tempType, lb, extent, &m_tripletType );
        MPI_Type_commit(&m_tripletType);
        MPI_Type_free(&tempType);
        //  gsInfo<<"rank: "<<gsMpi::worldRank()<<" constructor of gsConnectionHandler finished\n"<<std::flush;

    }
    static uPtr make(const std::vector<connectionPair> & pairs, MPI_Comm comm)
    {
        return memory::make_unique(new gsConnectionHandler(pairs, comm));
    }

    ~gsConnectionHandler()
    {
        delete[] m_requ;
        delete[] m_sRequ;
        MPI_Type_free(&m_tripletType);
    }

    void accumulateDistributedVector(gsMatrix<T> & inout) const;
    void distributeAccumulatedVector(gsMatrix<T> & inout) const;

    void accumulateSparseMatrix(gsSparseMatrix<T>& inout) const;

    void postAccumulate() const;
    void startAccumulate(const  gsMatrix<T> & input) const;
    void finishAccumulate(gsMatrix<T> & result) const;
    gsMpiComm getComm() const {return m_comm;}

protected:
    void init(const std::vector<connectionPair> & pairs);

    inline void updateInterfaces(const gsMatrix<T> & buff, int p, gsMatrix<T>& inout) const;

    /*
    struct pairComperator
    {
      explicit pairComperator(int i): n(i) { }
      inline bool operator()(const std::pair<index_t, int> & m) const { return m.first == n; }
    private:
      int n;
    };
*/
    MPI_Datatype getTripletType() const {return m_tripletType;}

protected:
    gsMpiComm m_comm;
    mutable std::vector<MPI_Datatype> m_sendTypes;
    MPI_Datatype m_tripletType;

    std::map<unsigned, std::set<procDof> > m_connPairs;

    size_t m_nNeighbors;
    gsSortedVector<int> m_neighbor;
    std::vector<std::vector<index_t> > m_buffToIdx;
    std::map<unsigned,unsigned> m_multiplicity;

    mutable std::vector<gsMatrix<T> > m_buffer;
    mutable MPI_Request* m_requ;
    mutable MPI_Request* m_sRequ;
};


class gsParallelGlobalLocalHandler;

template<typename T>
class gsParallelOperator : public gsDistributedOperator<T>
{
public:
    typedef memory::shared_ptr<gsParallelOperator> Ptr;
    typedef memory::unique_ptr<gsParallelOperator> uPtr;

    typedef std::pair<int,index_t> procDof;
    typedef std::pair<procDof,procDof> connectionPair;

public:

    gsParallelOperator(const std::vector<connectionPair> & pairsTestTrial, typename gsLinearOperator<T>::Ptr localOperator , MPI_Comm comm) : m_locOp(localOperator)
    {
        m_TestSpaceConnection = memory::make_shared(new gsConnectionHandler<T>(pairsTestTrial, comm));
        m_TrialSpaceConnection = m_TestSpaceConnection;
        m_accumulateAfter = false;
    }

    gsParallelOperator(const std::vector<connectionPair> & pairsTest,const std::vector<connectionPair> & pairsTrial, typename gsLinearOperator<T>::Ptr localOperator , MPI_Comm comm) : m_locOp(localOperator)
    {
        m_TestSpaceConnection = memory::make_shared(new gsConnectionHandler<T>(pairsTest, comm));
        m_TrialSpaceConnection = memory::make_shared(new gsConnectionHandler<T>(pairsTrial, comm));
        m_accumulateAfter = false;
    }

    gsParallelOperator(typename gsConnectionHandler<T>::Ptr connectionTestSpace,typename gsConnectionHandler<T>::Ptr connectionTrialSpace,  typename gsLinearOperator<T>::Ptr localOperator ) :
        m_locOp(localOperator), m_TestSpaceConnection(connectionTestSpace), m_TrialSpaceConnection(connectionTrialSpace)
    {
        m_accumulateAfter = false;
    }

    static uPtr make(const std::vector<connectionPair> & pairsTestTrial, typename gsLinearOperator<T>::Ptr localOperator , MPI_Comm comm)
    {
        return memory::make_unique(new gsParallelOperator(pairsTestTrial, localOperator,comm));
    }
    static uPtr make(const std::vector<connectionPair> & pairsTrial,const std::vector<connectionPair> & pairsTest, typename gsLinearOperator<T>::Ptr localOperator , MPI_Comm comm)
    {
        return memory::make_unique(new gsParallelOperator(pairsTrial, pairsTest,localOperator,comm));
    }
    static uPtr make(typename gsConnectionHandler<T>::Ptr connectionTrialSpace,typename gsConnectionHandler<T>::Ptr connectionTestSpace,  typename gsLinearOperator<T>::Ptr localOperator )
    {
        return memory::make_unique(new gsParallelOperator(connectionTrialSpace, connectionTestSpace,localOperator));
    }

    //TODO: Maybe this should go away...
    void accumulateAfterApplication(bool accumulate) {m_accumulateAfter = accumulate;}

    inline const typename gsConnectionHandler<T>::Ptr& getConnectionHandlerTestSpace() const {return m_TestSpaceConnection;}
    inline const typename gsConnectionHandler<T>::Ptr& getConnectionHandlerTrialSpace() const {return m_TrialSpaceConnection;}

    void postAccumulate() const
    {
        getConnectionHandlerTrialSpace()->postAccumulate();
    }
    void startAccumulate(const  gsMatrix<T> & input) const
    {
        getConnectionHandlerTrialSpace()->startAccumulate(input);
    }
    void finishAccumulate(gsMatrix<T> & result) const
    {
        getConnectionHandlerTrialSpace()->finishAccumulate(result);
    }
    void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
    {
        distributed = input;
        getConnectionHandlerTrialSpace()->distributeAccumulatedVector(distributed);
    }
    void accumulate(const gsMatrix<T>& input, gsMatrix<T> & accumulated) const
    {
        accumulated = input;
        getConnectionHandlerTrialSpace()->accumulateDistributedVector(accumulated);
    }

    void accumulate(gsSparseMatrix<T>& inout) const
    {
        GISMO_ASSERT(inout.rows() == inout.cols(), "This method is currently only implemented for square matrices");
        getConnectionHandlerTrialSpace()->accumulateSparseMatrix(inout);

    }


    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
    {
        result.setZero(m_locOp->rows(),input.cols());
        m_locOp->apply(input,result);
        if(m_accumulateAfter) m_TrialSpaceConnection->accumulateDistributedVector(result);
    }

    virtual index_t rows() const {return m_locOp->rows();}
    virtual index_t cols() const {return m_locOp->cols();}

    const typename gsLinearOperator<T>::Ptr getLinearOperator() const {return m_locOp;}

    gsMatrix<T> toDenseMatrix(gsParallelGlobalLocalHandler & handler) const;


protected:
    typename gsLinearOperator<T>::Ptr m_locOp;

    typename gsConnectionHandler<T>::Ptr m_TestSpaceConnection;
    typename gsConnectionHandler<T>::Ptr m_TrialSpaceConnection;
    bool m_accumulateAfter;
};






template<typename T>
class gsInterfaceConnection
{
protected:
    struct gsConnEntry
    {
        unsigned proc;
        unsigned glob;
        unsigned loc;
    };
    // moved to unsigned[3]

public:
    typedef memory::shared_ptr<gsInterfaceConnection> Ptr;
    typedef memory::unique_ptr<gsInterfaceConnection> uPtr;

    typedef std::pair<int,index_t> procDof;
    typedef std::pair<procDof,procDof> connectionPair;

public:

    gsInterfaceConnection(int rank): m_rank(rank) {}

    virtual ~gsInterfaceConnection()
    {
        MPI_Type_free(&m_connType);
    }

    void init()
    {
        //int blocklengths[3] = {1,1,1};
        MPI_Aint offsets[3];
        offsets[0] = offsetof(gsConnEntry, proc) ; //most likely going to be 0
        offsets[1] = offsetof(gsConnEntry, glob);
        offsets[2] = offsetof(gsConnEntry, loc);
        MPI_Datatype types[3];
        types[0] = MPITraits<unsigned>::getType();
        types[1] = MPITraits<unsigned>::getType();
        types[2] = MPITraits<unsigned>::getType();
        int lengths[3];
        lengths[0] = 1;
        lengths[1] = 1;
        lengths[2] = 1;
        MPI_Datatype tempType;
        MPI_Aint lb, extent;
        MPI_Type_create_struct(3, lengths, offsets, types,&tempType);

        MPI_Type_get_extent( tempType, &lb, &extent );
        MPI_Type_create_resized( tempType, lb, extent, &m_connType );

        // MPI_Type_contiguous(3,MPITraits<unsigned>::getType(),&m_connType);
        MPI_Type_commit(&m_connType);
        MPI_Type_free(&tempType);

        std::map<unsigned, std::set<procDof> > connMap;
        generateGlobalConnectionMap(connMap);
        initConnectionPair(connMap);
    }

    const std::vector<connectionPair>& getConnectionPairs() const {return m_conn;}

protected:

    virtual void generateGlobalConnectionMap(std::map<unsigned, std::set<procDof> >& connMap) = 0;

    void initConnectionPair(const std::map<unsigned, std::set<procDof> >& connMap)
    {
        //TODO: handle hanging Dirichlet-Neumann nodes

        //TODO: reserve m_conn;
        m_conn.reserve(connMap.size());

        for(std::map<unsigned,std::set<procDof> >::const_iterator it=connMap.begin(); it!=connMap.end(); ++it)
        {
            for(std::set<procDof>::iterator itA=(*it).second.begin(); itA!=(*it).second.end(); ++itA)
            {
                if(itA==(*it).second.end())
                    break;

                for(std::set<procDof>::iterator itB=itA; itB!=(*it).second.end(); ++itB)
                {
                    if(itB==itA)
                        continue;

                    if((itA->first == m_rank) ^ (itB->first == m_rank))
                        m_conn.push_back(connectionPair(*itA,*itB));
                }
            }
        }
    }

protected:
    int m_rank;
    std::vector<connectionPair> m_conn;

    MPI_Datatype m_connType;
};

class gsParallelGlobalLocalHandler
{
public:
    typedef memory::shared_ptr<gsParallelGlobalLocalHandler> Ptr;
    typedef memory::unique_ptr<gsParallelGlobalLocalHandler> uPtr;

    typedef std::pair<int,index_t> procDof;
    typedef std::pair<procDof,procDof> connectionPair;


public:
    gsParallelGlobalLocalHandler() {}

    gsParallelGlobalLocalHandler(gsDofMapper globMapper, std::vector<gsDofMapper> locMappers, MPI_Comm comm) : m_globMapper(globMapper) , m_locMappers(locMappers), m_comm(comm)
    {
        // gsInfo<<"rank: "<<gsMpi::worldRank()<<"; size= "<<m_locMappers.size()<<"\n";
    }

    static uPtr make(gsDofMapper globMapper, std::vector<gsDofMapper> locMappers, MPI_Comm comm)
    {
        return memory::make_unique(new gsParallelGlobalLocalHandler(globMapper, locMappers, comm));
    }

    template<typename U>
    void extractLocalVector(const gsMatrix<U> & globVector, gsMatrix<U>& locVector, int proc=-1) const
    {

        //  gsInfo<<"rank "<<gsMpi::worldRank()<<"\n new glob Vector: "<<globVec.transpose()<<"\n";
        /*
        for(int p =0; p<m_comm.size();++p)
        {
            m_comm.barrier();
            if(p==m_comm.rank())
            {
                if(proc ==-1) proc = m_comm.rank();
                gsDebug<< m_globMapper.asVector().transpose()<<"\n\n";
                gsDebug<< m_locMappers[proc].asVector().transpose()<<"\n\n";
                m_globMapper.print(gsDebug);
                m_locMappers[proc].print(gsDebug);


                gsInfo<<"Rank "<<gsMpi::worldRank()<<"\n";
                for(size_t i=0; i<m_locMappers[proc].mapSize();++i)
                    gsInfo<<" "<<m_globMapper.index(i,proc);
                gsInfo<<"\n";

                locVector.setZero(m_locMappers[proc].freeSize(),globVector.cols());
                gsInfo<<"Local Vector: "<<locVector.rows()<<" - globalVector: "<<globVector.rows()<<"\n"<<std::flush;
                for(size_t i=0; i<m_locMappers[proc].mapSize();++i)
                    if(m_locMappers[proc].is_free(i))
                        locVector.row(m_locMappers[proc].index(i))= globVector.row(m_globMapper.index(i,proc));
            }
        }
        */

        if(proc ==-1) proc = m_comm.rank();
        locVector.setZero(m_locMappers[proc].freeSize(),globVector.cols());
        for(size_t i=0; i<m_locMappers[proc].mapSize();++i)
            if(m_locMappers[proc].is_free(i))
                locVector.row(m_locMappers[proc].index(i))= globVector.row(m_globMapper.index(i,proc));
    }

    template<typename U>
    void addLocalVectorToGlobal(const gsMatrix<U> & locVector, gsMatrix<U>& globVector, int proc=-1) const
    {
        /*
        // debug variant
        for(int p =0; p<m_comm.size();++p)
        {
        //    m_comm.barrier();
            if(p==m_comm.rank())
            {
                if(proc ==-1) proc = m_comm.rank();
                gsDebugVar("Adding local vector");
                gsInfo<<"Rank: "<<gsMpi::worldRank()<<"; locVector: "<<locVector.rows()<<", globVector: "<<globVector.rows()<<"\n";

                std::vector<bool> visited(m_locMappers[proc].coupledSize(),false);
       //         gsInfo<<"loc:\n";m_locMappers[proc].print();
      //          gsInfo<<"glob:\n";m_globMapper.print();

                for(size_t i=0; i<m_locMappers[proc].mapSize();++i)
                {
                    if(m_locMappers[proc].is_coupled(i))
                        gsInfo<<"Rank: "<<gsMpi::worldRank()<<" ("<<i<<", cidx = "<<m_locMappers[proc].cindex(i)<<")\n"<<std::flush;
                    else if(m_locMappers[proc].is_free(i))
                        gsInfo<<"Rank: "<<gsMpi::worldRank()<<" ("<<i<<", fidx = "<<m_locMappers[proc].index(i)<<")\n"<<std::flush;

                    if(m_locMappers[proc].is_free(i) && !m_locMappers[proc].is_coupled(i))
                        globVector.row(m_globMapper.index(i,proc))+=locVector.row(m_locMappers[proc].index(i)) ;

                    if(m_locMappers[proc].is_coupled(i) && visited[m_locMappers[proc].cindex(i)] == false)
                    {
                        globVector.row(m_globMapper.index(i,proc))+=locVector.row(m_locMappers[proc].index(i)) ;
                        visited[m_locMappers[proc].cindex(i)].flip();
                    }
                }
            }
        }
        */

        //The m_locMapper are permuted anyway, so use tagged for extracting the coupled indices.
        if(proc ==-1) proc = m_comm.rank();
        std::vector<bool> visited(m_locMappers[proc].taggedSize(),false);
        for(size_t i=0; i<m_locMappers[proc].mapSize();++i)
        {
            if(m_locMappers[proc].is_free(i) && !m_locMappers[proc].is_tagged(i))
                globVector.row(m_globMapper.index(i,proc))+=locVector.row(m_locMappers[proc].index(i)) ;

            if(m_locMappers[proc].is_tagged(i) && visited[m_locMappers[proc].tindex(i)] == false)
            {
                globVector.row(m_globMapper.index(i,proc))+=locVector.row(m_locMappers[proc].index(i)) ;
                visited[m_locMappers[proc].tindex(i)].flip();
            }
        }
    }

    template<typename U>
    void buildGlobalVector(const gsMatrix<U> & locVector, gsMatrix<U>& globVector) const
    {
        globVector.setZero(m_globMapper.freeSize(),locVector.cols());
        addLocalVectorToGlobal(locVector,globVector,m_comm.rank());

        m_comm.sum(globVector.data(),(int)globVector.rows()*globVector.cols());
    }

    index_t localSize() const {return m_locMappers[m_comm.rank()].freeSize();}
    index_t globalSize() const {return m_globMapper.freeSize();}
    
    gsMpiComm getComm() const {return m_comm;}

protected:

    gsDofMapper m_globMapper;
    std::vector<gsDofMapper> m_locMappers;
    //gsDofMapper m_locMappers;
    gsMpiComm m_comm;

};



template<typename T,bool a, bool b> class gsPatchSubassambledLocalOperator;

template<typename T>
class gsPatchSubassembledTopology
{
    template<typename,bool a,bool b> friend class gsPatchSubassambledLocalOperator;

public:
    typedef memory::shared_ptr<gsPatchSubassembledTopology> Ptr;
    typedef memory::unique_ptr<gsPatchSubassembledTopology> uPtr;

    typedef std::pair<int,index_t> procDof;
    typedef std::pair<procDof,procDof> connectionPair;

public:

    gsPatchSubassembledTopology(const gsSortedVector<size_t>& myPatches,
                                     const std::vector<gsMultiBasis<T> >& bases, gsMatrix<bool> basesInteraction, const std::vector<std::vector<gsDofMapper> >& locMappers
                                     , const std::vector<gsDofMapper>& globalMapper)
        : m_myPatches(myPatches)
    {
        init(bases, basesInteraction, locMappers,globalMapper,bases.front().topology());
    }

    gsPatchSubassembledTopology(const gsSortedVector<size_t>& myPatches,
                                     const gsMultiBasis<T> & bases, const std::vector<gsDofMapper> & locMappers, const gsDofMapper& globalMapper)
        : m_myPatches(myPatches)
    {
        std::vector<gsMultiBasis<T> > VecBases(1,bases);
        std::vector<std::vector<gsDofMapper> > vecMapper(1,locMappers);
        std::vector<gsDofMapper> vecGlobalMapper(1,globalMapper);
        gsMatrix<bool> basesInteraction(1,1);
        basesInteraction(0,0) = true;
        init(VecBases, basesInteraction, vecMapper,vecGlobalMapper,bases.topology());
    }


    static uPtr make(const gsSortedVector<size_t>& myPatches,
                     const std::vector<gsMultiBasis<T> >& bases, gsMatrix<bool> basesInteraction, const std::vector<std::vector<gsDofMapper> >& locMappers,const std::vector<gsDofMapper>& globalMapper)
    {
        return memory::make_unique(new gsPatchSubassembledTopology(myPatches, bases, basesInteraction, locMappers,globalMapper));
    }

    static uPtr make( const gsSortedVector<size_t>& myPatches,
                     const gsMultiBasis<T> & bases, const std::vector<gsDofMapper> & locMappers,const gsDofMapper& globalMapper)
    {
        return memory::make_unique(new gsPatchSubassembledTopology(myPatches, bases, locMappers,globalMapper));
    }

    inline void extractLocalVector(const gsMatrix<T> & globVector, gsMatrix<T>& locVector, int patch) const
    {
        extractLocalVector_(globVector,locVector,patch,m_locMappers[patch],m_globMapper_patch);
    }
    inline void addLocalVectorToGlobal(const gsMatrix<T> & locVector, gsMatrix<T>& globVector, int patch) const
    {
        addLocalVectorToGlobal_(locVector,globVector,patch,m_locMappers[patch],m_globMapper_patch);
    }
    inline void setLocalVectorToGlobal(const gsMatrix<T> & locVector, gsMatrix<T>& globVector, int patch) const
    {
        setLocalVectorToGlobal_(locVector,globVector,patch,m_locMappers[patch],m_globMapper_patch);
    }

    const gsDofMapper & getLocalMapper() const {return m_globMapper; }
    const gsDofMapper & getLocalPatchMapper() const {return m_globMapper_patch; }
    gsDofMapper & getLocalMapper() {return m_globMapper; }
    gsDofMapper & getLocalPatchMapper() {return m_globMapper_patch; }

    ///reorders the globals dofs according to \a globalMapper.
    void reorderLike(const gsDofMapper & globalMapper)
    {
        gsDofMapper* patchMapper = &m_globMapper_patch;
        gsVector<index_t> perm = gsVector<index_t>::LinSpaced(patchMapper->freeSize(),0,patchMapper->freeSize()-1);

        //    gsInfo<<"rank "<<gsMpi::worldRank()<<" starting with filling permutation: "<<"\n"; gsInfo<<perm.transpose()<<"\n";
        std::vector<std::vector<std::pair<index_t,index_t> > > preImage(globalMapper.freeSize());
        for(index_t i=0; i<globalMapper.freeSize();++i)
            preImage[i].reserve(2);

        // gsDebug <<" " << "calculate preimage "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);
        for(size_t i=0; i<globalMapper.mapSize();++i)
        {
            size_t gl = globalMapper.mapIndex(i);
            if(globalMapper.is_free_index(gl))
            {
                size_t patch =0;
                for(; patch<globalMapper.numPatches();++patch)
                    if(globalMapper.offset(patch) > i)
                        break;
                patch--;
                preImage[gl].push_back(std::make_pair(patch, i - globalMapper.offset(patch) - globalMapper.firstIndex()) );
            }
        }

        // gsDebug <<" " << "fill perm "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);
        int substract = 0;
        for(index_t i=0; i<globalMapper.freeSize();++i)
        {
            GISMO_ASSERT(globalMapper.is_free_index(i), "the index is not free, something went wrong.");

            const std::vector<std::pair<index_t,index_t> >& preIm = preImage[i];
            bool isContained = false;
            for(size_t j = 0; j<preIm.size();++j)
            {
                if(m_myPatches.bContains(preIm[j].first))
                {
                    index_t idx1 = i - substract;
                    index_t idx2  = patchMapper->freeIndex(preIm[j].second,preIm[j].first);
                    if(idx1 !=idx2)
                        perm[idx2] = idx1;
                    isContained = true;
                    break;
                }
            }
            if(!isContained)
                substract++;
        }

        /*
        for(int p =0; p<gsMpi::worldSize();++p)
        {
            MPI_Barrier(gsMpi::worldComm());
            if(p==gsMpi::worldRank())
            {
        gsInfo<<"m_globMapper Before perm: \n";
        for(size_t i=0; i<m_globMapper.mapSize();++i)
            gsInfo<<m_globMapper.mapIndex(i)<<" , ";
        gsInfo<<"\n perm: \n"<<perm.transpose()<<"\n";
        */
        // gsDebug <<" " << "start permutation "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);
            m_globMapper_patch.markCoupledAsTagged();
            m_globMapper_patch.permuteFreeDofs(perm);
            m_globMapper.markCoupledAsTagged();
            m_globMapper.permuteFreeDofs(perm);

            //if(distributeInput)
            {
                //permute the multiplicity stuff
                std::map<size_t,size_t> multiplicityTemp;
                for(std::map<size_t,size_t>::const_iterator it = m_multiplicity.begin(); it!=m_multiplicity.end();++it)
                    multiplicityTemp[perm[it->first]]=it->second;
                m_multiplicity = multiplicityTemp;
            }

        // gsDebug <<" " << "finished permutation "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);

        /*
        gsInfo<<"\n\nm_globMapper After perm: \n";
        for(size_t i=0; i<m_globMapper.mapSize();++i)
            gsInfo<<m_globMapper.mapIndex(i)<<" , ";

        gsInfo<<"\n\nglobMapper: \n";
        for(size_t i=0; i<globalMapper.mapSize();++i)
            gsInfo<<globalMapper.mapIndex(i)<<" , ";
        gsInfo<<"\n\n";
            }
        }
        */
    }


    void reorder(const gsDofMapper & dofmapperRef, gsDofMapper & dofmapper) const
    {
        gsDofMapper* patchMapper = &dofmapper;
        gsVector<index_t> perm = gsVector<index_t>::LinSpaced(patchMapper->freeSize(),0,patchMapper->freeSize()-1);

        //    gsInfo<<"rank "<<gsMpi::worldRank()<<" starting with filling permutation: "<<"\n"; gsInfo<<perm.transpose()<<"\n";
        std::vector<std::vector<std::pair<index_t,index_t> > > preImage(dofmapperRef.freeSize());
        for(index_t i=0; i<dofmapperRef.freeSize();++i)
            preImage[i].reserve(2);

        // gsDebug <<" " << "calculate preimage "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);
        for(size_t i=0; i<dofmapperRef.mapSize();++i)
        {
            size_t gl = dofmapperRef.mapIndex(i);
            if(dofmapperRef.is_free_index(gl))
            {
                size_t patch =0;
                for(; patch<dofmapperRef.numPatches();++patch)
                    if(dofmapperRef.offset(patch) > i)
                        break;
                patch--;
                preImage[gl].push_back(std::make_pair(patch, i - dofmapperRef.offset(patch) - dofmapperRef.firstIndex()) );
            }
        }

        // gsDebug <<" " << "fill perm "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);
        int substract = 0;
        gsSortedVector<size_t>::const_iterator it;
        for(index_t i=0; i<dofmapperRef.freeSize();++i)
        {
            GISMO_ASSERT(dofmapperRef.is_free_index(i), "the index is not free, something went wrong.");

            const std::vector<std::pair<index_t,index_t> >& preIm = preImage[i];
            bool isContained = false;
            for(size_t j = 0; j<preIm.size();++j)
            {
                it =m_myPatches.find_it_or_fail(preIm[j].first);
                if(it!=m_myPatches.end())
                {
                    index_t idx1 = i - substract;
                    index_t idx2  = patchMapper->freeIndex(preIm[j].second,std::distance(m_myPatches.begin(),it));
                    if(idx1 !=idx2)
                        perm[idx2] = idx1;
                    isContained = true;
                    break;
                }
            }
            if(!isContained)
                substract++;
        }

        /*
        for(int p =0; p<gsMpi::worldSize();++p)
        {
            MPI_Barrier(gsMpi::worldComm());
            if(p==gsMpi::worldRank())
            {
        gsInfo<<"m_globMapper Before perm: \n";
        for(size_t i=0; i<m_globMapper.mapSize();++i)
            gsInfo<<m_globMapper.mapIndex(i)<<" , ";
        gsInfo<<"\n perm: \n"<<perm.transpose()<<"\n";
        */
        // gsDebug <<" " << "start permutation "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);
            patchMapper->markCoupledAsTagged();
            patchMapper->permuteFreeDofs(perm);

        // gsDebug <<" " << "finished permutation "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);

        /*
        gsInfo<<"\n\nm_globMapper After perm: \n";
        for(size_t i=0; i<m_globMapper.mapSize();++i)
            gsInfo<<m_globMapper.mapIndex(i)<<" , ";

        gsInfo<<"\n\nglobMapper: \n";
        for(size_t i=0; i<globalMapper.mapSize();++i)
            gsInfo<<globalMapper.mapIndex(i)<<" , ";
        gsInfo<<"\n\n";
            }
        }
        */
    }


protected:
    void inline init(const std::vector<gsMultiBasis<T> > & bases,const gsMatrix<bool>& basesInteraction, const std::vector<std::vector<gsDofMapper> >& locMappers, const std::vector<gsDofMapper>& globalMapper,const gsBoxTopology& m_patches)
    {
        init_(bases,basesInteraction,locMappers,globalMapper,m_patches,m_locMappers,m_globMapper,m_globMapper_patch);
        //if(distributeInput)
            calculateMultiplicity();
    }

    void init_(const std::vector<gsMultiBasis<T> > & bases,const gsMatrix<bool>& basesInteraction, const std::vector<std::vector<gsDofMapper> >& locMappers, const std::vector<gsDofMapper>& globalMapper,const gsBoxTopology& m_patches,
               std::vector<gsDofMapper >& m_locMappers, gsDofMapper& m_globMapper, gsDofMapper& m_globMapper_patch)
    {
        if(locMappers.size() == 1)
            m_locMappers = locMappers.front();
        else
        {
            m_locMappers.resize(m_patches.nBoxes());
            for(size_t npi = 0; npi<m_myPatches.size();++npi)
            {
                size_t np = m_myPatches[npi];

                size_t size=0;
                for(size_t c =0; c<locMappers.size();++c)
                    size+=locMappers[c][np].size();
                gsVector<index_t> s(1);
                s<<size;

                m_locMappers[np] = gsDofMapper(s);
                size =0;
                for(size_t c =0; c<locMappers.size();++c)
                {
                    for(int i=0; i<locMappers[c][np].size();++i)
                        if(locMappers[c][np].is_boundary(i))
                            m_locMappers[np].eliminateDof(size+i,0);
                    size+=locMappers[c][np].size();
                }
                m_locMappers[np].finalize();
            }
        }


        gsVector<index_t> size(1),sizes(m_patches.nBoxes());
        sizes.setZero();
        size(0) = 0;
        for(size_t i=0; i<m_myPatches.size();++i)
        {
            sizes[m_myPatches[i]] = m_locMappers[m_myPatches[i]].size();
            size[0]+=m_locMappers[m_myPatches[i]].size();
        }

        m_globMapper = gsDofMapper(size);
        m_globMapper_patch = gsDofMapper(sizes);


        /*
        gsMatrix<unsigned>  b1,b2;
        for(typename gsMultiPatch<T>::const_iiterator iit = m_patches.iBegin();iit!=m_patches.iEnd();++iit)
        {
            boundaryInterface bI = *iit;

            if(m_myPatches.bContains(bI.first().patch) && m_myPatches.bContains(bI.second().patch))
            {
                size_t i1 = *m_myPatches.find_ptr_or_fail(bI.first().patch);
                size_t i2 = *m_myPatches.find_ptr_or_fail(bI.second().patch);
                size_t sizec1 = 0;
                for(size_t c1 = 0; c1<bases.size();++c1)
                {
                    size_t sizec2 = 0;
                    for(size_t c2 = 0;c2<bases.size();++c2)
                    {
                        if(basesInteraction(c1,c2))
                        {
                            bases[c1].basis(bI.first().patch).matchWith(bI,bases[c2].basis(bI.second().patch),b1,b2);
                            b1 = (b1.array() + sizec1).matrix();
                            b2 = (b2.array() + sizec2).matrix();
                            m_globMapper_patch.matchDofs(bI.first().patch,b1,bI.second().patch,b2);
                            // gsInfo<<"merge: "+util::to_string(bI.first().patch)+" - "+util::to_string(b1.transpose())+" and \n "+util::to_string(bI.second().patch)+ " - "+util::to_string(b2.transpose())+"\n\n";
                            b1 = (b1.array()+sizes.head(i1).sum()).matrix();
                            b2 = (b2.array()+sizes.head(i2).sum()).matrix();
                            m_globMapper.matchDofs(0,b1,0,b2);

                        }
                        sizec2+=locMappers[c2][bI.second().patch].size();
                    }
                    sizec1+=locMappers[c1][bI.first().patch].size();
                }

            }

        }
        */
        std::vector<std::vector<std::vector<std::pair<index_t,index_t> > > >preImage(bases.size());
       // gsInfo<<"Global Mapper before preimage:  \ncoupled="<<globalMapper[0].coupledSize()<<"  \ntagged: "<<globalMapper[0].taggedSize()<<"\n\n";
        for(size_t c=0; c<bases.size();++c)
        {
            preImage[c].resize(globalMapper[c].freeSize());
            for(index_t i=0; i<globalMapper[c].freeSize();++i)
                preImage[c][i].reserve(2);

            // gsDebug <<" " << "calculate preimage "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);
            for(size_t i=0; i<globalMapper[c].mapSize();++i)
            {
                size_t gl = globalMapper[c].mapIndex(i);
                if(globalMapper[c].is_coupled_index(gl) || globalMapper[c].is_tagged_index(gl))
                {
                    size_t patch =0;
                    for(; patch<globalMapper[c].numPatches();++patch)
                        if(globalMapper[c].offset(patch) > i)
                            break;
                    patch--;
                    if(m_myPatches.bContains(patch))
                        preImage[c][gl].push_back(std::make_pair(patch, i - globalMapper[c].offset(patch) - globalMapper[c].firstIndex()) );
                }
            }
        }
        //gsInfo<<"preimage size: "<<preImage[0].size()<<"\n";

        //TODO: extend to multiple components
        for(size_t i=0; i<preImage[0].size();++i)
        {
           // gsInfo<<"size: "<<preImage[0][i].size()<<"\n";
            if(preImage[0][i].size()>1)
            {
                size_t i1 = *m_myPatches.find_ptr_or_fail(preImage[0][i][0].first);
                for(size_t j=1; j<preImage[0][i].size();++j)
                {
                    m_globMapper_patch.matchDof(preImage[0][i][0].first,preImage[0][i][0].second,preImage[0][i][j].first,preImage[0][i][j].second);
                    size_t i2 = *m_myPatches.find_ptr_or_fail(preImage[0][i][j].first);
                    m_globMapper.matchDof(0,preImage[0][i][0].second+sizes.head(i1).sum(),0,preImage[0][i][j].second+sizes.head(i2).sum());
                }
            }
        }

        //TODO: include the multiple components in a right way
        for(size_t idx = 0; idx<m_myPatches.size();++idx)
        {
            for(size_t c = 0; c<bases.size();++c)
            {
                index_t shift = sizes.head(m_myPatches[idx]).sum();
                for(int i =0; i<bases[c][m_myPatches[idx]].size();++i)
                    if(m_locMappers[m_myPatches[idx]].is_boundary(i))
                    {
                        m_globMapper_patch.eliminateDof(i,m_myPatches[idx]);
                        //gsInfo<<"rank: "+util::to_string(gsMpi::worldRank())+": eliminate: "+util::to_string(m_myPatches[idx])+" - "+util::to_string(i)+"\n";
                        m_globMapper.eliminateDof(i+shift,0);
                    }
            }
        }
        m_globMapper_patch.finalize();
        m_globMapper.finalize();
    }



    inline void extractLocalVector_(const gsMatrix<T> & globVector, gsMatrix<T>& locVector, int patch, const gsDofMapper & mapper, const gsDofMapper & globalMapper) const
    {
        locVector.setZero(mapper.freeSize(),globVector.cols());
        for(size_t i=0; i<mapper.mapSize();++i)
        {
            if(mapper.is_free(i))
                locVector.row(mapper.index(i))= globVector.row(globalMapper.index(i,patch));
        }



    }

    inline void addLocalVectorToGlobal_(const gsMatrix<T> & locVector, gsMatrix<T>& globVector, int patch, const gsDofMapper & mapper,  const gsDofMapper & globalMapper) const
    {
            for(size_t i=0; i<mapper.mapSize();++i)
            {
                if(mapper.is_free(i))
                    globVector.row(globalMapper.index(i,patch))+=locVector.row(mapper.index(i)) ;
            }

    }
    inline void setLocalVectorToGlobal_(const gsMatrix<T> & locVector, gsMatrix<T>& globVector, int patch, const gsDofMapper & mapper,  const gsDofMapper & globalMapper) const
    {
        //the local operators are stored distributed
            for(size_t i=0; i<mapper.mapSize();++i)
            {
                if(mapper.is_free(i))
                    globVector.row(globalMapper.index(i,patch))=locVector.row(mapper.index(i)) ;
            }
    }


protected:
    void calculateMultiplicity()
    {
        std::vector<size_t> sizes(m_globMapper_patch.numPatches());
        sizes.back() = m_globMapper_patch.mapSize() - m_globMapper_patch.offset(sizes.size()-1);
        for(size_t i=0; i<m_globMapper_patch.numPatches()-1;++i)
            sizes[i] = m_globMapper_patch.offset(i+1)-m_globMapper_patch.offset(i);

        for(size_t k = 0; k< sizes.size();++k)
            for(size_t i=0; i<sizes[k];++i)
            {
                if(m_globMapper_patch.is_coupled(i,k))
                {
                    size_t idx= m_globMapper_patch.index(i,k);
                    if(m_multiplicity.count(idx) > 0)
                        m_multiplicity[idx]++;
                    else
                        m_multiplicity[idx] = 1;
                }
            }
    }

protected:
    gsSortedVector<size_t> m_myPatches;

    gsDofMapper m_globMapper, m_globMapper_patch;
    std::vector<gsDofMapper> m_locMappers;

    std::map<size_t, size_t > m_multiplicity; //only needed for col dofs
};


template<typename T, bool preventOutputAccumulation = false, bool distributeInput = false>
class gsPatchSubassambledLocalOperator : public gsLinearOperator<T>
{

public:
    typedef memory::shared_ptr<gsPatchSubassambledLocalOperator> Ptr;
    typedef memory::unique_ptr<gsPatchSubassambledLocalOperator> uPtr;

    typedef std::pair<int,index_t> procDof;
    typedef std::pair<procDof,procDof> connectionPair;

public:

    gsPatchSubassambledLocalOperator(std::vector<typename gsLinearOperator<T>::Ptr> locOps,typename gsPatchSubassembledTopology<T>::Ptr topology)
        : m_locOps(give(locOps)), m_topologyRow(topology), m_topologyCol(topology){}

    gsPatchSubassambledLocalOperator(std::vector<typename gsLinearOperator<T>::Ptr> locOps,typename gsPatchSubassembledTopology<T>::Ptr topologyRow,typename gsPatchSubassembledTopology<T>::Ptr topologyCol)
        : m_locOps(give(locOps)), m_topologyRow(topologyRow), m_topologyCol(topologyCol){}

    static uPtr make(std::vector<typename gsLinearOperator<T>::Ptr> locOps,typename gsPatchSubassembledTopology<T>::Ptr topology)
    {
        return memory::make_unique(new gsPatchSubassambledLocalOperator(locOps, topology));
    }

    static uPtr make(std::vector<typename gsLinearOperator<T>::Ptr> locOps,typename gsPatchSubassembledTopology<T>::Ptr topologyRow,typename gsPatchSubassembledTopology<T>::Ptr topologyCol)
    {
        return memory::make_unique(new gsPatchSubassambledLocalOperator(locOps, topologyRow,topologyCol));
    }

    void apply(const gsMatrix<T> &input, gsMatrix<T> &x) const
    {
        x.setZero(rows(),input.cols());



        // gsInfo<<"myPatches: ";for(size_t i=0; i<m_myPatches.size();++i) gsInfo<<m_myPatches[i]<<", "; gsInfo<<"\n";
        if(distributeInput)
        {
            gsMatrix<T>& input_ =const_cast<gsMatrix<T> &>(input);
            for(std::map<size_t,size_t>::const_iterator it = m_topologyCol->m_multiplicity.begin(); it!=m_topologyCol->m_multiplicity.end();++it)
                input_.row(it->first) /= it->second;
        }


        for(size_t i=0; i<m_topologyRow->m_myPatches.size();++i)
        {
            m_topologyCol->extractLocalVector(input,m_locVec,m_topologyCol->m_myPatches[i]);
            m_locOps[m_topologyCol->m_myPatches[i]]->apply(m_locVec,m_locRes);
            if(preventOutputAccumulation)
                m_topologyRow->setLocalVectorToGlobal(m_locRes,x,m_topologyRow->m_myPatches[i]);
            else
                m_topologyRow->addLocalVectorToGlobal(m_locRes,x,m_topologyRow->m_myPatches[i]);
        }

        //repair input
        if(distributeInput)
        {
            gsMatrix<T>& input_ =const_cast<gsMatrix<T> &>(input);
            for(std::map<size_t,size_t>::const_iterator it = m_topologyCol->m_multiplicity.begin(); it!=m_topologyCol->m_multiplicity.end();++it)
                input_.row(it->first) *= it->second;
        }
    }

    virtual index_t rows() const {return m_topologyRow->getLocalMapper().freeSize();}
    virtual index_t cols() const {return m_topologyCol->getLocalMapper().freeSize();}

    typename gsPatchSubassembledTopology<T>::Ptr getTopologyRow() {return m_topologyRow;}
    typename gsPatchSubassembledTopology<T>::Ptr getTopologyCol() {return m_topologyCol;}


    /// size of the matrix is this->rows(), this->cols()
    /// It requires the input operators to be matrix operators of type gsSparseMatrix<T>
    gsSparseMatrix<T> buildSparseMatrix(const gsSortedVector<index_t>& dofsR, const gsSortedVector<index_t>& dofsC) const
    {
        size_t totalNNZ = 0;
        const gsSortedVector<size_t>& myPatches =  m_topologyRow->m_myPatches;
        std::vector<const gsMatrixOp<gsSparseMatrix<T> >* > mats(m_locOps.size());
        for(size_t i=0; i< myPatches.size();++i)
        {
            const gsMatrixOp<gsSparseMatrix<T> >* op= dynamic_cast<const gsMatrixOp<gsSparseMatrix<T> >* >(m_locOps[myPatches[i]].get());
            GISMO_ASSERT(op!=NULL, "the underlying operators must be of type gsMatrixOp<gsSparseMatrix<T> >");
            mats[myPatches[i]] = op;
            totalNNZ+= op->matrix().nonZeros();
        }
        
        //TODO: this is maybe a bit much
        gsSparseMatrix<T> result(this->rows(),this->cols());
      //  result.reservePerColumn((index_t)totalNNZ/this->cols());


        //  std::vector<std::map<unsigned,index_t> > patchDofsR(m_locOps.size());
        //  std::vector<std::map<unsigned,index_t> > patchDofsC(m_locOps.size());

        /*
        for(size_t k=0; k<m_myPatches.size();++k)
        {
            size_t patch = m_myPatches[k];
            patchDofsR[patch].reserve(dofsR.size());
            patchDofsC[patch].reserve(dofsC.size());
        }
        */
        /*
        std::vector<std::pair<index_t,index_t> > preIm;
        for(size_t i=0; i< dofsC.size();++i)
        {
            m_globMapper_patchCol.preImage(dofsC[i],preIm);
            for(size_t j=0; j<preIm.size();++j)
                patchDofsC[preIm[j].first][dofsC[i]] =preIm[j].second;
        }
        for(size_t i=0; i< dofsR.size();++i)
        {
            m_globMapper_patchRow.preImage(dofsR[i],preIm);
            for(size_t j=0; j<preIm.size();++j)
                patchDofsR[preIm[j].first][dofsR[i]] = preIm[j].second;
        }
        */
        std::vector<gsVector<index_t> > patchDofsR(myPatches.size());
        for(size_t i=0; i< myPatches.size();++i)
            patchDofsR[i] = give(m_topologyRow->m_locMappers[myPatches[i]].inverseAsVector());

        gsSparseEntries<T> entries;
        entries.reserve(totalNNZ);
        for(size_t k=0; k<myPatches.size();++k)
        {
            size_t patch = myPatches[k];
            for(size_t i=0; i<m_topologyCol->m_locMappers[patch].mapSize();++i)
            {
                index_t idxC = m_topologyCol->m_locMappers[patch].index(i);
                index_t glC= m_topologyCol->m_globMapper_patch.index(i,patch);
                if(m_topologyCol->m_locMappers[patch].is_free_index(idxC) && dofsC.bContains(glC))
                {
                    for(typename gsSparseMatrix<T>::InnerIterator it(mats[patch]->matrix(),idxC);it;++it)
                    {
                        //   std::map<unsigned,index_t>::iterator iit = patchDofsR[patch].find(it.row());
                        index_t idxR = patchDofsR[k][it.row()];
                        //                 gsInfo<<"actual row: "<<it.row()<<"\n"<<std::flush;

                        index_t glR= m_topologyRow->m_globMapper_patch.index(idxR,patch);
                        if(dofsR.bContains(glR))
                            entries.add(glR,glC,it.value());
                            //result.coeffRef(glR,glC)+= it.value();
                    }
                }
            }
        }
        result.setFrom(entries);
        result.makeCompressed();
        return result;
    }


    gsSparseMatrix<T> buildSparseMatrix(const gsSortedVector<index_t>& dofs) const
    {
        return buildSparseMatrix(dofs, dofs);
    }

    gsSparseMatrix<T> buildSparseMatrix(bool sumUp = true) const
    {
        size_t totalNNZ = 0;
        const gsSortedVector<size_t>& myPatches =  m_topologyRow->m_myPatches;
        std::vector<const gsMatrixOp<gsSparseMatrix<T> >* > mats(m_locOps.size());
        for(size_t i=0; i< myPatches.size();++i)
        {
            const gsMatrixOp<gsSparseMatrix<T> >* op= dynamic_cast<const gsMatrixOp<gsSparseMatrix<T> >* >(m_locOps[myPatches[i]].get());
            GISMO_ASSERT(op!=NULL, "the underlying operators must be of type gsMatrixOp<gsSparseMatrix<T> >");
            mats[myPatches[i]] = op;
            totalNNZ+= op->matrix().nonZeros();
        }
        std::vector<gsVector<index_t> > patchDofsR(myPatches.size());
        std::vector<gsVector<index_t> > patchDofsC(myPatches.size());
        for(size_t i=0; i< myPatches.size();++i)
        {
            patchDofsR[i] = give(m_topologyRow->m_locMappers[myPatches[i]].inverseAsVector());
            patchDofsC[i] = give(m_topologyCol->m_locMappers[myPatches[i]].inverseAsVector());
        }


        gsSparseMatrix<T> result(this->rows(),this->cols());

        gsSparseEntries<T> entries;
        entries.reserve(totalNNZ);
        for(size_t k=0; k<myPatches.size();++k)
        {
            size_t patch = myPatches[k];
            for (int j=0; j<mats[patch].matrix().outerSize(); ++j)
              for (typename gsSparseMatrix<T>::InnerIterator it(mats[patch].matrix(),j); it; ++it)
              {
                  entries.add(m_topologyRow->m_globMapper_patch.index(patchDofsR[patch][it.row()],patch),
                          m_topologyCol->m_globMapper_patch.index(patchDofsC[patch][it.col()],patch),it.value());
              }
        }
        if(sumUp)
            result.setFrom(entries);
        else
            result.setFromTriplets(entries.begin(),entries.end(),take_first<T>());
        result.makeCompressed();
        return result;

    }

    const std::vector<typename gsLinearOperator<T>::Ptr> & getLocalOps() const {return m_locOps;}

protected:
    std::vector<typename gsLinearOperator<T>::Ptr> m_locOps;
    typename gsPatchSubassembledTopology<T>::Ptr m_topologyRow;
    typename gsPatchSubassembledTopology<T>::Ptr m_topologyCol;

    mutable gsMatrix<T> m_locVec, m_locRes;

private:
    template <typename U>
    struct take_first {
        U operator() (const U& a, const U&b) { return a; }
    };
};



template<typename T>
class gsPatchInterfaceConnections : public gsInterfaceConnection<T>
{
    typedef gsInterfaceConnection<T> Base;
    typedef typename Base::gsConnEntry gsConnEntry;
public:

    typedef memory::shared_ptr<gsPatchInterfaceConnections> Ptr;
    typedef memory::unique_ptr<gsPatchInterfaceConnections> uPtr;

    typedef typename Base::procDof procDof;
    typedef typename Base::connectionPair connectionPair;

    gsPatchInterfaceConnections(const gsSortedVector<size_t>& myPatches,
                                const std::vector<gsMultiBasis<T> >& bases, gsMatrix<bool> basesInteraction, const std::vector<std::vector<gsDofMapper> >& locMappers,
                                const gsDofMapper& locOpPatchMappers,const gsDofMapper& globalMapper, const gsMpiComm& comm  )
        : Base(comm.rank()),  m_patches(bases.front().topology()), m_myPatches(myPatches), m_bases(bases), m_basesInteraction(basesInteraction), m_locMappers(locMappers), m_globMapper(globalMapper),
          m_locOpPatchMappers(locOpPatchMappers), m_comm(comm)
    {

        m_patchToProc.resize(m_patches.nBoxes());
        for(size_t npi = 0; npi<m_myPatches.size();++npi)
            m_patchToProc[m_myPatches[npi]] = m_comm.rank();
        m_comm.sum<int>(m_patchToProc.data(),(int)m_patchToProc.size());
    }

    gsPatchInterfaceConnections(const gsSortedVector<size_t>& myPatches,
                                const gsMultiBasis<T>& bases,const std::vector<gsDofMapper>& locMappers, const gsDofMapper& locOpPatchMappers, const gsDofMapper& globalMapper, const gsMpiComm& comm)
        : Base(comm.rank()),  m_patches(bases.topology()), m_myPatches(myPatches), m_globMapper(globalMapper), m_locOpPatchMappers(locOpPatchMappers), m_comm(comm)
    {
        m_basesInteraction.setZero(1,1);
        m_basesInteraction(0,0)=true;
        m_locMappers.push_back(locMappers);
        m_bases.push_back(bases);

        m_patchToProc.resize(m_patches.nBoxes());
        for(size_t npi = 0; npi<m_myPatches.size();++npi)
            m_patchToProc[m_myPatches[npi]] = m_comm.rank();
        m_comm.sum<int>(m_patchToProc.data(),(int)m_patchToProc.size());
    }

    gsDofMapper generateProcGlobalMapper() const
    {
        gsMatrix<index_t>  b1,b2;
        gsVector<index_t> size;
        size.setZero(m_comm.size());

        //    for(size_t p=0; p< m_comm.size();++p)
        //        size(p)+= procLocalMappers[p].mapSize();

        for(index_t np = 0; np< m_patches.nBoxes();++np)
            for(size_t c = 0; c<m_bases.size();++c)
                size(m_patchToProc[np])+= m_locMappers[c][np].size();

        gsDofMapper procDofMapper(size);

        std::vector<size_t> accSum(m_patches.nBoxes(),0);
        for(index_t np = 0; np< m_patches.nBoxes();++np)
            for(index_t i=0; i<np;++i)
                if(m_patchToProc[i]==m_patchToProc[np])
                {
                    for(size_t c = 0; c<m_bases.size();++c)
                        accSum[np]+=m_locMappers[c][i].size();
                }



        for(typename gsMultiPatch<T>::const_iiterator iit = m_patches.iBegin();iit!=m_patches.iEnd();++iit)
        {
            boundaryInterface bI = *iit;
            int p1 = m_patchToProc[bI.first().patch];
            int p2 = m_patchToProc[bI.second().patch];


            size_t sizec1 = 0;
            for(size_t c1 = 0; c1<m_bases.size();++c1)
            {
                size_t sizec2 = 0;
                for(size_t c2 = 0;c2<m_bases.size();++c2)
                {
                    if(m_basesInteraction(c1,c2))
                    {
                        m_bases[c1].basis(bI.first().patch).matchWith(bI,m_bases[c2].basis(bI.second().patch),b1,b2);
                        for(index_t i=0; i<b1.rows();++i)
                        {
                            index_t procL1 = b1(i,0)+accSum[bI.first().patch]+sizec1;
                            index_t procL2 = b2(i,0)+accSum[bI.second().patch]+sizec2;

                            procDofMapper.matchDof(p1,procL1,p2,procL2);

                            if(!m_locMappers[c1][bI.first().patch].is_free(b1(i,0)))
                                procDofMapper.eliminateDof(procL1,p1);
                            if(!m_locMappers[c2][bI.second().patch].is_free(b2(i,0)))
                                procDofMapper.eliminateDof(procL2,p2);
                        }
                    }
                    sizec2 += m_locMappers[c2][bI.second().patch].size();
                }
                sizec1 += m_locMappers[c1][bI.first().patch].size();
            }
        }


        std::vector<std::vector<std::pair<index_t,index_t> > >preImage(m_globMapper.freeSize());
        for(index_t i=0; i<m_globMapper.freeSize();++i)
            preImage[i].reserve(2);

        // gsDebug <<" " << "calculate preimage "<<"\n"<<std::flush; MPI_Barrier(MPI_COMM_WORLD);
        for(size_t i=0; i<m_globMapper.mapSize();++i)
        {
            size_t gl = m_globMapper.mapIndex(i);
            if(m_globMapper.is_coupled_index(gl) || m_globMapper.is_tagged_index(gl))
            {
                size_t patch =0;
                for(; patch<m_globMapper.numPatches();++patch)
                    if(m_globMapper.offset(patch) > i)
                        break;
                patch--;
                if(m_myPatches.bContains(patch))
                    preImage[gl].push_back(std::make_pair(patch, i - m_globMapper.offset(patch) - m_globMapper.firstIndex()) );
            }
        }


        //TODO: extend to multiple components
        for(size_t i=0; i<preImage.size();++i)
        {
            if(preImage[i].size()>1)
            {
                int p1 = m_patchToProc[preImage[i][0].first];

                for(size_t j=1; j<preImage[i].size();++j)
                {
                    int p2 = m_patchToProc[preImage[i][j].first];
                    procDofMapper.matchDof(p1,preImage[i][0].second+accSum[preImage[i][0].first],p2,preImage[i][j].second+accSum[preImage[i][j].first]);
                }
            }
        }




        for(index_t np = 0; np< m_patches.nBoxes();++np)
        {
            size_t sizec = 0;
            for(size_t c = 0; c<m_bases.size();++c)
            {
                for(int i=0; i<m_locMappers[c][np].size();++i)
                    if(!m_locMappers[c][np].is_free(i))
                        procDofMapper.eliminateDof(i+sizec+accSum[np],m_patchToProc[np]);
                sizec+=m_locMappers[c][np].size();
            }
        }
        /*
        for(int p=0; p<m_comm.size();++p)
            for(int i=0; i<procLocalMappers[p].mapSize();++i)
            {
                if(!procLocalMappers[p].is_free(i))
                    procDofMapper.eliminateDof(i,p);
            }
*/
        procDofMapper.finalize();
        gsVector<index_t> perm = generatePermutation(procDofMapper);
        procDofMapper.markCoupledAsTagged();
        procDofMapper.permuteFreeDofs(perm);

        return procDofMapper;
    }



protected:
    gsVector<index_t> generatePermutation(const gsDofMapper& procGlobalMapper)const
    {
        std::vector<size_t> accSum(m_patches.nBoxes(),0);
        for(index_t np = 0; np< m_patches.nBoxes();++np)
            for(index_t i=0; i<np;++i)
                if(m_patchToProc[i]==m_patchToProc[np])
                {
                    for(size_t c = 0; c<m_bases.size();++c)
                        accSum[np]+=m_locMappers[c][i].size();
                }

        gsVector<index_t> perm = gsVector<index_t>::LinSpaced(procGlobalMapper.freeSize(),0,procGlobalMapper.freeSize()-1);
        //    gsInfo<<"rank "<<gsMpi::worldRank()<<" starting with filling permutation: "<<"\n"; gsInfo<<perm.transpose()<<"\n";
        for(index_t np =0; np< m_patches.nBoxes();++np)
        {
            size_t sizec = 0;
            for(size_t c=0; c<m_bases.size();++c)
            {
                for(int i=0; i<m_locMappers[c][np].size();++i)
                {
                    if(m_locMappers[c][np].is_free(i))
                    {
                        index_t idx1 = m_globMapper.index(i+sizec,np);
                        index_t idx2  = procGlobalMapper.index(i+accSum[np]+sizec,m_patchToProc[np]);
                        if(idx1 !=idx2)
                            perm[idx2] = idx1;
                    }
                }
                sizec+= m_locMappers[c][np].size();
            }
        }
        //  gsInfo<<"Rank "<<gsMpi::worldRank()<<": "<<perm.transpose()<<"\n";

        return perm;
    }

    void generateGlobalConnectionMap(std::map<unsigned, std::set<procDof> >& connMap)
    {
        gsMatrix<index_t>  b1,b2,glob;

        size_t size = 0;
        for(size_t c = 0; c<m_bases.size();++c)
        {
            for(size_t i=0; i<m_myPatches.size();++i)
            {
                size+=m_locMappers[c][m_myPatches[i]].coupledSize();
            }
        }

        std::vector<gsConnEntry > locConnection;
        locConnection.reserve(2*size);
        std::vector<int> connSize(m_comm.size(),0);
        gsConnEntry connEntry;

        for(typename gsMultiPatch<T>::const_iiterator iit = m_patches.iBegin();iit!=m_patches.iEnd();++iit)
        {
            boundaryInterface bI = *iit;

            if(m_myPatches.bContains(bI.first().patch) | m_myPatches.bContains(bI.second().patch))
            {
                int proc1 = m_patchToProc[bI.first().patch];
                int proc2 = m_patchToProc[bI.second().patch];

                size_t sizec1 = 0;
                for(size_t c1 = 0; c1<m_bases.size();++c1)
                {
                    size_t sizec2 = 0;
                    for(size_t c2 = 0;c2<m_bases.size();++c2)
                    {
                        if(m_basesInteraction(c1,c2))
                        {
                            m_bases[c1].basis(bI.first().patch).matchWith(bI,m_bases[c2].basis(bI.second().patch),b1,b2);
                            m_globMapper.localToGlobal(b2,bI.second().patch,glob);
                            for(int i=0;i<glob.rows();i++)
                            {
                                if(m_globMapper.is_free_index(glob(i,0)))
                                {
                                    //glob proc loc

                                    connEntry.glob = glob(i,0);
                                    connEntry.proc = m_comm.rank();
                                    connEntry.loc = m_comm.rank() == proc1 ? m_locOpPatchMappers.index(b1(i,0)+sizec1,bI.first().patch) :
                                                                             m_locOpPatchMappers.index(b2(i,0)+sizec2,bI.second().patch);
                                    locConnection.push_back(connEntry);

                                    /*
                                    locConnection.push_back(glob(i,0));
                                    locConnection.push_back(m_comm.rank());
                                    locConnection.push_back(m_comm.rank() == proc1 ? m_locOpPatchMappers.index(b1(i,0)+sizec1,bI.first().patch) :
                                                                                     m_locOpPatchMappers.index(b2(i,0)+sizec2,bI.second().patch));
                                                                                     */
                                    if(proc1==proc2)
                                    {
                                        /*
                                        locConnection.push_back(glob(i,0));
                                        locConnection.push_back(m_comm.rank());
                                        locConnection.push_back(m_locOpPatchMappers.index(b2(i,0)+sizec2,bI.second().patch));
                                        */
                                        connEntry.loc = m_locOpPatchMappers.index(b2(i,0)+sizec2,bI.second().patch);
                                        locConnection.push_back(connEntry);
                                    }
                                }
                                else
                                {;} //TODO: Make the local mappers right, think about it
                            }
                        }
                        sizec2 += m_locMappers[c2][bI.second().patch].size();
                    }
                    sizec1 += m_locMappers[c1][bI.first().patch].size();
                }
            }
        }

        /*
        for(int p =0; p<m_comm.size();++p)
        {
            m_comm.barrier();
            if(p==m_comm.rank())
            {
                for(int i=0; i<locConnection.size();++i)
                {
                    const gsConnEntry &  connEntry = locConnection[i];
                    std::cout << "rank "<<gsMpi::worldRank()<<": "<<"("<<connEntry.glob<<","<<connEntry.proc<<","<<connEntry.loc<<")"<<"\n";
                }
            }
        }
        */

        connSize[m_comm.rank()]=locConnection.size();
        m_comm.sum(connSize.data(),connSize.size());
        size_t sum =0;
        std::vector<int> rDisp(m_comm.size(),0);
        for(index_t i=0; i<m_comm.size();++i)
        {
            sum+=connSize [i];
            if(i>0) rDisp[i] = rDisp[i-1]+connSize [i-1];
        }
        locConnection.reserve(sum);
//        int dummy[3];
        //locConnection.insert(locConnection.begin(),rDisp[m_comm.rank()],0); //insert enough empty space before the content
        locConnection.insert(locConnection.begin(),rDisp[m_comm.rank()],gsConnEntry()); //insert enough empty space before the content
        locConnection.resize(sum);//resize to the desired size

        //std::vector<int> sDisp(m_comm.size(),rDisp[m_comm.rank()]);
        //std::vector<int> sCount(m_comm.size(),connSize[m_comm.rank()]);

        MPI_Allgatherv(//locConnection.data(),sCount.data(),sDisp.data(),Base::m_connType,
                       MPI_IN_PLACE,        connSize[m_comm.rank()], Base::m_connType,
                locConnection.data(),connSize.data(),rDisp.data(),Base::m_connType,
                m_comm);

        for(size_t i=0; i<locConnection.size();++i)
        {
          const gsConnEntry& entry = locConnection[i];
          //  connMap[locConnection[i]].insert(procDof(locConnection[i+1],locConnection[i+2]));
          connMap[entry.glob].insert(procDof(entry.proc,entry.loc));
        }


        /*
        for(int p =0; p<m_comm.size();++p)
        {
            m_comm.barrier();
            if(p==m_comm.rank())
            {

                for(typename std::map<unsigned, std::set<procDof> >::const_iterator it = connMap.begin();
                    it != connMap.end(); ++it)
                {
                    std::cout << "rank "<<gsMpi::worldRank()<<": "<<it->first << " (";
                    for(typename std::set<procDof>::const_iterator sit = it->second.begin();sit != it->second.end();++sit)
                        gsInfo<<"("<<sit->first<<","<<sit->second<<"), ";
                    gsInfo<<")\n";
                }
            }
        }
        */
    }


protected:

    gsBoxTopology m_patches;
    const gsSortedVector<size_t>& m_myPatches;
    std::vector<gsMultiBasis<T> > m_bases;
    gsMatrix<bool> m_basesInteraction;
    std::vector<std::vector<gsDofMapper> > m_locMappers;
    const gsDofMapper & m_globMapper;
    const gsDofMapper& m_locOpPatchMappers;
    const gsMpiComm& m_comm;

    std::vector<int> m_patchToProc;


};


} // namespace gismo

#endif

#ifndef GISMO_BUILD_LIB
#include <gsIETI/gsParallelOperator.hpp>
#endif
