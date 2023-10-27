/** @file gsIETIAssemblerMPI.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2014-12-03
*/

#pragma once
#include <gsCore/gsConfig.h>


#ifdef GISMO_WITH_MPI

#include <gsUtils/gsStopwatch.h>

#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsBoxTopology.h>
#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETIUtils.h>

#include <gsMpi/gsMpi.h>
#include <gsMpi/gsMpiComm.h>

#include <gsUtils/gsSortedVector.h>



#include <functional>


namespace gismo {

template<typename T> class gsIETIJumpOperatorMPI;

struct gsIETIMPIOptions: gsIETIOptions
{
    gsIETIMPIOptions() : gsIETIOptions(), nSppHolder(1) {}
    int nSppHolder;
};

struct gsIETIInfoMPI
{
    unsigned lagrangeMultReduce;

    //Processors, not patches
    unsigned nNeigbours;
    gsSortedVector<size_t> procNeigbour;
    gsSortedVector<size_t> patchNeigbour;

    //indexed over patches
    //  std::vector<unsigned> startIdx;
    //  std::vector<unsigned> nLag2Neighb;
};

template<typename T>
class gsIETI_MPI_Base
{
protected:
    typedef typename gsIETIAssembler<T>::patchDof patchDof;
public:
    ~gsIETI_MPI_Base()
    {
        delete m_embeddingT_request;
        delete[] m_Spp_request;
        delete[] m_lagrangeRequest;
        delete m_Ibcast_request;
        if(m_options.nSppHolder!=m_comm.size() && m_options.nSppHolder!=1)
        {
            MPI_Comm_free(&m_embComm);
            if(isSppHolder())
                MPI_Comm_free(&m_embMasterComm);
        }
    }

    const gsIETIInfoMPI & getInfoMPI() const {return infoMPI;}

    bool isSppHolder() const {return m_isSppHolder;}

protected:
    gsIETI_MPI_Base(gsSortedVector<size_t> myPatches, MPI_Comm comm) :
        m_comm(comm),
        m_commEmbT(comm),
        m_commEmbTMaster(comm)
    {
        m_patchIdx = myPatches;
        m_Ibcast_request= NULL;
        m_Spp_request =NULL;
        m_lagrangeRequest=NULL;
    }

    virtual void setOptions(const gsIETIOptions & opt) {m_options = opt;}

    virtual void init_MPI(const gsIETIInfo& info, const std::vector<std::pair<patchDof,patchDof> >& lagrangeTable, std::vector<index_t>& freeElimLagr);

protected:
    void init_embT(const gsIETIInfo& info) const;

    void send_embT(const gsIETIInfo& info, size_t np, size_t nRhs) const;
    void send_embT(const gsIETIInfo& info, size_t nRhs) const;
    void finalize_embT(gsMatrix<T>& result) const;

    void init_Spp(const gsIETIInfo& info);
    void sendSPP() ;

    void finalizeSpp(const gsIETIInfo& info, const std::vector<std::vector<index_t> >& pDofsLoc2Glob, typename gsIETIAssembler<T>::sparseSPDfact*& LU_Spp);


public:
    void initEmbedding(const gsIETIInfo &info, gsMatrix<T> & xP)const;
    void initEmbedding_(const gsIETIInfo &info, gsMatrix<T> & xP, std::vector<gsMatrix<T> > &xPis)const;


    void accumulate(const gsMatrix<T>& input, gsMatrix<T>& result) const
    {
        //  gsInfo<<"rank: "<<m_comm.rank()<<" m_lagDataType: "<<m_lagDataType.size()<<"  lagSize: "<<m_lagrangeBuff.size()<<"\n";

        postAccumulate();

        // gsInfo<<"rank: "<<m_comm.rank()<<" posted wait: "<<"\n";

        startAccumulate(input);

        finishAccumulate(result);
    }

    void postAccumulate() const
    {
        for(size_t i = 0; i<infoMPI.nNeigbours;++i)
        {
            m_lagrangeBuff[i].setZero();
            int p = infoMPI.procNeigbour[i];
            MPI_Irecv(m_lagrangeBuff[i].data(),1,m_lagDataType[i],p,4,m_comm,&m_lagrangeRequest[i]);
        }
    }

    void startAccumulate(const  gsMatrix<T> & input) const
    {
        m_inputBuff = input;
        MPI_Request req;
        for(size_t i = 0; i<infoMPI.nNeigbours;++i)
        {
            int p = infoMPI.procNeigbour[i];
            MPI_Isend(m_inputBuff.data(),1,m_lagDataType[i],p,4,m_comm,&req);
            MPI_Request_free(&req);
        }
    }

    void finishAccumulate(gsMatrix<T> & result) const
    {
        result = m_inputBuff;
        int idx;
        for(size_t i = 0; i<infoMPI.nNeigbours;++i)
        {
            MPI_Waitany(infoMPI.nNeigbours,m_lagrangeRequest,&idx,MPI_STATUS_IGNORE);
            result+=m_lagrangeBuff[idx];
        }
    }

    void distribute(const gsMatrix<T>& input, gsMatrix<T>& output) const
    {
        output.setZero(input.rows(),input.cols());
        for(index_t i =0; i< input.rows();++i)
            output.row(i)=input.row(i)/m_mult[i];
    }

    void buildGlobalLagrangeMultiplier(const gsMatrix<T> & loc,gsMatrix<T>& glob) const
    {
        glob.setZero(m_originalLagrangeSize,loc.cols());
        for(size_t i=0; i<m_reducedLagrangeTable.size();++i)
            glob.row(m_lagLocToGlob[i])=loc.row(i);
        m_comm.sum(glob.data(),glob.rows()*glob.cols());

    }

    bool hasPatch(size_t np)const
    {
        return m_patch2proc[np] == m_comm.rank();
    }

    //Method to set the buffer size
    //For multible rhs the solver solves only for one rhs after the another -> resized to one column
    //For the construction of the solution the buffer must be resized again to its original size!
    void setBufferCols(index_t newCols) const
    {
           m_embT_buffer.setZero(m_embT_buffer.rows(),newCols);
           m_emb_buffer.setZero(m_emb_buffer.rows(),newCols);
           for(size_t i=0; i<m_lagrangeBuff.size();++i)
               m_lagrangeBuff[i].setZero(m_lagrangeBuff[i].rows(),newCols);
    }



public:
    std::vector<size_t> m_req2patch;
    std::vector<size_t> m_patch2req;
    gsSortedVector<size_t> m_patchIdx;
    std::vector<int> m_patch2proc;
    std::vector<gsSortedVector<size_t> > m_proc2patch;

    gsMpiComm m_comm;
    gsMpiComm m_commEmbT;
    gsMpiComm m_commEmbTMaster;

    gsIETIInfoMPI infoMPI;

protected:
    gsIETIOptions m_options;
    gsVector<int> m_nPDofs_group;
    gsVector<int> m_pDispl_group;
    gsVector<int> m_sppDispl;
    gsVector<int> m_nSppDofs;
    size_t m_oldPos;

    bool m_isSppHolder;
    int m_idxHolder;
    gsSortedVector<int> m_SppMaster;
    gsSortedVector<int> m_embGroup;

    std::vector<std::pair<patchDof,patchDof> > m_reducedLagrangeTable;
    std::vector<size_t> m_lagLocToGlob;
    std::vector<size_t> m_mult;
    size_t m_originalLagrangeSize;

  //  mutable std::vector<gsMatrix<T> > m_embT_buffer;
 //   mutable std::vector<gsMatrix<T> > m_Spp_buffer;

    mutable std::vector<gsMatrix<T> > m_lagrangeBuff;
    mutable gsMatrix<T> m_inputBuff;
    mutable gsMatrix<T> m_Spp_sendBuffer;
    mutable gsMatrix<T> m_Spp_recvBuffer;
    mutable gsMatrix<T> m_emb_buffer;
    mutable gsMatrix<T> m_embT_buffer2;
    mutable gsMatrix<T> m_embT_buffer;
    mutable bool m_finishedEmbT;

    mutable MPI_Request* m_Ibcast_request;
    mutable MPI_Request* m_embeddingT_request;
    MPI_Request* m_Spp_request;
    MPI_Request* m_lagrangeRequest;

    //one for each other processor (the one, which are no neighbors dont get any)
    //They work for both send and recv;
    std::vector<MPI_Datatype> m_lagDataType;

   // MPI_Comm m_graph_comm;
  //  MPI_Comm m_neigbours;
    MPI_Comm m_embComm;
    MPI_Comm m_embMasterComm;
  //  MPI_Group m_neighbourGoup;
  //  MPI_Group m_world_group;

};


/**
 * @brief gsIETIAssemblerMPI is the main assembling routine for the IETI method.
 * It assembles the stiffness matrices on each patch and calculates the LU factorization,
 * it also provides the bookkeeping for different splittings.

 */
template<typename T>
class gsIETIAssemblerMPI : public virtual gsIETIAssembler<T> , public virtual gsIETI_MPI_Base<T>
{
public:
    /// \brief  the method used for facorizing the matrices
    typedef typename gsIETIAssembler<T>::sparseLUfact sparseLUfact;
    typedef typename gsIETIAssembler<T>::sparseSPDfact sparseSPDfact;

    typedef typename gsIETIAssembler<T>::gsIETIAssembler Base;
    typedef typename gsIETI_MPI_Base<T>::gsIETI_MPI_Base MPI_Base;

    /// \brief  shortcut
    typedef typename gsIETIAssembler<T>::patchDof patchDof;

    friend class gsIETIJumpOperatorMPI<T>;

    typedef memory::shared_ptr<gsIETIAssemblerMPI> Ptr;
    typedef memory::unique_ptr<gsIETIAssemblerMPI> uPtr;


protected:
    /// \brief  shortcut for assembling matrices
    typedef Eigen::Triplet<T> Trip;

    virtual void printTime(gsStopwatch& time, std::string text)
    {
        if(m_comm.rank()==0)
        {
#pragma omp master
            {
                time.stop();
                gsInfo<<text<<time<<std::endl<<std::flush;
                time.restart();
            }
        }
    }

    virtual void printWarn(std::string text)
    {
        if(m_comm.rank()==0)
        {
#pragma omp master
            {
                gsWarn<<text<<std::endl;
            }
        }
    }

    virtual void print(std::string text)
    {
        if(m_comm.rank()==0)
        {
#pragma omp master
            {
                gsInfo<<text<<std::endl<<std::flush;
            }
        }
    }

    virtual void storeSpp(const gsMatrix<T>& in, gsMatrix<T>& out);

    virtual void extractPatch(size_t np, const gsMatrix<T>& rhs, gsMatrix<T>& rhsLocal) const;
public:

    /**
     * @brief gsIETIAssemblerMPI The constructor of the IETI assembler
     * @param patches the multi patch representation
     * @param bases the basis corresponding to the patch
     * @param bconditions the boundary conditions
     * @param rhs the right hand side as gsFunction<T>
     * @param isMinimalEnergy should the minimum energy method be used (default is true)
     * @param strat_ the strategy for selecting the primal variables (default is primalDofMethod::all)
     * @param dirStrategy strategy for incormporating the dirichlet boundary (default is dirichlet::elimination)
     */
    gsIETIAssemblerMPI(gsAssembler<T> &assembler,gsSortedVector<size_t> myPatches = gsSortedVector<size_t>(), MPI_Comm comm = gsMpi::worldComm());

    ~gsIETIAssemblerMPI()
    {

    }

    static gsOptionList defaultOptions();
    virtual void setOptions(const gsOptionList & opt);              ///< Set the options based on a gsOptionList

    virtual void init();

    virtual void setNewRhs(const gsMatrix<T>& rhs);
    virtual void giveAssembledMatrices(std::vector<gsSparseMatrix<T> >& matrices, const gsMatrix<T> &rhs);

    /**
     * @brief assemble assembles all required matrices and LU factorizations and the corresponding right hand sides.
     *
     * This functions calls:
     *      assembleKiiKibKbb(tripletList);
     *      assembleSpp();
     *      assembleRhs(rhs_loc);
     *
     * In case of Minimum energy also:
     *      assembleC();
     *      assembleLUofKC(tripletList);
     *
     * and otherwise:
     *      assembleKrrKrpKpp(tripletList);
     *
     */
    virtual void assemble(const gsMultiPatch<T>& curSol=gsMultiPatch<T>());

    /**
     * @brief processSolution takes the results of IETISolver and returns from the solution (xP,(x2)_i) the
     * solution on the whole domain as a gsMatrix<T>. The output is the same as calling the solver without IETI.
     * IETISolver::calulateSolution needs to be called before. In order to receive a gsField on the whole domain,
     * use the assembler specific constructSolution routine.
     *
     * @param uP primal component of the solution
     * @param u2 ICO minimalEnergy u2 is the solution on the interface and else the solution on the patch.
     *   (if isMinimumEnergy: dofsB[np] x nRhs ; else dofsR[np] x nRhs)
     * @return solVec the solution of the system in matrix form.
     */
    virtual void processSolution(const gsMatrix<T>& uP_, const std::vector<gsMatrix<T> >& u2, gsMatrix<T>& solVec) const;

    /**
     * @brief embedding transforms a representation of the solution vector from a primal and remaining to patch wise values on the interface
     * @param xP the global primal component ( totalPdofs x nRhs)
     * @param x2 the local remaining/dual component for each patch (if isMinimumEnergy: dofsB[np] x nRhs ; else dofsR[np] x nRhs)
     * @param w the patch wise values on the interface  (dofsB[np] x nRhs)
     */
    virtual void initEmbedding(const gsIETIInfo& info, gsMatrix<T>& xP, std::vector<gsMatrix<T> > & xPis)const;
    virtual void embedding(std::vector<gsMatrix<T> >  & xPi,std::vector<gsMatrix<T> >const & x2,std::vector<gsMatrix<T> >& w) const;
    virtual void embedding(gsMatrix<T>  & xPi, gsMatrix<T> const & x2, size_t np, gsMatrix<T>& w) const;
    virtual void embeddingTransPartial(const gsMatrix<T> & w,size_t np, gsMatrix<T> &  u2) const;

    bool hasFinalizedEmbT() const
    {
        return MPI_Base::m_finishedEmbT;
    }

    virtual void finalizeEmbeddingTrans(gsMatrix<T>& result,bool saveToRhsP=false) const
    {
        MPI_Base::finalize_embT(result);
        if(saveToRhsP && MPI_Base::isSppHolder())
        {
            const gsMatrix<T> & rhs_p= m_rhs_p;
            const_cast<gsMatrix<T>& >(rhs_p) = result;
        }
    }


    virtual void embeddingTrans(const std::vector<gsMatrix<T> >& w,gsMatrix<T>&  uP, std::vector<gsMatrix<T> > &  u2) const;

    virtual void assembleRhs(const std::vector<gsMatrix<T> >& rhs_loc, size_t np);

    virtual void assembleRhsFreeElim();

    virtual void combineToCommonSolution(gsMatrix<T>& solVec) const;

  //  const gsIETIMPIOptions & getOptionsMPI() const {return m_options;}
    /*
    void accumulateCollective(const gsMatrix<T>& input, gsMatrix<T>& result) const
    {
      //  gsInfo<<"rank: "<<m_comm.rank()<<" m_lagDataType: "<<m_lagDataType.size()<<"  lagSize: "<<m_lagrangeBuff.size()<<"\n";
        for(size_t i = 0; i<infoMPI.procNeigbour.size();++i)
        {
            int p = infoMPI.procNeigbourHood[i];
           gsInfo<<"rank: "<<m_comm.rank()<<" lag buff: "<<m_lagrangeBuff[i].rows()<<" group rank: "<<p<<"\n";
           MPI_Irecv(m_lagrangeBuff[i].data(),1,m_lagDataType[i],p,4,m_neigbours,&m_lagrangeRequest[i]);
        }


        gsInfo<<"rank: "<<m_comm.rank()<<" posted! "<<"\n";

        m_inputBuff = input;

        MPI_Request req;
        for(size_t i = 0; i<infoMPI.procNeigbour.size();++i)
        {
           int p = infoMPI.procNeigbourHood[i];
           MPI_Isend(m_inputBuff.data(),1,m_lagDataType[i],p,4,m_neigbours,&req);
           MPI_Request_free(&req);
        }

        gsInfo<<"rank: "<<m_comm.rank()<<" sent! "<<"\n";

        result = input;
        int idx;
        for(size_t i = 0; i<infoMPI.procNeigbour.size();++i)
        {
            gsInfo<<"Waiting\n";
           MPI_Waitany(infoMPI.nNeigbours,m_lagrangeRequest,&idx,MPI_STATUS_IGNORE);
           gsInfo<<"Waiting, got "<<infoMPI.procNeigbour[idx]<<"\n";
           result+=m_lagrangeBuff[idx];
        }

    }
*/


protected:
    /*
    struct lagrangeComp
    {
        bool operator() (std::pair<patchDof,patchDof> a,std::pair<patchDof,patchDof> b)
        {
            size_t pa = std::min(a.first.first,a.second.first);
            size_t pb = std::min(b.first.first,b.second.first);

            if(pa < pb)
              return true;
            else if(pa>pb)
              return false;
            else //pa == pb
            {
                size_t paa = std::max(a.first.first,a.second.first);
                size_t pbb = std::max(b.first.first,b.second.first);

                if(paa < pbb)
                  return true;
                else if(paa>pbb)
                  return false;
                else //paa == pbb
                {
                    //pa==pb

                    unsigned da = pa==a.first.first? a.first.second : a.second.second;
                    unsigned db = pb==b.first.first? b.first.second : b.second.second;

                    GISMO_ASSERT(da!=db, "the same Pair is compared");

                    if(da< db)
                        return true;
                    else if(da>db)
                        return false;
                }
            }
        }

    };
*/
public:


    gsIETIAssemblerMPI() : MPI_Base(gsSortedVector<size_t>(),gsMpi::worldComm()) {}

protected:
    /// \brief  An info Structure
    using gsIETIAssembler<T>::info;


    /// @brief pointer to the used assembler
    using gsIETIAssembler<T>::m_assembler;

    /// @brief The multipatch domain
    using gsIETIAssembler<T>::m_patches;

    /// \brief   The discretization bases corresponding to \a m_patches and to
    /// the number of solution fields that are to be computed
    /// m_basis[i]: The multi-basis for unknown i
    using gsIETIAssembler<T>::m_basis;

    /// \brief backup of the "full" boundary conditions
    using gsIETIAssembler<T>::m_bConditions;

    /// \brief  contains global information of the primal dofs to the system (for each component)
    using gsIETIAssembler<T>::m_primalDofMapper;

    /// \brief  global standart mapper represent the matched patches (for each component)
    using gsIETIAssembler<T>::m_stdMapper;

    /// \brief  local dof mapper, map a basis function of a patch to the real (where elimination happend) index (locally) and for each component .
    /// m_locDofsMapper[np] -> all mappers on patch np.
    using gsIETIAssembler<T>::m_locDofsMapper;

    /// \brief  contains information of coupled edges in 3D
    using gsIETIAssembler<T>::m_edgeMapper;

    ///The stored primal dof mappings dont have a dirichlet boundary, so you dont have to care.
    //(these maps works for uPi -> uP)

    /// \brief  local number of patch k to the global number (in uP).  length: (nPatches)
    using gsIETIAssembler<T>::m_pDofsLoc2Glob;

    //for all stored dofs it is assumed that the dirichlet boundary is not eliminated!!!!
    //use mapper.index(i) to get the real index in the matrix.

    /// \brief  contains the local numbers of Vertex Primal dofs for each patch
    using gsIETIAssembler<T>::m_primalVdofs;

    /// \brief  Contains the Edges for each patch in 3D
    using gsIETIAssembler<T>::m_primalEdges;

    /// \brief  contains the local number of the remaining dofs for each patch
    using gsIETIAssembler<T>::m_remDofs;

    /// \brief  contains the local numbers of interface dofs for each patch
    using gsIETIAssembler<T>::m_boundDofs;

    /// \brief  on a patch k contains for a remaining dof i true if it is also on the interface, i.e. is in m_boundDofs[k]
    using gsIETIAssembler<T>::m_remDofIsAlsoBoundDof;

    /// \brief  contains for each patch the the interface sides
    using gsIETIAssembler<T>::m_patchISides;

    /// \brief  boundary conditions for each patch
    using gsIETIAssembler<T>::m_bc_loc;

    /// \brief  maps the global (patchwise) dof to the local interface index or interior index
    using gsIETIAssembler<T>::m_glob2BoundInteriorIndex;

    /// \brief  contains true if a global(patchwise) index is a interface index
    using gsIETIAssembler<T>::m_globIsBoundIndex;

    /// \brief  maps a global index to the repective primal or remaining index (ONLY FOR noMinimalEnergy)
    using gsIETIAssembler<T>::m_glob2PrimalRemainingIndex;

    /// \brief  contains true if a global index is a primal one (ONLY FOR noMinimalEnergy)
    using gsIETIAssembler<T>::m_globIsPrimalIndex;

    /// \brief  lagrange mult -> the pair of patchDofs
    using gsIETIAssembler<T>::m_lagrangeTable;

    /// \brief  global interface dof -> merged nodes+patches
    using gsIETIAssembler<T>::m_globalConnectionTable;

    /// \brief  lagrange multipl which couple a eliminated and free dof
    using gsIETIAssembler<T>::m_freeElimLagr;
    //------------------------------------------------------

    /// \brief  The S-orth Local Basis of the Primal variable Space
    using gsIETIAssembler<T>::m_Phi;

    /// \brief  Local Lu factorization of [K C^T,C 0]
    using gsIETIAssembler<T>::m_LU_KC;

    /// \brief  Local Matrix representation of C
    using gsIETIAssembler<T>::m_C;

    /// \brief  Local LU factorization of Krr
    using gsIETIAssembler<T>::m_LU_Krr;

    /// \brief  Global Schur complement wrt. the primal Variables
    using gsIETIAssembler<T>::m_LU_Spp;

    /// \brief  Local Matrix Representation of Krp
    using gsIETIAssembler<T>::m_Krp;

    /// \brief  Local representations of Kpp
    using gsIETIAssembler<T>::m_Kpp;

    // Vector rhs
    /// \brief  Primal components of the rhs
    using gsIETIAssembler<T>::m_rhs_p;

    /// \brief  Remaining or dual components of the rhs
    using gsIETIAssembler<T>::m_rhs_d;

    /// \brief  Interior components of the rhs
    using gsIETIAssembler<T>::m_rhs_i;

    /// \brief  dirichlet dofs  (if available, e.g. if using nitsche)
    using gsIETIAssembler<T>::m_dirDofs;

    /// \brief  rhs vector containing the dirichlet values for coupled  eliminated and free dofs (for gsIETIJumpOperator)
    using gsIETIAssembler<T>::m_rhs_dir;

    /// \brief  Local LU factorization of Kii
    using gsIETIAssembler<T>::m_LU_Kii;

    /// \brief  Local matrix representation of Kib
    using gsIETIAssembler<T>::m_Kib;
    using gsIETIAssembler<T>::m_Kbi;

    /// \brief  Local matrix representation of Kbb
    using gsIETIAssembler<T>::m_Kbb;

    using gsIETIAssembler<T>::m_IETIoptions; //gsIETIMPIOptions m_options;

protected:
    /* auxiliary variables */
    using gsIETIAssembler<T>::m_scalingsKC;
    using gsIETIAssembler<T>::m_scalingsKii;
    using gsIETIAssembler<T>::m_scalingsKrr;
    using gsIETIAssembler<T>::orient;

public: //stuff from MPI
    using gsIETI_MPI_Base<T>::infoMPI;
    using gsIETI_MPI_Base<T>::m_comm;
    using gsIETI_MPI_Base<T>::m_req2patch;
    using gsIETI_MPI_Base<T>::m_patch2req;
    using gsIETI_MPI_Base<T>::m_patchIdx;
    using gsIETI_MPI_Base<T>::m_patch2proc;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIETIAssemblerMPI.hpp)
#endif

#endif
