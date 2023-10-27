/** @file gsIETIdGAssemblerMPI.h

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


//#include <map>

#include <gsAssembler/gsVisitorDg2.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsBoxTopology.h>
#include <gsUtils/gsSortedVector.h>

#include <gsIETI/gsIETIAssemblerMPI.h>
#include <gsIETI/gsIETIdGAssembler.h>
#include <gsIETI/gsIETIUtils.h>




namespace gismo {

template<typename T> class gsIETIJumpOperatorMPI;

template<typename T>
class gsIETIdG_MPI_Base : public virtual gsIETI_MPI_Base<T>
{
    typedef gsIETI_MPI_Base<T> Base;
    typedef typename gsIETI_MPI_Base<T>::patchDof patchDof;

public:
    ~gsIETIdG_MPI_Base()
    {
        delete[] m_ddof_request;
    }

protected:
    gsIETIdG_MPI_Base(gsSortedVector<size_t> myPatches, MPI_Comm comm) : gsIETI_MPI_Base<T>(myPatches, comm) {}

    void init_MPI(const gsIETIInfo& info,
                  const std::vector<std::pair<patchDof,patchDof> >& lagrangeTable,
                  std::vector<index_t>& freeElimLagr,
                  const std::vector<std::pair<patchDof,patchDof> >& elimExtraBasisConnection);

    void init_ddof();

    void send_ddof();

    bool finishOneDdof(int &p, gsMatrix<T>& ddof);

public:
    using gsIETI_MPI_Base<T>::m_req2patch;
    using gsIETI_MPI_Base<T>::m_patch2proc;
    using gsIETI_MPI_Base<T>::m_patch2req;
    using gsIETI_MPI_Base<T>::m_patchIdx;
    using gsIETI_MPI_Base<T>::m_comm;
    using gsIETI_MPI_Base<T>::infoMPI;

  //  using gsIETI_MPI_Base<T>::m_graph_comm;
  //  using gsIETI_MPI_Base<T>::m_neigbours;
  //  using gsIETI_MPI_Base<T>::m_neighbourGoup;
  //  using gsIETI_MPI_Base<T>::m_world_group;

protected:
    std::vector<std::pair<patchDof,patchDof> > m_ddof_connection;
    std::vector<std::vector<T> > m_ddof_shiftS;
    std::vector<std::vector<T> > m_ddof_shiftR;
protected:
    mutable std::vector<gsMatrix<T> > m_ddof_Sbuffer;
    mutable std::vector<gsMatrix<T> > m_ddof_Rbuffer;
    MPI_Request* m_ddof_request;
    int m_ddof_count;
};


/**
 * @brief gsIETIdGAssembler is the main assembling routine for the IETI method.
 * It assembles the stiffness matrices on each patch and calculates the LU factorization,
 * it also provides the bookkeeping for different splittings.
 *
 *The code is NOT working for domains, which consist of no interior basis functions. in such a case,
 *just do one refinement.

 */
template<typename T>
class gsIETIdGAssemblerMPI : public gsIETIdGAssembler<T>, public gsIETIAssemblerMPI<T>, public gsIETIdG_MPI_Base<T>
{
public:
    /// \brief  the method used for facorizing the matrices
    typedef typename gsIETIdGAssembler<T>::sparseLUfact sparseLUfact;
    typedef typename gsIETIdGAssembler<T>::sparseSPDfact sparseSPDfact;

    /// \brief  shortcut
    typedef typename gsIETIdGAssembler<T>::patchDof patchDof;

    friend class gsIETIJumpOperatorMPI<T>;
protected:
    /// \brief  shortcut for assembling matrices
    typedef Eigen::Triplet<real_t> Trip;

    typedef gsIETIAssemblerMPI<T> Base;
    typedef gsIETIdGAssembler<T> dG_Base;
    typedef typename gsIETIAssemblerMPI<T>::MPI_Base MPI_Base;
    typedef gsIETIdG_MPI_Base<T> MPIdG_Base;

    /**
     * @brief compCalc calculates the index in a consecutive number of the
     *  dofs of a given index \a idx on patch \a np for a component \a c.
     *  e.g. the second dof of the third component on patch 1 results in 300th dof.
     *  (np = 0, idx = 1, c = 2 -> return 299)
     * @param np the patch number
     * @param idx the index of the dof of the given component
     * @param c the component
     * @return  the consecutive index
     */
    inline unsigned compCalc(int np,unsigned idx, int c) const
    {
        return c*(m_basis.front()[np].size()+m_numExtraBasisTotal[np]) +idx;
    }
    /**
     * @brief compCalcBack given a consecutive dof number on a patch, the function
     * calculates the index w.r.t. to its component. In order to get the component
     * use getComp(..).
     * e.g. the 300th dof on patch 1 has index 2
     * (np =0, idx  = 299 -> return 1)
     * @param np the patch number
     * @param idx the consecutive index of the dof
     * @return  the index of the dof w.r.t. its component
     */
    inline unsigned compCalcBack(int np,unsigned idx) const
    {
        return idx%(m_basis.front()[np].size()+m_numExtraBasisTotal[np]);
    }
    /**
     * @brief getComp given a consecutive dof number on a patch, the function calculates
     * its component. In order to get the index use compCalcBack(..).
     * e.g. the 300th dof on patch 1 has component 3
     * (np = 0, idy = 299 -> return 2)
     * @param np the patch number
     * @param idx the consecutive index of the dof
     * @return  the component of the dof
     */
    inline int getComp(int np,unsigned idx) const
    {
        return idx/(m_basis.front()[np].size()+m_numExtraBasisTotal[np]);
    }


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



    virtual void assembleDgInterfaceContribution(std::vector<gsSparseMatrix<T> >& matrices,std::vector<gsMatrix<T> >& rhs) const;

    /*
    void prepareActives(const boundaryInterface & bi, gsMatrix<unsigned>& actives1,gsMatrix<unsigned>& actives2,
                        gsMatrix<unsigned>& activesExtra1,gsMatrix<unsigned>& activesExtra2) const;
*/
    virtual void combineToCommonSolution(gsMatrix<T>& solVec) const;

    void updateDDofs();

public:

    /**
     * @brief gsIETIdGAssembler The constructor of the IETI assembler
     * @param patches the multi patch representation
     * @param bases the basis corresponding to the patch
     * @param bconditions the boundary conditions
     * @param rhs the right hand side as gsFunction<T>
     * @param isMinimalEnergy should the minimum energy method be used (default is true)
     * @param strat_ the strategy for selecting the primal variables (default is primalDofMethod::all)
     * @param dirStrategy strategy for incormporating the dirichlet boundary (default is dirichlet::elimination)
     */
    gsIETIdGAssemblerMPI(gsAssembler<T> &assembler,gsSortedVector<size_t> myPatches = gsSortedVector<size_t>(), MPI_Comm comm = gsMpi::worldComm());

    static gsOptionList defaultOptions();
    virtual void setOptions(const gsOptionList & opt);              ///< Set the options based on a gsOptionList

    virtual void init();

    virtual void prepareDDofs();

    virtual void assemble(const gsMultiPatch<T>& curSol=gsMultiPatch<T>());

private:
    /*
    gsIETIdGAssemblerMPI() :
        gsIETIAssembler<T>(),
        gsIETI_MPI_Base<T>(),
        gsIETIdGAssembler<T>(),
        gsIETIAssemblerMPI<T>(),
        gsIETIdG_MPI_Base<T>()
    {}
    */

protected:

    /// @brief pointer to the used assembler
    using gsIETIAssemblerMPI<T>::m_assembler;

    /// @brief The multipatch domain
    using gsIETIAssemblerMPI<T>::m_patches;

    /// \brief   The discretization bases corresponding to \a m_patches and to
    /// the number of solution fields that are to be computed
    /// m_basis[i]: The multi-basis for unknown i
    using gsIETIAssemblerMPI<T>::m_basis;

    /// \brief backup of the "full" boundary conditions
    using gsIETIAssemblerMPI<T>::m_bConditions;

    /// \brief  An info Structure
    using gsIETIAssemblerMPI<T>::info;

    /// \brief  contains global information of the primal dofs to the system (for each component)
    using gsIETIAssemblerMPI<T>::m_primalDofMapper;

    /// \brief  global standart mapper represent the matched patches (for each component)
    using gsIETIAssemblerMPI<T>::m_stdMapper;

    /// \brief  local dof mapper, map a basis function of a patch to the real (where elimination happend) index (locally) and for each component .
    /// m_locDofsMapper[np] -> all mappers on patch np.
    using gsIETIAssemblerMPI<T>::m_locDofsMapper;

    /// \brief  contains information of coupled edges in 3D
    using gsIETIAssemblerMPI<T>::m_edgeMapper;

    ///The stored primal dof mappings dont have a dirichlet boundary, so you dont have to care.
    //(these maps works for uPi -> uP)

    /// \brief  local number of patch k to the global number (in uP).  length: (nPatches)
    using gsIETIAssemblerMPI<T>::m_pDofsLoc2Glob;

    //for all stored dofs it is assumed that the dirichlet boundary is not eliminated!!!!
    //use mapper.index(i) to get the real index in the matrix.

    /// \brief  contains the local numbers of Vertex Primal dofs for each patch
    using gsIETIAssemblerMPI<T>::m_primalVdofs;

    /// \brief  Contains the Edges for each patch in 3D
    using gsIETIAssemblerMPI<T>::m_primalEdges;

    /// \brief  contains the local number of the remaining dofs for each patch
    using gsIETIAssemblerMPI<T>::m_remDofs;

    /// \brief  contains the local numbers of interface dofs for each patch
    using gsIETIAssemblerMPI<T>::m_boundDofs;

    /// \brief  on a patch k contains for a remaining dof i true if it is also on the interface, i.e. is in m_boundDofs[k]
    using gsIETIAssemblerMPI<T>::m_remDofIsAlsoBoundDof;

    /// \brief  contains for each patch the the interface sides
    using gsIETIAssemblerMPI<T>::m_patchISides;

    /// \brief  boundary conditions for each patch
    using gsIETIAssemblerMPI<T>::m_bc_loc;

    /// \brief  maps the global (patchwise) dof to the local interface index or interior index
    using gsIETIAssemblerMPI<T>::m_glob2BoundInteriorIndex;

    /// \brief  contains true if a global(patchwise) index is a interface index
    using gsIETIAssemblerMPI<T>::m_globIsBoundIndex;

    /// \brief  maps a global index to the repective primal or remaining index (ONLY FOR noMinimalEnergy)
    using gsIETIAssemblerMPI<T>::m_glob2PrimalRemainingIndex;

    /// \brief  contains true if a global index is a primal one (ONLY FOR noMinimalEnergy)
    using gsIETIAssemblerMPI<T>::m_globIsPrimalIndex;

    /// \brief  lagrange mult -> the pair of patchDofs
    using gsIETIAssemblerMPI<T>::m_lagrangeTable;

    /// \brief  global interface dof -> merged nodes+patches
    using gsIETIAssemblerMPI<T>::m_globalConnectionTable;

    /// \brief  lagrange multipl which couple a eliminated and free dof
    using gsIETIAssemblerMPI<T>::m_freeElimLagr;
    //------------------------------------------------------

    /// \brief  The S-orth Local Basis of the Primal variable Space
    using gsIETIAssemblerMPI<T>::m_Phi;

    /// \brief  Local Lu factorization of [K C^T,C 0]
    using gsIETIAssemblerMPI<T>::m_LU_KC;

    /// \brief  Local Matrix representation of C
    using gsIETIAssemblerMPI<T>::m_C;

    /// \brief  Local LU factorization of Krr
    using gsIETIAssemblerMPI<T>::m_LU_Krr;

    /// \brief  Global Schur complement wrt. the primal Variables
    using gsIETIAssemblerMPI<T>::m_LU_Spp;

    /// \brief  Local Matrix Representation of Krp
    using gsIETIAssemblerMPI<T>::m_Krp;

    /// \brief  Local representations of Kpp
    using gsIETIAssemblerMPI<T>::m_Kpp;

    // Vector rhs
    /// \brief  Primal components of the rhs
    using gsIETIAssemblerMPI<T>::m_rhs_p;

    /// \brief  Remaining or dual components of the rhs
    using gsIETIAssemblerMPI<T>::m_rhs_d;

    /// \brief  Interior components of the rhs
    using gsIETIAssemblerMPI<T>::m_rhs_i;

    /// \brief  dirichlet dofs  (if available, e.g. if using nitsche)
    using gsIETIAssemblerMPI<T>::m_dirDofs;

    /// \brief  rhs vector containing the dirichlet values for coupled  eliminated and free dofs (for gsIETIJumpOperatorMPI)
    using gsIETIAssemblerMPI<T>::m_rhs_dir;

    /// \brief  Local LU factorization of Kii
    using gsIETIAssemblerMPI<T>::m_LU_Kii;

    /// \brief  Local matrix representation of Kib
    using gsIETIAssemblerMPI<T>::m_Kib;

    /// \brief  Local matrix representation of Kbb
    using gsIETIAssemblerMPI<T>::m_Kbb;

    using gsIETIAssemblerMPI<T>::m_IETIoptions;
    //----------------------------------------------------------
    // new Stuff
    //DG - specific stuff
    //ordering of the extra dofs:
    // [0,....,n-1]<- usual dofs on the patch
    // [n,...,ntot]<- extra dofs [n,...,n0-1,n0,...,n1-1,..., ntot] , where ni are the number of the neigboring patch corresponding to side.index()==i

    /// \brief contains the number of extra basis functions for a specific side of a patch
    /// m_numExtraBasis[patch][side.index()]. contains [4,5,.,7], i.e. sides, which dont have an interface have undefined value.
    using gsIETIdGAssembler<T>::m_numExtraBasis;

    using gsIETIdGAssembler<T>::m_numExtraBasisTotal;

    using gsIETIdGAssembler<T>::m_boundaryBasis;

    using gsIETIdGAssembler<T>::m_elimExtraBasisConnection;

protected:
    /* auxiliary variables */
    using gsIETIAssemblerMPI<T>::m_scalingsKC;
    using gsIETIAssemblerMPI<T>::m_scalingsKii;
    using gsIETIAssemblerMPI<T>::m_scalingsKrr;

    using gsIETIAssemblerMPI<T>::orient;

public: //stuff from MPI
    using gsIETIAssemblerMPI<T>::infoMPI;
    using gsIETIAssemblerMPI<T>::m_comm;
    using gsIETIAssemblerMPI<T>::m_req2patch;
    using gsIETIAssemblerMPI<T>::m_patch2req;
    using gsIETIAssemblerMPI<T>::m_patchIdx;
    using gsIETIAssemblerMPI<T>::m_patch2proc;

};



} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIETIdGAssemblerMPI.hpp)
#endif
//#include "gsIETIdGAssemblerMPI.hpp"

#endif
