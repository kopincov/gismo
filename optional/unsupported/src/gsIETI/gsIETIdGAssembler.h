/** @file gsIETIdGAssembler.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2014-12-03
*/



#pragma once

//#include <map>

#include <gsAssembler/gsVisitorDg2.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsBoxTopology.h>
#include <gsUtils/gsSortedVector.h>

#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETIUtils.h>


namespace gismo {

template<typename T> class gsIETIJumpOperator;


/**
 * @brief gsIETIAssembler is the main assembling routine for the IETI method.
 * It assembles the stiffness matrices on each patch and calculates the LU factorization,
 * it also provides the bookkeeping for different splittings.
 *
 *The code is NOT working for domain, which consists of no interior basis functions. in such a case,
 *just do one refinement.

 */
template<typename T>
class gsIETIdGAssembler : public virtual gsIETIAssembler<T>
{
public:
    /// \brief  the method used for facorizing the matrices
    typedef typename gsIETIAssembler<T>::sparseLUfact sparseLUfact;
    typedef typename gsIETIAssembler<T>::sparseSPDfact sparseSPDfact;

    /// \brief  shortcut
    typedef typename gsIETIAssembler<T>::patchDof patchDof;

    friend class gsIETIJumpOperator<T>;
protected:
    /// \brief  shortcut for assembling matrices
    typedef Eigen::Triplet<real_t> Trip;

    typedef gsIETIAssembler<T> Base;


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

    /**
     * @brief setInfo Populates the entries of the gsIETIInfo structure. Needs to be called after
         init_patchInterfaceSides();
         init_primalDofs();
         init_boundaryDofs();
         init_jumpOperatorData();

        and called before
         init_remainingDofs();

     */
    void setInfo();


    /**
     * @brief init_mappers initializes the mappers for the bookkeeping. Needs to be called before all other init methods.
     */
    void init_mappers();

    /**
     * @brief init_jumpOperatorData initializes the data for the jump operator B,
     * i.e. gsIETIJumpOperator. Needs to be called before setInfo(..);
     */
    void init_jumpOperatorData();

    /**
     * @brief init_patchInterfaceSides Initializes the sides for each patch lieing on an interface.
     * Needs to be called before setInfo(..);
     */
    void init_patchInterfaceSides();

    /**
     * @brief init_primalDofs Initializes the info for the primal variables, needs to be called after
     *       init_patchInterfaceSides();
     *
     *  and called before
     *      init_boundaryDofs(isMinimalEnergy);
     *      setInfo(bool isMinimalEnergy);
     *      init_remainingDofs();
     *
     *
     * Global Primal Dof are ordered in the following way:
     *
     *[c1,c2,c3,....,cN,f1,...,fM,e1,e2,...,eK]
     *
     * where ci are vertex values, fi are face averages and ei are edge averages
     *
     * in 2D there are no faces, just edges, but they are computed like the faces. (just interfaces)
     * (an edge is always a 2D object and a face always a 3D object, independent of the object dimension)
     *
     **/
    void init_primalDofs();

    /**
     * @brief init_remainingDofs Initializes the splitting for primal and remaining dofs. This is only needed in case of
     *  minimumEnergy is false.
     *  Has to be called after all init_* functions and setInfo.
     */
    void init_remainingDofs();

    /**
     * @brief init_boundaryDofs Initializes the splitting in boundary and remaining dofs. Has to be called after
     *    init_primalDofs()
     *  and called before
     *    setInfo(..)
     *    init_boundaryDofs();
     */
    void init_boundaryDofs();



    /**
     * @brief assembleC assembes the constraint matrix C for the minimal Energy method, which is used to enforce the primal constraints
     */
    void assembleC(size_t np);

    //DG- related functions

    int dgOffset(const int patch,const boxSide &side) const
    {
        //side =1 -> offset =0; (side ==0 -> internal dof!)
        int offset =0;
        for(int i=1; i< side.index();i++)
            offset+=m_numExtraBasis[patch][i];
        return offset;
    }


    index_t dgFindSideIndex(const patchSide& side,const index_t& localBasisNumber) const
    {
        //find the index of the basis function on the boundary w.r.t to the side

        for(int k=0; k<m_boundaryBasis[side.patch][side.index()].rows();k++)
            if((m_boundaryBasis[side.patch][side.index()])(k,0)==localBasisNumber)
                return k;
        return -1;

    }

    index_t dgFindCorrespondingExtraIndex(const patchSide& side,const patchSide& side2, const index_t& localBasisNumber2) const
    {
        unsigned n= m_basis.front().size(side.patch);
        int k= dgFindSideIndex(side2, localBasisNumber2); //optimize this!!!!!!! probably O(logN)
        if(k==-1)
            return -1;
        else
            return n+dgOffset(side.patch,side)+k;
    }

    int getIndexOfExtraSide(size_t np, const patchSide& side) const
    {
        std::map<patchSide,int>::const_iterator it = m_patchSideToExtraIndex[np].find(side);
        return it!=m_patchSideToExtraIndex[np].end()? it->second : -1;
    }

    virtual void assembleDgInterfaceContribution(std::vector<gsSparseMatrix<T> >& matrices,std::vector<gsMatrix<T> >& rhs) const;

    void prepareActives(const boundaryInterface & bi, gsMatrix<index_t>& actives1,gsMatrix<index_t>& actives2,
                        gsMatrix<index_t>& activesExtra1,gsMatrix<index_t>& activesExtra2) const;

    virtual void updateDDofs();

    virtual void setupMGforKC(size_t np);
    virtual void getSchur(gsMatrix<T> &schur, size_t np);

    virtual void setInterfaceConditions(std::vector<iFace::strategy>& iFaceStrategy) {m_iCoupling = iFaceStrategy;}
public:
    gsIETIdGAssembler() {}

public:

    /**
     * @brief gsIETIAssembler The constructor of the IETI assembler
     * @param patches the multi patch representation
     * @param bases the basis corresponding to the patch
     * @param bconditions the boundary conditions
     * @param rhs the right hand side as gsFunction<T>
     * @param isMinimalEnergy should the minimum energy method be used (default is true)
     * @param strat_ the strategy for selecting the primal variables (default is primalDofMethod::all)
     * @param dirStrategy strategy for incormporating the dirichlet boundary (default is dirichlet::elimination)
     */
    gsIETIdGAssembler(gsAssembler<T> &assembler);

    static gsOptionList defaultOptions();
    virtual void setOptions(const gsOptionList & opt);              ///< Set the options based on a gsOptionList

    /// \brief  destructor
    //~gsIETIdGAssembler()
    //{
    //    for(size_t np = 0; np<info.numberPatches;np++)
    //        freeAll(m_boundaryBasis[np]);
    //}

    /// \brief initializes the IETI-method
      virtual void init();

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
    void assemble(const gsMultiPatch<T>& curSol= gsMultiPatch<T>());

protected:

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

    /// \brief  An info Structure
    using gsIETIAssembler<T>::info;

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

    /// \brief  contains for each patch the the interface sides
    using gsIETIAssembler<T>::m_patchAveragesSides;

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

    /// \brief  Local matrix representation of Kbb
    using gsIETIAssembler<T>::m_Kbb;

    using gsIETIAssembler<T>::m_IETIoptions;

    // members for inexact version
    using gsIETIAssembler<T>::m_exC;
    using gsIETIAssembler<T>::m_timings;
    using gsIETIAssembler<T>::m_K;
    using gsIETIAssembler<T>::m_S;
    using gsIETIAssembler<T>::m_inexact_KC;
    //using gsIETIAssembler<T>::m_Kreg;

    //----------------------------------------------------------
    // new Stuff
    //DG - specific stuff
    //ordering of the extra dofs:
    // [0,....,n-1]<- usual dofs on the patch
    // [n,...,ntot]<- extra dofs [n,...,n0-1,n0,...,n1-1,..., ntot] , where ni are the number of the neigboring patch corresponding to side.index()==i

    /// \brief contains the number of extra basis functions for a specific side of a patch
    /// m_numExtraBasis[patch][side.index()]. contains [4,5,.,7], i.e. sides, which dont have an interface have undefined value.
    std::vector< std::vector<index_t> > m_numExtraBasis;

    std::vector<index_t> m_numExtraBasisTotal;

    std::vector< std::vector< gsMatrix<index_t> > > m_boundaryBasis;

    std::vector<std::pair<patchDof,patchDof> >  m_elimExtraBasisConnection;

    std::vector<iFace::strategy> m_iCoupling;

    std::vector<std::map<patchSide,int> > m_patchSideToExtraIndex;

protected:
    /* auxiliary variables */
    using gsIETIAssembler<T>::m_scalingsKC;
    using gsIETIAssembler<T>::m_scalingsKii;
    using gsIETIAssembler<T>::m_scalingsKrr;


    using gsIETIAssembler<T>::orient;
};



} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIETIdGAssembler.hpp)
#endif
