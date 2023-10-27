/** @file gsIETIAssembler.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2014-12-03
*/


#pragma once
#include <gsUtils/gsStopwatch.h>

#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsBoxTopology.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsPreconditioner.h>

#include <gsUtils/gsSortedVector.h>
#include <gsIETI/gsIETIUtils.h>

#include <gsSolver/gsPatchPreconditionersCreator2.h>
#include <gsSolver/gsConjugateGradient.h>
#include <gsSolver/gsMinimalResidual.h>
#include <gsSolver/gsBramblePasciakCG.h>
#include <gsSolver/gsGMRes.h>

#include <iomanip>
namespace gismo {

template<typename T> class gsIETIJumpOperator;

struct gsIETIOptions;

/**
 * @brief gsIETIAssembler is the main assembling routine for the IETI method.
 * It assembles the stiffness matrices on each patch and calculates the LU factorization,
 * it also provides the bookkeeping for different splittings.

 */
template<typename T>
class gsIETIAssembler
{
public:
    /// \brief  the method used for facorizing the matrices
#if defined(GISMO_WITH_PARDISO)
    typedef typename gsSparseSolver<T>::PardisoLU sparseLUfact; //Threadsave
    typedef typename gsSparseSolver<T>::PardisoLLT sparseSPDfact; //Threadsave
  //  typedef typename gsSparseSolver<T>::PardisoLU sparseSPDfact; //Threadsave
#elif defined(_OPENMP)
    typedef typename gsSparseSolver<T>::LU sparseLUfact; //Threadsave
    typedef typename gsSparseSolver<T>::LU sparseSPDfact; //Threadsave
#elif(defined(GISMO_WITH_SUPERLU))
    typedef typename gsSparseSolver<T>::SuperLU sparseLUfact; //Not Threadsave
    typedef typename gsSparseSolver<T>::SuperLU sparseSPDfact; //Not Threadsave
#else
    typedef typename gsSparseSolver<T>::LU sparseLUfact; //Threadsave
    typedef typename gsSparseSolver<T>::LU sparseSPDfact; //Threadsave
#endif

    /// \brief  shortcut
    typedef std::pair<size_t, index_t> patchDof;

    typedef memory::shared_ptr<gsIETIAssembler> Ptr;
    typedef memory::unique_ptr<gsIETIAssembler> uPtr;

    friend class gsIETIJumpOperator<T>;
protected:

    class nonSymBlockPreconditioner;

    class orthProjToDualsubspace;

    class orthProjToDualofDualsubspace;

    class permutation;

protected:
    /// \brief  shortcut for assembling matrices
    typedef Eigen::Triplet<T> Trip;

    /**
     * @brief getEdgesFromBoundary given a 3D patch side, this method calculates all 4 Edges to this side
     * @param side the given patchSide of a 3D cube
     * @param e the vector containing the 4 Edges
     */
    virtual void getEdgesFromBoundary(const patchSide& side ,std::vector<Edge>& e)
    {
        gsMatrix<index_t> bI = m_basis.front().basis(side.patch).boundary(side);
        typename gsBasis<T>::uPtr bbasis = m_basis.front().basis(side.patch).boundaryBasis(side);
        bool hasDirichlet;

        for(int i=1;i<5;i++)
        {
            boxSide side2D(i);
            gsMatrix<index_t> bbI = bbasis->boundary(side2D);

            int n = bbI.rows();
            hasDirichlet=false;

            // if(n>2)
            {
                for(int j=1;j<n-1;j++)
                {
                    if(!m_locDofsMapper[side.patch][0].is_free((bI)((bbI)(j,0),0)))
                    {
                        hasDirichlet = true;
                        break;
                    }
                }
                e.push_back(Edge(hasDirichlet,n,orient(side2D.index()-1,side.index()), side, side2D));
            }

        }
    }

    virtual void printTime(gsStopwatch& time, std::string text)
    {
#pragma omp master
        {
            time.stop();
            gsInfo<<text<<time<<std::endl<<std::flush;
            time.restart();
        }
    }

    virtual void printWarn(std::string text)
    {
#pragma omp master
        {
            gsWarn<<text<<std::endl;
        }
    }

    virtual void print(std::string text)
    {
#pragma omp master
        {
            gsInfo<<text<<std::endl;
        }
    }


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
    virtual inline unsigned compCalc(int np,unsigned idx, int c) const
    {
        return c*m_basis.front()[np].size() +idx;
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
    virtual inline unsigned compCalcBack(int np,unsigned idx) const
    {
        return idx%m_basis.front()[np].size();
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
    virtual inline int getComp(int np,unsigned idx) const
    {
        return idx/m_basis.front()[np].size();
    }

    /**
     * @brief setInfo Populates the entries of the gsIETIInfo structure. Needs to be called after
         init_patchInterfaceSides();
         init_primalDofs();
         init_boundaryDofs(isMinimalEnergy);
         init_jumpOperatorData();

        and called before
         init_remainingDofs();

     */
    virtual void setInfo();


    /**
     * @brief init_mappers initializes the mappers for the bookkeeping. Needs to be called before all other init methods.
     */
    virtual void init_mappers();

    /**
     * @brief init_jumpOperatorData initializes the data for the jump operator B,
     * i.e. gsIETIJumpOperator. Needs to be called before setInfo(..);
     */
    virtual void init_jumpOperatorData();

    /**
     * @brief init_patchInterfaceSides Initializes the sides for each patch lieing on an interface.
     * Needs to be called before setInfo(..);
     */
    virtual void init_patchInterfaceSides();

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
    virtual void init_primalDofs();

    /**
     * @brief init_remainingDofs Initializes the splitting for primal and remaining dofs. This is only needed in case of
     *  minimumEnergy is false.
     *  Has to be called after all init_* functions and setInfo.
     */
    virtual void init_remainingDofs();

    /**
     * @brief init_boundaryDofs Initializes the splitting in boundary and remaining dofs. Has to be called after
     *    init_primalDofs()
     *  and called before
     *    setInfo(..)
     *    init_boundaryDofs();
     *
     * @param isMinimumEnergy which algorithm is used
     */
    virtual void init_boundaryDofs();

    virtual void init_InitialGuess();

    /**
     * @brief assembleLUofKC assembles the LU factorization of the blockmatrix [K C^T;C 0].
     *
     * due to stability reasons for very jumping coefficients, we store the LU factorization of:
     * [ K , \alpha C^T; \alpha C^T, 0], hence, if you solve with a rhs, you should use:
     * [F1 ,\alpha F2].
     * @param tripletList a vector of triplets representing the matrix K, already ordered as
     *          [Kbb Kbi; Kib Kii].
     */
    virtual void assembleLUofKC(gsSparseMatrix<T>& matrix, size_t np);

    /**
     * @brief assembleKrrKrpKpp builds the submatrices Krp, Kpp from the given vector of Triplets and stores a LU factorization of Krr
     * @param tripletList the vector of triplets representing the matrix K, already ordered asn
     *          [Kpp Kpr; Krp Krr].
     */
    virtual void assembleKrrKrpKpp(gsSparseMatrix<T>& matrix, size_t np);

    /**
     * @brief assembleKiiKibKbb builds the submatrices  Kib, Kbb from the given vector of Triplets and stores a LU factorization of Kii.
     * Required for preconditioning for all cases and for processing the rhs and the solution in case of minimal energy.
     * @param tripletList a vector of triplets representing the matrix K, already ordered as
     *          [Kbb Kbi; Kib Kii].
     */
    virtual void assembleKiiKibKbb(gsSparseMatrix<T>& matrix, size_t np);

    /**
     * @brief assembleSpp Assembles the schur complement w.r.t. the primal variables.
     *  In case of minimum Energy, it is the sumbmatrix  w.r.t. the primal varibales of the schur complement S of K w.r.t. to b.
     *  In the other case it is calculated as Kpp - Kpr* Krr^-1 * Krp
     *
     *  Needs to be called after
     *      assembleC();
     *      assembleLUofKC(tripletList);
     *   or
     *      assembleKrrKrpKpp(tripletList);
     *
     */
    virtual void assembleSpp(std::vector<gsMatrix<T> >& spp_loc);
    virtual void assembleSppLoc(gsMatrix<T>& spp_loc, size_t np);
    virtual void storeSpp(const gsMatrix<T>& in, gsMatrix<T>& out);

    /**
     * @brief assembleC assembes the constraint matrix C for the minimal Energy method, which is used to enforce the primal constraints
     */
    virtual void assembleC(size_t np);

    /**
     * @brief assembleRhs Assembles the rhs accoring to the desired splitting and/or transformation to the schur complement.
     *   Needs to be called after ALL assemble routines
     *
     * @param rhs_loc the patch-local contributions of the rhs
     */
    virtual void assembleRhs(const std::vector<gsMatrix<T> >& rhs_loc, size_t np);

    /**
     * @brief assembleRhsFreeElim Assembles the right hand side for the lagrange multipliers, which connect a eliminated and free dof
     */
    virtual void assembleRhsFreeElim();

    /**
     * @brief getInterfaceAverageValue This functions calculates for one patchSide of the parameter domain the interface average
     * for all nonzero basis function on this side
     * @param side The given patch side (2D or 3D)
     * @param average a vector of length boundaryBasis(side)->size() ,containig for each entry the average.
     * @param c the component of the basis (for scalar problems c = 0)
     * Note that for basis on the boundary of side, the average is not computed, the corresponding entry is 0
     */
    virtual void getInterfaceAverageValue(const patchSide& side, gsMatrix<T>& average, int c = 0);

    /**
     * @brief getFaceEdgeAverageValue This functions calculates for one Edge of a 3D cube in the parameter domain the average
     * for all nonzero basis function on this edge. For Edge averages of a 2D domain see getInterfaceAverageValue(..).
     * @param edge the given edge (only for 3D)
     * @param average  a vector of length edge.number ,containig for each entry the average.
     * @param c the component of the basis (for scalar problems c = 0)
     * Note that for basis on the boundary of side, the average is not computed, the corresponding entry is 0
     */
    virtual void getFaceEdgeAverageValue(const Edge& edge,gsMatrix<T>& average, int c = 0,bool eliminateCorners=true);

    /**
     * @brief checkPrimalDofStrategy The function checks if the primal dof strategy is appropriate for the desired
     * problem and used algorithm. This function changes the member m_primalDofstrat.
     * @param isMinimalEnergy true if minimal energy algorithm is used
     * @param dim the dimension of the desired object
     */
    virtual void checkPrimalDofStrategy(bool isMinimalEnergy,int dim);

    /// \brief  populates the orient matrix
    virtual void InstantiateEdge()
    {
        orient.resize(4,7);

        orient.col(static_cast<int>(boundary::west))<< 4,6,1,9;
        orient.col(static_cast<int>(boundary::east))<< 5,7,2,10;
        orient.col(static_cast<int>(boundary::south))<< 4,5,0,8;
        orient.col(static_cast<int>(boundary::north))<< 6,7,3,11;
        orient.col(static_cast<int>(boundary::front))<<1,2,0,3;
        orient.col(static_cast<int>(boundary::back))<< 9,10,8,11;
    }

    virtual void makeReorderingBoundInt(gsSparseMatrix<T>& matrix, size_t np)
    {
        gsVector<index_t> permVector(matrix.cols());
        for(index_t i = 0; i< matrix.cols();i++)
        {
            permVector[i] = m_glob2BoundInteriorIndex[np][i];
            if(m_globIsBoundIndex[np][i] ==false)
                permVector[i]+= info.dofsB[np];//shift the interior ones
        }

        // Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> perm(permVector);
        if(m_IETIoptions.KCSolver==IETILocalSolver::FastDiagonalization || m_IETIoptions.KCSolver==IETILocalSolver::Multigrid)
            m_K[np]=memory::make_shared(new gsSparseMatrix<T>(matrix));

        m_permMat[np]= Eigen::PermutationMatrix<Dynamic,Dynamic,index_t>(permVector);
        m_permInvMat[np]= m_permMat[np].inverse();
        matrix=matrix.twistedBy( m_permMat[np]);

    }
    virtual void makeReorderingPrimalRem(gsSparseMatrix<T>& matrix, size_t np)
    {
        gsVector<index_t> permVector(matrix.cols());
        for(index_t i = 0; i< matrix.cols();i++)
        {
            permVector[i] = m_glob2PrimalRemainingIndex[np][i];
            if(m_globIsPrimalIndex[np][i] ==false)
                permVector[i]+= info.dofsP[np];//shift the remaining ones
        }

        Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> perm(permVector);
        if(m_IETIoptions.KrrSolver==IETILocalSolver::FastDiagonalization || m_IETIoptions.KrrSolver==IETILocalSolver::Multigrid)
            m_K[np]=memory::make_shared(new gsSparseMatrix<T>(matrix));
        matrix=matrix.twistedBy(perm);
    }

    virtual void assembleInit();
    virtual void reserveSpace(std::vector<gsMatrix<T> >& dirDofs, std::vector<gsMatrix<T> >& rhs_loc);
    virtual void assembleLocal(gsAssembler<T>* A, gsSparseMatrix<T>& matrix, gsMatrix<T>& rhs, gsMatrix<T> &dirDofs, size_t np, const gsMultiPatch<T> &curSol);

    void constructKC(gsSparseMatrix<T> & KC, const gsSparseMatrix<T> &C, IETILocalSolver::solvers type, size_t np);
    virtual void getSchur(gsMatrix<T>& schur,size_t np);
    virtual void setupMGforKC(size_t np);
    void getPrimalInitialGuess(const gsMatrix<T> & allSol, gsMatrix<T>& sol, index_t c, int np) const ;

    virtual void extractPatch(size_t np, const gsMatrix<T>& rhs, gsMatrix<T>& rhsLocal) const;
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
    gsIETIAssembler(gsAssembler<T> &assembler);


    /// \brief  destructor
    virtual ~gsIETIAssembler()
    {
        if(m_IETIoptions.KCSolver==IETILocalSolver::Direct)
            freeAll(m_LU_KC);
        if(m_IETIoptions.KrrSolver==IETILocalSolver::Direct)
            freeAll(m_LU_Krr);
        if(m_IETIoptions.KiiSolver==IETILocalSolver::Direct)
            freeAll(m_LU_Kii);

        if(m_LU_Spp!=NULL)
            delete m_LU_Spp;
    }

    static gsOptionList defaultOptions();
    virtual void setOptions(const gsOptionList & opt);              ///< Set the options based on a gsOptionList

    /// \brief  returns the info structure
    const gsIETIInfo & getInfo() const {return info;}

    gsIETIOptions & getOptions() {return m_IETIoptions; }
    const gsIETIOptions & getOptions() const {return m_IETIoptions; }


    void printTiming() {m_timings.print(gsInfo);}

    /**
     * @brief solveKC solves [K ,C^T;C , 0]* sol = rhs. In order to obtain numerical stability, parts of the matrix
     * are scaled, this implies also a scaling of the rhs. Note, that the rhs is changed during the function call,
     * but the change is reverted at the end. So in fact, rhs stays unchanged.
     * @param np the actual patch
     * @param rhs the right hand side
     * @param sol the solution
     */
     template<bool exact, bool basis>
     void solveKC(unsigned np, gsMatrix<T>& rhs, gsMatrix<T>& sol) const;

    /**
     * @brief solveKii solves Kii* sol = rhs. In order to obtain numerical stability, parts of the matrix
     * are scaled, this implies also a scaling of the rhs. Note, that the rhs is changed during the function call,
     * but the change is reverted at the end. So in fact, rhs stays unchanged.
     * @param np the actual patch
     * @param rhs the right hand side
     * @param sol the solution
     */
     template<bool exact>
     void solveKii(unsigned np, gsMatrix<T>& rhs, gsMatrix<T>& sol) const;

    /**
     * @brief solveKr solves Krr* sol = rhs. In order to obtain numerical stability, parts of the matrix
     * are scaled, this implies also a scaling of the rhs. Note, that the rhs is changed during the function call,
     * but the change is reverted at the end. So in fact, rhs stays unchanged.
     * @param np the actual patch
     * @param rhs the right hand side
     * @param sol the solution
     */
    void solveKrr(unsigned np, gsMatrix<T>& rhs, gsMatrix<T>& sol) const;

    /**
     * @brief solveSpp solves Spp* sol = rhs.
     * @param rhs the right hand side
     * @param sol the solution
     */
    virtual void solveSpp(const gsMatrix<T>& rhs, gsMatrix<T>& sol) const;

    /**
     * @brief Applys the patch local stiffness matrix, reordered for B and I.
     */
    virtual void applyK(size_t np,const gsMatrix<T>& u, gsMatrix<T> & sol  )const
    {
        sol.setZero(u.rows(),u.cols());

        sol.topRows(info.dofsB[np]) = *m_Kbb[np]*u.topRows(info.dofsB[np]) + *m_Kbi[np]*u.bottomRows(info.dofsI[np]);
        sol.bottomRows(info.dofsI[np]) = (*m_Kib[np])*u.topRows(info.dofsB[np]) + *m_Kii[np]*u.bottomRows(info.dofsI[np]);

    }
    virtual void applyDP(const gsMatrix<T>& xP , std::vector<gsMatrix<T> >& uD) const;
    virtual void applyDP(unsigned np, const gsMatrix<T>& xP , gsMatrix<T>& uD) const;
    virtual void applyPD(const std::vector<gsMatrix<T> >& uD, gsMatrix<T>& xP) const;
    virtual void applyPD(unsigned np, const gsMatrix<T>& xD , gsMatrix<T>& uP) const;

    /**
     * @brief Applys the stiffness matrix, reordered for B and I.
     */
    virtual void applySpp(const gsMatrix<T>& uP, gsMatrix<T>& solP)const
    {
        solP.setZero(uP.rows(),uP.cols());

        solP = m_Spp*uP;

    }

    /// \brief  returns the conraint matrix C for one patch
    virtual const gsSparseMatrix<T> & getC(int patch) const {GISMO_ASSERT(!m_IETIoptions.opt.askSwitch("NoMinimumEnergy"), "Not available for this algorithm, choose NoMinimumEnergy as false"); return m_C[patch];}

    /// \brief  return the matrix Krp for one patch
    virtual const typename gsSparseMatrix<T>::Ptr & getKrp(int patch) const {GISMO_ASSERT(m_IETIoptions.opt.askSwitch("NoMinimumEnergy"), "Not available for this algorithm, choose NoMinimumEnergy as true");return m_Krp[patch];}

    /// \brief  return the matrix Kib for one patch
    virtual const typename gsSparseMatrix<T>::Ptr & getKib(int patch) const {return m_Kib[patch];}
    virtual const typename gsSparseMatrix<T>::Ptr & getKbi(int patch) const {return m_Kbi[patch];}

    /// \brief  return the matrix Kbb for one patch
    virtual const typename gsSparseMatrix<T>::Ptr& getKbb(int patch) const {return m_Kbb[patch];}

    /// \brief  return the matrix Kbb
    virtual const std::vector<typename gsSparseMatrix<T>::Ptr >& getKbb() const {return m_Kbb;}

    /// \brief return the timings and iteration count for Multigrid
    const gsIETITimings & getTiming() const {return m_timings;}

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
    virtual void assemble(const gsMultiPatch<T>& curSol=gsMultiPatch<T>());

    virtual void giveAssembledMatrices(std::vector<gsSparseMatrix<T> > &matrices, const gsMatrix<T> &rhs);
    virtual void setNewRhs(const gsMatrix<T>& rhs);
    /**
     * @brief numberRhs returns the number of used rhs
     */
    virtual index_t numberRhs() const  { return info.nRhs; }

    /**
     * @brief systemSize returns the number of unknown in the final system, i.e. the size of the solution vector used in the PCG.
     */
    virtual unsigned systemSize()  const { return m_IETIoptions.opt.getSwitch("SaddlePoint") ? info.dofTotalB+info.dofTotalI+info.dofTotalP+info.lagrangeMult : info.lagrangeMult;}

    /**
     * @brief getDirRhs returns ICO nCoupledElimDofs() == true, a vector for the rhs of the system corresponding to B.
     * @return the rhs vector b.
     */
    virtual const gsMatrix<T>& getDirRhs() const {return m_rhs_dir;}

    /**
     * @brief nCoupledElimDofs checks if there are any lagrange multip. which couple a dirichlet and neumann boundary.
     * @return true if there are none.
     */
    virtual bool nCoupledElimDofs() const {return m_freeElimLagr.empty();}

    /**
     * @brief getRawRhs   gives you in case of minimal energy the corresponding components of fb - Kbi * Kii^-1 fi
     * and for not minimal energy fp, fr  (fp entries of f corresponding to the vertices)
     * @param rhs_p the primal component of the corresponding rhs
     * @param rhs_d the remaining component of the corresponding rhs
     */
    virtual void getRawRhs(gsMatrix<T> & rhs_p,std::vector<gsMatrix<T> >& rhs_d) const{
        rhs_p = m_rhs_p;
        rhs_d = m_rhs_d;
    }

    virtual void getRawRhs(gsMatrix<T> & rhs_p, std::vector<gsMatrix<T> >& rhs_b,std::vector<gsMatrix<T> >& rhs_i) const{
        rhs_p = m_rhs_p;
        rhs_b = m_rhs_d;
        rhs_i = m_rhs_i;
    }

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
    virtual void processSolution(const gsMatrix<T>& uP, const std::vector<gsMatrix<T> >& u2, gsMatrix<T>& solVec) const;
    virtual void processSolution(typename gsMatrix<T>::BlockView& view, gsMatrix<T>& solVec) const;

    /**
     * @brief embedding transforms a representation of the solution vector from a primal and remaining to patch wise values on the interface
     * @param xP the global primal component ( totalPdofs x nRhs)
     * @param x2 the local remaining/dual component for each patch (if isMinimumEnergy: dofsB[np] x nRhs ; else dofsR[np] x nRhs)
     * @param w the patch wise values on the interface  (dofsB[np] x nRhs)
     */
    virtual void embedding(gsMatrix<T> const & xP,std::vector<gsMatrix<T> >const & x2,std::vector<gsMatrix<T> >& w) const;
    virtual void embedding(gsMatrix<T> const & xP, gsMatrix<T> const & x2, size_t np, gsMatrix<T>& w) const;
    /**
     * @brief embeddingTrans transforms a representation of the solution vecor form a patch wise interface values to global primal component and a remaining/dual component
     * @param w the patch wise values on the interface  (dofsB[np] x nRhs)
     * @param uP the global primal component ( totalPdofs x nRhs)
     * @param u2 the local remaining/dual component for each patch (if isMinimumEnergy: dofsB[np] x nRhs ; else dofsR[np] x nRhs)
     */
    virtual void embeddingTrans(const std::vector<gsMatrix<T> >& w,gsMatrix<T>&  uP, std::vector<gsMatrix<T> > &  u2) const;
    virtual void embeddingTrans(const gsMatrix<T> & w,size_t np, gsMatrix<T>&  uP, gsMatrix<T> &  u2) const;
    /**
     * @brief assemblePrimal assembles the patch contributions of the primal variables to the global representation
     * @param wi patch local contribution (dofsP[patch] x nRhs)
     * @param patch the actual patch
     * @param uP the updated global primal variable vector (totalPdofs x nRhs)
     */
    virtual void inline assemblePrimal(gsMatrix<T> const & wi, int patch, gsMatrix<T>& uP) const
    {
#if defined(_OPENMP) && _OPENMP > 200505
        // this OpenMP command does not exist prior to OpenMP 3.0
        if(omp_get_level()!=0)
#else
        if(false)
#endif
        {
            for(int rhs=0; rhs<wi.cols();rhs++)
                for(index_t i=0; i<info.dofsP[patch];i++)
#pragma omp atomic
                    uP(m_pDofsLoc2Glob[patch][i],rhs) += wi(i,rhs);
        }
        else
            for(index_t i=0; i<info.dofsP[patch];i++)
                uP.row(m_pDofsLoc2Glob[patch][i]) += wi.row(i);
    }


    /**
     * @brief distributePrimal distributes the global primal variable vector to a certain patch
     * @param xP the global primal variable vector (totalPdofs x nRhs)
     * @param patch the actual patch
     * @param xPi patch local contribution (dofsP[patch] x nRhs)
     */
    virtual void inline distributePrimal(gsMatrix<T>const & xP,int patch,gsMatrix<T> & xPi) const
    {
        xPi.setZero(info.dofsP[patch],xP.cols());
        for(index_t i=0; i<info.dofsP[patch];i++) {
            xPi.row(i) += xP.row(m_pDofsLoc2Glob[patch][i]);
        }
    }

    void calcRhsNorm(std::vector<gsMatrix<T> >& rhs)
    {
        gsDofMapper stdMapper = m_assembler->system().colMapper(0);
        gsVector<T> rh(stdMapper.freeSize());
        rh.setZero();
        for(size_t np=0; np<info.numberPatches;np++)
        {

            int sz  =  m_basis[0][np].size();
            for (index_t i = 0;i < sz ; ++i)
            {
                int idx = m_locDofsMapper[np][0].index(i);
                int globI = stdMapper.index(i,np);

                if(stdMapper.is_free(i,np))
                {
                    rh(globI)+= rhs[np](idx,0);
                }
            }


        }
        //    gsDebugVar(rh.transpose());
        m_rhsNorm=rh.norm();
    }
    void calcMatrix(std::vector<gsSparseMatrix<T> >& matrices)
    {
        gsDofMapper stdMapper = m_assembler->system().colMapper(0);
        gsSparseMatrix<T> mat(m_assembler->numDofs(),m_assembler->numDofs());
        mat.reserve(m_assembler->numColNz());
        for(size_t np=0; np<info.numberPatches;np++)
        {

            int sz  =  m_basis[0][np].size();
            for (index_t i = 0;i < sz ; ++i)
            {
                int idxi = m_locDofsMapper[np][0].index(i);
                int globI = stdMapper.index(i,np);
                if(stdMapper.is_free(i,np))
                {
                    for (index_t j = 0;j < sz ; ++j)
                    {
                        int idxj = m_locDofsMapper[np][0].index(j);
                        int globJ = stdMapper.index(j,np);

                        if(stdMapper.is_free(j,np))
                            mat(globI,globJ)+= matrices[np](idxi,idxj);
                    }
                }
            }


        }
        gsDebugVar(mat.toDense());
    }

    T getRhsNorm() {return m_rhsNorm;}



    void projectToDualSubspaceA(gsMatrix<T>& u, size_t np) const
    {
        gsStopwatch time;
        u.topRows(m_PhiA[np].rows()) -= m_PhiA[np]*m_exC[np]*u.topRows(info.dofsB[np]+info.dofsI[np]);
        m_timings.projection[np]+=time.stop();
    }
    void projectToDualSubspace(gsMatrix<T>& u, size_t np) const
    {
        gsStopwatch time;
        u.topRows(m_Phi[np].rows()) -= m_Phi[np]*m_C[np]*u.topRows(info.dofsB[np]);
        m_timings.projection[np]+=time.stop();
    }



public:
    gsIETIAssembler() {}


protected:

    /// @brief pointer to the used assembler
    gsAssembler<T>* m_assembler;

    /// @brief The multipatch domain
    gsMultiPatch<T> m_patches;

    /// \brief   The discretization bases corresponding to \a m_patches and to
    /// the number of solution fields that are to be computed
    /// m_basis[i]: The multi-basis for unknown i
    std::vector< gsMultiBasis<T> > m_basis;

    /// \brief backup of the "full" boundary conditions
    gsBoundaryConditions<T> m_bConditions;

    /// \brief  An info Structure
    gsIETIInfo info;

    /// \brief  contains global information of the primal dofs to the system (for each component)
    std::vector<gsDofMapper> m_primalDofMapper;

    /// \brief  global standart mapper represent the matched patches (for each component)
    std::vector<gsDofMapper> m_stdMapper;

    /// \brief  local dof mapper, map a basis function of a patch to the real (where elimination happend) index (locally) and for each component .
    /// m_locDofsMapper[np] -> all mappers on patch np.
    std::vector<std::vector< gsDofMapper> > m_locDofsMapper;

    /// \brief  contains information of coupled edges in 3D
    gsDofMapper m_edgeMapper;

    ///The stored primal dof mappings dont have a dirichlet boundary, so you dont have to care.
    //(these maps works for uPi -> uP)

    /// \brief  local number of patch k to the global number (in uP).  length: (nPatches)
    std::vector< std::vector<index_t> > m_pDofsLoc2Glob;

    //for all stored dofs it is assumed that the dirichlet boundary is not eliminated!!!!
    //use mapper.index(i) to get the real index in the matrix.

    /// \brief  contains the local numbers of Vertex Primal dofs for each patch
    std::vector< gsSortedVector<index_t> > m_primalVdofs;

    /// \brief  Contains the Edges for each patch in 3D
    std::vector< gsSortedVector<Edge > > m_primalEdges;

    /// \brief  contains the local number of the remaining dofs for each patch
    std::vector< std::vector<index_t> > m_remDofs;

    /// \brief  contains the local numbers of interface dofs for each patch
    std::vector< gsSortedVector<index_t> > m_boundDofs;

    /// \brief  on a patch k contains for a remaining dof i true if it is also on the interface, i.e. is in m_boundDofs[k]
    std::vector<std::vector<bool> >m_remDofIsAlsoBoundDof;

    /// \brief  contains for each patch the the interface sides
    std::vector< std::vector<patchSide> > m_patchISides;

    /// \brief  contains for each patch the the interface sides
    std::vector< std::vector<patchSide> > m_patchAveragesSides;

    /// \brief  boundary conditions for each patch
    std::vector<gsBoundaryConditions<T> >  m_bc_loc;

    /// \brief  maps the global (patchwise) dof to the local interface index or interior index
    std::vector< std::vector<index_t> > m_glob2BoundInteriorIndex;

    /// \brief  contains true if a global(patchwise) index is a interface index
    std::vector< std::vector<bool> > m_globIsBoundIndex;

    /// \brief  maps a global index to the repective primal or remaining index (ONLY FOR noMinimalEnergy)
    std::vector< std::vector<index_t> > m_glob2PrimalRemainingIndex;

    /// \brief  contains true if a global index is a primal one (ONLY FOR noMinimalEnergy)
    std::vector< std::vector<bool> > m_globIsPrimalIndex;

    /// \brief  lagrange mult -> the pair of patchDofs
    std::vector<std::pair<patchDof,patchDof> > m_lagrangeTable;

    /// \brief  global interface dof -> merged nodes+patches
    std::map<index_t,std::set<patchDof> > m_globalConnectionTable;

    /// \brief  lagrange multipl which couple a eliminated and free dof
    std::vector<index_t> m_freeElimLagr;
    //------------------------------------------------------

    /// \brief  The S-orth Local Basis of the Primal variable Space
    std::vector<gsMatrix<T> > m_Phi; //in [B,I] or [B] order (dependent on saddle point or not)
    std::vector< gsMatrix<T> > m_PhiA; //in [normal] order

    /// \brief  Local Lu factorization of [K C^T,C 0]
    std::vector<sparseLUfact* >m_LU_KC;

    /// \brief  Local Matrix representation of C
    std::vector<gsSparseMatrix<T> > m_C;

    /// \brief  Local LU factorization of Krr
    std::vector<sparseSPDfact* >m_LU_Krr;

    /// \brief  Global Schur complement wrt. the primal Variables
    sparseSPDfact* m_LU_Spp;

    /// \brief  Local Matrix Representation of Krp
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Krp;

    /// \brief  Local representations of Kpp
    std::vector<gsMatrix<T> > m_Kpp;

    // Vector rhs
    /// \brief  Primal components of the rhs
    gsMatrix<T> m_rhs_p;

    /// \brief  Remaining or dual components of the rhs
    std::vector<gsMatrix<T> > m_rhs_d;

    /// \brief  Interior components of the rhs
    std::vector<gsMatrix<T> > m_rhs_i;

    /// \brief  dirichlet dofs  (if available, e.g. if using nitsche)
    std::vector<gsMatrix<T> > m_dirDofs;

    /// \brief  rhs vector containing the dirichlet values for coupled  eliminated and free dofs (for gsIETIJumpOperator)
    gsMatrix<T> m_rhs_dir;

    /// \brief  Local LU factorization of Kii
    std::vector<sparseSPDfact* >m_LU_Kii;

    /// \brief  Local matrix representation of Kib
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kib;
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kbi;


    /// \brief  Local matrix representation of Kbb
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kbb;
    //----------------------- inexact part
    gsSparseMatrix<T> m_Spp;

    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kii;//<- Ptr To save memory when using MG

    std::vector<typename gsLinearOperator<T>::Ptr> m_inexact_Kii_BS;
    std::vector<real_t>  m_maxEigKii;

    std::vector<typename gsPreconditionerOp<T>::Ptr> m_inexact_Kii_MG;
    //std::vector<std::vector< gsSparseMatrix<real_t, RowMajor> > > m_transferMatrices_Kii; //TODO: Keep them

    std::vector<Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > m_permMat;
    std::vector<Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > m_permInvMat;

    std::vector<gsBlockOp<>::Ptr > m_KC;
    std::vector<gsMatrix<T> > m_S;
    std::vector<gsBlockOp<>::Ptr> m_inexact_KC;
    std::vector<real_t> m_maxEigKC;

    std::vector<gsSparseMatrix<T> > m_exC;
    std::vector<typename gsSparseMatrix<T>::Ptr > m_K; //<- Ptr //multigrid needs the original matrix
  //  std::vector<gsSparseMatrix<T> > m_Kreg; //regularized matrix for pure neumann
  //  std::vector<std::vector< gsSparseMatrix<real_t, RowMajor> > > m_transferMatrices_KC; //TODO: Keep them

    std::vector<std::vector< Eigen::PermutationMatrix<Dynamic,Dynamic,index_t> > > m_primalPerm;
    std::vector<std::vector<boxCorner> > m_primalCorners;

    //std::vector<gsMinimalResidual<> > m_IterBasis;
    mutable std::vector<gsBramblePasciakCG<gsBPCG_Types::SchoeberlZulehner>::Ptr > m_KCBasis;
    mutable std::vector<gsConjugateGradient<>::Ptr> m_KCSolve;
    mutable std::vector<gsConjugateGradient<>::Ptr> m_KiiSolve;


protected:

    /// \brief  information of the edge number. (e.g. the 2. edge on the boundary left has edge number 8, -> orient(2,bounday::left.index() )
    gsMatrix<int> orient;

    std::vector<Eigen::IterScaling<gsSparseMatrix<T> > > m_scalingsKC;
    std::vector<Eigen::IterScaling<gsSparseMatrix<T> > > m_scalingsKii;
    std::vector<Eigen::IterScaling<gsSparseMatrix<T> > > m_scalingsKrr;

    bool m_isInit;
    T m_rhsNorm;

    gsIETIOptions m_IETIoptions;

    mutable gsIETITimings m_timings;
    bool m_useGivenMats;
    std::vector<gsSparseMatrix<T> > m_tempMatrix;
    std::vector<gsMatrix<T> > m_tempRhs;
};




template <typename T>
class gsIETIAssembler<T>::nonSymBlockPreconditioner : public gsLinearOperator<T>
{
public:
    /// Shared pointer for gsIdentityOp
    typedef memory::shared_ptr< nonSymBlockPreconditioner > Ptr;

    /// Unique pointer for gsIdentityOp
    typedef memory::unique_ptr< nonSymBlockPreconditioner > uPtr;

    /// Constructor taking the dimension of the identity operator
    nonSymBlockPreconditioner(typename gsLinearOperator<T>::Ptr precA, typename gsLinearOperator<T>::Ptr precS,typename gsLinearOperator<T>::Ptr B )
        : m_precA(precA), m_precS(precS), m_B(B)
    {

    }

    /// Make command returing a shared pointer
    static Ptr make(typename gsLinearOperator<T>::Ptr precA, typename gsLinearOperator<T>::Ptr precS, typename gsLinearOperator<T>::Ptr B )
    { return memory::make_shared( new nonSymBlockPreconditioner(precA,precS,B) ); }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        x.setZero(rows(),input.cols());
        m_precA->apply(input.topRows(m_precA->cols()),temp1);
        x.topRows(m_precA->rows())=temp1;
        m_B->apply(temp1,temp2);
        x.bottomRows(m_precS->rows())=temp2;
        m_precS->apply(input.bottomRows(m_precS->rows()),temp2);
        x.bottomRows(m_precS->rows())+=temp2;
    }

    index_t rows() const {return m_precA->rows()+m_B->rows();}

    index_t cols() const {return rows();}
private:
    typename gsLinearOperator<T>::Ptr m_precA;
    typename gsLinearOperator<T>::Ptr m_precS;
    typename gsLinearOperator<T>::Ptr m_B;
    mutable gsMatrix<T> temp1,temp2;
};

template <typename T>
class gsIETIAssembler<T>::orthProjToDualsubspace : public gsLinearOperator<>
{
public:
    /// Shared pointer for gsIdentityOp
    typedef memory::shared_ptr< orthProjToDualsubspace > Ptr;

    /// Unique pointer for gsIdentityOp
    typedef memory::unique_ptr< orthProjToDualsubspace > uPtr;

    /// Constructor taking the dimension of the identity operator
    orthProjToDualsubspace(const gsMatrix<T>& Phi, const gsSparseMatrix<T>& C, const gsIETIInfo & info, size_t np) : m_Phi(Phi), m_C(C), m_info(info), np(np){}

    /// Make command returing a shared pointer
    static Ptr make(const gsMatrix<T>& Phi, const gsSparseMatrix<T>& C, const gsIETIInfo & info, size_t np)
    { return memory::make_shared( new orthProjToDualsubspace(Phi,C,info,np) ); }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x=input;
        x.topRows(m_Phi.rows()) -= m_Phi*m_C*x.topRows(m_C.cols());
    }

    index_t rows() const {return m_Phi.rows();}

    index_t cols() const {return m_C.cols();}
private:
    const gsMatrix<T>& m_Phi;
    const gsSparseMatrix<T>& m_C;
    const gsIETIInfo & m_info;
    size_t np;

};

template <typename T>
class gsIETIAssembler<T>::orthProjToDualofDualsubspace : public gsLinearOperator<>
{
public:
    /// Shared pointer for gsIdentityOp
    typedef memory::shared_ptr< orthProjToDualofDualsubspace > Ptr;

    /// Unique pointer for gsIdentityOp
    typedef memory::unique_ptr< orthProjToDualofDualsubspace > uPtr;

    /// Constructor taking the dimension of the identity operator
    orthProjToDualofDualsubspace(const gsMatrix<T>& Phi, const gsSparseMatrix<T>& C, const gsIETIInfo & info, size_t np) : m_Phi(Phi), m_C(C), m_info(info), np(np){}

    /// Make command returing a shared pointer
    static Ptr make(const gsMatrix<T>& Phi, const gsSparseMatrix<T>& C, const gsIETIInfo & info, size_t np)
    { return memory::make_shared( new orthProjToDualofDualsubspace(Phi,C,info,np) ); }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x=input;
        x.topRows(m_C.cols()) -= m_C.transpose()*m_Phi.transpose()*x.topRows(m_Phi.rows());
    }

    index_t rows() const {return m_C.cols();}

    index_t cols() const {return m_Phi.rows();}
private:
    const gsMatrix<T>& m_Phi;
    const gsSparseMatrix<T>& m_C;
    const gsIETIInfo & m_info;
    size_t np;
};

template <typename T>
class gsIETIAssembler<T>::permutation : public gsLinearOperator<>
{
public:
    /// Shared pointer for gsIdentityOp
    typedef memory::shared_ptr< permutation > Ptr;

    /// Unique pointer for gsIdentityOp
    typedef memory::unique_ptr< permutation > uPtr;

    /// Constructor taking the dimension of the identity operator
    permutation(const Eigen::PermutationMatrix<Dynamic,Dynamic,index_t>& perm) : m_perm(perm) {}

    /// Make command returing a shared pointer
    static Ptr make(const Eigen::PermutationMatrix<Dynamic,Dynamic,index_t>& perm)
    { return memory::make_shared( new permutation(perm) ); }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x = m_perm*input;
    }

    index_t rows() const {return m_perm.rows();}

    index_t cols() const {return m_perm.cols();}
private:
    const Eigen::PermutationMatrix<Dynamic,Dynamic,index_t>& m_perm;
};



} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIETIAssembler.hpp)
#endif
