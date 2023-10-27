/** @file gsRecipeAssembler.h

    @brief generic assembler on multipatch geometries based on recipe assembler

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, F. Buchegger
*/
#pragma once

#include <gsRecipeAssembler/gsRecipe.h>
#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssembler/gsPDEOperators.h>
#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsDomainIterator.h>

#include <gsMapUtils/gsMapperUtils.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsPde/gsPde.h>
#include <vector>


namespace gismo{


struct gsIntegrationRule
{
    gsQuadRule<real_t>         rule;
    gsDomainIterator<real_t>*  subdomains;

    ~gsIntegrationRule() {delete subdomains;}
};


/**
 * @brief The gsRecipeAssembler
 *
 *        There are two pieces for Galerkin method: the discretization of the solution space
 *        and that of the test space. Each of those is represented by a gsDiscretization object.
 */

enum {
    target=-1
};

class GISMO_EXPORT gsRecipeAssembler
{
public:
    typedef std::vector<gsPhysicalSpace*>           Discretization;
    typedef std::vector<gsPhysicalSpace::spacePtr>  SpaceList;
    typedef Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,index_t> PermMatrix;
    // matrix types
    typedef gsSparseMatrix<real_t>                   SMatT;
    typedef gsMatrix<real_t>                         FMatT;
protected:
    // domain and space discretization
    const gsMultiPatch<real_t>        &m_domain;
    Discretization                     m_space;

    // mapping from patch spaces to global space
    gsWeightMapper<real_t>            *m_map;
    PermMatrix                         m_permutation;
    std::vector<index_t>               m_shiftsSource;
    std::vector<index_t>               m_shiftsTarget;

    // matrices
    gsSparseMatrix<real_t>             m_elimM;
    gsMatrix<real_t>                   m_elimRhs;

    gsSparseMatrix<real_t>             m_sysM;
    gsSparseMatrix<real_t>             m_rhsMod;
    gsMatrix<real_t>                   m_rhs;

    // eliminated dofs from global space
    std::vector<std::vector<index_t> > m_eliminatedDofs;
    std::vector<index_t>               m_eliminatedTarget;
protected:
    // constructors

    gsRecipeAssembler(const gsMultiPatch<real_t> &domain)
        : m_domain(domain), m_map(NULL)
    {}

    gsRecipeAssembler(const gsMultiPatch<real_t> &domain, const Discretization &discretization)
        : m_domain(domain), m_space(discretization), m_map(NULL)
    {}

public:
    virtual ~gsRecipeAssembler();

public:

    /**
     * @brief eliminateDofs
     *        set some degrees of freedom as eliminated
     * @param dofs
     * @param isTarget are the dofs provides as patch space indices or as global indices?
     */
    void eliminateDofs(const std::vector<index_t> &dofs, int spaceId);

    /**
     * @brief eliminateDofs
     *        set some degrees of freedom as free, by removing them from the list of eliminated dofs
     * @param dofs
     * @param isTarget are the dofs provides as patch space indices or as global indices?
     */
    void freeDofs(const std::vector<index_t> &other, int spaceId);

    /**
     * @brief assemble
     *        assemble the system
     */
    virtual void assemble ();


    /**
     * @brief getSystemMatrix
     * @return a const reference to the system matrix
     */
    virtual const gsSparseMatrix<real_t>& getSystemMatrix()     const {return m_sysM;}

    /**
     * @brief getRhsModMatrix
     * @return a const reference to the rhs modification matrix, that is the linear operator
     *         that maps a given value of eliminated dofs to a given modification of the rhs
     */
    virtual const gsSparseMatrix<real_t>& getRhsModMatrix()     const {return m_rhsMod;}

    /**
     * @brief getEliminatedMatrix
     * @return a const reference to the system for the eliminated dofs
     */
    virtual const gsSparseMatrix<real_t>& getEliminatedMatrix() const {return m_elimM;}

    /**
     * @brief getSystemRhs
     * @return a const reference to the rhs matrix
     */
    virtual const gsMatrix<real_t>&       getSystemRhs()        const {return m_rhs;}

    /**
     * @brief  getEliminatedRhs
     * @return a const reference to the rhs for the eliminated dofs system
     */
    virtual const gsMatrix<real_t>&       getEliminatedRhs()    const {return m_elimRhs;}

    /**
     * @brief testtest
     * @return the shift index of the for each space (eliminated DoFs are not removed)
     */
    virtual std::vector<index_t>    getShifts(){return m_shiftsTarget;}


    /**
     * @brief reconstructSolution
     *        provides the coefficients for the space corresponding to a given solution for both
     *        free and eliminated dofs
     * @param space  index of the space in the vector
     * @param sysSol
     * @param elimDof
     * @return
     */
    virtual gsMatrix<real_t>  reconstructSolution(index_t space, const gsMatrix<real_t> &sysSol, const gsMatrix<real_t> &elimDof ) const;

    /**
     * @brief getSysSize
     * @return the size of the system matrix
     */
    virtual index_t getSysSize () const {return m_map->getNrOfTargets()-getElimSize ();}

    /**
     * @brief getRhsDim
     * @return the size of the system matrix
     */
    virtual index_t getRhsDim () const {return 1;}

     /**
     * @brief getFreeLimit
     * @return the size of the system matrix
     */
    virtual index_t getFreeLimit () const  {return getSysSize();}

     /**
     * @brief getElimSize
     * @return the size of the system matrix
     */
    virtual index_t getElimSize () const {return m_eliminatedTarget.size();}

    /**
     * @brief setSpace
     *        set the discretization space for the assembler.
     *        This function also reset the vectors of eliminated dofs.
     * @param space
     */
    virtual void setSpace (const Discretization &space);

    virtual const Discretization & getSpace ()
    { return m_space; }


    virtual void reset();
protected:
    // initialization functions

    /**
     * @brief init
     *        initialize data  before assembling
     */
    virtual void init();


    /**
     * @brief initMappers
     *        construct the mapper from the physical space mappers
     *        and initialize m_shiftsTarget.
     */
    virtual void initMappers();

    /**
     * @brief remapEliminatedDOFS
     *        change mapper so that it respect the eliminated DOFS.
     *        After this function is called it is assumed that
     *
     *        getSysSize
     *        getFreeLimit
     *        getRhsDim
     *
     *        return the correct numbers
     */
    virtual void remapEliminatedDOFS();

    /**
     * @brief initSystemMatrices
     *        prepare the matrices, this is done after the mapper has been initialized
     *        and adapted for the eliminated dofs
     */
    virtual void initSystemMatrices();

    /**
     * @brief getSysEstimatedEntries
     * @return the estimated number of nonzero entries per row in the system matrix
     */
    virtual gsVector<index_t> getSysEstimatedEntries();

    /**
     * @brief getSysRhsModEstimatedEntries
     * @return the estimated number of nonzero entries per row in the system matrix
     */
    virtual gsVector<index_t> getSysRhsModEstimatedEntries();

    /**
     * @brief getElimEstimatedEntries
     * @return the estimated number of nonzero entries per row in the eliminated matrix
     */
    virtual gsVector<index_t> getElimEstimatedEntries();

    /**
     * @brief postProcess
     *        called after the assembling is completed
     */
    virtual void postProcess();


    /**
     * @brief getPatchRecipe
     * @param patch
     * @return
     */
    virtual gsRecipe<real_t>    getPatchRecipe     (index_t patch)  = 0;
    /**
     * @brief getBoundaryRecipe
     * @param bc
     * @return
     */
    virtual gsRecipe<real_t>    getBoundaryRecipe  (patchSide s)  = 0;


    /**
     * @brief getPatchIntegration
     *        returns the used integration rule on each patch, can be reimplemented by subclasses
     * @param patch
     * @return
     */
    virtual gsIntegrationRule getPatchIntegration     (index_t patch );

    /**
     * @brief getPatchIntegration
     *        returns the used integration rule for a boundary condition, can be reimplemented by subclasses
     * @param patch
     * @return
     */
    virtual gsIntegrationRule getBoundaryIntegration  (patchSide ps);

    /**
     * @brief getPatchSpaces
     *        provides the list of evaluators for each patch
     * @param patch
     * @return
     */
    virtual SpaceList getPatchSpaces (index_t patch );

    /**
     * @brief getPatchMaxDegree
     * @param patch
     * @return
     */
    virtual gsVector<index_t> getPatchMaxDegree(index_t patch);

    virtual index_t getPatchMaxOverlap(index_t patch)
    { return (getPatchMaxDegree(patch).array()*2+1).prod() * totalScalarSize(); }
    virtual index_t totalScalarSize()
    {
        index_t total=0;
        for (size_t s=0; s<m_space.size();++s)
            total+=m_space[s]->targetDim;
        return total;
    }
protected:
    // utility functions for getting writing rules

protected:
    // standard writers
    typedef gsMatAndRhsModWriter<SMatT,SMatT>        SysW;
    typedef gsMatAndRhsModWriter<FMatT,gsNullWriter<real_t> > RhsW;
    // multiplier writer
    typedef gsMultiplierWriter<SysW>                 SymW;
    typedef gsShiftWriter<gsMultiplierWriter<SysW> > MulW;
    // writers for the eliminated system
    typedef gsBoundaryWriter<SMatT>                  EliSysW;
    typedef gsBoundaryWriter<FMatT>                  EliRhsW;
    // writers with coefficients
    typedef gsCoeffWriter<SysW>                      SysWC;
    typedef gsCoeffWriter<SymW>                      SymWC;

    gsLocalToGlobalMapper<real_t> *getSysWriter();
    gsLocalToGlobalMapper<real_t> *getLagrangeMultiplierWriter();
    gsLocalToGlobalMapper<real_t> *getRhsWriter();
    gsLocalToGlobalMapper<real_t> *getEliSysWriter();
    gsLocalToGlobalMapper<real_t> *getEliRhsWriter();

    gsLocalToGlobalMapper<real_t> *getSysWriter(real_t coef);
    gsLocalToGlobalMapper<real_t> *getLagrangeMultiplierWriter(real_t coef);

protected:
    virtual void assembleBoundaries ();
    virtual void assembleVolumes    ();
};

}

