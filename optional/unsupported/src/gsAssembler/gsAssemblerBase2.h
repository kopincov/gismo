/** @file gsAssemblerBase.h

    @brief Provides generic assembler routines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <gsCore/gsBasisRefs.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsStdVectorRef.h>
#include <gsCore/gsAffineFunction.h> // needed by DG
#include <gsAssembler/gsVisitorDg2.h> // should be improved ASAP

#include <gsAssembler/gsAssemblerOptions.h>

namespace gismo
{

/** @brief The assembler class provides to generic routines for volume
 * and boundary integrals that are used for for matrix and rhs
 * generation
*/
template <class T>
class gsAssemblerBase2
{
private:
    typedef gsStdVectorRef<gsDofMapper> gsDofMappers;

public:


    /// @brief default constructor
    /// \note none of the data fields are inititalized, use
    /// additionally an appropriate initialize function
    gsAssemblerBase2()
    {}


    gsAssemblerBase2(const gsAssemblerBase2& other)
    {
        m_patches     = other.m_patches;
        m_bases       = other.m_bases;
        m_bConditions = other.m_bConditions;
    }


    /// @brief Constructor using a multipatch domain
    /// \note Rest of the data fields should be initialized in a
    /// derived constructor
    gsAssemblerBase2(const gsMultiPatch<T> & patches) :
        m_patches(patches)
    { }

    /// @brief Constructor using a multipatch domain, a
    /// vector of multibases and the boundary conditions.
    /// \note Rest of the data fields should be initialized in a
    /// derived constructor
    gsAssemblerBase2(const gsMultiPatch<T>                    & patches,
                     std::vector< gsMultiBasis<T> > const     & bases,
                     gsBoundaryConditions<T> const            & bconditions) :
        m_patches(patches),
        m_bases(bases),
        m_bConditions(bconditions)
    { }

    virtual ~gsAssemblerBase2()
    { }

public:
    /// @brief Intitialize function for, sets data fields
    /// using a multi-patch, a vector of multi-basis and boundary conditions.
    /// \note Rest of the data fields should be initialized in the
    /// derived function initializePdeSpecific() .
    void initialize(gsMultiPatch<T>                           patches,
                    gsStdVectorRef<gsMultiBasis<T> > const  & bases,
                    gsBoundaryConditions<T> const           & bconditions)
    {
        m_patches = give(patches);
        m_bases.clear();
        m_bases = bases;
        m_bConditions = bconditions;

        initializePdeSpecific();

    }

    /// @brief Intitialize function for, sets data fields
    /// using a multi-patch, a multi-basis and boundary conditions.
    /// \note Rest of the data fields should be initialized in the
    /// derived function initializePdeSpecific() .
    void initialize(gsMultiPatch<T>                 patches,
                    gsMultiBasis<T> const         & bases,
                    gsBoundaryConditions<T> const & bconditions)
    {
        m_patches = give(patches);
        m_bases.clear();
        m_bases.push_back(bases);
        m_bConditions = bconditions;

        initializePdeSpecific();

    }

    /// @brief Intitialize function for single patch assembling, sets data fields
    /// using a gsGeometry, a basis reference for each component (vector) and
    /// boundary conditions. Works for scalar and vector valued PDES
    /// \note Rest of the data fields should be initialized in the
    /// derived function initializePdeSpecific() .
    void initializeSinglePatch(const gsGeometry<T>           & patch,
                               const gsBasisRefs<T>          & basis,
                               gsBoundaryConditions<T> const & bconditions)
    {
        m_patches = gsMultiPatch<T>(patch);
        m_bConditions = bconditions;

        m_bases.clear();
        for(size_t c=0;c<basis.size();c++)
            m_bases.push_back(gsMultiBasis<T>(basis[c]));

        initializePdeSpecific();
    }

public:

    /// @brief Returns a minimal copy of the given Assembler, i.e. containing only the members
    /// which cannot initialized via initializePdeSpecific(). Creates a copy of those functions, which are not threadsave.
    /// This function is of particular importance when using the parallel version of IETI
    virtual gsAssemblerBase2* minimalClone() const {GISMO_NO_IMPLEMENTATION; }

public:

    /// @brief Generic assembly routine for volume or boundary integrals
    template<class ElementVisitor>
    void apply(ElementVisitor & visitor,
               index_t patchIndex = 0,
               boxSide side = boundary::none)
    {
        //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";
        
        const gsBasisRefs<T> bases(m_bases, patchIndex);
        const gsDofMappers mappers(m_dofMappers);
        
        gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
        gsVector<T> quWeights; // Temp variable for mapped weights
        unsigned evFlags(0);

        // Initialize
        visitor.initialize(bases, QuRule, evFlags);

        // Initialize geometry evaluator -- TODO: Initialize inside visitor
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_patches[patchIndex]));

        // Initialize domain element iterator -- using unknown 0
        typename gsBasis<T>::domainIter domIt = bases[0].makeDomainIterator(side);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

            // Perform required evaluations on the quadrature nodes
            visitor.evaluate(bases, /* *domIt,*/ *geoEval, quNodes);
            
            // Assemble on element
            visitor.assemble(*domIt, *geoEval, quWeights);
            
            // Push to global matrix and right-hand side vector
            visitor.localToGlobal(mappers, m_ddof, patchIndex, m_matrix, m_rhs);
        }
    }


    /// @brief Generic assembly routine for patch-interface integrals
    template<class InterfaceVisitor>
    void apply(InterfaceVisitor & visitor,
               const boundaryInterface & bi)
    {
        //gsDebug<<"Apply DG on "<< bi <<".\n";
        
        const gsDofMappers mappers(m_dofMappers);
        const gsAffineFunction<T> interfaceMap(m_patches.getMapForInterface(bi));

        const index_t patch1      = bi.first().patch;
        const index_t patch2      = bi.second().patch;
        const gsBasis<T> & B1 = m_bases[0][patch1];// (!) unknown 0
        const gsBasis<T> & B2 = m_bases[0][patch2];
        
        gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
        gsVector<T> quWeights;         // Mapped weights
        // Evaluation flags for the Geometry map
        unsigned evFlags(0);

        const int bSize1      = B1.numElements( bi.first() .side() );
        const int bSize2      = B2.numElements( bi.second().side() );
        const int ratio = bSize1 / bSize2;
        GISMO_ASSERT(bSize1 >= bSize2 && bSize1%bSize2==0,
                     "DG assumes nested interfaces. Got bSize1="<<
                     bSize1<<", bSize2="<<bSize2<<"." );
        
        // Initialize
        visitor.initialize(B1, B2, QuRule, evFlags);
        
        // Initialize geometry evaluators
        typename gsGeometryEvaluator<T>::uPtr geoEval1(getEvaluator(evFlags, m_patches[patch1]));
        typename gsGeometryEvaluator<T>::uPtr geoEval2(getEvaluator(evFlags, m_patches[patch2]));

        // Initialize domain element iterators
        typename gsBasis<T>::domainIter domIt1 = B1.makeDomainIterator( bi.first() .side() );
        typename gsBasis<T>::domainIter domIt2 = B2.makeDomainIterator( bi.second().side() );

        //typename gsBasis<T>::domainIter domIt = B2.makeDomainIterator(B2, bi);
        
        int count = 0;
        // iterate over all boundary grid cells on the "left"
        for (; domIt1->good(); domIt1->next() )
        {
            count++;
            // Get the element of the other side in domIter2
            //domIter1->adjacent( bi.orient, *domIter2 );
            
            // Compute the quadrature rule on both sides
            QuRule.mapTo( domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights);
            interfaceMap.eval_into(quNodes1,quNodes2);

            // Perform required evaluations on the quadrature nodes
            visitor.evaluate(B1, *geoEval1, B2, *geoEval2, quNodes1, quNodes2);

            // Assemble on element
            visitor.assemble(*domIt1,*domIt2, *geoEval1, *geoEval2, quWeights);
            
            // Push to global patch matrix (m_rhs is filled in place)
            visitor.localToGlobal(mappers,m_ddof, patch1, patch2, m_matrix, m_rhs);
            
            if ( count % ratio == 0 ) // next master element ?
            {
                domIt2->next();
            }

        }
    }

    /// @brief forces the Assembler to calculete the Dirichlet dofs
    /// without calling the assemble function. Use this function with care.
    /// This is especially usefull if one has a solution vector,
    /// but has called the assemble function, or just wants to calculate the dirichlet values.
    virtual void computeDirichletDofs() {GISMO_NO_IMPLEMENTATION;}

    virtual T penalty(int k) const {GISMO_UNUSED(k); GISMO_NO_IMPLEMENTATION;}

public:

    /// @brief Main assemble routine
    virtual void assemble() {GISMO_NO_IMPLEMENTATION}

    /// @brief Main non-linear assemble routine with input from
    /// current solution
    virtual void assemble(const gsMultiPatch<T> & curSolution)
    {GISMO_UNUSED(curSolution); GISMO_NO_IMPLEMENTATION}

    /// @brief Reconstruct solution from computed solution vector
    virtual void constructSolution(const gsMatrix<T>& solVector,
                                   gsMultiPatch<T>& result) const
    {GISMO_UNUSED(solVector); GISMO_UNUSED(result); GISMO_NO_IMPLEMENTATION}

    /// @brief Update solution by adding the computed solution vector
    /// to the current solution
    virtual void updateSolution(const gsMatrix<T>& solVector,
                                gsMultiPatch<T>& result) const
    {GISMO_UNUSED(solVector); GISMO_UNUSED(result); GISMO_NO_IMPLEMENTATION}


    /// @brief Returns the used-dG Visitor (should be an abstract class)
    virtual gsVisitorDg2<T> visitorDg(const boundaryInterface& bi) {GISMO_UNUSED(bi);GISMO_NO_IMPLEMENTATION;}

public:

    /// @brief Return the multipatch.
    const gsMultiPatch<T> & patches() const { return m_patches; }

    /// @brief Return the multi-basis
    const gsMultiBasis<T> & multiBasis(index_t k) const { return m_bases[k]; }

    /// @brief Returns the number of multi-bases
    size_t numMultiBasis() const {return m_bases.size(); }

    /// @brief Return the DOF mapper for unknown \em i.
    const gsDofMapper& dofMapper(unsigned i = 0) const     { return m_dofMappers[i]; }

    /// @brief Returns the number of dofMappers (corresponds to the number of components)
    size_t numDofMappers() const {return m_dofMappers.size();}

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const { return m_matrix; }

    /// @brief Returns true if only the lower triangular part of the
    /// matrix is computed (for symmetric problems)
    virtual bool isSymmertric() const { return false; }

    /// @brief Returns the Dirichlet values (if applicable)
    const gsMatrix<T> & dirValues() const { return m_ddof; }

    /// @brief Sets any Dirichlet values to homogeneous (if applicable)
    void homogenizeDirichlet() { m_ddof.setZero(); }

    /// @brief Returns the left-hand side vector(s)
    /// ( multiple right hand sides possible )
    const gsMatrix<T> & rhs() const { return m_rhs; }

    /// @brief Returns the number of (free) degrees of freedom
    int numDofs() const { return m_dofs; }
    
    /// @brief Returns the options of the assembler
    const gsAssemblerOptions & options() const { return m_options; }

    /// @brief Returns the boundary conditions (if applicable)
    const gsBoundaryConditions<T> & bConditions() const {return m_bConditions;}

protected:

    /// @brief Prototype for initializing PDE specific members. In order to provide
    /// working init functions, it SHOULD (must) be overridden in the derived class.
    virtual void initializePdeSpecific() {GISMO_NO_IMPLEMENTATION;}

protected:

    // *** Input data ***

    /// The multipatch domain
    gsMultiPatch<T> m_patches;

    /// The discretization bases corresponding to \a m_patches and to
    /// the number of solution fields that are to be computed
    /// m_bases[i]: The multi-basis for unknown i
    std::vector< gsMultiBasis<T> > m_bases;
    
    /// The Dof mapper is used to map patch-local DoFs to the global DoFs
    /// One for each unknown, one for each patch
    /// m_dofMappers[i]: DoF Mapper for unknown i
    std::vector<gsDofMapper>  m_dofMappers;
    
    /// Reference Quadrature rule
    gsQuadRule<T> QuRule;

    // *** Convenience members - not used by the gsAssemblerBase interface ***

    /// Dirichlet DoF fixed values (if applicable)
    gsMatrix<T> m_ddof; //-- not used in gsAssemblerBase

    /// Boundary conditions
    gsBoundaryConditions<T> m_bConditions; //-- not used in gsAssemblerBase

    /// Options
    gsAssemblerOptions m_options; //-- not used in gsAssemblerBase


    // *** Outputs ***
    

    /// Global matrix
    gsSparseMatrix<T> m_matrix;

    /// Right-hand side ( multiple right hand sides possible )
    gsMatrix<T>       m_rhs;


    // *** Information ***


    /// number of degrees of freedom (excluding eliminated etc)
    // to do: rename to m_matSize
    int m_dofs;

};


} // namespace gismo

