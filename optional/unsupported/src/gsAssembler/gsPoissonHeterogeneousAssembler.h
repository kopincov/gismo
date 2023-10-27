
/** @file gsPoissonHeterogeneousAssembler.h

    @brief Provides assembler for the Poisson equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#pragma once

#include <gsIO/gsOptionList.h>
#include <gsPde/gsPoissonHeterogeneousPde.h>
#include <gsAssembler/gsVisitorPoissonHeterogeneous.h> // Stiffness volume integrals and load vector
#include <gsAssembler/gsVisitorDg2.h>
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals
#include <gsAssembler/gsVisitorNitsche2.h> // Nitsche boundary integrals
#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsPoissonAssembler.h>
//#include <gsAssembler/gsAssemblerUtils.h>


namespace gismo
{

/** @brief
    Implementation of an (multiple right-hand side) Poisson solver.

    The Poisson equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).
*/
template <class T>
class gsPoissonHeterogeneousAssembler : public gsPoissonAssembler<T>
{
public:
    typedef gsPoissonAssembler<T> Base;

public:
    /// Default Constructor
    gsPoissonHeterogeneousAssembler() {}

public:

    gsPoissonHeterogeneousAssembler( const gsPoissonHeterogeneousPde<T> & pde,
                                     const gsMultiBasis<T>          & bases,
                                     dirichlet::strategy           dirStrategy,
                                     iFace::strategy               intStrategy = iFace::glue)
    {
        m_options.setInt("DirichletStrategy" , dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);
        m_mixedIFaceConditions.resize(pde.domain().nInterfaces(),intStrategy);

        gsAssembler<T>::initialize(pde, bases, m_options);
    }

    gsPoissonHeterogeneousAssembler( const gsPoissonHeterogeneousPde<T> & pde,
                                     const gsMultiBasis<T>          & bases,
                                     dirichlet::strategy           dirStrategy,
                                     std::vector<iFace::strategy>  &  intStrategy)
        : m_mixedIFaceConditions(intStrategy)
    {
        GISMO_ASSERT(intStrategy.size() == (size_t)pde.domain().nInterfaces(), "Number of interfaces does not match the given number of interface strategies." );

        m_options.setInt("DirichletStrategy" , dirStrategy);
        m_options.setInt("InterfaceStrategy", -1);

        gsAssembler<T>::initialize(pde, bases, m_options);
    }

    virtual gsAssembler<T>* clone() const
    {
        return new gsPoissonHeterogeneousAssembler<T>(*this);
    }

    virtual gsAssembler<T>* create() const
    {
        return new gsPoissonHeterogeneousAssembler<T>;
    }

public:
    virtual void refresh()
    {
        if(m_options.getInt("InterfaceStrategy")!=-1)
            Base::scalarProblemGalerkinRefresh();
        else
        {
        // Check for coherency
        GISMO_ASSERT(this->check(), "Incoherent data in Assembler");

        GISMO_ASSERT(1==m_bases.size(), "Expecting a single discrete space "
                                        "for standard scalar Galerkin");

        // 1. Obtain a map from basis functions to matrix columns and rows
        gsDofMapper mapper;
        m_bases.front().getMapper(
            (dirichlet::strategy)(m_options.getInt("DirichletStrategy")),
            (iFace::dg),
            this->pde().bc(), mapper, 0,false);

        int ii=0;
        for(typename gsMultiPatch<T>::iiterator it= m_pde_ptr->domain().iBegin();it!=m_pde_ptr->domain().iEnd();++it)
        {
            boundaryInterface bI = *it;
            if(m_mixedIFaceConditions[ii] == iFace::conforming)
                mapper.matchDofs(bI.first().patch,m_bases.front()[bI.first().patch].boundary(bI.first().side()),
                        bI.second().patch,m_bases.front()[bI.second().patch].boundary(bI.second().side()));
            ii++;
        }
        mapper.finalize();
        if ( 0 == mapper.freeSize() ) // Are there any interior dofs ?
            gsWarn << " No internal DOFs, zero sized system.\n";

        // 2. Create the sparse system
        m_system = gsSparseSystem<T>(mapper);//1,1
        }
    }

    /// Main assembly routine.
    virtual void assemble()
    {
        GISMO_ASSERT(m_system.initialized(),
                     "Sparse system is not initialized, call initialize() or refresh()");

        m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        Base::computeDirichletDofs();

        // Clean and reserve the sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.reserve(nz,m_pde_ptr->numRhs());

        // Assemble volume integrals
        Base::template push<gsVisitorPoissonHeterogeneous<T> >();

        // Enforce Neumann boundary conditions
        Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );

        // If requested, enforce Dirichlet boundary conditions by Nitsche's method
        if ( m_options.getInt("DirichletStrategy") == dirichlet::nitsche )
            Base::template push<gsVisitorNitsche2<T> >(m_pde_ptr->bc().dirichletSides());
        // If requested, enforce Dirichlet boundary conditions by diagonal penalization
        else if ( m_options.getInt("DirichletStrategy") == dirichlet::penalize )
            Base::penalizeDirichletDofs();

        // If we are in in dg (Discontinuous Galerkin) mode: add
        // interface contributions
        if ( m_options.getInt("InterfaceStrategy") == iFace::dg)
            push<gsVisitorDg2<T,1> >(m_pde_ptr->patches().interfaces());
        else if(m_options.getInt("InterfaceStrategy") == -1 )
        {
            const std::vector<boundaryInterface>& allI =  m_pde_ptr->patches().interfaces();
            std::vector<boundaryInterface> interfaces;
            for(size_t i=0; i<allI.size();++i)
                if(m_mixedIFaceConditions[i]==iFace::dg)
                    interfaces.push_back(allI[i]);
            push<gsVisitorDg2<T,1> >(interfaces);
        }
        // Assembly is done, compress the matrix
        Base::finalize();

        //  gsInfo<<"the matrix starts here: \n"<<matrices[np].toDense()<<"\n\n";

        /*
        int max=0;
        for(int i=0; i<m_system.matrix().cols();++i)
        {
            if(max<m_system.matrix().col(i).nonZeros())
                max = m_system.matrix().col(i).nonZeros();

        }
        gsInfo<<"Max number of entries: "<<max <<"\n"<<"usage: "<<(real_t)max/this->options().numColNz(m_bases[0][0])<<"\n";
        //*/
    }


    virtual T penalty(int k) const
    {
        const int deg = m_bases[0][k].maxDegree();
        return (deg + m_bases[0][k].dim()) * (deg + 1) * T(2.0);
    }


    /// @brief Returns the used-dG Visitor (should be an abstract class)
    virtual gsVisitorDg2<T,1> visitorDg(const boundaryInterface& bi)
    {
        //TODO: make this more abstract (for IETI)

        //Needs to be reworked when dG is properly incorporated
        return gsVisitorDg2<T,1>(*m_pde_ptr, bi.first(),bi.second());
        //gsVisitorDg2<T>(penalty(bI->first().patch),bi,ppde->getAlpha());
    }

    std::vector<iFace::strategy> getMixedInterfaceConditions() const
    {
        return m_mixedIFaceConditions;
    }
protected:
    // DG contributions
    /*
    void assembleDg()
    {
        for ( typename gsMultiPatch<T>::iiterator
              it = m_patches.iBegin(); it != m_patches.iEnd(); ++it )
        {
            const boundaryInterface & iFace =
                    ( m_bases[0][it->first() .patch].numElements(it->first() .side() ) <
                    m_bases[0][it->second().patch].numElements(it->second().side() ) ?
                        it->getInverse() : *it );

            const T pp = penalty( iFace.first().patch );

            gsVisitorDg2<T> dg(pp, iFace,m_alpha,1);

            this->apply(dg, iFace);
        }
    }
*/

public:
    // TODO: move this at some point in the future to gsAssembler
    template<class DGVisitor>
    void push(const std::vector<boundaryInterface> & bIs)
    {

        for(typename gsMultiPatch<T>::const_iiterator  bI= bIs.begin(); bI!= bIs.end();++bI)
        {
            const boundaryInterface & iFace =
                    ( m_bases[0][bI->first() .patch].numElements(bI->first() .side() ) <
                    m_bases[0][bI->second().patch].numElements(bI->second().side() ) ?
                        bI->getInverse() : *bI );
            DGVisitor visitor(*m_pde_ptr, iFace.first(), iFace.second());

            const gsAffineFunction<T> interfaceMap(m_pde_ptr->domain().getMapForInterface(iFace));

            const index_t patch1      = iFace.first().patch;
            const index_t patch2      = iFace.second().patch;
            const gsBasis<T> & B1 = m_bases[0][patch1];// (!) unknown 0
            const gsBasis<T> & B2 = m_bases[0][patch2];

            gsQuadRule<T> QuRule ; // Quadrature rule
            gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
            gsVector<T> quWeights;         // Mapped weights
            // Evaluation flags for the Geometry map
            unsigned evFlags(0);

            const int bSize1      = B1.numElements( iFace.first() .side() );
            const int bSize2      = B2.numElements( iFace.second().side() );
            const int ratio = bSize1 / bSize2;
            GISMO_ASSERT(bSize1 >= bSize2 && bSize1%bSize2==0,
                         "DG assumes nested interfaces. Got bSize1="<<
                         bSize1<<", bSize2="<<bSize2<<"." );

            // Initialize
            visitor.initialize(B1, B2, QuRule, evFlags);

            // Initialize geometry evaluators
            typename gsGeometryEvaluator<T>::uPtr geoEval1(getEvaluator(evFlags, m_pde_ptr->domain()[patch1]));
            typename gsGeometryEvaluator<T>::uPtr geoEval2(getEvaluator(evFlags, m_pde_ptr->domain()[patch2]));

            // Initialize domain element iterators
            typename gsBasis<T>::domainIter domIt1 = B1.makeDomainIterator( iFace.first() .side() );
            typename gsBasis<T>::domainIter domIt2 = B2.makeDomainIterator( iFace.second().side() );


            int count = 0;
            // iterate over all boundary grid cells on the "left"
            for (; domIt1->good(); domIt1->next() )
            {
                count++;
                // Get the element of the other side in domIter2
                //domIter1->adjacent( iFace.orient, *domIter2 );

                // Compute the quadrature rule on both sides
                QuRule.mapTo( domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights);
                interfaceMap.eval_into(quNodes1,quNodes2);

                // Perform required evaluations on the quadrature nodes
                visitor.evaluate(B1, *geoEval1, B2, *geoEval2, quNodes1, quNodes2);

                // Assemble on element
                visitor.assemble(*domIt1,*domIt2, *geoEval1, *geoEval2, quWeights);

                // Push to global patch matrix (m_rhs is filled in place)

                //this just works for poisson (second one is better) -- FIXME
                visitor.localToGlobal(m_system.dofMappers()[0],m_ddof[0], patch1, patch2, m_system.matrix(), m_system.rhs());

                //visitor.localToGlobal(patch1, patch2, m_ddof, m_system); // USE THIS: TODO: implement in SparseSystem

                if ( count % ratio == 0 ) // next master element ?
                {
                    domIt2->next();
                }

            }
        }
    }

protected:
    std::vector<iFace::strategy> m_mixedIFaceConditions;
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};



//////////////////////////////////////////////////
//////////////////////////////////////////////////



} // namespace gismo
