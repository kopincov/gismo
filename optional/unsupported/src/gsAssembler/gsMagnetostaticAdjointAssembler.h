
/** @file gsMagnetostaticAdjointAssembler.h

    @brief Provides assembler for the (now: linear) magnetostatic equation. Combines the expr. assembler and the heterogeneous Poisson assembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Schneckenleitner
*/

#pragma once

#include <gsIO/gsOptionList.h>
#include <gsPde/gsPoissonHeterogeneousPde.h>
#include <gsAssembler/gsVisitorPoissonHeterogeneous.h> // Stiffness volume integrals and load vector
#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsAssembler/gsExpressions.h>
#include <gsPde/gsMagnetostaticAdjointPde.h>

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
class gsMagnetostaticAdjointAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

    typedef gsExprAssembler<real_t>::space space;
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::variable variable;
    typedef gsExprAssembler<real_t>::solution solution;

//public:
    /// Default Constructor
    //gsMagnetostaticAssembler(): matCoeff(NULL), rhs(NULL), phi(NULL), G(NULL) {}

public:

    gsMagnetostaticAdjointAssembler( const gsMagnetostaticAdjointPde<T> & pde, const gsMultiBasis<T> & bases,
                              const gsOptionList opt, const std::vector<std::pair<unsigned, unsigned> > & interfaces):
            m_interfaces(interfaces), m_domain(pde.domain())
    {
        m_mixedIFaceConditions.resize(pde.domain().nInterfaces(), iFace::conforming);
        m_domain = pde.domain();
        Base::initialize(pde, bases, opt);
    }

    virtual gsAssembler<T>* clone() const
    {
        return new gsMagnetostaticAdjointAssembler<T>(*this);
    }


    virtual gsAssembler<T>* create() const
    {
        //this should be a default constructor, however with const refs this is not possible
        return new gsMagnetostaticAdjointAssembler<T>(*this);
    }


public:
    const gsMatrix<T> & getRhsFull()
    {
        return m_system.rhs();
    }

    virtual void refresh()
    {

        // Check for coherency
        GISMO_ASSERT(this->check(), "Incoherent data in Assembler");

        GISMO_ASSERT(1==m_bases.size(), "Expecting a single discrete space "
                                        "for standard scalar Galerkin");

        gsExprAssembler<T> exprAss(1,1);
        exprAss.setIntegrationElements(m_bases.front());


        typedef gsExprAssembler<real_t>::space space;

        space phi = exprAss.getSpace(m_bases.front());
        phi.setInterfaceCont(0);
        phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));

        exprAss.setOptions(m_options);

        exprAss.initSystem();

        // read the Dirichlet dofs and the mapper
        space u = exprAss.testSpace(0);

        // 2. Create the sparse system
        //m_system = gsSparseSystem<T>(mapper);//1,1
        m_ddof.resize(1);
        gsDofMapper mapper = u.mapper();
        if(!mapper.isFinalized())
            mapper.finalize();
        m_system = gsSparseSystem<T>(mapper);


    }

    /// Main assembly routine.
    virtual void assemble()
    {
        using namespace gismo::expr;
        
        GISMO_ASSERT(m_system.initialized(),
                     "Sparse system is not initialized, call initialize() or refresh()");

        typedef gsExprAssembler<real_t>::space space;
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::solution solution;

        gsExprAssembler<T> exprAss(1,1);
        exprAss.setOptions(m_options);
        exprAss.setIntegrationElements(m_bases.front());

        gsMagnetostaticAdjointPde<T>* ppde = dynamic_cast<gsMagnetostaticAdjointPde<T>* >(m_pde_ptr.get());
        gsMultiPatch<T> domain = ppde->domain();

        space phi = exprAss.getSpace(m_bases.front());
        geometryMap G = exprAss.getMap(domain);
        variable matCoeff = exprAss.getCoeff(*(ppde->getAlpha()), G);
        variable B_d = exprAss.getCoeff(*(ppde->getReference()), G); // method rhs() gives only the 0-th piece of the rhs; works well for IETI
        variable uh = exprAss.getCoeff(*(ppde->getStateSolution()));

        phi.setInterfaceCont(0);
        phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));

        exprAss.initSystem();

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        //Base::computeDirichletDofs();

        // Clean and reserve the sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.reserve(nz,m_pde_ptr->numRhs());

        // Assemble volume integrals (lhs, rhs)
        //Magnetostatic PDE
        exprAss.assemble( igrad(phi, G) * igrad(phi, G).tr() * matCoeff.val() * meas(G) );

        gsMultiPatch<T> InterfaceMP(m_domain);

        if (domain.nPatches() == 1)
        {
            for(std::vector<std::pair<unsigned, unsigned> >::iterator it = m_interfaces.begin(); it != m_interfaces.end(); it++)
            {
                if((*it).first == (domain.piece(0)).id())
                {

                    std::vector<boundaryInterface> cont;
                    const boundaryInterface* cBI = m_domain.findInterface(it->first,it->second);

                    if(cBI==NULL)
                        GISMO_ERROR("Interface "+util::to_string(it->first)+" - "+util::to_string(it->second)+" not valid.");

                    // A new boundaryInterface must be created in order to have the multipatch to be
                    boundaryInterface bI(patchSide(0,cBI->first().side()),patchSide(cBI->second().patch,cBI->second().side()),
                                         static_cast<short_t>(2));

                    cont.push_back(bI);
                    //gsInfo << "This is the index: " << composite << "\n";
                    //const boundaryInterface* sides = m_domain.findInterface((*it).first, (*it).second);
                    //cont.push_back(*sides);
                    //InterfaceMP.addInterface(&InterfaceMP.patch((*it).first), sides->first(), &InterfaceMP.patch((*it).second), sides->second());
                    exprAss.assembleRhsInterface( (igrad(phi, G) * (tv(G)/tv(G).norm())) * (-2 * (fjac(uh).tr() * (jac(G).ginv()) * (tv(G)/tv(G).norm()) - B_d)) * nv(G).norm(), cont );
                }
            }
        }
        else
        {
            createInterfaceTopology(InterfaceMP);
            exprAss.assembleRhsInterface( (igrad(phi, G) * (tv(G)/tv(G).norm())) * (-2 * (fjac(uh).tr() * (jac(G).ginv()) * (tv(G)/tv(G).norm()) - B_d)) * nv(G).norm(), InterfaceMP.interfaces() );
        }

        m_system.matrix() = exprAss.matrix();
        m_system.rhs() = exprAss.rhs();

        //gsInfo<<"mat: "<<m_system.matrix().toDense()<<"\n\n";

        // read the Dirichlet dofs and the mapper
        space u = exprAss.testSpace(0);
        m_ddof[0]= u.fixedPart();

        
        // Assembly is done, compress the matrix
        Base::finalize();
    }
    
    /// Main assembly routine.
    virtual void assembleFull()
    {
        using namespace gismo::expr;
        
        GISMO_ASSERT(m_system.initialized(),
                     "Sparse system is not initialized, call initialize() or refresh()");

        typedef gsExprAssembler<real_t>::space space;
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::solution solution;

        gsExprAssembler<T> exprAss(1,1);
        exprAss.setOptions(m_options);
        exprAss.setIntegrationElements(m_bases.front());

        gsMagnetostaticAdjointPde<T>* ppde = dynamic_cast<gsMagnetostaticAdjointPde<T>* >(m_pde_ptr.get());

        space phi = exprAss.getSpace(m_bases.front());
        geometryMap G = exprAss.getMap(ppde->domain());
        variable matCoeff = exprAss.getCoeff(*(ppde->getAlpha()), G);
        variable B_d = exprAss.getCoeff(*(ppde->getReference()), G); // method rhs() gives only the 0-th piece of the rhs; works well for IETI
        variable uh = exprAss.getCoeff(*(ppde->getStateSolution()));
        
        phi.setInterfaceCont(0);
        phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));

        exprAss.initSystem();

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        //Base::computeDirichletDofs();

        // Clean and reserve the sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.reserve(nz,m_pde_ptr->numRhs());

        // Assemble volume integrals (lhs, rhs)
        gsInfo<<"Before Assembling all!\n";
        //Magnetostatic PDE
        gsMultiPatch<T> InterfaceMP(m_domain);
        createInterfaceTopology(InterfaceMP);
        exprAss.assemble( igrad(phi, G) * igrad(phi, G).tr() * matCoeff.val() * meas(G) );
        exprAss.assembleRhsInterface( (igrad(phi, G) * (tv(G)/tv(G).norm())) * (-2 * (fjac(uh).tr() * (jac(G).ginv()) * (tv(G)/tv(G).norm()) - B_d)) * nv(G).norm(), InterfaceMP.interfaces() );

        m_system.matrix() = exprAss.matrix();
        m_system.rhs() = exprAss.rhs();

        // read the Dirichlet dofs and the mapper
        space u = exprAss.testSpace(0);
        m_ddof[0]= u.fixedPart();
        
        // Assembly is done, compress the matrix
        Base::finalize();
    }

    // Helper to create the multipatch domain only with the interfaces
    void createInterfaceTopology(gsMultiPatch<T> & mp)
    {
        mp.clearTopology();

        for(std::vector<std::pair<unsigned, unsigned> >::iterator it = m_interfaces.begin(); it != m_interfaces.end(); it++)
        {
            std::vector<boundaryInterface> cont;
            const boundaryInterface* cBI = m_domain.findInterface(it->first,it->second);

            if(cBI==NULL)
                GISMO_ERROR("Interface "+util::to_string(it->first)+" - "+util::to_string(it->second)+" not valid.");

            mp.addInterface(*cBI);
        }
    }

    /// Main assembly routine.
    virtual void constructSolution(gsMatrix<T> & solVector, gsMultiPatch<T> & mp)
    {
        using namespace gismo::expr;

        GISMO_ASSERT(m_system.initialized(),
                     "Sparse system is not initialized, call initialize() or refresh()");

        typedef gsExprAssembler<real_t>::space space;
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::variable variable;
        typedef gsExprAssembler<real_t>::solution solution;

        gsExprAssembler<T> exprAss(1,1);
        exprAss.setOptions(m_options);
        exprAss.setIntegrationElements(m_bases.front());

        gsPoissonHeterogeneousPde<T>* ppde = dynamic_cast<gsPoissonHeterogeneousPde<T>* >(m_pde_ptr.get());

        space phi = exprAss.getSpace(m_bases.front());
        exprAss.getMap(ppde->domain());

        phi.setInterfaceCont(0);
        phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));

        exprAss.initSystem();

        solution uGrad = exprAss.getSolution(phi, solVector);
        uGrad.extract( mp );

    }


protected:
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

    std::vector<std::pair<unsigned, unsigned> > m_interfaces;
    gsMultiPatch<T> m_domain;
        /*
    gsExprAssembler<real_t> m_ExprAss;
    space phi;
    geometryMap G;
    variable rhs;
    variable matCoeff;
    //*/
    std::vector<iFace::strategy> m_mixedIFaceConditions;
    gsPiecewiseFunction<T> m_rhs;

};



//////////////////////////////////////////////////
//////////////////////////////////////////////////



} // namespace gismo
