
/** @file gsMagnetostaticAssembler.h

    @brief Provides assembler for the (now: linear) magnetostatic equation. Combines the expr. assembler and the heterogeneous Poisson assembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Schneckenleitner
*/

#pragma once

#include <fstream>
#include <gsIO/gsOptionList.h>
#include <gsPde/gsMagnetostaticPde.h>
#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExpressions.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>
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
class gsMagnetostaticAssembler : public gsAssembler<T>
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

    gsMagnetostaticAssembler( const gsMagnetostaticPde<T> & pde, const gsMultiBasis<T> & bases, const gsOptionList opt, const gsPiecewiseFunction<> & f)
    {
        m_mixedIFaceConditions.resize(pde.domain().nInterfaces(), iFace::conforming);
        m_rhs = f;
        Base::initialize(pde, bases, opt);
    }

    virtual gsAssembler<T>* clone() const
    {
        return new gsMagnetostaticAssembler<T>(*this);
    }


    virtual gsAssembler<T>* create() const
    {
        //this should be a default constructor, however with const refs this is not possible
        return new gsMagnetostaticAssembler<T>(*this);
    }


public:
    const gsPiecewiseFunction<T>* getAlpha()
    {
        return m_alpha;
    }

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

        //phi.setInterfaceCont(0);
        //phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));

        gsExprAssembler<T> exprAss(1,1);
        exprAss.setIntegrationElements(m_bases.front());

        gsMagnetostaticPde<T>* ppde = dynamic_cast<gsMagnetostaticPde<T>* >(m_pde_ptr.get());

        typedef gsExprAssembler<real_t>::space space;
        typedef gsExprAssembler<real_t>::geometryMap geometryMap; // Need geometryMap for spacetime

        space phi = exprAss.getSpace(m_bases.front());
        phi.setInterfaceCont(0);
        phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));
        exprAss.getMap(ppde->domain());

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
        // Write data of interest to a file
        //std::ofstream file;
        //file.open("localRhs.txt", std::ios_base::app);

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

        gsMagnetostaticPde<T>* ppde = dynamic_cast<gsMagnetostaticPde<T>* >(m_pde_ptr.get());
        gsMultiPatch<T> dom = ppde->domain();


        space phi = exprAss.getSpace(m_bases.front());
        geometryMap G = exprAss.getMap(ppde->domain());
        variable matCoeff = exprAss.getCoeff(*(ppde->getAlpha()), G);
        variable rhs = exprAss.getCoeff(*(ppde->rhs()), G); // method rhs() gives only the 0-th piece of the rhs; works well for IETI
        
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
        exprAss.assemble( igrad(phi, G) * igrad(phi, G).tr() * matCoeff.val() * meas(G), igrad(phi, G) * rhs * meas(G) );

        // Test PDE to test the gsMagnetostaticAssembler against the gsPoissonHeterogeneousAssembler
        //exprAss.assemble( igrad(phi, G) * igrad(phi, G).tr() * matCoeff.val() * meas(G), phi * rhs * meas(G) );
        
        m_system.matrix() = exprAss.matrix();

        //if(file.is_open())
            //file << "Local matrix: " << dom[0].id() << "\n" << exprAss.matrix().toDense() << "\n";
            //file << "Local rhs: " << dom[0].id() << "\n" << exprAss.rhs() << "\n\n";

        //gsInfo << "Local matrix: " << dom[0].id() << "\n" << exprAss.matrix().toDense() << "\n";

        m_system.rhs() = exprAss.rhs();
        //gsInfo<<"mat: "<<m_system.matrix().toDense()<<"\n\n";

        // read the Dirichlet dofs and the mapper
        space u = exprAss.testSpace(0);
        m_ddof[0]= u.fixedPart();
        
        // Assembly is done, compress the matrix
        Base::finalize();

        // Close the file
        //file.close();

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

        gsMagnetostaticPde<T>* ppde = dynamic_cast<gsMagnetostaticPde<T>* >(m_pde_ptr.get());

        space phi = exprAss.getSpace(m_bases.front());
        geometryMap G = exprAss.getMap(ppde->domain());
        variable matCoeff = exprAss.getCoeff(*(ppde->getAlpha()), G);
        variable rhs = exprAss.getCoeff(m_rhs, G);
        
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
        exprAss.assemble( igrad(phi, G) * igrad(phi, G).tr() * matCoeff.val() * meas(G), igrad(phi, G) * rhs * meas(G) );

        m_system.matrix() = exprAss.matrix();
        m_system.rhs() = exprAss.rhs();

        // read the Dirichlet dofs and the mapper
        space u = exprAss.testSpace(0);
        m_ddof[0]= u.fixedPart();
        
        // Assembly is done, compress the matrix
        Base::finalize();
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

        gsMagnetostaticPde<T>* ppde = dynamic_cast<gsMagnetostaticPde<T>* >(m_pde_ptr.get());

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

    std::vector<iFace::strategy> m_mixedIFaceConditions;
    gsPiecewiseFunction<T> m_rhs;
    const gsPiecewiseFunction<T>* m_alpha;
    
};



//////////////////////////////////////////////////
//////////////////////////////////////////////////



} // namespace gismo
