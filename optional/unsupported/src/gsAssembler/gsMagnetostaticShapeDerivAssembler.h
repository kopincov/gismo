
/** @file gsMagnetostaticShapeDerivAssembler.h

    @brief Provides assembler for the (now: linear) shape derivate of the magnetostatic problem

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
#include <gsPde/gsMagnetostaticShapeDerivPde.h>

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
class gsMagnetostaticShapeDerivAssembler : public gsAssembler<T>
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

    gsMagnetostaticShapeDerivAssembler( const gsMagnetostaticShapeDerivPde<T> & pde, const gsMultiBasis<T> & bases,
                              const gsOptionList opt): m_domain(pde.domain())
    {
        m_mixedIFaceConditions.resize(pde.domain().nInterfaces(), iFace::conforming);
        Base::initialize(pde, bases, opt);
    }

    virtual gsAssembler<T>* clone() const
    {
        return new gsMagnetostaticShapeDerivAssembler<T>(*this);
    }


    virtual gsAssembler<T>* create() const
    {
        //this should be a default constructor, however with const refs this is not possible
        return new gsMagnetostaticShapeDerivAssembler<T>(*this);
    }

    const gsMatrix<T> & getRhsFull()
    {
        return m_system.rhs();
    }

    virtual void refresh()
    {
        // Check for coherency
        GISMO_ASSERT(this->check(), "Incoherent data in Assembler");

        //GISMO_ASSERT(2==m_bases.size(), "Expecting a single discrete space "
        //                                "for standard scalar Galerkin");

        //phi.setInterfaceCont(0);
        //phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));

        gsExprAssembler<T> exprAss(1,1);
        exprAss.setIntegrationElements(m_bases.front());

        typedef gsExprAssembler<real_t>::space space;

        space phi = exprAss.getSpace(m_bases.front(), 2);
        phi.setInterfaceCont(0);
        phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));

        exprAss.setOptions(m_options);

        exprAss.initSystem();

        // read the Dirichlet dofs and the mapper
        space u = exprAss.testSpace(0);

        // 2. Create the sparse system
        //m_system = gsSparseSystem<T>(mapper);//1,1
        m_ddof.resize(2);
        gsDofMapper mapper = u.mapper();

        std::vector<gsDofMapper> mappers;
        gsVector<index_t> dims(2);

        dims(0) = 1;
        dims(1) = 1;

        for(size_t d = 0; d < 2; d++) {
            mappers.push_back(mapper);
        }
        if(!mapper.isFinalized())
            mapper.finalize();

        m_system = gsSparseSystem<T>(mappers, dims);

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

        //Define the unit functions
        gsFunctionExpr<> unit1("1.0", "0.0", 2);
        gsFunctionExpr<> unit2("0.0", "1.0", 2);

        gsExprAssembler<T> exprAss(1,1);
        exprAss.setOptions(m_options);
        exprAss.setIntegrationElements(m_bases.front());

        gsMagnetostaticShapeDerivPde<T>* ppde = dynamic_cast<gsMagnetostaticShapeDerivPde<T>* >(m_pde_ptr.get());
        gsMultiPatch<T> domain = ppde->domain();

        space phi = exprAss.getSpace(m_bases.front(), 2);
        geometryMap G = exprAss.getMap(domain);
        variable scale = exprAss.getCoeff(*(ppde->getScaling()));
        variable matCoeff = exprAss.getCoeff(*(ppde->getAlpha()), G);
        variable magnetization = exprAss.getCoeff(*(ppde->getRhs()), G); // method rhs() gives only the 0-th piece of the rhs; works well for IETI
        variable uh = exprAss.getCoeff(*(ppde->getStateSolution()));
        variable ph = exprAss.getCoeff(*(ppde->getAdjointSolution()));

        phi.setInterfaceCont(0);
        phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));

        exprAss.initSystem();

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        //Base::computeDirichletDofs();

        // Clean and reserve the sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.reserve(nz,m_pde_ptr->numRhs());

        // Set the required data
        // Set the Ansatz space
      //  m_ExprAss.setOptions(m_options);
        //m_ExprAss.setIntegrationElements(m_bases[0]);

        // Initialize the system
     //   m_ExprAss.initSystem();

        // Assemble volume integrals (lhs, rhs)
        //Auxiliary problem for computing the shape derivative

        // rhs expression
        /*auto rhs_expr = ijac(phi, G) % (id(2) * (matCoeff.val() * (fjac(uh).tr() * jac(G).ginv() * (fjac(ph).tr() * jac(G).ginv()).tr())).val() ) +
                        ijac(phi, G) % (id(2) * ((-fjac(ph).tr()) * jac(G).ginv() * magnetization / 1.086).val() ) +
                        ijac(phi, G) % ((1/1.086) * magnetization * (fjac(ph).tr() * jac(G).ginv())).tr() +
                        ijac(phi, G) % ((-matCoeff.val()) * ( ((jac(G).ginv()).tr() * fjac(ph)) * (fjac(uh).tr() * jac(G).ginv()) )) +
                        ijac(phi, G) % ((-matCoeff.val()) * ( ((jac(G).ginv()).tr() * fjac(uh)) * (fjac(ph).tr() * jac(G).ginv()) ));
        */ // auto isn't allowed in C++98

        // Assemble the lhs and rhs, lhs depends on the choice of the Hilbert space
        exprAss.assemble( ((phi * phi.tr()) + (scale.val() * (ijac(phi, G) % ijac(phi, G))) ) * meas(G),
                          ijac(phi, G) % (id(2) * (matCoeff.val() * (fjac(uh).tr() * jac(G).ginv() * (fjac(ph).tr() * jac(G).ginv()).tr())).val() ) +
                              ijac(phi, G) % (id(2) * ((-fjac(ph).tr()) * jac(G).ginv() * magnetization / 1.086).val() ) +
                              ijac(phi, G) % ((1/1.086) * magnetization * (fjac(ph).tr() * jac(G).ginv())).tr() +
                              ijac(phi, G) % ((-matCoeff.val()) * ( ((jac(G).ginv()).tr() * fjac(ph)) * (fjac(uh).tr() * jac(G).ginv()) )) +
                              ijac(phi, G) % ((-matCoeff.val()) * ( ((jac(G).ginv()).tr() * fjac(uh)) * (fjac(ph).tr() * jac(G).ginv()) ))
                                  * meas(G) );

        // Test PDE to test the gsMagnetostaticAssembler against the gsPoissonHeterogeneousAssembler
        //exprAss.assemble( igrad(phi, G) * igrad(phi, G).tr() * matCoeff.val() * meas(G), phi * rhs * meas(G) );

        m_system.matrix() = exprAss.matrix();
        m_system.rhs() = exprAss.rhs();

        //gsInfo<<"mat: "<<m_system.matrix().toDense()<<"\n\n";

        // read the Dirichlet dofs and the mapper
        space u = exprAss.testSpace(0);
        m_ddof[0]= u.fixedPart().col(0);
        m_ddof[1]= u.fixedPart().col(1);
        /*
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
        */

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

    /// Main assembly routine.
    virtual void constructSolution(gsMatrix<T> & solVector, gsMultiPatch<T> & shapeGradient)
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

        gsMagnetostaticShapeDerivPde<T>* ppde = dynamic_cast<gsMagnetostaticShapeDerivPde<T>* >(m_pde_ptr.get());

        space phi = exprAss.getSpace(m_bases.front(), 2);
        geometryMap G = exprAss.getMap(ppde->domain());
        
        phi.setInterfaceCont(0);
        phi.addBc((m_pde_ptr->bc()).get("Dirichlet"));

        exprAss.initSystem();

        solution uGrad = exprAss.getSolution(phi, solVector);
        uGrad.extract( shapeGradient );

    }


  

protected:
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

    gsMultiPatch<T> m_domain;

    std::vector<iFace::strategy> m_mixedIFaceConditions;
    
};



//////////////////////////////////////////////////
//////////////////////////////////////////////////



} // namespace gismo
