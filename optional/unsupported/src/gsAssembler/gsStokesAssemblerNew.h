#ifndef GSSTOKESASSEMBLERNEW_H
#define GSSTOKESASSEMBLERNEW_H


#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsPde/gsStokesPde.h>

#include <gsAssembler/gsVisitorStokes.h> // Stiffness volume integrals
//#include <gsAssembler/gsVisitorNeumannMixed.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNeumann.h>

namespace gismo
{

/** @brief
    Implementation of an Stokes solver.

    The Stokes problem: \f{eqnarray*}{
                        -\nu\Delta\mathbf{u} - \nabla p &=&\mathbf{f} \quad \text{in} \quad \Omega,\\
                        \nabla \cdot \mathbf{u} &=& 0 \quad \text{in} \quad \Omega,\\
                        \mathbf{u} &=& \mathbf{g} \quad \text{on} \quad \Gamma_D,\\
                        \nu\nabla\mathbf{u}\cdot\mathbf{n} + p \mathbf{n} &=& \mathbf{h} \quad \text{on} \quad \Gamma_N.
                        \f}
    The Neumann condition \f$ \mathbf{h}\f$ is set by applying Neumann for the velocity.
    \note The sign pressure of the pressure is fliped to obtain a symmetric linear system!
    TODO add:
            the spaces for the weak formulation we using this solver
            inf-sup warning


    See gsStokesAssembler2 for general description of input.
    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bconditions is a gsBoundaryConditions object.
    \param[in] bases is a vector of gsBasis for each component
    of the unknown of \f$\mathbf{u}\f$.
    \param[in] rhs is the right-hand side (source term) of the Stokes problem,
    \f$\mathbf{f}\f$.

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::iFace::strategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).

*/
template <class T>
class gsStokesAssemblerNew : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:

    gsStokesAssemblerNew() {}

    /** @brief Main Constructor of the assembler object.

    \param[in] pde A boundary value Stokes problem
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
    */
    gsStokesAssemblerNew( const gsStokesPde<T>          & pde,
                          const std::vector<gsMultiBasis<T> > & bases,
                          dirichlet::strategy           dirStrategy,
                          //element type!
                          iFace::strategy               intStrategy = iFace::glue)
    {
        m_options = gsAssembler<>::defaultOptions();
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);
        m_options.setInt("DirichletValues"  , dirichlet::interpolation);
        
        m_elementType = 0; // fixme : move away from here
        Base::initialize(pde, bases, m_options);
    }

    virtual gsAssembler<T>* clone() const
    {
        return new gsStokesAssemblerNew<T>(*this);;
    }
    virtual gsAssembler<T>* create() const
    {
        return new gsStokesAssemblerNew<T>;
    }

public:
    // Refresh routine
    virtual void refresh()
    {
        // 1. Common v-basis ?--> copy mappers, no=RT
        // 2. Div conforming ? --> no=block diagonal


        //space e_1 B_i, e_2 B_j, e_3 B_k
        const gsStokesPde<T>* pde_ptr = static_cast<const gsStokesPde<T>*>(m_pde_ptr.get());
        const index_t tarDim = pde_ptr->rhs()->targetDim(); //Target dimension for velocity
        std::vector<gsDofMapper> mappers;
        gsVector<index_t> rowInd = gsVector<index_t>::Zero(tarDim+1);
        gsVector<index_t> colvar= gsVector<index_t>::Zero(tarDim+1);

        switch (m_elementType)
        {
        case 0: // eg. Taylor-Hood

            //One for velocity and one for pressure
            mappers.resize(2);

            // Conforming interfaces mapper for velocity components with BCs
            if ( m_options.getInt("DirichletStrategy") == dirichlet::elimination)
                m_bases.front().getMapper(true, m_pde_ptr->bc(), mappers.front() );
            else
                m_bases.front().getMapper(true, mappers.front() );

            // Conforming interfaces mapper for pressure
            m_bases.back().getMapper(true, mappers.back() );

            rowInd(tarDim)=1; //block 0 to block tarDim uses  dofMapper 0, block tarDim uses dofMapper 1

            colvar(tarDim)=1;; //same usage for basis functions as for dofMappers
            break;

        case 1: // eg. Raviart-Thomas

            // distinct basis for each velocity component
            mappers.resize(tarDim+1);

            // Conforming interfaces mapper for each velocity component with BCs
            if ( m_options.getInt("DirichletStrategy") == dirichlet::elimination)
                for ( index_t i = 0; i!= tarDim; ++i)
                    //Note: "-1": common BC's for all components of velocity
                    m_bases[i].getMapper(true, m_pde_ptr->bc(), -1, mappers[i] );
            else
                for ( index_t i = 0; i!= tarDim; ++i)
                    m_bases[i].getMapper(true, mappers[i] );

            // Conforming interfaces mapper for pressure
            m_bases.back().getMapper(true, mappers.back() );

            for ( index_t i = 0; i!= tarDim; ++i)
                rowInd(i) = i;
            rowInd(tarDim)=tarDim; //the i-th block uses the i-th dofMapper
            for ( index_t i = 0; i!= tarDim; ++i)
                colvar(i) = i;
            colvar(tarDim)=tarDim; //same usage for basis functions as for dofMappers

            break;

        default :
            GISMO_ERROR("Error while initializing for element "<< m_elementType);
            break;
        }

        m_system = gsSparseSystem<T>(mappers,rowInd,rowInd,colvar); //rows and colums uses the same mappers
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,.333); //TODO: Improve
        m_system.reserve(nz, this->pde().numRhs());

        m_ddof.resize(tarDim+1);

    }

    // Main assembly routine
    virtual void assemble()
    {
        GISMO_ASSERT(m_system.initialized(),
                     "Sparse system is not initialized, call initialize() or refresh()");

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        const gsStokesPde<T>* pde_ptr = static_cast<const gsStokesPde<T>*>(m_pde_ptr.get());
        const index_t tarDim = pde_ptr->rhs()->targetDim();
        for ( index_t i = 0; i!= tarDim+1; ++i)
            Base::computeDirichletDofs(i);

        // Clean the sparse system
        m_system.setZero();



        // Assemble volume integrals
        switch (m_elementType)
        {
        case 0: // Taylor-Hood
            Base::template push<gsVisitorStokes<T> >();
            break;
        case 1: // Div-conforming
            //assembleDC();
            break;
        default :
            GISMO_ERROR("Something went terribly wrong.");
            break;
        }

        // Enforce Neumann boundary conditions
        Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );

        // If requested, enforce Dirichlet boundary conditions by Nitsche's method
        // if ( m_options.dirStrategy == dirichlet::nitsche ) //TODO: not supported
        //    Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());

        // If requested, enforce Dirichlet boundary conditions by diagonal penalization
        if ( m_options.getInt("DirichletStrategy") == dirichlet::penalize )
            Base::penalizeDirichletDofs();


        // Assembly is done, compress the matrix
        Base::finalize();
    }



    void assembleNitsche()
    {
    }

protected:
    unsigned m_elementType;

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

};

}
#endif // GSSTOKESASSEMBLERNEW_H
