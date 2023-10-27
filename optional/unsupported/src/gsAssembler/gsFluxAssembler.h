#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsFluxPde.h>
#include <gsAssembler/gsVisitorFlux.h>
#include <gsAssembler/gsVisitorNitsche.h>
#include <gsAssembler/gsVisitorNeumann.h>

namespace gismo
{

/** @brief
    Implementation of an (multiple right-hand side) assembler for an equation
    of the following form.

     \f$-\text{div}(F(\nabla\mathbf{u}))=\mathbf{f} \f$

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system.

    \ingroup Assembler
*/
template <class T>
class gsFluxAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:

    gsFluxAssembler()
    { }

    /** @brief Main Constructor of the assembler object.

    \param[in] pde A boundary value Poisson problem
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
    */
    gsFluxAssembler( const gsFluxPde<T>                & pde,
                        const gsMultiBasis<T>          & bases,
                        dirichlet::strategy              dirStrategy,
                        iFace::strategy                  intStrategy = iFace::glue)
    {
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        Base::initialize(pde, bases, m_options);
    }

    virtual gsAssembler<T>* clone() const
    {
        return new gsFluxAssembler<T>(*this);
    }
    virtual gsAssembler<T>* create() const
    {
        return new gsFluxAssembler<T>;
    }

    virtual void refresh()
    {
        Base::scalarProblemGalerkinRefresh();
    }

    /**
     * @brief assemble selects as initial solution the zero function and calls
     * afterwards void assemble(const gsMultiPatch<T>& curSolution).
     */
    virtual void assemble()
    {

        // Reserve sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.reserve(nz, this->pde().numRhs());
        
        Base::computeDirichletDofs();

        //Create initial value. Start with zero function.
        gsMultiPatch<T>sol(m_pde_ptr->domain());
        for(size_t np=0; np<m_pde_ptr->domain().nPatches();++np)
        {
            gsMatrix<T> & coeffs = sol.patch(np).coefs();
            coeffs.setZero(coeffs.rows(),1);
        }
        //Call the nonlinear assembling routine with this starting value
        assemble(sol);


    }

    // Main assembly routine
    virtual void assemble(const gsMultiPatch<T>& curSolution)
    {
        if(m_ddof.size()==0)
            Base::computeDirichletDofs();

        m_system.setZero();
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.reserve(nz, this->pde().numRhs());
        const gsFluxPde<T>* fluxpde = static_cast<const gsFluxPde<T>*>(m_pde_ptr.get());
        fluxpde->setCurSolution(curSolution);

        // Visit elements
        Base::template push<gsVisitorFlux<T> >();

        // Enforce Neumann boundary conditions //TODO: this must be improved
        Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );

        if ( m_options.getInt("DirichletStrategy") == dirichlet::nitsche )
            Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());
        else if ( m_options.getInt("DirichletStrategy") == dirichlet::penalize )
            Base::penalizeDirichletDofs();

        if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
            gsWarn <<"DG option is ignored.\n";

        // Assembly is done, compress the matrix
        Base::finalize();
    }

    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() might return a lower diagonal
    /// matrix, if we exploit possible symmetry during assembly
    /// (check: m_matrix.symmetry() == true )
    Eigen::SparseSelfAdjointView< typename gsSparseMatrix<T>::Base, Lower> fullMatrix()
    {
        return m_system.matrix().template selfadjointView<Lower>();
    }

protected:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPoissonAssembler.hpp)
#endif

