/** @file gsGeneralizedPoissonAssembler.h

    @brief Provides base interface for -Laplace u + alpha u = f

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsAssembler/gsPoissonAssembler.h>
#include <gsAssembler/gsVisitorDg.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{    

// Local overwrite of gsPoissonPde
template<class T=real_t>
class gsGeneralizedPoissonPde : public gsPoissonPde<T>
{
public:    
    gsGeneralizedPoissonPde( ) { }

    /// Constructor
    gsGeneralizedPoissonPde(const gsMultiPatch<T>         &domain,
                            const gsBoundaryConditions<T> &bc,
                            const gsPiecewiseFunction<T>  &rhs,
                            T                             alpha,
                            const gsFunction<T>           *sol = NULL)
    : gsPoissonPde<T>(domain,bc,rhs,sol), m_alpha(alpha) {}
    
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
            os<<"Poisson's equation  -\u0394u + "<<m_alpha<<" u = f ,  with:\n";
            os<<"Source function f= "<< m_rhs <<".\n";
            return os; 
    }
    
    inline T alpha() const { return m_alpha; } 

protected:
    T m_alpha;
    using gsPoissonPde<T>::m_rhs;
}; // class gsGeneralizedPoissonPde


// Local overwrite gsVisitorPoisson
template <class T=real_t, bool paramCoef = false>
class gsVisitorGeneralizedPoisson : public gsVisitorPoisson<T,paramCoef>
{
public:
    gsVisitorGeneralizedPoisson(const gsPde<T> & pde) : gsVisitorPoisson<T,paramCoef>(pde)
    {
        mod_pde_ptr = static_cast<const gsGeneralizedPoissonPde<T>*>(&pde);
    }
    
    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>     & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);
           
            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, physGrad);
           
            localRhs.noalias() += weight * ( bVals.col(k) * rhsVals.col(k).transpose() ) ;
            localMat.noalias() += weight * ( physGrad.transpose() * physGrad + mod_pde_ptr->alpha() * bVals.col(k) * bVals.col(k).transpose() );
        }
    }
    
private:
    const gsGeneralizedPoissonPde<T>* mod_pde_ptr;
    using gsVisitorPoisson<T>::pde_ptr;
    using gsVisitorPoisson<T>::rhs_ptr;
    using gsVisitorPoisson<T>::basisData;
    using gsVisitorPoisson<T>::localRhs;
    using gsVisitorPoisson<T>::localMat;
    using gsVisitorPoisson<T>::physGrad;
    using gsVisitorPoisson<T>::rhsVals;
    using gsVisitorPoisson<T>::md;
};


/** @brief
    This does the same as gsPoissonAssembler, but it does not only assemble
    -Laplace u = f, but -Laplace u + alpha u = f

    \ingroup Assembler
*/
template <class T=real_t>
class gsGeneralizedPoissonAssembler : public gsPoissonAssembler<T>
{
public:
    typedef gsAssembler<T> Base;
    
    gsGeneralizedPoissonAssembler()
    { }

    /** @brief Main Constructor of the assembler object.

    \param[in] pde A boundary value Poisson problem
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] opt A set of options for the assembly process
    */
    gsGeneralizedPoissonAssembler( const gsGeneralizedPoissonPde<T>  & pde,
                                   const gsMultiBasis<T>             & bases)
    {
        Base::initialize(pde, bases, m_options);
    }
    
    /** @brief Main Constructor of the assembler object.

    \param[in] pde A boundary value Poisson problem
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
    */
    gsGeneralizedPoissonAssembler( const gsGeneralizedPoissonPde<T>   & pde,
                                   const gsMultiBasis<T>              & bases,
                                   dirichlet::strategy                dirStrategy,
                                   iFace::strategy                    intStrategy = iFace::glue )
    {
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        Base::initialize(pde, bases, m_options);
    }
    
    /** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] basis a multi-basis that contains patch-wise bases
    \param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.
    \param[in] alpha parameter alpha
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
    
    \ingroup Assembler
    */
    gsGeneralizedPoissonAssembler( gsMultiPatch<T> const         & patches,
                                   gsMultiBasis<T> const         & basis,
                                   gsBoundaryConditions<T> const & bconditions,
                                   const gsFunction<T>           & rhs,
                                   T                             alpha,
                                   dirichlet::strategy           dirStrategy,
                                   iFace::strategy               intStrategy)
    {
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        typename gsPde<T>::Ptr pde( new gsGeneralizedPoissonPde<T>(patches,bconditions,rhs,alpha) );
        Base::initialize(pde, basis, m_options);
    }
    
    /** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] basis a multi-basis that contains patch-wise bases
    \param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.
    \param[in] alpha parameter alpha
    \param[in] opt reference to a gsOptionList object
    
    \ingroup Assembler
    */
    gsGeneralizedPoissonAssembler( gsMultiPatch<T> const         & patches,
                                   gsMultiBasis<T> const         & basis,
                                   gsBoundaryConditions<T> const & bconditions,
                                   const gsFunction<T>           & rhs,
                                   T                             alpha = 0,
                                   const gsOptionList&           opt = gsAssembler<>::defaultOptions() )
    {
        m_options.update(opt);
        typename gsPde<T>::Ptr pde( new gsGeneralizedPoissonPde<T>(patches,bconditions,rhs,alpha) );
        Base::initialize(pde, basis, m_options);
    }
    
    void assemble()
    {
        GISMO_ASSERT(m_system.initialized(), 
                     "Sparse system is not initialized, call initialize() or refresh()");

        // Reserve sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.reserve(nz, this->pde().numRhs());

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        Base::computeDirichletDofs();
        // Clean the sparse system
        //m_system.setZero(); //cf. r7323
        // Assemble volume integrals
        Base::template push<gsVisitorGeneralizedPoisson<T> >();
        // Enforce Neumann boundary conditions
        Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );
        // If requested, enforce Dirichlet boundary conditions by Nitsche's method
        if ( m_options.getInt("DirichletStrategy") == dirichlet::nitsche )
           Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());
        // If requested, enforce Dirichlet boundary conditions by diagonal penalization
        else if ( m_options.getInt("DirichletStrategy") == dirichlet::penalize )
           Base::penalizeDirichletDofs();
        if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
           Base::template pushInterface<gsVisitorDg<T> >();

        // Assembly is done, compress the matrix
        Base::finalize();
    }
protected:
    using gsPoissonAssembler<T>::m_pde_ptr;
    using gsPoissonAssembler<T>::m_options;
    using gsPoissonAssembler<T>::m_system;
    using gsPoissonAssembler<T>::m_bases;
};


} // namespace gismo
