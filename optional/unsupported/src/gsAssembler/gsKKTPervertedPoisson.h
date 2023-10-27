/** @file gsKKTPervertedPoisson.h

    @brief Provides assembler and solver for the optimal control peverted Poisson.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#pragma once

#include <gsAssembler/gsAssemblerBase2.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsAssembler/gsVisitorKKTPervertedPoisson.h>
#include <gsAssembler/gsVisitorKKTPervertedPoissonBoundary.h>
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNeumannBiharmonic.h> // Neumann boundary integrals


namespace gismo
{

/** @brief
    Implementation of a KKT system for the perverted Poisson, see article:

    Robust preconditioners for PDE-constrained optimization with limited observations.
    K-A Mardal, B F Nielsen and M Nordaas

    UNDER DEVELOPMENT!!!
*/
template <class T>
class gsKKTPervertedPoisson : public gsAssemblerBase2<T>
{
public:
    typedef gsAssemblerBase2<T> Base;

public:

/** @brief Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a vector multi-basis that contains patch-wise
    bases for each unknown u, f and w.
    \param[in] bconditions contains the boundary conditions for the state u
    \param[in] g_des is the desired state.
*/
    gsKKTPervertedPoisson( gsMultiPatch<T>               const & patches,
                       std::vector<gsMultiBasis<T> > const & bases,
                       gsBoundaryConditions<T> const & bconditions,
                       gsBoundaryConditions<T> const & bconditions2,
                       std::vector<gsFunction<T> *> const & desired,
                       const T alpha,
                       dirichlet::strategy dstrategy )
    :  m_desired(desired),
       m_alpha(alpha),
       m_bConditions2(bconditions2)
    {
        m_options.dirStrategy = dstrategy;
        this->initialize(patches,bases,bconditions);
    }

    /// Sets the KKT assembler options
    void setOptions(const gsAssemblerOptions  & options)
    {
        // to do
    }

    /// Set the control parameter
    void setAlpha(T val)
    { 
        m_alpha = val;
    }

    //GISMO_ASSERT( m_patches.nPatches() == 1, "Only implemented for single patch");

    /// Main assembly routine
    void assemble();
    
    void assembleNitsche()
    {
        //Not developed/implemented!
    }
    
    /// Computes the Dirichlet DoF values by interpolation or projection
    void computeDirichletDofs()
    {
        // If we have a homogeneous Dirichlet problem fill boundary
        // DoFs with zeros ( ++ to update)
        if ( m_options.dirStrategy == dirichlet::none)
            m_ddof.setZero( m_dofMappers[1].boundarySize(), 1 );


        // If the Dirichlet strategy is elimination then precompute
        // Dirichlet dofs (m_dofMapper excludes these from the system)
        if ( m_options.dirStrategy == dirichlet::elimination)
        {
            computeDirichletDofsIntpl();//L2Proj
        }

    }

    /// Reconstruct solution field from computed solution vector
    gsField<T> * constructSolution(const gsMatrix<T> & solVector, int unk) const;


protected:

    //Wrapper for all PDE specific initializations
    void initializePdeSpecific()
    {
        // Initialize
        init();
    }

    void init()
    {
        //One for control, state and Lagrange multiplier
        m_dofMappers.resize(3);

        //Control: f
        m_bases.front().getMapper(true, m_dofMappers.front() );

        //State: u, Conforming interfaces mapper with BCs
        if ( m_options.dirStrategy == dirichlet::elimination)
            m_bases[1].getMapper(true, m_bConditions, m_dofMappers[1] );
        else
            m_bases[1].getMapper(true, m_dofMappers[1] );

        //Lagrange multiplier: w
        m_bases.back().getMapper(true, m_dofMappers.back() );

        // Set size of the system
        m_dofs = m_dofMappers.front().freeSize()  + m_dofMappers[1].freeSize() + m_dofMappers.back().freeSize();

    }


    void computeDirichletDofsIntpl()
    {
        const gsDofMapper & mapper = m_dofMappers[1];

        m_ddof.resize( mapper.boundarySize(), 1 ); //--mrhs

        // Iterate over all patch-sides with Dirichlet-boundary conditions
        for ( typename gsBoundaryConditions<T>::const_iterator
                  it = m_bConditions.dirichletBegin();
              it != m_bConditions.dirichletEnd(); ++it )
        {
            const int k   = it->patch();
            const gsBasis<T> & basis = (m_bases[1])[k];

            // Get dofs on this boundary
            gsMatrix<index_t> boundary = basis.boundary(it->side()) ;

            // If the condition is homogeneous then fill with zeros
            if ( it->isHomogeneous() )
            {
                for (index_t i=0; i!= boundary.size(); ++i)
                {
                    const int ii= mapper.bindex( (boundary)(i) , k );
                    m_ddof.row(ii).setZero();
                }
                continue;
            }

            // Get the side information
            int dir = it->side().direction( );
            index_t param = (it->side().parameter() ? 1 : 0);

            // Compute grid of points on the face ("face anchors")
            std::vector< gsVector<T> > rr;
            rr.reserve( this->patches().parDim() );

            for ( short_t i=0; i < this->patches().parDim(); ++i)
            {
                if ( i==dir )
                {
                    gsVector<T> b(1);
                    b[0] = ( basis.component(i).support() ) (0, param);
                    rr.push_back(b);
                }
                else
                {
                    rr.push_back( basis.component(i).anchors().transpose() );
                }
            }

            // Compute dirichlet values
            gsMatrix<T> fpts =
                it->function()->eval( m_patches[it->patch()].eval(  gsPointGrid<T>( rr ) ) );

            // Interpolate dirichlet boundary
            typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());
            typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
            const gsMatrix<T> & dVals =  geo->coefs();

            // Save corresponding boundary dofs
            for (index_t k2=0; k2!= boundary.size(); ++k2)
            {
                const int ii= mapper.bindex( (boundary)(k2) , it->patch() );
                m_ddof.row(ii) = dVals.row(k2);
            }
        }
    }

protected:

    // Desired state g
    std::vector<gsFunction<T> * > m_desired;

    // Parameter
    T m_alpha;

protected:

    /// Boundary conditions for state u
    gsBoundaryConditions<T> m_bConditions2;

    //moved to m_options!
    //dirichlet::strategy m_dirStrategy;
    // iFace::strategy m_intStrategy;

    // Members from gsAssemblerBase
    using Base::m_patches;
    using Base::m_bases;
    using Base::m_bConditions;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_matrix;
    using Base::m_rhs;
    using Base::m_dofs;

};


//////////////////////////////////////////////////
//////////////////////////////////////////////////

template<class T>
void gsKKTPervertedPoisson<T>::assemble()
{

    computeDirichletDofs();

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    // JS2: This need to be imporved for the KKT system
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_bases[0].dim(); ++i)
        nonZerosPerCol *= 4 * m_bases[0].maxDegree(i) + 1;
    
    m_matrix.resize(m_dofs, m_dofs); // Clean global matrix
    m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );
    
    // Resize the right-hand side vector
    m_rhs.setZero(m_dofs, 1 );
    

    // Initialize the element visitor
    gsVisitorKKTPervertedPoisson<T> KKTPP(m_alpha);
    for (size_t np=0; np < m_patches.nPatches(); ++np )
    {
        //Assemble on patch np
        this->apply(KKTPP, np);
    }

    /// Assemble boundary mass matrix (and rhs)
    gsVisitorKKTPervertedPoissonBoundary<T> MassB1(m_desired,boundary::west);
    this->apply(MassB1, 0, boundary::west);

    gsVisitorKKTPervertedPoissonBoundary<T> MassB2(m_desired,boundary::east);
    this->apply(MassB2, 0, boundary::east);

    gsVisitorKKTPervertedPoissonBoundary<T> MassB3(m_desired,boundary::south);
    this->apply(MassB3, 0, boundary::south);

    gsVisitorKKTPervertedPoissonBoundary<T> MassB4(m_desired,boundary::north);
    this->apply(MassB4, 0, boundary::north);


    
    // Enforce Neumann boundary conditions for state. if any
    /*for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions.neumannBegin();
          it != m_bConditions.neumannEnd(); ++it )
    {
        gsVisitorNeumann<T> neumann(*it->function(), it->side());
        this->apply(neumann, it->patch(), it->side() );
    }
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions2.neumannBegin();
          it != m_bConditions2.neumannEnd(); ++it )
    {
        gsVisitorNeumannBiharmonic<T> neumann(*it->function(), it->side());
        this->apply(neumann, it->patch(), it->side() );
    }*/

    // If requested, enforce Dirichlet boundary conditions
    // by Nitsche's method
    // JS2: NOT IMPLEMENTED
    if ( m_options.dirStrategy == dirichlet::nitsche )
    {
        GISMO_ERROR("The Nitsche method is not implemented for this problem");
        //assembleNitsche();
    }

    
    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}

template<class T>
gsField<T> *  gsKKTPervertedPoisson<T>::constructSolution(const gsMatrix<T>& solVector, int unk) const
//gsField<T> & result ) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    index_t shift = 0;
    for (index_t k = 0; k< unk; ++k)
        shift += m_dofMappers[k].freeSize();


    const gsDofMapper & mapper = m_dofMappers[unk];

    gsPiecewiseFunction<T> * sols = new gsPiecewiseFunction<T>;

    gsMatrix<T> coeffs;

    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {
        // Reconstruct solution coefficients on patch p
        const int sz  = m_bases[unk][p].size();
        coeffs.resize( sz, 1);

        for (index_t i = 0; i < sz; ++i)
        {
            if ( mapper.is_free(i, p) ) // DoF value is in the solVector
            {
                coeffs.row(i) = solVector.row( mapper.index(i, p) + shift);
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof.row( mapper.bindex(i, p) );
            }
        }

        sols->addPiecePointer( m_bases[unk][p].makeGeometry( give(coeffs) ) );
    }

    return new gsField<T>(m_patches, typename gsFunctionSet<T>::Ptr(sols), true);
}

} // namespace gismo



