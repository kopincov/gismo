/** @file gsStokesAssembler.h

    @brief Provides assembler and solver for the Poisson equation, incl. adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, J. Sogn
*/

#pragma once

#include <gsAssembler/gsAssemblerBase2.h>
#include <gsAssembler/gsAssemblerOptions.h>

#include <gsAssembler/gsVisitorStokes.h> // Stiffness volume integrals
//#include <gsAssembler/gsVisitorStokesDC.h> // div-conforming 

//#include <gsAssembler/gsVisitorNeumannMixed.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNeumann2.h> // Neumann boundary integrals

//#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals

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
class gsStokesAssembler2 : public gsAssemblerBase2<T>
{
public:
    typedef gsAssemblerBase2<T> Base;

public:

/** @brief Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a vector multi-basis that contains patch-wise bases for each unknown
    \param[in] bconditions contains the boundary conditions for the velocity
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.
*/
    gsStokesAssembler2( gsMultiPatch<T>               const & patches,
                       std::vector<gsMultiBasis<T> > const & bases,
                       gsBoundaryConditions<T> const & bconditions,
                       const gsFunction<T>           & rhs, 
                       dirichlet::strategy dstrategy )
    :  Base(patches), 
       m_rhsFun(&rhs),
       m_viscosity(1.0)
    {
        //TODO: override     initialize().
        // and use them in the constructor

        m_bConditions = bconditions;
        m_elementType = 0;
        m_options.dirStrategy=dstrategy;
        //to input: pressure, compute velocity space according to element type
        m_bases = bases;

        // Initialize
        init();

    }

    /// Sets the Stokes assembler options
    void setOptions(const gsAssemblerOptions  & options)
    {
        // to do
    }

    /// Set the kinematic viscosity (default is \f$ \nu = 1.0\f$)
    void setViscosity(T val)
    { 
        m_viscosity = val;
    }

    /// Main assembly routine
    void assemble();
    
    void assembleNeumann()
    {
        for ( typename gsBoundaryConditions<T>::const_iterator
                  it = m_bConditions.neumannBegin();
              it != m_bConditions.neumannEnd(); ++it )
        {
            gsVisitorNeumann2<T> neumann(*it->function(), it->side());
            
            // Note: it->unknown()
            this->apply(neumann, it->patch(), it->side() );
        }
    }

    void assembleNitsche()
    {
        /*
        for ( typename gsBoundaryConditions<T>::const_iterator
                  it = m_bConditions.dirichletBegin();
              it != m_bConditions.dirichletEnd(); ++it )
        {
            gsVisitorNitsche<T> nitsche(*it->function(), penalty(it->patch()), it->side());
            
            // Note: it->unknown()
            this->apply(nitsche, it->patch(), it->side() );
        }
        */
    }
    
    /// Computes the Dirichlet DoF values by interpolation
    void computeDirichletDofs();

    /// Reconstruct solution field from computed solution vector
    gsField<T> constructSolution(const gsMatrix<T> & solVector, int unk) const;

    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() returns a lower diagonal matrix,
    /// since we exploit symmetry during assembly.
    Eigen::SparseSelfAdjointView< typename gsSparseMatrix<T>::Base, Lower> fullMatrix()
    {
        return m_matrix.template selfadjointView<Lower>();
    }

protected:

    // TODO: initializes the pde specfic stuff
    //void initializePdeSpecific(){ }

    void init()
    {
        // 1. Common v-basis ?--> copy mappers, no=RT
        // 2. Div conforming ? --> no=block diagonal


        //space e_1 B_i, e_2 B_j, e_3 B_k

        const index_t tarDim = m_rhsFun->targetDim(); //Target dimension for velocity

        switch (m_elementType)
        {
        case 0: // eg. Taylor-Hood

            //One for velocity and one for pressure
            m_dofMappers.resize(2);

            // Conforming interfaces mapper for velocity components with BCs
            if ( m_options.dirStrategy == dirichlet::elimination)
                m_bases.front().getMapper(true, m_bConditions, m_dofMappers.front() );
            else
                m_bases.front().getMapper(true, m_dofMappers.front() );


            // Conforming interfaces mapper for pressure
            m_bases.back().getMapper(true, m_dofMappers.back() );

            // Set size of the system -----------
            m_dofs = tarDim * m_dofMappers.front().freeSize()  + m_dofMappers.back().freeSize();
            
            break;

        case 1: // eg. Raviart-Thomas

            // distinct basis for each velocity component
            m_dofMappers.resize(tarDim+1);

            // Conforming interfaces mapper for each velocity component with BCs
            if ( m_options.dirStrategy == dirichlet::elimination)
                for ( index_t i = 0; i!= tarDim; ++i) 
                    //Note: "-1": common BC's for all components of velocity
                    m_bases[i].getMapper(true, m_bConditions, -1, m_dofMappers[i] );
            else
                for ( index_t i = 0; i!= tarDim; ++i) 
                    m_bases[i].getMapper(true, m_dofMappers[i] );

            // Conforming interfaces mapper for pressure
            m_bases.back().getMapper(true, m_dofMappers.back() );
            
            // Set size of the system -----------
            m_dofs = m_dofMappers.back().freeSize();
            for ( index_t i = 0; i!= tarDim; ++i)
                m_dofs += m_dofMappers[i].freeSize();

            break;

        default :
            GISMO_ERROR("Error while initializing for element "<< m_elementType);
            break;
        }

    }


    void assembleTH()
    {
        // Initialize the element visitor
        gsVisitorStokes<T> stokes(*m_rhsFun, m_viscosity);
        for (size_t np=0; np < m_patches.nPatches(); ++np )
        {
            //Assemble on patch np
            this->apply(stokes, np);
        }
    }

/*
    void assembleDC()
    {
        // Initialize the element visitor
        gsVisitorStokesDC<T> stokes(*m_rhsFun, m_viscosity);
        for (size_t np=0; np < m_patches.nPatches(); ++np )
        {
            //Assemble on patch np
            this->apply(stokes, np);
        }
    }
//*/  
  
protected:

    // Right hand side function
    const gsFunction<T> * m_rhsFun;

    // Kinematic viscosity
    T m_viscosity;

    // Element type
    int m_elementType;

protected:

    // Geometrical transformation type
    //ValueTransformationType m_geoTrans;

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
void gsStokesAssembler2<T>::assemble()
{
    // If we have a homogeneous Dirichlet problem fill boundary
    // DoFs with zeros ( ++ to update)
    if ( m_options.dirStrategy == dirichlet::none)
        m_ddof.setZero( m_dofMappers[0].boundarySize(), m_rhsFun->targetDim() );
    
    // If the Dirichlet strategy is elimination then precompute
    // Dirichlet dofs (m_dofMapper excludes these from the system)
    if ( m_options.dirStrategy == dirichlet::elimination)
        computeDirichletDofs();
    
    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_bases.front().dim(); ++i)
        nonZerosPerCol *= 4 * m_bases.front().maxDegree(i) + 1;
    
    m_matrix.resize(m_dofs, m_dofs); // Clean global matrix
    m_matrix.reserve( gsVector<index_t>::Constant(m_dofs, nonZerosPerCol) );
    
    // Resize the right-hand side vector
    m_rhs.setZero(m_dofs, 1 );
    
    // Assemble volume integrals
    switch (m_elementType)
    {
    case 0: // Taylor-Hood
        assembleTH();
        break;
    case 1: // Div-conforming
        //assembleDC();
        break;
    default :
        GISMO_ERROR("Something went terribly wrong.");
        break;
    }

    // If requested, enforce velocity Dirichlet boundary conditions
    // by Nitsche's method
    if ( m_options.dirStrategy == dirichlet::nitsche )
        assembleNitsche();
    
    // Enforce Neumann boundary conditions for velocity. if any
    assembleNeumann();
    
    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}


template<class T>
void gsStokesAssembler2<T>::computeDirichletDofs()
{
    const gsDofMapper & mapper = m_dofMappers.front();
    
    m_ddof.resize( mapper.boundarySize(), m_rhsFun->targetDim() ); //--mrhs
    
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions.dirichletBegin();
          it != m_bConditions.dirichletEnd(); ++it )
    {
        const int unk = it->unknown();
        const int k   = it->patch();
        const gsBasis<T> & basis = (m_bases[unk])[k];

        // Get dofs on this boundary
        gsMatrix<index_t> boundary = basis.boundary(it->side()) ;

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t k2=0; k2!= boundary.size(); ++k2)
            {
                const int ii= mapper.bindex( (boundary)(k2) , it->patch() );
                m_ddof.row(ii).setZero();
            }
            continue;
        }

        // Get the side information
        int dir = it->side().direction();
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( this->patches().parDim() );

        for ( int i=0; i < this->patches().parDim(); ++i)
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
        typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors( fpts );
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t k3=0; k3!= boundary.size(); ++k3)
        {
            const int ii= mapper.bindex( (boundary)(k3) , it->patch() );
            m_ddof.row(ii) = dVals.row(k3);
        }
    }
}


template<class T>
gsField<T> gsStokesAssembler2<T>::constructSolution(const gsMatrix<T>& solVector, int unk) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");
    GISMO_ASSERT(m_dofs ==solVector.rows(), "Something went wrong, is solution vector valid?");

    const gsDofMapper & mapper = m_dofMappers[unk];

    gsPiecewiseFunction<T> * sols = new gsPiecewiseFunction<T>;
    
    const index_t td  = m_rhsFun->targetDim();
    const index_t dim = ( unk == 0 ? td : 1 );
    gsMatrix<T> coeffs;

    // Point to the correct entries of the solution vector
    gsAsConstMatrix<T> solV = ( unk == 0 ? 
        gsAsConstMatrix<T>(solVector.data(), m_dofMappers[0].freeSize(), dim)
        :
        gsAsConstMatrix<T>(solVector.data() + td*m_dofMappers[0].freeSize(), 
                           m_dofMappers[1].freeSize(), 1)
        );

    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {    
        // Reconstruct solution coefficients on patch p
        const int sz  = m_bases[unk][p].size();
        coeffs.resize( sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if ( mapper.is_free(i, p) ) // DoF value is in the solVector
            {
                coeffs.row(i) = solV.row( mapper.index(i, p) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                // Assuming Taylor-Hood 
                coeffs.row(i) = m_ddof.row( mapper.bindex(i, p) );
                //TODO: linearized (stacked) as a big column 
            }
        }
        
        sols->addPiecePointer( m_bases[unk][p].makeGeometry( coeffs ) );
    }

    return gsField<T>(m_patches, typename gsPiecewiseFunction<T>::Ptr(sols), true);
}

} // namespace gismo



