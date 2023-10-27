/** @file gsShellAssembler.hpp

    @brief Provides non linear elasticity system matrices for thin shells.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Goyal, A. Mantzaflaris
*/

#include <gsThinShell/gsShellAssembler.h>

#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsField.h>
#include <gsUtils/gsPointGrid.h>

// Element visitors
#include <gsThinShell/gsVisitorLinShell.h>
#include <gsThinShell/gsVisitorNonLinShell.h>
#include <gsThinShell/gsVisitorShellNeumann.h>
#include <gsThinShell/gsVisitorMassShell.h>

// ---
#include <gsCore/gsBoxTopology.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>


namespace gismo
{

template<class T>
gsShellAssembler<T>::gsShellAssembler( gsMultiPatch<T> const & patches,
                                       T thickness,
                                       T E_modulus,
                                       T poissons_ratio,
                                       gsBoundaryConditions<T> const & bconditions,
                                       const gsFunction<T> & surface_force,
                                       const clamped_t & clamped, 
                                       const gsPointLoads<T> & pLoads
                                       )
    : m_thickness(thickness),
    m_surfaceForce(&surface_force),
    m_clamped(clamped),
    m_pLoads(pLoads)

{
    m_options.dirStrategy = dirichlet::elimination;
    // Initialize material properties
    m_lambda = E_modulus * poissons_ratio / ( (1 + poissons_ratio)*(1-2*poissons_ratio)) ;
    m_mu     = E_modulus / (2*(1 + poissons_ratio)) ;
    
    this->initialize(patches,gsMultiBasis<T>(patches),bconditions);

}

template<class T>
gsShellAssembler<T>::gsShellAssembler( gsMultiPatch<T> const & patches,
                                       T thickness,
                                       T rho,
                                       T E_modulus,
                                       T poissons_ratio,
                                       gsBoundaryConditions<T> const & bconditions,
                                       const gsFunction<T> & surface_force,
                                       const clamped_t & clamped, 
                                       const gsPointLoads<T> & pLoads
                                       )
    : m_thickness(thickness),
    m_rho(rho),
    m_surfaceForce(&surface_force),
    m_clamped(clamped),
    m_pLoads(pLoads)

{
    m_options.dirStrategy = dirichlet::elimination;
    // Initialize material properties
    m_lambda = E_modulus * poissons_ratio / ( (1 + poissons_ratio)*(1-2*poissons_ratio)) ;
    m_mu     = E_modulus / (2*(1 + poissons_ratio)) ;
    
    this->initialize(patches,gsMultiBasis<T>(patches),bconditions);
    this->assembleMassMatrix();

}

template<class T>
void gsShellAssembler<T>::initializePdeSpecific()
{
    // z-clamped Dirichlet dofs for the side requested
    initMappers(m_clamped);

    // Determine system size
    m_dofs = m_dofMappers[0].freeSize()
           + m_dofMappers[1].freeSize()
           + m_dofMappers[2].freeSize();

    // Add shifts to global matrix indices to get the global dof
    // number for each unknown coordinate
    m_dofMappers[1].setShift( m_dofMappers[0].freeSize()                              );
    m_dofMappers[2].setShift( m_dofMappers[0].freeSize() + m_dofMappers[1].freeSize() );
    m_dofMappers[1].setBoundaryShift( m_dofMappers[0].boundarySize() );
    m_dofMappers[2].setBoundaryShift( m_dofMappers[0].boundarySize() +
                                      m_dofMappers[1].boundarySize());
}



template<class T>
void gsShellAssembler<T>::assembleNeumann()
{
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions.neumannBegin();
          it != m_bConditions.neumannEnd(); ++it )
    {
        gsVisitorShellNeumann<T> neumann(*it->function(), it->side());

        // Note: it->unknown()
        this->apply(neumann, it->patch(), it->side() );
    }
}

template<class T> 
void gsShellAssembler<T>::applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<index_t> acts,globalActs;

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        if ( m_pLoads[i].parametric )
        {
            m_bases.front().basis(m_pLoads[i].patch).active_into( m_pLoads[i].point, acts );
            m_bases.front().basis(m_pLoads[i].patch).eval_into  ( m_pLoads[i].point, bVals);
        }
        else
        {
            gsWarn<< "Point loads parametric for now.\n";
        }

        // translate patch-local indices to global dof indices
        for (size_t j = 0; j< 3; ++j)
        {
            if (m_pLoads[i].value[j] != 0.0)
            {
                m_dofMappers[j].localToGlobal(acts, m_pLoads[i].patch, globalActs);

                for (index_t k=0; k < globalActs.rows(); ++k)
                {
                    if (int(globalActs(k,0)) < m_dofs)
                        m_rhs(globalActs(k,0), 0) += bVals(k,0) * m_pLoads[i].value[j];
                }
            }
        }
    }
}

template<class T>
void gsShellAssembler<T>::assembleMassMatrix()
{
    //computeDirichletDofs();
    index_t numDirichlet = 0;
    for (int i = 0; i < 3; ++i)
        numDirichlet += m_dofMappers[i].boundarySize();
    m_ddof.setZero(numDirichlet, 1);

    if (m_dofs == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    // Pre-allocate non-zero elements for each column of the sparse matrix
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
        nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;

    m_matrix = gsSparseMatrix<T>(m_dofs, m_dofs); // Clean matrices
    m_matrix.reserve( gsVector<index_t>::Constant(m_dofs, nonZerosPerCol) );
        
    // Assemble mass matrix integrals
    gsVisitorMassShell<T> visitor(m_rho,m_thickness);
    for (size_t np = 0; np < m_patches.nPatches(); ++np)
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();
    
    // Store the mass matrix
    m_matrix.swap(m_massMatrix);
}

template<class T>
void gsShellAssembler<T>::initMappers(const std::vector<std::pair<patchSide,index_t> > & clamped)
{
    typedef std::vector<std::pair<patchSide,index_t> >::const_iterator cIterator;

    m_dofMappers.resize(3);
    // Note the last false: mappers are not finalized
    m_bases.front().getMapper(true, m_bConditions, 0, m_dofMappers[0], false ); // x -coord
    m_bases.front().getMapper(true, m_bConditions, 1, m_dofMappers[1], false ); // y -coord
    m_bases.front().getMapper(true, m_bConditions, 2, m_dofMappers[2], false); // z -coord
    

    for( cIterator it = clamped.begin(); it!=clamped.end(); ++it)
    {
        gsDofMapper & mapper  = m_dofMappers[it->second];
        const patchSide & cur = it->first;
        // Get boundary dofs
        gsMatrix<index_t> bDofs = m_bases[0][cur.patch].boundary(cur);
        
        // Cast to tensor b-spline basis
        const gsTensorBSplineBasis<2,T> * tp = 
            dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&m_bases[0][cur.patch]);
        
        if ( tp != NULL) // clamp adjacent dofs
        {
            const int str = tp->stride( cur.direction() );            
            if ( cur.parameter() )
            {
                for ( index_t k=0; k<bDofs.size(); ++k)
                    mapper.matchDof( cur.patch, (bDofs)(k,0),
                                     cur.patch, (bDofs)(k,0) - str );
            }
            else
            {
                for ( index_t k=0; k<bDofs.size(); ++k)
                    mapper.matchDof( cur.patch, (bDofs)(k,0),
                                     cur.patch, (bDofs)(k,0) + str );
            }
        }
        else
            gsWarn<<"Unable to apply clamped condition.\n";

    }

    // Need to finalize the mappers
    m_dofMappers[0].finalize();
    m_dofMappers[1].finalize();
    m_dofMappers[2].finalize();
}


template<class T>
void gsShellAssembler<T>::assemble()
{
    //computeDirichletDofs();
    index_t numDirichlet = 0;
    for (int i = 0; i < 3; ++i)
        numDirichlet += m_dofMappers[i].boundarySize();
    m_ddof.setZero(numDirichlet, 1);

    if (m_dofs == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    int nonZerosPerCol = 3;
    for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
        nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;

    m_matrix = gsSparseMatrix<T>(m_dofs, m_dofs); // Clean matrices
    m_matrix.reserve( gsVector<index_t>::Constant(m_dofs, nonZerosPerCol) );
        
    // Resize the load vector
    m_rhs.setZero(m_dofs, 1 );

    // Assemble volume stiffness and load vector integrals
    gsVisitorLinShell<T> visitor(m_thickness,m_lambda,m_mu,*m_surfaceForce);
    for (size_t np = 0; np < m_patches.nPatches(); ++np)
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    // Enforce Neumann boundary conditions
    assembleNeumann();
    
    // Apply any point-loads
    if ( m_pLoads.numLoads() != 0 )
        applyLoads();

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}

template<class T>
void gsShellAssembler<T>::assemble(const gsMultiPatch<T> & deformed)
{
    if ( m_ddof.size() == 0 )
    {
        assemble();
        return;
    }

    // Initialize the matrix and rhs vector
    m_matrix.setZero();
    
    // Resize the load vector
    m_rhs.setZero(m_dofs, 1 );

    gsVisitorNonLinShell<T> visitor(m_thickness,m_lambda,m_mu,
                                    *m_surfaceForce, deformed.patch(0));

    // Assemble volume stiffness and load vector integrals
    for (size_t np = 0; np < m_patches.nPatches(); ++np)
    {
        //visitor.setDeformed( deformed.patch(np) );

        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    // Enforce Neumann forces
    assembleNeumann();
    
    // Apply any point-loads
    if ( m_pLoads.numLoads() != 0 )
        applyLoads();

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}

template<class T> // AM: is non-homogeneous not called yet
void gsShellAssembler<T>::computeDirichletDofsIntpl()
{
    const gsDofMapper & mapper = m_dofMappers.front();
    
    m_ddof.resize( mapper.boundarySize(), 1 ); //--mrhs
    
    // Iterate over all patch-sides with Dirichlet-boundary conditions
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
            for (index_t i=0; i!= boundary.size(); ++i)
            {
                const index_t ii= mapper.bindex( (boundary)(i) , k );
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


template<class T>
void  gsShellAssembler<T>::updateSolution(const gsMatrix<T>& solVector, 
                                          gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    for (size_t p = 0; p < m_patches.nPatches(); ++p)
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();

        for (index_t j = 0; j < 3; ++j)
        {
            const gsDofMapper & mapper = m_dofMappers[j];
            for (index_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector
                {
                    coeffs(i,j) += solVector( mapper.index(i, p), 0);
                }
            }
        }
    }
}

template<class T>
void  gsShellAssembler<T>::constructSolution(const gsMatrix<T>& solVector, 
                                             gsMultiPatch<T>& result) const
{
    // The IETI method cannot pass this assertion, since assemble() is never called on the whole domain
    //GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");
    
    // The final solution is the deformed shell, therefore we add the
    // solVector to the undeformed coefficients
    result = m_patches;

    for (size_t p = 0; p < m_patches.nPatches(); ++p)
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size(); // m_patches[p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();

        for (index_t j = 0; j < 3; ++j) // For all components x, y, z
        {
            const gsDofMapper & mapper = m_dofMappers[j];// grab mapper for this component

            for (index_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector ?
                {
                    coeffs(i,j) += solVector( mapper.index(i, p), 0);
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    coeffs(i,j) += m_ddof(mapper.bindex(i, p), 0);
                }
            }
        }
    }
}


}// namespace gismo
