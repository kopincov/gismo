/** @file gsMasorny.hpp

    @brief Provides implementation for the Masonry problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):  Y. Xia, A. Mantzaflaris
*/

#include <gsSelfSuppSurf/gsVisitorNonLinLoad.h>
namespace gismo
{

template<class T> 
void gsMasonry<T>::assemble(const gsMultiPatch<T> & curSolution)
{
    // Reserve sparse system
    m_system.reserve(m_bases[0], m_options, this->pde().numRhs());
    //m_system.matrix().setZero();
    
    // Visit elements
    gsVisitorNonLinLoad<T> visitor(curSolution, *coeff_ptr, *extload_ptr);// Newton iteration Jacobian
    for (size_t np=0; np < m_pde_ptr->domain().nPatches(); ++np )
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    // Enforce Neumann boundary conditions
    Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );

    // Apply any point-loads
    if ( m_pLoads.numLoads() != 0 )
        applyLoads();
    
    // Assembly is done, compress the matrix
    m_system.matrix().makeCompressed();
}

template<class T> 
void gsMasonry<T>::assemble()
{
    Base::assemble();

    if ( m_pLoads.numLoads() != 0 )
    {
        // Apply any point-loads
        applyLoads();
        
        // Assembly is done, compress the matrix
        m_system.matrix().makeCompressed();
    }
}

template<class T> 
void gsMasonry<T>::assembleSystem(const gsMultiPatch<T> & curSolution)
{
    this->computeDirichletDofs();

    // Reserve sparse system
    m_system.reserve(m_bases[0], m_options, this->pde().numRhs());
    
    // Visit elements
    gsVisitorMasonryRhs<T> visitor(curSolution, *coeff_ptr, *extload_ptr);
    for (size_t np=0; np < m_pde_ptr->domain().nPatches(); ++np )
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    if ( m_pLoads.numLoads() != 0 )
    {
        // Apply any point-loads
        applyLoads();
        
        // Assembly is done, compress the matrix
        m_system.matrix().makeCompressed();
    }
}

template<class T> 
void gsMasonry<T>::applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<index_t> acts;

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        //m_patches
        if ( m_pLoads[i].parametric )
        {
            m_bases.front().basis(m_pLoads[i].patch).active_into( m_pLoads[i].point, acts );
            m_bases.front().basis(m_pLoads[i].patch).eval_into  ( m_pLoads[i].point, bVals);
        }
        else
        {
            gsWarn<< "Point loads parametric for now.\n";
            //m_pLoads[i].toParametric(); // to do
        }

        // translate patch-local indices to global dof indices
        const gsDofMapper & mapper = m_system.colMapper(0);
        mapper.localToGlobal(acts, m_pLoads[i].patch, acts);

        for (index_t k=0; k < acts.rows(); ++k)
        {            
            if ( mapper.is_free_index( acts(k,0) ) ) // interior node?
                m_system.rhs()(acts(k,0), 0) += bVals(k,0) * m_pLoads[i].value[2];
            else
            {
                gsWarn<< "Boundary Load ?\n";
            }
        }
    }

}


template<class T>
void  gsMasonry<T>::updateSolution(const gsMatrix<T>& solVector, 
    gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(this->numDofs() == m_system.rhs().rows(),
                 "Something went wrong, assemble() not called?");

    for (index_t p=0; p < m_pde_ptr->domain().nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();

        for (index_t j = 0; j < 1; ++j) //m_dim == 1 scalar problem !
        {
            const gsDofMapper & mapper = m_system.colMapper(j);
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

//template<class T>
//void gsMasonry<T>::CalGeoHesnMatr(gsMatrix<T> & GeoHesnMatr)
//{

//}




}// namespace gismo
