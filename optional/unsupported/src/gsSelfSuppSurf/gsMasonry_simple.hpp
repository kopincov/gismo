/** @file gsMasorny.hpp

    @brief Provides implementation for the Masonry probelm

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, Y. Xia
*/

//#include <gsAssembler/gsMasonry.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsField.h>
#include <gsUtils/gsPointGrid.h>

#include <gsSelfSuppSurf/gsVisitorNonLinLoad.h>

namespace gismo
{


template<class T> // assumes that assemble is called first
void gsMasonry_simple<T>::assemble(const gsMultiPatch<T> & curSolution)
{
    //m_matrix = gsSparseMatrix<T>(m_dofs, m_dofs); // Clean matrices
    //m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );
    m_matrix.setZero();

    // Resize the load vector
    m_rhs.setZero(m_dofs, m_rhsFun->targetDim() );

    // Visit elements
    gsVisitorNonLinLoad<T> visitor(curSolution);
    for (index_t np=0; np < m_patches.nPatches(); ++np )
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        visitor.setPatch(np);
        this->apply(visitor, np);
    }

    // Enforce Neumann boundary conditions
    this->assembleNeumann();

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}


template<class T>
void  gsMasonry_simple<T>::updateSolution(const gsMatrix<T>& solVector, 
                                               gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    for (index_t p=0; p < m_patches.nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();

        for (index_t j = 0; j < 1; ++j) //m_dim == 1 scalar problem !
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



}// namespace gismo
