/** @file gsFittingSystem.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsMatrix/gsMatrix.h>
#include <gsFitting/gsFittingUtilsGen.h>
#include <gsFitting/gsFittingSystem.h>

namespace gismo
{

template<class GenBasis, class T>
_gsFittingSystem<GenBasis, T>::
_gsFittingSystem(GenBasis& basis, unsigned dim_im,
                 bool split_dim)
: m_local_global(basis, dim_im), m_A_mat(), m_bc()
{
    init_system(dim_im, split_dim);
}

template<class GenBasis, class T>
void _gsFittingSystem<GenBasis, T>::
init_system(unsigned dim_im, bool split_dim)
{
    m_split_dim = split_dim;
    m_dim = dim_im;

    int total_num_basis = m_local_global.size_totalDOF();
    int num_basis = m_local_global.size_freeDof();
    if(! split_dim)
        num_basis *= m_dim;

    int s_sm = 1;
    if(m_split_dim)
        s_sm = m_dim;

    m_rhs.setZero(num_basis, s_sm);
    m_A_mat.resize(num_basis, num_basis);
    m_A_mat.setZero();

    int s_bc = total_num_basis - num_basis;
    // if s_bc == 0, then there is no boundary condition
    if(s_bc > 0)
        m_bc.setZero(total_num_basis - num_basis, m_dim);
    else
        m_bc.setZero(0,0);
    setDirichletDof();
}


template<class GenBasis, class T> void
_gsFittingSystem<GenBasis, T>::
add_rhs(int ind_dof, int patch, const gsMatrix<T>& loc_rhs)
{
    if(this->is_free(ind_dof, patch))
    {
        const int global_ind = m_local_global.localToGlobal
            (ind_dof, patch, m_split_dim);
        if(m_split_dim)
        {
            if(loc_rhs.cols() == 1)
                m_rhs.row(global_ind) += loc_rhs.transpose();
            else
                m_rhs.row(global_ind) += loc_rhs;
        }
        else {
            if(loc_rhs.cols() == 1)
            {
                for(unsigned d = 0;d < m_dim;d++)
                    m_rhs(global_ind + d,0) += loc_rhs(d, 0);
            }
            else
            {
                for(unsigned d = 0;d < m_dim;d++)
                    m_rhs(global_ind + d,0) += loc_rhs(0, d);
            }
        }
    }
}

/// Note that a little improvement is possible:
/// we do not need to compute the values of
/// the lines that are removed
template<class GenBasis, class T> void
_gsFittingSystem<GenBasis, T>::
add_matrix(int ind_i, int ind_j, int patch,
           const std::vector<T>& localA)
{
    if(this->is_free(ind_i, patch))
    {
        /*   GISMO_ASSERT(this->is_free(ind_i, patch), "The primal (index i) should be free (the lines that are not free are not considered!!"); */
        int global_i = m_local_global.localToGlobal
            (ind_i, patch, m_split_dim);
        int global_j;
        bool is_j_free = this->is_free(ind_j, patch);
        if(is_j_free)
            global_j = m_local_global.localToGlobal
                (ind_j, patch, m_split_dim);
        else
            global_j = m_local_global.bindex(ind_j, patch);

        /// loop on the dimensions
        int d1_2, s = 1;
        /// In the case where the dimensions can be split, we fill in
        /// the matrix only when d1 == d2  (block diagonal with identical blocks)
        if(! m_split_dim)
            s = m_dim;


        if(this->is_free(ind_j, patch))
        {
            for(int d1 = 0;d1 < s;d1++)
            {
                for(int d2 = 0;d2 < s;d2++)
                {
                    d1_2 = d1*m_dim + d2;
                    if((!m_split_dim || d1 == d2)
                       && localA[d1_2] != 0.)
                    {
                        m_entries_mat.add(global_i + d1, global_j + d2,
                                          localA[d1_2]);
                    }
                }
            }
        }
        else
        {
            for(int d1 = 0;d1 < s;d1++)
            {
                for(int d2 = 0;d2 < s;d2++)
                {
                    d1_2 = d1*m_dim + d2;
                    if(localA[d1_2] != 0.)
                    {
                        for(unsigned d3 = 0;d3 < m_dim;d3++)
                        {
                            if(m_split_dim)
                            {
                                m_rhs(global_i, d3)
                                    -= localA[d1_2]
                                    * m_bc(global_j, d3);
                            }
                            else if(d1 == d2)
                            {
                                m_rhs(global_i + d3, 0)
                                    -= localA[d1_2]
                                    * m_bc(global_j, d3);
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class GenBasis, class T> gsFunctionSet<T>*
_gsFittingSystem<GenBasis, T>::solve(bool print_messages)
{
    m_A_mat.setFrom(m_entries_mat);
    m_A_mat.makeCompressed();

    typename gsSparseSolver<T>::BiCGSTABILUT solver( m_A_mat );

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
        m_result.setZero(0, 0);
        return NULL;
    }

    gsMatrix<T> x;
    // Solves for many right hand side  columns
    x = solver.solve(m_rhs);
    if(print_messages)
        gsInfo << "Inverse computed" << std::endl;

    /// In case we could not split the dimensions for the resolution,
    /// we now split the dimensions in the vector x
    if(! m_split_dim)
    {
        int num_basis = m_local_global.size_freeDof();
        m_result.setZero(num_basis, m_dim);
        for(int i = 0;i < num_basis;i++)
        {
            for(unsigned d1 = 0;d1 < m_dim;d1++)
                m_result(i, d1) = x(i*m_dim + d1);
        }
    }
    else  m_result = x;

    gsFittingSystem<GenBasis, T>* _this
        = static_cast<gsFittingSystem<GenBasis, T>* >(this);
    return _this->makeGeometry();
}

template<class T> gsFunctionSet<T>*
gsFittingSystem<gsBasis<T>, T>::makeGeometry()
{
    gsBasis<T>& basis = m_local_global.basis();
    if(m_bc.cols() > 0)
    {
        const int numDOF = m_local_global.size_totalDOF();
        gsMatrix<> coefs(numDOF, m_dim);
        for (int i = 0; i != numDOF; i++)
        {
            if(m_local_global.is_free(i, 0))
            {
                const int globalI
                    = m_local_global.localToGlobal(i, 0, true);
                coefs.row(i) = m_result.row(globalI);
            }
            else
            {
                const int bIndex = m_local_global.bindex(i, 0);
                coefs.row(i) = m_bc.row(bIndex);
            }

        }
        return basis.makeGeometry( give(coefs) ).release();
    }
    else
        return basis.makeGeometry( give(m_result) ).release();
}

template<class T> gsFunctionSet<T>*
gsFittingSystem<gsMultiBasis<T>, T>::makeGeometry()
{
    gsMultiBasis<T>& basis = m_local_global.basis();
    std::vector< patchSide > boundaries
        = basis.topology().boundaries();
    std::vector< boundaryInterface > interfaces
        = basis.topology().interfaces();
    unsigned size = basis.nBases();
    typename gsMultiPatch<T>::PatchContainer newPatches;

    for(unsigned patch = 0;patch < size;patch++)
    {
        const int numBasisFun = basis[patch].size();
        gsMatrix<> local_coefs(numBasisFun, m_dim);
        for (int i = 0; i < numBasisFun; i++)
        {
            if(m_local_global.is_free(i, patch))
            {
                const int globalI = m_local_global.
                    localToGlobal(i, patch, true);
                local_coefs.row(i) = m_result.row(globalI);
            }
            else
            {
                GISMO_ASSERT(m_bc.size() > 0, "ERROR: All the degrees of freedom should be active if there is no boundary condition");
                const int globalI = m_local_global.bindex(i, patch);
                local_coefs.row(i) = m_bc.row(globalI);

            }

        }
        gsGeometry<T>* geom_loc = basis[patch].makeGeometry
            (give(local_coefs)).release();
        GISMO_ASSERT(geom_loc != NULL,
                     "Error during the construction of the geometry associated with the linear system");
        newPatches.push_back(geom_loc);
    }
    return new gsMultiPatch<T>(newPatches, boundaries,
                               interfaces);
}

template<class GenBasis, class T> void
_gsFittingSystem<GenBasis, T>::
setDirichletDof()
{
    gsBoundaryConditions<T>* bc = m_local_global.get_bc();
    if(bc != NULL)
    {
        int s;
        for ( typename gsBoundaryConditions<T>::const_iterator
                  iter = bc->dirichletBegin();
              iter != bc->dirichletEnd(); ++iter )
        {
            int patch = iter->patch();
            gsBasis<T>& basis = getBasisGen<T>
                (m_local_global.basis(), patch);
            const gsMatrix<index_t> ind_boundary
                = basis.boundary(iter->side());
            gsGeometry<T>* boundary
                = static_cast<gsGeometry<T>*>(iter->function().get());
            const gsMatrix<T> & coeffs = boundary->coefs();
            // Save corresponding boundary dofs

            s = ind_boundary.cols();
            for (index_t l = 0; l < s;l++)
            {
                const int ii = m_local_global.bindex
                    (ind_boundary.at(l), iter->patch());
                m_bc.row(ii) = coeffs.row(l);
            }

        }
    }
}


}// namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFittingSystem.hpp)
#endif
