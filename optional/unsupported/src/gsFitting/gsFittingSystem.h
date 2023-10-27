/** @file gsFittingSystem.h

    @brief Contains the class gsFittingSystem that contains the
    sparse matrix and the RHS

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once


#include <gsCore/gsForwardDeclarations.h>
#include <gsMatrix/gsMatrix.h>

namespace gismo
{

/**
   @brief
   Class containing the sparse matrix and the RHS for gsFitting
   (TODO use gsSparseSystem?)
**/
template<class GenBasis, class T=real_t>
class _gsFittingSystem
{

public:

    /// Constructors
    _gsFittingSystem(GenBasis& basis, unsigned dim,
                     bool split_dim);

    inline bool is_free(int ind, int patch)
    {
        return m_local_global.is_free(ind, patch);
    }

    inline bool split_dim(){  return m_split_dim;  }

    void actualize_mapper(bool repair)
    {   m_local_global.actualize_mapper(repair);   }

    void set_bc(gsBoundaryConditions<T>* bc)
    {  m_local_global.set_bc(bc);  }

    /// Solves the system and construct the
    /// geometry associated with the solution.
    gsFunctionSet<T>* solve(bool print_messages);

    void add_rhs(int ind_dof, int patch, const gsMatrix<T>& loc_rhs);
    void add_matrix(int ind_i, int ind_j, int patch,
                    const std::vector<T>& localA);

    void init_system(unsigned dim_im, bool split_dim);

    void setDirichletDof();


private:
    /// Construct the geometry of the result from the coefficients
    /// obtained by solving the linear system.
    gsFunctionSet<T>* makeGeometry();


protected:

    /// The dimension of the image
    /// (dimension of the space where the target domain is embedded)
    unsigned m_dim;

    /// Calls the mapper when necessary (multipatch)
    /// to obtain a global index from a local index
    gsLocalGlobal<GenBasis, T> m_local_global;

    /// If true, we solve the system independently on each dimension
    bool m_split_dim;

    /// The matrix associated with the linear system.
    /// Used in compute linear (and possibly for updating
    /// the Lagrange multipliers)
    gsSparseMatrix<T> m_A_mat;

    /// In case Dirichlet boundary conditions have been set,
    /// contains the coefficients set
    gsMatrix<T> m_bc;

    /// The entries used for the construction of the matrix
    gsSparseEntries<T> m_entries_mat;

    /// The right hand side
    gsMatrix<T> m_rhs;

    /// The result of the linear system
    gsMatrix<T> m_result;
}; /// class _gsFittingSystem


/**
   @brief
   Class containing the sparse matrix and the RHS for gsFitting
   (TODO use gsSparseSystem?)
**/
template<class GenBasis, class T=real_t>
class gsFittingSystem : public _gsFittingSystem<GenBasis, T>
{

public:
    typedef _gsFittingSystem<GenBasis, T> Base;

    /// Constructors
    gsFittingSystem(GenBasis& basis, unsigned dim_im,
                    bool split_dim)
    : Base::_gsFittingSystem(basis, dim_im, split_dim){  }

    gsFittingSystem(GenBasis& basis, unsigned dim_im,
                    bool split_dim,
                    gsBoundaryConditions<T>& bc)
    : Base::_gsFittingSystem(basis, dim_im, split_dim, bc){  }

}; /// class _gsFittingSystem


/**
   @brief
   Class containing the sparse matrix and the RHS for gsFitting
   (TODO use gsSparseSystem?)
**/
template<class T>
class gsFittingSystem<gsBasis<T>, T>
    : public _gsFittingSystem<gsBasis<T>, T>
{

public:
    typedef _gsFittingSystem<gsBasis<T>, T> Base;

    /// Constructors
    gsFittingSystem(gsBasis<T>& basis, unsigned dim_im,
                    bool split_dim)
    : Base::_gsFittingSystem(basis, dim_im, split_dim){  }
    gsFittingSystem(gsBasis<T>& basis, unsigned dim_im,
                    bool split_dim,
                    gsBoundaryConditions<T>& bc)
    : Base::_gsFittingSystem(basis, dim_im, split_dim, bc){  }

    gsFunctionSet<T>* makeGeometry();

protected:
    using Base::m_bc;
    using Base::m_dim;
    using Base::m_result;
    using Base::m_local_global;

}; /// class _gsFittingSystem


/**
   @brief
   Class containing the sparse matrix and the RHS for gsFitting
   (TODO use gsSparseSystem?)
**/
template<class T>
class gsFittingSystem<gsMultiBasis<T>, T>
    : public _gsFittingSystem<gsMultiBasis<T>, T>
{

public:
    typedef _gsFittingSystem<gsMultiBasis<T>, T> Base;

    /// Constructors
    gsFittingSystem(gsMultiBasis<T>& basis, unsigned dim_im,
                    bool split_dim)
    : Base::_gsFittingSystem(basis, dim_im, split_dim){  }
    gsFittingSystem(gsMultiBasis<T>& basis, unsigned dim_im,
                    bool split_dim,
                    gsBoundaryConditions<T>& bc)
    : Base::_gsFittingSystem(basis, dim_im, split_dim, bc){  }

    gsFunctionSet<T>* makeGeometry();

protected:
    using Base::m_bc;
    using Base::m_dim;
    using Base::m_local_global;
    using Base::m_result;

}; /// class _gsFittingSystem

}// namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFittingSystem.hpp)
#endif
