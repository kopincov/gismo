/** @file gsLocalGlobal.h

    @brief Contains the class gsLocalGlobal that permits to call
    the mapper whenever needed (case multipatch)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsDofMapper.h>

#include <gsFitting/gsFittingIterator.h>

namespace gismo
{
template<class GenBasis, class T> class gsLocalGlobal;

/**
   Simple class to call the mapper in case of a multipatch
*/
template<class GenBasis, class T>
class _gsLocalGlobal
{
public:
    _gsLocalGlobal(GenBasis& basis, int dim,
                   gsBoundaryConditions<T>* bc)
    : m_basis(basis)
    {
        m_dim = dim;
        m_bc = bc;
    }

    virtual ~_gsLocalGlobal(){  }

    inline int size_freeDof()
    {  return m_mapper.freeSize();  }

    inline int size_totalDOF()
    {  return m_mapper.size();  }

    inline bool is_free(int ind, int patch)
    {
        int global_ind = m_mapper.index(ind, patch);
        return m_mapper.is_free_index(global_ind);
    }

    inline int add_split_dimension(int global_ind,
                                   bool split_dim)
    {
        if(split_dim)
            return global_ind;
        else
            return m_dim*global_ind;
    }

    void set_bc(gsBoundaryConditions<T>* bc)
    {
        m_bc = bc;
        gsLocalGlobal<GenBasis, T>* _this
            = static_cast<gsLocalGlobal<GenBasis, T>*>(this);
        _this->actualize_mapper(false);
    }

    /// Returns the global index of a degree of freedom
    /// associated with the index at some patch
    int localToGlobal(int ind, int patch,
                      bool split_dim)
    {
        return this->add_split_dimension
            (m_mapper.index(ind, patch),
             split_dim);
    }

    int bindex(int ind, int patch)
    {   return m_mapper.bindex(ind, patch);   }

    GenBasis& basis()
    {  return m_basis;  }
    gsBoundaryConditions<T>* get_bc()
    {  return m_bc;  }

protected:
    /// The basis (used for adapting the mapping if necessary,
    /// i.e., if it has been modified)
    GenBasis& m_basis;

    /// The mapper used to compute the global index
    /// of a degree of freedom
    gsDofMapper m_mapper;

    /// The dimension of the image
    /// (used in case the dimensions are not split)
    int m_dim;

    /// The Dirichlet boundary conditions
    gsBoundaryConditions<T> *m_bc;
}; // class _gsLocalGlobal


template<class GenBasis, class T>
class gsLocalGlobal : public _gsLocalGlobal<GenBasis, T>
{
    typedef _gsLocalGlobal<GenBasis, T> Base;

public:
    gsLocalGlobal(GenBasis& basis, int dim,
                   gsBoundaryConditions<T>* bc = NULL) :
    Base::_gsLocalGlobal(basis, dim, bc) {  }


    inline int nPatches();

    void actualize_mapper(bool repair);
}; // class gsLocalGlobal


template<class T>
class gsLocalGlobal<gsBasis<T>, T>
    : public _gsLocalGlobal<gsBasis<T>, T>
{
    typedef _gsLocalGlobal<gsBasis<T>, T> Base;
    using Base::m_basis;
    using Base::m_mapper;
    using Base::m_bc;

public:
    gsLocalGlobal(gsBasis<T>& basis, int dim,
                  gsBoundaryConditions<T>* bc = NULL) :
    Base::_gsLocalGlobal(basis, dim, bc)
    { actualize_mapper(false);  }

    ~gsLocalGlobal(){  }

    inline int nPatches(){ return 1; }

    void actualize_mapper(bool repair)
    {
        if(m_bc == NULL) {
            m_mapper = gsDofMapper(m_basis);
            m_mapper.finalize();
        }
        else {
            gsMultiBasis<T> tmp(m_basis);
            tmp.getMapper(true, *m_bc, m_mapper);
        }
    }
}; // class gsLocalGlobal


template<class T>
class gsLocalGlobal<gsMultiBasis<T>, T>
    : public _gsLocalGlobal<gsMultiBasis<T>, T>
{
    typedef _gsLocalGlobal<gsMultiBasis<T>, T> Base;
    using Base::m_basis;
    using Base::m_mapper;
    using Base::m_bc;

public:
    gsLocalGlobal(gsMultiBasis<T>& basis, int dim,
                  gsBoundaryConditions<T>* bc = NULL) :
    Base::_gsLocalGlobal(basis, dim, bc)
    {  actualize_mapper(false); }

    ~gsLocalGlobal(){  }


    inline int nPatches(){ return m_basis.nPatches(); }

    void actualize_mapper(bool repair)
    {
        if(repair)
            m_basis.repairInterfaces
                (m_basis.topology().interfaces());
        if(m_bc == NULL)
            m_basis.getMapper(true, m_mapper);
        else
            m_basis.getMapper(true, *m_bc, m_mapper);
    }

}; // class gsLocalGlobal

}
