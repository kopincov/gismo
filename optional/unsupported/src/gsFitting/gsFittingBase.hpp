/** @file gsFittingBase.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/
#pragma once

#include <gsFitting/gsParameterLines.h>

#include <gsFitting/gsFittingBase.h>
#include <gsFitting/gsFittingIdConstr.hpp>
#include <gsFitting/gsFittingUtilsGen.h>
#include <gsFitting/gsFittingUtilsGen.hpp>

#include <gsFitting/gsFittingSystem.hpp>


namespace gismo
{

template<short_t d, class GenBasis, class T>
void gsFittingBase<d, GenBasis, T>::constructId()
{
    m_identity_constr.init_identity_mapping(m_dim, m_energy);
    m_identity = m_identity_constr.get_template();
}

template<short_t d, class GenBasis, class T> void
gsFittingBase<d, GenBasis, T>::
init_gsFittingBase(bool remove_result,
                   unsigned dim_im, gsFittingParam<T>& param)
{
    m_remove_result = remove_result;
    m_result_linear = NULL;
    m_current = NULL;
    m_deform_minim = param.deform_min;
    m_dim = dim_im;
    m_identity = param.p_identity;
    m_print_messages = param.print_messages;
    m_output = param.output;
    m_is_displacement = false;

//    m_bc = NULL;
}

template<short_t d, class GenBasis, class T>
gsFittingBase<d, GenBasis, T>::
gsFittingBase(unsigned dim_im, GenBasis &basis,
              bool remove_result,
              gsFittingParam<T>& param,
              const gsPointContainer<T>& pts_LS)
: m_energy(dim_im, basis, param, pts_LS),
  m_system(basis, dim_im, false), m_basis(basis),
  m_identity_constr(basis, false, param.output, 1)
{
    init_gsFittingBase(remove_result, dim_im, param);

    if(m_identity == NULL)
        m_identity_given = false;
}

template<short_t d, class GenBasis, class T>
gsFittingBase<d, GenBasis, T>::
gsFittingBase(unsigned dim_im, GenBasis &basis,
              bool remove_result,
              gsFittingParam<T>& param)
: m_energy(dim_im, basis, param,
           gsPointContainer<T>( nBasesGen<T>(basis))),
  m_system(basis, dim_im, false), m_basis(basis),
  m_identity_constr(basis, false, param.output, 1)
{
    init_gsFittingBase(remove_result, dim_im, param);

    if(m_identity == NULL)
        m_identity_given = false;
}

template<short_t d, class GenBasis, class T>
void gsFittingBase<d, GenBasis, T>::set_properties_integrands()
{
    quadrEner().set_properties_integrands();

    if(!isLinear())
    {
        if(m_identity == NULL)
            m_identity_given = false;
    }
}


template<short_t d, class GenBasis, class T>
void gsFittingBase<d, GenBasis, T>::compute()
{
    this->set_properties_integrands();
    GISMO_ASSERT(m_result_linear == NULL, "Result already computed");
    if(m_deform_minim && m_identity == NULL)
        constructId();
    computeLinear();
}

template<short_t d, class GenBasis, class T>
void gsFittingBase<d, GenBasis, T>::computeLinear()
{
    GISMO_ASSERT(m_result_linear == NULL,
                 "The result should be NULL before the computation");
    bool split_dim = canSplitDimension();

    m_system.init_system(m_dim, split_dim);

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b
    m_energy.assemble(m_system);

    m_result_linear = m_system.solve(m_print_messages);
}

template<short_t d, class GenBasis, class T>
void gsFittingBase<d, GenBasis, T>::
set_current_map(gsFunctionSet<T>* current)
{
    if(m_current != NULL && m_current != m_identity)
        delete m_current;
    m_current = current;

    energy().set_current_map(current);
}

template<short_t d, class GenBasis, class T>
void gsFittingBase<d, GenBasis, T>::
actualize_basis()
{
    m_system.actualize_mapper(true);
    m_energy.actualize_basis();
}

/*
template<short_t d, class GenBasis, class T>
void gsFittingBase<d, GenBasis, T>::
add_fixed_border(gsMultiPatch<T>& border,
                 gsMultiPatch<T>& temp_border,
                 std::vector<int>& ind)
{
    bool new_bc = (m_bc == NULL);
    if(new_bc)
        m_bc = new gsBoundaryConditions<T>();
    unsigned s = ind.size();
    std::vector<T> init(d-1);
    std::vector<T> end(d-1);
    std::fill(init.begin(), init.end(), 0.);
    std::fill(end.begin(), end.end(), 1.);

    int del = 1 + 2*d;
    for(unsigned i = 0;i < s;i += del)
    {
        boxSide s(ind[i+1]);
        int ind_patch = ind[i];
        gsBasis<T>& basis = getBasisGen<T>(m_basis, ind_patch);
        for(unsigned j = 0;j < d-1;j++)
        {
            init[j] = ind[i+2+j];
            end[j] = ind[i+1+d+j];
        }

        gsGeometry<T>& patch_geom = border.patch(ind[i+2]);
        gsGeometry<T>& patch_temp = temp_border.patch(ind[i+2]);
        gsGeometry<T>* new_geom
            = basisProjectionInterpolation<T>
            (s, basis, patch_geom, patch_temp, init, end);
        m_bc->addCondition(ind[i], ind[i+1],
                          condition_type::dirichlet, new_geom);
    }
    if(new_bc)
        m_system.set_bc(m_bc);
}

template<short_t d, class GenBasis, class T>
void gsFittingBase<d, GenBasis, T>::
add_fixed_border(std::vector<gsMatrix<T> >& pts,
                 std::vector<gsMatrix<T> >& param,
                 std::vector<int>& ind)
{

    bool new_bc = (m_bc == NULL);
    if(new_bc)
        m_bc = new gsBoundaryConditions<T>();
    int s = ind.size();
    std::vector<T> init(d-1);
    std::vector<T> end(d-1);
    std::fill(init.begin(), init.end(), 0.);
    std::fill(end.begin(), end.end(), 1.);


    int del = 1 + 2*d;
    for(unsigned i = 0;i < s;i += del)
    {
        boxSide s(ind[i+1]);
        int ind_patch = ind[i];
        gsBasis<T>& basis = getBasis(m_basis, ind_patch);
        for(unsigned j = 0;j < d-1;j++)
        {
            init[j] = ind[i+2+j];
            end[j] = ind[i+1+d+j];
        }

         gsGeometry<T>* new_geom
            = basisProjectionInterpolation<T>
             (s, m_basis, pts[i], param[i], init, end);
        m_bc.addCondition(ind[i], ind[i+1],
                          condition_type::dirichlet, new_geom);
    }
    if(new_bc)
        m_system.set_bc(m_bc);
        }*/

} // namespace gismo
