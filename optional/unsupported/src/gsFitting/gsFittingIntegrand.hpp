/** @file gsFittingIntegrand.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingIntegrand.h>
#include <gsFitting/gsFittingUtilsGen.h>

namespace gismo
{

template<short_t d, class T>
gsFittingIntegrand<d, T>::
gsFittingIntegrand(bool is_smoothing, bool split_dim,
                   bool isLinear)
{
    m_stride[0] = d; m_stride[1] = d*(d+1)/2;
    m_is_smoothing = is_smoothing;
    m_split_dim = split_dim;
    m_isLinear = isLinear;
    m_dim = 0;

    m_image = NULL;
    m_der1 = NULL;
    m_der2 = NULL;
    m_points = NULL;
    m_points_param_space = NULL;

    m_image_curr = NULL;
    m_der1_curr = NULL;
    m_der2_curr = NULL;
}

template<short_t d, class T>
void gsFittingIntegrand<d, T>::
setLocalValues(gsMatrix<T>* points,
               gsMatrix<T>* points_param_space,
               gsMatrix<T>* image,
               unsigned nActives, unsigned dim_im,
               gsMatrix<T>* der1,
               gsMatrix<T>* der2,
               unsigned patch,
               bool have_current_mapping)
{
    m_image = image;
    m_der1 = der1;
    m_der2 = der2;

    m_dim = dim_im;
    m_nActives = nActives;
    m_patch = patch;
    m_points = points;
    m_points_param_space = points_param_space;
    m_have_current_mapping = have_current_mapping;
    actualize();
}

template<short_t d, class T>
void gsFittingIntegrand<d, T>::
setLocalValuesCurr(gsMatrix<T>* image, gsMatrix<T>* der1,
                   gsMatrix<T>* der2)
{
    m_image_curr = image;
    m_der1_curr = der1;
    m_der2_curr = der2;
    m_have_current_mapping = true;
    actualizeCurr();
}

template<short_t d, class T>
void gsFittingIntegrand<d, T>::
setLocalValuesError(unsigned m_dim,
                    gsMatrix<T>* points,
                    gsMatrix<T>* points_param_space,
                    gsMatrix<T>* image,
                    gsMatrix<T>* der1,
                    gsMatrix<T>* der2,
                    unsigned patch)
{
    m_image_curr = image;
    m_der1_curr = der1;
    m_der2_curr = der2;
    m_patch = patch;
    m_points = points;
    m_points_param_space = points_param_space;
    this->actualize_energy();
}


template<short_t d, class GenGeom, class T>
void gsFittingIntegrandL2Dist<d, GenGeom, T>::actualize()
{
    if(m_f != NULL)
    {
        m_image_f.setZero(m_dimAmbient, m_points->cols());
        fitting_eval_into<T>(*m_f, *m_points_param_space,
                             m_image_f, m_patch);
    }
}

template<short_t d, class GenGeom, class T>
T gsFittingIntegrandL2Dist<d, GenGeom, T>::
compute_matrix(int base_i, int base_j, int ind_pt,
               int d1, int d2)
{
    if(d1 == d2)
        return m_coeff * (*m_image)(base_i,ind_pt)
            *(*m_image)(base_j,ind_pt);
    else
        return 0.;
}

template<short_t d, class GenGeom, class T>
void gsFittingIntegrandL2Dist<d, GenGeom, T>::
compute_RHS(int base_i, int ind_pt, gsMatrix<T>& res)
{
    if(m_f != NULL)
    {
        res = m_coeff * (*m_image)(base_i,ind_pt)*m_image_f.col(ind_pt);
        if(m_have_current_mapping)
        {
            res.noalias() -= m_coeff * (*m_image)(base_i,ind_pt)
                * m_image_curr->col(ind_pt);
        }
    }
    else
        res.setZero(m_dim, 1);
}

template<short_t d, class GenGeom, class T>
T gsFittingIntegrandL2Dist<d, GenGeom, T>::
compute_energy(int ind_pt)
{
    if(m_f == NULL)
        return 0.;
    gsVector<real_t> res = (*m_image_curr).col(ind_pt)
        - m_image_f.col(ind_pt);
    return res.squaredNorm();
}

template<short_t d, unsigned o, class T>
T gsFittingIntegrandLin<d, o, T>::
computeHess(int base_i, int base_j, int ind_pt)
{
    T res = 0; // temporary variable
    T val;
    for (unsigned s = 0; s < m_stride[1]; s++)
    {
        val = (*m_der2)(base_i * m_stride[1] + s, ind_pt) *
            (*m_der2)(base_j * m_stride[1] + s, ind_pt);
        if (s < d)
            res += val;  // d^2u N_i * d^2u N_j + ...
        else
            res += 2.*val; // Mixed derivatives
    }
    return res;
}

template<short_t d, unsigned o, class T>
T gsFittingIntegrandLin<d, o, T>::
computeGrad(int base_i, int base_j, int ind_pt)
{
    T res = 0; // temporary variable

    for (unsigned s = 0; s < m_stride[0]; s++)
    {
        res += (*m_der1)(base_i*m_stride[0] + s, ind_pt)
            * (*m_der1)(base_j*m_stride[0] + s, ind_pt);
    }
    return res;
}

template<short_t d, unsigned o, class T>
T gsFittingIntegrandLin<d, o, T>::
compute_matrix(int base_i, int base_j, int ind_pt, int d1, int d2)
{
    if(d1 == d2)
    {
        if(o < 2)
            return m_coeff_grad * computeGrad(base_i, base_j, ind_pt);
        else
            return m_coeff_grad * computeGrad(base_i, base_j, ind_pt)
                + m_coeff_hess * computeHess(base_i, base_j, ind_pt);
    }
    else   return 0.;
}


template<unsigned o, class T> void
gsFittingId<o, T>::actualize(int patch, gsMatrix<T>& points)
{
    if(m_map != NULL)
    {
        gsGeometry<T>* _template =
            dynamic_cast<gsGeometry<T>*>(m_map);
        if(_template != NULL)
        {
            fitting_eval_into<T>(*_template, points,
                                  m_image_id, patch);
            if(o > 0)
                fitting_deriv_into<T>(*_template, points,
                                      m_der1_id, patch);
        } else {
            gsMultiPatch<T>* patches =
                dynamic_cast<gsMultiPatch<T>*>(m_map);
            GISMO_ASSERT(patches != NULL, "ERROR the object is neither a gsGeometry nor a gsMultiPatch");
            fitting_eval_into<T>(*patches, points,
                                  m_image_id, patch);
            if(o > 0)
                fitting_deriv_into<T>(*patches, points,
                                      m_der1_id, patch);
        }
    }
}


template<short_t d, unsigned o, class T>
T gsFittingIntegrandLin<d, o, T>::energyGrad(int ind_pt)
{
    gsMatrix<T> del = (*m_der1_curr).col(ind_pt);
    if(! m_identity_map.is_null())
        del -= m_identity_map.der1_id().col(ind_pt);

    return del.squaredNorm();
}

template<short_t d, unsigned o, class T>
T gsFittingIntegrandLin<d, o, T>::energyHess(int ind_pt)
{
    T res = 0; // temporary variable
    T val;
    for(unsigned i = 0;i < m_dim;i++)
    {
        for (unsigned s = 0; s < m_stride[1]; s++)
        {
            val = math::pow((*m_der2_curr)(i*m_stride[1] + s, ind_pt), 2);
            if (s < d)
                res += val;   // d^2u N_i * d^2u N_j + ...
            else
                res += 2.*val;  // Mixed derivatives
        }
    }
    return res;
}

template<short_t d, unsigned o, class T>
void gsFittingIntegrandLin<d, o, T>::
RHS_grad(int base_i, int ind_pt, gsMatrix<T>& res)
{
    T val;
    gsMatrix<T> tmp(m_dim, 1);
    tmp.setZero();
    for(unsigned i = 0;i < m_dim;i++)
    {
        for(unsigned s = 0;s < m_stride[0];s++)
        {
            val = (*m_der1_curr)(i*m_stride[0] + s, ind_pt);
            if(! m_identity_map.is_null())
                val -= m_identity_map.der1_id()
                    (i*m_stride[0] + s, ind_pt);
            tmp(i,0) -= (*m_der1)(base_i*m_stride[0] + s, ind_pt)
                * val;
        }
    }
    res += tmp*m_coeff_grad;
}

template<short_t d, unsigned o, class T>
void gsFittingIntegrandLin<d, o, T>::
RHS_grad_id(int base_i, int ind_pt, gsMatrix<T>& res)
{
    gsMatrix<T> tmp(m_dim, 1);
    tmp.setZero();
    for(unsigned i = 0;i < m_dim;i++)
    {
        for(unsigned s = 0;s < m_stride[0];s++)
        {
            tmp(i,0) += (*m_der1)(base_i*m_stride[0] + s, ind_pt)
                * m_identity_map.der1_id()
                (i*m_stride[0] + s, ind_pt);
        }
    }
    res += tmp*m_coeff_grad;
}

template<short_t d, unsigned o, class T>
void gsFittingIntegrandLin<d, o, T>::
RHS_hess(int base_i, int ind_pt, gsMatrix<T>& res)
{
    T val;
    gsMatrix<T> tmp(m_dim, 1);
    tmp.setZero();
    for(unsigned i = 0;i < m_dim;i++)
    {
        for (unsigned s = 0; s < m_stride[1]; s++)
        {
            val = (*m_der2)(base_i*m_stride[1] + s, ind_pt)
                * (*m_der2_curr)(i*m_stride[1] + s, ind_pt);
            if (s < d)
                tmp(i,0) -= val;      // d^2u N_i * d^2u N_j + ...
            else
                tmp(i,0) -= 2.*val;   // Mixed derivatives
        }
    }
    res += tmp*m_coeff_hess;
}

template<short_t d, unsigned o, class T>
gsFittingIntegrandLin<d, o, T>* gsFittingIntegrandLin<d, o, T>::copy()
{
    return new gsFittingIntegrandLin<d, o, T>
        (m_coeff_linear, m_coeff_grad/m_coeff_linear,
         m_coeff_hess/m_coeff_linear, m_isLinear);
}

template<short_t d, class GenGeom, class T>
gsFittingIntegrandL2Dist<d, GenGeom, T>*
gsFittingIntegrandL2Dist<d, GenGeom, T>::copy()
{
    return new gsFittingIntegrandL2Dist
        <d, GenGeom, T>(m_dimAmbient, m_f);
}


template<short_t d, class T>
T gsFittingIntegrandL2Norm<d, T>::compute_energy(int ind_pt)
{
    return (*m_image_curr).col(ind_pt).squaredNorm();
}


} // namespace gismo
