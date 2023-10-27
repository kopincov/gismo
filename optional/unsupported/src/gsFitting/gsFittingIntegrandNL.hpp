/** @file gsFittingIntegrandNL.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsFitting/gsFittingIntegrand.h>
#include <gsFitting/gsFittingUtilsGen.h>

namespace gismo
{

template<short_t d, class T>
void gsFittingIntegrandOrth<d, T>::actualizeCurr()
{
    unsigned nBasis = m_nActives;
    unsigned nPoints = m_points->cols();
    unsigned nComp = d*(d-1)/2;
    unsigned ind2 = 0;
    unsigned ind = 0;
    T val, val1, val2;
    m_mixed_metric.setZero(nComp*nBasis*m_dim, nPoints);
    m_mixed_metric_diag.setZero(m_stride[0]*nBasis*m_dim, nPoints);
    computeMetric();

    for(unsigned basis = 0;basis < nBasis;basis++)
    {
        for(unsigned i = 0;i < m_stride[0];i++)
        {
            for(unsigned j = i+1;j < m_stride[0];j++)
            {
                for(unsigned dim = 0;dim < m_dim;dim++)
                {
                    for(unsigned pt = 0;pt < nPoints;pt++)
                    {
                        val1 = (*m_der1_curr)(dim*m_stride[0] + j, pt);
                        val2 = (*m_der1_curr)(dim*m_stride[0] + i, pt);

                        if(! m_identity_map.is_null())
                        {
                            val1 -= m_identity_map.der1_id()
                                (dim*m_stride[0] + j, pt);
                            val2 -= m_identity_map.der1_id()
                                (dim*m_stride[0] + i, pt);

                        }
                        m_mixed_metric(ind2, pt)
                            += (*m_der1)(basis*m_stride[0] + i, pt) * val1
                            + (*m_der1)(basis*m_stride[0] + j, pt) * val2;
                    }
                    ind2++;
                }
            }
        }

        for(unsigned i = 0;i < m_stride[0];i++)
        {
            for(unsigned dim = 0;dim < m_dim;dim++)
            {
                for(unsigned pt = 0;pt < nPoints;pt++)
                {
                    val = (*m_der1_curr)(dim*m_stride[0] + i, pt);
                    if(! m_identity_map.is_null())
                        val -= m_identity_map.der1_id()
                            (dim*m_stride[0] + i, pt);

                    m_mixed_metric_diag(ind, pt)
                        += (*m_der1)(basis*m_stride[0] + i, pt)
                        * val;
                }
                ind++;
            }
        }
    }
}



template<short_t d, class T>
void gsFittingIntegrandOrth<d, T>::computeMetric()
{
    unsigned nPoints = m_points->cols();
    unsigned nComp = d*(d-1)/2;
    unsigned ind = 0;
    T val, val1, val2;
    m_metric.setZero(nComp, nPoints);
    m_metric_diag.setZero(d, nPoints);

    for(unsigned i = 0;i < m_stride[0];i++)
    {
        for(unsigned j = i+1;j < m_stride[0];j++)
        {
            for(unsigned pt = 0;pt < nPoints;pt++)
            {
                for(unsigned dim = 0;dim < m_dim;dim++)
                {
                    val1 = (*m_der1_curr)(dim*m_stride[0] + i, pt);
                    val2 = (*m_der1_curr)(dim*m_stride[0] + j, pt);
                    if(! m_identity_map.is_null())
                    {
                        val1 -= m_identity_map.der1_id()
                            (dim*m_stride[0] + i, pt);
                        val2 -= m_identity_map.der1_id()
                            (dim*m_stride[0] + j, pt);
                    }
                    m_metric(ind, pt) += val1 * val2;
                }
            }
            ind++;
        }
    }
    for(unsigned i = 0;i < m_stride[0];i++)
    {
        for(unsigned pt = 0;pt < nPoints;pt++)
        {
            for(unsigned dim = 0;dim < m_dim;dim++)
            {
                val = (*m_der1_curr)(dim*m_stride[0] + i, pt);
                if(! m_identity_map.is_null())
                    val -= m_identity_map.der1_id()
                        (dim*m_stride[0] + i, pt);
                m_metric_diag(i, pt) += math::pow(val, 2);
            }
        }
    }
}

template<short_t d, class T>
T gsFittingIntegrandOrth<d, T>::compute_energy(int ind_pt)
{
    T res = 0.;
    unsigned s = d*(d-1)/2;
    for(unsigned i = 0;i < s;i++)
        res += math::pow(m_metric(i, ind_pt), 2);
    for(unsigned i = 0;i < m_stride[0];i++)
        res += math::pow(m_metric_diag(i, ind_pt), 2);
    return res;
}

template<short_t d, class T>
T gsFittingIntegrandOrth<d, T>::
compute_matrix(int base_i, int base_j, int ind_pt,
               int d1, int d2)
{
    unsigned der_ij = 0;
    T res = 0.;
    unsigned ind_i, ind_j;
    for(unsigned i = 0;i < m_stride[0];i++)
    {
        for(unsigned j = i+1;j < m_stride[0];j++)
        {
            ind_i = getPosMixed(base_i, der_ij);
            ind_j = getPosMixed(base_j, der_ij);
            res += m_mixed_metric(ind_i + d1, ind_pt)
                * m_mixed_metric(ind_j + d2, ind_pt);
            if(d1 == d2)
            {
                res += m_metric(der_ij, ind_pt) *
                    ( (*m_der1)(base_i*m_stride[0] + i, ind_pt)
                      * (*m_der1)(base_j*m_stride[0] + j, ind_pt)
                      +  (*m_der1)(base_i*m_stride[0] + j, ind_pt)
                      * (*m_der1)(base_j*m_stride[0] + i, ind_pt) );
            }
            der_ij++;
        }
        ind_i = getPosMixedDiag(base_i, i);
        ind_j = getPosMixedDiag(base_j, i);
        res += 4. * m_mixed_metric_diag(ind_i + d1, ind_pt)
            * m_mixed_metric_diag(ind_j + d2, ind_pt);
        if(d1 == d2)
        {
            res += 2. * m_metric_diag(i, ind_pt)
                * (*m_der1)(base_i*m_stride[0] + i, ind_pt)
                      * (*m_der1)(base_j*m_stride[0] + i, ind_pt);
        }
    }
    return m_coeff_orth * res;
}

template<short_t d, class T>
void gsFittingIntegrandOrth<d, T>::
compute_RHS(int base_i, int ind_pt, gsMatrix<T>& res)
{
    res.setZero(m_dim, 1);
    unsigned der_ij = 0;
    unsigned ind;
    for(unsigned i = 0;i < m_stride[0];i++)
    {
        for(unsigned j = i+1;j < m_stride[0];j++)
        {
            ind = getPosMixed(base_i, der_ij);
            for(unsigned dim = 0;dim < m_dim;dim++)
            {
                res(dim, 0) -= m_coeff_orth
                    * m_mixed_metric(ind + dim, ind_pt)
                    * m_metric(der_ij, ind_pt);
            }
            der_ij++;
        }

        ind = getPosMixedDiag(base_i, i);
        for(unsigned dim = 0;dim < m_dim;dim++)
        {
            res(dim, 0) -= 2. * m_coeff_orth
                * m_mixed_metric_diag(ind + dim, ind_pt)
                * m_metric_diag(i, ind_pt);
        }
    }
}


template<class T>
void gsFittingIntegrandW<T>::setNormals()
{
    unsigned nPoints = m_points->cols();

    m_Dx.clear(); m_Dy.clear();
    m_Dx.resize(nPoints); m_Dy.resize(nPoints);

    m_inverse_area.setZero(1, nPoints);
    m_normal.setZero(3, nPoints);
    m_squareNorm.setZero(1, nPoints);
    for(unsigned pt = 0;pt < nPoints;pt++)
    {
        m_Dx[pt].setZero();
        m_Dy[pt].setZero();
        for(unsigned dim = 0;dim < m_dim;dim++)
        {
            m_Dx[pt](dim) = (*m_der1_curr)(dim * m_stride[0], pt);
            m_Dy[pt](dim) = (*m_der1_curr)(dim * m_stride[0] + 1, pt);
        }
        m_normal.col(pt) = m_Dx[pt].cross(m_Dy[pt]);
        m_inverse_area(0, pt) = 1./m_normal.col(pt).norm();
        if(m_dim == 2 && m_normal(2, pt) < 0.)
            m_inverse_area(0, pt) *= -1;
        m_squareNorm(0, pt) = m_Dx[pt].squaredNorm()
            + m_Dy[pt].squaredNorm();
    }
}


template<class T>
void gsFittingIntegrandW<T>::actualizeCurr()
{
    unsigned nBasis = m_nActives;
    unsigned nPoints = m_points->cols();
    unsigned ind = 0;
    T Dx_basis, Dy_basis;
    gsMatrix<T> cross(3,1);

    setNormals();

    m_mixed_grad.setZero(nBasis * 3, nPoints);
    m_mixed_pv.setZero(nBasis*m_dim*3, nPoints);
    for(unsigned basis = 0;basis < nBasis;basis++)
    {
        for(unsigned dim = 0;dim < m_dim;dim++)
        {
            for(unsigned pt = 0;pt < nPoints;pt++)
            {
                cross.setZero();
                Dx_basis = (*m_der1)(basis*m_stride[0], pt);
                Dy_basis = (*m_der1)(basis*m_stride[0] + 1, pt);
                add_cross_prod(Dx_basis, dim, m_Dy[pt], cross);
                add_cross_prod(-Dy_basis, dim, m_Dx[pt], cross);
                for(unsigned dim2 = 0;dim2 < 3;dim2++)
                    m_mixed_pv(ind + dim2, pt) = cross(dim2);
            }
            ind += 3;
        }
        for(unsigned pt = 0;pt < nPoints;pt++)
        {
            for(unsigned dim2 = 0;dim2 < m_dim;dim2++)
            {
                m_mixed_grad(basis * 3 + dim2, pt)
                    = (*m_der1_curr)(dim2 * m_stride[0], pt)
                    * (*m_der1)(basis * m_stride[0], pt)
                    + (*m_der1_curr)(dim2 * m_stride[0] + 1, pt)
                    * (*m_der1)(basis * m_stride[0] + 1, pt);
            }
        }
    }
}


template<class T>
void gsFittingIntegrandW<T>::
compute_RHS(int base_i, int ind_pt, gsMatrix<T>& res)
{
    res.setZero(m_dim, 1);
    T val = m_coeff_winslow * m_squareNorm(0, ind_pt)
        * math::pow(m_inverse_area(0, ind_pt), 3);

    for(unsigned dim = 0;dim < m_dim;dim++)
    {
        res(dim, 0) -= 2. * m_coeff_winslow * m_mixed_grad
            (3 * base_i + dim, ind_pt)
            * m_inverse_area(0, ind_pt);

        res(dim, 0) += val * dot_prod_normal_mixedCross
            (base_i, dim, ind_pt);
    }
}


template<class T>
T gsFittingIntegrandW<T>::
compute_matrix(int base_i, int base_j, int ind_pt,
               int d1, int d2)
{
    T val;
    T res = 0.;
    if(d1 == d2)
    {
        val = (*m_der1)(base_i * m_stride[0], ind_pt)
            * (*m_der1)(base_j * m_stride[0], ind_pt)
            + (*m_der1)(base_i * m_stride[0] + 1, ind_pt)
            * (*m_der1)(base_j * m_stride[0] + 1, ind_pt);
        res += 2. * val * m_inverse_area(0, ind_pt);
    }

    val = 2. * m_mixed_grad(base_j * 3 + d2, ind_pt) *
        math::pow(m_inverse_area(0, ind_pt), 3);
    res -= val * dot_prod_normal_mixedCross(base_i, d1, ind_pt);

    /// same computation with i and j switched
    val = 2. * m_mixed_grad(base_i * 3 + d1, ind_pt) *
        math::pow(m_inverse_area(0, ind_pt), 3);
    res -= val * dot_prod_normal_mixedCross(base_j, d2, ind_pt);


    val = m_squareNorm(0, ind_pt)
        * math::pow(m_inverse_area(0, ind_pt), 3);
    res -= val * dot_prod_mixedCross(base_i, d1, base_j,
                                     d2, ind_pt);
    res -= val * dot_prod_normal_crossBasis(base_i, d1, base_j,
                                            d2, ind_pt);

    val = 6. * m_squareNorm(0, ind_pt)
        * math::pow(m_inverse_area(0, ind_pt), 5);
    res += val * dot_prod_normal_mixedCross(base_i, d1, ind_pt)
    * dot_prod_normal_mixedCross(base_j, d2, ind_pt);

    return m_coeff_winslow * res;
}

template<class T>
T gsFittingIntegrandW<T>::compute_energy(int ind_pt)
{
    if(m_inverse_area(0, ind_pt) < 0.)
    {
        gsWarn << "THE DETERMINANT IS NEGATIVE" << std::endl;
        return -1.;
    }
    return m_coeff_winslow * m_squareNorm(0, ind_pt)
        * m_inverse_area(0, ind_pt);
}

//////////// W3D ///////////////////

/// Computes the values of the numerator of the energy (e_1 and e_2)
template<class T>
void gsFittingIntegrandW3D<T>::computeValues()
{
    unsigned nPoints = m_points->cols();
    m_e1.setZero(nPoints, 1);
    m_e2.setZero(nPoints, 1);

    for(unsigned pt = 0;pt < nPoints;pt++)
    {
        m_e1(pt, 0) = m_metric(0, pt) * m_metric(3, pt)
            + m_metric(0, pt) * m_metric(5, pt)
            + m_metric(3, pt) * m_metric(5, pt);
        m_e2(pt, 0) = math::pow(m_metric(1, pt), 2)
            + math::pow(m_metric(2, pt), 2)
            + math::pow(m_metric(4, pt), 2);
    }
}


template<class T>
void gsFittingIntegrandW3D<T>::computeFirstDer()
{
    unsigned nBasis = m_nActives;
    unsigned nPoints = m_points->cols();
    int ind[2];
    int ind_cur;

    m_De1.setZero(3 * nBasis, nPoints);
    m_De2.setZero(3 * nBasis, nPoints);

    for(unsigned pt = 0;pt < nPoints;pt++)
    {
        for(unsigned b1 = 0;b1 < nBasis;b1++)
        {
            for(unsigned dim = 0;dim < 3;dim++)
            {
                ind_cur = 3 * b1 + dim;
                for(unsigned der = 0;der < 3;der++)
                {
                    setIndDiagMetric(der, ind);
                    m_De1(ind_cur, pt)
                        += 2. * (*m_der1_curr)(dim * m_stride[0] + der, pt)
                        * (*m_der1)(b1 * m_stride[0] + der, pt)
                        * ( m_metric(ind[0], pt) + m_metric(ind[1], pt) );

                    for(unsigned der2 = der + 1;der2 < 3;der2++)
                    {
                        m_De2(ind_cur, pt)
                            += 2. * ( (*m_der1_curr)(dim * m_stride[0] + der, pt)
                                     * (*m_der1)(b1 * m_stride[0] + der2, pt)
                                     + (*m_der1_curr)(dim * m_stride[0] + der2, pt)
                                     * (*m_der1)(b1 * m_stride[0] + der, pt) )
                            * valMetric(der, der2, pt);
                    }
                }
            }
        }
    }
}

template<class T>
void gsFittingIntegrandW3D<T>::computeSecondDer()
{
    unsigned nBasis = m_nActives;
    unsigned nPoints = m_points->cols();
    int ind[2];
    T val, val2;
    T b1_der1, b1_der2, b2_der1, b2_der2;
    int size = 9 * nBasis * nBasis;
    int pos;

    m_DDe1.setZero(size, nPoints);
    m_DDe2.setZero(size, nPoints);

    for(unsigned pt = 0;pt < nPoints;pt++)
    {
        for(unsigned b1 = 0;b1 < nBasis;b1++)
        {
            for(unsigned b2 = 0;b2 < nBasis;b2++)
            {
                for(unsigned der = 0;der < 3;der++)
                {
                    setIndDiagMetric(der, ind);
                    val = 2. * (*m_der1)(b1 * m_stride[0] + der, pt)
                        * (*m_der1)(b2 * m_stride[0] + der, pt)
                        * ( m_metric(ind[0], pt) + m_metric(ind[1], pt) );

                    val2 = 0.;
                    for(unsigned der2 = der + 1;der2 < 3;der2++)
                    {
                        val2 += 2. * ( (*m_der1)(b1 * m_stride[0] + der, pt)
                                       * (*m_der1)(b2 * m_stride[0] + der2, pt)
                                       + (*m_der1)(b1 * m_stride[0] + der2, pt)
                                       * (*m_der1)(b2 * m_stride[0] + der, pt) )
                            * valMetric(der, der2, pt);
                    }

                    for(unsigned dim = 0;dim < 3;dim++)
                    {
                        /// terms on the diagonal of the dimensions
                        pos = 9 * (b1 * nBasis + b2) + 3 * dim + dim;
                        m_DDe1(pos, pt) += val;
                        m_DDe2(pos, pt) += val2;

                        for(unsigned dim2 = 0;dim2 < 3;dim2++)
                        {
                            pos = 9 * (b1 * nBasis + b2) + 3 * dim + dim2;

                            for(unsigned der2 = der + 1;der2 < 3;der2++)
                            {
                                b1_der1 = (*m_der1_curr)(dim * m_stride[0] + der, pt)
                                    * (*m_der1)(b1 * m_stride[0] + der, pt);
                                b1_der2 = (*m_der1_curr)(dim * m_stride[0] + der2, pt)
                                    * (*m_der1)(b1 * m_stride[0] + der2, pt);
                                b2_der1 = (*m_der1_curr)(dim2 * m_stride[0] + der, pt)
                                    * (*m_der1)(b2 * m_stride[0] + der, pt);
                                b2_der2 = (*m_der1_curr)(dim2 * m_stride[0] + der2, pt)
                                    * (*m_der1)(b2 * m_stride[0] + der2, pt);
                                m_DDe1(pos, pt)
                                    += 2. * ( b1_der1 * b2_der2 + b2_der1 * b1_der2 );

                                /// The derivatives are crossed
                                b1_der2 = (*m_der1_curr)(dim * m_stride[0] + der, pt)
                                    * (*m_der1)(b1 * m_stride[0] + der2, pt);
                                b1_der1 = (*m_der1_curr)(dim * m_stride[0] + der2, pt)
                                    * (*m_der1)(b1 * m_stride[0] + der, pt);
                                b2_der2 = (*m_der1_curr)(dim2 * m_stride[0] + der, pt)
                                    * (*m_der1)(b2 * m_stride[0] + der2, pt);
                                b2_der1 = (*m_der1_curr)(dim2 * m_stride[0] + der2, pt)
                                    * (*m_der1)(b2 * m_stride[0] + der, pt);
                                m_DDe2(pos, pt)
                                    += 2. * (b1_der1 + b1_der2) * (b2_der1 + b2_der2);
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class T>
void gsFittingIntegrandW3D<T>::computeDerDeter()
{
    unsigned nPoints = m_points->cols();
    unsigned nBasis = m_nActives;
    T val;
    int dim_ij, ind, sign;
    int der3, dim3;
    int size_mixed_deter = math::pow(nBasis, 2) * 3;

    /// Will contain the derivatives of the bases
    gsMatrix<T> mat_der_b(nBasis, m_stride[0]);
    mat_der_b.setZero();

    gsMatrix<T> mat_det(m_stride[0], 3);
    mat_det.setZero();

    m_Ddeter.setZero(3 * nBasis, nPoints);
    m_DDdeter.setZero(size_mixed_deter, nPoints);

    for(unsigned pt = 0;pt < nPoints;pt++)
    {

        /// Set the derivatives of the current map and compute the determinant
        for(unsigned der = 0;der < m_stride[0];der++)
        {
            for(unsigned dim = 0;dim < m_dim;dim++)
            {
                mat_det(dim, der)
                    = (*m_der1_curr)(dim*m_stride[0] + der, pt);
            }

            for(unsigned b1 = 0;b1 < nBasis;b1++)
            {
                mat_der_b(b1, der) = (*m_der1)(b1*m_stride[0] + der, pt);
            }
        }

        /// Compute the mixed determinant
        for(unsigned b1 = 0;b1 < nBasis;b1++)
        {

            for(unsigned dim = 0;dim < 3;dim++)
            {
                for(unsigned der = 0;der < 3;der++)
                {
                    unsigned dimMj = (der + 2) % 3;
                    unsigned dimPj = (der + 1) % 3;
                    unsigned dimMi = (dim + 2) % 3;
                    unsigned dimPi = (dim + 1) % 3;

                    val = mat_det(dimMi, dimMj) * mat_det(dimPi, dimPj)
                        - mat_det(dimPi, dimMj) * mat_det(dimMi, dimPj);

                    m_Ddeter(3 * b1 + dim, pt)
                        += mat_der_b(b1, der) * val;

                    for(unsigned b2 = 0;b2 < nBasis;b2++)
                    {
                        for(unsigned dim2 = dim+1;dim2 < 3;dim2++)
                        {
                            if(dim == 0)
                                dim_ij = dim2 - 1;
                            else
                                dim_ij = 2;
                            dim3 = complementary(dim, dim2);
                            ind = 3 * (nBasis * b1 + b2) + dim_ij;
                            for(unsigned der2 = der + 1;der2 < 3;der2++)
                            {
                                der3 = complementary(der, der2);
                                sign = BaseNL::sgn(der, der2, der3)
                                    * BaseNL::sgn(dim, dim2, dim3);
                                m_DDdeter(ind, pt) +=
                                    ( mat_der_b(b1, der)
                                      * mat_der_b(b2, der2)
                                      - mat_der_b(b2, der)
                                      * mat_der_b(b1, der2) )
                                    * mat_det(der3, dim3)
                                    * sign;
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Here, the dimension is known: d=3
template<class T>
void gsFittingIntegrandW3D<T>::computeMetric()
{
    unsigned nPoints = m_points->cols();

    gsMatrix<T> mat_der(3, 3);
    int ind;

    m_metric.setZero(6, nPoints);
    m_inv_deter.setZero(nPoints, 1);
    m_deter.setZero(nPoints, 1);

    for(unsigned pt = 0;pt < nPoints;pt++)
    {

        /// Set the derivatives of the current map
        /// and compute the determinant
        for(unsigned der = 0;der < m_stride[0];der++)
        {
            for(unsigned dim = 0;dim < m_dim;dim++)
            {
                mat_der(dim, der)
                    = (*m_der1_curr)(dim*m_stride[0] + der, pt);
            }

        }
        m_deter(pt, 0) = mat_der.determinant();
        m_inv_deter(pt, 0) = 1.0 / m_deter(pt, 0);

        ind = 0;
        for(unsigned der = 0;der < m_stride[0];der++)
        {

            for(unsigned der2 = der;der2 < m_stride[0];der2++)
            {
                m_metric(ind, pt) = ( mat_der.col(der).transpose()
                                      * mat_der.col(der2) ).value();
                ind++;
            }
        }

    }
}

template<class T>
void gsFittingIntegrandW3D<T>::actualizeCurr()
{
    computeMetric();
    computeValues();
    computeFirstDer();
    computeSecondDer();
    computeDerDeter();
}

template<class T>
T gsFittingIntegrandW3D<T>::
compute_matrix(int base_i, int base_j, int ind_pt,
               int d1, int d2)
{
    unsigned nBasis = m_nActives;
    T res = 0.;
    int ind_DDe = 9 * (base_i * nBasis + base_j) + 3 * d1 + d2;
    int ind_i = 3 * base_i + d1;
    int ind_j = 3 * base_j + d2;
    T val4 = 0.;

    T val1 = ( m_DDe1(ind_DDe, ind_pt) - m_DDe2(ind_DDe, ind_pt) )
        * m_inv_deter(ind_pt, 0);
    res = val1;

    if(d1 != d2)
    {
        int ind_det, ind_der;
        int _d1, _d2;
        T sgn;
        if(d1 < d2)
        {
            sgn = 1.;
            _d1 = d1;
            _d2 = d2;
        }
        else
        {
            sgn = -1.;
            _d1 = d2;
            _d2 = d1;
        }
        if(_d1 == 0)
            ind_der = _d2 - 1;
        else
            ind_der = 2;
        ind_det = 3 * (nBasis * base_i + base_j) + ind_der;

        val4 = sgn * m_DDdeter(ind_det, ind_pt)
            * math::pow(m_inv_deter(ind_pt, 0), 2)
            * ( m_e1(ind_pt, 0) - m_e2(ind_pt, 0) );
        res -= val4;
    }
    T val2 = ( ( m_De1(ind_i, ind_pt) - m_De2(ind_i, ind_pt) )
               * m_Ddeter(ind_j, ind_pt)
               + ( m_De1(ind_j, ind_pt) - m_De2(ind_j, ind_pt) )
               * m_Ddeter(ind_i, ind_pt) )
        * math::pow(m_inv_deter(ind_pt, 0), 2);
    res -= val2;

    T val3 = 2. * math::pow(m_inv_deter(ind_pt, 0), 3)
        * m_Ddeter(ind_i, ind_pt) * m_Ddeter(ind_j, ind_pt)
        * ( m_e1(ind_pt, 0) - m_e2(ind_pt, 0) );
    res += val3;
    return res * m_coeff_winslow3D;
}

template<class T>
void gsFittingIntegrandW3D<T>::
compute_RHS(int base_i, int ind_pt, gsMatrix<T>& res)
{
    int ind_i;
    res.setZero(3,1);
    for(unsigned dim = 0;dim < 3;dim++)
    {
        ind_i = 3 * base_i + dim;
        res(dim, 0) += ( m_De2(ind_i, ind_pt)
                         - m_De1(ind_i, ind_pt) )
            * m_inv_deter(ind_pt, 0);
        res(dim, 0) += ( m_e1(ind_pt, 0)
                         - m_e2(ind_pt, 0) )
            * math::pow(m_inv_deter(ind_pt, 0), 2)
            * m_Ddeter(ind_i, ind_pt);
        res(dim, 0) *= m_coeff_winslow3D;
    }
}

} // namespace gismo
