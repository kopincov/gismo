/** @file gsFittingQuadrature.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsFitting/gsParameterLines.h>

#include <gsFitting/gsFittingIntegrand.h>
#include <gsFitting/gsFittingIntegrand.hpp>
#include <gsFitting/gsFittingIntegrandNL.hpp>
#include <gsFitting/gsFittingIdConstr.hpp>


namespace gismo
{

/// If m_split_dim is true, then entries_mat is associated with
/// a matrix of size N*N, with N the number of degrees of freedom
/// Otherwise, A_mat is a matrix of size (d*N)*(d*N), with d the dimension and
/// N the number of degrees of freedom
template<short_t d, class GenBasis, unsigned o, class T>
void gsFittingQuadrature<d, GenBasis, o, T>::
assembleSystem(gsFittingSystem<GenBasis, T>& system, bool isLinear)
{
    bool split_dim = system.split_dim();
// #define DBG_EXPORT_GAUSS_POINTS
    unsigned s = m_integrand.size();
    GISMO_ASSERT(s > 0, "Error: gsFittingQuadrature constructed without any integrand!!");

    gsVector<int> numNodes(d);
    for ( int i = 0; i!= d; ++i )
        numNodes[i] = m_size_quadr[i];
    gsGaussRule<T> QuRule( numNodes ); // Reference Quadrature rule
    gsMatrix<T> quNodes, quNodes_param;
    gsVector<T> quWeights;
    std::vector< gsMatrix<T> > localA;

    /// Dimensions used in case where we cannot split the dimensions
    unsigned _d, _s;

    /// Construct the "generic template iterator".
    gsGenGeomIterator<T>* domIt = generic_template_iterator();

    /// If we cannot split the dimensions, we need to
    /// fill m_dim*m_dim matrices at each cell
    if(split_dim)
    {
        /// We only need one matrix in the case where the dimensions
        /// can be considered independently
        localA.push_back(gsMatrix<T>());
        _d = 1;
        _s = 1;

    } else {
        _s = m_dim*m_dim;
        _d = m_dim;
        for(unsigned i = 0;i < _s;i++)
            localA.push_back(gsMatrix<T>());
    }

    gsMatrix<T> localB;
    gsMatrix<T> B_tmp;
    gsMatrix<T> centerPointDom;
    unsigned d1_2;
    int patch;

#ifdef DBG_EXPORT_GAUSS_POINTS
    /// DBG: export the points of the quadrature
    std::string file = "exp_gauss_points_" + util::to_string(d);
    std::vector< gsVector<T> > dbg_pts;
#endif

    GISMO_ENSURE(isLinear || m_current != NULL,
                 "minimize a nonlinear energy without current solution !!");

    for (; domIt->good(); domIt->next() )
    {
        patch = domIt->patch();

        /// Map the Quadrature rule to the element
        /// and compute basis derivatives
        compute_quadrature_points(*domIt, QuRule, quNodes_param,
                                  quNodes, quWeights);

        gsBasis<T>& cur_basis =
            domIt->currentDomainBasis(m_basis);
        centerPointDom = domIt->imageCenterPoint();
        compute_basis(quNodes, cur_basis, centerPointDom);

        if(m_current != NULL)
        {
            gsGeometry<T>& cur_geom =
                domIt->currentDomainGeom(m_basis, *m_current);
            compute_current(quNodes, cur_geom);
        }

        const index_t numActive = m_actives.rows();

        for(d1_2 = 0;d1_2 < _s;d1_2++)
            localA[d1_2].setZero(numActive, numActive);

        localB.setZero(m_dim, numActive);
        B_tmp.setZero(m_dim, 1);

        for (unsigned _integr = 0;_integr < s;_integr++)
        {
            m_integrand[_integr]
                ->setLocalValues(&quNodes, &quNodes_param,
                                 &m_image, numActive, m_dim, &m_der1,
                                 &m_der2, domIt->patch());
            if(m_current != NULL)
            {
                m_integrand[_integr]->setLocalValuesCurr
                    (&m_image_curr, &m_der1_curr, &m_der2_curr);
            }
        }

        /// perform the quadrature and compute the terms that
        /// must be added to the matrix and the RHS
        for (index_t k = 0; k < quWeights.rows(); ++k)
        {
#ifdef DBG_EXPORT_GAUSS_POINTS
            dbg_pts.push_back(quNodes.col(k));
#endif
            const T weight = quWeights[k];

            for (unsigned _integr = 0;_integr < s;_integr++)
            {
                for (index_t i = 0; i != numActive; ++i)
                {
                    m_integrand[_integr]
                        ->compute_RHS(i, k, B_tmp);
                    localB.col(i) += weight * B_tmp;
                    for (index_t j = 0; j != numActive; ++j)
                    {
                        /// _d is set to 1 whenever we can split
                        /// dimensions and d otherwise
                        for(unsigned d1 = 0;d1 < _d;d1++)
                        {
                            for(unsigned d2 = 0;d2 < _d;d2++)
                            {
                                d1_2 = d1*m_dim + d2;
                                localA[d1_2](i, j) +=
                                    weight * m_integrand[_integr]
                                    ->compute_matrix(i, j, k, d1, d2);
                            }
                        }
                    }
                }
            }
        }
        /// add the new terms to the matrix and to the RHS
        for (index_t i = 0; i != numActive; ++i)
        {
            system.add_rhs(m_actives(i, 0), patch, localB.col(i));

            for (index_t j = 0; j != numActive; ++j)
                add_matrix(system, i, j, m_actives(i, 0),
                           m_actives(j, 0), patch, localA);
        }
    }
    delete domIt;

    /// DBG
#ifdef DBG_EXPORT_GAUSS_POINTS
    gsMatrix<T> exp(m_dim, dbg_pts.size());
    for(unsigned i = 0;i < dbg_pts.size();i++)
        exp.col(i) = dbg_pts[i];
    gsWriteParaviewPoints(exp, file);
#endif
}

template<short_t d, class GenBasis, unsigned o, class T>
void gsFittingQuadrature<d, GenBasis, o, T>::
compute_basis(gsMatrix<T>& pts,
              gsBasis<T>& basis,
              gsMatrix<T>& centerPoint)
{
    basis.active_into(centerPoint, m_actives);
    basis.eval_into(pts, m_image);
    if(o > 0)
    {
        basis.deriv_into(pts, m_der1);
        if(o > 1)
            basis.deriv2_into(pts, m_der2);
    }
}

template<short_t d, class GenBasis, unsigned o, class T>
void gsFittingQuadrature<d, GenBasis, o, T>::
compute_current(gsMatrix<T>& pts,
               gsGeometry<T>& geom)
{
    geom.eval_into(pts, m_image_curr);
    if(o > 0)
    {
        geom.deriv_into(pts, m_der1_curr);
        if(o > 1)
            geom.deriv2_into(pts, m_der2_curr);
    }
}


template<short_t d, class GenBasis, unsigned o, class T>
void gsFittingQuadrature<d, GenBasis, o, T>::
set_properties_integrand()
{
    m_split_dim = true;
    m_isLinear = true;
    unsigned s = m_integrand.size();
    if(s > 0)
    {
        for(unsigned i = 0;i < s;i++)
        {
            m_isLinear = m_isLinear
                && m_integrand[i]->isLinear();
            m_split_dim = m_split_dim
                && m_integrand[i]->canSplitDimension();
        }
    }
}

template<short_t d, class GenBasis, unsigned o, class T>
void gsFittingQuadrature<d, GenBasis, o, T>::computeEnergies(int type)
{

    gsVector<int> numNodes(d);
    for ( int i = 0; i!= d; ++i )
        numNodes[i] = std::max(m_basis.degree(i), (short_t)2);//+1;
    gsGaussRule<T> QuRule( numNodes ); // Reference Quadrature rule
    gsMatrix<T> quNodes, quNodes_param;
    gsVector<T> quWeights;
    unsigned s = m_integrand.size();
    //  gsMatrix<T> centerPoint_param, centerPoint;

    /// Construct the "generic template iterator".
    gsGenGeomIterator<T>* domIt = generic_template_iterator();
    T error_cell, smoothingEner_cell, value;
    T _error_cell, _smoothingEner_cell;

    int ind, patch, prec_patch = 0;

    reset_error();

    GISMO_ENSURE(m_current != NULL, "compute the error with a result null !!");
    ind = 0;
    for (; domIt->good(); domIt->next() )
    {
        /// Map the Quadrature rule to the element
        /// and compute basis derivatives
        compute_quadrature_points(*domIt, QuRule, quNodes_param,
                                  quNodes, quWeights);
        gsGeometry<T>& cur_geom =
            domIt->currentDomainGeom(m_basis, *m_current);
        compute_current(quNodes, cur_geom);

        /// m_image, m_der1, m_der2 are computed by the integrands
        for (unsigned _integr = 0;_integr < s;_integr++)
        {
            m_integrand[_integr]
                ->setLocalValuesError(m_dim, &quNodes,
                                      &quNodes_param,
                                      &m_image_curr,
                                      &m_der1_curr, &m_der2_curr,
                                      domIt->patch());
        }

        /// perform the quadrature and compute the error
        error_cell = 0.;
        smoothingEner_cell = 0.;
        for (index_t k = 0; k < quWeights.rows(); ++k)
        {
            const T weight = quWeights[k];

            for (unsigned _integr = 0;_integr < s;_integr++)
            {
                value = weight * m_integrand[_integr]
                    ->compute_energy(k);
                if(value < 0.)
                {
                    GISMO_ASSERT(m_integrand[_integr]->is_smoothing(),
                        "-1 in an energy of type fitting. This should not happend!!");
                    m_maxEnergySmoothing = value;
                    m_minEnergySmoothing = value;
                    m_totalEnergySmoothing = value;
                    delete domIt;
                    return;
                }
                if(! m_integrand[_integr]->is_smoothing())
                    error_cell += value;
                else
                    smoothingEner_cell += value;
            }
        }
        _error_cell = set_local_value(type, error_cell, 1);
        _smoothingEner_cell = set_local_value(type,
                                              smoothingEner_cell, 0);

        patch = domIt->patch();
        if(prec_patch != patch)
            ind = 0;
        gsMatrix<T>& tmp_err = m_errors[patch];
        gsMatrix<T>& tmp_energy = m_energySmoothing[patch];
        tmp_err(ind,0) = _error_cell;
        tmp_energy(ind,0) = _smoothingEner_cell;
        prec_patch = patch;
        ind++;
    }
    delete domIt;
}


template<short_t d, class GenBasis, unsigned o, class T>
T gsFittingQuadrature<d, GenBasis, o, T>::
set_local_value(int type, T _err, bool is_error)
{
    T err = 0.;

    switch (type)
    {
    case 0:
        err = _err;
        break;
    case 1:
        err = sqrt(_err);
        break;
    default:
        gsWarn << "Unknown type in get_Error(errors, type)...\n";
        break;
    }

    if(is_error)
    {
        if ( err > m_maxError ) m_maxError = err;
        if ( err < m_minError ) m_minError = err;
        m_totalError += err;
    } else {
        if ( err > m_maxEnergySmoothing ) m_maxEnergySmoothing = err;
        if ( err < m_minEnergySmoothing ) m_minEnergySmoothing = err;
        m_totalEnergySmoothing += err;
    }
    return err;
}


template<short_t d, class GenBasis, unsigned o, class T>
void gsFittingQuadrature<d, GenBasis, o, T>::reset_error()
{
    gsGenGeomIterator<T>* domIt = generic_template_iterator();
    m_errors.clear();
    m_energySmoothing.clear();

    unsigned s = domIt->nPatches();
    for(unsigned i = 0;i < s;i++)
    {
        m_errors.push_back(gsMatrix<T>(domIt->numElement(i),1));
        m_energySmoothing.push_back
            (gsMatrix<T>(domIt->numElement(i),1));
    }
    m_maxError = 0;
    m_maxEnergySmoothing = 0;
    m_minError = std::numeric_limits<double>::max();
    m_minEnergySmoothing = std::numeric_limits<double>::max();
    m_totalEnergySmoothing = 0;
    m_totalError = 0;
    delete domIt;
}


template<short_t d, class GenBasis, unsigned o, class T>
bool gsFittingQuadrature<d, GenBasis, o, T>::
copy_non_smoothing_integrand(gsFittingQuadrature<d, GenBasis, o, T>& origin,
                             gsFittingQuadrature<d, GenBasis, o, T>& dest)
{
    unsigned s = origin.m_integrand.size();
    bool res = false;
    for(unsigned i =0;i < s;i++)
    {
        if(! origin.m_integrand[i]->is_smoothing())
        {
            dest.m_integrand.
                push_back(origin.m_integrand[i]->copy() );
            res = true;
        }
    }
    if(res)
        dest.set_properties_integrand();
    return res;
}

template<short_t d, class GenBasis, unsigned o, class T>
void gsFittingQuadrature<d, GenBasis, o, T>::
compute_quadrature_points(gsGenGeomIterator<T>& domIt,
                          gsGaussRule<T>& QuRule,
                          gsMatrix<T>& quNodes_param,
                          gsMatrix<T>& quNodes,
                          gsVector<T>& quWeights)
{
    QuRule.mapTo(domIt.lowerCorner(), domIt.upperCorner(),
                 quNodes_param, quWeights);
    domIt.eval_into(quNodes_param, quNodes);
}



template<short_t d, class GenBasis, unsigned o, class T>
gsGenGeomIterator<T>* gsFittingQuadrature<d, GenBasis, o, T>::
generic_template_iterator()
{
    if(m_domain == NULL)
    {
        return new gsGenGeomIterator<T>(m_basis);
    }
    else
    {
        gsGeometry<T>* geom_tmp
            = dynamic_cast<gsGeometry<T>*>(m_domain);
        if(geom_tmp != NULL)
            return new gsGenGeomIteratorSimplePatch<T>(*geom_tmp);
        else
        {
            gsMultiPatch<T>* patches_tmp
                = dynamic_cast<gsMultiPatch<T>*>(m_domain);
            if(patches_tmp != NULL)
                return new gsGenGeomIteratorMultiPatch
                    <T>(*patches_tmp);
            else
            {
                gsWarn << "Error during the quadrature (gsFitting)" << std::endl;
                return NULL;
            }
        }
    }
}

template<short_t d, class GenBasis, unsigned o, class T>
gsFittingQuadrature<d, GenBasis, o, T>&
gsFittingQuadrature<d, GenBasis, o, T>::
operator= ( const gsFittingQuadrature<d, GenBasis, o, T>& other )
{
    m_basis = other.m_basis;
    m_domain = other.m_domain;
    init_size_quadr(other.m_size_quadr);
    removeAll();
    m_integrand = other.m_integrand;
    m_split_dim = other.m_split_dim;
    m_isLinear = other.m_isLinear;
    m_current = other.m_current;
    m_dim = other.m_dim;
    return *this;
}

} // namespace gismo
