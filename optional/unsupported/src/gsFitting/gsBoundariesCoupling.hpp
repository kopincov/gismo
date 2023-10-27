/** @file gsBoundariesCoupling.hpp

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsBoundariesCoupling.h>
#include <gsFitting/gsFittingUtilsGen.hpp>

namespace gismo
{

template<class T> gsFittingBorder<T>::
gsFittingBorder(const gsFittingBorder<T>& source)
: m_points(source.m_points),
      m_strong(source.m_strong), m_type(source.m_type)
{
    init();
    m_geometry = source.m_geometry;
    if(m_type == type_border::trim_surf)
        m_trimmed = source.m_trimmed;
}

template<class T>
gsFittingBorderCoupling<T>::
gsFittingBorderCoupling(gsFittingBorder<T>& cond_template,
                        gsFittingBorder<T>& cond_target)
: m_cond_template(cond_template), m_cond_target(cond_target),
  m_type(type_coupling::unknownType)
{
    if(cond_target.type() == cond_template.type())
    {
        if(cond_target.type() == type_border::pts)
            m_type = type_coupling::PtsPts;
        else if(cond_target.type() == type_border::geom)
            m_type = type_coupling::GeomGeom;
        else if(cond_target.type() == type_border::trim_surf)
            m_type = type_coupling::TrimTrim;
    }
    else if(cond_template.type() == type_border::side)
    {
        if(cond_target.type() == type_border::geom)
            m_type = type_coupling::StrongGeom;
        else if(cond_target.type() == type_border::pts)
            m_type = type_coupling::StrongPts;
    }
    if(m_type == type_coupling::unknownType)
        gsWarn << "Coupling unknown" << std::endl;
}

template<class T> void
gsFittingBorder<T>::getBoundingBox(gsMatrix<T>& bbox)
{
    GISMO_ASSERT(m_type != type_border::side,
                 "Should not be called with side");
    if(m_type == type_border::pts)
        getBoundingBoxPts<T>(m_points, bbox);
    else if(m_type == type_border::geom)
        setBoundingBox<T>(*m_geometry, bbox);
    else if(m_type == type_border::trim_surf)
    {
        setBoundingBox<T>(*m_trimmed.getTP(), bbox);
    }
}

template<class T> index_t
gsFittingBorder<T>::dim_im()
{
    GISMO_ASSERT(m_type != type_border::side,
                 "Should not be called with side");
    if(m_type == type_border::pts)
        return m_points.rows();
    else if(m_type == type_border::geom)
        return m_geometry->geoDim();
    else if(m_type == type_border::trim_surf)
        return m_trimmed.geoDim();
    return 0;
}

template<class T> template<class GenBasis>
void gsFittingBorderCoupling<T>::
initialize(gsFittingParam<T>& param, GenBasis& basis)
{
    if(m_cond_template.type() == type_border::side)
        gsInfo << "Strong BC" << std::endl;
    else if(m_cond_target.type() == type_border::geom)
    {
        refineBoundary(m_cond_template.geometry(), basis);
        if( !param.continuous_fitting)
        {
            std::vector< gsMatrix<T> > tmp_pts, tmp_param;
            sample_boundary_points<T>(m_cond_target.geometry(),
                                      m_cond_template.geometry(),
                                      tmp_pts, tmp_param);
            m_cond_target.points() = tmp_pts[0];
            m_cond_template.points() = tmp_param[0];
        }
    }
    else if(m_cond_target.type() == type_border::trim_surf)
        initializeTrim(param, basis);
}

template<class T> template<class GenBasis>
void gsFittingBorderCoupling<T>::
associatePts(gsFittingEnergy<GenBasis, T>& energy,
             GenBasis& basis, index_t ind)
{
    energy.addPointsLS(m_cond_target.points(),
                     m_cond_template.points(), ind);
}


template<class T> template<short_t d, class GenBasis>
void gsFittingBorderCoupling<T>::
associateGeom(gsFittingEnergy<GenBasis, T>& energy,
              GenBasis& basis, index_t ind)
{
    addBCQuadrature<d, GenBasis, gsGeometry<T>, 2, T>
        (basis, m_cond_template.geometry(), m_cond_target.geometry(),
         energy.dim_im(), energy.quadrEner());
}

template<class T> template<class GenBasis>
void gsFittingBorderCoupling<T>::
associateTrim(gsFittingEnergy<GenBasis, T>& energy,
              GenBasis& basis)
{
    /// GISMO_ASSERT...    ensures that the initialization has been called


}

template<class T> template<short_t d, class GenBasis>
void gsFittingBorderCoupling<T>::
associate(gsFittingBase<d, GenBasis, T>& fitting, index_t ind)
{
    if(m_cond_template.type() == type_border::side)
        associateStrong<d, GenBasis>(fitting);
    else if(m_cond_target.type() == type_border::pts)
        associatePts<GenBasis>(fitting.energy(),
                               fitting.basis(), ind);
    else if(m_cond_target.type() == type_border::geom)
    {
        /// In case we have sampled the points
        if(m_cond_template.points().size() > 0)
            associatePts<GenBasis>(fitting.energy(),
                                   fitting.basis(), ind);
        else
            associateGeom<d, GenBasis>(fitting.energy(),
                                       fitting.basis(), ind);
    }
    else if(m_cond_target.type() == type_border::trim_surf)
        associateTrim<GenBasis>(fitting.energy(),
                                fitting.basis());
}


template<class T> template<short_t d, class GenBasis>
void gsFittingBorderCoupling<T>::
associateStrong(gsFittingBase<d, GenBasis, T>& fitting)
{

}

template<class T> template<class GenBasis> void
gsFittingBorderCoupling<T>::
initializeTrim(gsFittingParam<T>& param, GenBasis& basis)
{
    gsFittingParam<T> param_face(param);
    param_face.interiorKnots =
        param.interiorKnots_boundary;

    gsPlanarDomain<T>& domain_temp
        ( m_cond_template.trimmed().domain() );
    gsPlanarDomain<T>& domain_targ
        ( m_cond_target.trimmed().domain() );
    gsTemplateTargetFitting<T> mapping_domain(domain_temp,
                                              domain_targ,
                                              param_face);
    gsFunctionSet<T>* _res = mapping_domain.template
        computeMappingTrimming<2>();
    gsGeometry<T>* res = static_cast<gsGeometry<T>*>(_res);

    /// The parametrization is a composition of two mappings
    /// (res and the mapping from the planar domain to the 3D face)
    /// Hence, we only consider the sampling to simplify

    std::vector< gsMatrix<T> > points_loc(1);
    std::vector< gsMatrix<T> > param_loc(1);

    gsGeometry<T>& geom_targ
        (*m_cond_target.trimmed().getTP());
    gsGeometry<T>& geom_temp
        (*m_cond_template.trimmed().getTP());

    sample_boundary_points<T>(res->basis(), geom_targ, geom_temp,
                           points_loc, param_loc, res);
    m_cond_target.points() = points_loc[0];
    m_cond_template.points() = param_loc[0];
}


} /// namespace gismo
