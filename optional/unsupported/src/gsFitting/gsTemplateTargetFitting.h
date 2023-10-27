/** @file gsTemplateTargetFitting.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

namespace gismo
{


template<class T>
class gsTemplateTargetFitting
{
protected:
  /// The template geometry
  gsFittingGeometry<T> m_template;

  /// The target geometry
  gsFittingGeometry<T> m_target;

  /// The coupling between the boundaries of the template
  /// and the target geometries
  gsFittingGeometryCoupling<T> m_coupling;

  /// The parameters of the computation
  gsFittingParam<T>* m_param_fitting;

public:
  gsTemplateTargetFitting(gsFittingParam<T>& param)
  : m_template(), m_target(), m_coupling(m_template, m_target),
    m_param_fitting(&param)
  {  init();  }

  gsTemplateTargetFitting(std::vector<gsMatrix<T> >& pts,
                          std::vector<gsMatrix<T> >& param,
                          bool unique_patch,
                          gsFittingParam<T>& param_fitting)
  : m_template(param, unique_patch), m_target(pts, unique_patch),
    m_coupling(m_template, m_target), m_param_fitting(&param_fitting)
  { init(); }

  gsTemplateTargetFitting(gsMultiPatch<T>& geom,
                          gsMultiPatch<T>& templ,
                          bool unique_patch,
                          gsFittingParam<T>& param)
  : m_template(templ, unique_patch), m_target(geom, unique_patch),
    m_coupling(m_template, m_target), m_param_fitting(&param)
  {  init(); }

  gsTemplateTargetFitting(const gsTemplateTargetFitting<T>&
                          source)
  : m_template(source.m_template), m_target(source.m_target),
    m_coupling(source.m_coupling),
    m_param_fitting(source.m_param_fitting)
  {  init(); }


  ////////// The geometries have only one patch ////////////

  gsTemplateTargetFitting(gsMatrix<T>& pts, gsMatrix<T>& param,
                          gsFittingParam<T>& param_fitting)
  : m_template(param_fitting), m_target(pts),
    m_coupling(m_template, m_target),
    m_param_fitting(&param_fitting)
  {  init(); }

  gsTemplateTargetFitting(gsGeometry<T>& geom,
                          gsGeometry<T>& templ,
                          gsFittingParam<T>& param_fitting)
  : m_template(templ), m_target(geom),
    m_coupling(m_template, m_target),
    m_param_fitting(&param_fitting)
  {  init(); }

  gsTemplateTargetFitting(gsSolid<T>& geom,
                          gsSolid<T>& templ,
                          gsFittingParam<T>& param_fitting)
  : m_template(templ), m_target(geom),
    m_coupling(m_template, m_target),
    m_param_fitting(&param_fitting)
  {  init(); }

  gsTemplateTargetFitting(gsPlanarDomain<T>& geom,
                          gsPlanarDomain<T>& templ,
                          gsFittingParam<T>& param_fitting)
  : m_template(templ), m_target(geom),
    m_coupling(m_template, m_target),
    m_param_fitting(&param_fitting)
  {  init(); }

  gsTemplateTargetFitting<T>&
  operator=(const gsTemplateTargetFitting<T>& source)
  {
    m_template = gsFittingGeometry<T>(source.m_template);
    m_target = gsFittingGeometry<T>(source.m_target);
    m_param_fitting = source.m_param_fitting;
    m_coupling = gsFittingGeometryCoupling<T>(m_template, m_target);
    return *this;
  }

  index_t nPatches() {  return m_target.nPatches();  }

  index_t dim_im()  { return m_target.dim_im();  }

  index_t dim_dom()  { return m_template.dim_im();  }

  gsFittingParam<T>& fitting_param()
  {   return *m_param_fitting;  }

  template<unsigned d> gsFunctionSet<T>*
  computeMappingTrimming(gsMultiBasis<T>* basis = NULL);

  gsFunctionSet<T>*
  computeMappingTrimmingDyn(gsMultiBasis<T>* basis = NULL);

  template<class GenBasis> gsFunctionSet<T>* computeMapping(GenBasis& basis);

private:
  void init()
  {
    GISMO_ASSERT(m_template.nPatches() == m_target.nPatches(),
                 "Error in the data");
  }

  //////////////////////////////////////////////////////
};


} /// namespace gismo
