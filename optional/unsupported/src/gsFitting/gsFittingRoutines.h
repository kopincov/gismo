/** @file gsFittingRoutines.h

    @brief Contains the three routines fittingTargetDyn,
    fittingPointsDyn, and fittingSolid, that should be called to use
    gsFitting. For the two first routines, if the basis and the topology
    are not given, we suppose that the fitting is single patch
    and the basis constructed on the bounding box of the
    image of the template geometry. Otherwise, the basis
    is either kept unchanged or modified (while keeping the domain
    unchanged) depending on the parameters entered by the user.
    Finally, depending on the dimension, type of boundary,etc, the
    appropriate fitting routine is called.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingIterative.h>
#include <gsFitting/gsFittingHierar.h>
#include <gsFitting/gsFittingIntegrandNL.h>
#include <gsFitting/gsFittingUtilsSampling.h>
#include <gsFitting/gsFittingUtilsIO.h>
#include <gsFitting/gsFittingUtilsSolids.h>

#include <gsFitting/gsTemplateTargetFitting.h>

namespace gismo
{

/**
   Fits a mapping such that the image of the target
    geometry is the template geometry.
 */
template <class GenGeom, class T>
int fittingTargetDyn(gsFittingParam<T> &param,
                     GenGeom& geom, GenGeom& templ,
                     gsMultiPatch<T>* topology,
                     gsMultiBasis<T>* basis);


/** Fits a mapping such that the image of the points param are
    the points pts.
 */
template <class T>
int fittingPointsDyn(gsFittingParam<T> &fitting_param,
                     std::vector<gsMatrix<T> > &pts,
                     std::vector<gsMatrix<T> > &param,
                     gsMultiPatch<T>* topology,
                     gsMultiBasis<T>* basis);

/** Fits a mapping such that the image of the solid temp
    is the solid geom. TODO: use multipatch
 */
template<class T> int
fittingSolid(gsSolid<T>* temp, gsSolid<T>* geom,
             gsFittingParam<T>& fitting_param);

template<class T>
int fittingTrimming(gsMultiPatch<T>* topology,
                    gsMultiBasis<T>* basis,
                    gsTemplateTargetFitting<T>& data);

/***** ONLY THE FUNCTIONS ABOVE SHOULD BE USED ****/

/*
  Fits a mapping such that the image of the geometry
  templ is the geometry geom.
  GenGeom: the type of geometry of the boundary
  (gsGeometry or gsMultiPatch)
*/
template <short_t d, class GenGeom, class T>
int fittingTargetSinglePatch(gsFittingParam<T> &param,
                             GenGeom& geom, GenGeom& templ);
template <short_t d, class GenGeom, class T>
int fittingTargetMultiPatch(gsFittingParam<T> &param,
                            gsMultiBasis<T> &basis,
                            GenGeom& geom, GenGeom& templ);

/*
  Fits a mapping such that the image of the points
  param are the points pts.
  GenGeom: the type of geometry of the boundary
  (gsGeometry or gsMultiPatch)
*/
template <short_t d, class T>
int fittingPointsMultiPatch(gsFittingParam<T> &fitting_param,
                            gsMultiBasis<T>& basis,
                            std::vector<gsMatrix<T> > &pts,
                            std::vector<gsMatrix<T> > &param);
template <short_t d, class T>
int fittingPointsSinglePatch(gsFittingParam<T> &fitting_param,
                             std::vector<gsMatrix<T> > &pts,
                             std::vector<gsMatrix<T> > &param);

} /// namespace gismo
