/** @file gsFittingUtilsSolids.h

    @brief Contains routines permitting to construct a mapping from a
    given template solid to a given target solid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingUtilsInit.h>
#include <gsFitting/gsFittingBase.h>
#include <gsFitting/gsFittingUtilsIO.h>
#include <gsFitting/gsFittingConstr.h>

#include <gsUtils/gsMesh/gsBoundingBox.h>

#include <gsModeling/gsSolidElement.h>
#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsSolid.h>

#include <gsIO/gsReadFile.h>

#include <gsHSplines/gsTHBSplineBasis.h>


namespace gismo
{

/// Returns a spline basis having the same domain as surf
/// with nbInteriorKnots interior knots and a certain degree.
template<short_t d, class T> gsTensorBSplineBasis<d>
getBasisFromTemplate(gsGeometry<T>& surf, int nbInteriorKnots,
                     int degree);

/// We fit a gsSolid using gsFittingLinSyst keeping at each step
/// the parametrisation of the boundary (edges and faces) fixed
template<class T> gsFunctionSet<T>*
fitSolidFixedParam(gsSolid<T>* temp, gsSolid<T>* geom,
                   gsFittingParam<T>& fitting_param,
                   int nPointsSample = 100);

/// We fit a Face using gsFittingLinSyst keeping the parametrization
///  of the boundary (edges) fixed
template<class T> gsFunctionSet<T>*
fitFaceFixedParam(gsTrimSurface<T>* temp, gsTrimSurface<T>* geom,
                  gsBasis<T>& basis, gsFittingParam<T>& param,
                  int nPointsSample = 100);

/// Constructs a basis having a given bounding box.
template<short_t d, class T> gsTensorBSplineBasis<d>
getBasisFromBoundingBox(gsMatrix<T> bbox,
                        int nbInteriorKnots, int degree);
template<short_t d, class T> gsTensorBSplineBasis<d>
getBasisFromBoundingBox(gsBoundingBox<T>& bbox,
                        int nbInteriorKnots, int degree);

/*
  Routine used for the multipatch solid parametrization:
  applies the middle edge splitting to each face to set
  the target points (boundary points of the faces, i.e.,
  the edges)
 */
template<short_t d, class T> void
set_strong_bc_face(gsMultiPatch<T>& boundary,
                   gsMultiBasis<T>& basis_face,
                   gsDofMapper mapper);

} /// namespace gismo
