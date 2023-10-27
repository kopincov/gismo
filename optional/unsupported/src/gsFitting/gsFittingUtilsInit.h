/** @file gsFittingUtilsInit.h

    @brief Contains some routines used for initializing datas
    necessary to apply gsFitting

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

#include <gsFitting/gsFittingBase.h>

namespace gismo
{

/// Computes the bounding box of the points pts
template<class T> void
getBoundingBoxPts(gsMatrix<T>& pts, gsMatrix<T>& res);

/// Join the points and parameters into one matrix
template <class T>
void junctionPoints(std::vector< gsMatrix<T> > &pts,
                    std::vector< gsMatrix<T> > &param,
                    std::vector< gsMatrix<T> > &pts_res,
                    std::vector< gsMatrix<T> > &param_res,
                    gsMatrix<T>* bb = NULL);


/// Creates a copy of the basis of the multipatch mp.
/// If is_thb, then the resulting basis is a THB-Splines
template<short_t d, class T> void
basesFromMultiPatch(gsMultiPatch<T>& mp, gsMultiBasis<T>& res,
                    bool keep_basis, int degree = 2,
                    int num_internal_knots = 0,
                    bool is_thb = 1);
template<short_t d, class T> void
copyBases(gsMultiBasis<T>& mp, gsMultiBasis<T>& res,
          bool keep_basis,
          int degree = 2, int num_internal_knots = 0,
          bool is_thb = 1);
template<class T> void
basesFromMultiPatchDyn(gsMultiPatch<T>& mp, gsMultiBasis<T>& res,
                       bool keep_basis, int degree = 2,
                       int num_internal_knots = 0,
                       bool is_thb = 1);
template<class T> void
copyBasesDyn(gsMultiBasis<T>& mp, gsMultiBasis<T>& res,
             bool keep_basis, int degree = 2,
             int num_internal_knots = 0,
             bool is_thb = 1);

/**
 * Joins all the parameter and image points into one gsMatrix.
 * Removes moreover all the points that are not includedin
 * the bounding box bb (if given)
 */
template<class T> bool
isIncluded(const typename gsMatrix<T>::Column& vect,
           gsMatrix<T>& bb_param);

template<class T> gsBasis<T>
getBasisDyn(short_t d, std::vector< gsKnotVector<T> >& KV,
            bool use_refinement, gsBasis<T>& res);

} /// namespace gismo
