/** @file gsFittingUtilsSampling.h

    @brief This file contains sampling routines.
    The sampling is each time performed by using gsPointsSampling, i.e. we compute the points of the gaussian quadrature associated with some given basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingUtilsInit.h>

namespace gismo
{

/// Sample points from geom and templ geometries.
/// The domain in which basis and these two geometries are
/// defined must be the same
template<class T>
void sample_boundary_points(gsBasis<T>& basis,
                            gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            gsGeometry<T>* temp_reparam = NULL,
                            int* local_size_sample = NULL);
template<class T>
void sample_boundary_points(gsBasis<T>& basis,
                            gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int local_size_sample,
                            gsGeometry<T>* temp_reparam = NULL);
template<class T>
void sample_boundary_points(gsMultiBasis<T>& basis,
                            gsMultiPatch<T>& geom,
                            gsMultiPatch<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            gsMultiPatch<T>* temp_reparam = NULL,
                            int* local_size_sample = NULL);
template<class T>
void sample_boundary_points(gsMultiBasis<T>& basis,
                            gsMultiPatch<T>& geom,
                            gsMultiPatch<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int local_size_sample,
                            gsMultiPatch<T>* temp_reparam = NULL);

/// Sample points from geom and templ geometries.
/// The domain where these two geometries are
/// defined must be the same
template<class T>
void sample_boundary_points(gsMultiPatch<T>& geom,
                            gsMultiPatch<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            gsMultiPatch<T>* temp_reparam,
                            int local_size_sample);
template<class T>
void sample_boundary_points(gsMultiPatch<T>& geom,
                            gsMultiPatch<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            gsMultiPatch<T>* temp_reparam = NULL,
                            int* local_size_sample = NULL);
template<class T>
void sample_boundary_points(gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int local_size_sample,
                            gsGeometry<T>* temp_reparam = NULL);
template<class T>
void sample_boundary_points(gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            gsGeometry<T>* temp_reparam = NULL,
                            int* local_size_sample = NULL);

template<class T>
void sample_boundary_points(gsBasis<T>& basis,
                            gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int pos,
                            int* local_size_sample = NULL);
template<class T>
void sample_boundary_points(gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int pos, int* local_size_sample = NULL);


/// Refine the boundary until it has "as many degrees of freedom
/// as the domain along this boundary"
template<class T>
void refineBoundary(gsMultiPatch<T>& bound, gsBasis<T>& basisDom);
template<class T>
void refineBoundary(gsGeometry<T>& bound, gsBasis<T>& basisDom);

template<class T> void refineBoundary(gsMultiPatch<T>& bound,
                                      gsMultiBasis<T>& basisDom);
template<class T> void refineBoundary(gsGeometry<T>& bound,
                                      gsMultiBasis<T>& basisDom);



} /// namespace gismo
