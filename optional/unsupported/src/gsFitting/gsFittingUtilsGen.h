/** @file gsFittingUtilsGen.h

    @brief Small routines used for keeping genericity of the code.
    All these routines are implemented for single and multi patch.
    Some of these functions have been defined here in order to keep
    the stable part unchanged.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingQuadrature.h>
#include <gsFitting/gsFittingIntegrandNL.h>

namespace gismo
{

template<short_t d, class GenBasis, class T>
class gsFittingBase;

/// Returns the bounding box of a geometry.
/// In gsMultiPatch this function is defined,
/// while this function is not implemented in gsGeometry.
/// Called when the boundary of the template domain is given
/// in order to define the basis.
template<class T>
void setBoundingBox(gsGeometry<T>& geom, gsMatrix<T>& bb);
template<class T>
void setBoundingBox(gsMultiPatch<T>& geom, gsMatrix<T>& bb);

/// Add the Integrand associated with the L2 distance
/// to the template to the quadrature quadr
template<short_t d, class T> void
add_L2DistMapping(gsFittingQuadrature<d, gsBasis<T>, 2, T> &quad,
                  gsFunctionSet<T>* mapping);
template<short_t d, class T> void
add_L2DistMapping(gsFittingQuadrature<d, gsMultiBasis<T>, 2, T> &quad,
                  gsFunctionSet<T>* mapping);

/// Returns the first basis of b (trivial if b is a single patch basis)
template<class T> inline gsBasis<T>& firstBasis(gsMultiBasis<T>& b)
{ return b[0];  }
template<class T> inline gsBasis<T>& firstBasis(gsBasis<T>& b)
{ return b;  }

/// Returns the basis of b in some patch
/// (the patch is not considered if b is a single patch basis)
template<class T> inline gsBasis<T>& getBasisGen(gsMultiBasis<T>& b,
                                              int patch)
{ return b[patch];  }
template<class T> inline gsBasis<T>& getBasisGen(gsBasis<T>& b,
                                              int patch)
{    return b;  }

/// Returns the size of the basis
template<class T> inline int nBasesGen(gsMultiBasis<T>& b)
{ return b.nBases();  }
template<class T> inline int nBasesGen(gsBasis<T>& b)
{    return 1;  }

template<class T> inline int nPatchesGen(gsMultiPatch<T>& g)
{ return g.nPatches();  }
template<class T> inline int nPatchesGen(gsGeometry<T>& g)
{    return 1;  }

/// Evaluates a geometry at some points.
/// The basis is here to distinguish the case single
/// patch from the case multipatch.
/// The result is set in res.
template<class T> void
generic_eval_into(gsBasis<T>& basis, gsFunctionSet<T>& geometry,
                  std::vector< gsMatrix<T> > &param,
                  std::vector< gsMatrix<T> > &res);
template<class T> void
generic_eval_into(gsMultiBasis<T>& basis, gsFunctionSet<T>& patches,
                  std::vector< gsMatrix<T> > &param,
                  std::vector< gsMatrix<T> > &res);

/// Supposing that x contains the coefficients of a geometry and
/// supposing that these coefficients are defined on the same basis
/// than the geometry patches, we some these two geometries.
/// The basis is here to distinguish the case single
/// patch from the case multipatch.
template<class T>
void addMappingToCoeff(gsBasis<T>& basis,
                       gsFunctionSet<T> &patches, gsMatrix<T>& x);
template<class T>
void addMappingToCoeff(gsMultiBasis<T>& basis,
                       gsFunctionSet<T> &patches, gsMatrix<T>& x);

/// Supposing that geom_left and geom_right are geometries defined
/// on the same basis, we some these two geometries.
/// The basis is here to distinguish the case single
/// patch from the case multipatch.
template<class T>
void sumGenericGeometries(gsBasis<T>& basis,
                          gsFunctionSet<T> &geom_left,
                          gsFunctionSet<T> &geom_right,
                          T coeff = 1.);
template<class T>
void sumGenericGeometries(gsMultiBasis<T>& basis,
                          gsFunctionSet<T> &patches_left,
                          gsFunctionSet<T> &patches_right,
                          T coeff = 1.);

/// Casts the object geom into a gsGeometry
/// The basis is here to distinguish the case single
/// patch from the case multipatch.
template<class T> gsGeometry<T>&
getGeometry(gsBasis<T>& basis, gsFunctionSet<T>& geom);
template<class T> gsMultiPatch<T>&
getGeometry(gsMultiBasis<T>& basis, gsFunctionSet<T>& geom);

/// Returns the n^th geometry (n is not taken into account here)
/// The basis is here to distinguish the case single
/// patch from the case multipatch.
template<class T> gsGeometry<T>&
getNthGeometry(gsBasis<T>& basis,
               gsFunctionSet<T>& geom, int n);
template<class T> gsGeometry<T>&
getNthGeometry(gsMultiBasis<T>& basis,
               gsFunctionSet<T>& geom, int n);

/// Computes the value of some geometry at the points pos
template<class T>
void fitting_eval_into(gsGeometry<T>& geom, gsMatrix<T>& pos,
                       gsMatrix<T>& res, int ind_patch = 0);
template<class T>
void fitting_eval_into(gsMultiPatch<T>& patches, gsMatrix<T>& pos,
                       gsMatrix<T>& res, int ind_patch = 0);

/// Computes the derivatives of some geometry at the points pos
template<class T>
void fitting_deriv_into(gsGeometry<T>& geom, gsMatrix<T>& pos,
                       gsMatrix<T>& res, int ind_patch = 0);
template<class T>
void fitting_deriv_into(gsMultiPatch<T>& patches, gsMatrix<T>& pos,
                        gsMatrix<T>& res, int ind_patch = 0);


/// Returns the dimension of the domain of the basis
template<class T> int dimGenBasis(gsMultiBasis<T>& basis)
{  return basis[0].dim();}

template<class T> int dimGenBasis(gsBasis<T>& basis)
{  return basis.dim();  }

/// Copies a given basis
template<class T>
gsBasis<T> copyBasis(gsGeometry<T>& geom);

template<class T>
gsMultiBasis<T> copyBasis(gsMultiPatch<T>& patches);

/// Copies the coefficient of a geometry
template<class T> void
geometryToCoeffGen(gsBasis<T>& basis, gsFunctionSet<T>& func,
                   gsMatrix<T>& res,
                   gsLocalGlobal<gsBasis<T>, T>& local_global);

template<class T> void
geometryToCoeffGen(gsMultiBasis<T>& basis, gsFunctionSet<T>& func,
                   gsMatrix<T>& res,
                   gsLocalGlobal<gsMultiBasis<T>, T>& local_global);

/// Return false if the template and target are now empty
/// (all the bc are now imposed strongly)

template<class T> bool
split_multiPatch(gsFunctionSet<T>* templ,
                 gsFunctionSet<T>* target,
                 std::vector<int>& indices,
                 gsMultiPatch<T>& bc_target,
                 gsMultiPatch<T>& bc_template);
template<class T> bool
split_multiPatchDyn(gsFunctionSet<T>* templ,
                    gsFunctionSet<T>* target,
                    std::vector<int>& indices,
                    gsMultiPatch<T>& bc_target,
                    gsMultiPatch<T>& bc_template);

} /// namespace gismo
