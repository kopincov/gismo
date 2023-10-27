/** @file gsFittingUtilsGen.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingUtilsGen.h>

namespace gismo
{
template<class T>
gsBasis<T> copyBasis(gsGeometry<T>& geom)
{   return gsBasis<T>(geom.basis());  }

template<class T>
gsMultiBasis<T> copyBasis(gsMultiPatch<T>& patches)
{
    std::vector<gsBasis<T>*> liste_basis
        = patches.basesCopy();
    return gsMultiBasis<T>(liste_basis, patches);

}

template<class T> void
generic_eval_into(gsBasis<T>& basis, gsFunctionSet<T>& geometry,
                  std::vector< gsMatrix<T> > &param,
                  std::vector< gsMatrix<T> > &res)
{
    gsGeometry<T>& _geometry =
        static_cast<gsGeometry<T>&>(geometry);
    _geometry.eval_into(param[0], res[0]);
}

template<class T> void
generic_eval_into(gsMultiBasis<T>& basis, gsFunctionSet<T>& patches,
                  std::vector< gsMatrix<T> > &param,
                  std::vector< gsMatrix<T> > &res)
{
    gsMultiPatch<T>& _patches =
        static_cast<gsMultiPatch<T>&>(patches);
    unsigned s = _patches.nPatches();
    for(unsigned i = 0;i < s;i++)
    {
        _patches[i].eval_into(param[i], res[i]);
    }
}


template<class T>
void sumGenericGeometries(gsBasis<T>& basis,
                          gsFunctionSet<T> &geom_left,
                          gsFunctionSet<T> &geom_right,
                          T coeff)
{
    gsGeometry<T>& _geom_left =
        static_cast<gsGeometry<T>&>(geom_left);
    gsGeometry<T>& _geom_right =
        static_cast<gsGeometry<T>&>(geom_right);
    gsMatrix<T>& coeffs_left = _geom_left.coefs();
    gsMatrix<T>& coeffs_right = _geom_right.coefs();
    coeffs_left += coeffs_right * coeff;
}

template<class T>
void sumGenericGeometries(gsMultiBasis<T>& basis,
                          gsFunctionSet<T> &patches_left,
                          gsFunctionSet<T> &patches_right,
                          T coeff)
{
    gsMultiPatch<T>& _patches_left =
        static_cast<gsMultiPatch<T>&>(patches_left);
    gsMultiPatch<T>& _patches_right =
        static_cast<gsMultiPatch<T>&>(patches_right);
    unsigned s = _patches_left.nPatches();

    for(unsigned i = 0;i < s;i++)
    {
        gsMatrix<T>& coeffs_left = _patches_left.patch(i).coefs();
        gsMatrix<T>& coeffs_right = _patches_right.patch(i).coefs();
            coeffs_left += coeff * coeffs_right;
    }
}


template<class T> gsGeometry<T>&
getGeometry(gsBasis<T>& basis, gsFunctionSet<T>& geom)
{
    return static_cast<gsGeometry<T>&>(geom);
}

template<class T> gsMultiPatch<T>&
getGeometry(gsMultiBasis<T>& basis, gsFunctionSet<T>& geom)
{
    return static_cast<gsMultiPatch<T>&>(geom);
}

/// Returns the n^th geometry (n is not taken into account here)
template<class T> gsGeometry<T>&
getNthGeometry(gsBasis<T>& basis, gsFunctionSet<T>& geom, int n)
{
    return getGeometry(basis, geom);
}

/// Returns the n^th geometry
template<class T> gsGeometry<T>&
getNthGeometry(gsMultiBasis<T>& basis,
               gsFunctionSet<T>& geom, int n)
{
    return getGeometry<T>(basis, geom).patch(n);
}

template<class T>
void fitting_eval_into(gsGeometry<T>& geom, gsMatrix<T>& pos, gsMatrix<T>& res, int ind_patch)
{
    geom.eval_into(pos, res);
}

template<class T>
void fitting_eval_into(gsMultiPatch<T>& patches, gsMatrix<T>& pos, gsMatrix<T>& res, int ind_patch)
{
    patches[ind_patch].eval_into(pos, res);
}

template<class T>
void fitting_deriv_into(gsGeometry<T>& geom,
                       gsMatrix<T>& pos,
                       gsMatrix<T>& res,
                       int ind_patch)
{
    geom.deriv_into(pos, res);
}

template<class T>
void fitting_deriv_into(gsMultiPatch<T>& patches,
                       gsMatrix<T>& pos, gsMatrix<T>& res,
                       int ind_patch)
{
    patches[ind_patch].deriv_into(pos, res);
}


template<short_t d, class T> void
add_L2DistMapping(gsFittingQuadrature<d, gsBasis<T>, 2, T> &quad,
                  gsFunctionSet<T>* mapping)
{
    GISMO_ASSERT(mapping != NULL,
                 "ERROR: the mapping is not defined");
    gsGeometry<T>* _map =
        static_cast<gsGeometry<T>*>(mapping);
    quad.addIntegrand(new gsFittingIntegrandL2Dist
                      <d, gsGeometry<T>, T>(d, _map) );
}

template<short_t d, class T> void
add_L2DistMapping(gsFittingQuadrature<d, gsMultiBasis<T>, 2, T> &quad,
                  gsFunctionSet<T>* mapping)
{
    GISMO_ASSERT(mapping != NULL,
                 "ERROR: the mapping is not defined");
    gsMultiPatch<T>* _map = static_cast<gsMultiPatch<T>*>
        (mapping);
    quad.addIntegrand(new gsFittingIntegrandL2Dist
                      <d, gsMultiPatch<T>, T>(d, _map) );
}


template<class T> void
setBoundingBox(gsGeometry<T>& geom, gsMatrix<T>& bbox)
{
    boundingBox(geom, bbox);
}

template<class T> void
setBoundingBox(gsMultiPatch<T>& geom, gsMatrix<T>& bbox)
{
    geom.boundingBox(bbox);
}

template<class T> void
geometryToCoeffGen(gsBasis<T>& basis, gsFunctionSet<T>& func,
                   gsMatrix<T>& res,
                   gsLocalGlobal<gsBasis<T>, T>& local_global)
{
    gsGeometry<T>& geom = getGeometry(basis, func);
    res = geom.coefs().replicate(1, 1);
}

template<class T> void
geometryToCoeffGen(gsMultiBasis<T>& basis,
                   gsFunctionSet<T>& func,
                   gsMatrix<T>& res,
                   gsLocalGlobal<gsMultiBasis<T>, T>& local_global)
{
    gsMultiPatch<T>& geom = getGeometry(basis, func);
    unsigned nPatches = basis.nBases();
    unsigned size_tot = local_global.sizeDof();

    res.resize(size_tot, basis.domainDim());
    res.setZero();
    for(unsigned i = 0;i < nPatches;i++)
    {
        const unsigned numBasisFun = basis.size(i);
        for(unsigned j = 0;j < numBasisFun;j++)
        {
            const int globalI = local_global.localToGlobal
                (j, i, true);
            res.row(globalI) = geom.patch(i).coef(j);
        }
    }
}

template<class T> bool
split_multiPatch(gsBasis<T>& basis, gsFunctionSet<T>* templ,
                 gsFunctionSet<T>* target,
                 std::vector<int>& indices,
                 gsMultiPatch<T>& bc_target,
                 gsMultiPatch<T>& bc_template)
{
    GISMO_ENSURE(indices.size() == 0, "There should be no boundary conditions");
    return true;
}

template<class T> bool
split_multiPatch(gsFunctionSet<T>* templ,
                 gsFunctionSet<T>* target,
                 std::vector<int>& indices,
                 gsMultiPatch<T>& bc_target,
                 gsMultiPatch<T>& bc_template)

{
    gsMultiPatch<T>* _templ
        = static_cast<gsMultiPatch<T>*>(templ);
    gsMultiPatch<T>* _target
        = static_cast<gsMultiPatch<T>*>(target);

    typedef typename gsMultiPatch<T>::PatchContainer Container;
    Container cont_templ = _templ->patches();
    Container cont_target = _target->patches();
    unsigned s = _templ->nPatches();
    const std::vector<int> ind_bc = indices_bc(s, indices);
    int nb_bc = 0;
    int del = 1 + 2*s;
    for(unsigned i = 0;i < s;i++)
    {
        std::vector<int>::const_iterator pos
            = std::find (ind_bc.begin(),
                         ind_bc.end(), ind_bc[i]);
        if(pos != indices.end())
        {
            /// This patch is in the set of indices.
            /// We remove it from template and target
            /// geometries and we add it to new bc targets and
            /// geometries
            bc_target.addPatch(_target->patch(i));
            bc_template.addPatch(_templ->patch(i));

            cont_templ.erase(cont_templ.begin() + i);
            cont_target.erase(cont_target.begin() + i);

            indices[i * del] = nb_bc;
            nb_bc++;
            i--;
        }
    }
    if(_templ->nPatches() == 0)
    {
        _templ->clear();
        _target->clear();
        return false;
    }
    return true;
}


template<class T> bool
split_multiPatchDyn(gsFunctionSet<T>* templ,
                    gsFunctionSet<T>* target,
                    std::vector<int>& indices,
                    gsMultiPatch<T>& bc_target,
                    gsMultiPatch<T>& bc_template)
{
    if(dynamic_cast<gsMultiPatch<T>* >(templ) != NULL)
        return split_multiPatch(templ, target, indices,
                                bc_template, bc_template);
    return true;
}


}  /// namespace gismo
