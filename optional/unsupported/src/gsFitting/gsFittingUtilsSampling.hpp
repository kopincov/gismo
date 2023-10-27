/** @file gsFittingUtilsSampling.hpp


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingUtilsSampling.h>
#include <gsFitting/gsFittingUtilsInit.hpp>

namespace gismo
{

template<class T>
void sample_boundary_points(gsBasis<T>& basis,
                            gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector< gsMatrix<T> >& res_geom,
                            std::vector< gsMatrix<T> >& res_template,
                            gsGeometry<T>* temp_reparam,
                            int* local_size_sample)
{
    std::vector<gsMatrix<T> > param;
    gsPointsSampling<gsBasis<T>, T>
        sampler(basis, local_size_sample);
    unsigned total_size = sampler.total_size();
    unsigned d_geom = geom.targetDim();
    unsigned d_param = geom.domainDim();

    res_geom.clear();
    res_template.clear();

    param.push_back(gsMatrix<T>(d_param, total_size));
    sampler.points_sampling(param);

    if(temp_reparam != NULL)
    {
        std::vector<gsMatrix<T> > reparam;
        gsMatrix<T> bb = geom.basis().support();

        reparam.push_back(gsMatrix<T>(d_param, total_size));
        temp_reparam->eval_into(param[0], reparam[0]);

        /// We remove all the parameters with image by temp_reparam
        /// not included in the parameter domain of the geometry
        for(unsigned i = 0;i < total_size;i++)
        {
            if(! isIncluded<T>(reparam[0].col(i), bb))
            {
                reparam[0].removeCol(i);
                param[0].removeCol(i);
                total_size--;
                i--;
            }
        }
        res_geom.push_back(gsMatrix<T>
                           (d_geom, reparam[0].cols()));
        geom.eval_into(reparam[0], res_geom[0]);
    }
    else
    {
        res_geom.push_back(gsMatrix<T>(d_geom, total_size));
        geom.eval_into(param[0], res_geom[0]);
    }

    res_template.push_back(gsMatrix<T>(d_geom, param[0].cols()));
    templ.eval_into(param[0], res_template[0]);
}

template<class T>
void sample_boundary_points(gsBasis<T>& basis,
                            gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int pos,
                            int* local_size_sample)
{
    std::vector<gsMatrix<T> > param;
    gsPointsSampling<gsBasis<T>, T>
        sampler(basis, local_size_sample);
    unsigned total_size = sampler.total_size();
    unsigned d_geom = geom.targetDim();
    unsigned d_param = geom.domainDim();

    param.push_back(gsMatrix<T>(d_param, total_size));
    res_geom[pos] = gsMatrix<T>(d_geom, total_size);
    res_template[pos] = gsMatrix<T>(d_geom, total_size);

    sampler.points_sampling(param);
    geom.eval_into(param[pos], res_geom[pos]);
    templ.eval_into(param[pos], res_template[pos]);
}


template<class T>
void sample_boundary_points(gsBasis<T>& basis,
                            gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int local_size_sample,
                            gsGeometry<T>* temp_reparam)
{
    int* p_local_size = NULL;
    if(local_size_sample > 0)
    {
        int tab_local_size[4];
        unsigned dim = basis.domainDim();
        for(unsigned i = 0;i < dim;i++)
            tab_local_size[i] = local_size_sample;
        p_local_size = tab_local_size;
    }
    sample_boundary_points(basis, geom, templ, res_geom,
                           res_template, temp_reparam,
                           p_local_size);
}

/// The domain is supposed to be a single patch.
/// geom is a multipatch geometry
/// having image in this domain
/// res_geom and res_template are both vectors of size 1
template<class T>
void sample_boundary_points(gsMultiBasis<T>& basis,
                            gsMultiPatch<T>& geom,
                            gsMultiPatch<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            gsMultiPatch<T>* temp_reparam,
                            int* local_size_sample)
{
    GISMO_ENSURE(temp_reparam == NULL, "NOT IMPLEMENTED YET (TODO)");
    gsMultiPatch<T>& _geom =
        static_cast<gsMultiPatch<T>&>(geom);
    gsMultiPatch<T>& _templ =
        static_cast<gsMultiPatch<T>&>(templ);
    gsMultiPatch<T>* _reparam = NULL;
    if(temp_reparam != NULL)
        _reparam = static_cast<gsMultiPatch<T>*>
            (temp_reparam);


    std::vector< gsMatrix<T> > tmp_geom;
    std::vector< gsMatrix<T> > tmp_template, tmp_reparam;

    std::vector<gsMatrix<T> > param;
    gsPointsSampling<gsMultiBasis<T>, T>
        sampler(basis, local_size_sample);

    std::vector<index_t> sizes = sampler.sizes();
    unsigned s = sizes.size();
    unsigned total_size;
    unsigned d_geom = geom.targetDim();
    unsigned d_param = geom.domainDim();

    for(unsigned i = 0;i < s;i++)
       param.push_back(gsMatrix<T>(d_param, sizes[i]));

    sampler.points_sampling(param);

    total_size = 0;
    for(unsigned i = 0;i < s;i++)
    {
        /// If there is a reparameterization of the template
        /// we first sample this reparameterization
        tmp_geom.push_back(gsMatrix<T>(d_geom, sizes[i]));
        tmp_template.push_back(gsMatrix<T>(d_geom, sizes[i]));

        if(_reparam != NULL)
        {
            tmp_reparam.push_back(gsMatrix<T>(d_param, sizes[i]));
            _reparam[i].eval_into(param[i], tmp_reparam[i]);
            _geom[i].eval_into(tmp_reparam[i], tmp_geom[i]);

        } else
            _geom[i].eval_into(param[i], tmp_geom[i]);

        _templ[i].eval_into(param[i], tmp_template[i]);

        total_size += sizes[i];
    }

    res_geom.clear();
    res_template.clear();
    res_geom.push_back(gsMatrix<T>(d_geom, total_size));
    res_template.push_back(gsMatrix<T>(d_geom, total_size));

    s = tmp_geom.size();
    unsigned _j = 0;
    for(unsigned i = 0;i < s;i++)
    {
        for(int j = 0;j < tmp_geom[i].cols();j++)
        {
            res_geom[0].col(_j) = tmp_geom[i].col(j);
            res_template[0].col(_j) = tmp_template[i].col(j);
            _j++;
        }
    }
}


template<class T>
void sample_boundary_points(gsMultiBasis<T>& basis,
                            gsMultiPatch<T>& geom,
                            gsMultiPatch<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int local_size_sample,
                            gsMultiPatch<T>* temp_reparam)
{
    int* p_local_size = NULL;
    if(local_size_sample > 0)
    {
        int tab_local_size[4];
        unsigned dim = basis.domainDim();
        for(unsigned i = 0;i < dim;i++)
            tab_local_size[i] = local_size_sample;
        p_local_size = tab_local_size;
    }
    sample_boundary_points(basis, geom, templ, res_geom,
                           res_template, temp_reparam,
                           p_local_size);
}

/// Check that there exist some basis functions that are active
/// on the two corners.
bool checkActivesLocal(int                 nb_min,
                       gsMatrix<index_t>& actives_lowerCorner,
                       gsMatrix<index_t>& actives_upperCorner)
{
    index_t s_lower = actives_lowerCorner.rows();
    index_t s_upper = actives_upperCorner.rows();
    int nb = 0;
    for(index_t i = 0;i < s_upper;i++)
    {
        for(index_t j = 0;j < s_lower;j++)
        {
            if(actives_upperCorner.coeff(i,0)
               == actives_lowerCorner.coeff(j,0))
                nb++;
        }
    }
    return nb >= nb_min;
}



/// Check if the domain bound is refined enough
template<class T>
bool checkActives(gsGeometry<T>& bound, gsBasis<T>& basisDom)
{
    gsGenGeomIteratorSimplePatch<T> domIt(bound);

//    gsMatrix<T> centerPoint;
    gsMatrix<T> lowerCorner;
    gsMatrix<T> upperCorner;
//    gsMatrix<index_t> actives_centerPoint;
    gsMatrix<index_t> actives_lowerCorner;
    gsMatrix<index_t> actives_upperCorner;
    short_t dim = bound.geoDim();
    int nb_min = 1;
    for(short_t i = 0;i < dim;i++)
        nb_min *= basisDom.degree(i);

    for (; domIt.good(); domIt.next() )
    {
        /// Map the Quadrature rule to the element
        /// and compute basis derivatives
        gsBasis<T>& cur_basis =
            domIt.currentDomainBasis(basisDom);
        //    centerPoint = domIt.imageCenterPoint();
        lowerCorner = domIt.imageLowerCorner();
        upperCorner = domIt.imageUpperCorner();

        //    cur_basis.active_into(centerPoint, actives_centerPoint);
        cur_basis.active_into(lowerCorner, actives_lowerCorner);
        cur_basis.active_into(upperCorner, actives_upperCorner);

        /// We check if all the active basis functions in the center
        /// are also active in at least one of the corners
        if(!checkActivesLocal(nb_min, actives_lowerCorner,
                              actives_upperCorner))
            return false;
    }
    return true;
}


/// Refine the boundary until it has "as many degrees of freedom
/// as the domain along this boundary"
template<class T>
void refineBoundary(gsGeometry<T>& bound, gsBasis<T>& basisDom)
{
    while(! checkActives(bound, basisDom))
    {
        bound.uniformRefine();
    }
}


template<class T>
void refineBoundary(gsMultiPatch<T>& bound, gsBasis<T>& basisDom)
{
    unsigned s = bound.nPatches();
    for(unsigned i = 0;i < s;i++)
        refineBoundary(bound.patch(i), basisDom);
}

template<class T> void refineBoundary(gsMultiPatch<T>& bound,
                                      gsMultiBasis<T>& basisDom)
{
    unsigned s = bound.nPatches();
    for(unsigned i = 0;i < s;i++)
        refineBoundary(bound.patch(i), basisDom.basis(i));
}
template<class T> void refineBoundary(gsGeometry<T>& bound,
                                      gsMultiBasis<T>& basisDom)
{
    gsWarn << "Should not be here!! The case of a single patch boundary and a multipatch basis is not considered!!!" << std::endl;
}


template<class T>
void sample_boundary_points(gsMultiPatch<T>& geom,
                            gsMultiPatch<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            gsMultiPatch<T>* temp_reparam,
                            int* local_size_sample)
{
    gsMultiBasis<T> basis_tmp = copyBasis(templ);
    sample_boundary_points(basis_tmp, geom, templ,
                           res_geom, res_template,
                           temp_reparam,
                           local_size_sample);
}

template<class T>
void sample_boundary_points(gsMultiPatch<T>& geom,
                            gsMultiPatch<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int local_size_sample,
                            gsMultiPatch<T>* temp_reparam)
{
    int* p_local_size = NULL;
    if(local_size_sample > 0)
    {
        int tab_local_size[4];
        unsigned dim = geom.domainDim();
        for(unsigned i = 0;i < dim;i++)
            tab_local_size[i] = local_size_sample;
        p_local_size = tab_local_size;
    }
    sample_boundary_points(geom, templ, res_geom,
                           res_template, temp_reparam, p_local_size);
}

template<class T>
void sample_boundary_points(gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            gsGeometry<T>* temp_reparam,
                            int* local_size_sample)
{
    sample_boundary_points<T>(templ.basis(), geom, templ,
                              res_geom, res_template,
                              temp_reparam,
                              local_size_sample);
}


template<class T>
void sample_boundary_points(gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int pos, int* local_size_sample)
{
    sample_boundary_points<T>(templ.basis(), geom, templ,
                              res_geom, res_template, pos,
                              local_size_sample);

}

template<class T>
void sample_boundary_points(gsGeometry<T>& geom,
                            gsGeometry<T>& templ,
                            std::vector<gsMatrix<T> >& res_geom,
                            std::vector<gsMatrix<T> >& res_template,
                            int local_size_sample,
                            gsGeometry<T>* temp_reparam)
{
    int* p_local_size = NULL;
    if(local_size_sample > 0)
    {
        int tab_local_size[4];
        unsigned dim = geom.domainDim();
        for(unsigned i = 0;i < dim;i++)
            tab_local_size[i] = local_size_sample;
        p_local_size = tab_local_size;
    }
    sample_boundary_points(geom, templ, res_geom,
                           res_template, temp_reparam, p_local_size);
}



} /// namespace gismo
