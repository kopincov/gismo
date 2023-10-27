/** @file gsFittingUtilsInit.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingUtilsInit.h>

#include <gsFitting/gsFittingIterative.hpp>
#include <gsFitting/gsFittingUtilsIO.hpp>
#include <gsFitting/gsFittingUtilsSolids.hpp>

#include <gsFitting/gsFittingBase.h>
#include <gsFitting/gsFittingBase_.cpp>
#include <gsFitting/gsFittingBase.hpp>

namespace gismo
{

template <class T>
void junctionPoints(std::vector<gsMatrix<T> > &pts,
                    std::vector<gsMatrix<T> > &param,
                    std::vector<gsMatrix<T> > &pts_res,
                    std::vector<gsMatrix<T> > &param_res,
                    gsMatrix<T>* bb_param)
{
    int s_total = 0;
    unsigned nPatch = pts.size();
    int dim_pts = pts[0].rows();
    int dim_param = param[0].rows();
    unsigned s_mat;
    int ind;

    pts_res.clear();
    param_res.clear();
    for(unsigned i = 0;i < nPatch;i++)
        s_total += pts[i].cols();
    pts_res.push_back(gsMatrix<T>(dim_pts, s_total) );
    param_res.push_back(gsMatrix<T>(dim_param, s_total) );

    ind = 0;
    for(unsigned i = 0;i < nPatch;i++)
    {
        s_mat = pts[i].cols();
        for(unsigned j = 0;j < s_mat;j++)
        {
            if(bb_param == NULL
               || isIncluded<T>(param[i].col(j), *bb_param))
            {
                pts_res[0].col(ind) = pts[i].col(j);
                param_res[0].col(ind) = param[i].col(j);
                ind++;
            }
        }
    }
}


template<class T> bool
isIncluded(const typename gsMatrix<T>::Column& vect,
           gsMatrix<T>& bb)
{
    unsigned size = vect.rows();
    for(unsigned i = 0;i < size;i++)
    {
        if(vect(i, 0) > bb(i, 1) || vect(i, 0) < bb(i, 0))
            return false;
    }
    return true;
}


template <class T>
void junctionPoints(std::vector< std::vector< gsMatrix<T> > > &pts,
                    std::vector< std::vector< gsMatrix<T> > > &param,
                    std::vector< gsMatrix<T> > &pts_res,
                    std::vector< gsMatrix<T> > &param_res,
                    gsMatrix<T>* bb_param = NULL)
{
    int s_total = 0;
    unsigned nPatch = pts.size();
    int dim_pts = pts[0][0].rows();
    int dim_param = param[0][0].rows();
    unsigned s_mat;
    int ind;

    pts_res.clear();
    param_res.clear();
    for(unsigned i = 0;i < nPatch;i++)
        s_total += pts[i][0].cols();
    pts_res.push_back(gsMatrix<T>(dim_pts, s_total) );
    param_res.push_back(gsMatrix<T>(dim_param, s_total) );

    ind = 0;
    for(unsigned i = 0;i < nPatch;i++)
    {
        s_mat = pts[i][0].cols();
        for(unsigned j = 0;j < s_mat;j++)
        {
            if(bb_param == NULL || isIncluded<T>
               (param[i][0].col(j), *bb_param))
            {
                pts_res[0].col(ind) = pts[i][0].col(j);
                param_res[0].col(ind) = param[i][0].col(j);
                ind++;
            }
            else
                gsInfo << "DBG: ENTER" << std::endl;
        }
    }
}


template<class T> void
basesFromMultiPatchDyn(gsMultiPatch<T>& mp, gsMultiBasis<T>& result,
                       bool keep_basis, int degree,
                       int num_internal_knots, bool is_thb)
{
    int dim = mp.dim();
    if(dim == 2)
    {
        basesFromMultiPatch<2, T>(mp, result, keep_basis, degree,
                                         num_internal_knots, is_thb);
    }
    else
    {
        GISMO_ASSERT(dim == 3, "Only the cases d=2 and d=3 are considered");
        basesFromMultiPatch<3, T>(mp, result, degree,
                                  num_internal_knots, is_thb);

    }
}

template<class T> void
copyBasesDyn(gsMultiBasis<T>& basis,
             gsMultiBasis<T>& result, bool keep_basis,
             int degree, int num_internal_knots, bool is_thb)
{
    int dim = basis.dim();
    if(dim == 2)
    {
        copyBases<2, T>(basis, result, keep_basis, degree,
                        num_internal_knots, is_thb);
    }
    else
    {
        GISMO_ASSERT(dim == 3, "Only the cases d=2 and d=3 are considered");
        copyBases<3, T>(basis, result, keep_basis, degree,
                               num_internal_knots, is_thb);

    }
}


template<short_t d, class T> void
copyBases(gsMultiBasis<T>& bases, gsMultiBasis<T>& result,
          bool keep_basis, int degree, int num_internal_knots,
          bool is_thb)
{
    typename gsMultiBasis<T>::BasisContainer newBases;
    unsigned nPatches = bases.nBases();

    if(keep_basis)
    {
        typename gsMultiBasis<T>::BasisContainer newBases;
        newBases.clear();
        for (unsigned i = 0; i != nPatches; i++)
        {
            if(is_thb)
                newBases.push_back( new gsTHBSplineBasis<d>
                                    ( bases.basis(i) ) );
            else
                newBases.push_back
                    ( bases.basis(i).clone().release() );
        }
        result = gsMultiBasis<T>(newBases, bases.topology());
    }
    else
    {
        typename gsMultiBasis<T>::BasisContainer newBases;
        std::vector< gsKnotVector<T> > KV;
        for (unsigned i = 0; i != nPatches; i++)
        {
            const gsMatrix<T> param = bases.basis(i).support();
            KV.clear();
            for(unsigned dim = 0;dim < d;dim++)
            {
                KV.push_back(gsKnotVector<T> (param(dim, 0),
                                              param(dim, 1),
                                              num_internal_knots,
                                              degree + 1) );
            }
            if(is_thb)
            {
                gsTensorBSplineBasis<d> tmp(KV);
                newBases.push_back(new gsTHBSplineBasis<d>(tmp));
            }
            else
                newBases.push_back(new gsTensorBSplineBasis<d>(KV));
        }
        result = gsMultiBasis<T>(newBases, bases.topology());
    }

}


template<short_t d, class T> void
basesFromMultiPatch(gsMultiPatch<T>& mp, gsMultiBasis<T>& result,
                    bool keep_basis, int degree,
                    int num_internal_knots, bool is_thb)
{

    unsigned nPatches = mp.nPatches();
    std::vector< patchSide > boundaries = mp.boundaries();
    std::vector< boundaryInterface > interfaces = mp.interfaces();

    gsBoxTopology topology(d, nPatches, boundaries, interfaces);

    if(keep_basis)
    {
        if(is_thb)
        {
            typename gsMultiBasis<T>::BasisContainer newBases;
            newBases.clear();
            for (unsigned i = 0; i != nPatches; i++)
                newBases.push_back( new gsTHBSplineBasis<d>( mp.basis(i) ) );
            result = gsMultiBasis<T>(newBases, topology);
        }
        else
        {
            typename gsMultiBasis<T>::BasisContainer newBases
                = mp.basesCopy();
            result = gsMultiBasis<T>(newBases, topology);
        }
    }
    else
    {
        typename gsMultiBasis<T>::BasisContainer newBases;
        std::vector< gsKnotVector<T> > KV;
        for (unsigned i = 0; i != nPatches; i++)
        {
            const gsMatrix<T> param = mp.parameterRange(i);
            KV.clear();
            for(unsigned dim = 0;dim < d;dim++)
            {
                KV.push_back(gsKnotVector<T> (param(dim, 0),
                                              param(dim, 1),
                                              num_internal_knots,
                                              degree + 1) );
            }
            if(is_thb)
            {
                gsTensorBSplineBasis<d> tmp(KV);
                newBases.push_back(new gsTHBSplineBasis<d>(tmp));
            }
            else
                newBases.push_back(new gsTensorBSplineBasis<d>(KV));
        }
        result = gsMultiBasis<T>(newBases, topology);
    }

}

template<class T>
void getBoundingBoxPts(gsMatrix<T>& pts, gsMatrix<T>& res)
{
    unsigned size = pts.cols();
    unsigned dim = pts.rows();
    for(unsigned d = 0;d < dim;d++)
    {
        res(d,0) = std::numeric_limits<T>::max();
        res(d,1) = -std::numeric_limits<T>::max();
    }

    for(unsigned i = 0;i < size;i++)
    {
        for(unsigned d = 0;d < dim;d++)
        {
            if(pts(d, i) < res(d, 0))
                res(d, 0) = pts(d, i);
            if(pts(d, i) > res(d, 1))
                res(d, 1) = pts(d, i);
        }
    }
}

template<class T> void
boundingBox(gsGeometry<T>& geom, gsMatrix<T> & result)
{
    result.col(0) = geom.coefs().colwise().minCoeff().transpose();
    result.col(1) = geom.coefs().colwise().maxCoeff().transpose();
}


} /// namespace gismo
