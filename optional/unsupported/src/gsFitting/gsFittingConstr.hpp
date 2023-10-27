/** @file gsFittingConstr.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingConstr.h>

#include <gsFitting/gsFittingUtilsGen.hpp>
#include <gsFitting/gsFittingBase_.cpp>
#include <gsFitting/gsFittingBase.hpp>

#include <gsFitting/gsFittingIterative.hpp>
#include <gsFitting/gsFittingReinit.hpp>
#include <gsFitting/gsFittingUtilsSampling.h>
#include <gsFitting/gsFittingUtilsSampling.hpp>
#include <gsFitting/gsFittingEnergy.hpp>

#include <gsFitting/gsBoundariesCoupling.hpp>
#include <gsFitting/gsFittingGeometry.hpp>

namespace gismo
{

template<short_t d, class GenBasis,
         unsigned o, class T> void
addSmoothingQuadrature(GenBasis& basis, int dim_im,
                       gsFittingQuadrEnergy<GenBasis, T>& ener,
                       gsFittingParam<T>& param)
{
    std::vector< gsFittingQuadrature<d, GenBasis, o, T> >&
        vect_quadr(ener.template getQuadratures<d>());
    unsigned s = vect_quadr.size();
    GISMO_ASSERT(s == 0, "There already exist a smoothing?");

    vect_quadr.push_back(gsFittingQuadrature<d, GenBasis, o, T>(basis, dim_im));
    addSmoothingIntegrands<d, GenBasis, o, T>(param, vect_quadr[s]);
}


template<short_t d, class GenBasisDom,
         class GenGeomBC, unsigned o, class T> void
addBCQuadrature(GenBasisDom& basis, GenGeomBC& domain,
                GenGeomBC& geom, unsigned dim_im,
                gsFittingQuadrEnergy<GenBasisDom, T>& ener,
                int size_gauss_quadr)
{
    std::vector< gsFittingQuadrature<d-1, GenBasisDom, o, T> >&
        vect_quadr ( ener.template getQuadratures<d-1>() );
    unsigned s = vect_quadr.size();
    GISMO_ASSERT(s == 0, "There already exist a fitting?");

    int *_size_gauss_quadr = NULL;
    int tab_size_gauss_quadr[d-1];
    if(size_gauss_quadr > 0)
    {
        for(unsigned i = 0;i < d-1;i++)
            tab_size_gauss_quadr[i] = size_gauss_quadr;
        _size_gauss_quadr = tab_size_gauss_quadr;
    }

    vect_quadr.push_back(gsFittingQuadrature<d-1, GenBasisDom, o, T>(basis, dim_im, &domain, _size_gauss_quadr));

    vect_quadr[s].addIntegrand( new gsFittingIntegrandL2Dist
                                <d-1, GenGeomBC, T>(dim_im, &geom) );
}


template<class GenBasis, class T> void
add_winslow_smoothing(gsFittingQuadrature<2, GenBasis, 2, T>&
                      quadr, T coeff_NL, T coeff_winslow)
{
    quadr.addIntegrand( new gsFittingIntegrandW<T>
                        (coeff_NL, coeff_winslow) );
}

template<class GenBasis, class T> void
add_winslow_smoothing(gsFittingQuadrature<3, GenBasis, 2, T>&
                      quadr, T coeff_NL, T coeff_winslow)
{
    quadr.addIntegrand( new gsFittingIntegrandW3D<T>
                        (coeff_NL, coeff_winslow) );
}

template<short_t d, class GenBasis, unsigned o, class T> void
addSmoothingIntegrands(gsFittingParam<T>& param,
                       gsFittingQuadrature<d, GenBasis,
                       o, T>& quadr_ener)
{
    if(param.coeff_linear_global > 0.)
    {
        GISMO_ASSERT(param.coeff_linear_gradient > 0
                     || param.coeff_linear_hessian,
            "The global linear coefficient is positive while all linear coefficients are null. This might lead to an error.");
        quadr_ener.addIntegrand
            ( new gsFittingIntegrandLin<d, o, T>
              (param.coeff_linear_global,
               param.coeff_linear_gradient,
               param.coeff_linear_hessian, true) );
    }
    if(param.coeff_NL_global > 0.)
    {
        GISMO_ASSERT(param.coeff_NL_winslow > 0.
                     || param.coeff_NL_metric > 0.,
                     "The global nonlinear coefficient is positive while all nonlinear coefficients are null. This might lead to an error.");
        add_NL_smoothing<d, GenBasis, o, T>(quadr_ener, param);
    }
}


template<short_t d, class GenBasis, unsigned o, class T> void
add_NL_smoothing(gsFittingQuadrature<d, GenBasis, o, T>& quadr,
                 gsFittingParam<T>& param)
{
    if(param.coeff_NL_metric > 0)
    {
        quadr.addIntegrand( new gsFittingIntegrandOrth<d, T>
                            (param.coeff_NL_global,
                             param.coeff_NL_metric));
    }
    if(param.coeff_NL_winslow > 0)
        add_winslow_smoothing<GenBasis, T>
            (quadr, param.coeff_NL_global,
             param.coeff_NL_winslow);
}


template<short_t d, class GenBasis, unsigned o, class T>
gsFunctionSet<T>*
fitPoints(std::vector< gsMatrix<T> >& points,
          std::vector< gsMatrix<T> >& parameters,
          GenBasis& basis, gsFittingParam<T>& param,
          gsFunctionSet<T>* templ, gsFunctionSet<T>* geom)
{
    unsigned dim_im = points[0].rows();
    gsFittingInitializerLS<d, GenBasis, T>
        initializer(points, parameters, templ, geom);
    return fitData(dim_im, basis, param, initializer);
}

/// basisDom: the basis of the template domain
/// (which is different from the basis of the boundary conditions)
template<short_t d, class GenBasisDom, class GenGeomBC,
         unsigned o, class T> gsFunctionSet<T>*
fitGeometry(GenGeomBC* geom, GenGeomBC* temp,
            GenBasisDom& basisDom, gsFittingParam<T>& param)
{
    if(!param.continuous_fitting)
    {
        std::vector< gsMatrix<T> > vect_pts;
        std::vector< gsMatrix<T> > vect_param;
        refineBoundary(*temp, basisDom);

        int nbPatches = nPatchesGen<T>(*geom);
        std::vector<int> bc = indices_bc
            (nbPatches, param.fixed_boundaries);

        sample_boundary_points<T>(*geom, *temp,
                                  vect_pts, vect_param);

        return fitPoints<d, GenBasisDom, o, T>
            (vect_pts, vect_param, basisDom, param,
             temp, geom);

    } else {
        unsigned dim_im = geom->geoDim();
        gsFittingInitializerContinuous
            <d, GenBasisDom, GenGeomBC, o, T> initializer
            (geom, temp);
        return fitData(dim_im, basisDom, param, initializer);
    }
}


template<short_t d, class GenBasis, class T> gsFunctionSet<T>*
fitData(int dim_im, GenBasis& basis, gsFittingParam<T>& param,
        gsFittingInitializer<d, GenBasis, T>& initializer)
{
    bool remove_result = false;

    gsFittingReinit<d, GenBasis, T>
        reinit(param, param.threshold_reinit,
               param.prop_reinit);


    gsFittingHierarchical<d, GenBasis, T>* p_hIter = NULL;
    if(param.use_refinement)
    {
        p_hIter = new gsFittingHierarchical
            <d, GenBasis, T>(basis, param.threshold_refinement,
                             param.extension_refinement,
                             param.max_num_iter_refin);
    }
    gsFittingAdaptSmoo<d, GenBasis, T>
        subPb(param.coeff_reducing_smoothing,
              param.inf_smoothing_parameter,
              param.max_num_reducing,
              param.coeff_linear_global,
              param.coeff_NL_global,
              param.print_messages, p_hIter);
    gsFittingImproveAfterCV<d, GenBasis, T>* p_iter = NULL;
    if(param.use_reducing_smoo)
        p_iter = &subPb;
    else
        p_iter = p_hIter;
    gsFittingIterNL <d, GenBasis, T>
        fitting(dim_im, basis, remove_result,
                param, p_iter, &reinit);

    addSmoothingQuadrature<d, GenBasis, 2, T>
        (basis, dim_im, fitting.quadrEner(),
         param);

    if(param.coeff_tikhonov > 0.)
        fitting.add_regularization_term
        (param.coeff_tikhonov);

    initializer.init(fitting);
    fitting.compute();
    return fitting.result();
}


std::vector<int> indices_bc(int nPatches,
                            const std::vector<int>& indices)
{
    int del = 1 + 2*nPatches;
    std::vector<int> res;
    unsigned s = indices.size();
    for(unsigned i = 0;i < s;i += del)
    {
        std::vector<int>::iterator it
            = std::find (res.begin(), res.end(), indices[i]);
        if(it == res.end())
            res.push_back(indices[i]);
    }
    return res;
}



template<short_t d, class GenBasis, class T> gsFunctionSet<T>*
fitData(int dim_im, GenBasis& basis, gsFittingParam<T>& param,
        gsFittingGeometryCoupling<T>& coupling)
{
    bool remove_result = false;

    gsFittingReinit<d, GenBasis, T>
        reinit(param, param.threshold_reinit,
               param.prop_reinit);

    gsFittingHierarchical<d, GenBasis, T>* p_hIter = NULL;
    if(param.use_refinement)
    {
        p_hIter = new gsFittingHierarchical
            <d, GenBasis, T>(basis, param.threshold_refinement,
                             param.extension_refinement,
                             param.max_num_iter_refin);
    }
    gsFittingAdaptSmoo<d, GenBasis, T>
        subPb(param.coeff_reducing_smoothing,
              param.inf_smoothing_parameter,
              param.max_num_reducing,
              param.coeff_linear_global,
              param.coeff_NL_global,
              param.print_messages, p_hIter);
    gsFittingImproveAfterCV<d, GenBasis, T>* p_iter = NULL;
    if(param.use_reducing_smoo)
        p_iter = &subPb;
    else
        p_iter = p_hIter;
    gsFittingIterNL <d, GenBasis, T>
        fitting(dim_im, basis, remove_result,
                param, p_iter, &reinit);

    addSmoothingQuadrature<d, GenBasis, 2, T>
        (basis, dim_im, fitting.quadrEner(),
         param);

    if(param.coeff_tikhonov > 0.)
        fitting.add_regularization_term
        (param.coeff_tikhonov);

    /// Initialize the data that must be initialized
    coupling.template initialize<GenBasis>(param, basis);

    /// Sets the constraints to the system
    coupling.template associate<d, GenBasis>(fitting);
    fitting.compute();
    return fitting.result();
}

} /// namespace gismo
