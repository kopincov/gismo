/** @file gsFittingRoutines.hpp


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingRoutines.h>
#include <gsFitting/gsFittingUtilsInit.h>

#include <gsFitting/gsFittingIterative.hpp>
#include <gsFitting/gsFittingUtilsIO.hpp>
#include <gsFitting/gsFittingUtilsSolids.hpp>
#include <gsFitting/gsTemplateTargetFitting.hpp>

namespace gismo
{

template <short_t d, class GenGeom, class T>
int fittingTargetSinglePatch(gsFittingParam<T> &param,
                             GenGeom& geom, GenGeom& templ)
{
    gsMatrix<T> bbox;
    setBoundingBox(templ, bbox);
    std::vector< gsKnotVector<T> > KV;
    for(unsigned i = 0; i != d; i++)
    {
        KV.push_back(gsKnotVector<T> (bbox(i, 0), bbox(i, 1),
                                      param.interiorKnots,
                                      param.degree + 1 ));
    }
    gsBasis<T>* p_basis = NULL;
    gsTensorBSplineBasis<d> basis( KV );
    gsTHBSplineBasis<d> thb(basis);
    if(param.use_refinement)
        p_basis = &thb;
    else
        p_basis = &basis;

    gsFunctionSet<T> *result = fitGeometry
        <d, gsBasis<T>, GenGeom, 2, T>
        (&geom, &templ, *p_basis, param);
    GISMO_ENSURE(result != NULL, "Error during the construction of the parameterization");
    exportGeometry<gsGeometry<T> >(result, param.output, true);

    return 0;
}


template <short_t d, class GenGeom, class T>
int fittingTargetMultiPatch(gsFittingParam<T> &param,
                            gsMultiBasis<T> &basis,
                            GenGeom& geom, GenGeom& templ)
{
    gsFunctionSet<T> *result = fitGeometry
        <d, gsMultiBasis<T>, GenGeom, 2, T>
        (&geom, &templ, basis, param);
    GISMO_ENSURE(result != NULL, "Error during the construction of the parameterization");
    exportGeometry<gsMultiPatch<T> >(result, param.output, true);

    return 0;
}


template <short_t d, class T>
int fittingPointsMultiPatch(gsFittingParam<T> &fitting_param,
                            gsMultiBasis<T>& basis,
                            std::vector<gsMatrix<T> > &pts,
                            std::vector<gsMatrix<T> > &param)
{
    gsFunctionSet<T> *result = fitPoints<d, gsMultiBasis<T>, 2, T>
        (pts, param, basis, fitting_param, NULL, NULL);

    GISMO_ENSURE(result != NULL, "Error during the construction of the parameterization");
    exportGeometry<gsMultiPatch<T> >(result, fitting_param.output,
                                     true);

    return 0;
}

template <short_t d, class T>
int fittingPointsSinglePatch(gsFittingParam<T> &fitting_param,
                             std::vector<gsMatrix<T> > &pts,
                             std::vector<gsMatrix<T> > &param)
{
    gsFunctionSet<T> *result = NULL;

    gsMatrix<T> bbox(d, 2);
    getBoundingBox(param[0], bbox);
    gsBasis<T>* p_basis = NULL;

    std::vector< gsKnotVector<T> > KV;
    for(unsigned i = 0; i != d; i++)
    {
        KV.push_back(gsKnotVector<T> ( bbox(i, 0), bbox(i, 1),
                                       fitting_param.interiorKnots,
                                       fitting_param.degree + 1 ));
    }
    gsTensorBSplineBasis<d> basis( KV );
    gsTHBSplineBasis<d> thb(basis);
    if(fitting_param.use_refinement)
        p_basis = &thb;
    else
        p_basis = &basis;

    result = fitPoints<d, gsBasis<T>, 2, T>
        (pts, param, *p_basis, fitting_param, NULL, NULL);

    GISMO_ENSURE(result != NULL, "Error during the construction of the parameterization");
    exportGeometry<gsGeometry<T> >(result, fitting_param.output,
                                   true);

    return 0;
}


template<class T> int
fittingSolid(gsSolid<T>* temp, gsSolid<T>* geom,
             gsFittingParam<T>& fitting_param)
{
    unsigned nPointsSample = 20;
    gsFunctionSet<T>* result = fitSolidFixedParam
        (temp, geom, fitting_param, nPointsSample);

    GISMO_ENSURE(result != NULL, "Error during the construction of the parameterization");
    exportGeometry<gsGeometry<T> >(result, fitting_param.output,
                                   true, 3);
    return 0;
}

template <class GenGeom, class T>
int fittingTargetDyn(gsFittingParam<T> &param,
                     GenGeom& geom, GenGeom& templ,
                     gsMultiPatch<T>* topology, gsMultiBasis<T>* basis)
{
    int dim = templ.geoDim();
    gsMultiBasis<T> new_basis;

    if(topology != NULL)
    {
        basesFromMultiPatchDyn<T>(*topology, new_basis,
                                  param.keep_basis, param.degree,
                                  param.interiorKnots,
                                  param.use_refinement);
        basis = &new_basis;
    } else if(basis != NULL)
    {
        copyBasesDyn<T>(*basis, new_basis,
                        param.keep_basis, param.degree,
                        param.interiorKnots,
                        param.use_refinement);
        basis = &new_basis;
    }

    if(dim == 2)
    {
        if(basis == NULL)
            return fittingTargetSinglePatch
                <2, GenGeom, T>(param, geom, templ);
        else
        {
            return fittingTargetMultiPatch
                <2, GenGeom, T>(param, *basis, geom, templ);
        }

    } else
    {
        GISMO_ENSURE(dim == 3, "Only the cases d=2 and d=3 are considered");
        if(basis == NULL)
            return fittingTargetSinglePatch
                <3, GenGeom, T>(param, geom, templ);
        else
            return fittingTargetMultiPatch
                <3, GenGeom, T>(param, *basis, geom, templ);
    }
}


template <class T>
int fittingPointsDyn(gsFittingParam<T> &fitting_param,
                     std::vector< gsMatrix<T> > &pts,
                     std::vector< gsMatrix<T> > &param,
                     gsMultiPatch<T>* topology,
                     gsMultiBasis<T>* basis)
{
    GISMO_ASSERT(param.size() > 0, "Error: no points given");
    gsMultiBasis<T> new_basis;
    if(topology != NULL)
    {
        basesFromMultiPatchDyn<T>
            (*topology, new_basis, fitting_param.keep_basis,
             fitting_param.degree, fitting_param.interiorKnots,
             fitting_param.use_refinement);
        basis = &new_basis;
    } else if(basis != NULL)
    {
        if(fitting_param.use_refinement
           || ! fitting_param.keep_basis)
        {
            copyBasesDyn<T>(*basis, new_basis,
                            fitting_param.keep_basis,
                            fitting_param.degree,
                            fitting_param.interiorKnots,
                            fitting_param.use_refinement);
            basis = &new_basis;
        }
    }

    if(basis != NULL)
    {
        GISMO_ENSURE(basis->nBases() == pts.size(),
                     "Number of patches different in the points and in the multibasis");
    }
    int dim = param[0].rows();
    if(dim == 2)
    {
        if(basis == NULL)
            return fittingPointsSinglePatch
                <2, T>(fitting_param, pts, param);
        else
            return fittingPointsMultiPatch
                <2, T>(fitting_param, *basis, pts, param);

    }
    else
    {
        GISMO_ENSURE(dim == 3, "Only the cases d=2 and d=3 are considered");
        if(basis == NULL)
            return fittingPointsSinglePatch
                <3, T>(fitting_param, pts, param);
        else
            return fittingPointsMultiPatch
                <3, T>(fitting_param, *basis, pts, param);
    }
}


template<class T>
int fittingTrimming(gsMultiPatch<T>* topology,
                    gsMultiBasis<T>* basis,
                    gsTemplateTargetFitting<T>& data)
{
    gsFunctionSet<T> *result = NULL;
    gsFittingParam<T>& fitting_param(data.fitting_param());
    bool single_patch = true;

    GISMO_ASSERT(data.nPatches() > 0, "Error: no points given");
    gsMultiBasis<T> new_basis;
    if(topology != NULL)
    {
        single_patch = false;
        basesFromMultiPatchDyn<T>
            (*topology, new_basis, fitting_param.keep_basis,
             fitting_param.degree, fitting_param.interiorKnots,
             fitting_param.use_refinement);
        basis = &new_basis;
    } else if(basis != NULL)
    {
        single_patch = false;
        if(fitting_param.use_refinement
           || ! fitting_param.keep_basis)
        {
            copyBasesDyn<T>(*basis, new_basis,
                            fitting_param.keep_basis,
                            fitting_param.degree,
                            fitting_param.interiorKnots,
                            fitting_param.use_refinement);
            basis = &new_basis;
        }
    }
    result = data.computeMappingTrimmingDyn(basis);
    if(single_patch)
        exportGeometry< gsGeometry<T> >
            (result, fitting_param.output, true);
    else
        exportGeometry< gsMultiPatch<T> >
            (result, fitting_param.output, true);
    return 0;
}

} /// namespace gismo
