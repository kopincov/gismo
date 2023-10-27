/** @file gsFittingGeometry.hpp

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#include <gsFitting/gsFittingGeometry.h>

namespace gismo
{

template<class T>
gsFittingGeometry<T>::
gsFittingGeometry(gsSolid<T>& solid)
: m_boundaries()
{
    GISMO_ASSERT(solid.nVolumes() == 1, "We only consider solids that have only one volume");

    m_boundaries.push_back
        ( std::vector< gsFittingBorder<T> >() );
    unsigned s = solid.nHalfFaces();
    for(unsigned i = 0;i < s;i++)
        m_boundaries[0].push_back
            (gsFittingBorder<T>(*solid.getHalfFaceFromID(i)->surf));
}

template<class T>
gsFittingGeometry<T>::
gsFittingGeometry(gsMultiPatch<T>& geom, bool unique_patch)
: m_boundaries()
{
    unsigned s = geom.size();
    if(! unique_patch)
    {
        for(unsigned i = 0;i < s;i++)
        {
            m_boundaries.push_back
                ( std::vector< gsFittingBorder<T> >() );
            m_boundaries[i].push_back
                (gsFittingBorder<T> (geom.patch(i) ) );
        }
    }
    else
    {
        m_boundaries.push_back(std::vector< gsFittingBorder<T> >() );
        for(unsigned i = 0;i < s;i++)
        {
            m_boundaries[0].push_back
                ( gsFittingBorder<T> (geom.patch(i) ) );
        }
    }
}

template<class T> gsFittingGeometry<T>::
gsFittingGeometry(std::vector<gsMatrix<T> >& pts, bool unique_patch)
: m_boundaries()
{
    unsigned s = pts.size();
    if(! unique_patch)
    {
        for(unsigned i = 0;i < s;i++)
        {
            m_boundaries.push_back
                ( std::vector< gsFittingBorder<T> >() );
            m_boundaries[i].push_back
                (gsFittingBorder<T> (pts[i]) );
        }
    }
    else
    {
        m_boundaries.push_back(std::vector< gsFittingBorder<T> >() );
        for(unsigned i = 0;i < s;i++)
        {
            m_boundaries[0].push_back
                ( gsFittingBorder<T> ( pts[i] ) );
        }
    }
}

template<class T> void
gsFittingGeometry<T>::getBoundingBox(gsMatrix<T>& res)
{

    GISMO_ASSERT(nPatches() == 1,
                 "Bounding boxes can only be computed when the number of patches is one!");
    if(m_bbox.size() > 0)
        res = m_bbox[0];
    unsigned s = m_boundaries[0].size();
    if(s > 0)
    {
        int dim = m_boundaries[0][0].dim_im();
        res.setZero(dim, 2);
        gsMatrix<T> bbox(dim, 2);
        m_boundaries[0][0].getBoundingBox(res);
        for(unsigned i = 1;i < s;i++)
        {
            m_boundaries[0][i].getBoundingBox(bbox);
            for(int d = 0;d < dim;d++)
            {
                if(bbox(d, 0) < res(d, 0))
                    res(d, 0) = bbox(d, 0);
                if(bbox(d, 1) > res(d, 1))
                    res(d, 1) = bbox(d, 1);
            }
        }
    }
    else
        res.setZero(3, 2);
}


template<class T> gsFittingGeometry<T>::
gsFittingGeometry(gsPlanarDomain<T>& domain)
: m_boundaries()
{
    unsigned numLoops = domain.numLoops();
    unsigned nCurves;
    m_boundaries.push_back
        ( std::vector< gsFittingBorder<T> >() );
    for(unsigned i = 0;i < numLoops;i++)
    {
        gsCurveLoop<T>& loop( domain.loop(i) );
        nCurves = loop.numCurves();
        for(unsigned j = 0;j < nCurves;j++)
            m_boundaries[0].push_back
                ( gsFittingBorder<T>(loop.curve(j)) );
    }
}

template<class T> gsFittingGeometryCoupling<T>::
gsFittingGeometryCoupling(gsFittingGeometry<T>& templ,
                          gsFittingGeometry<T>& target)
: m_template(&templ), m_target(&target), m_coupling()
{
    unsigned nPatches = templ.nPatches();
    unsigned nBorder;

    for(unsigned i = 0;i < nPatches;i++)
    {
        nBorder = target.nBorder(i);
        m_coupling.push_back(ContainerPatch());
        for(unsigned j = 0;j < nBorder;j++)
            m_coupling[i].push_back(
                gsFittingBorderCoupling<T>
                (templ.border(i, j), target.border(i, j)) );
    }
}

template<class T> template<unsigned d>
gsTensorBSplineBasis<d, T> gsFittingGeometryCoupling<T>::
getBasisTrimming(gsFittingParam<T>& param_fitting)
{
    gsMatrix<T> bbox;
    m_template->getBoundingBox(bbox);
    std::vector< gsKnotVector<T> > KV;
    for(unsigned i = 0; i != d; i++)
        KV.push_back(gsKnotVector<T> (bbox(i, 0), bbox(i, 1),
                                      param_fitting.interiorKnots,
                                      param_fitting.degree + 1));

    gsTensorBSplineBasis<d> basis( KV );
    return basis;
}


} /// namespace gismo
