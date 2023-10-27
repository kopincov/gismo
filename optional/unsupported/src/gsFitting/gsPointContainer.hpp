/** @file gsPointContainer.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsPointContainer.h>
#include <gsFitting/gsFittingUtilsSampling.hpp>


namespace gismo
{


template<class T>
gsPointContainer<T>::
gsPointContainer(const gsPointContainer<T>& it)
: m_points(it.m_points), m_params(it.m_params)
{
    init_total_points();
    m_geometry = it.m_geometry;
    m_template = it.m_template;
    m_local_size = it.m_local_size;
}

template<class T>
gsPointContainer<T>::
gsPointContainer(index_t nPatches) : m_points(nPatches), m_params(nPatches)
{
    m_local_size = NULL;
    m_template = NULL;
    m_geometry = NULL;
    init_total_points();
}


template<class T>
gsPointContainer<T>::
gsPointContainer(std::vector< gsMatrix<T> >& pts,
                 std::vector< gsMatrix<T> >& params,
                 gsFunctionSet<T>* _template,
                 gsFunctionSet<T>* geometry,
                 index_t *local_size)
: m_points(pts), m_params(params)
{
    init_total_points();
    m_template = _template;
    m_geometry = geometry;
    m_local_size = local_size;
}

/// Reset the points of the container
template<class T> void gsPointContainer<T>::
resetPoints(std::vector< gsMatrix<T> >& pts,
            std::vector< gsMatrix<T> >& param)
{
    m_points = std::vector< gsMatrix<T> >(pts);
    m_params = std::vector< gsMatrix<T> >(param);
    init_total_points();
}

template<class T> void gsPointContainer<T>::
resetPoints(gsMatrix<T>& pts, gsMatrix<T>& param)
{
    m_points.clear();
    m_params.clear();
    m_points.push_back(pts);
    m_params.push_back(param);
    init_total_points();
}

template<class T> void gsPointContainer<T>::
addPoints(std::vector< gsMatrix<T> >& pts,
          std::vector< gsMatrix<T> >& param)
{
    unsigned s = pts.size();
    GISMO_ASSERT(s == nPatches(), "Wrong number of patches");
    for(unsigned i = 0;i < s;i++)
        addPoints(pts[i], param[i], i);
    init_total_points();
}

template<class T> void gsPointContainer<T>::
addPoints(gsMatrix<T>& pts, gsMatrix<T>& param, index_t ind)
{
    GISMO_ASSERT(ind < nPatches(), "Wrong index");
    unsigned s = pts.cols();
    GISMO_ASSERT((index_t) s == param.cols(),
                 "The number of points does not coincide!!");

    gsMatrix<T>& cur_pts(m_points[ind]);
    gsMatrix<T>& cur_params(m_params[ind]);

    unsigned nPts = cur_pts.cols();

    cur_pts.conservativeResize(pts.rows(), nPts + s);
    cur_params.conservativeResize(param.rows(), nPts + s);

    for(unsigned i = 0;i < s;i++)
    {
        cur_pts.col(i + nPts) = pts.col(i);
        cur_params.col(i + nPts) = param.col(i);
    }
    init_total_points();
}


/// When the basis is refined, this function must be called to
/// resample the points
template<class T> bool gsPointContainer<T>::
actualize_degrees_freedom(gsMultiBasis<T>& new_basis)
{
    if(m_geometry == NULL || m_template == NULL)
        return false;
    gsMultiPatch<T>* patches =
        dynamic_cast<gsMultiPatch<T>*>(m_geometry);

    gsMultiPatch<T>* _template =
        static_cast<gsMultiPatch<T>*>(m_template);
    GISMO_ASSERT(patches != NULL,
                 "ERROR the object should be a gsMultiPatch");
    GISMO_ASSERT(_template != NULL,
                 "ERROR the object should be a gsMultiPatch");
    size_t size_temp = _template->nPatches();

    GISMO_ASSERT(new_basis.nBases() == size_temp,
                 "The number of patches should be the same");

    m_points.clear();
    m_params.clear();
    m_points.resize(size_temp);
    m_params.resize(size_temp);

    for(size_t i = 0;i < size_temp;i++)
    {
        refineBoundary(_template->patch(i), new_basis.basis(i));
        sample_boundary_points<T>(patches->patch(i), _template->patch(i),
                                  m_points, m_params, i,
                                  m_local_size);
    }
    init_total_points();
    return true;
}


/// When the basis is refined, this function must be called to
/// resample the points
template<class T> bool gsPointContainer<T>::
actualize_degrees_freedom(gsBasis<T>& new_basis)
{
    if(m_geometry == NULL || m_template == NULL)
        return false;
    gsGeometry<T>* geometry =
        dynamic_cast<gsGeometry<T>*>(m_geometry);

    if(geometry != NULL)
    {
        gsGeometry<T>* _template =
            static_cast<gsGeometry<T>*>(m_template);
        refineBoundary(*_template, new_basis);
        sample_boundary_points<T>(_template->basis(),
                                  *geometry, *_template,
                                  m_points, m_params);
    } else {
        gsMultiPatch<T>* patch =
            dynamic_cast<gsMultiPatch<T>*>(m_geometry);
        GISMO_ASSERT(patch != NULL, "ERROR the object is neither a gsGeometry nor a gsMultiPatch");

        gsMultiPatch<T>* _template =
            static_cast<gsMultiPatch<T>*>(m_template);
        refineBoundary(*_template, new_basis);
        sample_boundary_points<T>(*patch, *_template,
                                  m_points, m_params, NULL);
    }
    init_total_points();
    return true;
}

/// Initialize the attributes containing the number of points
template<class T>
void gsPointContainer<T>::init_total_points()
{
    unsigned total = 0;
    m_size = m_points.size();
    unsigned _size = m_size;
    m_global_index.clear();
    m_global_index.push_back(total);

    for(unsigned i = 0;i < _size;i++)
    {
        total += m_points[i].cols();
        m_global_index.push_back(total);
    }
}

template<class T> void
gsPointContainer<T>::
eval_into(gsGeometry<T>& geometry,
          std::vector< gsMatrix<T> >& res)
{
    index_t size = m_global_index[1];
    short_t d = geometry.targetDim();
    res.push_back(gsMatrix<T>(d, size));
    geometry.eval_into(m_params[0], res[0]);
}

template<class T> void
gsPointContainer<T>::
eval_into(gsBasis<T>& basis, gsFunctionSet<T>& geometry,
          std::vector< gsMatrix<T> >& res)
{
    gsGeometry<T>& _geometry =
        static_cast<gsGeometry<T>&>(geometry);
    this->eval_into(_geometry, res);
}

template<class T> void
gsPointContainer<T>::
eval_into(gsMultiPatch<T>& patches,
          std::vector< gsMatrix<T> >& res)
{
    index_t size;
    short_t d = patches.targetDim();
    for(index_t i = 0;i < m_size;i++)
    {
        size = m_global_index[i+1] - m_global_index[i];
        res.push_back(gsMatrix<T>(d, size));
        patches[i].eval_into(m_params[i], res[i]);
    }
}

template<class T> void
gsPointContainer<T>::
eval_into(gsMultiBasis<T>& basis,
          gsFunctionSet<T>& patches,
          std::vector< gsMatrix<T> >& res)
{
    gsMultiPatch<T>& _patches =
        static_cast<gsMultiPatch<T>&>(patches);
    this->eval_into(_patches, res);
}

} // namespace gismo
