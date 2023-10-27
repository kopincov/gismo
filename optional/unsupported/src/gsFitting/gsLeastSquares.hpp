/** @file gsLeastSquares.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsLeastSquares.h>
#include <gsFitting/gsPointContainer.hpp>
#include <gsFitting/gsFittingUtilsGen.hpp>


namespace gismo
{

template<class GenBasis, class T>
gsLeastSquares<GenBasis, T>::
gsLeastSquares(GenBasis& basis)
: gsPointContainer<T>::gsPointContainer(nBasesGen<T>(basis)),
m_basis(basis), m_output("")
{
    m_current_map = NULL;
    init_m_error();
    m_export_points = false;
}

template<class GenBasis, class T>
gsLeastSquares<GenBasis, T>::
gsLeastSquares(GenBasis& basis, const gsPointContainer<T>& points,
                std::string& output, bool export_points)
: gsPointContainer<T>::gsPointContainer(points), m_basis(basis), m_output(output)
{
    m_current_map = NULL;
    m_export_points = export_points;
    init_m_error();
}

template<class GenBasis, class T>
gsLeastSquares<GenBasis, T>::
gsLeastSquares(GenBasis& basis,
               std::vector< gsMatrix<T> >& pts,
               std::vector< gsMatrix<T> >& params,
               std::string& output, bool export_points,
               gsFunctionSet<T>* _template,
               gsFunctionSet<T>* geometry,
               index_t *local_size) : gsPointContainer<T>::
gsPointContainer(pts, params, _template, geometry, local_size),
m_basis(basis), m_output(output)
{
    m_export_points = export_points;
    init_m_error();
    m_current_map = NULL;
}


template<class GenBasis, class T>
gsLeastSquares<GenBasis, T>::
gsLeastSquares(const gsLeastSquares<GenBasis, T>& ls)
: gsPointContainer<T>::gsPointContainer(ls),
m_basis(ls.m_basis), m_output(ls.m_output)
{
    init_m_error();
    m_current_map = ls.m_current_map;
    m_export_points = ls.m_export_points;
}


template<class GenBasis, class T> void
gsLeastSquares<GenBasis, T>::init_m_error()
{
    m_errors.clear();
    for(int i = 0;i < m_size;i++)
        m_errors.push_back(gsMatrix<T>(this->size(i),1));
}

template<class GenBasis, class T>
void gsLeastSquares<GenBasis, T>::computeErrors(index_t type)
{
    gsMatrix<T> val_i;
    gsPointIterator<T> it(*this);
    std::vector< gsMatrix<T> > image;
    T err, _err;
    m_minError = std::numeric_limits<double>::max();
    m_maxError = 0.;
    m_totalError = 0.;

    GISMO_ASSERT(m_current_map != NULL,
                 "The current map should be given to compute the error");

    computeImage(m_current_map, image);

    for(;it.good();it.next())
    {
        const typename gsMatrix<T>::Column curr_point
            = it.currPoint();
        const typename gsMatrix<T>::Column curr_image
            = it.getPointUsingIterator(image);
        _err = (curr_point.col(0) -
                curr_image.col(0)).squaredNorm();
        err = get_local_error(type, _err);
        setError(it, err);
    }
}


template<class GenBasis, class T>
bool gsLeastSquares<GenBasis, T>::
actualize_degrees_freedomLS()
{
    bool refined = this->actualize_degrees_freedom(m_basis);
    if(refined)
        init_m_error();
    return refined;
}

/// Fills the matrix with the least square method
template<class GenBasis, class T>
void gsLeastSquares<GenBasis, T>::
assemble(gsFittingSystem<GenBasis, T>& system)
{
    //for computing the value of the basis function
    gsMatrix<T> value;
    gsMatrix<index_t> actives;

    gsPointIterator<T> it(*this);

    std::vector< gsMatrix<T> > results;

    int dim = this->dim_im();
    gsMatrix<T> vect_diff(dim, 1);
    unsigned _d = 1;

    gsSparseEntries<T> entries_colloc_mat;

    int patch;
    bool split_dim = system.split_dim();

    if(! split_dim)
        _d = dim;

    std::vector<T> localA(_d * _d);
    std::fill(localA.begin(), localA.end(), 0);

    if(m_export_points)
    {
        std::string _path = m_output + "_pts_";
        std::string _path2 = m_output + "_image_";
        unsigned s = m_points.size();
        for(unsigned i = 0;i < s;i++)
        {
            std::string path = _path + "patch_"
                + util::to_string(i);
            std::string path2 = _path2 + "patch_"
                + util::to_string(i);
            gsWriteParaviewPoints(m_points[i], path2);
            gsWriteParaviewPoints(m_params[i], path);
        }
    }

    if(m_current_map != NULL)
        this->eval_into(m_basis, *m_current_map, results);

    for(;it.good();it.next())
    {
        patch = it.patch();
        const typename gsMatrix<T>::Column
            curr_param = it.currParam();
        const gsMatrix<T>& tmp = curr_param;
        /// computing the values of the basis
        /// functions at the current point

        gsBasis<T>& cur_basis = it.currentBasis(m_basis);

        vect_diff = it.currPoint();

        /// Case where we consider the deformation:
        /// we subtract the identity mapping
        if(m_current_map != NULL)
            vect_diff -= it.getPointUsingIterator(results);
        add_component_SM(vect_diff, it);

        cur_basis.eval_into(tmp, value);

        // which functions have been computed i.e. which are active
        cur_basis.active_into(tmp, actives);

        const index_t numActive = actives.rows();
        for (index_t i = 0; i != numActive; ++i)
        {
            if(value(i,0) != 0.)
            {
                system.add_rhs( actives(i, 0), patch,
                                value(i, 0)
                                * vect_diff.transpose() );

                for (index_t j = i; j != numActive; ++j)
                {
                    if(value(j,0) != 0.)
                    {
                        for(unsigned d1_2 = 0;d1_2 < _d;d1_2++)
                        {
                            localA[d1_2 * _d + d1_2]
                                = value(i, 0) * value(j, 0);
                        }
                        system.add_matrix(actives(i, 0),
                                          actives(j, 0),
                                          patch, localA);
                        if(j != i)
                            system.add_matrix(actives(j, 0),
                                              actives(i, 0),
                                              patch, localA);

                    }
                }
            }
        }
    }

    /*

      int numBasis = m_B.rows();
      if(!split_dim)
      numBasis /= dim;
      m_colloc_mat.resize(numBasis, this->size_total());
      m_colloc_mat.setZero();
      m_colloc_mat.setFrom(entries_colloc_mat);
      m_colloc_mat.makeCompressed();  */

}

/*
template<class GenBasis, class T> void
gsLeastSquares<GenBasis, T>::
add_colloc_matrix(T value, int local_ind,
                  gsPointIterator<T>& it,
                  gsSparseEntries<T>& entries)
{
    int i_colloc
        = localToGlobal(local_ind, true, it);
    int global_ind = it.globalIndex();
    entries.add(i_colloc, global_ind, value);
    } */


/// Reset the points of the container
template<class GenBasis, class T> void
gsLeastSquares<GenBasis, T>::
resetPointsLS(std::vector<gsMatrix<T> >& pts,
              std::vector<gsMatrix<T> >& param)
{
    this->resetPoints(pts, param);
    resetLMult();
    init_m_error();
}

/// Reset the points of the container
template<class GenBasis, class T> void
gsLeastSquares<GenBasis, T>::
resetPointsLS(gsMatrix<T>& pts, gsMatrix<T>& param)
{
    this->resetPoints(pts, param);
    resetLMult();
    init_m_error();
}

/// Reset the points of the container
template<class GenBasis, class T> void
gsLeastSquares<GenBasis, T>::
resetPointsLS(gsPointContainer<T>& cont)
{
    this->resetPoints(cont.points(), cont.params());
    resetLMult();
    init_m_error();
}

template<class GenBasis, class T>
void gsLeastSquares<GenBasis, T>::
addPointsLS(std::vector<gsMatrix<T> >& pts,
            std::vector<gsMatrix<T> >& param)
{
    this->addPoints(pts, param);
    resetLMult();
    init_m_error();
}

template<class GenBasis, class T>
void gsLeastSquares<GenBasis, T>::
addPointsLS(gsMatrix<T>& pts, gsMatrix<T>& param,
            index_t ind)
{
    this->addPoints(pts, param, ind);
    resetLMult();
    init_m_error();

}

template<class GenBasis, class T>
void gsLeastSquares<GenBasis, T>::
addPointsLS(gsPointContainer<T>& cont)
{
    this->addPoints(cont.points(), cont.params());
    resetLMult();
    init_m_error();
}


template<class GenBasis, class T>
T gsLeastSquares<GenBasis, T>::get_local_error(index_t type, T _err)
{
    T err = 0.;
    switch (type)
    {
    case 0:
        err = _err;
        break;
    case 1:
        err = sqrt(_err);
        break;
    default:
        gsWarn << "Unknown type in get_Error(errors, type)...\n";
        break;
    }

    if ( err > m_maxError ) m_maxError = err;
    if ( err < m_minError ) m_minError = err;
    m_totalError += err;
    return err;
}
template<class GenBasis, class T>
void gsLeastSquares<GenBasis, T>::
computeImage(gsFunctionSet<T>* map,
             std::vector<gsMatrix<T> >& res)
{
    this->eval_into(getGeometry<T>(m_basis, *map), res);
}

} // namespace gismo
