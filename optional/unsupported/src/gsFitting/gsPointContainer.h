/** @file gsPointContainer.h

    @brief Contains the class gsPointContainer. This class contains
    the points and their associated parameters used for least
    squares method

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsFittingIterator.h>

#include <gsCore/gsLinearAlgebra.h>
#include <gsMatrix/gsMatrix.h>

namespace gismo
{


/* @brief
   Class gsPointContainer: contains the points q_i and parameters p_i
   to be fitted when constructing a (multipatch) BSpline S
   satisfying S(p_i) = q_i. These points are possibly associated to
   different patches of the BSpline.
   In the case where these points are obtained by sampling, the template
   and target domains (containing respectively the points q_i and p_i)
   can be used to resample the points in case the basis is refined.
   TODO: in case the points have not been obtained by sampling,
      we should resample using a (linear, quadratic?) interpolation
*/
template<class T>
class gsPointContainer
{
public:
    /// Empty constructor
    gsPointContainer(index_t nPatches);

    /// Constructor: we copy the points and the parameters.
    /// The template and the geometry can also be set in order
    /// to be able to resample if the basis is refined.
    gsPointContainer(std::vector< gsMatrix<T> >& pts,
                     std::vector< gsMatrix<T> >& params,
                     gsFunctionSet<T>* _template = NULL,
                     gsFunctionSet<T>* geometry = NULL,
                                        index_t *local_size = NULL);

    /// Copy a point container
    gsPointContainer(const gsPointContainer<T>& it);


    /// Reset the points of the container
    void resetPoints(std::vector<gsMatrix<T> >& pts,
                     std::vector<gsMatrix<T> >& param);
    void resetPoints(gsMatrix<T>& pts, gsMatrix<T>& param);

    void addPoints(std::vector<gsMatrix<T> >& pts,
                   std::vector<gsMatrix<T> >& param);

    void addPoints(gsMatrix<T>& pts, gsMatrix<T>& param,
                   index_t ind);

    /// Reset the geometries associated with the points
    void resetGeom(gsFunctionSet<T>* _template = NULL,
                     gsFunctionSet<T>* geometry = NULL)
    {  m_template = _template;   m_geometry = geometry;  }

    /// In the case where the basis is refined,
    /// we refine the boundary and resample the points
    /// to take into account this new precision.
    bool actualize_degrees_freedom(gsBasis<T>& new_basis);
    bool actualize_degrees_freedom(gsMultiBasis<T>& new_basis);

    /// Returns the global index associated with a point
    /// of a given patch
    inline                    index_t globalIndex(                   index_t patch,                    index_t pos)
    {
        return m_global_index[patch] + pos;
    }

    /// Return the table containing the cumulated number of points
    inline std::vector<index_t> tab_global_indices()
    { return m_global_index; }


    /// Returns the dimension of the space if there is
    /// at least one point in the set and -1 otherwise
    inline                    index_t dimensionGeom()
    {
        if(m_global_index[m_size] > 0) {
            return m_points[0].rows();
        }
        else {
            return -1;
        }
    }

    /// Compute the number of points in each patches.
    /// The result is set in m_global_index
    void init_total_points();

    /// Returns the points (in the target domain)
    inline std::vector< gsMatrix<T> >& points(){ return m_points; }

    /// Returns the parameters (in the template domain)
    inline std::vector< gsMatrix<T> >& params(){ return m_params; }


    /// Return the parameter at the current position
    inline typename gsMatrix<T>::Column param(                   index_t patch,                    index_t pos)
    {
        GISMO_ASSERT(patch < nPatches() && pos < size(patch),
                     "ERROR: index outside of the range");
        return m_params[patch].col(pos);
    }

    /// Returns the point at the current position if m_is_good
    inline typename gsMatrix<T>::Column point(                   index_t patch,                    index_t pos)
    {
        GISMO_ASSERT(patch < nPatches() && pos < size(patch),
                     "ERROR: index outside of the range");
        return m_points[patch].col(pos);
    }

    inline                    index_t dim_im()
    {
        GISMO_ASSERT(nPatches() > 0, "ERROR: no points");
        return m_points[0].rows();

    }

    /// Returns the total number of points
    inline index_t size_total(){ return m_global_index[m_size];   }

    /// Returns the number of points in a given patch
    inline index_t size(                   index_t patch)
    { return m_points[patch].cols();   }

    /// Returns the number of patches
    inline index_t nPatches(){ return m_size;   }

    /// Evals the geometry at each points of the container
    /// The result is set in res
    void eval_into(gsBasis<T>& basis, gsFunctionSet<T>& geometry,
                   std::vector< gsMatrix<T> >& res);
    void eval_into(gsGeometry<T>& geometry,
                   std::vector< gsMatrix<T> >& res);

    void eval_into(gsMultiBasis<T>& basis,
                   gsFunctionSet<T>& patches,
                   std::vector< gsMatrix<T> >& res);
    void eval_into(gsMultiPatch<T>& patches,
                   std::vector< gsMatrix<T> >& res);

    bool isEmpty()
    {   return m_points.size() == 0;  }

protected:
    /// the points of the point cloud (in the target domain)
    std::vector<gsMatrix<T> > m_points;
    /// parameters associated to the points of the point cloud
    /// (in the template domain)
    std::vector<gsMatrix<T> > m_params;

    /// the number of patches
    index_t m_size;
    /// the accumulated total number of points
    std::vector<index_t> m_global_index;

    /// Case where the points to be fitted have been obtained
    /// by a sampling. Number of points to be chosen in each
    /// interval for the resampling
    int *m_local_size;

    /// Case where the points to be fitted have been obtained
    /// by a sampling. A pointor to the target domain
    gsFunctionSet<T>* m_geometry;

    /// Case where the points to be fitted have been obtained
    /// by a sampling. A pointor to the template domain
    gsFunctionSet<T>* m_template;
};


}
