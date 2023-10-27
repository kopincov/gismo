/** @file gsLeastSquares.h

    @brief Provides the class gsLeastSquares used for
    constructing the matrix and the RHS associated with the
    least squares method

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsFittingIterator.h>
#include <gsFitting/gsLocalGlobal.h>

namespace gismo
{

/*
  @brief Class gsLeastSquares: we assemble here the matrix and the RHS
  associated with the least squares method. In the case where a
  nonlinear energy is minimized, this energy can be "developped" near a
  current position (solution).
  Note that the boundary can be multipatch and the domain single patch. In that case all the
  points are included in the same element of the vector.
*/
template<class GenBasis, class T = real_t>
class gsLeastSquares : public gsPointContainer<T>
{

public:
    using gsPointContainer<T>::m_points;
    using gsPointContainer<T>::m_params;

    typedef gsPointContainer<T> Base;
    /// constructors

    /// Empty constructor
    gsLeastSquares(GenBasis& basis);

    /// Copy the points contained in the point container
    gsLeastSquares(GenBasis& basis, const gsPointContainer<T>& points,
                   std::string& output, bool export_points);
    /// Copy the points. The target and the template domains can be
    /// set in order to resample each time the basis is modified
    /// (see gsPointContainer for more details)
    gsLeastSquares(GenBasis& basis,
                   std::vector< gsMatrix<T> >& pts,
                   std::vector< gsMatrix<T> >& params,
                   std::string& output, bool export_points,
                   gsFunctionSet<T>* _template = NULL,
                   gsFunctionSet<T>* geometry = NULL,
                   index_t *local_size = NULL);

    /// Copy a gsLeastSquares
    gsLeastSquares(const gsLeastSquares<GenBasis, T>& it);

    /// Destructor
    virtual ~gsLeastSquares(){  }

    /// Add the elements associated with the
    /// least squares method to a metrix and a RHS.
    /// Note that the matrix is represented here by a
    /// gsSparseEntries to speed up the construction.
    void assemble(gsFittingSystem<GenBasis, T>& system);


    /// In case where the basis is refined,
    /// we refine the template and target geometries
    /// and resample the points
    /// to take into account this new precision.
    /// Return true if the points have been resampled
    virtual bool actualize_degrees_freedomLS();

    /// Reset the points of the container
    void resetPointsLS(std::vector<gsMatrix<T> >& pts,
                       std::vector<gsMatrix<T> >& param);
    void resetPointsLS(gsMatrix<T>& pts, gsMatrix<T>& param);
    void resetPointsLS(gsPointContainer<T>& cont);

    /// Adds some points of the container
    void addPointsLS(std::vector<gsMatrix<T> >& pts,
                       std::vector<gsMatrix<T> >& param);
    void addPointsLS(gsMatrix<T>& pts, gsMatrix<T>& param,
                     index_t ind);
    void addPointsLS(gsPointContainer<T>& cont);

    /// In case we use Lagrange multiplier,
    /// adds these Lagrange multiplier to the RHS
    virtual void add_component_SM(gsMatrix<T>& vect,
                                  gsPointIterator<T>& it)
    { }

    /// In case we rescale the energy, this function is called
    /// to rescale the Lagrange multipliers
    virtual void decreaseLMult(T coeff)
    { }

    virtual void resetLMult()
    { }


    /// Sets the current position
    void set_current_map(gsFunctionSet<T>* current_map)
    {  m_current_map = current_map;  }

    /// Returns the number of patches
    index_t nPatches(){  return m_points.size();   }

    /// Returns the maximal, minimal and total error of the
    /// current mapping respectively.
    /// !!! computeErrors must be called before these functions
    /// are called
    inline T maxError(){ return m_maxError;   }
    inline T minError(){ return m_minError;   }
    inline T totalError(){ return m_totalError;   }

    /// Returns the errors of the current mapping.
    /// !!! computeErrors must be called before this function
    /// is called
    inline std::vector< gsMatrix<T> >& errors()
    {  return m_errors;   }

    /// Returns the errors of the current mapping in the patch i.
    /// !!! computeErrors must be called before this function
    /// is called
    inline gsMatrix<T>& errors(unsigned i)
    {  return m_errors[i];   }

    /// Returns the errors of the current mapping in the
    /// patch i at the position pos.
    /// !!! computeErrors must be called before this function
    /// is called
    inline T error(index_t patch, index_t pos)
    {   return m_errors[patch](pos, 0);   }

    /// Returns the errors of the current mapping at the
    /// position given by the iterator it.
    /// !!! computeErrors must be called before this function
    /// is called
    inline T error(gsPointIterator<T>& it)
    {  return error(it.patch(), it.pos());  }


    /// Sets the errors of the current mapping in the
    /// patch i at the position pos.
    inline void setError(index_t patch, index_t pos, T err)
    {   m_errors[patch](pos, 0) = err;   }

    /// Sets the errors of the current mapping at the
    /// position given by the iterator it.
    inline void setError(gsPointIterator<T>&it, T err)
    {  setError(it.patch(), it.pos(), err);  }


    /// Resets the errors to zero
    void init_m_error();

    virtual void adaptMultipliers(gsFunctionSet<T>* sol,
                                  gsFunctionSet<T>* proj
                                  /*, gsSparseMatrix<T>& mat*/ ) {  }

    /// Computes the error of the current mapping at each point
    /// of the container.
    void computeErrors(index_t type = 0);

    void computeImage(gsFunctionSet<T>* map,
                      std::vector<gsMatrix<T> >& res);

private:
    /// Computes the error according to the type of norm
    /// that has been specified
    T get_local_error(index_t type, T _err);

    /*
    /// Add some entry to the collocation matrix
    /// (NOT USED ANYMORE)
    void add_colloc_matrix(T value, index_t local_ind,
                           gsPointIterator<T>& it,
                           gsSparseEntries<T>& entries); */

protected:
    using Base::m_size;

    /// A vector containing the error of the current
    /// mapping at each point
    std::vector< gsMatrix<T> > m_errors;

    /// Maximal, minimal and total error of the current
    /// mapping at each point
    T m_maxError;
    T m_minError;
    T m_totalError;

    /// The basis used for the construction of the matrix
    /// and of the RHS
    GenBasis& m_basis;

    /// The current position used for:
    ///  - developping the system in case one of the energy is nonlinear
    ///  - computing the error
    gsFunctionSet<T>* m_current_map;

    /// The path used for exporting files
    std::string m_output;

    /// Export all the points of the container before constructing
    /// the matrix if true
    bool m_export_points;

    /// The collocation matrix associated with the points of the
    /// container and the basis m_basis
    gsSparseMatrix<T> m_colloc_mat;

}; // class gsLeastSquares

}
