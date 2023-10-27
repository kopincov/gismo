/** @file gsLeastSquaresLMult.h

    @brief Contains the class gsLeastSquaresLMult extending gsLeastSquares.
    This class contains Lagrange multipliers associated with the points
    to be fitted

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsFittingIterator.h>

namespace gismo
{

/*
  @brief class gsLeastSquaresLMult. Permits to compute the
  Lagrange multipliers associated with the points to be fitted.
*/
template<class GenBasis, class T = real_t>
class gsLeastSquaresLMult : public gsLeastSquares<GenBasis, T>
{
    typedef gsLeastSquares<GenBasis, T> Base;

    using gsPointContainer<T>::m_points;
    using gsPointContainer<T>::m_params;
    using gsPointContainer<T>::m_size;
    using Base::m_colloc_mat;
    using Base::m_basis;

public:
    /// constructors

    /// Empty constructor
    gsLeastSquaresLMult(GenBasis& basis)
    : Base::gsLeastSquares(basis)
    {   initLMult();   }

    /// Copy the points contained in the point container
    gsLeastSquaresLMult(GenBasis& basis,
                        const gsPointContainer<T>& points,
                        std::string& output,
                        bool export_points, T coeff_displ)
    : Base::gsLeastSquares(basis, points, output,
                           export_points)
    {
        m_coeff_displacement = coeff_displ;
        initLMult();
    }

    /// Copy the points. The target and the template domains can be
    /// set in order to resample each time the basis is modified
    /// (see gsPointContainer for more details)
    gsLeastSquaresLMult(GenBasis& basis,
                        std::vector< gsMatrix<T> >& pts,
                        std::vector< gsMatrix<T> >& params,
                        std::string& output, bool export_points,
                        T coeff_displ,
                        gsFunctionSet<T>* _template = NULL,
                        gsFunctionSet<T>* geometry = NULL,
                        int *local_size = NULL)
    : Base::gsLeastSquares(basis, pts, params, output,
                           export_points, _template,
                           geometry, local_size)
    {
        initLMult();
        m_coeff_displacement = coeff_displ;
    }

    /// Copy a _gsLeastSquares
    gsLeastSquaresLMult(const gsLeastSquaresLMult<GenBasis, T>& it)
    : Base::gsLeastSquares(it)
    {
        initLMult();
        m_coeff_displacement = it.m_coeff_displacement;
    }

    /// Destructor
    virtual ~gsLeastSquaresLMult(){  }

    /// Called when the basis is modified.
    /// Here, we reset the Lagrange multipliers
    /// (should be improved)
    bool actualize_degrees_freedomLS()
    {
        bool resample =
            Base::actualize_degrees_freedomLS();
        if(resample)
            initLMult();
        return resample;
    }

    /// TODO: Can definitly be improved.
    /// The new Lagrange multiplier should be
    /// copied from Lagrange multipliers before the refinement
    void initLMult()
    {
        m_LMulti.clear();
        for(int i = 0;i < m_size;i++)
        {
            m_LMulti.push_back(gsMatrix<T>
                               (this->dim_im(), this->size(i)));
            m_LMulti[i].setZero();
        }
    }

    void resetLMult()
    {
        initLMult();
    }

    /// In case the scale of the energies change, we call this
    /// function to rescale the Lagrange multipliers as well
    void decreaseLMult(T val);

    /// Change the Lagrange multipliers
    /// sol: the new solution before projection
    /// proj: the new solution after projection
    void adaptMultipliers(gsFunctionSet<T>* sol,
                          gsFunctionSet<T>* proj
                          /*, gsSparseMatrix<T>& mat*/ );

    /// Adds the Lagrange multipliers to the right hand side
    void add_component_SM(gsMatrix<T>& vect,
                          gsPointIterator<T>& it);

private:
    /// Computes the error according to the type of norm
    /// that has been specified
    T get_local_error(int type, T _err);

protected:
    /// Langrange multipliers associated to each point
    std::vector<gsMatrix<T> > m_LMulti;

    /// Used for Uzawa algorithm.
    /// The Lagrange multipliers are adapted at each
    /// step by taking the difference between the position
    /// before and after the projection times this coefficient
    T m_coeff_displacement;

}; // class gsLeastSquaresLMult

}
