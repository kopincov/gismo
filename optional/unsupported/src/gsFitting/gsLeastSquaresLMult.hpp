/** @file gsLeastSquaresLMult.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/

#pragma once

#include <gsFitting/gsLeastSquares.hpp>
#include <gsFitting/gsPointContainer.hpp>
#include <gsFitting/gsLeastSquaresLMult.h>


namespace gismo
{

template<class GenBasis, class T> void
gsLeastSquaresLMult<GenBasis, T>::
adaptMultipliers(gsFunctionSet<T>* sol,
                 gsFunctionSet<T>* proj
                 /*, gsSparseMatrix<T>& mat*/ )
{
    bool use_mat = false;
    std::vector<gsMatrix<T> > image_proj, image_sol;
    std::vector<gsMatrix<T> > image_func_diff;
    gsPointIterator<T> it(*this);

    if(! use_mat)
    {
        gsPointIterator<T> it(*this);
        gsPointIterator<T> it2(*this);
        this->computeImage(sol, image_sol);
        this->computeImage(proj, image_proj);
    }
    /*  else
    {
        gsMatrix<T> coeff_diff, coeff_sol, coeff_proj;
        gsMatrix<T> coeff_im;
        geometryToCoeffGen<T>(m_basis, *sol, coeff_sol,
                              m_local_global);
        geometryToCoeffGen<T>(m_basis, *proj, coeff_proj,
                              m_local_global);
        coeff_diff = coeff_proj - coeff_sol;
        coeff_im = mat * coeff_diff;

        gsFunctionSet<T>* func_diff =
            makeGeometryGen(m_basis, coeff_im, m_local_global);

        this->computeImage(func_diff, image_func_diff);
        }*/


    for(;it.good();it.next())
    {
        typename gsMatrix<T>::Column curr_mult
            = it.getPointUsingIterator(m_LMulti);
        if(use_mat)
        {
            const typename gsMatrix<T>::Column curr_diff
                = it.getPointUsingIterator(image_func_diff);

            curr_mult -= m_coeff_displacement * curr_diff;
        }
        else
        {
            const typename gsMatrix<T>::Column curr_proj
                = it.getPointUsingIterator(image_proj);
            const typename gsMatrix<T>::Column curr_sol
                = it.getPointUsingIterator(image_sol);
            curr_mult -= m_coeff_displacement * ( curr_proj - curr_sol );
        }
    }
}


template<class GenBasis, class T> void
gsLeastSquaresLMult<GenBasis, T>::
decreaseLMult(T coeff)
{
    gsPointIterator<T> it(*this);

    for(;it.good();it.next())
    {
        typename gsMatrix<T>::Column curr_mult
            = it.getPointUsingIterator(m_LMulti);
        curr_mult *= coeff;
    }
}

template<class GenBasis, class T> void
gsLeastSquaresLMult<GenBasis, T>::
add_component_SM(gsMatrix<T>& vect, gsPointIterator<T>& it)
{
    const typename gsMatrix<T>::Column curr_mult
        = it.getPointUsingIterator(m_LMulti);
    vect -= curr_mult;
}


} // namespace gismo
