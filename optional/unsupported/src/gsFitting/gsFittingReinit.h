/** @file gsFittingReinit.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingIterative.h>

namespace gismo {

/**
   \brief class gsFittingReinit: in case the norm of the displacement of the first iteration is too important, we compute a new initialization of the nonlinear process. This is obtained by computing a simpler parametrization, i.e., with less deformation (difference between the initialization and the expected parametrization).
   To this purpose, we interpolate the target points with the current
   initialization. Taking the resulting points as new target points, we obtain the expected simpler problem. Supposing that one find a solution of this problem, we use it as a new initialization for the current process.

   NOT USED ANYMORE
*/
template <short_t d, class GenBasis, class T>
class gsFittingReinit :
        public gsFittingImproveAfterCV<d, GenBasis, T>
{

public:
    /**
       \brief
       Main constructor
    */
    gsFittingReinit(gsFittingParam<T>& param,
                    T threshold, T prop = 0.5)
    : gsFittingImproveAfterCV<d, GenBasis, T>::gsFittingImproveAfterCV(),
    m_fitting_param(param)
    {
        m_fitting_param.output = param.output + "_reinit";
        //     m_fitting_param.export_iterations = false;
        m_fitting_param.export_iterations = true;
        m_fitting_param.print_messages = true;
        m_fitting_param.use_refinement = false;
        m_prop = prop;
        m_threshold = threshold;
        //    m_fitting_param.export_points = true;

        m_fitting_param.ind = param.ind + 1;
        //   m_fitting_param.max_num_iter_NL = 1;
    }

    bool iterAdaptFitting(gsFittingIter<d, GenBasis, T>& fitting,
                          int iter);

    void interpolatePoints(gsLeastSquares<GenBasis, T>& target,
                           gsFunctionSet<T>* init,
                           std::vector<gsMatrix<T> >& pts,
                           std::vector<gsMatrix<T> >& param);

public:
    /// The coefficient used for the interpolation
    /// of the boundary
    T m_prop;

    /// If the first displacement is above this threshold,
    /// we recompute the initialization
    T m_threshold;

    gsFittingParam<T> m_fitting_param;
};

};// namespace gismo
