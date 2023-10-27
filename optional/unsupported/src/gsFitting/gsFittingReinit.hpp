/** @file gsFittingLocalReinit.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingReinit.h>
#include <gsFitting/gsFittingIterative.hpp>

namespace gismo {

template <short_t d, class GenBasis, class T>
bool gsFittingReinit<d, GenBasis, T>
::iterAdaptFitting(gsFittingIter<d, GenBasis, T>& fitting,
                   int iter)
{
    if(iter == 0 && fitting.normDisplacement() > m_threshold)
    {
        gsFunctionSet<T>* new_init = NULL;
        if(m_fitting_param.ind == 0)
        {
            /// We compute a first initialization by solving
            /// a linear problem

            /*   new_init = fitPoints<d, GenBasis, 2, T>
                (fitting.leastSquares().points(),
                 fitting.leastSquares().params(),
                 m_fitting_param, NULL, NULL);*/
        }
        else
        {
            gsInfo << std::endl <<
                "----------- Recompute the initialization --------------"
                   << std::endl << std::endl;
            std::vector<gsMatrix<T> > pts;
            std::vector<gsMatrix<T> > param;
            GISMO_ENSURE(! fitting.leastSquares().isEmpty(),
                         "Only implemented for least squares method for the moment");
            interpolatePoints(fitting.leastSquares(),
                              fitting.initialPos(),
                              pts, param);
            m_fitting_param.p_initialization = fitting.initialPos();
            fitPoints<d, GenBasis, 2, T>
                (pts, param, fitting.basis(),
                 m_fitting_param, NULL, NULL);
            gsInfo << std::endl <<
                "----------- End of the computation of the initialization --------------"
                   << std::endl << std::endl;
        }
        GISMO_ENSURE(new_init != NULL,
                     "Error during the computation of the new initialization");
        fitting.reinitialize(new_init);
        return true;
    }
    return false;
}



template <short_t d, class GenBasis, class T>
void gsFittingReinit<d, GenBasis, T>
::interpolatePoints(gsLeastSquares<GenBasis, T>& target,
                    gsFunctionSet<T>* init,
                    std::vector<gsMatrix<T> >& pts,
                    std::vector<gsMatrix<T> >& param)
{
    param = target.params();
    target.computeImage(init, pts);

    gsPointIterator<T> it(target);

    for(;it.good();it.next())
    {
        typename gsMatrix<T>::Column curr_target
            = it.currPoint();
        typename gsMatrix<T>::Column curr_id
            = it.getPointUsingIterator(pts);
        curr_id = m_prop * curr_id
            + (1. - m_prop) * curr_target;
    }
}


};// namespace gismo
