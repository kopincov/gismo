/** @file gsBoundariesSplit.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingLinSyst.h>
#include <gsFitting/gsFittingQuadrature.h>
#include <gsFitting/gsFittingIntegrandNL.h>

namespace gismo
{

/**
   We suppose that all the curves split are on the boundary,
   i.e., they can after the splitting be strongly imposed
*/
template<class T>
class gsFittingBCSplit
{

protected:
    struct type{
        T init;
        T end;
        int patch;
        bool strong;
        int side;
        bool sign_orien;

        type(T _init, T _end, int _patch, bool _strong,
             int _side, bool _sign_orien)
        {
            init = _init;
            end = _end;
            patch = _patch;
            strong = _strong;
            side = _side;
            sign_orien = _sign_orien;
        }
    };

    gsFittingBC<T>& m_origin;
    std::vector<type> m_splitting;
    std::vector< gsFittingBC<T> > m_new_bc;


public:
    gsFittingBCSplit(gsFittingBC<T>& origin,
                     int side, bool sign_orien,
                     bool strong = true)
    : m_origin(origin),
      m_splitting(type(0., 1., 0, strong, side, sign_orien)),
      m_new_bc()   {}
    gsFittingBCSplit(gsFittingBC<T>& origin,
                     const std::vector<type>& splitting)
    : m_origin(origin), m_splitting(splitting), m_new_bc()  {}

    /// Split the bc for the purpose of constructing a parametrization
    /// of the domain delimited by this bc
    template<class GenBasis>
    void split_bc(GenBasis& basis);

};  /// gsFittingBCSplit

}
