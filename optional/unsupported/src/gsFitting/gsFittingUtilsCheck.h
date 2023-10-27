/** @file gsFittingUtilsCheck.h

    @brief Files containing all the routines permitting to check the regularity of some mapping.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

namespace gismo
{

/*
  Returns true if the determinant of the gradient of func is positive at
  each corner of the domain.
  Here basis permits to determine if func is single or multi-patch.
 */
template<class T>
bool positiveDet(gsBasis<T>& basis,
                 gsFunctionSet<T>& func, T tol = 0.00001);
template<class T>
bool positiveDet(gsMultiBasis<T>& basis,
                 gsFunctionSet<T>& func, T tol = 0.00001);
template<class T>
bool positiveDet(gsMultiPatch<T>& func, T tol = 0.00001);
template<class T>
bool positiveDet(gsGeometry<T>& func, T tol = 0.00001);


} /// namespace gismo
