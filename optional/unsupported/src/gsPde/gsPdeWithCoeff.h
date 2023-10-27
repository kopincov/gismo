
/** @file gsPdeWithCoeff.h

    @brief Provides an interfaces to classes of type gsPde that is used by IETI

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/


#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{

template<class T> struct gsPdeWithCoeff { virtual T getCoeffForIETI(unsigned np) const = 0; };

}
