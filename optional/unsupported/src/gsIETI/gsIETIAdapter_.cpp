/** gsIETIAdapter_.cpp

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s):
    Created on:  2018-06-05
*/


//#include  <gsIETI/gsIETIAdapter.h>
#include  <gsIETI/gsIETIAdapter.hpp>

#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>
#endif

namespace gismo {

CLASS_TEMPLATE_INST IETIAdapter<real_t>;


#ifdef GISMO_WITH_MPI

CLASS_TEMPLATE_INST IETIAdapterMPI<real_t>;

#endif
} // namespace gismo
