/** gsParallelCG.cpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsConfig.h>
#ifdef GISMO_WITH_MPI
//#include  <gsIETI/gsParallelCG.h>
#include <gsIETI/gsParallelCG.hpp>

namespace gismo {


CLASS_TEMPLATE_INST gsParallelCG<real_t>;


} // namespace gismo
#endif
