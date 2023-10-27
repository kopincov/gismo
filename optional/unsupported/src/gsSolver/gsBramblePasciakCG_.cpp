/** @file gsBramblePascialCG_.cpp

    @brief General Conjugate gradient solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/
#include <gsSolver/gsBramblePasciakCG.h>
#include <gsSolver/gsBramblePasciakCG.hpp>
#include <gsCore/gsTemplateTools.h>

namespace gismo
{

CLASS_TEMPLATE_INST gsBramblePasciakCG<gsBPCG_Types::BP,real_t>;
CLASS_TEMPLATE_INST gsBramblePasciakCG<gsBPCG_Types::BP_Schur,real_t>;
CLASS_TEMPLATE_INST gsBramblePasciakCG<gsBPCG_Types::SchoeberlZulehner,real_t>;


} // end namespace gismo


