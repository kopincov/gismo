/** @file gsPDEOperators_.cpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsCore/gsTemplateTools.h>
#include <gsRecipeAssembler/gsPDEOperators.hpp>

namespace gismo {

CLASS_TEMPLATE_INST gsBilinearOp<real_t>;
CLASS_TEMPLATE_INST gsLaplaceOp<real_t>;
CLASS_TEMPLATE_INST gsLaplaceLaplaceOp<real_t>;
CLASS_TEMPLATE_INST gsGradGradOp<real_t>;
CLASS_TEMPLATE_INST gsGenericSecondOrderOp<real_t>;
CLASS_TEMPLATE_INST gsL2ScalarOp<real_t>;
CLASS_TEMPLATE_INST gsL2TestOp<real_t>;
CLASS_TEMPLATE_INST gsL2TestVecOp<real_t>;
CLASS_TEMPLATE_INST gsDivergenceOp<real_t>;
CLASS_TEMPLATE_INST gsGradientOp<real_t>;
CLASS_TEMPLATE_INST gsBoundaryL2TestOp<real_t>;
CLASS_TEMPLATE_INST gsBoundaryL2TestVecOp<real_t>;
CLASS_TEMPLATE_INST gsBoundaryL2ScalarOp<real_t>;
CLASS_TEMPLATE_INST gsBoundaryNormalDerValueOp<real_t>;
CLASS_TEMPLATE_INST gsBoundaryNormalDerNormalDerOp<real_t>;
CLASS_TEMPLATE_INST gsBoundaryNormalDerTestOp<real_t>;
CLASS_TEMPLATE_INST gsBoundaryNormalDerTestNormalOp<real_t>;
CLASS_TEMPLATE_INST gsBoundaryNormalDerTestVecOp<real_t>;
CLASS_TEMPLATE_INST gsLinElastOp<real_t>;

} // namespace gismo
