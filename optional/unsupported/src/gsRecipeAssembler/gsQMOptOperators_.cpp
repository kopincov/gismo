/** @file gsQMOptOperators_.cpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
**/

#include <gsCore/gsTemplateTools.h>
#include <gsRecipeAssembler/gsQMOptOperators.hpp>

namespace gismo {

CLASS_TEMPLATE_INST gsQualityMeasureOp<real_t>;
CLASS_TEMPLATE_INST gsQualitySystemMatrixOp<real_t>;
CLASS_TEMPLATE_INST gsQualityRHSOp<real_t>;
CLASS_TEMPLATE_INST gsQualityMeasureValueOp<real_t>;
CLASS_TEMPLATE_INST gsQualityConstraintOp<real_t>;

} // namespace gismo
