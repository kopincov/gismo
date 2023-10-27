/** @file gsPhysicalFuncSet_.cpp

    @brief instantiation of gsFunctionSet

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsTransformedFuncSet.h>

namespace gismo {

CLASS_TEMPLATE_INST gsTransformedFuncSetImpl<real_t, gsTransformPureIdentity<real_t> >;
CLASS_TEMPLATE_INST gsTransformedFuncSetImpl<real_t, gsTransformPureGradConforming<real_t> >;
CLASS_TEMPLATE_INST gsTransformedFuncSetImpl<real_t, gsTransformPureDivConforming<real_t> >;
//CLASS_TEMPLATE_INST gsTransformedFuncSetImpl<real_t, gsTransformPureRestriction<real_t> >;

}
