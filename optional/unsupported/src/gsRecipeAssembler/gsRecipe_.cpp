/** gsRecipe_.cpp.in

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s):
    Created on:  2014-10-28
*/

#include <gsCore/gsTemplateTools.h>
#include <gsRecipeAssembler/gsRecipe.hpp>

namespace gismo {
    CLASS_TEMPLATE_INST gsRecipeIngredient<real_t>;
    CLASS_TEMPLATE_INST gsRecipe<real_t>;

} // namespace gismo
