#include <gsCore/gsTemplateTools.h>

#include <gsBem/gsBemLaplace.h>
#include <gsBem/gsBemLaplace.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsBemLaplace<real_t>;
//CLASS_TEMPLATE_INST gsGreenFunction2d<real_t>;

}; // namespace gismo
