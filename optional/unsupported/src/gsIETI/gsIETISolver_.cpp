

#include <gsCore/gsTemplateTools.h>

#include <gsIETI/gsIETISolver.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsIETIJumpOperator<real_t>;
  CLASS_TEMPLATE_INST gsIETISolver<real_t>;
  CLASS_TEMPLATE_INST gsInexactIETISolver<real_t>;
}

