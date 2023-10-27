
#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsConfig.h>
#ifdef GISMO_WITH_MPI

#include <gsIETI/gsIETISolverMPI.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsIETIJumpOperatorMPI<real_t>;
  CLASS_TEMPLATE_INST gsIETISolverMPI<real_t>;
}

#endif
