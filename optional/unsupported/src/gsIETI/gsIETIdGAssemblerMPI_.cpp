

#include <gsCore/gsConfig.h>
#include <gsCore/gsTemplateTools.h>
#ifdef GISMO_WITH_MPI

#include <gsIETI/gsIETIdGAssemblerMPI.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsIETIdG_MPI_Base<real_t>;
  CLASS_TEMPLATE_INST gsIETIdGAssemblerMPI<real_t>;
}

#endif

