
#include <gsCore/gsConfig.h>
#include <gsCore/gsTemplateTools.h>

#ifdef GISMO_WITH_MPI

#include <gsIETI/gsIETIAssemblerMPI.hpp>

namespace gismo
{
  CLASS_TEMPLATE_INST gsIETI_MPI_Base<real_t>;
  CLASS_TEMPLATE_INST gsIETIAssemblerMPI<real_t>;
}

#endif

