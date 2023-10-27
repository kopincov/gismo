
#include <gsCore/gsTemplateTools.h>

#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETIAssembler.hpp>

namespace gismo
{

  CLASS_TEMPLATE_INST gsIETIAssembler<real_t>;

  TEMPLATE_INST  void gsIETIAssembler<real_t>::solveKii<true>(unsigned , gsMatrix<real_t>& , gsMatrix<real_t>& ) const;
  TEMPLATE_INST  void gsIETIAssembler<real_t>::solveKii<false>(unsigned , gsMatrix<real_t>& , gsMatrix<real_t>& ) const;

  TEMPLATE_INST  void gsIETIAssembler<real_t>::solveKC<true, true>(unsigned , gsMatrix<real_t>& , gsMatrix<real_t>& ) const;
  TEMPLATE_INST  void gsIETIAssembler<real_t>::solveKC<false, true>(unsigned , gsMatrix<real_t>& , gsMatrix<real_t>& ) const;
  TEMPLATE_INST  void gsIETIAssembler<real_t>::solveKC<true, false>(unsigned , gsMatrix<real_t>& , gsMatrix<real_t>& ) const;
  TEMPLATE_INST  void gsIETIAssembler<real_t>::solveKC<false, false>(unsigned , gsMatrix<real_t>& , gsMatrix<real_t>& ) const;
}
