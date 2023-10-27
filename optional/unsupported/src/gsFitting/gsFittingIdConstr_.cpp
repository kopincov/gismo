
#include <gsCore/gsTemplateTools.h>

#include  <gsFitting/gsFittingIdConstr.h>

namespace gismo {

CLASS_TEMPLATE_INST gsFittingIdConstr
<2, gsMultiBasis<real_t>, real_t>;
CLASS_TEMPLATE_INST gsFittingIdConstr
<2, gsBasis<real_t>, real_t>;
CLASS_TEMPLATE_INST gsFittingIdConstr
<3, gsMultiBasis<real_t>, real_t>;
CLASS_TEMPLATE_INST gsFittingIdConstr
<3, gsBasis<real_t>, real_t>;
CLASS_TEMPLATE_INST gsFittingIdConstr
<4, gsMultiBasis<real_t>, real_t>;
CLASS_TEMPLATE_INST gsFittingIdConstr
<4, gsBasis<real_t>, real_t>;

} // namespace gismo
