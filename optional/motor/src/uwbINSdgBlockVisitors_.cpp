
#include <gsCore/gsForwardDeclarations.h>

#include "uwbINSdgBlockVisitors.h"
#include "uwbINSdgBlockVisitors.hpp"

namespace gismo
{

CLASS_TEMPLATE_INST uwbINSdgBlockVisitor<real_t>;
CLASS_TEMPLATE_INST uwbINSdgBlockAVisitor<real_t>;
CLASS_TEMPLATE_INST uwbINSdgBlockBVisitor<real_t>;
CLASS_TEMPLATE_INST uwbINSdgBlockNVisitor<real_t>;
CLASS_TEMPLATE_INST uwbINSdgBlockPVisitor<real_t>;
CLASS_TEMPLATE_INST uwbINSdgBlockApVisitor<real_t>;
CLASS_TEMPLATE_INST uwbINSiFaceSizeVisitor<real_t>;

} // namespace gismo
