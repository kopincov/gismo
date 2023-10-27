#include <gsCore/gsTemplateTools.h>

#include <gsBoxSplines/gsBoxSplineBasis.h>
#include <gsBoxSplines/gsBoxSplineBasis.hpp>


namespace gismo
{

CLASS_TEMPLATE_INST gsBoxSplineBasis<2,real_t>;
    CLASS_TEMPLATE_INST gsBoxSplineBasis<3,real_t>;

}
