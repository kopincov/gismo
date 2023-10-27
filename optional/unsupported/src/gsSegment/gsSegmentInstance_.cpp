
#include <gsCore/gsTemplateTools.h>

#include <gsSegment/gsMVInterpolation.h>
#include <gsSegment/gsMVInterpolation.hpp>
#include <gsSegment/gsVolumeSegment.h>
#include <gsSegment/gsVolumeSegment.hpp>

namespace gismo
{

#define T real_t
#define uZ unsigned
#define Z int


CLASS_TEMPLATE_INST gsMVInterpolation<T>;
CLASS_TEMPLATE_INST gsMVInterpolationComponent<T>;
CLASS_TEMPLATE_INST gsVolumeSegment<T>;


#undef T
#undef uZ
#undef Z

}
