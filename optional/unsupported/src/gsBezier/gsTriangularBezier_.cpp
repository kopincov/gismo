#include <gsCore/gsTemplateTools.h>

#include <gsBezier/gsTriangularBezier.h>
#include <gsBezier/gsTriangularBezierXML.h>

namespace gismo
{


    CLASS_TEMPLATE_INST internal::gsXml< gsTriangularBezier<2,real_t> >;
}
