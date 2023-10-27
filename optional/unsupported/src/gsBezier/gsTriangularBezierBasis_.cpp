#include <gsCore/gsTemplateTools.h>

#include <gsBezier/gsTriangularBezierBasis.h>
#include <gsBezier/gsTriangularBezierBasisXML.h>

namespace gismo
{
   //CLASS_TEMPLATE_INST gsXml< gsTriangularBezier<2,real_t> >;

   CLASS_TEMPLATE_INST internal::gsXml< gsTriangularBezierBasis<2,real_t> >;
}
