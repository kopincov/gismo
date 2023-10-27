#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsBasisEvaluator.h>
#include <gsCore/gsBasisEvaluator.hpp>

#define T real_t

namespace gismo
{
// Quote from the C++ Standard:
// 'An explicit instantiation that names a class template
// specialization is also an explicit instantiation of the same kind
// (declaration or deﬁnition) of each of its members (not including
// members inherited from base classes) that has not been previously
// explicitly specialized in the translation unit containing the
// explicit instantiation, except as described below.'

TEMPLATE_INST
gsBasisEvaluator<T> * makeBasisEvaluator ( const gsBasis<T> &basis, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans);
TEMPLATE_INST
gsBasisEvaluator<T> * makeBasisEvaluator ( const std::vector<gsBasis<T> *> &basis, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans);
TEMPLATE_INST
gsBasisEvaluator<T> * makeBasisEvaluator ( const gsBasis<T> &basis, index_t shift, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans );
TEMPLATE_INST
gsBasisEvaluator<T> * makeBasisEvaluator ( const std::vector<gsBasis<T> *> &basis, const std::vector<index_t>& shifts, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans );
}

#undef T
