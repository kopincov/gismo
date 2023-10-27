

#include <gsCore/gsTemplateTools.h>
#include <gsRemappedBasis/gsSelector.h>
#include <gsRemappedBasis/tree.h>
#include <gsRemappedBasis/tree.hpp>

namespace gismo
{

TEMPLATE_INST std::ostream& operator<< (std::ostream &out, const BinaryTree<NodeData> &tree);

}
