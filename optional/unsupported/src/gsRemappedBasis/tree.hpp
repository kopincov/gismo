/** @file tree.hpp

    @brief This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    Created on:  2015-01-09
*/


#include  <gsRemappedBasis/tree.h>


namespace gismo {

template <typename Data>
std::ostream& operator<< (std::ostream &out, const BinaryTree<Data> &tree)
{
    out<<"\nRoot "<<tree.m_root;
    out<<"\nFree "<<tree.m_freeList;
    out<<"\n\n";

    for(size_t iter=0; iter<tree.m_node.size(); ++iter)
    {
        out<<"NodeId :"<<iter<<"\n";
        tree.m_node[iter].print(out)<<"\n";
    }
    return out;
}

} // namespace gismo
