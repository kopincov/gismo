/** @file gsMapFactoryCombine.cpp

    @brief Provides implementation of gsMapFactoryCombine class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/

#include <gsMapUtils/gsMapFactoryCombine.h>
#include <gsMapUtils/gsMapperUtils.h>

namespace gismo
{

gsWeightMapper<real_t> * gsMapFactoryCombine::makeMapper() const
{
    std::vector<gsWeightMapper<real_t>*> maps(m_children.size());
    for(size_t i=0;i<maps.size();++i)
        maps[i]=m_children[i]->makeMapper();
    gsWeightMapper<real_t>* mapper = combineMappers(maps, m_shifts, m_needShifts);
    return mapper;
}

}
