/** @file gsMapFactoryMultiBasis.cpp

    @brief Provides implementation of gsMapFactoryMultiBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/

#include <gsMapUtils/gsMapFactoryMultiBasis.h>
#include <gsMapUtils/gsMapFactoryDofMapper.h>

namespace gismo
{

gsWeightMapper<real_t> * gsMapFactoryMultiBasis::makeMapper() const
{
    gsDofMapper dofMap;
    m_basis.getMapper(m_conforming,dofMap);
    gsMapFactoryDofMapper mapper(dofMap,true);
    return mapper.makeMapper();
}

}
