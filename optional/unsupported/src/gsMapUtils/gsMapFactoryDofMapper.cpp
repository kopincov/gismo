/** @file gsMapFactoryDofMapper.cpp

    @brief Provides implementation of gsMapFactoryDofMapper class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsMapUtils/gsMapFactoryDofMapper.h>

namespace gismo
{

gsWeightMapper<real_t> * gsMapFactoryDofMapper::makeMapper() const
{
    gsWeightMapper<real_t> * mapper;
    index_t m_globals=m_dofMapper.size();
    index_t m_locals = m_dofMapper.mapSize();
    mapper = new gsWeightMapper<real_t>(m_locals,m_globals);
    index_t local=0;
    for(size_t patch=0;patch<m_dofMapper.numPatches();++patch)
    {
        const size_t patchSize= (patch!=m_dofMapper.numPatches()-1) ?
                    m_dofMapper.offset(patch+1)-m_dofMapper.offset(patch) :
                    m_locals-m_dofMapper.offset(m_dofMapper.numPatches()-1);
        for(size_t i = 0;i<patchSize;++i)
        {
            const index_t global = m_dofMapper.index(i,patch);
            if( m_keepEliminated || !m_dofMapper.is_boundary_index(global) )
                mapper->setEntry(local,global,1);
            local++;
        }
    }
    return mapper;
}

}
