/** @file gsMapFactoryDofMapper.h

    @brief Construct a gsWeightMapper from a gsDofMapper, optionally maps eliminated Dofs to the end.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsMapUtils/gsMapFactory.h>
#include <gsCore/gsDofMapper.h>

namespace gismo
{

class GISMO_EXPORT gsMapFactoryDofMapper : public gsMapFactory
{
public:

    /// Constructor providing a dofMapper and a flag if eliminated dofs should be kept
    gsMapFactoryDofMapper(const gsDofMapper& dofMapper,bool keepEliminated) :
        m_dofMapper(dofMapper),m_keepEliminated(keepEliminated){}

    virtual gsWeightMapper<real_t> * makeMapper() const;

private:
    const gsDofMapper& m_dofMapper;
    const bool m_keepEliminated;
};

}
