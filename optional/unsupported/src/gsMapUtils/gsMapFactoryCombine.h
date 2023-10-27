/** @file gsMapFactoryCombine.h

    @brief Constructs a gsWeightMapper by combining those constructed by a vector of gsMapFactory

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/

#pragma once

#include <gsMapUtils/gsMapFactory.h>

namespace gismo
{

class GISMO_EXPORT gsMapFactoryCombine : public gsMapFactory
{
public:

    /// Constructor using a path to input file (full path is expected)
    gsMapFactoryCombine( const std::vector<gsMapFactory*> &children, bool needShifts=true)
        : m_children(children),m_needShifts(needShifts)
    {}

    virtual gsWeightMapper<real_t> * makeMapper() const;
    virtual const std::vector<index_t> &getShifts() {return m_shifts;}
private:
    const std::vector<gsMapFactory*> m_children;
    mutable std::vector<index_t>     m_shifts;
    bool                             m_needShifts;
};

}
