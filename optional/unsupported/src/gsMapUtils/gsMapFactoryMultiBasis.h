/** @file gsMapFactoryMultiBasis.h

    @brief Constructs a gsWeightMapper for a MultiBasis using C^0 gluing.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/

#pragma once

#include <gsCore/gsMultiBasis.h>
#include <gsMapUtils/gsMapFactory.h>

namespace gismo
{

class GISMO_EXPORT gsMapFactoryMultiBasis : public gsMapFactory
{
public:

    /// Constructor using a path to input file (full path is expected)
    gsMapFactoryMultiBasis( const gsMultiBasis<real_t> &basis , bool conforming=true) :
        m_basis(basis), m_conforming(conforming)
    {}

    virtual gsWeightMapper<real_t> * makeMapper() const;

private:
    const gsMultiBasis<real_t> &m_basis;
    bool                        m_conforming;
};

}
