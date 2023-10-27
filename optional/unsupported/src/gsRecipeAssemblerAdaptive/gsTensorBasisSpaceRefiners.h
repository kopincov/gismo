/** @file gsTensorBasisSpaceRefiners.h

    @brief space refiners based on gsCompositeHBasis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/


#pragma once

#include <gsAssembler/gsAdaptiveRefUtils.h>

#include <gsRecipeAssemblerAdaptive/gsSpaceRefiner.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsCore/gsMultiBasis.h>

namespace gismo {

class gsMultiBasisRefiner : public gsSpaceRefiner
{
protected:
    const   gsMultiPatch<>     &m_domain;
    mutable gsMultiBasis<>      m_basis;
public:
    gsMultiBasisRefiner(
            const gsMultiPatch<>           &domain,
            const gsMultiBasis<>           &basis
            )
        : m_domain(domain), m_basis(basis)
    {}

    virtual std::vector<gsPhysicalSpace*> getSpaces()const
    {
        std::vector<gsPhysicalSpace*> result;
        result.push_back(new gsPhysicalSpaceScalar (m_basis,m_domain,INVERSE_COMPOSITION));
        return result;
    }

    virtual void updateSpaces (const gsMatrix<real_t>& ) // markedCells)
    {
        m_basis.uniformRefine();
    }

    const gsMultiBasis<> & getBasis() {return m_basis;}
};






}

