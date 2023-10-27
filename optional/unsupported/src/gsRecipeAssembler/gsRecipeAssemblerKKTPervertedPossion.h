/** @file gsRecipeAssemblerKKTPervertedPoisson.h

    @brief assemblers for a Perverted StokesEquation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, J. Sogn
**/

#pragma once

#include <gsPde/gsPoissonPde.h>
#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsMapUtils/gsL2GMapper.h>

namespace gismo {



enum {
    controlSpace = 0,
    stateSpace = 2,
    multiplierSpace = 1
};


class GISMO_EXPORT gsRecipeAssemblerKKTPervertedPoisson : public gsRecipeAssembler
{
protected:
    const gsPoissonPde<real_t>              &m_pde;
    gsBoundaryConditions<real_t>              m_desired;
    const real_t                           m_alpha;
public:

    gsRecipeAssemblerKKTPervertedPoisson(const gsPoissonPde<real_t> &pde, const gsBoundaryConditions<real_t> &desired, real_t m_alpha = 1.0);

public:
    const gsPoissonPde<real_t> &pde () const
    {return m_pde;}
public:
    real_t getAplha() const {return m_alpha;}


public:
    virtual void init()
    {
        collectEliminatedDofs  ();
        gsRecipeAssembler::init();
    }

    virtual void collectEliminatedDofs  ();

protected:
    // recipe providing functions
    virtual gsRecipe<real_t>    getPatchRecipe      (index_t patch);

    virtual gsRecipe<real_t>    getBoundaryRecipe   (patchSide ps);

    virtual gsIntegrationRule   getPatchIntegration (index_t patch );
};


}

