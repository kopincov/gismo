/** @file gsRecipeAssemblerPervertedStokes.h

    @brief assemblers for a Perverted StokesEquation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, J. Sogn
**/

#pragma once

#include <gsPde/gsPervertedStokesPde.h>
#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsMapUtils/gsL2GMapper.h>

namespace gismo {


enum method
{
    RaviartThomas,
    TaylorHood,
};

enum {
    velocityState = 0,
    pressureState = 1,
    velocityMultiplier = 2,
    pressureMultiplier = 3
};

GISMO_EXPORT std::vector<gsPhysicalSpace*> constructTHSpaces (
        const std::vector<gsBasis<real_t> *>       &bases,
        const gsMultiPatch<real_t>                 &geo,
        //    const gsMapFactory                         &factory,
        std::vector<std::vector<gsBasis<>*> >      *outBasis=NULL
        );

GISMO_EXPORT std::vector<gsPhysicalSpace*> constructPervertedTHSpaces (
    const std::vector<gsBasis<real_t> *>       &bases,
    const gsMultiPatch<real_t>                 &geo,
    //    const gsMapFactory                         &factory,
    std::vector<std::vector<gsBasis<>*> >      *outBasis,
    index_t                              numRefine = 2,
    index_t                              increasePolydegree = 1
    );


class GISMO_EXPORT gsRecipeAssemblerPervertedStokes : public gsRecipeAssembler
{
protected:
    const gsPervertedStokesPde<real_t>              &m_pde;
    //bool                                    m_zeroAverage;
    //bool                                    m_zeroAvgSet;
public:

    gsRecipeAssemblerPervertedStokes(const gsPervertedStokesPde<real_t> &pde);
public:
    const gsPervertedStokesPde<real_t> &pde () const
    {return m_pde;}
    virtual index_t getRhsDim() const
    { return m_pde.numRhs(); }
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

    virtual gsRecipe<real_t>    getDirichletRecipe  (const boundary_condition<real_t> &bc);

    virtual gsRecipe<real_t>    getNeumannRecipe    (const boundary_condition<real_t> &bc);

    virtual gsIntegrationRule   getPatchIntegration (index_t patch );

};


}

