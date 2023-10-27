/** @file gsRecipeAssemblerStokes.h

    @brief assemblers for StokesEquation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsPde/gsStokesPde.h>
#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsMapUtils/gsL2GMapper.h>

namespace gismo {


enum method
{
    RaviartThomas,
    TaylorHood,
    SubGrid
};

enum {
    velocity = 0,
    pressure = 1
};

GISMO_EXPORT std::vector<gsPhysicalSpace*> constructTHSpaces (
        const std::vector<gsBasis<real_t> *>       &bases,
        const gsMultiPatch<real_t>                 &geo,
        //    const gsMapFactory                         &factory,
        std::vector<std::vector<gsBasis<>*> >      *outBasis=NULL
        );


class GISMO_EXPORT gsRecipeAssemblerStokes : public gsRecipeAssembler
{
protected:
    const gsStokesPde<real_t>              &m_pde;
    bool                                    m_zeroAverage;
    bool                                    m_zeroAvgSet;
public:

    gsRecipeAssemblerStokes( const gsStokesPde<real_t> &pde);
public:
    void setZeroAverage(bool zero)
    { m_zeroAvgSet=true; m_zeroAverage=zero; }
    bool getZeroAverage() const
    { return m_zeroAverage; }
    const gsStokesPde<real_t> &pde () const
    {return m_pde;}


    // functions specifying sizes accordingly to zero average
    // and number of rhs
    virtual index_t getSysSize() const
    {
        return m_zeroAverage ? gsRecipeAssembler::getSysSize()+1 : gsRecipeAssembler::getSysSize();
    }
    virtual index_t getFreeLimit() const
    {
        return gsRecipeAssembler::getSysSize();
    }
    virtual index_t getRhsDim() const
    { return m_pde.numRhs(); }
public:
    virtual void init()
    {
        collectEliminatedDofs  ();
        if (!m_zeroAvgSet)
        {
            if (m_eliminatedDofs[0].size() || m_eliminatedTarget.size())
                m_zeroAverage=false;
            else
                m_zeroAverage=true;
        }

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

