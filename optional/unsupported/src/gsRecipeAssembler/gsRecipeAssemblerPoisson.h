/** @file gsRecipeAssemblerPoisson.h

    @brief assemblers for second order elliptic problems

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsPde/gsPoissonPde.h>
#include <gsMapUtils/gsL2GMapper.h>

namespace gismo {


class GISMO_EXPORT gsRecipeAssemblerPoisson : public gsRecipeAssembler
{
protected:
    const gsPoissonPde<real_t> &m_pde;
    bool                        m_zeroAverage;
    bool                        m_zeroAvgSet;
    dirichlet::strategy         m_strategy;
    std::vector<boundary_condition<real_t> > m_boundaryCond;
public:
    gsRecipeAssemblerPoisson( const gsPoissonPde<real_t> &pde);

public:

    void setZeroAverage(bool zero)
    {
        m_zeroAvgSet=true;
        m_zeroAverage=zero;
    }
    bool getZeroAverage() const
    {
        return m_zeroAverage;
    }

    void setDirichletStrategy(dirichlet::strategy strategy)
    {
        m_strategy=strategy;
    }

    // functions specifying sizes accordingly to zero average
    // and number of rhs
    virtual index_t getSysSize() const
    {
        return m_zeroAverage ? gsRecipeAssembler::getSysSize()+1 : gsRecipeAssembler::getSysSize();
    }

    virtual index_t getRhsDim() const
    {
        return m_pde.numRhs();
    }
protected:
    virtual index_t getFreeLimit() const
    {
        return gsRecipeAssembler::getSysSize();
    }
public:
    virtual void init()
    {
        collectEliminatedDofs  ();
        if (!m_zeroAvgSet)
        {
            if ( m_strategy != dirichlet::elimination
               || m_eliminatedDofs[0].size() || m_eliminatedTarget.size())
                m_zeroAverage=false;
            else
                m_zeroAverage=true;
        }
        gsRecipeAssembler::init();
    }


    virtual void collectEliminatedDofs  ();

protected:
    // recipe providing functions
    virtual gsRecipe<real_t>    getPatchRecipe     (index_t patch);

    virtual gsRecipe<real_t>    getBoundaryRecipe  (patchSide ps);

    virtual gsRecipe<real_t>    getDirichletRecipe (const boundary_condition<real_t> &bc);

    virtual gsRecipe<real_t>    getNitscheRecipe   (const boundary_condition<real_t> &bc);

    virtual gsRecipe<real_t>    getNeumannRecipe   (const boundary_condition<real_t> &bc);

    virtual gsRecipe<real_t>    getRobinRecipe     (const boundary_condition<real_t> &bc);

private:
    // utility functions
    real_t nitche (index_t patch);

    real_t hparam (index_t patch);
};


}

