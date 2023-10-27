/** @file gsRecipeAssemblerLinElast.h

    @brief assemblers for linear elasticity PDEs.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
**/

#pragma once

#include <gsPde/gsLinearElasticityPde.h>
#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsMapUtils/gsL2GMapper.h>

namespace gismo {

// In most parts copied from gsRecipeAssemblerStokes.
class GISMO_EXPORT gsRecipeAssemblerLinElast : public gsRecipeAssembler
{
protected:
    const gsLinearElasticityPde<real_t>     &m_pde;
    bool                                    m_zeroAverage;
    bool                                    m_zeroAvgSet;

    dirichlet::strategy         m_strategy;

public:
    gsRecipeAssemblerLinElast( const gsLinearElasticityPde<real_t> &pde)
        : gsRecipeAssembler(pde.domain()), m_pde(pde), m_zeroAvgSet(false)
    {
        m_space.resize(1);
        m_space[0]=NULL;

        m_lambda = pde.m_lambda;
        m_mu = pde.m_mu;
        m_YoungsModulus = pde.m_YoungsModulus;
        m_PoissonsRatio = pde.m_PoissonsRatio;
    }

    ~gsRecipeAssemblerLinElast(){}


public:
    void setZeroAverage(bool zero)
    { m_zeroAvgSet=true; m_zeroAverage=zero; }

    bool getZeroAverage() const
    { return m_zeroAverage; }

    const gsLinearElasticityPde<real_t> &pde () const
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
    {
        return m_pde.numRhs();
    }
protected:
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

    void setDirichletStrategy(dirichlet::strategy strategy)
    {
        m_strategy=strategy;
    }

    virtual void collectEliminatedDofs  ();

    real_t m_YoungsModulus;
    real_t m_PoissonsRatio;
    real_t m_mu;
    real_t m_lambda;

protected:
    // recipe providing functions
    virtual gsRecipe<real_t>    getPatchRecipe      (index_t patch);

    virtual gsRecipe<real_t>    getBoundaryRecipe   (patchSide ps);

    virtual gsRecipe<real_t>    getDirichletRecipe  (const boundary_condition<real_t> &bc);

    virtual gsRecipe<real_t>    getNeumannRecipe    (const boundary_condition<real_t> &bc);

    // overrides function from gsRecipeAssembler for over-integration
    //virtual gsIntegrationRule   getPatchIntegration (index_t patch );

};


}

