/** @file gsRecipeAssemblerBiharmonicSogn.h

    @brief assembler for 4'th order problems, assembles: (\Delta u, \Delta v)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
**/

#pragma once

#include <gsCore/gsConstantFunction.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsPde/gsBiharmonicPde.h>
#include <gsMapUtils/gsL2GMapper.h>

namespace gismo {


class GISMO_EXPORT gsRecipeAssemblerBiharmonicSogn : public gsRecipeAssembler
{
protected:
    const gsBiharmonicPde<real_t> &m_pde;
    bool                        m_zeroAverage;
    bool                        m_zeroAvgSet;
    dirichlet::strategy         m_strategy;
    std::vector<boundary_condition<real_t> > m_boundaryCond;
public:
    gsRecipeAssemblerBiharmonicSogn( const gsBiharmonicPde<real_t> &pde);

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

    virtual gsRecipeIngredient<real_t>  getDirichletFirstSysIngredient(const boundary_condition<real_t> &bc);
    virtual gsRecipeIngredient<real_t>  getDirichletFirstRhsIngredient(const boundary_condition<real_t> &bc);

    //Not implemented
    //virtual gsRecipeIngredient<real_t> getDirichletSecondRecipe(const boundary_condition<real_t> &bc);

    //virtual gsRecipeIngredient<real_t>    getNitscheRecipe  (const boundary_condition<real_t> &bc);

    virtual gsRecipeIngredient<real_t>  getNeumannFirstIngredient  (const boundary_condition<real_t> &bc);
    virtual gsRecipeIngredient<real_t>  getNeumannSecondIngredient  (const boundary_condition<real_t> &bc);


};


}

