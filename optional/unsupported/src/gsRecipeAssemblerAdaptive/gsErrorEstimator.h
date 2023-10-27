/** @file gsErrorEstimator.h

    @brief DESCRIPTION

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/



#pragma once

#include <gsAssembler/gsAdaptiveRefUtils.h>
#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>

namespace gismo {

class gsErrorEstimator
{
protected:
    gsMatrix<real_t> m_total;
    gsMatrix<real_t> m_local;
    bool m_upToDateFlag;
public:
    gsErrorEstimator() : m_upToDateFlag(false) { }
    virtual ~gsErrorEstimator() {}

    void estimate(
        const std::vector<gsPhysicalSpace*>  space,
        const std::vector<gsMatrix<real_t> > &coefs
        )
    {
        computeEstimate(space,coefs);
        m_upToDateFlag=true;
    }

    const gsMatrix<real_t> & getTotalErrorEstimate() const
    {
        if(!m_upToDateFlag)
            GISMO_ERROR("Estimate requested before the call to estimate.");
        return m_total;
    }

    const gsMatrix<real_t> & getLocalErrorEstimate() const
    {
        if(!m_upToDateFlag)
            GISMO_ERROR("Estimate requested before the call to estimate.");
        return m_local;
    }

    virtual void reset()
    {
        m_upToDateFlag=false;
    }

    virtual void computeEstimate(
        const std::vector<gsPhysicalSpace*>  &space,
        const std::vector<gsMatrix<real_t> > &coefs
        )   =0;
};

class gsErrorEstimatorPerCellExact : public gsErrorEstimator , public gsRecipeDistance
{
public:
    gsErrorEstimatorPerCellExact(const gsMultiPatch<real_t> &domain,const std::vector<gsFunction<real_t>*> functions)
        : gsRecipeDistance(domain,functions)
    {
        m_norms.clear();
        gsMatrix<> norms(1,2);
        norms<<0,2;
        for (size_t sp=0; sp<m_ref.size(); ++sp)
            m_norms.push_back(norms); // default to L2 norm
        m_storePerCell=true;
        m_storePerFunc=false;
    }

    virtual void computeEstimate(
        const std::vector<gsPhysicalSpace*>  &space,
        const std::vector<gsMatrix<real_t> > &coefs
        )
    {
        initInternalData(space, coefs);
        gsRecipeDistance::assemble();
        m_local = gsRecipeDistance::getPerCell();
        m_total = m_NormM;
    }
};

}
