/** @file gsRecipeAssemblerQMOpt.h

    @brief provides an assembler for optimizing quality measures,
    inherited from the gsRecipeAssembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssembler/gsQMOptOperators.h>
#include <gsMapUtils/gsL2GMapper.h>

namespace gismo{

//========================================= ASSEMBLER ==========================================//

class GISMO_EXPORT gsRecipeAssemblerQMOpt2D : public gsRecipeAssembler
{
public:
    const gsQualityMeasureWeights          &m_weights;
          gsMatrix<real_t>                 *m_qualityMeasureVals;
          gsMatrix<real_t>                 *m_areas;
    const bool                              m_assembleSystemMat;
    const bool                              m_assembleRhs;
    const bool                              m_assembleQualityMeasureVals;
    const bool                              m_assembleAreas;
    const bool                              m_assembleQualityDeriv;
    const bool                              m_assembleFitting;
          bool                              m_fixBoundaries;

    gsRecipeAssemblerQMOpt2D(const gsMultiPatch<real_t>                    &geo,
                                 const std::vector<gsBasis<real_t>*>       &spaces,
                                 const gsWeightMapper<real_t>              &mapper,
                                 const gsQualityMeasureWeights             &weights,

                                 const bool                                 assembleSystemMat,
                                 const bool                                 assembleRhs,
                                 const bool                                 assembleQualityMeasureVals,
                                 const bool                                 assembleAreas,
                                 const bool                                 assembleQualityDeriv,
                                 const bool                                 assembleFitting
                                 );

    gsRecipeAssemblerQMOpt2D(const gsMultiPatch<real_t>                    &geo,
                                 const std::vector<gsPhysicalSpace*>        phySpaces,
                                 const gsQualityMeasureWeights             &weights,

                                 const bool                                 assembleSystemMat,
                                 const bool                                 assembleRhs,
                                 const bool                                 assembleQualityMeasureVals,
                                 const bool                                 assembleAreas,
                                 const bool                                 assembleQualityDeriv,
                                 const bool                                 assembleFitting
                                 );

    ~gsRecipeAssemblerQMOpt2D()
    { freeAll(m_space); }

    virtual void init();

    virtual void fixBoundaries  ();

    virtual void initSystemMatrices();

    gsRecipe<real_t>    getPatchRecipe     (index_t patch);

    gsRecipe<real_t>    getBoundaryRecipe  (patchSide ) // ps)
    { return gsRecipe<real_t>(); }

    bool getQualityMeasureVals(gsMatrix<real_t> & qualityMeasureVals) const;

    bool getAreas(gsMatrix<real_t> & areas) const;

    bool getQualityDerivs(gsMatrix<real_t> & qualityDerivs) const;

protected:

    /**
     * @brief postProcess
     *        called after the assembling is completed
     */
    void postProcess();

    void assembleFittingMatrices();

};
}

