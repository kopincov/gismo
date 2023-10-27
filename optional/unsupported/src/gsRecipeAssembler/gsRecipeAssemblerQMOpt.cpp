/** @file gsRecipeAssemblerQMOpt.cpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssembler/gsQMOptOperators.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsRecipeAssembler/gsRecipeAssemblerQMOpt.h>

namespace gismo{

gsRecipeAssemblerQMOpt2D::gsRecipeAssemblerQMOpt2D(const gsMultiPatch<real_t>                &geo,
                                                   const std::vector<gsBasis<real_t>*>       &bases,
                                                   const gsWeightMapper<real_t>              &mapper,
                                                   const gsQualityMeasureWeights             &weights,

                                                   const bool                                 assembleSystemMat,
                                                   const bool                                 assembleRhs,
                                                   const bool                                 assembleQualityMeasureVals,
                                                   const bool                                 assembleAreas,
                                                   const bool                                 assembleQualityDeriv,
                                                   const bool                                 assembleFitting
                                                   )
    : gsRecipeAssembler(geo),m_weights(weights),m_qualityMeasureVals(NULL),m_areas(NULL),
      m_assembleSystemMat(assembleSystemMat),m_assembleRhs(assembleRhs),
      m_assembleQualityMeasureVals(assembleQualityMeasureVals),m_assembleAreas(assembleAreas),
      m_assembleQualityDeriv(assembleQualityDeriv), m_assembleFitting(assembleFitting),
      m_fixBoundaries(false)
{
    gsPhysicalSpaceScalar physicalSpace = gsPhysicalSpaceScalar(bases,geo,NO_TRANSFORMATION,mapper);
    std::vector<gsPhysicalSpaceScalar*> phySpaces;
    phySpaces.push_back(&physicalSpace);
    phySpaces.push_back(&physicalSpace);
    gsPhysicalSpace* physicalSpaceVector = new gsPhysicalSpaceVector(phySpaces);
    std::vector<gsPhysicalSpace* > vectorSpaces;
    vectorSpaces.push_back(physicalSpaceVector);
    setSpace(vectorSpaces);

    if(m_assembleQualityMeasureVals)
        m_qualityMeasureVals=new gsMatrix<real_t>();
    if(m_assembleAreas)
        m_areas=new gsMatrix<real_t>();
}

gsRecipeAssemblerQMOpt2D::gsRecipeAssemblerQMOpt2D(const gsMultiPatch<real_t>                &geo,
                                                   const std::vector<gsPhysicalSpace*> phySpaces,
                                                   const gsQualityMeasureWeights             &weights,

                                                   const bool                                 assembleSystemMat,
                                                   const bool                                 assembleRhs,
                                                   const bool                                 assembleQualityMeasureVals,
                                                   const bool                                 assembleAreas,
                                                   const bool                                 assembleQualityDeriv,
                                                   const bool                                 assembleFitting
                                                   )
    : gsRecipeAssembler(geo),m_weights(weights),m_qualityMeasureVals(NULL),m_areas(NULL),
      m_assembleSystemMat(assembleSystemMat),m_assembleRhs(assembleRhs),
      m_assembleQualityMeasureVals(assembleQualityMeasureVals),m_assembleAreas(assembleAreas),
      m_assembleQualityDeriv(assembleQualityDeriv), m_assembleFitting(assembleFitting),
      m_fixBoundaries(false)
{
    setSpace(phySpaces);

    if(m_assembleQualityMeasureVals)
        m_qualityMeasureVals=new gsMatrix<real_t>();
    if(m_assembleAreas)
        m_areas=new gsMatrix<real_t>();
}

void gsRecipeAssemblerQMOpt2D::fixBoundaries  ()
{
    m_fixBoundaries=true;
}

void gsRecipeAssemblerQMOpt2D::init()
{
    if(m_fixBoundaries)
    {
        typedef std::vector< patchSide >::const_iterator const_biterator;
        std::vector<index_t> boundaries;
        for(const_biterator it = m_domain.bBegin();it!=m_domain.bEnd();++it)
        {
            for(size_t i = 0;i<m_space.size();++i)
            {
                std::vector<index_t> boundaries_it = m_space[i]->boundaryDofs(*it,0);
                boundaries.insert(boundaries.end(), boundaries_it.begin(), boundaries_it.end());
            }
        }
        sort( boundaries.begin(), boundaries.end() );
        boundaries.erase( unique( boundaries.begin(), boundaries.end() ), boundaries.end() );
        eliminateDofs(boundaries,0);
    }
    gsRecipeAssembler::init();
}

void gsRecipeAssemblerQMOpt2D::initSystemMatrices()
{
    gsRecipeAssembler::initSystemMatrices();
    if(m_assembleQualityMeasureVals)
        m_qualityMeasureVals->setZero(m_weights.getSize()+1,1);
    if(m_assembleAreas)
        m_areas->setZero(m_weights.getSize()+1,m_domain.nPatches());
}

gsRecipe<real_t>    gsRecipeAssemblerQMOpt2D::getPatchRecipe     (index_t patch)
{
    gsRecipe<real_t>            patch_recipe(0);
    gsRecipeIngredient<real_t>  ingr;

    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);

    if(m_assembleSystemMat)
    {
        ingr.setOperator(new gsQualitySystemMatrixOp<real_t>(m_weights,patch));
        ingr.setRule(getSysWriter());
        patch_recipe.add(ingr);
    }

    if(m_assembleRhs || m_assembleQualityDeriv)
    {
        ingr.setOperator(new gsQualityRHSOp<real_t>(m_weights,patch));
        ingr.setRule(getRhsWriter());
        patch_recipe.add(ingr);
    }

    if(m_assembleQualityMeasureVals)
    {
        ingr.setOperator(new gsQualityMeasureValueOp<real_t>(m_weights,patch));
        ingr.setRule(new gsL2GPlain<gsMatrix<real_t> >(*m_qualityMeasureVals));
        patch_recipe.add(ingr);
    }

    if(m_assembleAreas)
    {
        gsQualityMeasureWeights weights(m_domain.nPatches());
        weights.m_area=1;
        ingr.setOperator(new gsQualityMeasureValueOp<real_t>(weights,patch));
        ingr.setRule(new gsL2GPlain<gsMatrix<real_t> >(*m_areas,0,patch));
        patch_recipe.add(ingr);
    }

    return patch_recipe;
}

bool gsRecipeAssemblerQMOpt2D::getQualityMeasureVals(gsMatrix<real_t> & qualityMeasureVals) const
{
    if(!m_assembleQualityMeasureVals)
        return false;
    qualityMeasureVals=*m_qualityMeasureVals;
    return true;
}

bool gsRecipeAssemblerQMOpt2D::getAreas(gsMatrix<real_t> & areas) const
{
    if(!m_assembleAreas)
        return false;
    areas.resize(m_areas->cols(),1);
    for(index_t i = 0;i<areas.rows();++i)
        areas(i,0)=m_areas->coeff(9,i);
    return true;
}

bool gsRecipeAssemblerQMOpt2D::getQualityDerivs(gsMatrix<real_t> & qualityDerivs) const
{
    if(!m_assembleQualityDeriv)
        return false;
    qualityDerivs=(m_rhs*-2);
    return true;
}

void gsRecipeAssemblerQMOpt2D::postProcess()
{
    if(m_assembleFitting)
        assembleFittingMatrices();
    gsRecipeAssembler::postProcess();
}

void gsRecipeAssemblerQMOpt2D::assembleFittingMatrices()
{
    GISMO_ASSERT(m_weights.m_approxPoints!=0,"only get here if we should do fitting.");
    gsSparseMatrix<real_t> temp_sysMat;
    gsMatrix<real_t> temp_rhs;
    gsLocalToGlobalMapper<real_t>* sysW = getSysWriter();
    gsLocalToGlobalMapper<real_t>* rhsW = getRhsWriter();
    gsMatrix<index_t> dummy;
    for(size_t patchId = 0;patchId<m_domain.nPatches();++patchId)
    {
        if(m_weights.m_pars[patchId].rows()==0||m_weights.m_pars[patchId].cols()==0)
            continue;
        const real_t weight = m_weights.m_approxPoints;
        gsMatrix<real_t> curPoints = m_weights.m_points[patchId]-m_domain.patch(patchId).eval(m_weights.m_pars[patchId]);
        const gsMatrix<real_t>& paramVals = m_weights.m_pars[patchId];
        curPoints.transposeInPlace();
        const index_t num_points = curPoints.rows();

        //for computing the value of the basis function
        gsPhysicalSpace::spacePtr evaluator = m_space[0]->getPatchSpace(patchId);
        gsFuncData<real_t> val(NEED_VALUE|NEED_ACTIVE);

        for(index_t k = 0; k < num_points; k++)
        {
            const gsMatrix<real_t> curr_point = paramVals.col(k);
            evaluator->compute(curr_point,val);

            //computing the values of the basis functions at the current point
            const gsMatrix<real_t>& value=val.values[0];

            // which functions have been computed i.e. which are active
            const gsMatrix<index_t>& actives=val.actives;

            temp_sysMat.resize(actives.rows(),actives.rows());
            temp_sysMat.setZero();
            temp_rhs.setZero(actives.rows(),1);

            const index_t numActive = actives.rows()/2;

            for (index_t i = 0; i != numActive; ++i)
            {
                temp_rhs(i,0)           += weight * ( value(2*i, 0) * curPoints(k,0) );
                temp_rhs(numActive+i,0) += weight * ( value(2*i, 0) * curPoints(k,1) );
                for (index_t j = 0; j != numActive; ++j)
                {
                    temp_sysMat(i, j)                     += weight * ( value(2*i, 0) * value(2*j, 0) );
                    temp_sysMat(numActive+i, numActive+j) += weight * ( value(2*i, 0) * value(2*j, 0) );
                }
            }

            sysW->store(actives,actives,temp_sysMat);
            rhsW->store(actives,dummy,temp_rhs);
        }
    }
}
}


