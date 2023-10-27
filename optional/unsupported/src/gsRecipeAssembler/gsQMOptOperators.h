/** @file gsQMOptOperators.h

    @brief Provides operators for the recipe assembler for optimizing certain
    quadratic functionals to improve parameterizations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsAssembler/gsSeminormH2.h>
#include <gsRecipeAssembler/gsPDEOperators.h>

namespace gismo
{


//========================================= QUALITY WEIGHTS ==========================================//

struct GISMO_EXPORT gsQualityMeasureWeights{
    gsQualityMeasureWeights(unsigned nPatches)
    {
        m_length = 0;
        m_orthogonality = 0;
        m_uniformity = 0;
        m_skewness = 0;
        m_eccentricity = 0;
        m_area = 0;
        m_selfIntersect = 0;
        m_areaInverse = 0;

        m_epsilon = 0;
        m_areas.setOnes(nPatches,1);
        m_weightArea=false;

        m_approxPoints = 0;
        gsMatrix<real_t> p;
        for(unsigned i=0;i<nPatches;++i)
        {
            m_points.push_back(p);
            m_pars.push_back(p);
        }
    }

public:
    real_t m_orthogonality;
    real_t m_skewness;
    real_t m_eccentricity;
    real_t m_uniformity;
    real_t m_length;
    real_t m_area;
    real_t m_selfIntersect;
    real_t m_areaInverse;

    real_t m_epsilon;
    gsMatrix<real_t> m_areas;
    bool m_weightArea;

    real_t m_approxPoints;
    std::vector<gsMatrix<real_t> > m_points;
    std::vector<gsMatrix<real_t> > m_pars;

    size_t getSize() const { return 9; }
};

//============================================ OPERATORS =============================================//

template<typename T>
class gsQualityMeasureOp : public gsBilinearOp<T>
{
protected:
    const gsQualityMeasureWeights m_weights;
    const unsigned m_patchId;

public:
    gsQualityMeasureOp(const gsQualityMeasureWeights weights,unsigned patchId)
        : m_weights(weights), m_patchId(patchId)
    {  }

public:
    virtual unsigned  testSpaceNeeds() const {return NEED_VALUE | NEED_GRAD | NEED_2ND_DER;}
    virtual unsigned  unknownSpaceNeeds() const {return 0;}
    virtual unsigned  geometryNeeds()    const {return  NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;}

protected:
    /// Look at class gsQualityMeasure
    void getLen(const gsPointFuncData<T>     &testSpace,
                std::vector<gsVector<T> >& Jlen) const;

    /// Look at class gsQualityMeasure
    void getUni(const gsPointFuncData<T>     &testSpace,
                std::vector<gsVector<T> >& Jlen) const;

    /// Look at class gsQualityMeasure
    void getOrt(const gsPointFuncData<T>     &testSpace,
                const gsPointMapData<T>  &geoEval,
                gsVector<T>& Jort) const;

    /// Look at class gsQualityMeasure
    void getArea(const gsPointFuncData<T>     &testSpace,
                 const gsPointMapData<T>  &geoEval,
                 gsVector<T>& Jarea) const;

    /// Look at class gsQualityMeasure
    void getSke(const gsPointFuncData<T>     &testSpace,
                const gsPointMapData<T>  &geoEval,
                gsVector<T>& Jske) const;

    /// Look at class gsQualityMeasure
    void getEcc(const gsPointFuncData<T>     &testSpace,
                const gsPointMapData<T>  &geoEval,
                std::vector<gsVector<T> >& Jecc) const;

    /// Look at class gsQualityMeasure
    void getSelfInt(const gsPointFuncData<T>     &testSpace,
                    const gsPointMapData<T>  &geoEval,
                    std::vector<gsVector<T> >& JselfInt) const;

    /// Look at class gsQualityMeasure
    void getAreaInv(const gsPointFuncData<T>     &testSpace,
                    const gsPointMapData<T>  &geoEval,
                    gsVector<T>& JareaInv) const;
};

template<typename T>
class gsQualitySystemMatrixOp : public gsQualityMeasureOp<T>
{
    using gsQualityMeasureOp<T>::m_weights;
    using gsQualityMeasureOp<T>::m_patchId;
public:
    gsQualitySystemMatrixOp(const gsQualityMeasureWeights weights,unsigned patchId)
        : gsQualityMeasureOp<T>(weights,patchId)
    {  }

public:
    virtual void  pointEval (
                const gsPointFuncData<T>  &testSpace,
                const gsPointFuncData<T>  &unknownSpace,
                const gsPointMapData <T>  &geoEval,
                gsRecipeAccumulator<T> result
                ) const;
};

template<typename T>
class gsQualityRHSOp : public gsQualityMeasureOp<T>
{
    using gsQualityMeasureOp<T>::m_weights;
    using gsQualityMeasureOp<T>::m_patchId;
public:
    gsQualityRHSOp(const gsQualityMeasureWeights weights,unsigned patchId)
        : gsQualityMeasureOp<T>(weights,patchId)
    {  }

public:
    virtual void  pointEval (
                const gsPointFuncData<T>  &testSpace,
                const gsPointFuncData<T>  &unknownSpace,
                const gsPointMapData <T>  &geoEval,
                gsRecipeAccumulator<T>     result
                ) const;
    virtual void      outputSize (unsigned & /*r*/,
                                  unsigned & c) const
    { c=1; }

};

template<typename T>
class gsQualityMeasureValueOp : public gsBilinearOp<T>
{
protected:
    const gsQualityMeasureWeights m_weights;
    const unsigned m_patchId;

public:
    gsQualityMeasureValueOp(const gsQualityMeasureWeights weights,unsigned patchId)
        : m_weights(weights),m_patchId(patchId)
    {  }

public:
    virtual unsigned  testSpaceNeeds() const {return NEED_VALUE | NEED_GRAD | NEED_2ND_DER;}
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE | NEED_GRAD | NEED_2ND_DER;}
    virtual unsigned  geometryNeeds()    const {return  NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;}

public:
    virtual void  pointEval (
                const gsPointFuncData<T>  &testSpace,
                const gsPointFuncData<T>  &unknownSpace,
                const gsPointMapData <T>  &geoEval,
                gsRecipeAccumulator<T>     result
                ) const;
    virtual void      outputSize (unsigned &r, unsigned &c) const
    {
        r=m_weights.getSize()+1;
        c=1;
    }
};

template<typename T>
class gsQualityConstraintOp : public gsQualityMeasureOp<T>
{
    using gsQualityMeasureOp<T>::m_weights;
    using gsQualityMeasureOp<T>::m_patchId;
public:
    gsQualityConstraintOp(const gsQualityMeasureWeights weights,unsigned patchId)
        : gsQualityMeasureOp<T>(weights,patchId)
    {  }

public:
    virtual void  pointEval (
                const gsPointFuncData<T>  &testSpace,
                const gsPointFuncData<T>  &unknownSpace,
                const gsPointMapData <T>  &geoEval,
                    gsRecipeAccumulator<T>     result
                ) const;
    virtual void      outputSize (unsigned & /*r*/,
                                  unsigned & c) const
    { c=1; }
};

} // end namespace gismo

