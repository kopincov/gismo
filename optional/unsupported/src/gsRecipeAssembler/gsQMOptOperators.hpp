/** @file gsQMOptOperators.hpp

    @brief Implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsRecipeAssembler/gsPDEOperators.h>
#include <gsRecipeAssembler/gsQMOptOperators.h>

namespace gismo {


template<typename T>
void gsQualityMeasureOp<T>::getLen(const gsPointFuncData<T>     &testSpace,
                                   std::vector<gsVector<T> >& Jlen) const
{
    GISMO_ASSERT( testSpace.info().first==2,"just implemented in 2D");

    Jlen.clear();
    gsVector<T> JlenMat;
    JlenMat.setZero(testSpace.actives().rows());
    Jlen.push_back(JlenMat);
    Jlen.push_back(JlenMat);
    Jlen.push_back(JlenMat);
    Jlen.push_back(JlenMat);
    for(int i = 0;i<JlenMat.rows();++i)
    {
        Jlen[0](i)=testSpace.jacobian(i)(0,0);
        Jlen[1](i)=testSpace.jacobian(i)(1,0);
        Jlen[2](i)=testSpace.jacobian(i)(0,1);
        Jlen[3](i)=testSpace.jacobian(i)(1,1);
    }
}

template<typename T>
void gsQualityMeasureOp<T>::getUni(const gsPointFuncData<T>     &testSpace,
                                   std::vector<gsVector<T> >& Jlen) const
{
    GISMO_ASSERT( testSpace.info().first==2,"just implemented in 2D");

    Jlen.clear();
    gsVector<T> JlenMat(testSpace.actives().rows());
    Jlen.push_back(JlenMat);
    Jlen.push_back(JlenMat);
    Jlen.push_back(JlenMat);
    Jlen.push_back(JlenMat);
    Jlen.push_back(JlenMat);
    Jlen.push_back(JlenMat);
    gsVector<T>basisDeriv2U,basisDeriv2V,basisDerivUV;
    gsMatrix<T> deriv2InPoint;
    for(int i = 0;i<JlenMat.rows();++i)
    {
        deriv2InPoint = testSpace.deriv2().col(i);
        // TODO check
        // this is very suspicious as the format is 1xx 1yy 1xy 2xx 2yy 2xy
        // and the following gives
        // 1xx 1xy 2yy
        // 1yy 2xx 2xy
        deriv2InPoint.resize(2,3);
        basisDeriv2U = deriv2InPoint.col(0);
        basisDeriv2V = deriv2InPoint.col(1);
        basisDerivUV = deriv2InPoint.col(2);
        Jlen[0](i)=basisDeriv2U(0);
        Jlen[1](i)=basisDeriv2U(1);
        Jlen[2](i)=basisDeriv2V(0);
        Jlen[3](i)=basisDeriv2V(1);
        Jlen[4](i)=basisDerivUV(0);
        Jlen[5](i)=basisDerivUV(1);
    }
}

template<typename T>
void gsQualityMeasureOp<T>::getOrt(const gsPointFuncData<T>     &testSpace,
                                   const gsPointMapData<T>  &geoEval,
                                   gsVector<T>& Jort) const
{
    GISMO_ASSERT( testSpace.info().first==2,"just implemented in 2D");
    Jort.setZero(testSpace.actives().rows());
    for(int i = 0;i<Jort.rows();++i)
    {
        Jort(i)=
                geoEval.jacobian().col(0).dot(testSpace.jacobian(i).col(1))
                +geoEval.jacobian().col(1).dot(testSpace.jacobian(i).col(0));
    }
}

template<typename T>
void gsQualityMeasureOp<T>::getArea(const gsPointFuncData<T> &testSpace,
                                    const gsPointMapData<T>  &geoEval,
                                    gsVector<T>& Jarea) const
{
    GISMO_ASSERT( testSpace.info().first==2,"just implemented in 2D");
    Jarea.setZero(testSpace.actives().rows());
    const gsMatrix<T>geomDeriv = geoEval.jacobian();
    const gsVector<T>geomDerivU = geomDeriv.col(0);
    const gsVector<T>geomDerivV = geomDeriv.col(1);
    gsVector<T> basisDerivU,basisDerivV;
    gsMatrix<T> derivInPoint;
    for(int i = 0;i<Jarea.rows();++i)
    {
        derivInPoint = testSpace.jacobian(i);
        basisDerivU = derivInPoint.col(0);
        basisDerivV = derivInPoint.col(1);
        Jarea(i)=basisDerivU(0)*geomDerivV(1)+geomDerivU(0)*basisDerivV(1)-
                basisDerivU(1)*geomDerivV(0)-geomDerivU(1)*basisDerivV(0);
    }
}

template<typename T>
void gsQualityMeasureOp<T>::getSke(const gsPointFuncData<T>     &testSpace,
                                   const gsPointMapData<T>  &geoEval,
                                   gsVector<T>& Jske) const
{
    GISMO_ASSERT( testSpace.info().first==2,"just implemented in 2D");
    Jske.setZero(testSpace.actives().rows());
    const gsMatrix<T>geomDeriv = geoEval.jacobian();
    const gsVector<T>geomDerivU = geomDeriv.col(0);
    const gsVector<T>geomDerivV = geomDeriv.col(1);
    gsVector<T>basisDerivU,basisDerivV;
    gsMatrix<T> derivInPoint;
    T gugv,gugu,gvgv;
    for(int i = 0;i<Jske.rows();++i)
    {
        derivInPoint = testSpace.jacobian(i);
        basisDerivU = derivInPoint.col(0);
        basisDerivV = derivInPoint.col(1);
        gugv=geomDerivU.dot(geomDerivV);
        gugu=geomDerivU.dot(geomDerivU);
        gvgv=geomDerivV.dot(geomDerivV);
        Jske(i)=(2*gugv*(basisDerivU.dot(geomDerivV)+basisDerivV.dot(geomDerivU))*gugu*gvgv)/math::pow(gugu*gvgv,2)-
                (gugv*gugv)*(2*basisDerivU.dot(geomDerivU)*gvgv+gugu*2*basisDerivV.dot(geomDerivV))/math::pow(gugu*gvgv,2);
    }
}

template<typename T>
void gsQualityMeasureOp<T>::getEcc(const gsPointFuncData<T>     &testSpace,
                                   const gsPointMapData<T>  &geoEval,
                                   std::vector<gsVector<T> >& Jecc) const
{
    const short_t parDim = testSpace.info().first;
    GISMO_ASSERT( parDim==2,"just implemented in 2D");
    Jecc.clear();
    gsVector<T> JeccMat(testSpace.actives().rows());
    Jecc.push_back(JeccMat);
    Jecc.push_back(JeccMat);
    Jecc.push_back(JeccMat);
    Jecc.push_back(JeccMat);
    const gsMatrix<T>geomDeriv = geoEval.jacobian();
    const gsVector<T>geomDerivU = geomDeriv.col(0);
    const gsVector<T>geomDerivV = geomDeriv.col(1);
    gsMatrix<T>geomDeriv2 = geoEval.deriv2();
    geomDeriv2.resize(2,3);
    const gsVector<T>geomDeriv2U = geomDeriv2.col(0);
    const gsVector<T>geomDeriv2V = geomDeriv2.col(1);
    gsVector<T>basisDerivU,basisDerivV,basisDeriv2U,basisDeriv2V;
    gsMatrix<T> derivInPoint,deriv2InPoint;
    T gugu,gvgv,guguu,gvgvv;
    for(int i = 0;i<JeccMat.rows();++i)
    {
        derivInPoint = testSpace.jacobian(i);
        basisDerivU = derivInPoint.col(0);
        basisDerivV = derivInPoint.col(1);
        deriv2InPoint = testSpace.deriv2().middleRows(i*(3*parDim),parDim*3);
        deriv2InPoint.resize(2,3);
        basisDeriv2U = deriv2InPoint.col(0);
        basisDeriv2V = deriv2InPoint.col(1);
        guguu=geomDerivU.dot(geomDeriv2U);
        gvgvv=geomDerivV.dot(geomDeriv2V);
        gugu=geomDerivU.dot(geomDerivU);
        gvgv=geomDerivV.dot(geomDerivV);
        Jecc[0](i)=((basisDerivU.dot(geomDeriv2U)+basisDeriv2U.dot(geomDerivU))*gugu-guguu*2*basisDerivU.dot(geomDerivU))/math::pow(gugu,2);
        Jecc[1](i)=((basisDerivV.dot(geomDeriv2V)+basisDeriv2V.dot(geomDerivV))*gvgv-gvgvv*2*basisDerivV.dot(geomDerivV))/math::pow(gvgv,2);
    }
}

template<typename T>
void gsQualityMeasureOp<T>::getSelfInt(const gsPointFuncData<T>     &testSpace,
                                       const gsPointMapData<T>  &geoEval,
                                       std::vector<gsVector<T> >& JselfInt) const
{
    GISMO_ASSERT( testSpace.info().first==2,"just implemented in 2D");
    JselfInt.clear();
    gsVector<T> JselfIntMat(testSpace.actives().rows());
    JselfInt.push_back(JselfIntMat);
    JselfInt.push_back(JselfIntMat);
    const gsMatrix<T>geomDeriv  = geoEval.jacobian();
    const gsVector<T>geomDerivU = geomDeriv.col(0);
    const gsVector<T>geomDerivV = geomDeriv.col(1);
    gsMatrix<T> derivInPoint;
    gsVector<T>basisDerivU,basisDerivV;
    for(int i = 0;i<JselfIntMat.rows();++i)
    {
        derivInPoint = testSpace.jacobian(i);
        basisDerivU = derivInPoint.col(0);
        basisDerivV = derivInPoint.col(1);
        JselfInt[0](i)=(-2*basisDerivU.dot(geomDerivU))/math::pow(geomDerivU.dot(geomDerivU)+m_weights.m_epsilon,2);
        JselfInt[1](i)=(-2*basisDerivV.dot(geomDerivV))/math::pow(geomDerivV.dot(geomDerivV)+m_weights.m_epsilon,2);
    }
}

template<typename T>
void gsQualityMeasureOp<T>::getAreaInv(const gsPointFuncData<T>     &testSpace,
                                       const gsPointMapData<T>  &geoEval,
                                       gsVector<T>& JareaInv) const
{
    GISMO_ASSERT( testSpace.info().first==2,"just implemented in 2D");
    JareaInv.resize(testSpace.actives().rows());
    const gsMatrix<T>geomDeriv = geoEval.jacobian();
    const gsVector<T>geomDerivU = geomDeriv.col(0);
    const gsVector<T>geomDerivV = geomDeriv.col(1);
    gsVector<T>basisDerivU,basisDerivV;
    gsMatrix<T> derivInPoint;
    for(int i = 0;i<JareaInv.rows();++i)
    {
        derivInPoint = testSpace.jacobian(i);
        basisDerivU = derivInPoint.col(0);
        basisDerivV = derivInPoint.col(1);
        JareaInv(i)=(-1)*math::pow(geomDerivU(0)*geomDerivV(1)-geomDerivU(1)*geomDerivV(0),-2)*
                (basisDerivU(0)*geomDerivV(1)+geomDerivU(0)*basisDerivV(1)-
                 basisDerivU(1)*geomDerivV(0)-geomDerivU(1)*basisDerivV(0));
    }
}



template<typename T>
void  gsQualitySystemMatrixOp<T>::pointEval (
        const gsPointFuncData<T>     &testSpace,
        const gsPointFuncData<T>     &/*unknownSpace*/,
        const gsPointMapData<T>  &geoEval,
        gsRecipeAccumulator<T> result
        ) const
{
    GISMO_ASSERT( testSpace.info().first==2,"just implemented in 2D");
    // returns local matrix
    // nr actives * nr actives
    // unknown  : cols
    // testfunctions : rows
    T patchArea=m_weights.m_areas(m_patchId,0);
    T sqrtPatchArea=math::pow(m_weights.m_areas(m_patchId,0),0.5);
    if(m_weights.m_length>0)
    {
        std::vector<gsVector<T> > Jlen;
        gsQualityMeasureOp<T>::getLen(testSpace,Jlen);
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_length*Jlen[0]*Jlen[0].transpose();
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_length*Jlen[1]*Jlen[1].transpose();
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_length*Jlen[2]*Jlen[2].transpose();
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_length*Jlen[3]*Jlen[3].transpose();
    }
    if(m_weights.m_uniformity>0)
    {
        std::vector<gsVector<T> > Jlen;
        gsQualityMeasureOp<T>::getUni(testSpace,Jlen);
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_uniformity*Jlen[0]*Jlen[0].transpose();
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_uniformity*Jlen[1]*Jlen[1].transpose();
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_uniformity*Jlen[2]*Jlen[2].transpose();
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_uniformity*Jlen[3]*Jlen[3].transpose();
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*2*m_weights.m_uniformity*Jlen[4]*Jlen[4].transpose();
        result+=(1.0/geoEval.measure())*1.0/sqrtPatchArea*2*m_weights.m_uniformity*Jlen[5]*Jlen[5].transpose();
    }
    if(m_weights.m_orthogonality>0)
    {
        gsVector<T> Jort;
        gsQualityMeasureOp<T>::getOrt(testSpace,geoEval,Jort);
        result+=(1.0/geoEval.measure())*(1.0/patchArea)*m_weights.m_orthogonality*Jort*Jort.transpose();
    }
    if(m_weights.m_area>0)
    {
        gsVector<T> Jarea;
        gsQualityMeasureOp<T>::getArea(testSpace,geoEval,Jarea);
        result+=(1.0/geoEval.measure())*(1.0/patchArea)*m_weights.m_area*Jarea*Jarea.transpose();
    }
    if(m_weights.m_skewness>0)
    {
        gsVector<T> Jske;
        gsQualityMeasureOp<T>::getSke(testSpace,geoEval,Jske);
        result+=(1.0/geoEval.measure())*(1.0/patchArea)*m_weights.m_skewness*Jske*Jske.transpose();
    }
    if(m_weights.m_eccentricity>0)
    {
        std::vector<gsVector<T> > Jecc;
        gsQualityMeasureOp<T>::getEcc(testSpace,geoEval,Jecc);
        result+=(1.0/geoEval.measure())*1.0/patchArea*m_weights.m_eccentricity*Jecc[0]*Jecc[0].transpose();
        result+=(1.0/geoEval.measure())*1.0/patchArea*m_weights.m_eccentricity*Jecc[1]*Jecc[1].transpose();
    }
    if(m_weights.m_selfIntersect>0)
    {
        std::vector<gsVector<T> > JselfInt;
        gsQualityMeasureOp<T>::getSelfInt(testSpace,geoEval,JselfInt);
        result+=(1.0/geoEval.measure())*1.0/patchArea*m_weights.m_selfIntersect*JselfInt[0]*JselfInt[0].transpose();
        result+=(1.0/geoEval.measure())*1.0/patchArea*m_weights.m_selfIntersect*JselfInt[1]*JselfInt[1].transpose();
    }
    if(m_weights.m_areaInverse>0)
    {
        gsVector<T> JareaInv;
        gsQualityMeasureOp<T>::getAreaInv(testSpace,geoEval,JareaInv);
        result+=(1.0/geoEval.measure())*(1.0/sqrtPatchArea)*m_weights.m_areaInverse*JareaInv*JareaInv.transpose();
    }
}



template<typename T>
void  gsQualityRHSOp<T>::pointEval (
        const gsPointFuncData<T>     &testSpace,
        const gsPointFuncData<T>     &/*unknownSpace*/,
        const gsPointMapData<T>  &geoEval,
        gsRecipeAccumulator<T> result
        ) const
{
    GISMO_ASSERT( testSpace.info().first==2,"just implemented in 2D");
    // returns local matrix
    // nr actives * nr actives
    // unknown  : cols
    // testfunctions : rows
    //

    T patchArea=m_weights.m_areas(m_patchId,0);
    T sqrtPatchArea=math::pow(m_weights.m_areas(m_patchId,0),0.5);
    if(m_weights.m_length>0)
    {
        std::vector<gsVector<T> > Jlen;
        gsQualityMeasureOp<T>::getLen(testSpace,Jlen);
        const gsMatrix<T>geomDeriv = geoEval.jacobian();
        const gsVector<T>geomDerivU = geomDeriv.col(0);
        const gsVector<T>geomDerivV = geomDeriv.col(1);
        result-=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_length*geomDerivU(0)*Jlen[0];
        result-=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_length*geomDerivU(1)*Jlen[1];
        result-=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_length*geomDerivV(0)*Jlen[2];
        result-=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_length*geomDerivV(1)*Jlen[3];
    }
    if(m_weights.m_uniformity>0)
    {
        std::vector<gsVector<T> > Jlen;
        gsQualityMeasureOp<T>::getUni(testSpace,Jlen);
        gsMatrix<T>geomDeriv2 = geoEval.deriv2();
        geomDeriv2.resize(2,3);
        const gsVector<T>geomDeriv2U = geomDeriv2.col(0);
        const gsVector<T>geomDeriv2V = geomDeriv2.col(1);
        const gsVector<T>geomDerivUV = geomDeriv2.col(2);
        result-=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_uniformity*geomDeriv2U(0)*Jlen[0];
        result-=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_uniformity*geomDeriv2U(1)*Jlen[1];
        result-=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_uniformity*geomDeriv2V(0)*Jlen[2];
        result-=(1.0/geoEval.measure())*1.0/sqrtPatchArea*m_weights.m_uniformity*geomDeriv2V(1)*Jlen[3];
        result-=(1.0/geoEval.measure())*2*1.0/sqrtPatchArea*m_weights.m_uniformity*geomDerivUV(0)*Jlen[4];
        result-=(1.0/geoEval.measure())*2*1.0/sqrtPatchArea*m_weights.m_uniformity*geomDerivUV(1)*Jlen[5];
    }
    if(m_weights.m_orthogonality>0)
    {
        gsVector<T> Jort;
        gsQualityMeasureOp<T>::getOrt(testSpace,geoEval,Jort);
        const gsMatrix<T>geomDeriv = geoEval.jacobian();
        const gsVector<T>geomDerivU = geomDeriv.col(0);
        const gsVector<T>geomDerivV = geomDeriv.col(1);
        result-=(1.0/geoEval.measure())*(1.0/patchArea)*m_weights.m_orthogonality*(geomDerivU.dot(geomDerivV))*Jort;
    }
    if(m_weights.m_area>0)
    {
        gsVector<T> Jarea;
        gsQualityMeasureOp<T>::getArea(testSpace,geoEval,Jarea);
        const gsMatrix<T>geomDeriv = geoEval.jacobian();
        const gsVector<T>geomDerivU = geomDeriv.col(0);
        const gsVector<T>geomDerivV = geomDeriv.col(1);
        result-=(1.0/geoEval.measure())*(1.0/patchArea)*m_weights.m_area*(geomDerivU(0)*geomDerivV(1)-geomDerivU(1)*geomDerivV(0))*Jarea;
    }
    if(m_weights.m_skewness>0)
    {
        gsVector<T> Jske;
        gsQualityMeasureOp<T>::getSke(testSpace,geoEval,Jske);
        const gsMatrix<T>geomDeriv = geoEval.jacobian();
        const gsVector<T>geomDerivU = geomDeriv.col(0);
        const gsVector<T>geomDerivV = geomDeriv.col(1);
        const T gugv = geomDerivU.dot(geomDerivV);
        const T gugu = geomDerivU.dot(geomDerivU);
        const T gvgv = geomDerivV.dot(geomDerivV);
        result-=(1.0/geoEval.measure())*(1.0/patchArea)*m_weights.m_skewness*(math::pow(gugv,2)/(gugu*gvgv))*Jske;
    }
    if(m_weights.m_eccentricity>0)
    {
        std::vector<gsVector<T> > Jecc;
        gsQualityMeasureOp<T>::getEcc(testSpace,geoEval,Jecc);
        const gsMatrix<T>geomDeriv = geoEval.jacobian();
        const gsVector<T>geomDerivU = geomDeriv.col(0);
        const gsVector<T>geomDerivV = geomDeriv.col(1);
        gsMatrix<T>geomDeriv2 = geoEval.deriv2();
        geomDeriv2.resize(2,3);
        const gsVector<T>geomDeriv2U = geomDeriv2.col(0);
        const gsVector<T>geomDeriv2V = geomDeriv2.col(1);
        const T gugu = geomDerivU.dot(geomDerivU);
        const T gvgv = geomDerivV.dot(geomDerivV);
        const T guguu = geomDerivU.dot(geomDeriv2U);
        const T gvgvv = geomDerivV.dot(geomDeriv2V);
        result-=(1.0/geoEval.measure())*1.0/patchArea*m_weights.m_eccentricity*(guguu/gugu)*Jecc[0];
        result-=(1.0/geoEval.measure())*1.0/patchArea*m_weights.m_eccentricity*(gvgvv/gvgv)*Jecc[1];
    }
    if(m_weights.m_selfIntersect>0)
    {
        std::vector<gsVector<T> > JselfInt;
        gsQualityMeasureOp<T>::getSelfInt(testSpace,geoEval,JselfInt);
        const gsMatrix<T>geomDeriv = geoEval.jacobian();
        const gsVector<T>geomDerivU = geomDeriv.col(0);
        const gsVector<T>geomDerivV = geomDeriv.col(1);
        result-=(1.0/geoEval.measure())*1.0/patchArea*m_weights.m_selfIntersect*
                (1/math::pow(geomDerivU.dot(geomDerivU)+m_weights.m_epsilon,2))*JselfInt[0];
        result-=(1.0/geoEval.measure())*1.0/patchArea*m_weights.m_selfIntersect*
                (1/math::pow(geomDerivV.dot(geomDerivV)+m_weights.m_epsilon,2))*JselfInt[1];
    }
    if(m_weights.m_areaInverse>0)
    {
        gsVector<T> JareaInv;
        gsQualityMeasureOp<T>::getAreaInv(testSpace,geoEval,JareaInv);
        const gsMatrix<T>geomDeriv = geoEval.jacobian();
        const gsVector<T>geomDerivU = geomDeriv.col(0);
        const gsVector<T>geomDerivV = geomDeriv.col(1);
        result-=(1.0/geoEval.measure())*(1.0/sqrtPatchArea)*m_weights.m_areaInverse/(geomDerivU(0)*geomDerivV(1)-geomDerivU(1)*geomDerivV(0))*JareaInv;
    }
}


template<typename T>
void  gsQualityMeasureValueOp<T>::pointEval (
        const gsPointFuncData<T>     &/*testSpace*/,
        const gsPointFuncData<T>     &/*unknownSpace*/,
        const gsPointMapData<T>  &geoEval,
        gsRecipeAccumulator<T> result
        ) const
{
    gsMatrix<T> locMat;
    locMat.setZero(m_weights.getSize()+1,1);
    // first deriv:
    const gsMatrix<T>geomDeriv  = geoEval.jacobian();
    const gsVector<T>geomDerivU = geomDeriv.col(0);
    const gsVector<T>geomDerivV = geomDeriv.col(1);
    // second deriv:
    gsMatrix<T>geomDeriv2 = geoEval.deriv2();
    geomDeriv2.resize(2,3);
    const gsVector<T>geomDeriv2U = geomDeriv2.col(0);
    const gsVector<T>geomDeriv2V = geomDeriv2.col(1);
    const gsVector<T>geomDerivUV = geomDeriv2.col(2);
    // helper values:
    const T gugv = geomDerivU.dot(geomDerivV);
    const T gugu = geomDerivU.dot(geomDerivU);
    const T gvgv = geomDerivV.dot(geomDerivV);
    const T guguu = geomDerivU.dot(geomDeriv2U);
    const T gvgvv = geomDerivV.dot(geomDeriv2V);

    T patchArea=m_weights.m_areas(m_patchId,0);
    T sqrtPatchArea=math::pow(m_weights.m_areas(m_patchId,0),0.5);
    // compute vals:
    // length:
    locMat(1,0)+=math::pow(geomDerivU(0),2)+math::pow(geomDerivU(1),2)+
            math::pow(geomDerivV(0),2)+math::pow(geomDerivV(1),2);
    // uniformity:
    locMat(2,0)+=math::pow(geomDeriv2U(0),2)+math::pow(geomDeriv2U(1),2)+
            math::pow(geomDeriv2V(0),2)+math::pow(geomDeriv2V(1),2)+
            2*math::pow(geomDerivUV(0),2)+2*math::pow(geomDerivUV(1),2);
    // orthognality:
    locMat(3,0)+=math::pow(geomDerivU.dot(geomDerivV),2);
    // area-squared:
    locMat(4,0)+=math::pow(geomDerivU(0)*geomDerivV(1)-geomDerivU(1)*geomDerivV(0),2);
    // skewness:
    locMat(5,0)+=math::pow(math::pow(gugv,2)/(gugu*gvgv),2);
    // eccentricity:
    locMat(6,0)+=math::pow(guguu/gugu,2)+math::pow(gvgvv/gvgv,2);
    // self-intersections:
    locMat(7,0)+=1/math::pow(gugu+m_weights.m_epsilon,2)+
            1/math::pow(gvgv+m_weights.m_epsilon,2);
    // area-inverse:
    locMat(8,0)+=math::pow(locMat(4,0),-1);
    //area
    locMat(9,0)+=math::abs(geomDerivU(0)*geomDerivV(1)-geomDerivU(1)*geomDerivV(0));
    // total quality measures + weights:
    locMat(0,0)+=1.0/sqrtPatchArea*m_weights.m_length*locMat(1,0)+1.0/sqrtPatchArea*m_weights.m_uniformity*locMat(2,0)+
            (1.0/patchArea)*m_weights.m_orthogonality*locMat(3,0)+(1.0/patchArea)*m_weights.m_area*locMat(4,0)+
            (1.0/patchArea)*m_weights.m_skewness*locMat(5,0)+1.0/patchArea*m_weights.m_eccentricity*locMat(6,0)+
            1.0/patchArea*m_weights.m_selfIntersect*locMat(7,0)+(1.0/sqrtPatchArea)*m_weights.m_areaInverse*locMat(8,0);
    result+=(1.0/geoEval.measure())*locMat;
}


template<typename T>
void   gsQualityConstraintOp<T>::pointEval (
        const gsPointFuncData<T>     &testSpace,
        const gsPointFuncData<T>     &/*unknownSpace*/,
        const gsPointMapData<T>  &geoEval,
        gsRecipeAccumulator<T> result
        ) const
{
    gsVector<T> Jarea;
    gsQualityMeasureOp<T>::getArea(testSpace,geoEval,Jarea);
    for(int i =0;i<Jarea.rows();++i)
        Jarea(i,0)=Jarea(i,0)<0?-Jarea(i,0):Jarea(i,0);
    result+=(1.0/geoEval.measure())*Jarea;
}





} // end namespace gismo
