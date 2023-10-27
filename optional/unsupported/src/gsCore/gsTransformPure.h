/**
 * @file gsTransformPure.hpp
 *
 * @brief Implementation of the mapping of functions between the parametric
 *   and the physical domain.
 *
 *   This file is part of the G+Smo library.
 *
 *   This Source Code Form is subject to the terms of the Mozilla Public
 *   License, v. 2.0. If a copy of the MPL was not distributed with this
 *   file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *   Author(s): A. Bressan, A. Manzaflaris
**/

#pragma once

#include <gsCore/gsFuncData.h>
#include <gsCore/gsEvalUtils.hpp>

namespace gismo {


/*
 * COMMON SCHEMA OF A PURE TRANSFORMATION OBJECT
template <typename T>
class gsTransformPure // No transformation
{
public:
    static unsigned getFuncFlags  (unsigned flags) {return flags;}
    static unsigned getGeoFlags (unsigned flags)   {return 0;}

    static void transform   (const gsMapData<T> &geoData,
                                    const gsFunctionSet<T> &parmFunc,
                                    gsFuncData<T> &phisFunc)
    {
        ...
    }

};
*/

template <typename T>
class gsTransformPureIdentity // No transformation
{
public:
    static unsigned getFuncFlags  (unsigned flags) {return flags;}
    static unsigned getGeoFlags   (unsigned ) {return 0;}
    static void transform   (const gsMapData<T>  &geoData,
                             const gsFunctionSet<T> &FF,
                             gsFuncData<T> &phisFunc)
    {
        FF.compute(geoData.points,phisFunc);
    }
};


template <typename T>
class gsTransformPureGradConforming
{
public:
    static unsigned getFuncFlags (unsigned flags)
    {
        unsigned result=flags& ~(NEED_CURL|NEED_DIV|NEED_LAPLACIAN);
        result |= flags&(NEED_LAPLACIAN|NEED_DERIV2)     ? NEED_DERIV|NEED_DERIV2 : 0;
        result |= flags&(NEED_CURL|NEED_DIV)             ? NEED_DERIV  : 0;
        return result;
    }
    static unsigned getGeoFlags  (unsigned flags)
    {
        unsigned result=0;
        result |= flags&NEED_VALUE     ? 0 : 0;
        result |= flags&NEED_DERIV     ? NEED_GRAD_TRANSFORM : 0;
        result |= flags&NEED_DERIV2    ? NEED_GRAD_TRANSFORM|NEED_DERIV2 : 0;
        result |= flags&NEED_DIV       ? NEED_GRAD_TRANSFORM : 0;
        result |= flags&NEED_CURL      ? NEED_GRAD_TRANSFORM : 0;
        result |= flags&NEED_LAPLACIAN ? NEED_GRAD_TRANSFORM|NEED_DERIV2 : 0;
        return result;
    }

    static void transform   (const gsMapData<T>  &geoData,
                             const gsFunctionSet<T> &FF,
                             gsFuncData<T> &phisFunc)
    {
        gsFuncData<T> parmFunc(getFuncFlags(phisFunc.flags));
        parmFunc.patchId=phisFunc.patchId;
        FF.compute(geoData.points,parmFunc);

        const unsigned flags=phisFunc.flags;
        phisFunc.dim.first  = geoData.dim.second;
        phisFunc.dim.second = parmFunc.dim.second;

        if (flags&NEED_ACTIVE)
            phisFunc.actives=parmFunc.actives;
        phisFunc.values.resize(parmFunc.values.size());
        if (flags&NEED_VALUE)
            phisFunc.values[0]=parmFunc.values[0];
        if (flags&(NEED_DERIV|NEED_DIV|NEED_CURL) )
        {
            transformDeriv(geoData,parmFunc,phisFunc);
            if (flags&NEED_DIV)
                convertValue<T>::derivToDiv(phisFunc.values[1],phisFunc.dim,phisFunc.divs);
            if (flags&NEED_CURL)
                convertValue<T>::derivToCurl(phisFunc.values[1],phisFunc.dim,phisFunc.curls);
        }
        if (flags&(NEED_DERIV2|NEED_LAPLACIAN) )
        {
            transformDeriv2(geoData,parmFunc,phisFunc);
            if (flags&NEED_LAPLACIAN)
                convertValue<T>::deriv2ToLaplacian(phisFunc.values[2],phisFunc.dim,phisFunc.laplacians);
        }
    }

    static inline void transformDeriv(const gsMapData<T>  &geoData,
                                      const gsFuncData<T> &parmFunc,
                                      gsFuncData<T> &phisFunc)
    {
        const index_t numPoints=geoData.points.cols();
        const index_t paramSize=parmFunc.dim.first;
        const index_t physiSize=phisFunc.dim.first;

        const index_t derivNum =parmFunc.values[1].rows()/paramSize;

        phisFunc.values[1].resize(physiSize*derivNum,numPoints);
        for (index_t p=0; p<numPoints; ++p)
        {
            phisFunc.values[1].reshapeCol(p,physiSize,derivNum)=
                    geoData.fundForm(p)*parmFunc.values[1].reshapeCol(p,paramSize,derivNum);
        }
    }

    static inline void transformDeriv2(const gsMapData<T>  &geoData,
                                       const gsFuncData<T> &prmFunc,
                                       gsFuncData<T> &phsFunc)
    {
        const index_t numPoints = geoData.points.cols();
        const index_t prmDomDim = prmFunc.dim.first;
        const index_t phsDomDim = phsFunc.dim.first;

        const index_t prmSecDerSize = (prmDomDim*(prmDomDim+1))/2;
        const index_t phsSecDerSize = (phsDomDim*(phsDomDim+1))/2;

        const index_t numBlocks = prmFunc.values[2].rows()/prmSecDerSize;

        phsFunc.values[2].resize(phsSecDerSize*numBlocks,numPoints);

        gsMatrix<T> cmpHessian(prmDomDim,prmDomDim);        // hessian of a component of the function
        gsMatrix<T> geoHessianMod(phsSecDerSize,phsDomDim); // hessian of a component of the parametrization

        for (index_t p=0; p<numPoints; ++p)
        {
            // First part: J^-T H J^-1
            for (index_t b=0; b<numBlocks;++b)
            {
                convertValue<T>::deriv2ToHessianSingle(
                            prmFunc.values[2].block(b*prmSecDerSize,p,prmSecDerSize,1),
                        cmpHessian,prmDomDim);
                convertValue<T>::hessianToDeriv2Single(
                            (geoData.fundForm(p)*cmpHessian*geoData.fundForm(p).transpose()).eval(),
                            phsFunc.values[2].block(b*phsSecDerSize,p,phsSecDerSize,1), phsDomDim);
            }
            // Second part: J^-T G J^-1 DDG J^-1
            for (index_t c=0;c<phsDomDim; ++c)
            {
                convertValue<T>::deriv2ToHessianSingle(
                            geoData.values[2].block(c*prmSecDerSize,p,prmSecDerSize,1),
                        cmpHessian,prmDomDim);
                convertValue<T>::hessianToDeriv2Single(
                            (geoData.fundForm(p)*cmpHessian*geoData.fundForm(p).transpose()).eval(),
                            geoHessianMod.col(c),
                            phsDomDim);
            }
            phsFunc.values[2].reshapeCol(p,phsSecDerSize,numBlocks).noalias() -=
                    geoHessianMod*geoData.fundForm(p)
                    *(prmFunc.values[1].reshapeCol(p,prmDomDim,numBlocks));
        }
    }

};


template <typename T>
class gsTransformPureDivConforming
{
public:
    static unsigned getFuncFlags (unsigned flags)
    {
        unsigned result = flags&(NEED_VALUE|NEED_ACTIVE|NEED_DERIV|NEED_DIV|SAME_ELEMENT);
        result |= flags&NEED_CURL   ? NEED_DERIV : 0;
        if (flags&NEED_DERIV2 || flags&NEED_LAPLACIAN)
            GISMO_ERROR("cannot computer 3rd derivative");
        return result;
    }
    static unsigned getGeoFlags  (unsigned flags)
    {
        unsigned result=0;
        result |= flags&NEED_VALUE  ? NEED_GRAD_TRANSFORM|NEED_MEASURE : 0;

        flags  |= flags&NEED_CURL   ? NEED_DERIV : 0;
        result |= flags&NEED_DERIV  ? NEED_GRAD_TRANSFORM|NEED_DERIV2 : 0;

        if (flags&NEED_DERIV2 || flags&NEED_LAPLACIAN) GISMO_ERROR("cannot computer 3rd derivative");
        return result;
    }

    static void transform   (const gsMapData<T>  &geoData,
                             const gsFunctionSet<T> &FF,
                             gsFuncData<T> &phisFunc)
    {
        gsFuncData<T> parmFunc(getFuncFlags(phisFunc.flags));
        parmFunc.patchId=phisFunc.patchId;
        FF.compute(geoData.points,parmFunc);

        const unsigned flags=phisFunc.flags;
        phisFunc.dim.first  = geoData.dim.second;
        phisFunc.dim.second = parmFunc.dim.second;

        phisFunc.values.resize(parmFunc.values.size());
        if (flags&NEED_ACTIVE)
            phisFunc.actives=parmFunc.actives;

        if (flags&NEED_VALUE)
            transformValues(geoData,parmFunc,phisFunc);
        if (flags&(NEED_DERIV|NEED_CURL) )
        {
            transformDeriv(geoData,parmFunc,phisFunc);
            if (flags&NEED_CURL)
                convertValue<T>::derivToCurl(phisFunc.values[1],phisFunc.dim,phisFunc.curls);
        }
        if (flags&NEED_DIV)
            transformDiv(geoData,parmFunc,phisFunc); // optimized without computing the full deriv transformation

        if (flags&NEED_DERIV2)
        {
            transformDeriv2(geoData,parmFunc,phisFunc);
            if (flags&NEED_LAPLACIAN)
                convertValue<T>::deriv2ToLaplacian(phisFunc.values[2],phisFunc.dim,phisFunc.laplacians);
        }
    }

    static inline void transformValues(const gsMapData<T>  &geoData,
                                       const gsFuncData<T> &parmFunc,
                                       gsFuncData<T> &phisFunc)
    {
        const index_t paramSize=parmFunc.dim.first;
        const index_t physiSize=phisFunc.dim.first;
        const index_t numValues = parmFunc.values[0].rows()/paramSize;
        const index_t numPoints = geoData.points.cols();

        phisFunc.values[0].resize(numValues*physiSize,numPoints);
        for (index_t p=0; p< numPoints;++p)
        {
            phisFunc.values[0].reshapeCol(p,physiSize,numValues)=
                    geoData.jacobian(p)*parmFunc.values[0].reshapeCol(p,paramSize,numValues)/geoData.measures(0,p)
                    *(geoData.dim.second==geoData.dim.first && geoData.jacobian(p).determinant()<0 ? -1 :1 )
                    ;
        }
    }

    static inline void transformDeriv(const gsMapData<T>  &geoData,
                                      const gsFuncData<T> &prmFunc,
                                      gsFuncData<T> &phsFunc)
    {
        // In Einstein notation this is
        // F(i)    = JG(i,l) w(l) det^-1
        // DF(i,j) = HG(i,l,k) J^-1(k,j) w(l) det^-1
        //           + JG(i,l) Dw(l,k) J^-1(k,j) det^-1
        //           - JG(i,l) w(l) Ddet(j) det^-2
        // where
        //  HG(i,l,k) is the derivative of the i-th component of G in dir k and l
        //  JG(i,l)   is the derivative of the i-th component of G in dir l
        //  Dw(l,k)   is the derivative of the l-th component of w in dir k

        // The code below is less clear as it must de-encode the Hessian of the
        // geometry and must loop over the active in order to be able to see
        // their parametric derivative as a matrix.

        typedef gsAsMatrix<T,-1,-1> matrixView;
        typedef gsAsConstMatrix<T,-1,-1> matrixViewC;
        typedef Eigen::Transpose<typename matrixViewC::Base> matrixViewCT;

        const index_t numPts = geoData.points.cols();
        const index_t prmDim = geoData.dim.first;
        const index_t phsDim = geoData.dim.second;
        const index_t secDer = ((prmDim+1)*prmDim)/2;
        const index_t numAct = prmFunc.values[1].rows()/prmDim/prmDim;


        gsMatrix<T>   HG  (prmDim*prmDim,phsDim);
        gsMatrix<T>   Ddet(prmDim,1);
        real_t        det;

        phsFunc.values[1].resize(phsDim*phsDim*numAct,numPts);

        for (index_t p=0; p<numPts;++p)
        {
            matrixViewC   W(prmFunc.values[0].col(p).data(), prmDim, numAct);
            matrixView   DF(phsFunc.values[1].col(p).data(), phsDim*phsDim, numAct);

            matrixViewCT JG=matrixViewC(geoData.values[1].col(p).data(), phsDim, prmDim).transpose();
            matrixViewC  JI(geoData.fundForms.col(p).data(), phsDim, prmDim); // inverse of JG transposed

            // init HG
            for (index_t i=0; i<prmDim;++i)
                convertValue<T>::deriv2ToHessianSingle(geoData.values[2].block(secDer*i,p,secDer,1),matrixView(HG.col(i).data(),prmDim,prmDim),prmDim);

            // compute det
            det=geoData.measures(0,p); // should add orientation data to geoData

            // compute Ddet
            if (prmDim==phsDim)
            {
                Ddet.setZero(prmDim,1);
                gsMatrix<T> tmp=JG;
                for (int d=0;d<prmDim;++d)
                {
                    for (int i=0; i<prmDim;++i)
                    {
                        for (int j=0;j<prmDim;++j)
                            tmp(j,i)=matrixView(HG.col(j).data(),prmDim,prmDim)(d,i);
                        Ddet(d)+=tmp.determinant();
                        for (int j=0;j<prmDim;++j) // strange that Eigen refuses  tmp.col(i)=JG.col(I)
                            tmp(j,i)=JG(j,i);
                    }
                }
            }
            else
                GISMO_ERROR("The implementation of the derivative of the determinant is missing for domainDim!=targetDim");


            for (index_t i=0;i<phsDim;++i)
            {
                gsMatrix<T> tmp(1,prmDim); // tmp=JG.row(i);
                // TODO find out why using JG.row(i) fails
                for (int j=0;j<prmDim;++j)
                    tmp(0,j)=JG(i,j);
                // doing a j block at a time
                DF.block(i*phsDim,0,phsDim,numAct) =
                        (JI*(matrixView(HG.col(i).data(),prmDim,prmDim)-Ddet*tmp/det)).eval()*W/det
                        ;
                // but we need to multiply DW from both sides so we must (for the moment) do a loop over the actives
                for (index_t a=0; a<numAct; ++a)
                {
                    matrixViewC DWT(prmFunc.values[1].col(p).data()+a*prmDim*prmDim, prmDim, prmDim);
                    DF.block(i*phsDim,a,phsDim,1)+= JI * DWT * tmp.transpose() / det;
                }
            }
        }
    }

    static inline void transformDeriv2(const gsMapData<T>  &geoData,
                                       const gsFuncData<T> &parmFunc,
                                       gsFuncData<T> &phisFunc)
    {
        GISMO_UNUSED(geoData); GISMO_UNUSED(parmFunc); GISMO_UNUSED(phisFunc);
        GISMO_ERROR("transformDeriv2 NOT IMPLEMENTED for div conforming");
    }
    static inline void transformDiv(const gsMapData<T>  &geoData,
                                    const gsFuncData<T> &parmFunc,
                                    gsFuncData<T> &phisFunc)
    {
        phisFunc.divs= parmFunc.divs*geoData.measures.asVector().asDiagonal().inverse();
    }
};


/**
 * To restrict a function defined on \f$\mathbb{R}^n\f$ to a lower dimensional
 * manifold parametrized by the geometry.
 *
 * This transformation project the derivatives on the tangent plane so that
 * the orthogonal component of the derivatives is NULL.
 *
 * The use is to compute norms on surfaces.
 */
template <typename T>
class gsTransformPureRestriction
{
public:
    static unsigned getFuncFlags (unsigned flags)
    {
        unsigned result=flags& ~(NEED_CURL|NEED_DIV|NEED_LAPLACIAN);
        result |= flags&NEED_CURL         ? NEED_DERIV  : 0;
        result |= flags&NEED_DERIV        ? NEED_DERIV  : 0;
        result |= flags&NEED_LAPLACIAN    ? NEED_DERIV2 : 0;
        return result;
    }
    static unsigned getGeoFlags  (unsigned flags)
    {
        unsigned result=0;
        result =  (flags& ~(SAME_ELEMENT)? NEED_VALUE : 0);
        result |= (flags& (NEED_DIV|NEED_CURL|NEED_LAPLACIAN|NEED_DERIV|NEED_DERIV2)) ? NEED_DERIV|NEED_GRAD_TRANSFORM : 0;
        return result;
    }

    static void transform   (const gsMapData<T>     &geoData,
                             const gsFunctionSet<T> &global,
                             gsFuncData<T>       &restricted)
    {
        gsFuncData<T> original(getFuncFlags(restricted.flags));
        original.patchId=restricted.patchId;
        global.compute(geoData.values[0],original);

        const unsigned flags=restricted.flags;
        restricted.dim.first=geoData.dim.second;
        restricted.dim.second=original.dim.second;

        if (flags&NEED_ACTIVE)
            restricted.actives=original.actives;

        restricted.values.resize(original.values.size());
        if (flags&NEED_VALUE)
            restricted.values[0]=original.values[0];
        if (flags&NEED_DERIV)
        {
            transformDeriv(geoData,original,restricted);
            if (flags&NEED_DIV)
                convertValue<T>::derivToDiv(restricted.values[1],restricted.dim, restricted.divs);
            if (flags&NEED_CURL)
                convertValue<T>::derivToCurl(restricted.values[1],restricted.dim,restricted.curls);
        }

        if (flags&NEED_DERIV2)
        {
            transformDeriv2(geoData,original,restricted);
            if (flags&NEED_LAPLACIAN)
                convertValue<T>::deriv2ToLaplacian(restricted.values[2],restricted.dim,restricted.laplacians);
        }
    }

    static inline void transformDeriv(const gsMapData<T>  &geoData,
                                      const gsFuncData<T> &original,
                                      gsFuncData<T>       &restricted)
    {
        const index_t numPoints=geoData.points.cols();
        const index_t derivSize=original.dim.first;
        const index_t derivNum =original.values[1].rows()/derivSize;

        restricted.values[1].resizeLike(original.values[1]);
        for (index_t p=0; p<numPoints; ++p)
        {
            restricted.values[1].reshapeCol(p,derivSize,derivNum)=
                    geoData.fundForm(p)*geoData.jacobian(p).transpose()*original.values[1].reshapeCol(p,derivSize,derivNum);
        }
    }

    static inline void transformDeriv2(const gsMapData<T>  &geoData,
                                       const gsFuncData<T> &original,
                                       gsFuncData<T>       &restricted)
    {

        const index_t numPoints=geoData.points.cols();
        const index_t domainDim=original.dim.first;
        const index_t derivSize=original.dim.first;
        const index_t numBlocks=original.values[2].rows()/derivSize;

        restricted.values[1].resizeLike(original.values[2]);

        gsMatrix<T> temp;
        gsMatrix<T> cmpHessian;

        for (index_t p=0; p<numPoints; ++p)
        {
            temp=geoData.fundForm(p)*geoData.jacobian(p); // projection matrix
            for (index_t b=0; b<numBlocks;++b)
            {
                convertValue<T>::deriv2ToHessianSingle(
                            original.values[2].block(b*derivSize,p,derivSize,1),
                        cmpHessian,domainDim);
                convertValue<T>::hessianToDeriv2Single(
                            temp.transpose()*cmpHessian*temp,
                            restricted.values[2].block(b*derivSize,p,derivSize,1), domainDim);
            }
        }
    }
};



} // namespace gismo
