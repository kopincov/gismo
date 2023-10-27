/** @file gsDistanceOperators.h

    declarations for distance operators

    Here there is an interface and some implementations that represent
    the act of adding the element-wise contribution of a matrix to
    the global object.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsRecipeAssembler/gsRecipe.h>
#include <gsRecipeAssembler/gsPDEOperators.h>


namespace gismo {

/**
 * @brief The gsDistanceOperator class
 *
 * base class to compute a distance, this is a temporary solution as the proper
 * way is to move coefficients in the evaluator.
 */
class GISMO_EXPORT gsDistanceOperator : public gsBilinearOp<real_t>
{
protected:
    const gsMatrix<real_t>    &m_coef;
    const gsFunction<real_t>  *m_func;

public:
    gsDistanceOperator(const gsMatrix<real_t> &coef, const gsFunction<real_t> *func=NULL)
        :   m_coef(coef), m_func(func)
    {}
    virtual unsigned  unknownSpaceNeeds() const = 0;
    virtual unsigned  geometryNeeds()      const = 0;

    virtual unsigned  testSpaceNeeds() const {return 0;}
    virtual void      outputSize (unsigned &r, unsigned &c) const { r=1; c=m_coef.cols();}
protected:
    // utility function that returns the transposed local coefs
    inline gsMatrix<real_t> getLocCoef(const gsMatrix<index_t> &active) const
    {
        return getLocCoef(active,m_coef);
    }
    static inline gsMatrix<real_t> getLocCoef(const gsMatrix<index_t> &active, const gsMatrix<real_t> &coefs)
    {
        const index_t nRows=active.rows();
        const index_t nCols=coefs.cols();
        gsMatrix<real_t> result(nRows,nCols);
        for(index_t r=0;r<nRows;++r)
            result.row(r)=coefs.row(active(r,0)).transpose();
        return result;
    }


    template <typename T>
    inline gsMatrix<real_t> computeLocValues2 (const T &val, const gsMatrix<index_t> &active) const
    {
        const gsMatrix<real_t> lCoefs=getLocCoef(active,m_coef);
        return val*lCoefs;
    }

    template <typename T>
    inline gsMatrix<real_t> computeLocValues2 (const T &val, const gsMatrix<index_t> &active, const gsMatrix<real_t> &coefs) const
    {
        const gsMatrix<real_t> lCoefs=getLocCoef(active,coefs);
        return val*lCoefs;
    }

};


/**
 * @brief The gsDistL2 class
 * compute square of the L2 norm
 */
class GISMO_EXPORT gsDistL2 : public gsDistanceOperator
{
public:
    gsDistL2(const gsMatrix<real_t> &coef, const gsFunction<real_t> *func=NULL)
        : gsDistanceOperator(coef, func)
    {}
    virtual void pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t> result
                ) const;
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE;}
    virtual unsigned  geometryNeeds()      const {return NEED_VALUE ;}
};

/**
 * @brief The gsDistH1 class
 * compute square of the H1 seminorm
 */
class GISMO_EXPORT gsDistH1 : public gsDistanceOperator
{
public:
    gsDistH1(const gsMatrix<real_t> &coef, const gsFunction<real_t> *func=NULL)
        : gsDistanceOperator(coef, func)
    {}
    virtual void pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t> result
                ) const;
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE | NEED_GRAD;}
    virtual unsigned  geometryNeeds()      const {return NEED_VALUE ;}
};

/**
 * @brief The gsDistH2 class
 * compute square of the H2 seminorm
 */
class GISMO_EXPORT gsDistH2 : public gsDistanceOperator
{
public:
    gsDistH2(const gsMatrix<real_t> &coef, const gsFunction<real_t> *func=NULL)
        : gsDistanceOperator(coef, func)
    {}
    virtual void pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t> result
                ) const;
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE | NEED_DERIV2;}
    virtual unsigned  geometryNeeds()      const {return NEED_VALUE ;}
};

/**
 * @brief The gsDistW0p class
 * compute p power of the W0p norm
 */
class GISMO_EXPORT gsDistW0p : public gsDistanceOperator
{
protected:
    real_t m_p;
public:
    gsDistW0p(real_t p,const gsMatrix<real_t> &coef, const gsFunction<real_t> *func=NULL)
        : gsDistanceOperator(coef, func), m_p(p)
    {}
    virtual void pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t> result
                ) const;
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE;}
    virtual unsigned  geometryNeeds()      const {return NEED_VALUE ;}
};

/**
 * @brief The gsDistW1p class
 * compute p power of the W1p seminorm
 */
class GISMO_EXPORT gsDistW1p : public gsDistanceOperator
{
protected:
    real_t m_p;
public:
    gsDistW1p(real_t p,const gsMatrix<real_t> &coef, const gsFunction<real_t> *func=NULL)
        : gsDistanceOperator(coef, func), m_p(p)
    {}
    virtual void pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t> result
                ) const;
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE | NEED_GRAD;}
    virtual unsigned  geometryNeeds()      const {return NEED_VALUE ;}
};



class GISMO_DEPRECATED GISMO_EXPORT gsSquaredL2DistOp : public gsBilinearOp<real_t>
{
private:
    gsMatrix<real_t>   &m_u_coef;
    gsFunction<real_t> *m_func;
public:
    gsSquaredL2DistOp (gsMatrix<real_t> &u_coef, gsFunction<real_t> *func)
        : m_u_coef(u_coef), m_func(func)
    {
    }
    virtual void pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t> result
                ) const;
    virtual unsigned  testSpaceNeeds() const {return 0;}
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE;}
    virtual unsigned  geometryNeeds()      const {return NEED_VALUE;}
    virtual void      outputSize (unsigned &r, unsigned &c) const { r=1; c=m_u_coef.cols();}
};


class GISMO_DEPRECATED GISMO_EXPORT gsSquaredH1DistOp : public gsBilinearOp<real_t>
{
private:
    const gsMatrix<real_t>   &m_u_coef;
    const gsFunction<real_t> &m_func;
    const gsFunction<real_t> &m_func_der;

public:
    gsSquaredH1DistOp (const gsMatrix<real_t> &u_coef, const gsFunction<real_t> &func, const gsFunction<real_t> &der)
        : m_u_coef(u_coef), m_func(func), m_func_der(der)
    {
    }
    virtual void pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t> result
                ) const;
    virtual unsigned  testSpaceNeeds() const {return 0;}
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE | NEED_GRAD;}
    virtual unsigned  geometryNeeds()      const {return NEED_VALUE ;}
    virtual void      outputSize (unsigned &r, unsigned &c) const { r=2; c=m_u_coef.cols();}
};

}
