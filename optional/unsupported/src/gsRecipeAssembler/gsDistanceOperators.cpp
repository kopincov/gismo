/** @file gsDistanceOperators.cpp

    implementations for distance operators

    Here there is an interface and some implementations that represent
    the act of adding the element-wise contribution of a matrix to
    the global object.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsRecipeAssembler/gsDistanceOperators.h>
#include <gsRecipeAssembler/gsRecipe.h>
#include <gsCore/gsFunction.h>

namespace gismo {



void  gsDistL2::pointEval (
                const gsPointFuncData<real_t> & /*testSpace*/,
                const gsPointFuncData<real_t> & unknownSpace,
                const gsPointMapData <real_t> & geoEval,
                gsRecipeAccumulator<real_t>     result
        ) const
{
    gsMatrix<real_t> valS=computeLocValues2(unknownSpace.value(),unknownSpace.actives());
    gsMatrix<real_t> valF;
    if (m_func)
        m_func->eval_into(geoEval.value(),valF);
    else
        valF.setZero(valS.rows(),valS.cols());
    // hack for multiple scalar equations at once
    if(valS.rows()==1 && valF.rows()==valS.cols())
        result += (valS-valF.transpose()).colwise().squaredNorm();
    else
        result += (valS-valF).colwise().squaredNorm();
}



void  gsDistH1::pointEval (
                const gsPointFuncData<real_t> & /*testSpace*/,
                const gsPointFuncData<real_t> & unknownSpace,
                const gsPointMapData <real_t> & geoEval,
                gsRecipeAccumulator<real_t>     result
        ) const
{
    gsMatrix<real_t> valS=computeLocValues2(unknownSpace.deriv(),unknownSpace.actives());
    gsMatrix<real_t> valF;
    valF.setZero(valS.rows(),valS.cols());
    if (m_func)
    {
        m_func->deriv_into(geoEval.value(),valF);
        valF.conservativeResize(valS.rows(),valS.cols());
    }
    result += (valS-valF).colwise().squaredNorm();
}



void  gsDistH2::pointEval (
                const gsPointFuncData<real_t> & /*testSpace*/,
                const gsPointFuncData<real_t> & unknownSpace,
                const gsPointMapData <real_t> & geoEval,
                gsRecipeAccumulator<real_t>     result
        ) const
{
    gsMatrix<real_t> valS=computeLocValues2(unknownSpace.deriv2(),unknownSpace.actives());
    gsMatrix<real_t> valF;
    valF.setZero(valS.rows(),valS.cols());
    if (m_func)
    {
        m_func->deriv2_into(geoEval.value(),valF);
        valF.conservativeResize(valS.rows(),valS.cols());
    }
    result += (valS-valF).colwise().squaredNorm();
}



void  gsDistW0p::pointEval (
                const gsPointFuncData<real_t> & /*testSpace*/,
                const gsPointFuncData<real_t> & unknownSpace,
                const gsPointMapData <real_t> & geoEval,
                gsRecipeAccumulator<real_t>     result
        ) const
{
    gsMatrix<real_t> valS=computeLocValues2(unknownSpace.value(),unknownSpace.actives());
    gsMatrix<real_t> valF;
    if (m_func)
        m_func->eval_into(geoEval.value(),valF);
    else
        valF.setZero(valS.rows(),valS.cols());
    valS-=valF;
    for (index_t r=0; r<valS.rows();++r)
        for (index_t c=0;c<valS.cols();++c)
            valS(r,c)=math::pow(valS(r,c),m_p);
    // hack for multiple scalar equations at once
    if(valS.rows()==1 && valF.rows()==valS.cols())
        result += (valS-valF.transpose()).colwise().squaredNorm();
    else
        result += valS.colwise().sum();
}


void  gsDistW1p::pointEval (
                const gsPointFuncData<real_t> & /*testSpace*/,
                const gsPointFuncData<real_t> & unknownSpace,
                const gsPointMapData <real_t> & geoEval,
                gsRecipeAccumulator<real_t>     result
        ) const
{
    gsMatrix<real_t> valS=computeLocValues2(unknownSpace.deriv(),unknownSpace.actives());
    gsMatrix<real_t> valF;
    if (m_func)
        m_func->deriv_into(geoEval.value(),valF);
    else
        valF.setZero(valS.rows(),valS.cols());
    valS-=valF;
    for (index_t r=0; r<valS.rows();++r)
        for (index_t c=0;c<valS.cols();++c)
            valS(r,c)=math::pow(valS(r,c),m_p);
    result += valS.colwise().sum();
}



void gsSquaredL2DistOp::pointEval (
                const gsPointFuncData<real_t> & /*testSpace*/,
                const gsPointFuncData<real_t> & unknownSpace,
                const gsPointMapData <real_t> & geoEval,
                gsRecipeAccumulator<real_t>     result
        ) const
{
    const index_t au = unknownSpace.actives().rows();
    const index_t tdim = unknownSpace.info().second;
    gsMatrix<real_t> value(tdim,m_u_coef.cols());
    gsMatrix<real_t> value_F(tdim,m_u_coef.cols());
    if (m_func)
        m_func->eval_into(geoEval.value(),value_F);
    else value_F.setZero(tdim,m_u_coef.cols());

    value.setZero(tdim,m_u_coef.cols());
    for (index_t i = 0; i < au; ++i)
        value+=unknownSpace.value().middleRows(i*tdim,tdim)*m_u_coef.row(unknownSpace.actives()(i));
    result += (value-value_F).colwise().squaredNorm();
}


void gsSquaredH1DistOp::pointEval (
                const gsPointFuncData<real_t> & /*testSpace*/,
                const gsPointFuncData<real_t> & unknownSpace,
                const gsPointMapData <real_t> & geoEval,
                gsRecipeAccumulator<real_t>     result
        ) const
{
    const index_t au = unknownSpace.actives().rows();
    index_t tdim = unknownSpace.info().second;
    gsMatrix<real_t> value(tdim,m_u_coef.cols());
    value.setZero(tdim,m_u_coef.cols());
    for (index_t i = 0; i < au; ++i)
        value+=unknownSpace.value().middleRows(i*tdim,tdim)*m_u_coef.row(unknownSpace.actives()(i));
    value-=m_func.eval(geoEval.value());
    result.row(0) += value.colwise().squaredNorm();

    tdim = unknownSpace.info().second*unknownSpace.info().first;
    value.resize(tdim,m_u_coef.cols());
    value.setZero(tdim,m_u_coef.cols());
    for (index_t i = 0; i < au; ++i)
        value+=unknownSpace.deriv().middleRows(i*tdim,tdim)*m_u_coef.row(unknownSpace.actives()(i));
    value-=m_func_der.eval(geoEval.value());
    result.row(1) += value.colwise().squaredNorm();
}

} // namespace gismo
