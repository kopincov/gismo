/** @file gsPDEOperators.hpp

    implementations for pde operators

    Here there is an interface and some implementations that represent
    the act of adding the element-wise contribution of a matrix to
    the global object.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsRecipeAssembler/gsPDEOperators.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometryEvaluator.h>
#include <gsCore/gsFunction.h>

namespace gismo {



template<typename T>
void  gsLaplaceLaplaceOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & unknownSpace,
        const gsPointMapData<T>  & /*geoEval*/,
        gsRecipeAccumulator<T>     result
        ) const
{
    result += testSpace.laplacian().transpose()*unknownSpace.laplacian();
}

template<typename T>
void gsLaplaceOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & unknownSpace,
        const gsPointMapData<T>  & /*geoEval*/,
        gsRecipeAccumulator<T>     result
        ) const
{
    result += testSpace.value().transpose()* unknownSpace.laplacian() ;
}



template<typename T>
void  gsGradGradOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & unknownSpace,
        const gsPointMapData <T> & /*geoEval*/,
        gsRecipeAccumulator<T>     result
        ) const
{
    result += testSpace.deriv().transpose()*unknownSpace.deriv();
}

template<typename T>
void  gsGenericSecondOrderOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & unknownSpace,
        const gsPointMapData <T> & geoEval,
        gsRecipeAccumulator<T>     result
        ) const
{
    gsMatrix<real_t> tmp;
    const index_t tarDim= unknownSpace.info().second;
    const index_t domDim= unknownSpace.info().first;

    if (m_A!=NULL)
    {
        m_A->eval_into(geoEval.value(),tmp);
        result += testSpace.deriv().transpose()*gsAsConstMatrix<real_t>(tmp.data(),tarDim*domDim,tarDim*domDim)*unknownSpace.deriv();
    }
    if (m_b!=NULL)
    {
        m_b->eval_into(geoEval.value(),tmp);
        result += testSpace.deriv().transpose()*gsAsConstMatrix<real_t>(tmp.data(),tarDim*domDim,tarDim)*unknownSpace.value();
    }
    if (m_c!=NULL)
    {
        m_c->eval_into(geoEval.value(),tmp);
        result += testSpace.value().transpose()*gsAsConstMatrix<real_t>(tmp.data(),tarDim,tarDim)*unknownSpace.value();
    }
}


template<typename T>
void gsL2ScalarOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & unknownSpace,
        const gsPointMapData <T> & /*geoEval*/,
        gsRecipeAccumulator<T>     result
        ) const
{
    result +=  testSpace.value().transpose()*unknownSpace.value();
}


template <typename T>
void gsDivergenceOp<T>::pointEval (
            const gsPointFuncData<T> & testSpace,
            const gsPointFuncData<T> & unknownSpace,
            const gsPointMapData <T> & /*geoEval*/,
            gsRecipeAccumulator<T>     result
            ) const
    { result += testSpace.value().transpose()*unknownSpace.div(); }



template<typename T>
void gsGradientOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & unknownSpace,
        const gsPointMapData <T> & /*geoEval*/,
        gsRecipeAccumulator<T>     result
        ) const
{
    result += testSpace.value().transpose()*unknownSpace.deriv();
}



template<typename T>
void  gsL2TestOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & /*unknownSpace*/,
        const gsPointMapData <T> & geoEval,
        gsRecipeAccumulator<T>     result
        ) const
{
    gsMatrix<T> temp, temp2;
    testFunc->eval_into(geoEval.value(), temp);
    //TODO use function sets
    const int blockSize= testSpace.info().second;
    const int numTest=temp.rows()/blockSize;
    temp2.resize(testSpace.value().cols(),numTest);
    for (int i=0; i< numTest;++i)
        temp2.col(i)=testSpace.value().transpose()*temp.middleRows(blockSize*i,blockSize);
    result+=temp2;
}



template<typename T>
void  gsL2TestVecOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & unknownSpace,
        const gsPointMapData <T> & geoEval,
        gsRecipeAccumulator<T>     result
        ) const
{
    static bool warned=false;
    if (!warned)
    {
        gsWarn<<"gsL2TestVecOp is going to be DEPRECATED, use gsL2TestOp\n";
        warned=true;
    }
    gsL2TestOp<T> fake(*testFunc);
    fake.pointEval(testSpace,unknownSpace,geoEval, result);
}

template<typename T>
void  gsBoundaryL2TestOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & /*unknownSpace*/,
        const gsPointMapData <T> & geoEval,
        gsRecipeAccumulator<T>     result
        ) const
{
    testFunc->eval_into(geoEval.value(),temp);
    temp.resize(testSpace.info().second,temp.size()/testSpace.info().second);
    result += testSpace.value().transpose()*temp;
}



template<typename T>
void gsBoundaryL2TestVecOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & unknownSpace,
        const gsPointMapData <T> & geoEval,
        gsRecipeAccumulator<T>     result
        ) const
{
    static bool warned=false;
    if (!warned)
    {
        gsWarn<<"gsBoundaryL2TestVecOp is going to be DEPRECATED, use gsBoundaryL2TestOp\n";
        warned=true;
    }
    gsBoundaryL2TestOp<T> fake(*testFunc);
    fake.pointEval(testSpace,unknownSpace,geoEval,result);
}


template<typename T>
void gsBoundaryL2ScalarOp<T>::pointEval (
        const gsPointFuncData<T> & testSpace,
        const gsPointFuncData<T> & unknownSpace,
        const gsPointMapData <T> & /*geoEval*/,
        gsRecipeAccumulator<T>     result
        ) const
{
    result += testSpace.value().transpose()*unknownSpace.value();
}



template<typename T>
void gsBoundaryNormalDerValueOp<T>::pointEval (
        const gsPointFuncData<T>  &testSpace,
        const gsPointFuncData<T>  &unknownSpace,
        const gsPointMapData <T>  &geoEval,
        gsRecipeAccumulator<T>  result
        ) const
{
    gsMatrix<T> tmp = (unknownSpace.jacobians()*geoEval.outNormal());
    tmp.resize(unknownSpace.info().second, unknownSpace.deriv().cols());
    result += testSpace.value().transpose()*tmp;
}



template<typename T>
void gsBoundaryNormalDerNormalDerOp<T>::pointEval (
        const gsPointFuncData<T>  &testSpace,
        const gsPointFuncData<T>  &unknownSpace,
        const gsPointMapData <T>  &geoEval,
        gsRecipeAccumulator<T>  result
        ) const
{
    typename gsPointMapData<T>::constColumn oN = geoEval.outNormal();
    if (testSpace.info().second==1)
        result += ((testSpace.jacobians()*oN)*(unknownSpace.jacobians()*oN).transpose());
    else
    {
        // the for loop avoids temporaries
        // unluckily it is not elegant
        const index_t u_activeNum=unknownSpace.jacobians().rows()/unknownSpace.info().second;
        const index_t t_activeNum=testSpace.jacobians().rows()/testSpace.info().second;

        for (index_t u_a=0; u_a<u_activeNum;++u_a)
            for (index_t t_a=0; t_a<t_activeNum;++t_a)
            {
            result(t_a,u_a)+= (testSpace.jacobian(t_a)*oN).dot((unknownSpace.jacobian(u_a)*oN));
            }
    }
}


template<typename T>
void gsBoundaryNormalDerTestOp<T>::pointEval (
        const gsPointFuncData<T>  &testSpace,
        const gsPointFuncData<T>  &/*unknownSpace*/,
        const gsPointMapData <T>  &geoEval,
        gsRecipeAccumulator<T>  result
        ) const
{
    testFunc->eval_into(geoEval.value(),temp);
    gsMatrix<T> tmp=(testSpace.jacobians()*geoEval.outNormal());
    const index_t tarD=testSpace.info().second;
    // TODO use function sets and destroy this hack
    if (testSpace.info().second==1)
        result += tmp*temp.transpose();
    else
    {
        tmp.resize(tarD,tmp.size()/tarD);
        result += tmp.transpose()*temp;
    }
}

template<typename T>
void gsBoundaryNormalDerTestNormalOp<T>::pointEval (
        const gsPointFuncData<T>  &testSpace,
        const gsPointFuncData<T>  &/*unknownSpace*/,
        const gsPointMapData <T>  &geoEval,
        gsRecipeAccumulator<T>  result
        ) const
{
    testFunc->eval_into(geoEval.value(),temp);
    gsMatrix<T> tmp=(testSpace.jacobians()*geoEval.outNormal());
    GISMO_ASSERT(testSpace.info().second == 1, "Not implementes for target dimension higher then one!");
    // TODO use function sets and destroy this hack
    result += tmp*(geoEval.outNormal().transpose()*temp).transpose();
}

template<typename T>
void  gsBoundaryNormalDerTestVecOp<T>::pointEval (
        const gsPointFuncData<T>  &testSpace,
        const gsPointFuncData<T>  &/*unknownSpace*/,
        const gsPointMapData <T>  &geoEval,
        gsRecipeAccumulator<T>  result
        ) const
{
    const index_t tarD=testSpace.info().second;
    testFunc->eval_into(geoEval.value(),temp);
    gsMatrix<T> tmp=testSpace.jacobians()*geoEval.outNormal();
    tmp.resize(tarD,tmp.size()/tarD);
    result += (temp.transpose()*tmp).transpose();
}

template<typename T>
void  gsLinElastOp<T>::pointEval (
        const gsPointFuncData<T>  &testSpace,
        const gsPointFuncData<T>  &unknownSpace,
        const gsPointMapData <T>  &/*geoEval*/,
        gsRecipeAccumulator<T>  result
        ) const
{
    index_t au=unknownSpace.value().cols();
    index_t at=testSpace.value().cols();
    for (index_t i=0; i<at;++i)
        for (index_t j=0; j<au;++j)
        {
            result(i,j)+=
                    (
                        m_mu*
                        (
                            testSpace.jacobian(i).array()
                            *(unknownSpace.jacobian(j)+unknownSpace.jacobian(j).transpose()).array()
                            ).sum()
                        +m_lambda*  ( testSpace.div().col(i).transpose()*unknownSpace.div().col(j) ).value()
                        );
        }
}


} // namespace gismo

