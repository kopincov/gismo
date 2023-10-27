/** @file gsQuasiInterpolate.hpp

    @brief Different Quasi-Interpolation Schemes based on the article
    "Spline methods (Lyche Morken)"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Haberleitner
*/


#pragma once

#include<gsCore/gsLinearAlgebra.h>
#include<gsCore/gsFunction.h>
#include<gsNurbs/gsBSpline.h>
#include<gsUtils/gsCombinatorics.h>

namespace gismo {

/*
template<typename T>
void qiCwiseData(const gsTensorBSplineBasis<T,2> & tbsp,
                 const gsVector<index_t,2> & ind,
                 std::vector<gsMatrix<T> > & qiNodes,
                 std::vector<gsMatrix<T> > & qiWeights)
{

}

template<typename T>
void gsQuasiInterpolate<T>::
void compute(const gsTensorBSplineBasis<T,2> & tbsp,
             const gsFunction<T> & fun,
             gsTensorBSpline<T,2> & res)
{
    gsMatrix<T> fval, coefs(tbsp.size(), fun.targetDim());
    
    std::vector<gsMatrix<T> > qiNodes(2), qiWeights(2);

    gsGridIterator<index_t,CUBE> it( tbsp.sizeCwise() );
    for(index_t c = 0; it; ++it, ++c) // For all coefficients
    {// todo: split into boundary CPs and interior CPs: rule is uniform in the interior
        
        // Get the QI data for the current coefficient
        qiCwiseData(tbsp, *it, quNodes, quWeights); // (border cases may be less points)

        // Apply tensor-product quasi-interpolation
        gsGridIterator<T,CWISE> pt(quNodes);        
        coefs.row(c).setZero();
        for(; pt; ++pt)
        {
            fun.eval_into(*pt, fval);
            const T wgt = qiWeights[0].at(pt.index(0)) * qiWeights[1].at(pt.index(1));
            coefs.row(c) += wgt * fval.transpose();
        }
    }
    res = gsTensorBSpline<T,2>(tbsp,give(coefs));
}
//*/

template<typename T>
void gsQuasiInterpolate<T>::Taylor(const gsBSplineBasis<T> &b, const gsFunction<T> &fun, const int &r, gsBSpline<T> &result)
{
    const gsKnotVector<T> & kv = b.knots();
    int deg = b.degree();
    gsMatrix<T> xj = b.anchors();

    int n = xj.size();
    int dim = fun.targetDim();
    gsMatrix<T> coefs(n,dim);

    std::vector<gsMatrix<T> > derivs;
    fun.evalAllDers_into(xj, r, derivs);


    gsMatrix<T> val;
    std::vector<T> knots;
    for(int j=0; j<n; j++)
    {
        val.setZero(1,dim);
        knots.clear();
        for(int i=j+1; i<=j+deg; i++)
            knots.push_back(kv[i]);

        for(int k=0; k<=r; k++) // (r+1) nodes
        {
            const T factor1 = derivProd(knots, deg-k, xj(j)); //coeff
            for(int i=0; i<dim; i++)
            {
                const T factor2 = derivs[k](i,j); //node
                val(i) += std::pow(-1.0,k) * factor1 * factor2;
            }
        }
        val /= factorial(deg);
        coefs.row(j) = val;
    }
    result = gsBSpline<T>(b, give(coefs));
}



template<typename T>
void gsQuasiInterpolate<T>::Schoenberg(const gsBSplineBasis<T> &b, const gsFunction<T> &fun, gsBSpline<T> &result)
{
    const gsMatrix<T> xj = b.anchors(); // 1 node per CP
    // coef == 1

    gsMatrix<T> coefs;
    fun.eval_into(xj, coefs);
    //result is is d x n but coefficients are n x d
    coefs.transposeInPlace();

    result = gsBSpline<T>(b, give(coefs));
}



template<typename T>
void gsQuasiInterpolate<T>::EvalBased(const gsBSplineBasis<T> &b, const gsFunction<T> &fun, const bool specialCase, gsBSpline<T> &result)
{
    const gsKnotVector<T> & kv = b.knots();
    const int n = b.size();

    gsMatrix<T> coefs(n, fun.targetDim());

    gsMatrix<T> knots(1,kv.size());
    for(unsigned int i=0; i<kv.size(); i++)
        knots(i) = kv[i];

    gsMatrix<T> TmpCoefs;

    int type = 0;
    if(specialCase)
    {
        type = b.degree();
        GISMO_ASSERT( (type == 1 || type == 2 || type == 3),
                      "quasiInterpolateEvalBased is implemented for special cases of deg 1, 2 or 3!");
    }

    switch(type)
    {
    case 1: //piecewise linear (section 8.2.1 of "Spline methods (Lyche Morken)")
    {
        fun.eval_into(knots, TmpCoefs);
        for(int i=0; i<n; i++)
            coefs.row(i) = TmpCoefs.col(i+1).transpose();
        break;
    }
    case 2: //3-point quadratic (section 8.2.2 of "Spline methods (Lyche Morken)")
    {
        fun.eval_into(knots, TmpCoefs);
        gsMatrix<T> knotsAvg(1, kv.size()-1);
        for(unsigned int i=0; i<kv.size()-1; i++)
            knotsAvg(i) = (kv[i]+kv[i+1])/2;
        gsMatrix<T> TmpCoefsAvg;
        fun.eval_into(knotsAvg, TmpCoefsAvg);

        coefs.row(0) = TmpCoefs.col(0);
        for(int i=1; i<n-1; i++)
        {
            // formula: (-a + 4b - c)/2;
            coefs.row(i).noalias() = 
                ( - TmpCoefs.col(i+1)
                  + 4 * TmpCoefsAvg.col(i+1)
                  - TmpCoefs.col(i+2) ) / 2;
        }
        coefs.row(n-1) = TmpCoefs.col(n);
        break;
    }
    case 3: //5-point cubic (section 8.2.3 of "Spline methods (Lyche Morken)")
    {
        fun.eval_into(knots, TmpCoefs);
        gsMatrix<T> knotsAvg(1, kv.size()-1);
        for(unsigned int i=0; i<kv.size()-1; i++)
            knotsAvg(i) = (kv[i]+kv[i+1])/2;
        gsMatrix<T> TmpCoefsAvg;
        fun.eval_into(knotsAvg, TmpCoefsAvg);

        coefs.row(0) = TmpCoefs.col(3).transpose();

        // formula: (- 5a + 40b - 24c + 8d - e)/18
        coefs.row(1).noalias() =
                ( -  5 * TmpCoefs.col(3)
                  + 40 * TmpCoefsAvg.col(3)
                  - 24 * TmpCoefs.col(4)
                  +  8 * TmpCoefsAvg.col(4)
                  - TmpCoefs.col(5) ) / 18;

        for(int i=2; i<n-2; i++)
        {
            // formula: (a - 8b + 20c - 8d + e)/6
            coefs.row(i).noalias() =
                    ( TmpCoefs.col(i+1)
                      -  8 * TmpCoefsAvg.col(i+1)
                      + 20 * TmpCoefs.col(i+2)
                      -  8 * TmpCoefsAvg.col(i+2)
                      + TmpCoefs.col(i+3) ) / 6;
        }

        // formula: (- a + 8b - 24c + 40d - 5e)/18
        coefs.row(n-2).noalias() =
                ( - TmpCoefs.col(n-2)
                  +  8 * TmpCoefsAvg.col(n-2)
                  - 24 * TmpCoefs.col(n-1)
                  + 40 * TmpCoefsAvg.col(n-1)
                  -  5 * TmpCoefs.col(n) ) / 18;

        coefs.row(n-1) = TmpCoefs.col(n);

        break;
    }
    default: //if none of the special cases, (Theorem 8.7 and Lemma 9.7 of "Spline methods (Lyche Morken)")
        gsMatrix<T> xik, weights;
        for(int i=0; i<n; i++)
        {
            //look for the greatest subinterval to chose the interpolation points from
            int gsi = greatestSubInterval(kv, i, i+kv.degree());

            //compute equally distributed points in greatest subinterval
            distributePoints(kv[gsi], kv[gsi+1], kv.degree()+1, xik);

            //compute the factors omega_{ik}
            computeWeights(xik, kv, i+1, weights);

            //compute the coefficients lambda_i of the quasi-interpolant
            coefs.row(i) = computeControlPoints(weights, fun, xik);
        }
        break;
    }
    //return the quasi-interpolant
    result = gsBSpline<T>(b, give(coefs));
}


template<typename T>
T gsQuasiInterpolate<T>::derivProd(const std::vector<T> &zeros, const int &order, const T &x)
{
    if(order == 0) // value
        return (x - gsAsConstMatrix<T,1>(zeros).array()).prod();

    std::vector<T> tmpZeros;
    
    if(order == 1) // first derivative
    {
        const index_t n1 = zeros.size() - 1;
        tmpZeros = zeros;
        T val = 0;
        for(typename std::vector<T>::iterator it = tmpZeros.begin(); it!=tmpZeros.end(); ++it)
        {
            std::iter_swap(it, tmpZeros.end()-1);
            val += (x - gsAsConstMatrix<T,1>(tmpZeros,1,n1).array()).prod(); // eval product
            std::iter_swap(it, tmpZeros.end()-1);
        }
        return val;
    }

    // Reccursion for higher order derivatives
    const int n = zeros.size();
    T val = 0;
    for(int i=0; i!=n; i++)
    {
        tmpZeros = zeros;
        tmpZeros.erase(tmpZeros.begin() + i);
        val += derivProd(tmpZeros, order-1, x);
    }
    return val;
}


template<typename T>
void gsQuasiInterpolate<T>::distributePoints(T a, T b, int n, gsMatrix<T> &points)
{
    points.resize(1,n);
    for(int k=0; k<n; k++)
        points.at(k) = a + (T)k/(n-1) * (b-a);
}


template<typename T>
void gsQuasiInterpolate<T>::computeWeights(const gsMatrix<T> &points, const gsKnotVector<T> &knots, const int &pos, gsMatrix<T> &weights)
{
    const int deg = knots.degree();
    weights.resize(1,deg+1);

    gsMatrix<T> pointsReduced(1,deg);
    gsVector<int> indices(deg);
    for(int k=0; k<deg+1; k++)
    {
        T constant = (T)factorial(deg);
        
        //get the list of points without 'points(k)' and
        //multiply the values of the numerator together, to get the total constant
        for(int i=0; i<k; i++)
        {
            pointsReduced(i) = points(i);
            constant *= points(k) - points(i);
        }
        for(int i=k+1; i<deg+1; i++)
        {
            pointsReduced(i-1) = points(i);
            constant *= points(k) - points(i);
        }

        //get all permutations of the list
        for(int i=0; i<deg; i++)
            indices[i]=i;

        T sum = 0;
        do
        {
            //compute the product of values for this permutation
            T factor = 1;
            for(int i=0; i<deg; i++)
            {
                factor *= (knots[pos+indices[i]] - pointsReduced(i));
            }
            sum += factor; //sum up the products for each permutation
        }
        while(std::next_permutation(indices.data(), indices.data()+deg));

        weights(k) = sum / constant;
    }
}


template<typename T>
gsMatrix<T> gsQuasiInterpolate<T>::computeControlPoints(const gsMatrix<T> &weights, const gsFunction<T> &fun, const gsMatrix<T> &xik)
{
    gsMatrix<T> funValues;
    fun.eval_into(xik, funValues);
    return (weights.asDiagonal() * funValues.transpose()).colwise().sum();
}


template<typename T>
int gsQuasiInterpolate<T>::greatestSubInterval(const gsKnotVector<T> &knots, const int &posStart, const int &posEnd)       //ToDo: move to gsKnotVector
{
    const int diff = posEnd-posStart;
    T maxDist=0.0;
    int maxInd = posStart;
    for(int i=1; i<diff+1; i++)
    {
        const T dist = knots[posStart+i+1] - knots[posStart+i];
        if(dist > maxDist)
        {
            maxDist = dist;
            maxInd = posStart+i;
        }
    }
    return maxInd;
}


} // gismo
