/** @file gsNormL2M.h

    @brief Computes the L2 norm of M = pI + H'*E(phi)H

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): K. Rafetseder
*/

#include<gsAssembler/gsNorm.h>

#pragma once

namespace gismo
{

template <class T>
class gsNormL2M : public gsNorm<T>
{
    friend  class gsNorm<T>;

public:

    gsNormL2M(const gsField<T> & _field1,
              const gsFunction<T> & _func2,
              bool _f2param = false)
        : gsNorm<T>(_field1,_func2), f2param(_f2param)
    {
        
    }

    gsNormL2M(const gsField<T> & _field1,
              const gsFunction<T> & _func2,
              const gsFunction<T> & _dfunc2,
              bool _f2param = false)
        : gsNorm<T>(_field1,_func2), f2param(_f2param)
    {

    }


    gsNormL2M(const gsField<T> & _field1)
        : gsNorm<T>(_field1), f2param(false)
    { }

public:
    
    T compute(bool storeElWise = false)
    {
        this->apply(*this,storeElWise);
        return m_value;
    }


protected:

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        // Setup Quadrature
        const unsigned d = basis.dim();
        gsVector<index_t> numQuadNodes( d );
        for (unsigned i = 0; i < d; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here
        
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE|NEED_VALUE|NEED_GRAD_TRANSFORM;
    }
    
    // Evaluate on element.
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                  const gsFunction<T>    & _func1,
                  const gsFunction<T>    & _func2,
                  gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        _func1.eval_into(quNodes, f1vals); // matrix with numQuNodes columns and 3 rows for p, phi1, phi2
        _func1.deriv_into(quNodes, f1ders); // matrix with numQuNodes columns and 3*2 rows for the gradient of p, phi1, phi2
        // get the gradients to columns
        //f1ders.resize(3*quNodes.rows(), quNodes.cols() );

        // Evaluate second function (defined on physical domain)
        geoEval.evaluateAt(quNodes);
        if(!f2param)
            _func2.eval_into(geoEval.values(), f2vals);
        else
        {
            _func2.eval_into(quNodes, f2vals); // matrix with numQuNodes columns and 3 rows for p, phi1, phi2
            _func2.deriv_into(quNodes, f2ders); // matrix with numQuNodes columns and 3*2 rows for the gradient of p, phi1, phi2
        }

    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {
        T sum(0.0);

        gsMatrix<T,2,2> Curl, Curl2;
        gsMatrix<T,2,2> symCurl, symCurl2;

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute physical gradients at k as a 2 (dim) x 3 (p, phi1, phi2) matrix
            geoEval.transformGradients(k, f1ders, f1pders);
            if(f2param)
                geoEval.transformGradients(k, f2ders, f2pders);
            
            const T weight = quWeights[k] *  geoEval.measure(k);

            // compute symCurl
            Curl(0,0) = f1pders(1,1);
            Curl(0,1) = -f1pders(0,1);

            Curl(1,0) = f1pders(1,2);
            Curl(1,1) = -f1pders(0,2);

            symCurl = 0.5*(Curl + Curl.transpose());

            // fill  vector M = [M_00, M_01, M_10, M_11] with corresponding matrix entries
            Mvals.resize(4,1);
            Mvals(0,0) = f1vals(0,k) + symCurl (0,0);
            Mvals(1,0) = symCurl (0,1);
            Mvals(2,0) = symCurl (1,0);
            Mvals(3,0) = f1vals(0,k) + symCurl (1,1);

            if(!f2param)
                sum += weight * (Mvals - f2vals.col(k)).squaredNorm();
            else
            {
                // compute symCurl
                //jac2.row(0) = f2pders.col(1).transpose();
                //jac2.row(1) = f2pders.col(2).transpose();

                Curl2(0,0) = f2pders(1,1);
                Curl2(0,1) = -f2pders(0,1);

                Curl2(1,0) = f2pders(1,2);
                Curl2(1,1) = -f2pders(0,2);

                symCurl2 = 0.5*(Curl2 + Curl2.transpose());

                // fill  vector M = [M_00, M_01, M_10, M_11] with corresponding matrix entries
                M2vals.resize(4,1);
                M2vals(0,0) = f2vals(0,k) + symCurl2(0,0);
                M2vals(1,0) = symCurl2 (0,1);
                M2vals(2,0) = symCurl2 (1,0);
                M2vals(3,0) = f2vals(0,k) +symCurl2 (1,1);

                sum += weight * (Mvals - M2vals).squaredNorm();
            }

        }
        accumulated += sum;
        return sum;
    }
    
    inline T takeRoot(const T v) { return math::sqrt(v);}

private:
    // first derivative of func2:
    //const gsFunction<T> * dfunc2; // If this is NULL a numerical approximation will be used

    using gsNorm<T>::m_value;
    using gsNorm<T>::m_elWise;
    //using gsNorm<T>::func2;

    gsMatrix<T> f1vals, f1ders, f2vals, f2ders, f1pders, f2pders, Mvals, M2vals;

    bool f2param;
};


} // namespace gismo

