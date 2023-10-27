/** @file gsDgNorm.h

    @brief Computes the H1 norm.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Moore
*/

// Under progresss

#include<gsAssembler/gsNorm.h>

#pragma once

namespace gismo
{

/** @brief The gsDGnorm class provides the functionality
 * to calculate the DG-norm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T>
class gsDgNorm : public gsNorm<T>
{
    friend  class gsNorm<T>;

public:

    gsDgNorm(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             bool _f2param = false) 
    : gsNorm<T>(_field1,_func2), dfunc2(NULL), f2param(_f2param)
    { 
        
    }

    gsDgNorm(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             const gsFunction<T> & _dfunc2,
             bool _f2param = false)
    : gsNorm<T>(_field1,_func2), dfunc2(&_dfunc2), f2param(_f2param)
    {

    }


    gsDgNorm(const gsField<T> & _field1) 
    : gsNorm<T>(_field1), dfunc2(NULL), f2param(false)
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
    void evaluate(gsGeometryEvaluator<T> & geoEval1,
                  gsGeometryEvaluator<T> & geoEval2,
                  const gsGeometry<T>    & func11,
                  const gsGeometry<T>    & func12,
                  const gsFunction<T>    & func21,
                  const gsFunction<T>    & func22,
                  const boundaryInterface & bi, // interface
                  const T mu,
                  gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        func11.deriv_into(quNodes, f11ders);
        func12.deriv_into(quNodes, f12ders);

        // Evaluate second function (defined of physical domain)
        geoEval1.evaluateAt(quNodes);
        geoEval2.evaluateAt(quNodes);
        if(dfunc2==NULL)
        {
            func21.deriv_into(geoEval1.values(), f21ders);
        }
        else
            dfunc2->eval_into(geoEval1.values(), f21ders);

        // ** Evaluate function v
        //gsMatrix<T> f2val = func2Param ? func2.deriv(quNodes)
        //: func2.eval( geoEval->values() );      
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element1, 
                     gsDomainIterator<T>    & element2,
                     gsGeometryEvaluator<T> & geoEval1,
                     gsGeometryEvaluator<T> & geoEval2,
                     gsVector<T> const      & quWeights)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Transform the gradients
            geoEval1.transformGradients(k, f1ders, f1pders);
            //if ( f2Param )
            //f2ders.col(k)=geoEval.gradTransforms().block(0, k*d,d,d) * f2ders.col(k);// to do: generalize
            
            const T weight = quWeights[k] *  geoEval.measure(k);
            sum += weight * (f1pders - f2ders.col(k)).squaredNorm();
        }
        return sum;
    }
    
private:
    // first derivative of func2:
    const gsFunction<T> * dfunc2; // If this is NULL a numerical approximation will be used

    using gsNorm<T>::m_value;
    using gsNorm<T>::m_elWise;

    gsMatrix<T> f1ders, f2ders;
    gsMatrix<T> f1pders;

    bool f2param;// not used yet
};


} // namespace gismo

