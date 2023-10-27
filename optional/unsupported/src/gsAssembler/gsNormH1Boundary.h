/** @file gsNormH1Boundary.h

    @brief Computes the H1 norm only on the boudary (surface integral).

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Moore
*/

#include <gsAssembler/gsNorm.h>
#pragma once



namespace gismo
{

/** @brief The gsNormH1Boundary class provides the functionality
 * to calculate the H1 - norm on a given boudary (surface integral)
 * between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T>
class gsNormH1Boundary : public gsNorm<T>
{
    friend class gsNorm<T>;
        
public:

    gsNormH1Boundary(const gsField<T> & _field1,
                     const gsFunction<T> & _func2,
                     boxSide s,
                     bool _f2param = false)
    : gsNorm<T>(_field1,_func2), f2param(_f2param), side(s)
    { }
    
    gsNormH1Boundary(const gsField<T> & _field1, boxSide s)
    : gsNorm<T>(_field1), f2param(false), side(s)
    { }

public:
    
    T compute(bool storeElWise = false)
    {
        this->apply(*this,storeElWise, side);
        return this->m_value;
    }


protected:

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        // Setup Quadrature
        const unsigned d = basis.dim();
        const int dir = side.direction();
        gsVector<index_t> numQuadNodes( d );
        for (unsigned i = 0; i < d; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;
    }
    
    // Evaluate on element.
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                         const gsFunction<T>    & _func1,
                         const gsFunction<T>    & _func2,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        _func1.deriv_into(quNodes, f1ders);
        // get the gradients to columns
        f1ders.resize(quNodes.rows(), quNodes.cols() );     
        
        // Evaluate second function (defined of physical domain)
        // Compute geometry related values
        geoEval.evaluateAt(quNodes);    
        _func2.deriv_into(geoEval.values(), f2ders);
        // get the gradients to columns
        f2ders.resize(quNodes.rows(), quNodes.cols() );
    }
    
    // assemble on element
    inline T compute(gsDomainIterator<T>    & element, 
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            //const T d = element.dim();
            // Transform the gradients
            geoEval.transformGradients(k, f1ders, f1pders);

            // Compute the unit normal  
            geoEval.outerNormal(k, side, unormal);
            const T weight = quWeights[k] * unormal.norm();
                 
            // f2ders : N X 1
            sum += weight * (f1pders - f2ders.col(k)).squaredNorm() ;
        }
        accumulated += sum;
        return sum;
    }
    
private:
    
    gsMatrix<T> f1ders, f2ders;
    gsMatrix<T> f1pders;
    gsVector<T> unormal;
    
    bool f2param;

protected:

    boxSide side;

};


} // namespace gismo

