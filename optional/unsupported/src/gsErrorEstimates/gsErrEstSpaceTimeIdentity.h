/** @file gsSpaceTimeErrorIdentity.h

    @brief Computes the main contribution of the space-time error identity ||f + Delta_x v - v_t||^2_Q$,
    where the full error identity is defined as

    EI = || grad_x (v - u_0)||^2_{\Sigma_0} ||f + Delta_x v - v_t||^2_Q.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Matculevich
*/

#include<gsAssembler/gsNorm.h>

#pragma once

namespace gismo
{

/** @brief The gsSpaceTimeErrorIdentity class provides the functionality
 * the main contribution of the space-time error identity
 * ||f + Delta_x v - v_t||^2_Q$.
 *
 * \ingroup Assembler
*/
    template <class T>
    class gsErrEstSpaceTimeIdentity : public gsNorm<T>
    {
        friend  class gsNorm<T>;
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;

    public:

        gsErrEstSpaceTimeIdentity(const gsField<T> & _approx,
                                  const gsFunction<T> & _rhs,
                                  bool _f2param = false)
                : gsNorm<T>(_approx,_rhs), dfunc2(NULL), f2param(_f2param)
        {

        }
        

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
            evFlags = NEED_MEASURE| NEED_VALUE| NEED_JACOBIAN | NEED_2ND_DER | NEED_GRAD_TRANSFORM;
        }

        // Evaluate on element.
        void evaluate(gsGeometryEvaluator<T> & geoEval,
                      const gsFunction<T>    & approx,
                      const gsFunction<T>    & rhs,
                      gsMatrix<T>            & quNodes)
        {
            // Evaluate first function
            approx.deriv_into(quNodes, approxGrad);
            approx.deriv2_into(quNodes, approxDeriv2);

            // Evaluate second function 
            geoEval.evaluateAt(quNodes);

            rhs.eval_into( geoEval.values(), rhsVals );

        }

        // assemble on element
        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & accumulated)
        {
            T sum(0.0);
            const unsigned d = element.dim();

            // An auxiliary eye matrix to calculate the laplacian
            gsMatrix<> ones = gsMatrix<>::Identity(d-1, 1);
            ones.setOnes();

            gsMatrix<T> approxSpaceLaplace;


            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                const T weight = quWeights[k] *  geoEval.measure(k);

                if (weight) {
                    // Transform the gradients
                    // ph_approxGrad : N X dim
                    // ph_approxGradSpace, ph_approxGradTime : N X (dim-1)
                    geoEval.transformGradients(k, approxGrad, ph_approxGrad);
                    geoEval.transformDeriv2Hgrad(k, approxGrad, approxDeriv2, ph_approxDeriv2);
                    ph_approxSpaceLaplace = ph_approxDeriv2.leftCols(d - 1);

                    Rows ph_approxGradTime = ph_approxGrad.bottomRows(1);
                    approxSpaceLaplace = ph_approxSpaceLaplace * ones;

                    sum += weight * math::pow(rhsVals(0, k) + approxSpaceLaplace(0, 0) - ph_approxGradTime(0, 0), 2.0);
                }
            }

            accumulated += sum;
            return sum;
        }


        inline T takeRoot(const T v) { return math::sqrt(v);}

    private:
        const gsFunction<T> * dfunc2; // If this is NULL a numerical approximation will be used

        using gsNorm<T>::m_value;
        using gsNorm<T>::m_elWise;

        gsMatrix<T> approxGrad, approxDeriv2, rhsVals;
        gsMatrix<T> ph_approxGrad, ph_approxSpaceLaplace, ph_approxDeriv2;
        gsVector<T> unormal;

        bool f2param;// not used yet
    };
} // namespace gismo
