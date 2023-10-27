/** @file gsSpaceTimeSpaceGradSliceNorm.h

    @brief Computes the space-time norm || grad_x w ||^2_{\Sigma_T}

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/
#pragma once

#include<gsAssembler/gsNorm.h>

namespace gismo
{
    /** @brief The gsSpaceTimeDtDeltaKNorm class provides the functionality
     * to calculate the norm || grad_x w ||^2_{\Sigma_T} on the top of the cylinder.
     *
     * \ingroup Assembler
    */

    template <class T>
    class gsSpaceTimeSpaceGradSliceNorm : public gsNorm<T>
    {
        friend class gsNorm<T>;
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;
        typedef gsNorm<T> Base;

    public:

        gsSpaceTimeSpaceGradSliceNorm(const gsField<T> & _field1,
                                const gsFunction<T> & _func2,
                                const std::vector<boundaryInterface>& sliceInterfaces,
                                const std::vector<patchSide>& sides,
                                bool _f2param = false)
                :gsNorm<T>(_field1,_func2), m_sliceInterfaces(sliceInterfaces), m_sides(sides), f2param(_f2param)
        {}

        gsSpaceTimeSpaceGradSliceNorm(const gsField<T> & _field1,
                                const gsFunction<T> & _func2,
                                const gsFunction<T> & _dfunc2,
                                const std::vector<boundaryInterface>& sliceInterfaces,
                                const std::vector<patchSide>& sides,
                                bool _f2param = false)
                : gsNorm<T>(_field1,_func2), m_sliceInterfaces(sliceInterfaces), m_sides(sides), f2param(_f2param)
        {}


        gsSpaceTimeSpaceGradSliceNorm(const gsField<T> & _field1,const std::vector<boundaryInterface>& sliceInterfaces,
                                const std::vector<patchSide>& sides)
                : gsNorm<T>(_field1), m_sliceInterfaces(sliceInterfaces), m_sides(sides), f2param(false)
        {}

    public:

        T compute(bool storeElWise = false)
        {
            m_value = T(0.0);
            side = * m_sides.begin(); // get the side from the vector of m_sides

            for (typename gsMultiPatch<T>::const_biterator bit = patchesPtr->bBegin(); bit != patchesPtr->bEnd(); ++bit){
                if (bit->side() == side) {
                    side = bit->side();
                    this->apply1(*this, storeElWise, bit->patch, side);
                }
            }
            return takeRoot(this->m_value);
        }
        inline T takeRoot(const T v) { return math::sqrt(v);}

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
            evFlags = NEED_MEASURE| NEED_VALUE | NEED_GRAD_TRANSFORM | NEED_JACOBIAN;
        }

        // Evaluate on element.
        inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                             const gsFunction<T>    & _func1,
                             const gsFunction<T>    & _func2,
                             gsMatrix<T>            & quNodes)
        {
            // Compute geometry related values
            geoEval.evaluateAt(quNodes);

            // Evaluate first function
            _func1.deriv_into(quNodes, f1grads);

            // Evaluate second function (defined of physical domain)
            _func2.deriv_into(f2param? quNodes : geoEval.values(), f2grads);

        }

        // assemble on element
        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & accumulated)
        {
            T sum(0.0);
            const int d = element.dim();

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                const T weight = quWeights[k] * unormal.norm();

                if (weight) {
                    geoEval.transformGradients(k, f1grads, ph_f1grads);

                    Rows ph_f1gradsSpace = ph_f1grads.topRows(d - 1);
                    Rows ph_f2gradsSpace = f2grads.topRows(d - 1);

                    geoEval.outerNormal(k, side, unormal);

                    sum += weight * ((ph_f1gradsSpace - ph_f2gradsSpace.col(k)).squaredNorm());
                }
            }
            accumulated += sum;
            return sum;
        }

    private:

        const std::vector<boundaryInterface>& m_sliceInterfaces;
        const std::vector<patchSide>& m_sides;

        using Base::m_value;
        using Base::m_elWise;
        using Base::patchesPtr;

        gsMatrix<T> f1grads, f2grads, ph_f1grads, ph_f2grads;

        gsVector<T> unormal;

        using gsNorm<T>::field1;
        using gsNorm<T>::func2;

        bool f2param;
        boxSide side;
    };
} // namespace gismo
