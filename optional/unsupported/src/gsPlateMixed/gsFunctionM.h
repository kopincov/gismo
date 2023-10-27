/** @file gsFunctionM.h

    @brief Function to compute the bending moments M

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/

#pragma once

#include <gsCore/gsFunction.h>

namespace gismo
{

template<class T>
class gsFunctionM : public gsFunction<T>
{
public:
    typedef gsFunction<T> Base;

    /// Shared pointer for gsFunctionM
    typedef memory::shared_ptr< gsFunctionM > Ptr;

    /// Unique pointer for gsFunctionM
    typedef memory::unique_ptr< gsFunctionM > uPtr;


    using Base::support;
    using Base::domainDim;
    using Base::targetDim;
public:
    gsFunctionM(const gsGeometry<T>& geometry,const gsGeometry<T>& p, const gsGeometry<T>& phi):
        m_geometry(geometry), m_p(p), m_phi(phi)
    {
        // Voigt notation to Matrix
        VoigtVtoM.col(0) << 0,1,0;
        VoigtVtoM.col(1) << 0,1,1;

        // Matrix to Voigt notation
        VoigtMtoV [0][0]=0;
        VoigtMtoV [1][1]=1;
        VoigtMtoV [0][1]=2;
        VoigtMtoV [1][0]=2;

    }


    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsFunctionM, clone)

    /*
      Member functions with non-virtual implementations
      (override the _into versions in derived classes).
    */

    /** \brief Evaluate the function at points \a u into \a result.
     *
     * Let \em n be the dimension of the source space ( n = domainDim() ).\n
     * Let \em m be the dimension of the image/target space ( m = targetDim() ).\n
     * Let \em N denote the number of evaluation points.
     *
     * \param[in] u gsMatrix of size <em>n</em> x <em>N</em>, where each
     * column of \em u represents one evaluation point.
     * \param[out] result gsMatrix of size <em>m</em> x <em>N</em>, where each
     * column of \em u represents the result of the function at the
     * respective valuation point.
     */
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m.setZero(3, u.cols());
        mHat.setZero(3, u.cols());

        m_p.eval_into(u, pVals);
        m_phi.deriv_into(u, phiGrads);

        m.row(0) = pVals + phiGrads.row(1);
        m.row(1) = pVals - phiGrads.row(2);
        m.row(2) = 0.5* (-phiGrads.row(0) + phiGrads.row(3));


        unsigned evFlags =NEED_JACOBIAN | NEED_MEASURE | NEED_NORMAL;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_geometry));
        geoEval->evaluateAt(u);


        for(index_t N=0; N!=u.cols();++N)
        {
            // definition of M
            m.col(N) /= geoEval->measure(N);

            // transform into local cartesian coordinate system

            // geometry quantities
            geoEval->normal(N,normal);
            normal.normalize();
            F.leftCols(2) = geoEval->jacobian(N);
            F.col(2)      = normal;
            FInvT = F.inverse().transpose();

            //localCartCoord.col(0) = FInvT.col(0);
            //localCartCoord.col(1) = F.col(1);
            localCartCoord.col(0) = F.col(0);
            localCartCoord.col(1) = FInvT.col(1);
            localCartCoord.col(0).normalize();
            localCartCoord.col(1).normalize();

            // loop over components M_11, M_22, M_12
            for(index_t comp =0; comp!=3; ++comp)
            {
                index_t alpha = VoigtVtoM(comp,0);
                index_t beta = VoigtVtoM(comp,1);

                for(index_t sigma = 0; sigma!=2; ++sigma)
                    for(index_t tau = 0; tau!=2; ++tau)
                        mHat(comp,N) += m(VoigtMtoV[sigma][tau],N)* localCartCoord.col(alpha).dot(F.col(sigma)) * localCartCoord.col(beta).dot(F.col(tau));

            }
        }

        result = mHat;
    }

    virtual short_t domainDim () const
    {
        return 2;
    }

private:

    const gsGeometry<T>& m_geometry;
    const gsGeometry<T>& m_p;
    const gsGeometry<T>& m_phi;

    mutable gsVector<T> normal; // Normal to the shell centerline
    mutable gsMatrix<T,3,3> F; // deformation gradient (a_1, a_2, a_3)
    mutable gsMatrix<T,3,3> FInvT; // (a^1, a^2, a^3)
    mutable gsMatrix<T,3,3> localCartCoord; // (e^1, e^2, e^3)

    mutable gsMatrix<> pVals;
    mutable gsMatrix<> phiGrads;
    mutable gsMatrix<> m;
    mutable gsMatrix<> mHat;

    // Voigt notation
    index_t VoigtMtoV [2][2];
    gsMatrix<index_t, 3, 2> VoigtVtoM;


}; // class gsFunction


} // namespace gismo
