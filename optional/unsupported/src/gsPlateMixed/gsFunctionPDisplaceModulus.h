/** @file gsFunctionPDisplaceModulus.h

    @brief Function to compute the physical displacement

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
class gsFunctionPDisplaceModulus : public gsFunction<T>
{
public:
    typedef gsFunction<T> Base;

    /// Shared pointer for gsFunctionPDisplaceModulus
    typedef memory::shared_ptr< gsFunctionPDisplaceModulus > Ptr;

    /// Unique pointer for gsFunctionPDisplaceModulus
    typedef memory::unique_ptr< gsFunctionPDisplaceModulus > uPtr;


    using Base::support;
    using Base::domainDim;
    using Base::targetDim;
public:
    gsFunctionPDisplaceModulus(const gsGeometry<T>& geometry,const gsGeometry<T>& u_co):
        m_geometry(geometry), m_u_co(u_co)
    {}

    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsFunctionPDisplaceModulus, clone)

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
        unsigned evFlags = NEED_JACOBIAN | NEED_MEASURE | NEED_NORMAL;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_geometry));
        geoEval->evaluateAt(u);
        m_u_co.eval_into(u, result);

        for(index_t N=0; N!=u.cols(); ++N)
        {

            geoEval->normal(N,normal);
            normal.normalize();
            F.leftCols(2) = geoEval->jacobian(N);
            F.col(2)      = normal;
            FInvT = F.inverse().transpose();

            // compute physical displacement from covariant components of displacement
            // check for singular point of geometry
            if(math::abs(geoEval->jacDet(N)) > 1e-14)
                result.col(N) = FInvT * result.col(N);
            else
                gsInfo<<"Geometry is evaluated at singular point! \n";

            result(0,N) = math::abs(result(0,N));
            result(1,N) = math::abs(result(1,N));
            result(2,N) = math::abs(result(2,N));

        }

        //result.row(2) *= 1e5;
    }

    virtual short_t domainDim () const
    {
        return 2;
    }

private:

    const gsGeometry<T>& m_geometry;
    const gsGeometry<T>& m_u_co;

    mutable gsVector<T> normal; // Normal to the shell centerline
    mutable gsMatrix<T,3,3> F; // deformation gradient (a_1, a_2, a_3)
    mutable gsMatrix<T,3,3> FInvT; // (a^1, a^2, a^3)


}; // class gsFunction


} // namespace gismo
