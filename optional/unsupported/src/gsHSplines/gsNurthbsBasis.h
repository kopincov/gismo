/** @file gsNurthbsBasis.h

    @brief Provides declaration of TensorNurbsBasis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

#pragma once

#include <gsCore/gsRationalBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsHSplines/gsNurthbs.h>

namespace gismo
{

// forward declaration
template<short_t d, class T> class gsNurthbs;

/** \brief
    A tensor product Non-Uniform Rational THB-spline (NURTHBS) basis.

    This is the rational version of gsTHBSplineBasis.\n
    See gsNurthbs for the associated gsGeometry.

    \tparam d dimension of the parameter domain
    \tparam T coefficient type

    \ingroup basis
    \ingroup Nurbs
*/
template<short_t d, class T>
class gsNurthbsBasis :
        public gsRationalBasis< gsTHBSplineBasis<d,T> >
{
    typedef gsRationalBasis<gsTHBSplineBasis<d,T> > Basis;

    /// Shared pointer for gsNurthbsBasis
    typedef memory::shared_ptr< gsNurthbsBasis > Ptr;

    /// Unique pointer for gsNurthbsBasis
    typedef memory::unique_ptr< gsNurthbsBasis > uPtr;

public:

    GISMO_CLONE_FUNCTION(gsNurthbsBasis)

    /// Construct NURTHBS basis by a THBSplineBasis plus weights
    gsNurthbsBasis( gsTHBSplineBasis<d,T> *bs, const gsMatrix<T> & w) :
    Basis( bs, w )  { }

    ~gsNurthbsBasis() { }  //destructor

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "Rational THB-spline basis. Underlying basis:" << std::endl;
        os << *m_src << std::endl << std::endl;
        return os;
    }

    /// @brief Create a gsGeometry of proper type for this basis with the
    /// given coefficient matrix.
    virtual memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const
    {
        return memory::unique_ptr<gsGeometry<T> >(new gsNurthbs<d,T>( *this , give(coefs) ));
    }

    GISMO_BASIS_ACCESSORS

    using gsRationalBasis< gsTHBSplineBasis<d,T> >::m_src;
    using gsRationalBasis< gsTHBSplineBasis<d,T> >::m_weights;
};




} //namespace
