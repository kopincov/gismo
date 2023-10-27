/** @file gsNurthbs.h

    @brief Represents a tensor-product NURBS patch

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsCore/gsRationalBasis.h>
#include <gsHSplines/gsNurthbsBasis.h>

namespace gismo
{

// forward declaration
template<short_t d, class T> class gsNurthbsBasis;

/** \brief
    A tensor product Non-Uniform Rational THB-spline function
    (NURTHBS) of parametric dimension \em d, with arbitrary target
    dimension.

    This is the geometry type associated with gsNurthbsBasis.

    \tparam d the parametric dimension of the tensor product
    \tparam T coefficient type

    \ingroup geometry
    \ingroup Nurbs
*/

template<short_t d, class T>
class gsNurthbs : public gsGeoTraits<d,T>::GeometryBase
{

public:

    typedef typename gsGeoTraits<d,T>::GeometryBase Base;
    typedef T Scalar_t;
    typedef gsNurthbsBasis<d,T>   Basis;

    /// Shared pointer for gsNurthbs
    typedef memory::shared_ptr< gsNurthbs > Ptr;

    /// Unique pointer for gsNurthbs
    typedef memory::unique_ptr< gsNurthbs > uPtr;

    /// Default empty constructor
    gsNurthbs() : Base() { }

    /// Constructor with gsNurThbsBasis and coefficients.
    gsNurthbs( const Basis & basis, gsMatrix<T> coefs ) :
    Base( basis, give(coefs) ) { }

    GISMO_CLONE_FUNCTION(gsNurthbs)

    GISMO_BASIS_ACCESSORS
};

} //namespace
