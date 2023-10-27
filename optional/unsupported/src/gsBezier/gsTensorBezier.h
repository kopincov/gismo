/** @file gsTensorBezier.h

    @brief Provides declaration of a tensor Bezier patch

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsBezier/gsTensorBernsteinBasis.h>

namespace gismo {


template <class T> class gsBernsteinBasis;
template <short_t d, class T> class gsTensorBezier;


/** \brief
    A tensor product of \em d piecewise Bezier functions, with arbitrary target dimension.

    This is the geometry type associated with gsTensorBernsteinBasis.
    
    \tparam d the parametric dimension of the tensor product
    \tparam T coefficient type

    \ingroup geometry
*/

template<short_t d, class T>
class gsTensorBezier : public gsGeoTraits<d,T>::GeometryBase
{ 

public:
    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    typedef gsTensorBernsteinBasis<d,T> Basis;

    /// Shared pointer for gsTensorBezier
    typedef memory::shared_ptr< gsTensorBezier > Ptr;

    /// Unique pointer for gsTensorBezier
    typedef memory::unique_ptr< gsTensorBezier > uPtr;

    typedef T Scalar_t;

public: 

    // TO DO : temporary, REMOVE
    gsTensorBezier( const gsTensorBasis<d,gsBernsteinBasis<T> > * basis, const gsMatrix<T> * coefs ) 
    { 
        this->m_basis =  new Basis(*basis) ;
        this->m_coefs = *coefs ;    
    }


    gsTensorBezier( const Basis * basis, const gsMatrix<T> * coefs )
    : Base (basis,coefs) { }

    gsTensorBezier( const Basis & basis, const gsMatrix<T> & coefs )
    : Base (basis,coefs) { }

    GISMO_CLONE_FUNCTION(gsTensorBezier)

    GISMO_BASIS_ACCESSORS
    
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const { os<<"gsTensorBezier\n"; return os; };

};

} // namespace gismo
