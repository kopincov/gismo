/** @file gsBernsteinBasis.h

    @brief Provides declaration of the gsTensorBernsteinBasis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsTensor/gsTensorBasis.h>
#include <gsBezier/gsBernsteinBasis.h>

namespace gismo
{


/** 
    @brief Class for a tensor product Bernstein basis

    \param T coefficient type
    \param d dimension of the parameter domain

    \ingroup basis
*/

  
template<short_t d, class T>
class gsTensorBernsteinBasis : public gsTensorBasis<d,T>  
{

public: 
    /// Base type
    typedef gsTensorBasis< d,T > Base;

    typedef gsBernsteinBasis<T> Family_t;

    /// Coordinate basis type
    typedef gsBernsteinBasis<T> Basis_t;

    /// Shared pointer for gsTensorBernsteinBasis
    typedef memory::shared_ptr< gsTensorBernsteinBasis > Ptr;

    /// Unique pointer for gsTensorBernsteinBasis
    typedef memory::unique_ptr< gsTensorBernsteinBasis > uPtr;

    /// Coefficient type
    typedef T Scalar_t;

    /// Associated geometry type
    typedef gsTensorBezier<d,T> GeometryType;

    /// Associated Boundary basis type
    typedef typename gsBernsteinTraits<d-1,T>::Basis BoundaryBasisType;

    using typename Base::iterator;
    using typename Base::const_iterator;

public:

    // Constructors forwarded from the base class
    gsTensorBernsteinBasis() : Base() { };

    gsTensorBernsteinBasis( Basis_t* x,  Basis_t*  y) : Base(x,y) { };

    gsTensorBernsteinBasis( Basis_t* x,  Basis_t* y, Basis_t* z ) : Base(x,y,z) { };

    gsTensorBernsteinBasis( std::vector<Basis_t* > const & bb ) : Base(bb) { };



public:

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "TensorBernsteinBasis of dimension " << this->dim()<< ", size "<< this->size() <<".";
        return os;
    }

    GISMO_CLONE_FUNCTION(gsTensorBernsteinBasis)
    
    GISMO_MAKE_GEOMETRY_NEW

};


} // namespace gismo
