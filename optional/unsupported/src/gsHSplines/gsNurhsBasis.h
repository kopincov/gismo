/** @file gsNurhsBasis.h

    @brief Provides declaration of Non-uniform Rational Hierarchical B-Spline basis class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

#pragma once

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHBSpline.h>


namespace gismo
{
    /** 
     * \brief
     * A non-uniform rational hierarchical B-spline basis of parametric dimension \em d.
     *
     * See \cite Kraft1997 for the theory behind this kind of basis.
     * 
     * \tparam d the dimension of the parameter domain
     * \tparam T coefficient type
     *
     * \ingroup basis
     * \ingroup HSplines
    */ 
    
template<short_t d, class T>
class gsNurhsBasis : public gsRationalBasis< gsHBSplineBasis<d,T> >
{
public:

    /// Source basis type
    typedef gsHBSplineBasis<d,T> Src_t;


public:

    gsNurhsBasis( const Src_t & basis ) : Base(basis) { }

    gsNurhsBasis( Src_t* basis, gsMatrix<T> w ) : Base(basis, give(w)) { }

};


} // end namespace gismo


