/** @file gsDeCasteljau.hpp

    @brief Provides definitions of the gsBernsteinBasis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsNurbs/gsKnotVector.h>


namespace gismo {

/// Executes DeCausteljau's algorithm for Bezier evaluation on points
/// u, breaks, degree deg, coefficient matrix coefs
template<class T>
void gsDeCasteljauEval( 
    const gsMatrix<T> &u,
    const gsKnotVector<T> & breaks,
    int deg,
    const gsMatrix<T> & coefs,
    gsMatrix<T>& result
    )
{
    assert( u.rows() == 1 ) ;
    result.resize( coefs.cols(), u.cols() ) ;
    gsMatrix<T> points(deg+1, coefs.cols() );

    for ( index_t j=0; j< u.cols(); j++ ) // for all points (entries of u)
	{
	    // Get to the correct curve piece
        const int ind = breaks.iFind( u(0,j) ) - breaks.begin();

	    // Transform interval
	    const T u0 = ( u(0,j) - breaks[ind])/(breaks[ind+1]-breaks[ind]) ;
	    const T u1 = T(1.0) - u0 ;
	    // Initialize points with the coefficients
	    points.noalias() = coefs.middleRows( deg*ind, deg+1 );

	    // Compute triangular scheme
	    for ( int k=deg; k!= 0; --k ) 
		for ( int i=0; i!=k; i++ ) 
		    points.row(i) = u1 * points.row(i) + u0 * points.row(i+1);    

	    // Same loop as above for triangular scheme
	    // for ( int k=1; k<= deg; k++ ) 
	    // 	for ( int i=0; i<=deg-k; i++ ) 
	    // 	    points.row(i) = u1 * points.row(i) + u0 * points.row(i+1);

	    // Save result
	    result.col(j).noalias() = points.row(0);
	}
}

/// Executes Horner's scheme for Bezier evaluation on points
/// u, breaks, degree deg, coefficient matrix coefs
/// TO DO :compare with gsDeCasteljauEval
template<class T>
void gsHornerBezier( 
    const gsMatrix<T> &u,
    const gsKnotVector<T> & breaks,
    int deg,
    const gsMatrix<T> & coefs,
    gsMatrix<T>& result
    )
{
    assert( u.rows() == 1 ) ;
    result.resize( coefs.cols(), u.cols() ) ;

    gsMatrix<T> tmp(1,coefs.cols() );
    T utmp;
    int binc; // binomial coefficient

    for ( index_t j=0; j< u.cols(); j++ ) // for all points (entries of u)
	{
        const int ind = breaks.iFind( u(0,j) ) - breaks.begin();
	    const T u0 = ( u(0,j) - breaks[ind])/(breaks[ind+1]-breaks[ind]) ;
	    const T u1 = T(1.0) - u0 ;

	    tmp.noalias() = u1 * coefs.row(deg*ind) ;
            utmp = 1.0;
            binc = 1;
	    for ( int k=1; k!=deg; ++k ) 
            {
                utmp *= u0;
                binc *= (deg-k+1)/k;
                tmp = (tmp + utmp * binc * coefs.row( deg*ind+k ) ) * u1;           
            }
            result.col(j).noalias() = tmp + utmp* u0 * coefs.row( deg*(ind+1) );
	}
}

/// Executes DeCausteljau's algorithm to subdivide a bezier curve at value val
template<class T, class Mat> 
void gsDeCasteljauSubdivide( 
    gsKnotVector<T> & breaks,
    const int deg,
    Mat & coefs,
    const T val,
    const bool update_breaks = true
    )
{
    Mat points;
    const int np= coefs.rows()-1;

    // Get to the correct curve piece
    const int ind = breaks.iFind( val ) - breaks.begin();

    // Initialize points with the coefficients
    // to do: remove points and compute inplace on coefs
    points = coefs.middleRows( deg*ind, deg+1 );

    // Make room for new coefficients
    coefs.conservativeResize( coefs.rows() + deg, coefs.cols() );

    // shift control points that are not affected
    for( index_t i = np; i>= deg*(ind+1); --i )
        coefs.row(i+deg) = coefs.row(i);
    
    coefs.middleRows( deg*(ind+1), deg+1 ) = coefs.middleRows( deg*(ind), deg+1 );

    // Transform interval
    const T u0 = ( val - breaks[ind])/(breaks[ind+1]-breaks[ind]) ;
    const T u1 = T(1.0) - u0 ;
    
    // Compute triangular scheme
    int s = deg*ind;
    for ( int k=deg; k!= 0; --k ) 
    {
        coefs.row(s++) = points.row(0);// Collect left 
        for ( int i=0; i!=k; i++ ) 
            points.row(i) = u1 * points.row(i) + u0 * points.row(i+1);    
    }
 
    coefs.middleRows( deg*(ind+1), deg+1 ) = points;

    // Update breaks
    if ( update_breaks )
        breaks.insert(val);    
}

}; //end namespace gismo
