// defined over a Lagrange basis

#pragma once

#include <ostream>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>
#include "gsLagrangeBasis.h"


namespace gismo
{

/** \brief
    The geometry class of a Lagrange Polyomial curve

    This is the geometry type associated with gsLagrangeBasis.

    \tparam T coefficient type

    \ingroup geometry
*/

template<class T>
class gsLagrangePoly : public gsGeoTraits<1,T>::GeometryBase
{

public:
    typedef gsLagrangeBasis<T> Basis;
    typedef typename gsGeoTraits<1,T>::GeometryBase Base;

    /// Shared pointer for gsLagrangePoly
    typedef memory::shared_ptr< gsLagrangePoly > Ptr;

    /// Unique pointer for gsLagrangePoly
    typedef memory::unique_ptr< gsLagrangePoly > uPtr;

    /// Default empty constructor
    gsLagrangePoly() : Base() { }

    /// Construct B-Spline by basis and coefficient matrix
    gsLagrangePoly( const Basis & basis, gsMatrix<T> coefs ) :
    Base( basis, give(coefs)) { }

    /// Construct Lagrange curve by degree, coefficient matrix and domain
    gsLagrangePoly( const unsigned & p, const gsMatrix<T> & coefs, const T & u0=0, const T & u1= 1) : Base()
    {
        gsLagrangeBasis<T> * lagrange_basis = new Basis(u0,u1,p-1);
        Base::m_coefs=coefs;
        Base::m_basis=lagrange_basis;
        GISMO_ASSERT( lagrange_basis->size() == this->m_coefs.rows(),
                      "The coefficient matrix of the geometry (rows="<<this->m_coefs.rows()<<") does not match the number of basis functions in its basis("<< lagrange_basis->size() <<").");
    }

    /// Construct a Lagrange curve that resembles the \a part -th part of a bezier sequence
    /// by interpolating at degree+1 values
    gsLagrangePoly( const gsBezier<T> & bezier_curve, int part) : Base()
    {
        gsBernsteinBasis<T> basis = bezier_curve.basis();
        gsKnotVector<T> * bezier_knots = basis.domain();
        T start = (*bezier_knots)[part];
        T end = (*bezier_knots)[part+1];
        gsLagrangeBasis<T> * lagrange_basis = new gsLagrangeBasis<T>(0,1, basis.degree()-1);
        gsMatrix<T> lagrange_coeffs;
        const std::vector<T> * lagrange_breaks = lagrange_basis->get_m_breaks();
        gsMatrix<T> u(1,lagrange_breaks->size());
        for(unsigned i = 0;i<lagrange_breaks->size();i++)
        {
            u(0,i)=start+(end-start)*lagrange_breaks->at(i);
        }
        bezier_curve.eval_into(u,lagrange_coeffs);
        Base::m_coefs=lagrange_coeffs.transpose();
        Base::m_basis=lagrange_basis;
        GISMO_ASSERT( lagrange_basis->size() == this->m_coefs.rows(),
                      "The coefficient matrix of the geometry (rows="<<this->m_coefs.rows()<<") does not match the number of basis functions in its basis("<< lagrange_basis->size() <<").");
    }

    
public:

    GISMO_BASIS_ACCESSORS
        
    GISMO_CLONE_FUNCTION(gsLagrangePoly)

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "Lagrange curve of degree "<<
            this->basis().degree()<< ", over [" << domainStart()
         << " " << domainEnd() << "] with breaks ";
        for ( typename std::vector<T>::const_iterator itr=
              this->basis().get_m_breaks()->begin(); itr != this->basis().get_m_breaks()->end(); ++itr )
            os << *itr << ", ";
        return os;
    }


//////////////////////////////////////////////////
// Additional members for univariate B-Splines
//////////////////////////////////////////////////

    /// Returns the starting value of the domain of the basis
    T domainStart() const { return this->basis().get_m_start(); }

    /// Returns the end value of the domain of the basis
    T domainEnd() const { return this->basis().get_m_end(); }

    /// gives back a BezierCurve, that resembles this Lagrange Curve in the parametric
    /// domain [0,1]. If the parametric domain of the Lagrange Curve is different from
    /// [0,1], the method reparameterizeToZeroOne() can be called to change this.
    gsBezier<T> * transformToBezier()
    {
        gsMatrix<T> transformMat;
        this->basis().getTransformationLagrangeBezier(transformMat);
        gsMatrix<T> newCoefs = transformMat*this->coefs();
        unsigned deg = this->basis().degree();
        gsBezier<T> * bezCurve = new gsBezier<T>(deg,newCoefs);
        return bezCurve;
    }

    /// reparameterize this curve to [0,1]
    void reparameterizeToZeroOne()
    {
        this->basis().reparameterizeToZeroOne();
    }

// Data members
private:

}; // class gsLagrangePoly



}; // namespace gismo
