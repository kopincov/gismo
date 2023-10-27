/** @file gsBernsteinBasis.h

    @brief Provides declaration of the gsBernsteinBasis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsBezier/gsBezier.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

// forward declaration
template <class T> class gsBezier;

/// \brief Traits for BernsteinBasis in more dimensions
template<short_t d, class T>
struct gsBernsteinTraits
{
    typedef gsTensorBernsteinBasis<d,T>  Basis;
    //typedef gsRatTensorBernstein<d,T>    RatBasis;
    typedef gsTensorBezier<d,T>          Geometry;
    //typedef gsRatTensorBezier<d,T>       RatGeometry;
};
template<class T>
struct gsBernsteinTraits<1,T>
{
    typedef gsBernsteinBasis<T>          Basis;
    //typedef gsRatBernsteinBasis<T>       RatBasis;
    typedef gsBezier<T>                  Geometry;
    //typedef gsRatBezier<T>               RatGeometry;
};
template<class T>
struct gsBernsteinTraits<0,T>
{
    typedef gsConstantBasis<T>           Basis;
    typedef gsConstantBasis<T>           RatBasis;
    typedef gsConstantFunction<T>        Geometry;
    typedef gsConstantFunction<T>        RatGeometry;
};


/** @brief
    Univariate piecewise-Bernstein basis.

    \tparam T coefficient type

    \ingroup basis
*/
  
template<class T>
class gsBernsteinBasis : public gsBasis<T>
{
public:
    typedef gsBasis<T> Base;

    typedef gsBernsteinBasis<T> Self_t;

    typedef gsBernsteinBasis<T> Family_t;

    /// Coefficient type
    typedef T Scalar_t;

    /// Associated geometry type
    typedef gsBezier<T> GeometryType;

    /// Associated Boundary basis type
    typedef gsBernsteinBasis<T> BoundaryBasisType;

    /// Dimension of the parameter domain
    static const short_t Dim = 1;

    /// Shared pointer for gsBernsteinBasis
    typedef memory::shared_ptr< gsBernsteinBasis > Ptr;

    /// Unique pointer for gsBernsteinBasis
    typedef memory::unique_ptr< gsBernsteinBasis > uPtr;

    static Ptr Make ( const gsKnotVector<T> & KV, const int & p)
    { return Ptr( new gsBernsteinBasis(KV,p) ); };

public:

    /// Default empty constructor
    gsBernsteinBasis()  : Base() { };

    /// Construct Bernstein basis along the knots KV and degree p
    gsBernsteinBasis ( const gsKnotVector<T> & KV, const int & p) :
    m_p(p), m_breaks(KV.unique(), 0)
    {
        if(p != KV.degree()  )
            std::cout << "gsBernsteinBasis Warning: Knots deg="<< KV.degree()<< " different than "<< p <<"\n";
    };
    
    gsBernsteinBasis<T>(T const& u0, T const& u1, int const& p, unsigned const & interior= 0)
    { 
        m_breaks= gsKnotVector<T>(u0,u1,interior);
        m_p = p;
    };

    gsBernsteinBasis<T>( gsBSplineBasis<T> const & bsb)
    { 
        m_breaks= gsKnotVector<T>(bsb.knots().unique(), 0);
        m_p = bsb.degree();
    };

    ~gsBernsteinBasis() { }; //destructor

public:

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    GISMO_CLONE_FUNCTION(gsBernsteinBasis)

    GISMO_MAKE_GEOMETRY_NEW

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        typename gsKnotVector<T>::const_iterator itr;
        os << "Bernstein Basis: deg=" << this->degree()
           << ", size="<< this->size() 
           << ", domain= [";
        for ( itr= m_breaks.begin(); itr != m_breaks.end(); ++itr )
            os << *itr << " ";
        os << "].\n";
        return os; };
    
    short_t domainDim() const { return Dim; }

    /// Returns the number of basis functions in the basis
    index_t size() const { return m_p * (m_breaks.size()-1) +1; }

    /// Returns the number of elements.
    size_t numElements() const { return m_breaks.size() - 1; }

    size_t numElements(boxSide const & s) const
    { return 1; }
    
    // Look at gsBasis class for a description
    gsBernsteinBasis<T>& component(short_t i) const;

    // Look at gsBasis class for a description
    void anchors_into(gsMatrix<T> & result) const;

    // Look at gsBasis class for a description
    void connectivity(const gsMatrix<T> & nodes, 
                      gsMesh<T> & mesh) const;

    // Look at gsBasis class for a description
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;

    /// Returns the indices of the basis functions that touch the domain boundary
    gsMatrix<index_t> boundary( ) const ;

    /// Returns the indices of the basis functions that touch the domain boundary
    gsMatrix<index_t> boundary(boxSide const & s, index_t offset ) const;

    /// Returns the boundary basis for side s
    GISMO_UPTR_FUNCTION_DEC(gsBernsteinBasis<T>, boundaryBasis, boxSide const &)

    // Look at gsBasis class for a description
    gsMatrix<T> support() const ;

    // Look at gsBasis class for a description
    gsMatrix<T> support(const index_t & i) const ;

    // Look at gsBasis class for a description
    // Adapted from Algorithm A2.2 from 'The NURBS BOOK' pg70.
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;


    // Look at gsBasis class for a description
    virtual void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalFunc_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, 
                               gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, 
                    gsMatrix<T>& result) const ;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    /// Evaluates the laplacians of non-zero basis functions at (the columns of) u.
    gsMatrix<T> laplacian(const gsMatrix<T> & u ) const ;
   
    /// Check the BernsteinBasis for consistency
    bool check() const
    {   // TO DO
        return true;
    };
    

//////////////////////////////////////////////////
// Additional members for univariate Bernstein basis
//////////////////////////////////////////////////

    /// Evaluates the non-zero basis functions and their
    /// first k derivatives at value u.

    /// Adapted from Algorithm A2.3 from 'The NURBS BOOK' pg72.
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n,
                                  std::vector<gsMatrix<T> >& result) const;

    // Look at gsBasis class for a description
    short_t degree(short_t i = 0) const
    { 
        GISMO_ASSERT(i==0,"Asked for degree(i) in 1D basis.");
        return m_p; 
    }

    /// Sets the degree of the basis as a reference
    inline void setDegree(short_t const & i) { m_p=i; };

    /// Returns the order of the basis
    inline unsigned order() const { GISMO_ASSERT(m_p >= -1, "order can't be negative"); return m_p+1; };

    /// Returns the starting value of the domain of the basis
    T domainStart() const { return m_breaks.first(); };

    /// Returns the ending value of the domain of the basis
    T domainEnd() const { return m_breaks.last(); };

    /// Returns the index of the first active (ie. non-zero) basis function at point u
    /// Takes into account non-clamped knots.
    inline unsigned firstActive(const T & u) const { 
        return m_p * (m_breaks.iFind(u) - m_breaks.begin());
    };

    /** \brief Number of active basis functions at an arbitrary parameter value.
     *
     * This assumes that this number doesn't change for different parameters.
     */
    inline unsigned numActive() const { return m_p + 1; }

    /// Returns the index of the first active (ie. non-zero) basis
    /// function at all columns (points) of u
    inline gsMatrix<index_t,1> * firstActive(const gsMatrix<T,1> & u) const
    {
        gsMatrix<index_t,1> * res = new gsMatrix<index_t,1>(1, u.cols());

        for( index_t i = 0; i < u.cols(); i++ )
            (*res)(0,i) = m_breaks.iFind( u(0,i) ) - m_breaks.begin();

        res->array() *= m_p;
        return res;
        // Don't forget to delete the res once you're done with it.
    }

    /// Returns the knot vector of the basis
    gsKnotVector<T> * domain() const { return const_cast<gsKnotVector<T>*>(&m_breaks); };

    /// Insert a knot
    void insertKnot(T const knot)
    { m_breaks.insert(knot); };
    
    void insertKnot_withCoefs(T const knot, gsMatrix<T> & coefs);


    /** Refine uniformly the basis by adding \a numKnots
     *  knots between every two distinct knots.
     */
    void uniformRefine(int numKnots = 1, int mul=1)
    { m_breaks.uniformRefine(numKnots); }

    void uniformRefine_withCoefs(gsMatrix<T> & coefs, int numKnots = 1, int mul=1);

    void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1, int mul=1);

    /// Apply k-refinement to the basis i times
    void uniform_k_refine(int const & i = 1) 
    { 
        m_p += i; 
        m_breaks.degreeElevate(i);
        m_breaks.uniformRefine();
    }
  
    void degreeElevate(short_t const & i = 1, short_t const dir = -1)
    { m_p+=i; }

    inline int trueSize(){ return size(); }

    // Look at gsBasis class for a description
    size_t elementIndex(const gsVector<T> & u ) const
    { return m_breaks.iFind(u(0,0)) - m_breaks.begin(); }

    // Same as gsBasis::elementIndex but argument is a value instead of a vector
    size_t elementIndex(T u ) const
    { return m_breaks.iFind(u) - m_breaks.begin(); }

// Data members
private:

    // Degree
    short_t m_p;
    // Knot vector
    gsKnotVector<T> m_breaks;


}; // class gsBernsteinBasis


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBernsteinBasis.hpp)
#endif

 
