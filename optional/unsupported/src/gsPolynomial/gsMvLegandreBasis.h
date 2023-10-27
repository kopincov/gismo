%%%%%
#pragma once
#include <cmath>

//#include <gismo.h>
#include <gsCore/gsForwardDeclarations.h>
//#include "gsMvLegendreBasis.hpp"
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

// forward declaration
template< class T>  class gsMvLegendreBasis;
template <class T>  class gsMvLegendrePoly;

/** @brief
    A multivariate Legendre basis.

    \tparam T coefficient type

    \ingroup basis
*/
template<class T>
class gsMvLegendreBasis : public gsBasis<T>
{
public:
    typedef gsBasis<T> Base;

    /// Coefficient type
    typedef T Scalar_t;

    /// Associated geometry type
    typedef gsMvLegendrePoly<T> GeometryType;

    /// Associated Boundary basis type
    typedef gsMvLegendreBasis<T> BoundaryBasisType;

    /// Dimension of the parameter domain
    static const int Dim = 1;

    /// Shared pointer for gsMvLegendreBasis
    typedef memory::shared_ptr< gsMvLegendreBasis > Ptr;

    /// Unique pointer for gsMvLegendreBasis
    typedef memory::unique_ptr< gsMvLegendreBasis > uPtr;

    static Ptr Make ( const gsVector<T> & vec, const T & start, const T & end)  // TODO: remove?
    { return Ptr( new gsMvLegendreBasis(vec,start,end) ); }

public:
    /// Construct MvLegendre basis along the parameter pars and Interval [start,end]
    explicit gsMvLegendreBasis (const index_t deg, const T & start = 1, const T & end = 1) :
    m_ubasis(deg,start, m_end)
    { }

public:

    /// Returns the dimension \em d of the parameter space.
    int dim() const { return Dim; }

    /// Returns the indices of active (non zero) basis functions at
    /// points (columns of) u, as a list of indices, in result
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;

    /// Returns the number of active (nonzero) basis functions at points \a u in \a result.
    void numActive_into(const gsMatrix<T> & u, gsVector<index_t>& result) const
    {
        result.setConstant(1, u.cols(), m_p+1);
    }

    /// Returns a bounding box for the basis' domain
    gsMatrix<T> support() const
    {
        gsMatrix<T> res(1,2);
        res << m_start , m_end ;
        return res;
    }

    /// Returns a bounding box for the basis' domain
    gsMatrix<T> support(const index_t & i) const
    {
        return support();
    }

    /// Evaluates the non-zero basis functions at value u.
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// Evaluates the (partial) derivatives of non-zero basis functions at (the columns of) u.
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// Evaluates the (partial) derivatives of the nonzero basis
    /// functions at points \a u into \a result.
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    GISMO_CLONE_FUNCTION(gsMvLegendreBasis)

    GISMO_MAKE_GEOMETRY_NEW

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "MvLegendre Basis: deg=" << this->degree()
           << ", size="<< this->size()
           << ", interval= [" << m_start << "," << m_end << "].\n";
        return os;
    }

    /// Returns the number of basis functions in the basis
    int size() const { return m_p+1; }

    /// Return the gsDomain which represents the parameter domain of
    /// this basis. Currently unused.
    gsDomain<T> * domain() const
    { GISMO_NO_IMPLEMENTATION }

    /// Returns the polynomial degree.
    short_t degree() const
    {
        return m_p;
    }

private:

    inline T _getA(const index_t j)
    {
        return (T)(2*j-1)/j;
        return math::sqrt((T)(2*j-1)) * math::sqrt((T)(2*j+1)) / (j*math::sqrt((T)3));
    }

    inline T _getB(const index_t j)
    {
        return 0;
    }
    
    inline T _getC(const index_t j)
    {
        return (T)(j-1)/j;
        return _getA(j) / _getA(j-1);
    }

private:
    gsLegendreBasis<T> m_ubasis;
}; // class gsMvLegendreBasis

} // namespace gismo


//include the hpp file that contains the file definitions
#include "gsMvLegendreBasis.hpp"

