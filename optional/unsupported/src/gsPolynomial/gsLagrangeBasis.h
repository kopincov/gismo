
#pragma once
#include <cmath>

//#include <gismo.h>
#include <gsCore/gsForwardDeclarations.h>
//#include "gsLagrangeBasis.hpp"
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

// forward declaration
template< class T>  class gsLagrangeBasis;
template <class T>  class gsLagrangePoly;

/*template<short_t d, class T>
  struct gsTraits<gsLagrangeBasis<T>, d>
  {
  typedef gsLagrangeBasis<T>   TensorBasisType;
  typedef gsLagrange<T>        TensorGeometryType;
  typedef gsLagrange<T>        TensorBoundaryType;
  };

/// Traits for LagrangeBasis in 1 dimension: specialization for d=1
template<class T>
struct gsTraits<gsLagrangeBasis<T>,1>
{
typedef gsLagrangeBasis<T>   TensorBasisType;
typedef gsLagrange<T>        TensorGeometryType;
typedef gsLagrange<T>        TensorBoundaryType;
};*/

/** @brief
    A univariate Lagrange basis.

    \tparam T coefficient type

    \ingroup basis
*/

template<class T>
class gsLagrangeBasis : public gsBasis<T>
{
public:
    typedef gsBasis<T> Base;

    /// Coefficient type
    typedef T Scalar_t;

    /// Associated geometry type
    typedef gsLagrangePoly<T> GeometryType;

    /// Associated Boundary basis type
    typedef gsLagrangeBasis<T> BoundaryBasisType;

    /// Dimension of the parameter domain
    static const short_t Dim = 1;

    /// Shared pointer for gsLagrangeBasis
    typedef memory::shared_ptr< gsLagrangeBasis > Ptr;

    /// Unique pointer for gsLagrangeBasis
    typedef memory::unique_ptr< gsLagrangeBasis > uPtr;

    static Ptr Make ( const gsVector<T> & vec, const T & start, const T & end)
    { return Ptr( new gsLagrangeBasis(vec,start,end) ); }

public:

    /// Default empty constructor
    gsLagrangeBasis()  : Base() { }

    /// Construct Lagrange basis along the parameter pars and Interval [start,end]
    gsLagrangeBasis ( const std::vector<T> & pars, const T & start, const T & end ) :
    m_p(pars.size()-1), m_breaks(pars), m_start(start), m_end(end)
    {
        check();
    }

    /// Construct Lagrange basis with equidistant breaks in the Interval [start,end]
    gsLagrangeBasis ( const T & start, const T & end, int amount_of_inner_breaks ) :
    m_p(amount_of_inner_breaks+1), m_start(start), m_end(end)
    {
        m_breaks.push_back(start);
        for(int i =0;i<amount_of_inner_breaks;i++)
        {
            m_breaks.push_back(
                start+( (i+1) * (end-start) / (amount_of_inner_breaks+1) ) );
        }
        m_breaks.push_back(end);
        check();
    }

    ~gsLagrangeBasis() { } //destructor

public:

    //////////////////////////////////////////////////
    // Virtual member functions required by the base class
    //////////////////////////////////////////////////

    /// Returns the dimension \em d of the parameter space.
    short_t dim() const { return Dim; }

    /// Returns the anchors (greville points) of the basis
    void anchors_into(gsMatrix<T>& result) const;

    /// Returns the anchor point for member \a i of the basis.
    void anchor_into(index_t i, gsMatrix<T>& result) const
    {
        result.resize(1,1);
        result(0,0)=m_breaks[i];
    }

    /// Returns the indices of active (non zero) basis functions at points (columns of) u, as a list of indices, in result
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;

    /// Returns the number of active (nonzero) basis functions at points \a u in \a result.
    void numActive_into(const gsMatrix<T> & u, gsVector<index_t>& result) const
    {
        result.resize(1,u.cols());
        for(short_t i = 0; i<u.cols();++i)
        {
            result(0,i)=m_p+1;
        }
    }

    /// Returns the indices of the basis functions that touch the domain boundary
    gsMatrix<index_t> boundary( ) const
    {
        gsMatrix<index_t> * res = new gsMatrix<index_t>(m_breaks.size(),1);
        for(short_t i = 0; i<=m_p;++i)
        {
            (*res)(i,0)=i;
        }
        return *res;
    }

    /// Returns the indices of the basis functions that touch the domain boundary
    gsMatrix<index_t> boundary(boundary::side const & s ) const { return boundary(); }

    /// Returns the boundary basis for side s.
    GISMO_UPTR_FUNCTION_DEC(gsBasis<T>, boundaryBasis, boundary::side const &)

    /// Returns a bounding box for the basis' domain
    gsMatrix<T> support() const
    {
        gsMatrix<T> * res = new gsMatrix<T>(1,2);
        *res << m_start , m_end ;
        return *res ;
    }


    /// Returns a bounding box for the basis' domain
    gsMatrix<T> support(const index_t & i) const
    {
        return support();
    }

    /// \brief Returns an interval that contains the parameter values in direction \ dir.
    ///
    /// Returns a 1x2 matrix, containing the two endpoints of the interval.
    gsMatrix<T> supportInterval(index_t dir) const;

    /// Evaluates the non-zero basis functions at value u.
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// Evaluates i-th basis functions at value u.
    void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
    {
        evalDerSingle_into(i,u,0,result);
    }

    /// Evaluates the (partial) derivatives of non-zero basis functions at (the columns of) u.
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// Evaluates the (partial)derivatives of the i-th basis function at (the columns of) u.
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        evalDerSingle_into(i,u,1,result);
    }

    /// Evaluates the (partial) derivatives of the nonzero basis functions at points \a u into \a result.
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// @brief Evaluate the (partial) derivatives of the \a i-th basis function
    /// at points \a u into \a result.
    void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        evalDerSingle_into(i,u,2,result);
    }

    /// @brief Evaluate the nonzero basis functions and their derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;

    /// @brief Evaluate the basis function \a i and its derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDersSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;

    /// @brief Evaluate the (partial) derivative(s) of order \a n the \a i-th basis function
    /// at points \a u into \a result.
    void evalDerSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;

    GISMO_CLONE_FUNCTION(gsLagrangeBasis)

    GISMO_MAKE_GEOMETRY_NEW
    
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "Lagrange Basis: deg=" << this->degree()
           << ", size="<< this->size()
           << ", parameters= [";
        for ( unsigned i = 0; i<m_breaks.size(); ++i )
            os << m_breaks[i] << " ";
        os << "], interval= [" << m_start << "," << m_end << "].\n";
        return os;
    }


    /// Prints the object as a string with extended details.
    std::string detail() const
    {
        // By default just uses print(..)
        std::ostringstream os;
        print(os);
        return os.str();
    }

    //////////////////////////////////////////////////
    // Virtual member that may be implemented or not by the derived class
    //////////////////////////////////////////////////

    /// Returns the number of basis functions in the basis
    short_t size() const { return m_p+1; }

    /// Returns the same basis for i=0, error otherwise
    gsLagrangeBasis<T>& component(index_t i) const
    {
        if ( i == 0 )
            return const_cast<gsLagrangeBasis&>(*this);
        else
            throw std::runtime_error("gsLagrangeBasis has only one component");
    }

    /// Return the gsDomain which represents the parameter domain of
    /// this basis. Currently unused.
    gsDomain<T> * domain() const
    { GISMO_NO_IMPLEMENTATION }

    /// Returns the polynomial degree.
    short_t maxDegree() const
    {
        return degree();
    }

    /// Returns the polynomial degree.
    short_t minDegree() const
    {
        return degree();
    }

    /// Returns the polynomial degree.
    short_t degree() const
    {
        return size()-1;
    }

    /// @brief Applies interpolation given the parameter values \a pts and
    /// values \a vals. May be reimplemented in derived classes.
    gsGeometry<T> * interpolateParameters(gsMatrix<T> const& pts,
                                          gsMatrix<T> const& vals) const
    { GISMO_NO_IMPLEMENTATION }

    //MEMBER FUNCTIONS (SHOULD BE PRIVATE)

    /// Stores a Matrix in result, which can be multiplied with the coefficients to
    /// get the coefficient of the Bezier representation.
    void getTransformationLagrangeBezier(gsMatrix<T> & result) const
    {
        gsMatrix<T> L_to_M;
        _getTransformationLagrangeMonomial(L_to_M);
        gsMatrix<T> M_to_B;
        _getTransformationMonomialBezier(M_to_B);
        result=M_to_B*L_to_M;
    }

    /// Check the LagrangeBasis for consistency
    bool check() const
    {
        bool consistent = true;
        // Checks if the breaks vector is strictly monotonic increasing.
        for(unsigned i=0;i<m_breaks.size()-1;i++)
            if(m_breaks[i]>=m_breaks[i+1])
            {
                GISMO_ERROR("gsLagrangeBasis Error: not strictly monoton increasing");
                consistent=false;
            }

        // Checks if the interval has zero length
        if(m_start==m_end)
        {
            GISMO_ERROR("gsLagrangeBasis Error: interval = 0");
            consistent=false;
        }
        return consistent;
    }

    /// Changes the basis so the curve is defined in the interval [0,1]
    void reparameterizeToZeroOne();

    // getters:

    /// Returns the degree of this basis.
    short_t get_m_p() const
    {
        return m_p;
    }

    /// Returns the breaks vector of this basis.
    const std::vector<T> * get_m_breaks() const
    {
        return & m_breaks;
    }

    /// Returns the start of the parameter interval.
    T get_m_start() const
    {
        return m_start;
    }

    /// Returns the end of the parameter interval.
    T get_m_end() const
    {
        return m_end;
    }

private:

    /// Stores a Matrix in result, which can be multiplied with the coefficients to
    /// get the coefficient of the Monomial representation.
    void _getTransformationLagrangeMonomial(gsMatrix<T> & result) const;

    /// Stores a Matrix in result, which can be multiplied with the coefficients to
    /// get the coefficient of the Bezier representation from Monomial coefficients.
    void _getTransformationMonomialBezier(gsMatrix<T> & result) const;

    /// takes a vector and changes it in place so it resembles the next point in the
    /// sequence. gives all the possibilities of vec.size() elements out of \a end.
    /// example with vec.size() = 3 and end = 4 (starting at [0,1,2]):
    /// [0,1,2] => [0,1,3] => [0,1,4] => [0,2,3] => [0,3,4] => [1,2,3] => [1,2,4]
    /// => [1,3,4] => [2,3,4]
    /// returns true for new element, or false if the end of sequence is reached.
    bool _nextPoint(std::vector<int> & vec, int end) const;

    /// Returns the denominator of the ith basis function.
    T _getFactor(int i) const
    {
        T prod = 1;
        for(short_t j = 0; j<=m_p;j++)
            if(i != j)
                prod*=(m_breaks[i]-m_breaks[j]);
        return prod;
    }

    // Data members
private:

    // Degree
    short_t const m_p;
    // param vector
    std::vector<T> m_breaks;
    // Interval
    T m_start;
    T m_end;




}; // class gsLagrangeBasis

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////

//include the hpp file that contains the file definitions
#include "gsLagrangeBasis.hpp"

//#ifndef GISMO_HEADERS_ONLY
//#include GISMO_HPP_HEADER(gsLagrangeBasis.hpp)
//#endif





