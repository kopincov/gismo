
#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsBoundary.h>

#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

  /** 
      @brief Class representing a Box Spline Basis, derived from gsBasis

      Template parameter
      \param T is the coefficient type

      \param d is the dimension of the domain (d=2 for surfaces, d=3 for solids)
  */
  
template<short_t d, class T>
class gsBoxSplineBasis : public gsBasis<T> // derived from gsBasis
{
public:

  typedef gsBasis<T> Base;

  /// Coefficient type
  typedef T Scalar_t;

  /// Associated geometry type
  //typedef gsGenericGeometry< gsBoxSplineBasis > GeometryType;

  /// Associated Boundary basis type
  typedef gsBoxSplineBasis<d-1,T> BoundaryBasisType;

  /// Dimension of the parameter domain
  static const int Dim = d; // d = 2 for surfaces, d = 3 for solids

  /// Shared pointer for gsBoxSplineBasis
  typedef memory::shared_ptr< gsBoxSplineBasis > Ptr;

  /// Unique pointer for gsBoxSplineBasis
  typedef memory::unique_ptr< gsBoxSplineBasis > uPtr;

public:

    /// Default empty constructor
    gsBoxSplineBasis() : Base() { };
    
    /**
     * \brief Constructs a box spline with respect to a direction matrix of ALL directions
     *
     * \param X a matrix of all directions
     *
      */
    gsBoxSplineBasis(const gsMatrix<T,d> & X)
        {
            assert( X.rows() == d ) ;

            m_X = getX(X);
            m_n = getn(X);

            check();

        }

    /**
     * \brief Constructs a box spline with respect to a direction matrix of DISTINCT directions
     *        and a multiplicity vector
     *
     * \param X a matrix of distinct directions
     * \param n a vector of multiplicities of the distinct directions
     *
      */
    gsBoxSplineBasis(const gsMatrix<T,d>& X, const gsVector<int>& n) : m_X(X), m_n(n)
        { 
            assert ( m_X.rows() == d ) ;
            assert ( m_X.cols() == m_n.rows() );

            check();
        }

    /**
     * \brief Constructs a box spline with respect to a direction matrix with main directions
     *       (1,0), (0,1), (1,1), (-1,1); ONLY FOR d = 2!
     *
     * \param i multiplicity of direction (1,0)
     * \param j multiplicity of direction (0,1)
     * \param k multiplicity of direction (1,1) (default value = 0)
     * \param l multiplicity of direction (-1,1) (default value = 0)
     *
      */
    gsBoxSplineBasis(const int i, const int j, const int k = 0, const int l = 0)
        { 
        assert ( d == 2 );

        int ncol = static_cast<int>(i!=0) + static_cast<int>(j!=0) + static_cast<int>(k!=0) + static_cast<int>(l!=0);
        m_n.resize(ncol);        
        //m_n(0) = i;
        //m_n(1) = j;

        m_X.resize(d,ncol);
        m_X.setZero();

        //for (unsigned p = 0; p < d; p++) m_X(p,p) = 1;

        int c = 0;
        if (i > 0) {
            m_n(c) = i;
            m_X(0,c++) = 1;
        }

        if (j > 0) {
            m_n(c) = j;
            m_X(1,c++) = 1;
        }

        if (k > 0) {
            m_n(c) = k;
            m_X.col(c++).setOnes();
        }

        if (l > 0) {
            m_n(c) = l;
            m_X.col(c) << -1,1;
        }

        /*

        m_fullX.resize(d, m_n.sum() );
        m_fullX.setOnes();
        for ( int s=0; s < i; s++)   m_fullX(1,s)=0;

        //m_fullX.leftColumns(i) = m_X.col(0).replicate(i); // Eigen block operations
        for ( int s=0; s < j; s++)   m_fullX(0,i+s)=0;
        //m_fullX.middleColumns(i) = m_X.col(0).replicate(i);
        for ( int s=0; s < l; s++)   m_fullX(0,i+j+k+s)=-1;
        */

	    check();

	}

    // TO DO: add domain info (n,m,l)

    GISMO_CLONE_FUNCTION(gsBoxSplineBasis)

public:

  void swap(gsBoxSplineBasis& other)
    {
      m_X.swap   (other.m_X   );
    }

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

  short_t domainDim() const { return Dim; }

  /// Returns the number of basis functions in the basis
  index_t size() const { return 1; }

  /// Returns the bounding box of the support of the default basis function
  gsMatrix<T> support() const ;

  /// Returns the bounding box of the support of the i-th basis function
  gsMatrix<T> support(const index_t & i) const ;

  /// Returns the boundary basis for side s
  //GISMO_UPTR_FUNCTION_DEC(gsBoxSplineBasis, boundaryBasis, boundary::side const &)

  /// Returns a bounding box for the basis' domain
  gsMatrix<T> * parameterRange() const ;

  /// Evaluates the non-zero basis functions at value u.
  virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

  /// Evaluates i-th basis functions at value u.  
  virtual void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

  /// Evaluates i-th basis functions at value u.
  virtual gsMatrix<T> evalresult(const gsMatrix<T> & u) const;

  /// Evaluates derivative of basis functions at value u.
  virtual gsMatrix<T> derivativeresult(const gsMatrix<T> & u, const gsVector<int> & r) const;

  /// Evaluates the (partial) derivatives of non-zero basis functions at (the columns of) u.
  void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

  /// Evaluates the (partial)derivatives of the i-th basis function at (the columns of) u.
  void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

  /// Evaluates the (partial) derivatives of a B-Spline given by coefs at (the columns of) u.
  void deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

  /// Evaluates the (partial) second derivatives of a B-Spline given by coefs at (the columns of) u.
  gsMatrix<T> deriv2(const gsMatrix<T> & u ) const ;

  // TODO: replace by GISMO_MAKE_GEOMETRY_NEW when GeometryType exists
  memory::unique_ptr<gsGeometry<T> > makeGeometry(gsMatrix<T> coefs ) const;

  /// Returns the matrix of DISTINCT directions out of full matrix (1st constructor)
  /// \param X direction matrix of all directions
  gsMatrix<T,d> getX(const gsMatrix<T,d> & X) const ;

  /// Returns the vector of multiplicities out of full matrix (1st constructor)
  /// \param X direction matrix of all directions
  gsVector<int> getn(const gsMatrix<T,d> & X) const ;

  /// Returns the matrix of ALL directions
  /// \param X direction matrix of distinct directions
  /// \param n vector of multiplicities
  gsMatrix<T,d> getFullX(const gsMatrix<T,d> & X, const gsVector<int> & n) const ;

  /// Check the BoxSplineBasis for consistency
  bool check() const
    {
        if ( m_X.colPivHouseholderQr().rank() < (index_t)(d) )
	    {
        gsWarn<< "gsBoxSplineBasis Error: Matrix doesn't have full rank:\n"<< m_X;
		return false;
	    }
	return true;
    }

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const;

//////////////////////////////////////////////////
// Additional members for univariate Spline basis
//////////////////////////////////////////////////


  /// Returns the degree of the basis
  //short_t degree() const { return m_p; };
  //using Base::degree;

  /// Sets the degree of the basis as a reference
  //inline void set_degree(int const & i) { m_knots.set_degree(i); m_p=i; };

  /// Returns the order of the basis
  //inline unsigned order() const { return m_p+1; };

  /// Returns the starting value of the domain of the basis
  //T domainStart() const { return m_knots[m_p]; };

  /// Returns the ending value of the domain of the basis
  //T domainEnd() const { return m_knots[this->size()]; };

  /// Returns the index of the first active (ie. non-zero) basis function at point u
  /// Takes into account non-clamped knots.
  //inline unsigned firstActive(T u) const { 
  //   return m_knots.findspan(u)-m_p; 
  //};

  /** \brief Number of active basis functions at an arbitrary parameter value.
   *
   * This assumes that this number doesn't change for different parameters.
   */
  inline unsigned numActive() const { return 1; }

  /// Returns the index of the first active (ie. non-zero) basis
  /// function at all columns (points) of u
  //inline gsMatrix<unsigned,1> * firstActive(const gsMatrix<T,1> & u) const { 
  //  gsMatrix<unsigned,1> * a = m_knots.findspan(u) ;
  //  *a -= gsMatrix<unsigned,1>::Constant(1,u.cols(),m_p )  ; 
  //  return a;
  //};

  /// Returns the knot vector of the basis
  //gsDomain * domain() const { return NULL};

  /** Refine uniformly the basis by adding \a numKnots
   *  knots between every two distinct knots.
   */
  void uniformRefine(int numKnots = 1, int mul = 1) { }

  /// Refine the basis uniformly and perform knot refinement for the given coefficient vector
  void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul = 1);

  /// Refine the basis uniformly and produce a sparse matrix which maps coarse coefficient vectors to refined ones
  void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1, int mul = 1);

// Internal computations
private:

  /*
   *static T _boxSplineRecBruteForce( const gsMatrix<T> & X, const gsMatrix<T> & u );
*/

  /**
   * @brief _boxSplineRec recursion function in the evaluation of a box spline
   * @param n vector of multiplicities updated in the recursion step
   * @param m current position in the recursion tree (for computing the delayed translation)
   * @param N vector of normals to the hyperplanes
   * @param Y matrix for computing the least squares approximation
   * @param t vector of coefficients from the least squares approximation
   * @param J row vector of ones for matrix computations
   * @param u vector of points where the spline is to be evaluated
   * @return returns the vector of values of a box spline at points u
   */
  gsMatrix<T> _boxSplineRec( const gsVector<int> & n,
                             const gsVector<int> & m,
                             const gsMatrix<T,d> & N,
                             const gsMatrix<T> & Y,
                             const gsMatrix<T> & t,
                             const gsMatrix<T,1> & J,
                             const gsMatrix<T> & u ) const;

  /**
   * @brief _computeBoxNormals computes the normals to all hyperplanes in m_X (by Kobbelt, 1996)
   *        working for d = 2 and d = 3
   * @param r number of rows to be selected before the base case is reached
   * @param k next direction in m_X to be considered for selection
   * @param M bitvector indicating the selected direction vectors in m_X
   * @return returns the matrix of normal vectors to the hyperplanes
   */
  gsMatrix<T,d> _computeBoxNormals(const int & r,
                                   const int & k,
                                   const gsVector<T> & M) const;

  //void expandDirection(gsMatrix<T,d> & XX);
    
    
// Data members
private:

  // Direction matrix with NO multiple directions
  gsMatrix<T,d> m_X ;

  // Vector of multiplicities of the directions with respect to m_X
  gsVector<int> m_n ;
  
  // Zonotope given as scalars corresponding to (columns of) the direction matrix?
  gsVector<int> m_zonotope ;

  // Domain given by knots, to do: add diagonals?
  gsKnotVector<T> m_knots[d];

  // Bezier matrix on reference triangle (numActivexmonomials)

}; // class gsBoxSplineBasis


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBoxSplineBasis.hpp)
#endif

