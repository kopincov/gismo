
#pragma once

#include <iostream>
#include <gsCore/gsBasis.h>


namespace gismo
{

/** \brief
    An univariate monomial polynomial basis.
    If the degree is \a p the basis is given by:
    \[ 1, x, x^2, ..., x^p \]
    The basis functions are numbered, starting from zero, as stated above.

    \tparam T coefficient type

    \ingroup basis
*/

template<class T>
class gsMonomialBasis : public gsBasis<T>
{
public:
    /// Shared pointer for gsMonomialBasis
    typedef memory::shared_ptr< gsMonomialBasis > Ptr;

    /// Unique pointer for gsMonomialBasis
    typedef memory::unique_ptr< gsMonomialBasis > uPtr;

    /// Dimension of the parameter domain
    static const int Dim = 1;

    /// Default constructor, which sets the degree to 0.
    gsMonomialBasis()  : m_p(0)
    { }

    /// \brief Constructs a monomial basis
    /// \param  degree degree of the monomial basis, which has to be >=0
    explicit gsMonomialBasis(int degree);
    
    // Pure virtual member functions required by the base class
    /// Returns the dimension of the parameter domain
    virtual short_t domainDim() const;

    GISMO_CLONE_FUNCTION(gsMonomialBasis)

    virtual memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const;
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const;
    /// Returns the numer of basis functions.
    virtual index_t size() const;


    // Virtual member functions of the base class which need to be redefined

    /// Returns the degree of the basis
    virtual short_t degree() const { return m_p; }

    virtual short_t setDegree() const { return m_p; } // TODO: really?

    /// Returns the degree of the basis (non-virtual)
    inline short_t deg() const { return m_p; }

    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
    {
        result.resize(m_p+1,u.cols());
        result.colwise() = 
            gsVector<index_t>::LinSpaced(m_p+1, 0, m_p); //.template cast<unsigned>();
    }
    
    /// \brief Evaluates the basis functions at values \a u.
    /// \param      u evaluation values as <em>1 x m</em> matrix.
    /// \param[out] result is a <em>size x m</em> matrix, where each row corresponds to a basis function
    /// and each column corresponds to an evaluation value.
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// Evaluates the \a i-th basis function at values \a u.
    virtual void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// \brief Evaluates the derivative of the basis functions at values \a u.
    /// \param      u evaluation values as <em>1 x m</em> matrix.
    /// \param[out] result is a <em>size x m</em> matrix, where each row corresponds to a basis function
    /// and each column corresponds to an evaluation value.
    virtual void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// Evaluates the derivatives of the \a i-th basis function at values \a u.
    virtual void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// \brief Evaluates the second derivative of the basis functions at values \a u.
    /// \param u    evaluation values as <em>1 x m</em> matrix.
    /// \param[out] result is a <em>size x m</em> matrix, where each row corresponds to a basis function
    /// and each column corresponds to an evaluation value.
    virtual void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// Evaluates the second derivatives of the \a i-th basis function at values \a u.
    virtual void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// \brief Evaluates all derivatives up to order \a n of the basis functions at values \a u.
    /// \param u    evaluation values as <em>1 x m</em> matrix.
    /// \param n    derivative order
    /// \param[out] result is a <em>size*(n+1) x m</em> matrix, where each column corresponds to an evaluation value
    /// and the first \a size rows include the function values, the second \a size rows include the first derivatives and so on up to order \a n.
    /// Each of the \a size rows corresponds to a basis function.
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;
    /// \brief Evaluates all derivatives up to order \a n of the i-th basis function at values \a u.
    /// \param u    evaluation values as <em>1 x m</em> matrix.
    /// \param n    derivative order
    /// \param[out] result is a <em>(n+1) x m</em> matrix, where each column corresponds to an evaluation value
    /// and each row correspond to a derivative starting from 0 up to order \a n.
    virtual void evalAllDersSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;

    /// \brief Evaluates the derivative of order \a n of the i-th basis function at values \a u.
    /// \param u    evaluation values as <em>1 x m</em> matrix.
    /// \param n    derivative order
    /// \param[out] result is a <em>1 x m</em> matrix, where each column corresponds to an evaluation value.
    virtual void evalDerSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;

    /// \brief Returns (a bounding box for) the domain of the whole basis.
    /// Returns a <em>1 x 2</em> matrix, containing the two diagonally extreme
    /// corners of a hypercube.
    virtual gsMatrix<T> support() const
    {
        gsMatrix<real_t> support(1,2);
        support<<0, 1;
        return support;
    }

    /// \brief Returns (a bounding box for) the support of the i-th basis function.
    /// Returns a <em>1 x 2</em> matrix, containing the two diagonally extreme
    /// corners of a hypercube.
    virtual gsMatrix<T> support(const index_t & ) const
    {
        gsMatrix<real_t> support(1,2);
        support<<0, 1;
        return support;
    }

    short_t degree(short_t i) const { return m_p; }

    /** \brief Evaluate the function described by \a coefs at points \a u,
     * i.e., evaluates a linear combination of coefs x BasisFunctions, into \a result.
     *
     * This function overrides the default implementation given in the base class.
     * The so called Horner scheme is used for the evaluation.
     *
     * \param u     evaluation values as <em>1 x m</em> matrix
     * \param coefs <em>size x n</em> coefficient matrix describing the geometry in this basis (\em n is the dimension of the coefficients)
     * \param[out] result  a matrix of size <em>n x m</em> with one function value as a column vector
     *              per evaluation point
     */
    virtual void evalFunc_into(const gsMatrix<T> & u,
                               const gsMatrix<T> & coefs,
                               gsMatrix<T>& result) const;



// Data members
private:
    /// Degree
    short_t m_p;

};  // class gsMonomialBasis

}   // namespce gismo


#include <gsPolynomial/gsMonomialBasis.hpp>
