
/** @file gsMvMonomialBasis.h

    @brief Provides dense multivariate monomial basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBasisFun.h>

#include <gsCore/gsDebug.h>

#include <gsCore/gsBoundary.h>

#include <gsCore/gsBasis.h>

#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

template<short_t d,class T>
class gsMvMonomialBasis : public gsBasis<T>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    /// Shared pointer for gsMvMonomialBasis
    typedef memory::shared_ptr< gsMvMonomialBasis > Ptr;

    /// Unique pointer for gsMvMonomialBasis
    typedef memory::unique_ptr< gsMvMonomialBasis > uPtr;

    typedef T Scalar_t;

    static const bool IsRational = false;

    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;

    /// Dimension of the parameter domain
    static const short_t Dim = d;

private:

    typedef Eigen::internal::variable_if_dynamic<unsigned,d> dimType;
    dimType m_d;
    short_t m_degree;
    std::vector<gsVector<index_t> > m_compositions;

public:
    gsMvMonomialBasis() : m_d(-1), m_degree(-1) { }

    gsMvMonomialBasis(unsigned _d, unsigned p) :  m_d(_d), m_degree(p)
    {
        getCompositions(m_compositions);
    }

    // enable_if<d"=-1>
    gsMvMonomialBasis(unsigned p) :  m_d(d), m_degree(p)
    {
        getCompositions(m_compositions);
    }

    /// Destructor
    ~gsMvMonomialBasis() { };

public:

    // Look at gsBasis class for a description
    short_t domainDim() const;

    // Look at gsBasis class for a description
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
    {
        const int sz = size();
        result.resize( sz, u.cols() );
        for ( int i = 0; i< sz; ++i )
            result.row(i).setConstant(i);// globally active
    }

    // Look at gsBasis class for a description
    gsMatrix<T> support() const;

    // Look at gsBasis class for a description
    gsMatrix<T> support(const index_t & i) const;

    // Look at gsBasis class for a description
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    // Look at gsBasis class for a description
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    // Look at gsBasis class for a description
    void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    GISMO_CLONE_FUNCTION(gsMvMonomialBasis)

    // Look at gsBasis class for a description
    memory::unique_ptr<gsGeometry<T> > makeGeometry(gsMatrix<T> coefs ) const
    { GISMO_NO_IMPLEMENTATION }

    // Look at gsBasis class for a description
    domainIter makeDomainIterator() const
    { GISMO_NO_IMPLEMENTATION }

    // Look at gsBasis class for a description
    std::ostream &print(std::ostream &os) const
    {
        os <<"Multivariate monomial basis basis of degree "
           << m_degree <<" and "<<m_d.value()<<" variables\n";
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

    /*
      Member functions that may be implemented or not in the derived class
    */

    /// The number of basis functions in this basis.
    index_t size() const;

    short_t totalDegree() const
    { return m_degree; }

    short_t degree(short_t i = 0) const
    { return m_degree; }

    int degreeOf(index_t k) const
    { return m_compositions[k].sum(); }

private:
    //create compositions for basis functions for m_degree
    void getCompositions( std::vector<gsVector<index_t> > & compos)const;

    // create compositions for basis functions for arbitrary degree
    void getCompositions( std::vector<gsVector<index_t> > & compos, index_t degree)const;

    //rt derivative
    //void rThDerivSingle(index_t r,index_t i, const gsMatrix<T> & dir, const gsMatrix<T> & u, gsMatrix<T>& result)const;

    ///find p in \a compos and return the index
    int findIndex(gsVector<index_t> const p, std::vector<gsVector<index_t> >const compos )const;

}; // class gsMvMonomialBasis


template<short_t d,class T>
short_t gsMvMonomialBasis<d,T>::domainDim() const { return m_d.value(); }

template<short_t d,class T>
index_t gsMvMonomialBasis<d,T>::size() const{
    return  m_compositions.size();
}

template<short_t d,class T>
gsMatrix<T> gsMvMonomialBasis<d,T>::support(const index_t & i) const
{return gsMatrix<T>();}

template<short_t d,class T>
gsMatrix<T> gsMvMonomialBasis<d,T>:: support() const
{return gsMatrix<T>();}

template<short_t d,class T>
void gsMvMonomialBasis<d,T>::evalSingle_into(index_t i,
                                             const gsMatrix<T> & u,
                                             gsMatrix<T>& result) const
{
    gsVector<T>  point;
    result.setZero(1,u.cols() );
    for(int j = 0; j < u.cols(); j++)
    {
        point =  u.col(j);
        T val = 1;
        for(int k = 0; k < m_compositions[i].size(); k++)
        {
            val = val * math::pow(point[k], static_cast<int>(m_compositions[i][k]));
        }
        result.at(j) = val;
    }
}

template<short_t d,class T>
void gsMvMonomialBasis<d,T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    result.setZero(m_compositions.size(), u.cols());
    gsMatrix<T> single_result;
    for(unsigned i = 0; i < m_compositions.size(); i++)
    {
        evalSingle_into(i, u, single_result);
        result.row(i) = single_result.row(0);
    }
}

template<short_t d,class T>
void gsMvMonomialBasis<d,T>::derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
GISMO_NO_IMPLEMENTATION
}

template<short_t d,class T>
void gsMvMonomialBasis<d,T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    result.setZero(m_compositions.size()*this->dim(), u.cols());
    gsMatrix<T> single_result;
    for(unsigned i = 0; i < m_compositions.size(); i++)
    {
        derivSingle_into(i, u, single_result);
        for(int j = 0; j < single_result.rows();j++)
        {
            result.row(m_d.value()*i+j) = single_result.row(j); // Map
        }
    }
}

template<short_t d,class T>
void gsMvMonomialBasis<d,T>::deriv2Single_into(index_t i, const gsMatrix<T> & u,
                                               gsMatrix<T>& result ) const
{
    result.setZero(m_d.value(), u.cols());
    gsVector<T>  point;

    for(index_t j = 0; j < u.cols(); ++j)
    {
        point =  u.col(j);
        for(index_t l = 0; l < m_compositions[i].size(); ++l)
        {
            T val = 1;
            for(index_t k = 0; k < m_compositions[i].size(); ++k)
            {
                const unsigned ex = m_compositions[i][k];
                if (ex<1 ) { val=0; break; }
                if (ex==1) { continue; }
                val = val * ex * math::pow(point[k], static_cast<int>(ex-1));
            }
            result(l,j) = val;
        }
    }
}

template<short_t d,class T>
void gsMvMonomialBasis<d,T>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    result.setZero(m_compositions.size()*((d * (d + 1)) / 2), u.cols());
    gsMatrix<T> single_result;
    for(unsigned i = 0; i < m_compositions.size(); i++)
    {
        deriv2Single_into(i, u, single_result);
        for(int j = 0; j < single_result.rows();j++){
            result.row(single_result.rows()*i+j) = single_result.row(j);
        }
    }
}

template<short_t d,class T>
void gsMvMonomialBasis<d,T>::getCompositions(std::vector<gsVector<index_t> > & compos)const{
    this->getCompositions(compos,m_degree);
}


namespace
{
bool comp_deg(const gsVector<index_t> & a, const gsVector<index_t>& b)
{ return a.sum()<b.sum(); }
}

template<short_t d,class T>
void gsMvMonomialBasis<d,T>::getCompositions(std::vector<gsVector<index_t> > & compos, index_t degree)const
{
    compos.reserve(numCompositions(degree,m_d.value()+1));
    gsVector<index_t> RR;
    firstComposition(degree,m_d.value()+1,RR);
    do
    {
        compos.push_back(RR.tail(m_d.value()));
    } while ( nextComposition(RR) );

    std::sort(compos.begin(), compos.end(), comp_deg );
}

template<short_t d,class T>
int gsMvMonomialBasis<d,T>::findIndex(gsVector<index_t> const p,
                                      std::vector<gsVector<index_t> >const compos )const
{
    //std::find
    for(index_t i = 0; i < compos.size();i++)
    {
        if(p==compos[i])return i;
    }
    return -1;
}

}; // namespace gismo


// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsMvMonomialBasis.hpp)
// #endif
