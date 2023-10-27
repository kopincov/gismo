#pragma once

#include "gsMonomialBasis.h"
#include <math.h>

namespace gismo
{

template <class T>
gsMonomialBasis<T> :: gsMonomialBasis(int degree)
{
    if(degree>=0)
        m_p=degree;
    else
        std::cout<<"Error: The degree has to be >=0"<<std::endl;
}

template <class T>
short_t gsMonomialBasis<T> :: domainDim() const
{
    return 1;   // Since we consider univariate basis
}

template <class T>
std::ostream & gsMonomialBasis<T> :: print(std::ostream &os) const
{
    os<<"gsMonomialBasis"<<std::endl;
    os<<"degree: "<<m_p<<std::endl;
    os<<"size: "<<size()<<std::endl;
    return os;
}

template <class T>
index_t gsMonomialBasis<T> :: size() const
{
    return m_p+1;
}

template <class T>
memory::unique_ptr<gsGeometry<T> > gsMonomialBasis<T> :: makeGeometry(gsMatrix<T> ) const
{
    GISMO_NO_IMPLEMENTATION;
}

template <class T>
void gsMonomialBasis<T> :: eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    result.resize(size(), u.cols());

    // 0th basis function is constant 1
    result.row(0).setOnes();

    for(int p=1; p<=m_p; p++)  // For all monomials up to m_p
        for(int eval_point=0; eval_point<u.cols(); eval_point++)    // For all values (columns of u)
            result(p,eval_point)=result(p-1,eval_point)*u(0,eval_point);

}

template <class T>
void gsMonomialBasis<T> :: evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    GISMO_ENSURE(static_cast<int>(i)<=m_p,"Error: There is no i-th basis function in the basis."<<std::endl<<
                 "Number of basis functions in the basis is "<<size()<<".");
        
    result.resize(1,u.cols());
    
    for(int eval_point=0; eval_point<u.cols(); eval_point++)    // For all values (columns of u)
        result(0,eval_point) = math::pow(u(0,eval_point),static_cast<int>(i));
}

template <class T>
void gsMonomialBasis<T> :: deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    result.resize(size(), u.cols());

    // Derivative 0th basis function is constant 0
    result.row(0).setZero();
    // Derivative 1st basis function is constant 1
    if(m_p>=1)
        result.row(1).setOnes();

    for(int p=2; p<=m_p; p++)  // For all monomials up to m_p
        for(int eval_point=0; eval_point<u.cols(); eval_point++)    // For all values (columns of u)
            result(p,eval_point)=(real_t)p/(p-1)*result(p-1,eval_point)*u(0,eval_point);
}

template <class T>
void gsMonomialBasis<T> :: derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    GISMO_ENSURE(static_cast<int>(i)<=m_p,"Error: There is no i-th basis function in the basis."<<std::endl<<
                 "Number of basis functions in the basis is "<<size()<<".");

    result.resize(1,u.cols());

    if(i==0)
        result.setZero(1,u.cols());
    else
        for(int eval_point=0; eval_point<u.cols(); eval_point++)    // For all values (columns of u)
            result(0,eval_point)=i* math::pow(u(0,eval_point),static_cast<int>(i)-1);

}

template <class T>
void gsMonomialBasis<T> :: deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    result.resize(size(), u.cols());

    // 2nd derivative of 0th and 1st basis function is constant 0
    result.row(0).setZero();
    if(m_p>=1)
        result.row(1).setZero();
    // 2nd derivative of 2nd basis function is constant 2
    if(m_p>=2)
        result.row(2).setConstant(2);

    for(int p=3; p<=m_p; p++)  // For all monomials up to m_p
        for(int eval_point=0; eval_point<u.cols(); eval_point++)    // For all values (columns of u)
            result(p,eval_point)=(real_t)p/(p-2)*result(p-1,eval_point)*u(0,eval_point);
}

template <class T>
void gsMonomialBasis<T> :: deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    GISMO_ENSURE(static_cast<int>(i)<=m_p,"Error: There is no i-th basis function in the basis."<<std::endl<<
                 "Number of basis functions in the basis is "<<size()<<".");

    result.resize(1,u.cols());

    if(i==0 || i==1)
        result.setZero(1,u.cols());
    else
        for(int eval_point=0; eval_point<u.cols(); eval_point++)
            result(0,eval_point)=i*(i-1) * math::pow(u(0,eval_point),static_cast<int>(i)-2);

}

template <class T>
void gsMonomialBasis<T> :: evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
{
    result.resize(size()*(n+1),u.cols());

    for (int deriv_order=0; deriv_order<=n; deriv_order++)
    {
        // Basis functions whose derivative is already zero
        for(int p=0; p<deriv_order && p<=m_p; p++)
            result.row(deriv_order*size()+p).setZero();

        // Basis function whose derivative is constant
        if(deriv_order<=m_p)
            result.row(deriv_order*size()+deriv_order).setConstant(factorial(deriv_order));

        // Remaining basis functions whose derivative is not constant
        for(int p=deriv_order+1; p<=m_p; p++)
            for(int eval_point=0; eval_point<u.cols(); eval_point++)
                result(deriv_order*size()+p,eval_point)=(real_t)p/(p-deriv_order)*result(deriv_order*size()+p-1,eval_point)*u(0,eval_point);
    }
}

template <class T>
void gsMonomialBasis<T> :: evalAllDersSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
{
    GISMO_ENSURE(static_cast<int>(i)<=m_p,"Error: There is no i-th basis function in the basis."<<std::endl<<
                 "Number of basis functions in the basis is "<<size()<<".");

    result.resize(n+1, u.cols());

    // Compute function values
    for(int eval_point=0; eval_point<u.cols(); eval_point++)
        result(0,eval_point) = math::pow(u(0,eval_point),static_cast<int>(i));

    // Compute derivatives up to order n
    for (int deriv_order=1; deriv_order<=n; deriv_order++)
    {
        if(static_cast<int>(i)<deriv_order)
            result.row(deriv_order).setZero();
        else
            for(int eval_point=0; eval_point<u.cols(); eval_point++)
                result(deriv_order,eval_point)=(i+1-deriv_order)*result(deriv_order-1,eval_point)/u(0,eval_point);
    }

}

template <class T>
void gsMonomialBasis<T> :: evalDerSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
{
    GISMO_ENSURE(static_cast<int>(i)<=m_p,"Error: There is no i-th basis function in the basis."<<std::endl<<
                 "Number of basis functions in the basis is "<<size()<<".");

    result.resize(1, u.cols());

    // Compute factor resulting from derivation
    real_t factor=1;
    for(int j=0; j<=n-1; j++)
        factor*=i-j;

    if((int)i<n)
        result.setZero();
    else
        for(int eval_point=0; eval_point<u.cols(); eval_point++)
            result(0,eval_point)=factor * math::pow(u(0,eval_point),static_cast<int>(i)-n);

}


template <class T>
void gsMonomialBasis<T> :: evalFunc_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const
{
    result.resize(coefs.cols(),u.cols());
    int n_coefs=coefs.rows();

    // Loop over the evaluation points
    for(int eval_point=0; eval_point<u.cols(); eval_point++)
    {
        // Horner Scheme
        result.col(eval_point)=coefs.row(n_coefs-1).transpose();
        for(int i=n_coefs-2; i>=0; i--)
            result.col(eval_point)=result.col(eval_point)*u.col(eval_point)+coefs.row(i).transpose();
    }
}



}

