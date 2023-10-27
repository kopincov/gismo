/** @file gsLinearCombinationOfFunctionsFunction.h

    @brief Provides declaration of the gsLinearCombinationOfFunctionsFunction class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunction.h>
#include <vector>

namespace gismo
{

/** 
    @brief Class defining a linear combination of functions: sum_i a[i]*f[i]

    \tparam T value type

    \ingroup function
    \ingroup Core
*/

template <class T=real_t>
class gsLinearCombinationOfFunctionsFunction : public gsFunction<T>
{
public:
    /// Shared pointer for gsLinearCombinationOfFunctionsFunction
    typedef memory::shared_ptr< gsLinearCombinationOfFunctionsFunction<T> > Ptr;

    /// Auto pointer for gsFunctionExpr
    typedef memory::unique_ptr< gsLinearCombinationOfFunctionsFunction<T> > uPtr;

    /// Pointer to gsFunction
    typedef typename gsFunction<T>::Ptr BasePtr;

private:
    gsLinearCombinationOfFunctionsFunction() { }

public:

    explicit gsLinearCombinationOfFunctionsFunction(std::vector<T> a, const std::vector<BasePtr>& f)
        : m_a(a), m_f(f)
    {
        GISMO_ASSERT( a.size() == f.size(), "The domain dims must match." );
        for ( unsigned i = 1; i<m_f.size(); ++i )
        {
            GISMO_ASSERT( m_f[i]->domainDim() == m_f[0]->domainDim(), "The domain dims must match." );
            GISMO_ASSERT( m_f[i]->targetDim() == m_f[0]->targetDim(), "The target dims must match." );
        }
    }

    explicit gsLinearCombinationOfFunctionsFunction(T a1, const BasePtr& f1)
        : m_a(1), m_f(1)
    {
        m_a[0] = a1;
        m_f[0] = f1;
    }

    explicit gsLinearCombinationOfFunctionsFunction(T a1, const BasePtr& f1, T a2, const BasePtr& f2)
        : m_a(2), m_f(2)
    {
        GISMO_ASSERT( f1->domainDim() == f2->domainDim(), "The domain dims must match." );
        GISMO_ASSERT( f1->targetDim() == f2->targetDim(), "The target dims must match." );
        m_a[0] = a1;
        m_a[1] = a2;
        m_f[0] = f1;
        m_f[1] = f2;
    }

    static Ptr make(std::vector<T> a, const std::vector<BasePtr>& f)    { return memory::make_shared( new gsLinearCombinationOfFunctionsFunction(a,f) );         }
    static Ptr make(T a1, const BasePtr& f1)                            { return memory::make_shared( new gsLinearCombinationOfFunctionsFunction(a1,f1) );       }
    static Ptr make(T a1, const BasePtr& f1, T a2, const BasePtr& f2)   { return memory::make_shared( new gsLinearCombinationOfFunctionsFunction(a1,f1,a2,f2) ); }

    GISMO_CLONE_FUNCTION(gsLinearCombinationOfFunctionsFunction)

    virtual const gsLinearCombinationOfFunctionsFunction & piece(const index_t k) const { return *this; }

    // Documentation in gsFunction class
    virtual short_t domainDim() const   { return m_f[0]->domainDim(); }

    // Documentation in gsFunction class
    virtual short_t targetDim() const   { return m_f[0]->targetDim(); }

    const gsVector<T> & value() const
    {
        gsVector<T> result = m_a[0] * m_f[0]->value();
        for ( unsigned i = 1; i<m_f.size(); ++i )
            result += m_a[i] * m_f[i]->value();
        return result;
    }

    T value(size_t i) const
    {
        T result = m_a[0] * m_f[0]->value(i);
        for ( unsigned i = 1; i<m_f.size(); ++i )
            result += m_a[i] * m_f[i]->value(i);
        return result;
    }

    // Documentation in gsFunction class
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == domainDim(), "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< domainDim());

        m_f[0]->eval_into(u,result);
        result *= m_a[0];
        gsMatrix<T> tmp(result);
        for ( unsigned i = 1; i<m_f.size(); ++i )
        {
             m_f[i]->eval_into(u,tmp);
             result.noalias() += m_a[i] * tmp;
        }
    }

    // Documentation in gsFunction class
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == domainDim(), "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< domainDim());

        m_f[0]->deriv_into(u,result);
        result *= m_a[0];
        gsMatrix<T> tmp(result);
        for ( unsigned i = 1; i<m_f.size(); ++i )
        {
             m_f[i]->deriv_into(u,tmp);
             result.noalias() += m_a[i] * tmp;
        }
    }

    // Documentation in gsFunction class
    virtual void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == domainDim(), "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< domainDim());

        m_f[0]->deriv2_into(u,result);
        result *= m_a[0];
        gsMatrix<T> tmp(result);
        for ( unsigned i = 1; i<m_f.size(); ++i )
        {
             m_f[i]->deriv2_into(u,tmp);
             result.noalias() += m_a[i] * tmp;
        }
    }

    // Documentation in gsFunction class
    virtual std::ostream &print(std::ostream &os) const
    {
        bool plus = false;
        for ( unsigned i = 0; i<m_f.size(); ++i )
        {
            if( m_a[i] != 0. )
            {
                if( plus ) os << "+";
                if( m_a[i] != 1. ) os << m_a[i] << "*";
                os << *m_f[i];
                plus = true;
            }
        }
        return os; 
    }

private:

    std::vector<T> m_a;
    std::vector<BasePtr> m_f;

};

}
