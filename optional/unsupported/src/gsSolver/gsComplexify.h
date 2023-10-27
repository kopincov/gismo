/** @file gsComplexify.h

    @brief Wrapper classes for Linear operators that are used in combination with complex numbers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#pragma once

#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/**
 * @brief The gsComplexify class
 *
 *  Wraps a linear operator, which is defined on real numbers T, as a linear operator, which is only defined on the corresponding complex number.
 *  This class is usefull, when you want to use an operator on real numbers in combinations with operators defined on complex number,
 *  e.g,  gsBlockOp, gsSumOp...
 *
 *  Given a complex number z = x + iy and an operator A, this class returs: Az = Ax + iAy
 */
template<typename T>
class gsComplexify : public gsLinearOperator<std::complex<T> >
{
public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<gsComplexify> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<gsComplexify> uPtr;


    gsComplexify(typename gsLinearOperator<T>::Ptr op) : m_op(op)
    {
        m_tempOut.setZero(m_op->rows(),1);
    }

    static  uPtr make(typename gsLinearOperator<T>::Ptr op) {return memory::make_unique(new gsComplexify(op));}

    virtual ~gsComplexify() {}

    /**
     * @brief apply the operator on the input vector and store the result in x
     * @param input Input vector
     * @param x     result vector
     */
    virtual void apply(const gsMatrix<std::complex<T> > & input, gsMatrix<std::complex<T> > & x) const
    {
        x.setZero(input.rows(),input.cols());
        if((input.imag().array().abs()>1.e-15).any())
        {
            //m_tempIn = input.imag();
            m_op->apply(input.imag(),m_tempOut);
            x.imag().swap(m_tempOut);
        }
        //m_tempIn = input.real();
        m_op->apply(input.real(), m_tempOut);
        x.real().swap(m_tempOut.real());
    }

    /// Returns the number of rows of the operator
    virtual index_t rows() const {return m_op->rows();}

    /// Returns the number of columns of the operator
    virtual index_t cols() const {return m_op->cols();}

private:
    typename gsLinearOperator<T>::Ptr m_op;
    mutable gsMatrix<T> m_tempOut;
};


/**
 * @brief The gsDeComplexify class
 *
 *  Wraps a linear operator, which is defined on complex numbers, i.e.,  std::complex<T>, as an linear operator
 *  which is only defined on the underlying real type T. This class is usefull, when you know that you linear
 *  returns only real numbers for real inputs, but is composed out operators based on complex numbers,
 *  e.g., used in gsBlockOp, gsSumOp... .
 *
 *  Given a real number x and an operator A, this class returs: Ax = A(x+i0).real()
 */
template<typename T>
class gsDeComplexify : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<gsDeComplexify> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<gsDeComplexify> uPtr;


    gsDeComplexify(typename gsLinearOperator<std::complex<T> >::Ptr op) : m_op(op)
    {
        m_temp.setZero(m_op->rows(),1);
    }

    static uPtr make(typename gsLinearOperator<std::complex<T> >::Ptr op) {return memory::make_unique(new gsDeComplexify(op));}

    virtual ~gsDeComplexify() {}

    /**
     * @brief apply the operator on the input vector and store the result in x
     * @param input Input vector
     * @param x     result vector
     */
    virtual void apply(const gsMatrix<T > & input, gsMatrix<T> & x) const
    {
        m_op->apply(input, m_temp);
        x = m_temp.real();
    }

    /// Returns the number of rows of the operator
    virtual index_t rows() const {return m_op->rows();}

    /// Returns the number of columns of the operator
    virtual index_t cols() const {return m_op->cols();}

private:
    typename gsLinearOperator<std::complex<T> >::Ptr m_op;
    mutable gsMatrix<std::complex<T> > m_temp;
};

}


#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>

#include <gsIETI/gsDistributedOperator.h>

namespace gismo
{

template<typename T>
class gsDistributedComplexify : public gsDistributedOperator<std::complex<T> >
{
public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<gsDistributedComplexify> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<gsDistributedComplexify> uPtr;


    gsDistributedComplexify(typename gsDistributedOperator<T>::Ptr op) : m_op(op)
    {
        m_tempOut.setZero(m_op->rows(),1);
    }

    static  uPtr make(typename gsDistributedOperator<T>::Ptr op) {return memory::make_unique(new gsDistributedComplexify(op));}

    virtual ~gsDistributedComplexify() {}

    virtual void apply(const gsMatrix<std::complex<T> > & input, gsMatrix<std::complex<T> > & x) const
    {
        x.setZero(input.rows(),input.cols());
        if((input.imag().array().abs()>1.e-15).any())
        {
            m_op->apply(input.imag(),m_tempOut);
            x.imag().swap(m_tempOut);
        }
        m_op->apply(input.real(), m_tempOut);
        x.real().swap(m_tempOut.real());
    }

    /// Returns the number of rows of the operator
    virtual index_t rows() const {return m_op->rows();}

    /// Returns the number of columns of the operator
    virtual index_t cols() const {return m_op->cols();}



    virtual void distribute(const gsMatrix<std::complex<T> > & input, gsMatrix<std::complex<T> > & distributed) const
    {
        distributed.setZero(input.rows(),input.cols());
        if((input.imag().array().abs()>1.e-15).any())
        {
            m_op->distribute(input.real(),m_tempOut);
            distributed.imag().swap(m_tempOut);
        }
        m_op->distribute(input.real(), m_tempOut);
        distributed.real().swap(m_tempOut.real());
    }

    virtual void accumulate(const gsMatrix<std::complex<T> > &input, gsMatrix<std::complex<T> > &accumulated) const
    {
        accumulated.setZero(input.rows(),input.cols());
        if((input.imag().array().abs()>1.e-15).any())
        {
            m_op->accumulate(input.real(),m_tempOut);
            accumulated.imag().swap(m_tempOut);
        }
        m_op->accumulate(input.real(), m_tempOut);
        accumulated.real().swap(m_tempOut.real());
    }

    virtual void postAccumulate() const
    {
        m_op->postAccumulate();
    }

    virtual void startAccumulate(const  gsMatrix<std::complex<T> > & input) const
    {
        // GISMO_ERROR("This function cannot properly handly post/start/finish accumulate, use accumulate() instead.");
        m_tempOut.setZero(input.rows(),2*input.cols());
        if((input.imag().array().abs()>1.e-15).any())
            m_tempOut.rightCols(input.cols()) = input.imag();

        m_tempOut.leftCols(input.cols())=input.real();

        m_op->startAccumulate(m_tempOut);
    }

    virtual void finishAccumulate(gsMatrix<std::complex<T> > & result) const
    {
        m_op->finishAccumulate(m_tempOut);
        result.resize(m_tempOut.rows(),m_tempOut.cols()/2);
        result.real()=m_tempOut.leftCols(m_tempOut.cols()/2);
        result.imag()=m_tempOut.rightCols(m_tempOut.cols()/2);
    }

private:
    typename gsDistributedOperator<T>::Ptr m_op;
    mutable gsMatrix<T> m_tempOut;
};

template<typename T>
class gsDistributedDeComplexify : public gsDistributedOperator<T>
{
public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<gsDistributedDeComplexify> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<gsDistributedDeComplexify> uPtr;


    gsDistributedDeComplexify(typename gsDistributedOperator<std::complex<T> >::Ptr op) : m_op(op)
    {
        m_temp.setZero(m_op->rows(),1);
    }

    static uPtr make(typename gsDistributedOperator<std::complex<T> >::Ptr op) {return memory::make_unique(new gsDistributedDeComplexify(op));}

    virtual ~gsDistributedDeComplexify() {}

    virtual void apply(const gsMatrix<T > & input, gsMatrix<T> & x) const
    {
        m_op->apply(input, m_temp);
        x = m_temp.real();
    }

    /// Returns the number of rows of the operator
    virtual index_t rows() const {return m_op->rows();}

    /// Returns the number of columns of the operator
    virtual index_t cols() const {return m_op->cols();}


    virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
    {
        m_op->distribute(input,m_temp);
        distributed.swap(m_temp.real());
    }

    virtual void postAccumulate() const
    {
        m_op->postAccumulate();
    }

    virtual void startAccumulate(const  gsMatrix<T> & input) const
    {
        m_op->startAccumulate(input);
    }

    virtual void finishAccumulate(gsMatrix<T> & result) const
    {
        m_op->finishAccumulate(m_temp);
        result.swap(m_temp.real());
    }


private:
    typename gsDistributedOperator<std::complex<T> >::Ptr m_op;
    mutable gsMatrix<std::complex<T> > m_temp;
};


}

#endif
