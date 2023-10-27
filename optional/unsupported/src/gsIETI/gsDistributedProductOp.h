/** @file gsDistributedProductOp.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
 **/

#pragma once
#include <gsCore/gsConfig.h>


#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>
#include <gsIETI/gsDistributedOperator.h>
#include <gsSolver/gsProductOp.h>

namespace gismo {


template<class T>
class gsDistributedProductOp : public gsDistributedOperator<T>
{
public:
    /// Shared pointer for gsScaledOp
    typedef memory::shared_ptr<gsDistributedProductOp> Ptr;

    /// Unique pointer for gsScaledOp
    typedef memory::unique_ptr<gsDistributedProductOp> uPtr;

    /// Shared pointer for gsLinearOperator
    typedef typename gsDistributedProductOp<T>::Ptr BasePtr;

    /// Empty constructor. To be filled with addOperator()
    gsDistributedProductOp() : m_op(0) {}

    /// Constructor taking a vector of Linear Operators
    gsDistributedProductOp(std::vector<BasePtr> ops)
        : m_op(give(ops)) {}

    /// Constructor taking two Linear Operators
    gsDistributedProductOp(BasePtr op0, BasePtr op1)
        : m_op(op0, op1) {}

    /// Constructor taking three Linear Operators
    gsDistributedProductOp(BasePtr op0, BasePtr op1, BasePtr op2 )
        : m_op(op0, op1, op2) {}


    /// Make command returning a smart pointer
    static uPtr make()
    { return uPtr( new gsDistributedProductOp() ); }

    /// Make command returning a smart pointer
    static uPtr make(std::vector<BasePtr> ops)
    { return uPtr( new gsDistributedProductOp(give(ops)) ); }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1)
    { return uPtr( new gsDistributedProductOp(give(op0),give(op1)) ); }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1, BasePtr op2)
    { return uPtr( new gsDistributedProductOp(give(op0),give(op1),give(op2)) ); }

    virtual ~gsDistributedProductOp() {}

    /// Add another operator
    void addOperator(BasePtr op)
    {
        m_op.addOperator(op);
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        m_op.apply(input,x);
    }

    ///Returns the number of rows in the preconditioner
    index_t rows() const {return m_op.rows();}

    ///Returns the number of columns in the preconditioner
    index_t cols() const {return m_op.cols();}


    virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
    {
         m_op.getOps().back()->distribute(input,distributed);
    }

    virtual void postAccumulate() const
    {
         m_op.getOps().back()->postAccumulate();
    }
    virtual void startAccumulate(const  gsMatrix<T> & input) const
    {
        m_op.getOps().back()->startAccumulate(input);
    }
    virtual void finishAccumulate(gsMatrix<T> & result) const
    {
        m_op.getOps().back()->finishAccumulate(result);
    }


private:
    gsProductOp<T> m_op;
};


}

#endif
