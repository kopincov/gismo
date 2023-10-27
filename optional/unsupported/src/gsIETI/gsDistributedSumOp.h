/** @file gsDistributedSumOp.h

    @brief Sum of Distributed Operators


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2018-05-14
*/


#pragma once
#include <gsCore/gsConfig.h>


#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>
#include <gsIETI/gsDistributedOperator.h>
#include <gsSolver/gsSumOp.h>


namespace gismo {

template<typename T>
class gsDistributedSumOp : public gsDistributedOperator<T>
{
    typedef typename gsDistributedOperator<T>::Ptr BasePtr;

public:
    typedef memory::shared_ptr<gsDistributedSumOp> Ptr;
    typedef memory::unique_ptr<gsDistributedSumOp> uPtr;


public:
    /// Empty constructor. To be filled with addOperator()
    gsDistributedSumOp() : m_op(0) {}

    /// Constructor taking a vector of Linear Operators
    gsDistributedSumOp(std::vector<BasePtr> ops)
        : m_op(give(ops)) {}

    /// Constructor taking two Linear Operators
    gsDistributedSumOp(BasePtr op0, BasePtr op1)
        : m_op(op0, op1) {}

    /// Constructor taking three Linear Operators
    gsDistributedSumOp(BasePtr op0, BasePtr op1, BasePtr op2 )
        : m_op(op0, op1, op2) {}


    /// Make command returning a smart pointer
    static uPtr make()
    { return uPtr( new gsDistributedSumOp() ); }

    /// Make command returning a smart pointer
    static uPtr make(std::vector<BasePtr> ops)
    { return uPtr( new gsDistributedSumOp(give(ops)) ); }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1)
    { return uPtr( new gsDistributedSumOp(give(op0),give(op1)) ); }

    /// Make command returning a smart pointer
    static uPtr make(BasePtr op0, BasePtr op1, BasePtr op2)
    { return uPtr( new gsDistributedSumOp(give(op0),give(op1),give(op2)) ); }

    virtual ~gsDistributedSumOp() {}

    /// Add another operator
    void addOperator(BasePtr op)
    {
        m_op.addOperator(op);
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        m_op.apply(input,x);
    }

    index_t rows() const
    {
        return m_op.rows();
    }

    index_t cols() const
    {
        return m_op.cols();
    }

    /// Return a vector of shared pointers to all operators
    const std::vector<typename gsLinearOperator<T>::Ptr >& getOps() const { return m_op.getOps(); }

    const gsSumOp<T>& getSumOp() const { return m_op; }
    gsSumOp<T>& getSumOp() { return m_op; }



    virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
    {
         const std::vector<typename gsLinearOperator<T>::Ptr >& ops = getOps();
         GISMO_ASSERT(!ops.empty(), "distribute does not work for 0 operators");
         static_cast<gsDistributedOperator<T>*>(ops.front().get())->distribute(input,distributed);
    }

    virtual void postAccumulate() const
    {
        const std::vector<typename gsLinearOperator<T>::Ptr >& ops = getOps();
        GISMO_ASSERT(!ops.empty(), "distribute does not work for 0 operators");
        static_cast<gsDistributedOperator<T>*>(ops.front().get())->postAccumulate();
    }
    virtual void startAccumulate(const  gsMatrix<T> & input) const
    {
        const std::vector<typename gsLinearOperator<T>::Ptr >& ops = getOps();
        GISMO_ASSERT(!ops.empty(), "distribute does not work for 0 operators");
        static_cast<gsDistributedOperator<T>*>(ops.front().get())->startAccumulate(input);
    }
    virtual void finishAccumulate(gsMatrix<T> & result) const
    {
        const std::vector<typename gsLinearOperator<T>::Ptr >& ops = getOps();
        GISMO_ASSERT(!ops.empty(), "distribute does not work for 0 operators");
        static_cast<gsDistributedOperator<T>*>(ops.front().get())->finishAccumulate(result);
    }

    private:

    gsSumOp<T> m_op;
};







} // namespace gismo
#endif


