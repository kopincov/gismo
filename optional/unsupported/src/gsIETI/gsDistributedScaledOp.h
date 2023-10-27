/** @file gsDistributedScaledOp.h

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
#include <gsSolver/gsLinearOperator.h>

namespace gismo {


template<class T>
class gsDistributedScaledOp : public gsDistributedOperator<T>
{
public:
    /// Shared pointer for gsScaledOp
    typedef memory::shared_ptr<gsDistributedScaledOp> Ptr;

    /// Unique pointer for gsScaledOp
    typedef memory::unique_ptr<gsDistributedScaledOp> uPtr;

    /// Shared pointer for gsLinearOperator
    typedef typename gsDistributedOperator<T>::Ptr BasePtr;

    /// Constructor taking a shared pointer to a linear operator and a scalar
    gsDistributedScaledOp(BasePtr op, T scalar = 1) : m_op(give(op)), m_scalar(scalar)    {}

    /// Make function returning a smart pointer
    static uPtr make(BasePtr op, T scalar = 1)
    { return uPtr( new gsDistributedScaledOp(give(op), scalar) ); }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        m_op->apply(input, x);
        x *= m_scalar;
    }

    ///Returns the number of rows in the preconditioner
    index_t rows() const {return m_op->rows();}

    ///Returns the number of columns in the preconditioner
    index_t cols() const {return m_op->cols();}


    virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
    {
         m_op->distribute(input,distributed);
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
        m_op->finishAccumulate(result);
    }


private:
    const BasePtr m_op;
    const T m_scalar;
};


}

#endif
