/** @file gsParallelPreconditionerOp.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    Created on: 2018-03-07
*/

#include <gsCore/gsConfig.h>
#pragma once

#ifdef GISMO_WITH_MPI
#include <gsIETI/gsDistributedOperator.h>
#include <gsIETI/gsParallelOperator.h>
#include <gsSolver/gsPreconditioner.h>


namespace gismo {


template<typename T>
class gsParallelPreconditionerOp : public gsDistributedOperator<T>
{
public:
    typedef memory::shared_ptr<gsParallelPreconditionerOp> Ptr;
    typedef memory::unique_ptr<gsParallelPreconditionerOp> uPtr;

public:

    gsParallelPreconditionerOp(typename gsPreconditionerOp<T>::Ptr precOp) : m_precOp(precOp)
    {
        m_distOp = dynamic_cast<gsDistributedOperator<T>* >(m_precOp->underlyingOp().get());
        GISMO_ASSERT(m_distOp != NULL, "the underlying operator has to be from type gsDistributedOperator");
    }

    static gsParallelPreconditionerOp::uPtr make(typename gsPreconditionerOp<T>::Ptr precOp)
    {
        return memory::make_unique(new gsParallelPreconditionerOp<T>(precOp));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        m_precOp->apply(input,x);
    }

    void postAccumulate() const
    {
        m_distOp->postAccumulate();
    }

    void startAccumulate(const  gsMatrix<T> & input) const
    {
        m_distOp->startAccumulate(input);
    }

    void finishAccumulate(gsMatrix<T> & result) const
    {
        m_distOp->finishAccumulate(result);
    }

    void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
    {
        m_distOp->distribute(input,distributed);
    }

    void accumulate(const gsMatrix<T>& input, gsMatrix<T> & accumulated) const
    {
        m_distOp->accumulate(input,accumulated);
    }

    index_t rows() const { return m_precOp->rows(); }
    index_t cols() const { return m_precOp->cols(); }


private:
    typename gsPreconditionerOp<T>::Ptr m_precOp;
    gsDistributedOperator<T>* m_distOp;
};



} // namespace gismo

#endif
