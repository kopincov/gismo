/** @file gsTimedPreconditioner.h

    @brief Allows to check timing of gsPreconditionerOp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsMultiGrid/gsForwardDeclarations.h>
#include <gsSolver/gsPreconditioner.h>
#include <gsUtils/gsStopwatch.h>

namespace gismo
{

/// @brief Wrapper class for a gsPreconditionerOp which measures the time
///
/// \ingroup Solver
template<class T=real_T>
class gsTimedPreconditionerOp : public gsPreconditionerOp<T>
{
public:

    typedef memory::shared_ptr<gsTimedPreconditionerOp> Ptr;
    typedef memory::unique_ptr<gsTimedPreconditionerOp> uPtr;
    typedef typename gsPreconditionerOp<T>::Ptr BasePtr;

    /// Constructor taking the name and a shared pointer to a smoother
    gsTimedPreconditionerOp(const std::string& name, const BasePtr& op, bool print_statistics_when_destructed=true)
        : m_op(op), m_calls(0), m_timing(0), m_destr(print_statistics_when_destructed) {}

    /// Make command returning a smart pointer
    static uPtr make(const std::string& name, const BasePtr& op, bool print_statistics_when_destructed=true)
    { return uPtr( new gsTimedPreconditionerOp(name,op,print_statistics_when_destructed) ); }

    /// Destructor. Prints the statistics
    virtual ~gsTimedPreconditionerOp()
    {
        if ( m_destr ) printStatistics();
    }

    virtual void step(const gsMatrix<>& f, gsMatrix<>& x) const
    {
        gsStopwatch time;
        m_op->step(f, x);
        m_timing += time.stop();
        ++m_calls;
        m_cols += f.cols();
    }

    virtual void stepT(const gsMatrix<>& f, gsMatrix<>& x) const
    {
        gsStopwatch time;
        m_op->stepT(f, x);
        m_timing += time.stop();
        ++m_calls;
        m_cols += f.cols();
    }

    real_t getTime() const { return m_timing; }  ///< Get the time spent in the operator
    index_t getCalls() const { return m_calls; } ///< Get the number of calls

    /// Print statistics
    void printStatistics() const
    {
        gsInfo << "Preconditioner " << m_name << " was called " << m_calls
            << " times (for " << m_cols << " right-hand-side columns in total) and used " << m_timing << " secs.\n";
        m_destr = false;
    }

    index_t rows() const { return m_op->rows(); }
    index_t cols() const { return m_op->cols(); }
    gsLinearOperator<real_t>::Ptr underlyingOp() const { return m_op->underlyingOp(); }

private:

    BasePtr m_op;
    mutable index_t m_calls;
    mutable index_t m_cols;
    mutable real_t m_timing;
    mutable bool m_destr;
};

} // namespace gismo
