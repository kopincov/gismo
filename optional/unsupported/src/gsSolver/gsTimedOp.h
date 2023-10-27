/** @file gsTimedOp.h

    @brief Provides a method to count the number of calls and runtimes of an operator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/
#pragma once

#include <gsCore/gsExport.h>
#include <gsUtils/gsStopwatch.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Wrapper class for a gsLinerarOperator which measures the time
///
/// \ingroup Solver
template<class T = real_t>
class gsTimedOp : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsTimedOp
    typedef memory::shared_ptr<gsTimedOp> Ptr;

    /// Unique pointer for gsTimedOp
    typedef memory::unique_ptr<gsTimedOp> uPtr;

    /// Base class
    typedef typename gsLinearOperator<T>::Ptr BasePtr;

    /// Constructor taking the name and a shared pointer to an operator
    gsTimedOp(const std::string& name, const BasePtr& op, bool print_statistics_when_destructed=true)
        : m_op(op), m_calls(0), m_cols(0), m_timing(0), m_name(name), m_destr(print_statistics_when_destructed) {}

    /// Make command returning a smart pointer
    static uPtr make(const std::string& name, const BasePtr& op, bool print_statistics_when_destructed=true)
    { return memory::make_unique( new gsTimedOp(name,op,print_statistics_when_destructed) ); }

    /// Destructor. Prints the statistics
    virtual ~gsTimedOp()
    {
        if ( m_destr ) printStatistics();
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
    {
        gsStopwatch time;
        m_op->apply(input, result);
        m_timing += time.stop();
        ++m_calls;
        m_cols += input.cols();
    }

    real_t getTime() const { return m_timing; }  ///< Get the time spent in the operator
    index_t getCalls() const { return m_calls; } ///< Get the number of calls

    /// Print statistics
    void printStatistics() const
    {
        gsInfo << "Operator " << m_name << " was called " << m_calls
            << " times (for " << m_cols << " right-hand-side columns in total) and used " << m_timing << " secs.\n";
        m_destr = false;
    }

    virtual index_t rows() const    { return m_op->rows(); }
    virtual index_t cols() const    { return m_op->cols(); }

private:

    BasePtr m_op;
    mutable index_t m_calls;
    mutable index_t m_cols;
    mutable real_t m_timing;
    std::string m_name;
    mutable bool m_destr;
};


}
