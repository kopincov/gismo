/** @file gsDistributedOperator.h

    @brief Interface class for gsLinerOperators, which are aware of MPI


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2018-02-22
*/


#pragma once
#include <gsCore/gsConfig.h>


#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>


namespace gismo {

template<typename T>
class gsDistributedOperator : public gsLinearOperator<T>
{
public:
    typedef memory::shared_ptr<gsDistributedOperator> Ptr;
    typedef memory::unique_ptr<gsDistributedOperator> uPtr;

    virtual ~gsDistributedOperator() {}

    virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const = 0;

    virtual void accumulate(const gsMatrix<T>& input, gsMatrix<T> & accumulated) const
    {
        postAccumulate();
        accumulated = input;
        startAccumulate(input);
        finishAccumulate(accumulated);
    }

    virtual void postAccumulate() const = 0;
    virtual void startAccumulate(const  gsMatrix<T> & input) const  = 0;
    virtual void finishAccumulate(gsMatrix<T> & result) const = 0;
};







} // namespace gismo
#endif


