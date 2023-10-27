/** @file gsPreconditionerFactory.h

    @brief DESCRIPTION

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/



#pragma once
#include <gsSolver/gsLinearOperator.h>

namespace gismo {

/// @brief Abstract class for stop criteria for iterative solvers.
///
/// \ingroup Solver
class gsPreconditionerFactory
{
public:
    gsPreconditionerFactory()
        : rightPrec(NULL), leftPrec(NULL)
    {
    }
    virtual ~gsPreconditionerFactory(){}

    virtual void compute(const gsLinearOperator &mat) const = 0;

    mutable gsLinearOperator *rightPrec;
    mutable gsLinearOperator *leftPrec;
};


class gsStaticPreconditionerFactory : public gsPreconditionerFactory
{
public:
    gsStaticPreconditionerFactory(gsLinearOperator *R=NULL, gsLinearOperator *L=NULL)
    {
        rightPrec = R;
        leftPrec  = L;
    }

    ~gsStaticPreconditionerFactory()
    {
        if (rightPrec) delete rightPrec;
        if (leftPrec)  delete leftPrec;
    }

    virtual void compute(const gsLinearOperator &mat) const
    {
        if(!rightPrec) rightPrec = new gsIdentityOp(mat.cols());
        if(!leftPrec)  leftPrec  = new gsIdentityOp(mat.rows());
    }
};



}

