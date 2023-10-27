/** @file gsPowerIteration.h

    @brief Provides a method to count the number of calls and runtimes of an operator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/
#pragma once

#include <gsCore/gsExport.h>
#include <gsSolver/gsLinearOperator.h>
#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>
#include <gsMpi/gsMpiComm.h>
#include <gsIETI/gsDistributedOperator.h>

#endif
namespace gismo
{

    template <typename T>
    inline T powerIteration( const typename gsLinearOperator<T>::Ptr& A, const typename gsLinearOperator<T>::Ptr& B, index_t N = 10)
    {
        GISMO_ASSERT( A->rows() == A->cols() && B->rows() == B->cols() && A->rows() == B->rows(), "Dimensions do not argee");

        gsMatrix<T> x, y, z;
        z.setRandom(A->rows(),1);

        for (index_t i = 0; i < N; ++i)
        {
            z.swap(x);
            x /= math::sqrt( x.col(0).dot(x.col(0)) );
            A->apply(x, y);
            B->apply(y, z);
        }

        return y.col(0).dot(z.col(0)) / x.col(0).dot(y.col(0));

    }

#ifdef GISMO_WITH_MPI

    // A: Acc -> Dist, B: Dist -> Acc  (like matrix and Prec)
    template <typename T>
    inline T powerIteration_MPI( const typename  gsDistributedOperator<T>::Ptr& A, const typename gsDistributedOperator<T>::Ptr& B, gsMpiComm comm, index_t N = 10)
    {
        GISMO_ASSERT( A->rows() == A->cols() && B->rows() == B->cols() && A->rows() == B->rows(), "Dimensions do not argee");

        gsMatrix<T> x_acc,x_,y_, z_,z_acc;
        z_.setRandom(A->rows(),1);
        A->accumulate(z_,z_acc);

        T scalarP;
        for (index_t i = 0; i < N; ++i)
        {
            z_acc.swap(x_acc);
            A->distribute(x_acc,x_);
            scalarP = x_acc.col(0).dot(x_.col(0));
            x_acc /= math::sqrt( comm.sum(scalarP) );
            A->apply(x_acc, y_);
            B->apply(y_, z_acc);
        }

        T scalars[2];
        scalars[0] =  y_.col(0).dot(z_acc.col(0));
        scalars[1] =  x_acc.col(0).dot(y_.col(0));
        comm.sum(scalars,2);

        return scalars[0] / scalars[1];

    }
    template <typename T>
    inline T powerIteration_MPI( const typename  gsDistributedOperator<T>::Ptr& A, gsMpiComm comm, index_t N = 10)
    {
        GISMO_ASSERT( A->rows() == A->cols(), "Dimensions do not argee");

        gsMatrix<T> x,x_, y,y_;
        y_.setRandom(A->rows(),1);
        A->accumulate(y_,y);

        T scalarP;
        for (index_t i = 0; i < N; ++i)
        {
            y.swap(x);
            y_.swap(x_);
            scalarP = x.col(0).dot(x_.col(0));
            x /= math::sqrt( comm.sum(scalarP) );
            A->apply(x, y_);
            A->accumulate(y_,y);
        }

        T scalars[2];
        scalars[0] =  y.col(0).dot(y_.col(0));
        scalars[1] =  x.col(0).dot(y_.col(0));
        comm.sum(scalars,2);

        return scalars[0] / scalars[1];

    }

#endif

}

