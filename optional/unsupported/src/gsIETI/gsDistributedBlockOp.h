/** @file gsDistributedBlockOp.h

    @brief Simple class create a block preconditioner structure.

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
#include <gsSolver/gsBlockOp.h>

namespace gismo
{

template<class T>
class gsDistributedBlockOp : public gsDistributedOperator<T>
{
public:

    /// Shared pointer for gsBlockOp
    typedef memory::shared_ptr< gsDistributedBlockOp<T> > Ptr;

    /// Unique pointer for gsBlockOp
    typedef memory::unique_ptr< gsDistributedBlockOp<T> > uPtr;

    /// Base class
    typedef memory::shared_ptr< gsDistributedOperator<T> > BasePtr;


    /// Constructor. Takes the number of blocks (nRows, nCols). Provide the contents of the blocks with addOperator
    gsDistributedBlockOp(index_t nRows, index_t nCols): m_op(nRows,nCols) { }

    /// Make function returning a smart pointer
    static uPtr make(index_t nRows, index_t nCols)
    { return memory::make_unique( new gsDistributedBlockOp(nRows,nCols) ); }


    /**
     * @brief Add a linear Operator \f$C_{ij}\f$ to the block structure
     * @param row row position in the block operator
     * @param col column position in the block operator
     * @param op shared pointer to the operator
     */
    void addOperator(index_t row, index_t col, const BasePtr& op)
    {
        m_op.addOperator(row,col,op);
    }

    /**
    * @brief Returns the pointer to a linear operator of a specific block (if existent), else a zero pointer
    * @param row row position in the block operator
    * @param col column position in the block operator
    */
    const typename gsLinearOperator<T>::Ptr & getOperator(index_t row, index_t col) const
    {
        return m_op.getOperator(row,col);
    }

    /**
     * @brief Apply the correct segment of the input vector on the preconditioners in the block structure and store the result.
     * @param input  Input vector
     * @param result Result vector
     */
    void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
    {
        m_op.apply(input,result);
    }

    /// Number of row blocks
    index_t rowBlocks() const {return m_op.rowBlocks();}
    /// Number of col blocks
    index_t colBlocks() const {return m_op.colBlocks();}

    index_t rows() const {return m_op.rows();}
    index_t cols() const {return m_op.cols() ;}


    virtual void accumulate(const gsMatrix<T>& input, gsMatrix<T> & accumulated) const
    {
        gsDistributedOperator<T>* op;
        index_t iter = 0;
        gsMatrix<T> result;
        accumulated.resize(input.rows(),input.cols());

        for(index_t r=0; r<m_op.rowBlocks();++r)
        {
            op = NULL;
            //Find the first non-zero operator
            for(index_t c=0; c<m_op.colBlocks();++c)
            {
                if(getOperator(r,c))
                {
                    op = static_cast<gsDistributedOperator<T>*>(getOperator(r,c).get());
                    break;
                }
            }
            GISMO_ASSERT(op!=0, "no operator found");
            op->accumulate(input.block(iter,0,op->rows(),input.cols()),result);
            accumulated.block(iter,0,op->rows(),input.cols()) = result;
            iter+=op->rows();
        }
    }

    //Distributes the image, more optimizations possible
    virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
    {
        gsDistributedOperator<T>* op;
        index_t iter = 0;
        gsMatrix<T> result;
        distributed.resize(input.rows(),input.cols());

        for(index_t r=0; r<m_op.rowBlocks();++r)
        {
            op = NULL;
            //Find the first non-zero operator
            for(index_t c=0; c<m_op.colBlocks();++c)
            {
                if(getOperator(r,c))
                {
                    op = static_cast<gsDistributedOperator<T>*>(getOperator(r,c).get());
                    break;
                }
            }
            GISMO_ASSERT(op!=0, "no operator found");
            op->distribute(input.block(iter,0,op->rows(),input.cols()),result);
            distributed.block(iter,0,op->rows(),input.cols()) = result;
            iter+=op->rows();
        }
    }

    //BUG: if the same parallel Operator is used several times, the send and recv ops will conflict!
    virtual void postAccumulate() const
    {
        gsDistributedOperator<T>* op;
        for(index_t r=0; r<m_op.rowBlocks();++r)
        {
            op = NULL;
            //Find the first non-zero operator
            for(index_t c=0; c<m_op.colBlocks();++c)
            {
                if(getOperator(r,c))
                {
                    op = static_cast<gsDistributedOperator<T>*>(getOperator(r,c).get());
                    break;
                }
            }
            GISMO_ASSERT(op!=0, "no operator found");
            op->postAccumulate();
        }
    }

    //BUG: if the same parallel Operator is used several times, the send and recv ops will conflict!
    virtual void startAccumulate(const  gsMatrix<T> & input) const
    {
        gsDistributedOperator<T>* op;
        index_t iter = 0;

        for(index_t r=0; r<m_op.rowBlocks();++r)
        {
            op = NULL;
            //Find the first non-zero operator
            for(index_t c=0; c<m_op.colBlocks();++c)
            {
                if(getOperator(r,c))
                {
                    op = static_cast<gsDistributedOperator<T>*>(getOperator(r,c).get());
                    break;
                }
            }
            GISMO_ASSERT(op!=0, "no operator found");
            op->startAccumulate(input.block(iter,0,op->rows(),input.cols()));
            iter+=op->rows();
        }
    }

    //BUG: if the same parallel Operator is used several times, the send and recv ops will conflict!
    virtual void finishAccumulate(gsMatrix<T> & result) const
    {
        gsDistributedOperator<T>* op;
        index_t iter = 0;
        gsMatrix<T> out;
        for(index_t r=0; r<m_op.rowBlocks();++r)
        {
            op = NULL;
            //Find the first non-zero operator
            for(index_t c=0; c<m_op.colBlocks();++c)
            {
                if(getOperator(r,c))
                {
                    op = static_cast<gsDistributedOperator<T>*>(getOperator(r,c).get());
                    break;
                }
            }
            GISMO_ASSERT(op!=0, "no operator found");
            op->finishAccumulate(out);
            if(r==0) result.resize(rows(),out.cols());
            result.block(iter,0,op->rows(),out.cols()) = out;
            iter+=op->rows();
        }
    }

private:

    gsBlockOp<T> m_op;
};
}

#endif
