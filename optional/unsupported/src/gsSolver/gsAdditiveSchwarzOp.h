/** @file gsAdditiveSchwarzOp.h

    @brief Provides the application of gsLinearOperators to subspaces 

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/
#pragma once

#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Class for representing an additive Schwarz type operator
///
/// The operator the sum of objects of type gsLinearOperator, pre- and post-multiplied with the given sparse transfer matrix
///
/// \ingroup Solver
template<typename T = real_t>
class gsAdditiveSchwarzOp GISMO_FINAL : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsAdditiveSchwarzOp
    typedef memory::shared_ptr<gsAdditiveSchwarzOp> Ptr;

    /// Unique pointer for gsAdditiveSchwarzOp
    typedef memory::unique_ptr<gsAdditiveSchwarzOp> uPtr;
    
    /// Base class
    typedef typename gsLinearOperator<T>::Ptr BasePtr;
    
    /// Matrix type
    typedef gsSparseMatrix<T> SparseMatr;

    /// Empty constructor. To be filled with addSubspace()
    gsAdditiveSchwarzOp() : m_sz(0), m_transfs(0), m_ops(0) {}
    
    /// Constructor taking vectors of transfer matrices and local operators
    gsAdditiveSchwarzOp( const std::vector<SparseMatr>& transfs, const std::vector<BasePtr>& ops )
    : m_sz(ops.size()), m_transfs(transfs), m_ops(ops)
    {
        for( index_t i=0; i<m_sz; ++i )
        {
            GISMO_ASSERT ( transfs[i].cols() == ops[i]->rows() && ops[i]->rows() == ops[i]->cols(), "Dimensions for operator do not fit." );
            GISMO_ASSERT ( transfs[i].rows() == transfs[0].rows(), "Dimensions of operators do not fit to each other." );
        }   
    }
    
    /// Make command returning a smart pointer to an empty operator. To be filled with addSubspace()
    static uPtr make() { return memory::make_unique( new gsAdditiveSchwarzOp() ); }
    
    /// Make command returning a smart pointer taking vectors of transfer matrices and local operators
    static uPtr make( const std::vector<SparseMatr>& transfs, const std::vector<BasePtr>& ops ) 
    { return memory::make_unique( new gsAdditiveSchwarzOp( transfs, ops ) ); }
   
    /// Add a subspace, specified by a transfer matrix, and a local operator for that space
    void addSubspace( const SparseMatr& transf, const BasePtr& op )
    {
        GISMO_ASSERT ( transf.cols() == op->rows() && op->rows() == op->cols(), "Dimensions for operator do not fit." );
        GISMO_ASSERT ( m_sz == 0 || transf.rows() == m_transfs[0].rows(), "Dimensions of operators do not fit to each other." );
        
        ++m_sz;
        m_transfs.push_back( transf );
        m_ops.push_back( op );
    }
    
    virtual void apply( const gsMatrix<T> & input, gsMatrix<T> & result ) const
    {
        // Here, we could make a permanently allocated vector
        gsMatrix<T> tmp;

        GISMO_ASSERT ( m_sz>0, "gsAdditiveSchwarzOp does not work for 0 operators." );

        for (index_t i=0; i<m_sz; ++i)
        {
            m_ops[i]->apply(m_transfs[i].transpose()*input,tmp);
            if ( i==0 )
                result  = m_transfs[0] * tmp;
            else
                result += m_transfs[i] * tmp;
        }
    }

    virtual index_t rows() const
    {
        GISMO_ASSERT ( m_sz>0, "gsAdditiveSchwarzOp does not work for 0 operators." );
        return m_transfs[0].rows();
    }
    
    virtual index_t cols() const
    {
        GISMO_ASSERT ( m_sz>0, "gsAdditiveSchwarzOp does not work for 0 operators." );
        return m_transfs[0].rows();
    }

private:

    index_t m_sz;
    std::vector<SparseMatr> m_transfs;
    std::vector<BasePtr> m_ops;
    
};

}

