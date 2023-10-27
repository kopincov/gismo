/** @file gsMultiplicativeOp.hpp

    @brief Allows to set up multiplicative Schwarz type preconditioners

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, J. Sogn
*/

namespace gismo
{

template<typename T>
void gsMultiplicativeOp<T>::step(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_underlying->rows() == x.rows() && x.rows() == f.rows()
               && m_underlying->cols() == m_underlying->rows() && x.cols() == f.cols(),
               "Dimensions do not match: "
               <<m_underlying->rows()<<"=="<<x.rows()<<"&&"<<x.rows()<<"=="<<f.rows()
               <<"&&"<<m_underlying->cols()<<"=="<<m_underlying->rows()<<"&&"<<x.cols()<<"=="<<f.cols());

    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    const size_t sz = m_transfers.size();
    for (size_t i = 0; i < sz; ++i)
    {
        // We hope that the following multiplications are fast due to RowMajor
        m_res[i] = m_transfers[i].transpose() * f - m_upds[i] * x;
        m_ops[i]->apply( m_res[i], m_p[i] );
        x += m_transfers[i] * m_p[i];
    }

}



template<typename T>
void gsMultiplicativeOp<T>::stepT(const gsMatrix<T>& f, gsMatrix<T>& x) const
{
    GISMO_ASSERT( m_underlying->rows() == x.rows() && x.rows() == f.rows()
               && m_underlying->cols() == m_underlying->rows() && x.cols() == f.cols(),
                 "Dimensions do not match.");

    GISMO_ASSERT( f.cols() == 1, "This operator is only implemented for a single right-hand side." );

    const size_t sz = m_transfers.size();
    for (size_t i = sz-1; i != (size_t)-1; --i)
    {
        // We hope that the following multiplications are fast due to RowMajor
        m_res[i] = m_transfers[i].transpose() * f - m_upds[i] * x;
        m_ops[i]->apply( m_res[i], m_p[i] );
        x += m_transfers[i] * m_p[i];
    }    

}


} // namespace gismo
