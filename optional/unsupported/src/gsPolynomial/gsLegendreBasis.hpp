
#pragma once

namespace gismo
{

template<class T>
void gsLegendreBasis<T>::active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
{
    result.resize(m_p+1,u.cols());
    result.colwise() = gsVector<index_t>::LinSpaced(m_p+1, 0, m_p);
}

template<class T>
void gsLegendreBasis<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    // P_i = (A*x + B) * P_{i-1} - C * P_{i-2}
    result.resize(m_p+1,u.cols());
    result.row(0).setOnes();

    if ( 0!=m_p )
        result.row(1) = u;

    if (m_p>1) // Bonnet’s recursion
    {
        for (int i=2; i<=m_p; ++i)
        {
            result.row(i) =
                _getA(i) * u.array() * result.row(i-1).array() -
                _getC(i)             * result.row(i-2).array();            
        }
    }
}

template<class T>
void gsLegendreBasis<T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    // P_i' =  A * (P_{i-1} + x * P_{i-1}') + B * P_{i-1}' - C * P_{i-2}'
    result.resize(m_p+1,u.cols());
    result.row(0).setZeros();
    
    if ( 0!=m_p )
        result.row(1).setOnes();

    if (m_p>1)
    {
        gsMatrix<T> val;
        eval_into(u, val);

        for (int i=2; i<=m_p; ++i)
        {
            // todo: inline val
            
            result.row(i) =
                _getA(i) * ( val.row(i-1).array() + u.array() * result.row(i-1).array() ) - 
                _getC(i) * result.row(i-2).array();
        }
    }
}

template<class T>
void gsLegendreBasis<T>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    GISMO_NO_IMPLEMENTATION
}

} // namespace gismo

