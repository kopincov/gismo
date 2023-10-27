
#pragma once

namespace gismo
{

template<class T>
void gsLagrangeBasis<T>::anchors_into(gsMatrix<T>& result) const
{
    result.resize(1,m_p+1);
    gsMatrix<> res(1,1);
    for(int j = 0; j<=m_p;++j)
    {
        anchor_into(j,res);
        result(0,j)=res(0,0);
    }
}

template<class T>
void gsLagrangeBasis<T>::active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
{
    result.resize(m_p+1,u.cols());
    result.colwise() = gsVector<index_t>::LinSpaced(m_p+1, 0, m_p);
}

template<class T>
gsBasis<T> * gsLagrangeBasis<T>::boundaryBasis_impl(boundary::side const &) const
{
    return new gsLagrangeBasis<T>(this->m_breaks,this->m_start,this->m_end);
}

template<class T>
gsMatrix<T> gsLagrangeBasis<T>::supportInterval(index_t dir) const
{
    gsMatrix<T> res(1,2);
    if(dir==0)
    {
        (*res)(0,0)=m_start;
        (*res)(0,1)=m_end;
    }
    else
    {
        (*res)(0,1)=m_start;
        (*res)(0,0)=m_end;
    }
    return res;
}

template<class T>
void gsLagrangeBasis<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    result.resize(m_p+1,u.cols());
    for (int i=0;i<=m_p;++i)
    {
        gsMatrix<> res(1,u.cols());
        evalSingle_into(i,u,res);
        for(int j=0;j<u.cols();++j)
        {
            result(i,j)=res(0,j);
        }
    }
}

template<class T>
void gsLagrangeBasis<T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    result.resize(m_p+1,u.cols());
    for (int i=0;i<=m_p;++i)
    {
        gsMatrix<> res(1,u.cols());
        derivSingle_into(i,u,res);
        for(int j=0;j<u.cols();++j)
        {
            result(i,j)=res(0,j);
        }
    }
}

template<class T>
void gsLagrangeBasis<T>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    result.resize(m_p+1,u.cols());
    for(int i=0;i<=m_p;i++)
    {
        gsMatrix<> res(1,u.cols());
        deriv2Single_into(i,u,res);
        for(int j=0;j<u.cols();j++)
        {
            result(i,j)=res(0,j);
        }
    }
}

template<class T>
void gsLagrangeBasis<T>::evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
{
    result.resize((n+1)*(m_p + 2),u.cols());
    for (int i=0;i<=m_p+1;++i)
    {
        gsMatrix<> res(n+1,u.cols());
        evalAllDersSingle_into(i,u,n,res);
        for(int j=0;j<u.cols();++j)
        {
            for(int k=0;k<n+1;k++)
            {
                result(i*(n+1)+k,j)=res(1,j);
            }
        }
    }
}

template<class T>
void gsLagrangeBasis<T>::evalAllDersSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
{
    result.resize(n+1, u.cols() );
    result.setZero();
    for(int order = 0;order<n+1;order++)
    {
        gsMatrix<> res(1,u.cols());
        evalDerSingle_into(i,u,order,res);
        for(int j=0;j<u.cols();j++)
            result(order,j)=res(0,j);
    }
}

template<class T>
void gsLagrangeBasis<T>::evalDerSingle_into(index_t i, const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
{
    result.resize(1, u.cols() );
    result.setZero();
    if(n<this->degree())
    {
        T fact = _getFactor(i);
        std::vector<int> vec(m_p-n);
        for(int k = 0; k<=m_p-n-1;k++) vec[k]=k;

        do
        {
            std::vector<real_t> mults(u.cols());
            for(int j=0;j<u.cols();j++)
            {
                mults[j]=1;
                for(int k = 0; k<=m_p-n-1;k++)
                {
                    int index = vec[k];
                    if(index>=static_cast<int>(i))
                        index++;
                    mults[j]*=(u(0,j)-m_breaks[index]);
                }
                result(0,j)+=mults[j];
            }
        }
        while(_nextPoint(vec, m_p));
        for(int j=0;j<u.cols();j++)
        {
            if(n>1)
                result(0,j)/=(fact/std::pow(2,n-1));
            else
                result(0,j)/=fact;
        }
    }
    else if(n==this->degree())
    {
        T fact = _getFactor(i);
        for(int j=0;j<u.cols();j++)
        {
            if(n>1)
                result(0,j)=m_p*std::pow(2,n-1)/fact;
            else
                result(0,j)=m_p/fact;
        }
    }
    else
    {
        for(int j=0;j<u.cols();j++)
            result(0,j)=0;
    }
}

template<class T>
void gsLagrangeBasis<T>::reparameterizeToZeroOne()
{
    T old_length = m_end-m_start;
    for(unsigned i=0;i<m_breaks.size();i++)
    {
        m_breaks[i]=(m_breaks[i]-m_start)/old_length;
    }
    m_start = 0;
    m_end = 1;
}

template<class T>
void gsLagrangeBasis<T>::_getTransformationLagrangeMonomial(gsMatrix<T> & result) const
{
    int basis_size = m_breaks.size();
    result.resize(basis_size, basis_size);
    result.setZero();
    for(int i = 0; i < basis_size; i++)
    {
        T fact = _getFactor(i);
        for(int j = 0; j < basis_size-1; j++)
        {
            std::vector<int> vec(m_p-j);
            for(int k = 0; k<m_p-j;k++)
            {
                vec[k]=k;
            }
            real_t sum = 0;
            do
            {
                real_t prod=1;
                for(int k = 0; k<m_p-j;k++)
                {
                    int index = vec[k];
                    if(index>=i)
                        index++;
                    prod*=m_breaks[index];
                }
                sum+=prod;
            }
            while(_nextPoint(vec,basis_size-1));
            result(i,j)=sum;
        }
        result(i,basis_size-1)=1;
        int sign = 1;
        for(int j = basis_size-1; j >= 0; j--)
        {
            result(i,j) *= sign / fact;
            sign *= -1;
        }
    }
    result.transposeInPlace(); // transpose for control points
}

template<class T>
void gsLagrangeBasis<T>::_getTransformationMonomialBezier(gsMatrix<T> & result) const
{
    std::vector<T> breaks = m_breaks;
    int basis_size = breaks.size();
    result.resize(basis_size, basis_size);
    result.setZero();
    int n=this->degree();
    for(int j = 0; j < basis_size; j++)
    {
        for(int i = 0; i <= j; i++)
        {
            result(i,j)=1./binomial(n,j)*binomial(n-i,j-i);
        }
    }
    result.transposeInPlace(); // transpose for control points
}

template<class T>
bool gsLagrangeBasis<T>::_nextPoint(std::vector<int> & vec, int end) const
{
    int n=vec.size(), changed=-1, i = 0;
    for(int j = n-1 ; j>=0 ; j--)
    {
        if(vec[j]+1<end-i){
            vec[j]++;
            changed=j;
            break;
        }
        i++;
    }
    if (changed == -1)
        return false;
    for(int j=changed;j<n-1;++j)
    {
        vec[j+1]=vec[j]+1;
    }
    return true;
}

} // namespace gismo

