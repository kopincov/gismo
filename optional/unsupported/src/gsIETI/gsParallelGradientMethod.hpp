
#include <gsIETI/gsParallelGradientMethod.h>

namespace gismo
{

template<typename T>
bool gsParallelGradientMethod<T>::initIteration( const VectorType& rhs, VectorType& x )
{
    if (Base::initIteration(rhs,x))
        return true;

    int n = m_mat->cols();
    int m = 1;                                                          // == rhs.cols();
    m_tmp.resize(n,m);

    //since rhs is distributed, it may happen that it is not properly distributed, i.e, some entries have large values,
    //which disapper when the vector is accumulated. This cancellation phenomenon leads to instabilities. One way to overcome this
    //is to accumulate the vector and distribute it in order to average the values. Here, we anyway need the accumulated rhs for the
    //norm, hence, we get the distributed rhs for free (no additional communication)
    VectorType rhs_acc,rhs_dist;
    accumulate(rhs,rhs_acc,*m_mat);
    m_rhs_norm = math::sqrt(dot(rhs_acc,rhs));
    distribute(rhs_acc, rhs_dist,*m_mat);

    m_mat->apply(x,m_tmp);

    m_res = rhs - m_tmp;                                                // initial residual


    accumulate(m_res,m_res_acc,*m_mat);

    m_error = math::sqrt(dot(m_res,m_res_acc)) / m_rhs_norm;
    if (m_error < m_tol)
        return true;


    return false;
}


template<typename T>
bool gsParallelGradientMethod<T>::step( VectorType& x )
{
    m_precond->apply(m_res, m_update_acc);

    // Update iterate
    x += m_damping * m_update_acc;

    // Update residual
    m_mat->apply(m_update_acc,m_tmp);
    m_res -= m_damping * m_tmp;
    accumulate(m_res,m_res_acc,*m_mat);

    m_error = math::sqrt(dot(m_res,m_res_acc)) / m_rhs_norm;
    if (m_error < m_tol)
        return true;

    return false;
}


template<typename T>
T gsParallelGradientMethod<T>::dot(const VectorType &a, const VectorType &b)
{
    T res= a.col(0).dot(b.col(0));
    return m_comm.sum(res);
}

}
