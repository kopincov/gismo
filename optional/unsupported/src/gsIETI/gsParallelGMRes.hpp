#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI
#include <gsIETI/gsParallelGMRes.h>

namespace gismo {

template<typename T>
T gsParallelGMRes<T>::dot(const VectorType &a, const VectorType &b)
{
    T res= a.col(0).dot(b.col(0));
    return m_comm.sum(res);
}

template<typename T>
void gsParallelGMRes<T>::dot(const VectorType &a, const VectorType &b, const VectorType &c, const VectorType &d, T* res)
{
    res[0] = a.col(0).dot(b.col(0));
    res[1] = c.col(0).dot(d.col(0));
    m_comm.sum(res,2);
}

template<typename T>
bool gsParallelGMRes<T>::initIteration( const VectorType& rhs, VectorType& x )
{
    if (Base::initIteration(rhs,x))
        return true;

    //since rhs is distributed, it may happen that it is not properly distributed, i.e, some entries have large values,
    //which disapper when the vector is accumulated. This cancellation phenomenon leads to instabilities. One way to overcome this
    //is to accumulate the vector and distribute it in order to average the values. Here, we anyway need the accumulated rhs for the
    //norm, hence, we get the distributed rhs for free (no additional communication)
    VectorType rhs_acc, rhs_dist;
    accumulate(rhs,rhs_acc,*m_mat);
    distribute(rhs_acc,rhs_dist,*m_mat);

    m_mat->apply(x,tmp);
    tmp = rhs_dist - tmp;

    m_precond->apply(tmp, residual_acc);
    distribute(residual_acc,residual,*m_precond);

   // beta = residual.norm(); // This is  ||r||


    //real_t rhsNorm2 = m_rhs_norm*m_rhs_norm;
    T res[2];
    dot(residual_acc,residual,rhs,rhs_acc,res);
    beta = math::sqrt(res[0]);
    m_rhs_norm = math::sqrt(res[1]);
    m_error = beta/m_rhs_norm;
    if(m_error < m_tol)
        return true;

    v_acc.resize(m_max_iters+1,gsMatrix<T>::Zero(x.rows(),1));
    v_acc[0] = give(residual_acc/beta);
    g.setZero(m_max_iters+1,1);
    g(0,0) = beta;
    Omega = gsMatrix<T>::Identity(2, 2);
    Omega_prev = gsMatrix<T>::Identity(1, 1);

    H.setZero(m_max_iters+1,m_max_iters);
    //h_tmp.setZero(m_max_iters,1);

    return false;
}

template<typename T>
void gsParallelGMRes<T>::finalizeIteration( VectorType& x )
{
   // gsInfo<<"rank: "<<m_comm.rank()<<"  doing now finalizing GMRes \n";
    //Remove last row of H and g
    //H.resize(m_num_iter,m_num_iter);
  //  H = H_prev.block(0,0,m_num_iter,m_num_iter);
 //   g_tmp.resize(m_num_iter,1);
 //   g_tmp = g.block(0,0,m_num_iter,1);

    //Solve H*y = g;
    //solveUpperTriangular(H, g_tmp);
    y = H.block(0,0,m_num_iter,m_num_iter).template triangularView<Eigen::Upper>().solve(g.block(0,0,m_num_iter,1));

    //Create the matrix from the column matrix in v.
    gsMatrix<T> V(m_mat->rows(),m_num_iter);
    for (index_t k = 0; k< m_num_iter; ++k)
    {
        V.col(k) = v_acc[k];
    }
    //Update solution
    x += V*y;

    // cleanup temporaries
    tmp.clear();
    g.clear();
    g_tmp.clear();
    h_tmp.clear();
    y.clear();
    w.clear();
    residual.clear();
    H_prev.clear();
    H.clear();
    Omega.clear();
    Omega_prev.clear();
    Omega_tmp.clear();
    Omega_prev_tmp.clear();
    v_acc.clear();
}

template<typename T>
bool gsParallelGMRes<T>::step( VectorType& x )
{


    GISMO_UNUSED(x); // The iterate x is never updated! Use finalizeIteration to obtain x.
    const index_t k = m_num_iter-1;
     //   gsInfo<<"rank: "<<m_comm.rank()<<" \t start a new GMRes iterations ("<<k<<")"<<"\n";
    /*
    H.setZero(k+2,k+1);
    h_tmp.setZero(k+2,1);

    if (k != 0)
    {
        H.block(0,0,k+1,k) = H_prev;
    }
    */
    gsMatrix<T>& w_acc = v_acc[k+1];
    Omega = gsMatrix<T>::Identity(k+2, k+2);
    m_mat->apply(v_acc[k],tmp);
    m_precond->apply(tmp, w_acc);
    distribute(w_acc,w,*m_precond);

    for (index_t i = 0; i< k+1; ++i)
    {
        //h_tmp(i,0) = (w.transpose()*v[i]).value(); //Typo h_l,k
        H(i,k) = dot(w,v_acc[i]); //v_acc is acc, w is dist
        w_acc = w_acc - H(i,k)*v_acc[i]; //w_acc,v_acc is acc,
        distribute(w_acc,w,*m_precond);
    }
    H(k+1,k) = math::sqrt(dot(w_acc,w)); // needs w_acc, w_dist

  //  if (math::abs(h_tmp(k+1,0)) < 1e-16) //If exact solution
  //      return true;
    v_acc[k+1]/=H(k+1,k);

    //v_acc.push_back(give(w_acc/H(k+1,k)));

    // Omega_prev is k+1 x k+1 (the extra rigfht lower entry would be anyway 1)
    H.block(0,k,k+1,1) = Omega_prev*H.block(0,k,k+1,1);
    //H.block(0,k,k+2,1) = h_tmp;

    //Find coef in rotation matrix
    T sk = H(k+1,k)/(math::sqrt(H(k,k)*H(k,k) + H(k+1,k)*H(k+1,k)));
    T ck = H(k,  k)/(math::sqrt(H(k,k)*H(k,k) + H(k+1,k)*H(k+1,k)));
    Omega(k,k)   = ck; Omega(k,k+1)   = sk;
    Omega(k+1,k) =-sk; Omega(k+1,k+1) = ck;

    //Rotate H and g
    H.topLeftCorner(k+2,k+1) = Omega*H.topLeftCorner(k+2,k+1);
   // H_prev = H;
    /*
    g_tmp.setZero(k+2,1);
    if (k != 0)
        g_tmp.block(0,0,k+1,1) = g;
    else
        g_tmp = g;
        */
    g.topRows(k+2) = Omega * g.topRows(k+2);

    T residualNorm2 = g(k+1,0)*g(k+1,0);
    m_error = math::sqrt(residualNorm2) / m_rhs_norm;
    if(m_error < m_tol)
        return true;

    //Resize rotation product

    Omega_prev_tmp.setZero(k+2,k+2);
    Omega_prev_tmp.block(0,0,k+1,k+1) = Omega_prev;
    Omega_prev_tmp(k+1,k+1) = 1;
    /*
    Omega_tmp.setZero(k+3,k+3);
    Omega_tmp.block(0,0,k+2,k+2) = Omega;
    Omega_tmp(k+2,k+2) = 1;
    */

    Omega_prev = Omega*Omega_prev_tmp;
    return false;
}


}

#endif
