#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI
#include <gsIETI/gsParallelCG.h>
#include <gsSolver/gsLanczosMatrix.h>
namespace gismo {

template<class T>
bool gsParallelCG<T>::initIteration(const VectorType& rhs, VectorType& x0)
{
    GISMO_ASSERT(rhs.cols()== 1, "Implemented only for single column right hand side matrix");


    if(m_calcEigenvals)
    {
        delta.clear();
        delta.resize(1,0);
        delta.reserve(m_max_iters / 3);

        gamma.clear();
        gamma.reserve(m_max_iters / 3);

        m_eigsAreCalculated =  true;
    }

    if (Base::initIteration(rhs,x0))
        return true;


    int n = m_mat->cols();
    int m = 1; // == rhs.cols();

    z.resize(n,m);
    z_acc.resize(n,m);
    tmp.resize(n,m);
    p.resize(n,m);
    p_acc.resize(n,m);

    //since rhs is distributed, it may happen that it is not properly distributed, i.e, some entries have large values,
    //which disapper when the vector is accumulated. This cancellation phenomenon leads to instabilities. One way to overcome this
    //is to accumulate the vector and distribute it in order to average the values. Here, we anyway need the accumulated rhs for the
    //norm, hence, we get the distributed rhs for free (no additional communication)
    VectorType rhs_acc, rhs_dist;
    accumulate(rhs,rhs_acc,*m_mat);
    m_rhs_norm = math::sqrt(dot(rhs_acc,rhs));
    distribute(rhs_acc, rhs_dist,*m_mat);


    // gsInfo<<"rank: "<<m_comm.rank()<<" size n="<<n<<"\n";
    //   gsInfo<<"rank: "<<m_comm.rank()<<" x0="<<x0.transpose()<<"\n";

    //This triggers the MPI communications for the operator

    m_mat->apply(x0,residual);  //apply the system matrix

    //  gsInfo<<"rank: "<<m_comm.rank()<<" tmp2: "<<tmp2.transpose()<<"\n";
    //  gsInfo<<"rank: "<<m_comm.rank()<<" rhs: "<<rhs.transpose()<<"\n";

    residual = rhs_dist - residual; //initial residual

    // gsInfo<<"rank: "<<m_comm.rank()<<" before precond: \n";
    //This triggers the MPI communications for the Preconditioner

    // gsInfo<<"rank: "<<m_comm.rank()<<" res: "<<residual.transpose()<<"\n";
    // gsInfo<<"rank: "<<m_comm.rank()<<" res_acc: "<<residual_acc.transpose()<<"\n";

    // postAccumulate(*m_precond);
    m_precond->apply(residual, p_acc);      //initial search direction
    //  gsInfo<<"rank: "<<m_comm.rank()<<" p: "<<p.transpose()<<"\n";

    //startAccumulate(p,*m_precond);

    //The dot function triggers the MPI_ALLREDUCE to calculate the global scalar product
    accumulate(residual,residual_acc, *m_mat);

    T res[2];
    dot(residual,p_acc,residual, residual_acc, res);
    absNew = Eigen::numext::real(res[0]);
    rhsNorm2 =Eigen::numext::real(res[1]);

    m_error = math::sqrt(rhsNorm2) / m_rhs_norm;
    //   absNew = Eigen::numext::real(dot(p,residual_acc));  // the square of the absolute value of r scaled by invM
    //   rhsNorm2 = Eigen::numext::real(dot(residual, residual_acc));

    //  gsInfo<<"rank: "<<m_comm.rank()<<" absNew: "<<absNew<<"  rhsNorm2: "<<rhsNorm2 <<"\n ";
    if (m_error <= m_tol)
        return true;

    // residualNorm2 = 0;
    // threshold = m_tol*m_tol*rhsNorm2;

    //gsInfo<<"Tol res init: "<<rhsNorm2<<std::endl;
    //gsInfo<<"Tol energy init : "<<absNew<<std::endl;
    //gsInfo<<"Tol res: "<<threshold<<std::endl;
    //gsInfo<<"Tol energy: "<<m_tol*m_tol*absNew<<std::endl;


    //   finishAccumulate(p_acc,*m_precond);
    //   gsInfo<<"rank: "<<m_comm.rank()<<" p_acc: "<<p_acc.transpose()<<"\n";

    return false;
}

template<class T>
bool gsParallelCG<T>::step( VectorType& x )
{
    //This triggers the MPI communications for the operaor
    // postAccumulate(*m_mat);
    // gsInfo<<"rank: "<<m_comm.rank()<<" p: "<<p_acc.transpose()<<"\n";
    m_mat->apply(p_acc,tmp); //apply system matrix

    //The dot function triggers the MPI_ALLREDUCE to calculate the global scalar product
    T alpha = absNew / dot(p_acc,tmp);   // the amount we travel on dir
    if(m_calcEigenvals)
        delta.back()+=(1./alpha);


    // gsInfo<<"rank: "<<m_comm.rank()<<" tmp: "<<tmp.transpose()<<"\n";
    // gsInfo<<"rank: "<<m_comm.rank()<<" alpha: "<<alpha<<"  dot(t,temp): "<<dot(p_acc,tmp) <<"\n ";

    residual -= alpha * tmp;              // update residual
    //  startAccumulate(residual,*m_mat);

    x += alpha * p_acc;                       // update solution


    //  gsInfo<<"rank: "<<m_comm.rank()<<" residual: "<<residual.transpose()<<"\n";

    //  gsInfo<<"rank: "<<m_comm.rank()<<" absNew: "<<absNew<<"  rhsNorm2: "<<residualNorm2 <<"\n ";

    //  postAccumulate(*m_precond);  //< for z
    m_precond->apply(residual, z_acc);          // approximately solve for "A z = residual"

    //   startAccumulate(z,*m_precond);
    //  gsInfo<<"rank: "<<m_comm.rank()<<" z: "<<z.transpose()<<"\n";

    T absOld = absNew;

    //The dot function triggers the MPI_ALLREDUCE to calculate the global scalar product
    accumulate(residual,residual_acc,*m_mat);
    //  finishAccumulate(residual_acc,*m_mat);
    T res[2];
    dot(residual,z_acc,residual, residual_acc, res);
    absNew = Eigen::numext::real(res[0]);
    residualNorm2 =Eigen::numext::real(res[1]);
    /*
    absNew = Eigen::numext::real(dot(residual_acc,z));     // update the absolute value of r
    residualNorm2 =Eigen::numext::real(dot(residual, residual_acc));
*/


    T beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
    //   finishAccumulate(z_acc,*m_precond);
    //The dot function triggers the MPI_ALLREDUCE to calculate the global scalar product
    m_error = math::sqrt(residualNorm2)/m_rhs_norm;
    if( m_error < m_tol)
        return true;

    p_acc = z_acc + beta * p_acc;                             // update search direction

    // gsInfo<<"rank: "<<m_comm.rank()<<" beta: "<<beta<<"  dot(residual,p): "<<absNew <<"\n ";
    // gsInfo<<"rank: "<<m_comm.rank()<<" p: "<<p.transpose()<<"\n";

    if(m_calcEigenvals)
    {
        gamma.push_back(-math::sqrt(beta)/alpha);
        delta.push_back(beta/alpha);
    }
    return false;
}

template<class T>
T gsParallelCG<T>::dot(const VectorType &a, const VectorType &b)
{
    T res= a.col(0).dot(b.col(0));
    return m_comm.sum(res);
}

template<typename T>
void gsParallelCG<T>::dot(const VectorType &a, const VectorType &b,
                          const VectorType &c, const VectorType &d, T* res)
{
    res[0] = a.col(0).dot(b.col(0));
    res[1] = c.col(0).dot(d.col(0));
    m_comm.sum(res,2);
}

/*
template<class T>
T gsParallelCG<T>::norm(const VectorType &a, const gsLinearOperator& precond)
{
    VectorType b = a;
    m_comm.sum(b.data(),b.cols()*b.rows());
  //  gsInfo<<"rank: "<<m_comm.rank()<<" b: "<<b.transpose()<<"\n";
    return b.squaredNorm();
}
*/

template<class T>
T gsParallelCG<T>::getConditionNumber()
{
    GISMO_ASSERT(m_eigsAreCalculated,"No data for eigenvalues was collected, call setCalcEigenvalues(true) and solve with an arbitrary right hand side");
    gsLanczosMatrix<T> L(gamma,delta);

    return L.maxEigenvalue()/L.minEigenvalue();
}

template<class T>
void gsParallelCG<T>::getEigenvalues(gsMatrix<T>& eigs )
{
    GISMO_ASSERT(m_eigsAreCalculated,"No data for eigenvalues was collected, call setCalcEigenvalues(true) and solve with an arbitrary right hand side");

    gsLanczosMatrix<T> LM(gamma,delta);

    gsSparseMatrix<T> L = LM.matrix();
    //LM.matrixForm(L);
    //  gsInfo<<"Lancos:\n "<<L.toDense()<<"\n";
    //there is probably a better option...
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,Dynamic,Dynamic> > eigensolver(L.toDense());

    eigs = eigensolver.eigenvalues();
}



}

#endif
