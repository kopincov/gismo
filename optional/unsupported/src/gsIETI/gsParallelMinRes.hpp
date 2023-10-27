#include <gsCore/gsConfig.h>
#ifdef GISMO_WITH_MPI
#include <gsIETI/gsParallelMinRes.h>

namespace gismo
{

template<class T>
bool gsParallelMinRes<T>::initIteration( const VectorType& rhs, VectorType& x )
{
    if (Base::initIteration(rhs,x))
        return true;

    int n = m_mat->cols();
    int m = 1; // = rhs.cols();

    //since rhs is distributed, it may happen that it is not properly distributed, i.e, some entries have large values,
    //which disapper when the vector is accumulated. This cancellation phenomenon leads to instabilities. One way to overcome this
    //is to accumulate the vector and distribute it in order to average the values. Here, we anyway need the accumulated rhs for the
    //norm, hence, we get the distributed rhs for free (no additional communication)
    VectorType rhs_acc, rhs_dist;
    accumulate(rhs,rhs_acc,*m_mat);
    //gsInfo<<"rhs: "<<rhs.minCoeff()<<" , "<<rhs.maxCoeff()<<" -- rhs_acc: "<<rhs_acc.minCoeff()<<" , "<<rhs_acc.maxCoeff()<<"\n";
    distribute(rhs_acc,rhs_dist,*m_mat);

    vPrev.setZero(n,m); vNew.setZero(n,m);v_acc.setZero(n,m);
    wPrev_acc.setZero(n,m); w_acc.setZero(n,m); wNew_acc.setZero(n,m);
    if (!m_inexact_residual) { AwPrev.setZero(n,m); Aw.setZero(n,m); AwNew.setZero(n,m); }

    m_mat->apply(x,negResidual);
    negResidual -= rhs_dist;

    v = -negResidual;
    accumulate(v,v_acc,*m_mat);
    m_precond->apply(v, z_acc);


    //  gsDebugVar(rhs.transpose());
    //  gsDebugVar(x.transpose());
    //  gsDebugVar(negResidual.transpose());
    //  gsDebugVar(z_acc.transpose());

    T res[3];
    dot(z_acc,v,v,v_acc,rhs,rhs_acc,res);
    m_rhs_norm = math::sqrt(res[2]);
    if(res[2] <=0) {
         x.setZero(rhs.rows(),rhs.cols()); // for sure zero is a solution
         m_error = 0.;
         return true; // iteration is finished
    }

    m_error = math::sqrt(res[1])/m_rhs_norm;
    if (res[1]<=0 || m_error < m_tol)
        return true;
    //if(m_comm.rank()==0) gsInfo<<"\tit - error: "<<m_num_iter<<" - "<<m_error<<" ("<<m_rhs_norm<<")\n"<<std::flush;

    if(m_comm.rank()==0) if(res[0]<=0) gsInfo<<"WARNING: dot(z_acc,v) = "<<res[0]<<"\n"<<std::flush;

    gammaPrev = 1;
    gamma = math::sqrt(res[0]);
    gammaNew = 1;
    eta = gamma;
    sPrev = 0; s = 0; sNew = 0;
    cPrev = 1; c = 1; cNew = 1;

    return false;
}

template<class T>
bool gsParallelMinRes<T>::step( VectorType& x )
{

    z_acc /= gamma;
    m_mat->apply(z_acc,Az);
    T delta = dot(z_acc,Az);

    vNew = Az - (delta/gamma)*v - (gamma/gammaPrev)*vPrev;
    m_precond->apply(vNew, zNew_acc);

    T t = dot(zNew_acc,vNew);
    gammaNew = math::sqrt(t);

    if(m_comm.rank()==0) if(t<=0) gsInfo<<"WARNING: dot(zNew,vNew) = "<<t<<"\n"<<std::flush;

    const T a0 = c*delta - cPrev*s*gamma;
    const T a1 = math::sqrt(a0*a0 + gammaNew*gammaNew);
    const T a2 = s*delta + cPrev*c*gamma;
    const T a3 = sPrev*gamma;

    /*
    if((isnan(t) || isnan(gammaNew) || isnan(a0) || isnan(a1) || isnan(a2) || isnan(a3) || isnan(cNew)  || isnan(sNew) || isnan(eta)))
      return true;
    */

    cNew = a0/a1;
    sNew = gammaNew/a1;
    wNew_acc = (z_acc - a3*wPrev_acc - a2*w_acc)/a1;
    if (!m_inexact_residual)
        AwNew = (Az - a3*AwPrev - a2*Aw)/a1;
    x += cNew*eta*wNew_acc;
    if (!m_inexact_residual)
        negResidual += cNew*eta*AwNew;

    if (m_inexact_residual)
        m_error =m_error* sNew; // see https://eigen.tuxfamily.org/dox-devel/unsupported/MINRES_8h_source.html
    else
    {
        accumulate(negResidual,v_acc,*m_mat); //reuse the vector v_acc
        T tt = dot(negResidual,v_acc);
        if(tt<=0){ gsInfo<<"WARNING: dot(negResidual,v_acc) = "<<tt<<"\n"<<std::flush; return true;}
        m_error = math::sqrt(tt) / m_rhs_norm;
    }
    //if(m_comm.rank()==0)gsInfo<<"\tit - error: "<<m_num_iter<<" - "<<m_error<<"\n"<<std::flush;

    eta = -sNew*eta;

    /*
    if(gsMpi::worldRank()==0 && (isnan(t) || isnan(gammaNew) || isnan(a0) || isnan(a1) || isnan(a2) || isnan(a3) || isnan(cNew)  || isnan(sNew) || isnan(eta)))
    {
        gsInfo<<"WARNING: nan appeared!\n";
	gsInfo<<"z_acc:\n"<<z_acc.transpose()<<"\n";
	gsInfo<<"Az:\n"<<Az.transpose()<<"\n";
	gsInfo<<"zNew_acc:\n"<<zNew_acc.transpose()<<"\n";
	gsInfo<<"vNew:\n"<<vNew.transpose()<<"\n";
	gsInfo<<"v:\n"<<v.transpose()<<"\n";
	gsInfo<<"vPrev:\n"<<vPrev.transpose()<<"\n";
	gsInfo<<"wNew_acc:\n"<<wNew_acc.transpose()<<"\n";
	gsInfo<<"wPrev_acc:\n"<<wPrev_acc.transpose()<<"\n";
	gsInfo<<"w_acc:\n"<<w_acc.transpose()<<"\n";
	gsInfo<<"delta:\n"<<delta<<"\n";
	gsInfo<<"gamma:\n"<<gamma<<"\n";
	gsInfo<<"t: "<<t<<"\n";
	gsInfo<<"c: "<<c<<"\n";
	gsInfo<<"s: "<<s<<"\n";
	gsInfo<<"cPrev: "<<cPrev<<"\n";
	gsInfo<<"sPrev: "<<sPrev<<"\n";
	gsInfo<<"a0: "<<a0<<"\n";
	gsInfo<<"a1: "<<a1<<"\n";
	gsInfo<<"a2: "<<a2<<"\n";
	gsInfo<<"a3: "<<a3<<"\n";
	gsInfo<<"cNew: "<<cNew<<"\n";
	gsInfo<<"sNew: "<<sNew<<"\n";
	gsInfo<<"eta: "<<eta<<"\n";
	gsInfo<<"m_error: "<<m_error<<"\n";
    }
    */

    // Test for convergence
    if (m_error < m_tol)
        return true;

    //Update variables
    vPrev.swap(v); v.swap(vNew);     // for us the same as: vPrev = v; v = vNew;
    wPrev_acc.swap(w_acc); w_acc.swap(wNew_acc);     // for us the same as: wPrev = w; w = wNew;
    if (!m_inexact_residual)
    { AwPrev.swap(Aw); Aw.swap(AwNew); } // for us the same as: AwPrev = Aw; Aw = AwNew;
    z_acc.swap(zNew_acc);                    // for us the same as: z = zNew;
    gammaPrev = gamma; gamma = gammaNew;
    sPrev = s; s = sNew;
    cPrev = c; c = cNew;



    return false;
}

template<class T>
void gsParallelMinRes<T>::finalizeIteration( VectorType& x )
{
    GISMO_UNUSED(x);
    // cleanup temporaries
    negResidual.clear();
    vPrev.clear(); v.clear(); vNew.clear(); v_acc.clear();
    wPrev_acc.clear(); w_acc.clear(); wNew_acc.clear();
    AwPrev.clear(); Aw.clear(); AwNew.clear();
    zNew_acc.clear(); z_acc.clear(); Az.clear();
}


template<class T>
T gsParallelMinRes<T>::dot(const VectorType &a, const VectorType &b)
{
    T res= a.col(0).dot(b.col(0));
    return m_comm.sum(res);
}

template<class T>
void gsParallelMinRes<T>::dot(const VectorType &a, const VectorType &b,
                              const VectorType &c, const VectorType &d, T* res)
{
    res[0] = a.col(0).dot(b.col(0));
    res[1] = c.col(0).dot(d.col(0));
    m_comm.sum(res,2);
}

template<class T>
void gsParallelMinRes<T>::dot(const VectorType &a, const VectorType &b,
                              const VectorType &c, const VectorType &d,
                              const VectorType &e, const VectorType &f, T* res)
{
    res[0] = a.col(0).dot(b.col(0));
    res[1] = c.col(0).dot(d.col(0));
    res[2] = e.col(0).dot(f.col(0));
    m_comm.sum(res,3);
}


}
#endif
