/// @file gsOptParamConditionNumber.h
/// Author: Elisabeth Pilgerstorfer
/// Opt...
#pragma once

#include <iostream>
#include <vector>

#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsBSpline.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{
/**
   @brief Class performing ... 
*/

template<class T>
class gsOptParamConditionNumber
{
private:
    typedef std::vector<T> Vector;
    
    /// default constructor
    gsOptParamConditionNumber(){ }
    
public:
    /// constructor
    gsOptParamConditionNumber(gsGeometry<T> const & input, int max_iter = 40) : 
        m_input(input), m_max_iter(max_iter)
        {
            m_output = NULL;
            //m_bound = computebound1D2(knots, current, 10, 2);
        }
    
    /// Destructor
    ~gsOptParamConditionNumber() 
    {
        if ( m_output )
            delete m_output;
    }
    
    
    /// Runs the optimizer
    void run();
    
    gsGeometry<T> & inputGeometry() { return m_input; }
    
    const gsGeometry<T> & result() { return *m_output; }
    
    T lastComputedBound() { return m_bound; }

    void setNumIterations(int max_iter) 
    {
        m_max_iter = max_iter;
    }
    
protected:
    
    /// pointer to the input B-spline object
    const gsGeometry<T> & m_input;

    /// pointer to the output B-spline object
    gsGeometry<T> * m_output;

    /// last computed bound value
    T m_bound;

    /// Options
    int m_max_iter;
    
private:

    /// Runs the optimizer for 1D
    void run1D(gsBSpline<T> const & bsp);

    /// Runs the optimizer for nD
    template <unsigned d>
    void runnD(gsTensorBSpline<d,T> const & bsp)
    { 
        //..            
    }

    // computes the lp-norm of a vector
    // input: vector, T number p
    // output T number that is the lp-norm
    T discretelpnorm(const Vector& vec,T p);

// computes the bound of the condition number in the 1D case, (integral of G' with gauss quadrature)
// input: knot vector, matrix with control points, T numbers p1 and p2 for continuous and discrete lp-norms
// output: T number that is the bound of the condition number of the stiffness matrix
    T computebound1D2(const gsKnotVector<T>& KV,gsMatrix<T> &coefs0, T p1, T p2);

// computes the bound of the condition number in the 2D case
// input: knot vectors KV1, KV2, matrix with control points coefs2D
// output: T number that is the bound of the condition number of the stiffness matrix in the 2D case
    T computebound2d(const gsKnotVector<T>& KV1,const gsKnotVector<T>& KV2,gsMatrix<T> &coefs2D);

// computes the gradient of the bound function
// input: knot vector KV, control points coefs0, T numbers p1 and p2 for continuous and discrete lp-norms
// output: vector of Ts that is the gradient of the bound-function
    Vector computegradbound2(const gsKnotVector<T>& KV,gsMatrix<T> &coefs0,T p1,T p2);

// computes the gradient of the 2D bound function
// input: knot vectors KV1 and KV2, matirx with control points coefs2D
// output: vector with gradient (if gradient needed as matrix, use function transvecmatgradbound
    Vector computegradbound2d(const gsKnotVector<T>& KV1,const gsKnotVector<T>& KV2,gsMatrix<T> &coefs2D);

// computes one step of the gradient method to compute the bound
// input: knot vector KV, control points coefs0, T numbers p1 and p2 for continuous and discrete lp-norms
// output: matrix with new control points
    void loopstep1D(const gsKnotVector<T>& KV,gsMatrix<T> &coefs0, T p1, T p2);


// computes one step of the gradient method to compute the bound
// input: knot vectors KV1 and KV2, control points coefs2D, vector coefschange to specify which control points can be changed
// output: new control points
    gsMatrix<T> loopstep2d(const gsKnotVector<T>& KV1, const gsKnotVector<T>& KV2, 
                           gsMatrix<T> &coefs2D,std::vector<int> coefschange);


    int coefsinrange(gsMatrix<T> &coefs0);

// computes the maximum eigenvalue of a symmetric 2x2 matrix with real coefficients
// input: matrix Nmat
// output: T number that is the maximum eigenvalue
    T lambdamax(gsMatrix<T>& Nmat);

// computes the inner product of two vectors vec1 and vec2
// input: two vectors vec1, vec2
// output: T number that is the inner product
    T innerprod(const Vector& vec1,const Vector& vec2);

// transforms the vector of the gradient of the bound into a matrix with 2 columns (since we are in 2D)
// input: vector with gradient of bound
// output: matrix with gradient of bound
    gsMatrix<T> transvecmatgradbound(Vector& gradboundvec);

// transforms the matrix of the gradient of the bound into a vector
// input: matrix with gradient of bound
// output: vector with gradient of bound
    Vector transmatvecgradbound(gsMatrix<T>& gradboundmat);

private:

}; // class gsOptParamConditionNumber


//////////////////////////////////////////////////
//////////////////////////////////////////////////

template<class T>
void gsOptParamConditionNumber<T>::run()
{
    const gsBSpline<T> * bsp = dynamic_cast<const gsBSpline<T>*>(&m_input);
    if ( bsp )
    {
        run1D(*bsp);
        return;
    }
    
    const gsTensorBSpline<2,T> * bsp2 = dynamic_cast<const gsTensorBSpline<2,T>*>(&m_input);
    if ( bsp2 )
    {
        runnD<2>(*bsp2);
        return;
    }
    
    const gsTensorBSpline<3,T> * bsp3 = dynamic_cast<const gsTensorBSpline<3,T>*>(&m_input);
    if ( bsp3 )
    {
        runnD<3>(*bsp3);
        return;
    }

    GISMO_ERROR("Error, geometry type cannot be optimized.");
}


template<class T>
void gsOptParamConditionNumber<T>::run1D(gsBSpline<T> const & bsp)
{
    const gsKnotVector<T> &  knots  = bsp.knots();
    //int numcpt = bsp.basis().size();

    //Vector gradientfct, gradientfctnew;
    //gradientfct = computegradbound2(knots,bcoefs,10,2);
    //gradientfctnew = computegradbound2(knots,bcoefs,10,2);

    gsMatrix<T> current  = bsp.coefs();

    //int aa;
    //aa=coefsinrange(bcoefs);
    
    for (int i=0; i!=m_max_iter; ++i)
    {
        //gradientfct = computegradbound2(knots, current, 10, 2);
        
        loopstep1D(knots, current, 10, 2);

        //gradientfctnew = computegradbound2(knots, current, 10, 2);
        //aa=coefsinrange(bcoefs);
    }

    m_bound = computebound1D2(knots,current, 10, 2);
    
    if ( m_output )
    {
        m_output->coefs() = current;
    }
    else
    {
        m_output = new gsBSpline<T>(knots,current);
    }
}


template<class T>
T gsOptParamConditionNumber<T>::discretelpnorm(const Vector& vec,T p)
    {
    T lpvalue=0;
    for (size_t i=0;i!=vec.size();++i)
    {
        lpvalue+=(math::pow(vec[i],p));
    }
    return (math::pow(lpvalue,1/p));
    }


template<class T>
T gsOptParamConditionNumber<T>::computebound1D2(const gsKnotVector<T>& KV,gsMatrix<T> &coefs0, T p1, T p2)
    {
    gsBSpline<T> geom_0(KV, coefs0);

    // Evaluation points
    T x1=-1./3.*sqrt(5.+2.*sqrt(10./7.));
    T x2=-1./3.*sqrt(5.-2.*sqrt(10./7.));
    T x3=0;
    T x4=1./3.*sqrt(5.-2.*sqrt(10./7.));
    T x5=1./3.*sqrt(5.+2.*sqrt(10./7.));

    // weights
    Vector weights;
    T w1=(322.-13.*sqrt(70.))/900.;
    T w2=(322.+13.*sqrt(70.))/900.;
    T w3=128./225.;
    T w4=(322.+13.*sqrt(70.))/900.;
    T w5=(322.-13.*sqrt(70.))/900.;
    weights.push_back(w1);
    weights.push_back(w2);
    weights.push_back(w3);
    weights.push_back(w4);
    weights.push_back(w5);

    Vector integrals;
    for (size_t i=2;i!=KV.size()-3;++i)
        {
        T a=KV[i];
        T b=KV[i+1];
        gsMatrix<T> paramvals(1,5);
        paramvals << (b-a)/2*x1+(a+b)/2, (b-a)/2*x2+(a+b)/2, (b-a)/2*x3+(a+b)/2, (b-a)/2*x4+(a+b)/2, (b-a)/2*x5+(a+b)/2;
        gsMatrix<T> results_inside(1,5);
        geom_0.jacobian_into(paramvals,results_inside);
        T sumj=0;
        for (int j=0;j!=5;++j)
            {
            sumj=sumj+(weights[j]*(math::pow(1/((b-a)*(b-a)*results_inside(0,j)*results_inside(0,j)),p1)));
            }
        sumj=sumj*(b-a)/2;
        integrals.push_back((math::pow(sumj,1/p1)));
        }
     return discretelpnorm(integrals,p2);
    }

template<class T>
T gsOptParamConditionNumber<T>::innerprod(const Vector& vec1,const Vector& vec2)
    {
    T result=0.;
    for(size_t i=0;i!=vec1.size();++i)
        result =result+vec1[i]*vec2[i];
    return result;
    }

template<class T>
typename gsOptParamConditionNumber<T>::Vector
gsOptParamConditionNumber<T>::computegradbound2(const gsKnotVector<T>& KV,gsMatrix<T> &coefs0,T p1,T p2)
    {
    int numcpt=coefs0.size();
    T delta=0.0001;
    std::vector<gsMatrix<T> > coefsdeltaplus;
    std::vector<gsMatrix<T> > coefsdeltaminus;
    for (int i=0; i!=numcpt;++i)
    {
    gsMatrix<T> coefsdeltainsideplus(numcpt,1);
    gsMatrix<T> coefsdeltainsideminus(numcpt,1);
    coefsdeltainsideplus=coefs0;
    coefsdeltainsideminus=coefs0;
    coefsdeltainsideplus(i,0)=coefsdeltainsideplus(i,0)+delta;
    coefsdeltainsideminus(i,0)=coefsdeltainsideminus(i,0)-delta;
    coefsdeltaplus.push_back(coefsdeltainsideplus);
    coefsdeltaminus.push_back(coefsdeltainsideminus);
    }

    Vector gradbound;
    gradbound.push_back(0);
    for (int i=1; i!=numcpt-1;++i)
    {
    gradbound.push_back((computebound1D2(KV,coefsdeltaplus[i],p1,p2)-computebound1D2(KV,coefsdeltaminus[i],p1,p2))/(2*delta));
    }
    gradbound.push_back(0);
    return gradbound;
    }


template<class T>
void gsOptParamConditionNumber<T>::loopstep1D(const gsKnotVector<T>& KV, gsMatrix<T> &coefs0, T p1, T p2)
   {
    int numcpt=coefs0.size();
    Vector gradbound=computegradbound2(KV,coefs0,p1,p2);

    gsMatrix<T> coefsgradbound(numcpt,1);
    for (int i=0;i!=numcpt;++i)
    {
    coefsgradbound(i,0)=gradbound[i];
    }

    gsMatrix<T> coefsges(numcpt,1);

    T sigma=0.5;
    int j=0;
    do
    {
    j=j+1;
    coefsges=coefs0-(math::pow(sigma,j))*coefsgradbound;
    }
    while (computebound1D2(KV,coefsges,p1,p2) > computebound1D2(KV,coefs0,p1,p2)-0.5*(math::pow(sigma,j))*innerprod(gradbound,gradbound));

    coefs0 -= math::pow(sigma,j)*coefsgradbound;
   }

template<class T>
int gsOptParamConditionNumber<T>::coefsinrange(gsMatrix<T> &coefs0)
{
    int a=0;
    for (int i=0;i!=coefs0.size();++i)
    {
        if ((coefs0(i)>=0) && (coefs0(i)<=1))
        {
        }
        else
        {a=a+1;}
    }
    return a;
}

template<class T>
T gsOptParamConditionNumber<T>::lambdamax(gsMatrix<T>& Nmat)
{
    return (Nmat(0,0)+Nmat(1,1)+sqrt(4*Nmat(1,0)*Nmat(1,0)+(Nmat(0,0)-Nmat(1,1))*(Nmat(0,0)-Nmat(1,1))))/2;
}


template<class T>
T gsOptParamConditionNumber<T>::computebound2d(const gsKnotVector<T>& KV1,const gsKnotVector<T>& KV2,gsMatrix<T> &coefs2D)
{
    gsTensorBSpline<2> *T_geom2D = new gsTensorBSpline<2>( KV1, KV2, coefs2D);
    // Evaluation points
    Vector xeval;
    T x1=-1./3.*sqrt(5.+2.*sqrt(10./7.));
    T x2=-1./3.*sqrt(5.-2.*sqrt(10./7.));
    T x3=0;
    T x4=1./3.*sqrt(5.-2.*sqrt(10./7.));
    T x5=1./3.*sqrt(5.+2.*sqrt(10./7.));
    xeval.push_back(x1);
    xeval.push_back(x2);
    xeval.push_back(x3);
    xeval.push_back(x4);
    xeval.push_back(x5);

    // weights
    Vector weights;
    T w1=(322.-13.*sqrt(70.))/900.;
    T w2=(322.+13.*sqrt(70.))/900.;
    T w3=128./225.;
    T w4=(322.+13.*sqrt(70.))/900.;
    T w5=(322.-13.*sqrt(70.))/900.;
    weights.push_back(w1);
    weights.push_back(w2);
    weights.push_back(w3);
    weights.push_back(w4);
    weights.push_back(w5);

    T p1=20;
    T absdet1;
    gsMatrix<T> intermedresult1;
    gsMatrix<T> matinvtransp1;
    gsMatrix<T> Nmat1;
    T lambdaN1;

    Vector vecges;
    Vector vecgeszaehler;
    gsMatrix<T> paramvals2d2(2,1);
    gsMatrix<T> results_inside2d2;
    T boundvalins;
    T boundvalouts;
    T boundvalinszaehler;
    T boundvaloutszaehler;
    Vector boundvalsxi2;
    Vector boundvalsxi2zaehler;
    T a1,a2,b1,b2;

    for(int ii=2;ii!=KV2.size()-3;++ii) // loop for y-direction
    {
    a2=KV2[ii];
    b2=KV2[ii+1];
    for(int kk=2;kk!=KV1.size()-3;++kk) // loop for x-direction
        {
        a1=KV1[kk];
        b1=KV1[kk+1];
        boundvalouts=0;
        boundvaloutszaehler=0;
        for (int i=0;i!=5;++i)
            {
              boundvalins=0;
              boundvalinszaehler=0;
            for (int j=0;j!=5;++j)
              {
                paramvals2d2 << (b1-a1)/2*xeval[i]+(a1+b1)/2,(b2-a2)/2*xeval[j]+(a2+b2)/2;
                T_geom2D->jacobian_into(paramvals2d2,results_inside2d2);
                absdet1=math::abs(results_inside2d2.determinant());
                intermedresult1=absdet1*results_inside2d2.inverse();
                matinvtransp1=(results_inside2d2.inverse()).transpose();
                Nmat1=intermedresult1*matinvtransp1;
                lambdaN1=lambdamax(Nmat1);
                boundvalins=boundvalins+(b2-a2)/2*weights[j]*math::pow(1/(math::abs(results_inside2d2.determinant())*(b1-a1)*(b2-a2)),p1);
                boundvalinszaehler=boundvalinszaehler+(b2-a2)/2*weights[j]*math::pow(lambdaN1*((b1-a1)*(b1-a1)+(b2-a2)*(b2-a2))/((b1-a1)*(b2-a2)),p1);
               }
              boundvalsxi2.push_back(boundvalins);
              boundvalsxi2zaehler.push_back(boundvalinszaehler);
              boundvalouts=boundvalouts+weights[i]*(b1-a1)/2*boundvalins;
              boundvaloutszaehler=boundvaloutszaehler+weights[i]*(b1-a1)/2*boundvalinszaehler;
             }
        vecges.push_back(math::pow(boundvalouts,1/p1));
        vecgeszaehler.push_back(math::pow(boundvaloutszaehler,1/p1));
        }
     }
  return discretelpnorm(vecgeszaehler,20)*discretelpnorm(vecges,20);
  }


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// computes the gradient of the 2D bound function
// input: knot vectors KV1 and KV2, matirx with control points coefs2D
// output: vector with gradient (if gradient needed as matrix, use function transvecmatgradbound
template<class T>
typename gsOptParamConditionNumber<T>::Vector gsOptParamConditionNumber<T>::computegradbound2d(const gsKnotVector<T>& KV1,const gsKnotVector<T>& KV2,gsMatrix<T> &coefs2D)
{
    T delta=0.0001;
    std::vector<gsMatrix<T> > coefsdeltaplus;
    std::vector<gsMatrix<T> > coefsdeltaminus;
    for (int i=0; i!=coefs2D.rows();++i)
    {
        for(int j=0;j!=coefs2D.cols();++j)
        {
            gsMatrix<T> coefsdeltainsideplus;
            gsMatrix<T> coefsdeltainsideminus;
            coefsdeltainsideplus=coefs2D;
            coefsdeltainsideminus=coefs2D;
            coefsdeltainsideplus(i,j)=coefsdeltainsideplus(i,j)+delta;
            coefsdeltainsideminus(i,j)=coefsdeltainsideminus(i,j)-delta;
            coefsdeltaplus.push_back(coefsdeltainsideplus);
            coefsdeltaminus.push_back(coefsdeltainsideminus);
            
        }
    }
    Vector gradbound;
    T toDivideByDelta;
    for (size_t i=0; i!=coefsdeltaplus.size();++i)
    {
        toDivideByDelta = (computebound2d(KV1,KV2,coefsdeltaplus[i])-computebound2d(KV1,KV2,coefsdeltaminus[i]));
        gradbound.push_back( toDivideByDelta /(2*delta));
    }
    return gradbound;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// transforms the vector of the gradient of the bound into a matrix with 2 columns (since we are in 2D)
// input: vector with gradient of bound
// output: matrix with gradient of bound
template<class T>
gsMatrix<T> gsOptParamConditionNumber<T>::transvecmatgradbound(Vector& gradboundvec)
   {
    int rows=gradboundvec.size()/2;
    gsMatrix<T> gradboundmat(rows,2);
    for (int i=0;i!=rows;++i)
    {
    gradboundmat(i,0)=gradboundvec[2*i];
    gradboundmat(i,1)=gradboundvec[2*i+1];
    }
    return gradboundmat;
   }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// transforms the matrix of the gradient of the bound into a vector
// input: matrix with gradient of bound
// output: vector with gradient of bound
template<class T>
typename gsOptParamConditionNumber<T>::Vector
gsOptParamConditionNumber<T>::transmatvecgradbound(gsMatrix<T>& gradboundmat)
   {
    gradboundmat.transpose();
     Vector gradboundvec;
     for(int i=0;i!=gradboundmat.size();++i)
     gradboundvec.push_back(gradboundmat(i));
     gradboundmat.transpose();

     return gradboundvec;
    }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// computes one step of the gradient method to compute the bound
// input: knot vectors KV1 and KV2, control points coefs2D, vector coefschange to specify which control points can be changed
// output: new control points
template<class T> gsMatrix<T> 
gsOptParamConditionNumber<T>::loopstep2d(const gsKnotVector<T>& KV1, const gsKnotVector<T>& KV2, 
                                     gsMatrix<T> &coefs2D, std::vector<int> coefschange)
   {
    Vector gradboundold=computegradbound2d(KV1,KV2,coefs2D);
    Vector gradbound;

    for(size_t i=0;i!=gradboundold.size();++i)
      {
        if(coefschange[i]==0)
        {gradbound.push_back(0);}
        else
        {
            gradbound.push_back(gradboundold[i]);
        }
      }

    gsMatrix<T> coefsgradbound;
    coefsgradbound=transvecmatgradbound(gradbound);
    gsMatrix<T> coefsges;
    Vector coefsgesvecold;
    Vector coefsgesvec;

    T sigma=0.5;
    int j=0;
    do
    {

    j=j+1;
    coefsges=coefs2D-(math::pow(sigma,j))*coefsgradbound;
    coefsgesvecold=transmatvecgradbound(coefsges);
      for(size_t i=0;i!=coefschange.size();++i)
      {
        if(coefschange[i]==0)
        {coefsgesvec.push_back(0);}
        else
        {coefsgesvec.push_back(coefsgesvecold[i]);}
      }
    }

    while (computebound2d(KV1,KV2,coefsges) > computebound2d(KV1,KV2,coefs2D)-0.5*(math::pow(sigma,j))*innerprod(gradbound,gradbound) );

    gsMatrix<T> coefsgesnew;
    coefsgesnew=coefs2D-math::pow(sigma,j)*coefsgradbound;
    return coefsgesnew;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template<class T>
std::ostream & operator <<( std::ostream &os, gsOptParamConditionNumber<T> const & cf)
{
    os <<"Parametrization optimizer.. \n";
    return os;
}

}; // namespace gismo
