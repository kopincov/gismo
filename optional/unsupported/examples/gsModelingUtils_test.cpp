//#define NDEBUG
#define debug 1
#include <iostream>

#include <gismo.h>
#include <gsModeling/gsModelingUtils.hpp>


using namespace gismo;


/// Multiply each element of a colume to each row of a matrix
template <class T>
gsMatrix<T> rowProduct(gsMatrix<T> const & col, gsMatrix<T> const & mat)
{
    GISMO_ASSERT( col.cols()==1 && col.rows()== mat.rows(), "Bad dimensions" );

    return col.asDiagonal() *  mat ;
}

int main()
{
  gsMatrix<> A(2, 2);
  gsMatrix<> C(1, 2);
  gsVector<> d(1);
  
  A << 2, 0, 0, 1;
  C << 1, 1;
  d << 1;
  
  gsInfo << "\nA = \n" << A << "\n\nC = \n" << C << "\n\nd = \n" << d << "\n";
  
  gsVector<> soln = criticalPointOfQuadratic(A, C, d);
  gsInfo << "\nThe critical point of X^real_t A X subject to C X = d is\n" << soln << "\n";
  
  gsMatrix<> bv(2,1); bv<<1,1;
  gsMatrix<> dv(1,1); dv<<1;
  gsVector<> solopt = optQuadratic<>(A, bv, C, dv);
  gsInfo << "\nThe critical point of X^real_t A X + bX subject to C X = d is\n" << solopt << "\n";
  solopt = optQuadratic<real_t>(A, bv,2,A,bv,1, C, dv);
  gsInfo << "\nThe critical point of X^real_t A X + bX subject to C X = d is\n" << solopt << "\n";
  
  
  
  //.......................  
  
  gsVector3d<> vec1; vec1 << -2,-2,0;
  gsVector3d<> vec2; vec2 << 2,2,0;
  gsInfo<<"\n"<< " Angle: " << conditionedAngle<>(vec1,vec2);
  
  gsVector3d<> normal; normal << -1,0,0;
  gsInfo<<"\n"<< "Conditioned Angle: " << conditionedAngle<>(vec1,vec2,normal);
  
  
  gsInfo<<"\n"<<"==================== removeCol ===================="<<"\n";
  gsMatrix<> Mrc(2,6); Mrc << 0,1,2,3,4,5,0,1,2,3,4,5;
  //removeCol<>(Mrc,3);
  //gsInfo<<"\n"<<"Mrc after removement of col 1: "<<"\n"<< Mrc <<"\n";
  //removeCol<>(Mrc,1, 3);gsInfo<<"\n"<<"Mrc after removement of one end: "<<"\n"<< Mrc <<"\n";
  removeCol<>(Mrc,2, 3);gsInfo<<"\n"<<"Mrc after removement of two ends: "<<"\n"<< Mrc <<"\n";
  
  gsInfo<<"\n"<< " ------------------- test addConstraints --------------------";    
  gsMatrix<> CC,dd;
  gsMatrix<> C1(1,2); C1 << 0,0;
  gsMatrix<> d1(1,1); d1 << 0;
  gsMatrix<> C2(2,2); C2 << 1,1,2,2;
  gsMatrix<> d2(2,1); d2 << 1,2;
  addConstraints<>(C1,d1,C2,d2,CC,dd);
  gsInfo<<"\n"<<"Assembled constraints: "<<"\n"<<CC<<"\n"<<dd;
  gsInfo<<"\n"<<" C2 row-multiplies with d2: "<<"\n"<< rowProduct<>(d2,C2);
  
  
  
  gsInfo<<"\n"<< " ------------------- test kroneckerProduct --------------------";    
  gsMatrix<> m1(1,2); m1<< 1,2;
  gsMatrix<> m2(2,2); m2<< 1,1,1,1;
  gsInfo<<"\n"<<"Kronecker product: "<<"\n"<< m1.kron(m2) <<"\n";
  
  gsInfo<<"\n"<< " ------------------- test knotVector == + innerProduct --------------------";
  int p1 = 1, p2 = 1;
  gsKnotVector<> knotv1(0,1,1,p1+1);
  gsKnotVector<> knotv2(0,1,1,p2+1);
  gsInfo<<"\n"<<"The two knot vectors are the same, true or false: \n"<< (knotv1==knotv2) <<"\n";
  
  gsBSplineBasis<> bs1(knotv1), bs2(knotv2);
  gsMatrix<> * K = innerProduct(bs1, bs2);
  gsInfo<<"\n"<<"innerProduct: \n"<< *K <<"\n";
  gsInfo<<"\n"<<"sum of innerProduct (should be 1): "<< K->sum() <<"\n";
  
  delete K;
  K= innerProduct1(bs1, bs2);
  gsInfo<<"\n"<<"innerProduct1: \n"<< *K <<"\n";
  delete K;
  K= innerProduct2(bs1, bs2);
  gsInfo<<"\n"<<"innerProduct2: \n"<< *K <<"\n";
  delete K;
   
  gsInfo<< "\n" << "==================== gsCurveLoop ====================" << "\n";
  // 	First make 4 spline curves which form the considered trimming loop
  gsMatrix<> tcp0(2,2); tcp0 << 0, 0, 1, 0;
  gsBSpline<> * tcurve0 = new gsBSpline<>( 0,1,0,1, give(tcp0) );
  gsMatrix<> tcp1(2,2); tcp1 << 1, 0, 1 , 1;
  gsBSpline<> * tcurve1 = new gsBSpline<>( 0,1,0,1, give(tcp1) );
  gsMatrix<> tcp2(2,2); tcp2 << 1, 1, 0 , 1;
  gsBSpline<> * tcurve2 = new gsBSpline<>( 0,1,0,1, give(tcp2) );
  gsMatrix<> tcp3(2,2); tcp3 << 0, 1, 0 , 0;
  gsBSpline<> * tcurve3 = new gsBSpline<>( 0,1,0,1, give(tcp3) );
  
  // 	Then construct a curve loop
  gsCurveLoop<> * tloop = new gsCurveLoop<>();
  tloop->insertCurve( tcurve0 );tloop->insertCurve( tcurve1 );
  tloop->insertCurve( tcurve2 );tloop->insertCurve( tcurve3 ); 
  gsInfo << "Sampled points: \n"<< tloop->sample(3) <<"\n";
  gsInfo << "Sampled points with one end: \n"<< tloop->sample(3,1) <<"\n";
  gsInfo << "Sampled points with no end: \n"<< tloop->sample(3,0) <<"\n";
  delete tloop;
  
  gsInfo<< "\n" << "==================== shared_ptr ====================" << "\n";
  memory::shared_ptr<int> pp1(new int);
  *pp1 = 1;
  gsInfo << *pp1;  
  gsInfo << "use_count: " << pp1.use_count() << '\n'; 
  
  
  gsInfo<< "\n" << "==================== Eigen ====================" << "\n";  
  gsMatrix<> a(3,1); a<<2,0,0;
  gsMatrix<> b(3,1); b<<0,2,0;
  gsMatrix<> c = (gsVector3d<>(a)).cross(gsVector3d<>(b));
  gsInfo<<"\n"<< "Norm: "<<"\n"<< a.norm();
  gsInfo<<"\n"<< "Normalized a: "<<"\n"<< a.normalized();
  gsInfo<<"\n"<< "Cross: "<<"\n"<< c;
  gsInfo<<"\n"<< "Size of the cross: "<< c.rows() <<" x "<<c.cols();
  
  
  gsInfo<<"\n";return 0;
}






