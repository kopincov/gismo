#include <string>
#include "gismo_unittest.h"
 
TEST(gsThbs_test_derivatives)
{
  std::string filename = "thbs_16.xml";
  gsFileData<>  data( filename );
  gsTHBSpline<2>::uPtr hbs = data.getFirst< gsTHBSpline<2> >();
  
//  gsTHBSplineBasis<2>  & HB = hbs->basis();
  gsMatrix<> para  = hbs->parameterRange();
  para(0,0)= 0.1;
  para(1,0)= 0.1;
  para(0,1)= 0.9;
  para(1,1)= 0.9;
  gsVector<> c0 = para.col(0);
  gsVector<> c1 = para.col(1);
  gsMatrix<> pts = uniformPointGrid(c0,c1, 8) ;
  
  // Compute first derivatives by a finite difference implementation and
  // by the gsTHBSplineBasis implementation and confront the results
  gsMatrix<> der1, hbs_der1;
  hbs->gsFunction<>::deriv_into(pts,der1); // finite differences
  hbs->deriv_into(pts,hbs_der1);           // implementation in gsTHBSplineBasis
  
  CHECK_MATRIX_CLOSE(der1, hbs_der1, (real_t)0.000001);

  gsMatrix<> der2, hbs_der2;
  hbs->gsFunction<>::deriv2_into(pts,der2);
  hbs->deriv2_into(pts,hbs_der2);
  // gsInfo<<"Second derivative by finite difference method "<<"\n"<<der2<<"\n";
  // gsInfo<<"Second derivative by THB spline derivative"  <<"\n"<<hbs_der2<<"\n";
  
  CHECK_MATRIX_CLOSE(der2, hbs_der2, (real_t)0.001);
}
