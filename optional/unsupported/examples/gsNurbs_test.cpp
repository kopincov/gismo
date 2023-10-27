
#include <gismo.h>
#include <gismo_dev.h>

#include <iostream>

using namespace gismo;

template<class gType> 
bool checkDynamicGeometryType(const gsGeometry<> *geo1, const gType *geo2, const std::string basisname, const std::string geoname)
{
    const bool ok = ( NULL != dynamic_cast<const gType *>(geo1) );
    gsInfo << "Type "<<(ok ? "" : " NOT ")<<"OK, "<<basisname<<"/"<<geoname<<".\n";
    delete geo1;
    return ok;
}


bool checkExactInterpolation(const gsBasis<> & basis, const gsFunction<> & f)
{
    gsGeometry<>::uPtr intp = basis.interpolateAtAnchors( f.eval(basis.anchors()) );

    gsMatrix<> paramRange = basis.support();
    const real_t maxDist = computeMaximumDistance<real_t>(f, *intp, paramRange.col(0), paramRange.col(1));
    gsInfo << "error: " << maxDist;
    
    if (maxDist > math::limits::epsilon()*1000)
    {
        gsInfo << " FAIL" << "\n";
        return false;
    }
    else
    {
        gsInfo << " OK" << "\n";
        return true;
    }
}


int main()
{
  bool passed = true;

  // *************************************************
  // Knot vector
  // *************************************************
  gsInfo<< "\n------------ Knot vectors tests ------------" << "\n";

  std::vector<real_t>::iterator it;

  gsKnotVector<> KV (-1, 0, 3,3, 1 ) ;
  gsInfo<< "A knot vector:"<< KV << "\n";
  KV.addConstant(1);
  gsInfo<< "Shifted by 1:"<< KV << "\n";

  std::vector<real_t> breaks = KV.unique() ;
  gsInfo << "Breaks:";
  for (it=breaks.begin(); it!=breaks.end(); ++it)
    gsInfo << " " << *it;
  gsInfo << "\n";

  gsKnotVector<> KV2( breaks, 2, 0 ) ;
  gsInfo<< "Knot vector with double knots: "<< KV2 << "\n";


  // *************************************************
  // Knot insertion
  // *************************************************
  {
  gsInfo<< "\n------------ Knot insertion test ------------" << "\n";

  gsKnotVector<> KV3(0, 1, 4, 4);
  gsBSplineBasis<> bsp3(KV3);
  gsInfo << "KV3 = " << KV3 << "\n";
  const int n = bsp3.size();    // == KV3.size() - KV3.degree() - 1

  gsMatrix<> eye = gsMatrix<>::Identity(n, n);
  gsMatrix<> col_before = eye.col(3);
  gsBSpline<> curve_before(KV3, col_before);

  gsInfo << "Coefficients before:" << "\n" << eye << "\n";
  gsInfo << "Inserting knot 0.5 ..." << "\n";
  gsBoehm(KV3, eye, (real_t)(0.5));
  gsInfo << "Coefficients after:" << "\n" << eye << "\n";

  gsMatrix<> col_after = eye.col(3);
  gsBSpline<> curve_after(KV3, give(col_after));

  gsMatrix<> paramRange = curve_before.parameterRange();
  const real_t maxDist = computeMaximumDistance<real_t>(curve_before, curve_after, paramRange.col(0), paramRange.col(1));

  gsInfo << "Estimated maximum distance = " << maxDist << "\n";
  }

  // *************************************************
  // B-Spline basis & B-Splines
  // *************************************************
  gsInfo<< "\n------------ B-Splines tests ------------" << "\n";

  // some evaluation points
  gsMatrix<real_t,1> u(1,5) ;

  //u << 0.05, 0.12 , 0.75  ;
  u << 0, 0.4, 0.5 , 0.98, 1  ;
  
  gsBSplineBasis<>::Ptr bsp = gsBSplineBasis<>::make(KV);
  // Sharing the same basis
  gsBasis<>::Ptr bsp2 = bsp;
  gsInfo<< "A shared "<< *bsp << "\n";

  gsInfo<< "The index of first active basis function on "<< u.at(0)<<": "
      << bsp->firstActive( u.at(0)) << "\n";

  gsMatrix<> grev;
  bsp->anchors_into(grev);
  gsInfo<< "The Greville points: "<< "\n";
  gsInfo<<  grev << "\n";

  gsMatrix<>  tmp = grev.transpose() ;
  gsBSpline<> B ( *bsp, tmp ) ;

  gsMatrix<> w( B.basis().size(), 1) ; 
  w.setOnes();
  gsNurbs<> rB (KV, give(w), give(tmp)) ;

  gsInfo<< "The "<< B << "\n";

  gsInfo<< "Evaluate the B-Spline at "<<u <<" : ";
  gsInfo<<  B.eval( u )  <<"\n";
  gsInfo<< "Evaluate the NURBS at "<<u <<" : ";
  gsInfo<<  rB.eval( u )  <<"\n";

  gsInfo<< "Evaluate non-zero basis functions on "<< u <<": \n";
  gsInfo<<  bsp->eval( u )   <<"\n";
  int k = bsp->firstActive(u(0,1));
  for ( int i = 0; i<= bsp->degree(); ++i) 
      gsInfo<< "Eval basis function "<< k +i <<": "<< bsp->evalSingle(k + i,u) << "\n";
   
  gsInfo<< "Is the point "<< u.col(3).transpose() <<" on curve: "<< B.contains(u.col(3) ) << "\n";

  gsInfo<< "Evaluate non-zero NURBS basis functions on "<< u <<": \n";
  gsInfo<<  rB.basis().eval( u )   <<"\n";


  gsInfo<< "B-Spline derivative : \n";
  gsInfo<<  bsp->deriv( u )   <<"\n";
  gsInfo<< "Nurbs derivative : \n";
  gsInfo<<  rB.basis().deriv( u )   <<"\n";
  gsInfo<< "B-Spline Jacobian : \n";
  gsInfo<<  B.jacobian( u.col(1) )  <<"\n";
  gsInfo<< "NURBS Jacobian : \n";
  gsInfo<<  rB.jacobian( u.col(1) )  <<"\n";

  gsInfo<< "Degree elevation: TO DO \n";
  //bsp->elevate();

  //gsInfo<<  * bsp->evalAllDers( u, 3 )  <<"\n";
    
  gsInfo<< "Sample uniformly distributed (in parameter domain) curve points: " << B.sample(3) <<"\n";

  // *************************************************
  // Interpolation
  // *************************************************
  {
  gsInfo<< "\n------------ Interpolation test ------------" << "\n";

  gsKnotVector<> KV3(0, 1, 4, 4);
  gsKnotVector<> KV_ref(0, 1, 9, 4);

  gsBSplineBasis<> bsp3(KV3);
  gsInfo << "KV3 = " << KV3 << "\n";
  gsInfo << "KV_ref = " << KV_ref << "\n";

  const int n = KV3.size() - KV3.degree() - 1;
  gsMatrix<> coeffs3( gsMatrix<>::Base::Zero(n, 1) );
  coeffs3( coeffs3.size() / 2 ) = 1.0;
  gsBSpline<> curve_before(KV3, coeffs3);
  gsBSplineBasis<> bsp_ref(KV_ref);
  gsMatrix<> an_pts = bsp_ref.anchors();
  gsGeometry<>::uPtr curve_after = bsp_ref.interpolateData(curve_before.eval(an_pts),an_pts);

  gsInfo << "Coefficients before:" << "\n" << curve_before.coefs().transpose() << "\n";
  gsInfo << "Interpolated coefficients:" << "\n" << curve_after->coefs().transpose() << "\n";

  gsMatrix<> paramRange = curve_before.parameterRange();
  const real_t maxDist = computeMaximumDistance<real_t>(curve_before, *curve_after, paramRange.col(0), paramRange.col(1));
  gsInfo << "Estimated maximum distance = " << maxDist << "\n";
  }

  // *************************************************
  // Tensor Product basis
  // *************************************************
  gsInfo<< "\n------------ Tensor B-Splines tests ------------" << "\n";

  // gsTensorBasis takes ownership of component bases! Therefore we have to clone them.
  gsTensorBSplineBasis<3>::Ptr tbasis =
    gsTensorBSplineBasis<3>::Ptr (new gsTensorBSplineBasis<3>(bsp->clone().release(), bsp->clone().release(), bsp->clone().release()));

  gsInfo<< "The "<< * tbasis << "\n";

  gsInfo<< "Size 0: "<< tbasis->size(0) << "\n";
  gsInfo<< "Size 1: "<< tbasis->size(1) << "\n";


  //gsInfo<< "The Greville points: \n"<< "\n";
  tbasis->anchors_into(grev) ;
  //gsInfo<< grev->transpose() << "\n";

  gsMatrix<> a(3,3);
  a<< .4, .1, 1, 
      .4, .7, 0,
      .4, .2, 0;

  gsInfo<<"Points: \n"<< a <<"\n";
  gsInfo<< "Evaluate Tensor-product basis (each column is the non-zero basis functions. Numbering runs first over first direction, then second, then ...): \n";
  gsMatrix<> e = tbasis->eval( a );
//  gsInfo<< *e  << "\n";
  gsInfo<< "Check the sum of every column: "<< e.colwise().sum() << "\n";

  gsInfo<< "Evaluate Tensor-product (partial) derivatives of non-zero basis functions (each column containts the gradients the non-zero basis functions. Numbering runs first over first direction, then second, then ...): \n";
//  gsInfo<< tbasis->deriv( a ) << "\n";

//  gsInfo<< "Index of active basis functions at the points : \n"<< * tbasis->active(a)  <<"\n";

  gsInfo << "Testing interpolation for 2D tensor product B-splines... ";
  passed = checkExactInterpolation(gsTensorBSplineBasis<2,real_t>(KV, KV), gsFunctionExpr<>("x^2 + x*y^2",2));

  gsInfo << "Testing interpolation for 3D tensor product B-splines... ";
  passed = checkExactInterpolation(*tbasis, gsFunctionExpr<>("-x + y^2 + x*z + z^2",3)) && passed;

  tmp = grev.transpose();
  gsTensorBSpline<3> tbsp( *tbasis, tmp ) ;
  gsInfo<< "The "<< tbsp  << "\n";
  gsInfo<< "Evaluate Tensor-product B-Spline" << ": \n";
  gsInfo<< tbsp.eval( a ) << "\n";
  gsInfo<< "Boundary South "<< *tbsp.boundary(boundary::south).get()  << "\n";

  // *************************************************
  // Bernstein basis & Bezier curves
  // *************************************************
  gsInfo<< "\n------------ Bernstein bases ------------" << "\n";

  gsBernsteinBasis<> bb(0.0,1.0,3,3) ;

  gsInfo<< "The "<< bb << "\n";
  gsInfo<< "Evaluate Bernstein polynomials on "<< u <<": \n";
  gsInfo<<  bb.eval( u )   <<"\n";
  gsInfo<< "Evaluate Bernstein derivatives on "<< u <<": \n";
  gsInfo<<  bb.deriv( u )   <<"\n";

  bb.anchors_into(grev) ;
  gsInfo<< "The Greville points: "<< "\n";
  gsInfo<<  grev << "\n";

  grev.transposeInPlace();
  gsBezier<> bez (bb, grev ) ;

  gsInfo<< "Evaluate identity Bezier curve at "<< u <<": ";
  gsInfo<<  bez.eval ( u )  <<"\n";
  gsInfo<< "Evaluate derivative at "<< u <<": ";
  gsInfo<<  bez.deriv( u )  <<"\n";

  bez.uniformRefine();
  gsInfo<< "Refined "<< bez << "\n";
  gsInfo<< "Evaluate refined identity Bezier curve at "<< u <<": ";
  gsInfo<<  bez.eval ( u )  <<"\n";

  // *************************************************
  // Tensor Bezier
  // *************************************************
  gsInfo<< "\n------------ Tensor Bezier tests ------------" << "\n";

  gsTensorBernsteinBasis<3> tbern(bb.clone().release(), bb.clone().release(), bb.clone().release());
  gsInfo<< "The "<< tbern << "\n";
  tbern.anchors_into(grev) ;
  gsInfo<< "Evaluate 3D Bernstein polynomials on "<< a <<": \n";
  e = tbern.eval( a ) ;
  gsInfo<< "Check the sum of every column: "<< e.colwise().sum() << "\n";


  gsInfo << "Interpolating: "<< tbsp << "\n";
  gsGeometry<>::uPtr tbsp_intpl = tbsp.basis().interpolateAtAnchors( tbsp.eval(tbsp.basis().anchors()) ); //todo Is it necessary?
  //gsInfo << "Coefficients before:\n" << tbsp2.coefs().transpose() << "\n";
  //gsInfo << "Interpolated coefficients:\n" << tbsp_intpl->coefs().transpose() << "\n";

  // *************************************************
  // NURBS
  // *************************************************
  gsInfo<< "\n------------ Rational / NURBS bases ------------" << "\n";
  
  gsNurbsBasis<>  nrbbasis( KV ) ;
  //nrbbasis->weight(2) = 20;
  gsInfo<< "1. "<<  nrbbasis  <<"\n";
  gsInfo<<"Evaluate NURBS basis: \n" <<  nrbbasis.eval(u)  <<"\n";

  // Test NURBS
  nrbbasis.anchors_into(grev) ;
  grev.transposeInPlace();
  gsNurbs<> N ( nrbbasis, grev ) ;
  gsInfo<<  "The "<< N  <<"\n";
  gsInfo<< "Evaluate the NURBS at "<<u <<" : ";
  gsInfo<< N.eval( u ) <<"\n";


  gsTensorNurbsBasis<3> ratbasis( *tbasis ) ;

  gsInfo<< "4. "<< ratbasis  <<"\n";

  gsMatrix<> tnrbC = tbsp.coefs();
  gsTensorNurbs<3> tnrbt( KV, KV, KV, give(tnrbC) );
  gsInfo<< "5. "<<  tnrbt  <<"\n";
  tnrbt.eval(a) ;

  gsInfo << "Testing interpolation for tensor product NURBS... ";
  passed = checkExactInterpolation(*tbasis, gsFunctionExpr<>("-x + y^2 + x*z + z^2",3)) && passed;

// TO DO   tnrbt->basis().uniformRefine();

  {
    gsNurbs<>::uPtr circle = gsNurbsCreator<>::NurbsCircle();
    gsNurbs<>::uPtr circle2 = circle->clone();
    circle2->uniformRefine();
    gsInfo << "NURBS circle: " << circle->coefsSize() << " control points; refined: " << circle2->coefsSize() << " control points. ";

    gsMatrix<> paramRange = circle->parameterRange();
    const real_t dist = computeMaximumDistance<real_t>(*circle, *circle2, paramRange.col(0), paramRange.col(1));

    gsInfo << "Maximum distance: " << dist << "\n";
    if (dist >math::limits::epsilon()*20) { passed = false; gsInfo << "FAIL" << "\n"; }

  }

  gsInfo << "\n------------ Types of geometry ------------" << "\n";

  if (!checkDynamicGeometryType(bez.basis().makeGeometry( bez.coefs()).release(), &bez, "gsBernsteinBasis", "gsBezier"))
      passed = false;
  if (!checkDynamicGeometryType(bsp->makeGeometry( B.coefs()).release(), &B, "gsBSplineBasis", "gsBSpline"))
      passed = false;
  if (!checkDynamicGeometryType(nrbbasis.makeGeometry(grev).release(), &N, "gsNurbsBasis", "gsNurbs"))
      passed = false;
  if (!checkDynamicGeometryType(tbasis->makeGeometry( tbsp.coefs()).release(), &tbsp, "gsTensorBSplineBasis", "gsTensorBSpline"))
      passed = false;
  if (!checkDynamicGeometryType(ratbasis.makeGeometry( tnrbt.coefs() ).release(), &tnrbt, "gsTensorNurbsBasis", "gsTensorNurbs"))
      passed = false;

  gsInfo << "\n";

  return passed ? 0 : 1;
}
