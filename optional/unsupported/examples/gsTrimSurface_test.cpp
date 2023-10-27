#define debug 1
#include <iostream>

#include <algorithm> // for copy
#include <iterator> // for ostream_iterator

#include <gismo.h>
#include <gismo_dev.h>




using namespace gismo;


int main(int argc, char *argv[])
{

if (1)
  {
  /// Master NURBS surface: degree 2, Bezier patch
  gsKnotVector<> Kv(0,2,1,4);
  gsMatrix<> C(25,3);
  C<< 0.124215, 0.0620467, -0.909323,
    0.990116, -0.112269, -0.0161917,
    1.99052, -0.213628, 1.12389,
    2.99999, -0.176554, -0.0225916,
    3.93311, -0.232422, -0.736062,
    -0.0212573, 1.14758, -0.383038,
    1, 1, 0,
    2, 1, 0,
    3, 1, 0,
    4.0773, 1.08316, -0.416151,
    -0.0554789, 2.0257, -0.138956, 
    1, 2, 0, 
    1.68928, 1.82642, 3.25086,
    3, 2, 0,
    3.94707, 2.11565, -0.357777,
    -0.112924, 2.99649, -0.0987988,
    1, 3, 0,
    2, 3, 0,
    3, 3, 0,
    4.0241, 3.11986, -0.280967,
    -0.209772, 4.06133, -0.303981,
    0.90333, 4.17869, -0.250397,
    1.98218, 4.13487, -0.119872,
    3.10997, 4.14081, -0.364305,
    4.09176, 4.18805, -0.43237;

  gsTensorBSpline<2,real_t>::Ptr master(new gsTensorBSpline<2, real_t>( Kv, Kv, give(C) ));
  //gsInfo<<"\n"<< "Control points of the base surface: \n"<< master->coefs()<<"\n";
  
  /// A trimming curve
  gsBSpline<> * trim_curve = gsNurbsCreator<>::BSplineFatCircle( 0.9 ).release();
  trim_curve->coefs().col(0).array() += 1; 
  trim_curve->coefs().col(1).array() += 1; 

  /// Planar domain by an outer boundary given by a curve
  gsPlanarDomain<> * domain = new gsPlanarDomain<>(trim_curve);
  // rotate the domain !!!

  gsTrimSurface<> ts(master, domain );
 
  // Print the trim surface
  gsInfo<< "I am a "<< ts <<"\n";

  // std::vector<real_t> x = gsBSplineRoot( *trim_curve , 0, 0.4);
  // gsAsMatrix<> xx (x) ;
  // gsInfo<< "Roots: ";;
  // gsInfo << xx << "\n";
  // gsInfo<< "Curve points:\n"<< trim_curve->eval(xx)  <<"\n";

  //gsMesh<> * pm = domain->toMesh(10);
  gsMesh<>::uPtr m = ts.toMesh(10);
  gsInfo << *m << "\n";
  gsWriteParaview( *m , "output" );
  //gsWriteParaview( *pm , "outputp" );
  
  // Call paraview on exit
  char cmd[100];
  strcpy(cmd,"paraview output.vtp\0");
  strcat(cmd," &");  
  
  //return system(cmd);  
  //return 0;
  };

  
  gsInfo<<"\n"<<"==================== Constructing a trimmed surface ===================="<<"\n";   
  //  First make 4 spline curves which form the considered trimming loop
  gsMatrix<> tcp0(2,2); tcp0 << 0, 0, 1, 0;
  gsBSpline<> * tcurve0 = new gsBSpline<>( 0,1,0,1, give(tcp0) );  
  gsMatrix<> tcp1(2,2); tcp1 << 1, 0, 1 , 1;
  gsBSpline<> * tcurve1 = new gsBSpline<>( 0,1,0,1, give(tcp1) );  
  gsMatrix<> tcp2(2,2); tcp2 << 1, 1, 0 , 1;
  gsBSpline<> * tcurve2 = new gsBSpline<>( 0,1,0,1, give(tcp2) );   
  gsMatrix<> tcp3(2,2); tcp3 << 0, 1, 0 , 0;
  gsBSpline<> * tcurve3 = new gsBSpline<>( 0,1,0,1, give(tcp3) );    
  
  // Then construct a curve loop
  gsCurveLoop<> * tloop = new gsCurveLoop<>();
  tloop->insertCurve( tcurve0 );tloop->insertCurve( tcurve1 );
  tloop->insertCurve( tcurve2 );tloop->insertCurve( tcurve3 );  
  gsPlanarDomain<> * dom= new gsPlanarDomain<>(tloop);
  
  // Define the master NURBS surface: degree 2, Bezier patch
  gsTensorBSpline<2,real_t>::Ptr surf = gsNurbsCreator<>::BSplineSquare();
  
  // Set up the trimmed surface
  gsTrimSurface<> tsurf(surf, dom );
  
  // Print the trim surface
  gsInfo<<"\n"<<"\n";
  //gsInfo<< "I am a "<< tsurf <<"\n";      
  
  
  ///-----------------------------------------------------------------------
  /// --- A shortcut to construct the above trimmed surface
  gsMatrix<> corner = gsMatrix<>(4,3);
    
  corner << 0, 0, 0,
	     1, 0, 0,  
	     1., 1., 0,  
	     0+ .0, 1- .0, 0; 
	     
  gsTrimSurface<> tsurf1(corner, 3, 3, 2);
  
  gsInfo<<"\n"<<"==================== function members ===================="<<"\n"; 
  gsInfo<<"\n"<<"Coords of the vertex 0 of loop 0: \n"<< tsurf1.vertexCoord(0,0);
  
  //gsInfo << "\n"<<"\n"<< " This is a: " << tsurf1<< "\n";
  
  real_t angle1, angle2, angle, angle1T, angle2T, angleT;

  tsurf1.cuttingAngles(0,2, &angle, &angle1, &angle2);
  tsurf1.cuttingAngles(2,0, &angleT, &angle1T, &angle2T);
  gsInfo << "\n" << "Angle: " << angle << " Angle 1: " << angle1 << " Angle 2: " << angle2;
  gsInfo << "\n" << "AngleT: " << angleT << " Angle 1T: " << angle1T << " Angle 2T: " << angle2T;
    
  gsVolumeSegment<> vs;
  real_t cost;
  cost = vs.costAuxiliaryEdge(angle, angle1, angle2, angleT, angle1T, angle2T);
  gsInfo << "\n" << " cost " << cost;
  // For checking, remove when done
   gsInfo << "\n";
  gsMatrix<> cj = tsurf1.derivatives(0);
  gsInfo << "\n" << "corner Jacobian at 0: " << cj << "\n";
  gsInfo << "\n" << "next Tangent coefs: " << tsurf1.UnitTangentCoefs_next(0, cj) << "\n";
  
  gsInfo << "\n" << "prev Tangent coefs: " << tsurf1.UnitTangentCoefs_prev(0, cj) << "\n";
  gsInfo << "\n" << "Tangent of bisectingTrimCurve at 0: " << (tsurf1.TangentCoefs_bisect(0))[0]*cj.col(0) + (tsurf1.TangentCoefs_bisect(0))[1]*cj.col(1) << "\n";
  cj = tsurf1.derivatives(2);
  gsInfo << "\n" << "corner Jacobian at 2: " << cj << "\n";
  gsInfo << "\n" << "Tangent of bisectingTrimCurve at 2: " << (tsurf1.TangentCoefs_bisect(2))[0]*cj.col(0) + (tsurf1.TangentCoefs_bisect(2))[1]*cj.col(1) << "\n";   
  
  gsBSpline<> cc = tsurf1.cuttingCurve(0,2);
  gsInfo << "\n" << "Auxiliary spline edge: " << cc << "\n";
  gsInfo << "\n" << "with control points:\n " << cc.coefs() << "\n";
  
  
  gsInfo << "\n";
  gsVector<> point(2);
  point << 0,0;
  gsInfo << "\n" << "Normal: " << tsurf1.unitNormal(point);
  gsInfo << "\n";
  //for (std::vector< gsGeometry<>* >::iterator it=loopv1.begin();it!=loopv1.end()-1;++it) {delete *it;};
  //delete tdom1;
  

  gsMatrix<> pts(2,1); pts<<0,0;
  gsMatrix<> nm1 = tsurf1.unitNormal(pts);
  gsInfo<<"\n"<<"Normals at: "<<"\n"<< pts <<"\n"<<" is: "<<"\n"<< nm1 <<"\n"; 
  
  gsMatrix<> nm = tsurf1.trimCurTangents(0,0,3);
  gsInfo<<"\n"<<"Tangents along the first trimming curve: "<<"\n"<< nm <<"\n";
  
  gsInfo <<"\n"; return 0;

}
