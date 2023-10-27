/** @file gsEvaluation_test.cpp

    @brief Tests the evaluation functionality of G+Smo function objects

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <iostream>

using namespace gismo;

// Test functions, defined at the end of this file
int run_test ( gsGeometry<> & geo,   
	       gsMatrix<> const & pts, 
	       bool const & vb = false);

int run_test ( gsBasis<> & basis,
	       gsMatrix<> const & pts,
	       bool const & vb = false);

// Main function
int main(int argc, char *argv[])
{
    std::string filename("");
    bool verbose = false;
    
    gsCmdLine cmd("This test verifies the functionality of evaluation members for G+Smo gsGeometry and gsBasis classes.");
    cmd.addSwitch("verbose", "Prints extended output on the screen", verbose);
    cmd.addPlainString("filename", "File containing the input geometry", filename);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    
    gsInfo<< "This test verifies the functionality of evaluation members of the gsGeometry and gsBasis classes.\n";
    
    gsGeometry<>::uPtr geo;
    gsBasis<>::uPtr basis;
    int result = 0;// 0 == success
  
    if ( ! filename.empty() )
    {
        gsFileData<> fdata(filename);

        geo = fdata.getAnyFirst<gsGeometry<> >();
        basis = fdata.getAnyFirst<gsBasis<> >();

        if (geo) // check if we've got data
        {
            gsInfo<< "  Got "<< *geo << "\n";
            gsMatrix<> para  = geo->parameterRange();
            para.col(0).array() += 0.01;// avoid endpoints
            para.col(1).array() -= 0.01;
            gsVector<> c0 = para.col(0) ;
            gsVector<> c1 = para.col(1) ;
            gsMatrix<> pts = uniformPointGrid( c0, c1 , 16 ) ;
            result += run_test( *geo, pts, verbose );
            result += run_test( geo->basis(), pts, verbose );
        }

        if (basis)
        {
            gsInfo<< "  Got "<< *basis << "\n";
            gsMatrix<> para  = basis->support();
            para.col(0).array() += 0.01;// avoid endpoints
            para.col(1).array() -= 0.01;
            gsVector<> c0 = para.col(0) ;
            gsVector<> c1 = para.col(1) ;
            gsMatrix<> pts = uniformPointGrid( c0, c1 , 16 ) ;
            result += run_test( *basis, pts, verbose );
        }

        return result;
    }
    else
    {
        geo = gsNurbsCreator<>::BSplineFatCircle();// BSpline Curve
        gsMatrix<> para  = geo->parameterRange();
        para.col(0).array() += 0.02;// avoid endpoints
        para.col(1).array() -= 0.01;
        gsVector<> c0 = para.col(0) ;
        gsVector<> c1 = para.col(1) ;
        gsMatrix<> pts = uniformPointGrid( c0, c1 , 3);//geo->basis().size() ) ;
        result += run_test( *geo, pts, verbose );
        result += run_test( geo->basis(), pts, verbose );
    
        geo = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(2));// Planar TensorBSpline patch
        //geo->uniformRefine(2);
        para  = geo->parameterRange();
        para.col(0).array() += 0.01;// avoid endpoints
        para.col(1).array() -= 0.01;
        c0 = para.col(0) ;
        c1 = para.col(1) ;
        pts = uniformPointGrid( c0, c1 , 3);//geo->basis().size() ) ;
        result += run_test( *geo, pts, verbose );
        result += run_test( geo->basis(), pts, verbose );
    
        geo = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(2)); // 3D TensorBSpline patch
        //geo->uniformRefine(2);
        para  = geo->parameterRange();
        para.col(0).array() += 0.01;// avoid endpoints
        para.col(1).array() -= 0.01;
        c0 = para.col(0) ;
        c1 = para.col(1) ;
        pts = uniformPointGrid( c0, c1 ,3);// geo->basis().size() ) ;
        result += run_test( *geo, pts, verbose );
        
        geo = gsReadFile<>("thbs_01.xml");  // 2D hier. truncated patch
        //geo->uniformRefine(2);
        para  = geo->parameterRange();
        para.col(0).array() += 0.01;// avoid endpoints
        para.col(1).array() -= 0.01;
        c0 = para.col(0) ;
        c1 = para.col(1) ;
        pts = uniformPointGrid( c0, c1 , 3);//geo->basis().size() / 2 ) ;
        result += run_test( *geo, pts, verbose );
        result += run_test( geo->basis(), pts, verbose );
    }
  
  // Result=0 : success, else failure

  if (result)
      gsInfo << "\n" << result << " tests failed." << "\n";
  else
      gsInfo << "\nAll tests ok." << "\n";
  return result;
}


int run_test( gsGeometry<> &  geo, const gsMatrix<> & pts, bool const & vb )
{
  gsMatrix<> m1, m2, m3;
  int result = 0;
  real_t error;
  
  gsInfo<<"\n***********  Running for "<< geo ;
  
  gsInfo<<"\n* Parametric dimension: "<< geo.parDim()         <<"\n" ;
  gsInfo<<"* Parametric range:\n"   << geo.parameterRange() <<"\n"  ;
  gsInfo<<"* Spatial dimension: "   << geo.geoDim()         <<"\n"    ;
  
  if (vb)
    gsInfo<<" Evaluation points:\n" << pts  <<"\n";
  else
    gsInfo<<" #Evaluation points: " << pts.cols()  <<"\n";
  
  // For the next lines to work, the minimum requirement is to
  // implement: gsBasis::active_into(u,result). 
  
  // For the next line to work, implement:
  // gsBasis::eval_into(u,result) and
  geo.eval_into(pts,m1);
  if (vb)
    gsInfo<<" Result of evaluation:\n"       << m1  <<"\n";
  
  // For the next line to work, implement:
  // gsBasis::deriv_into(u,result).
  geo.deriv_into(pts, m1 );
  geo.gsFunction<>::deriv_into(pts, m2 ); // This is using finite differences
  if (vb) {
    gsInfo<<" Result of first derivative:\n" << m1 <<"\n";
    gsInfo<<" Result of first derivative by finite differences:\n" << m2 <<"\n";
    gsInfo<<" Difference:\n"<<(m1-m2)<<"\n";
  }
  error = (m1 - m2).array().abs().maxCoeff();
  gsInfo<<" Distance: "<< error <<"\n";
  if (error > 1e-4)
      ++result;
  
  // For the next line to work, implement:
  // gsBasis::deriv2_into(u,result).    
  geo.deriv2_into(pts, m1 );
  geo.gsFunction<>::deriv2_into(pts, m2 ); // This is using finite differences
//  if (vb) {
    gsInfo<<" Result of second derivative:\n" << m1 <<"\n";
    gsInfo<<" Result of second derivative by finite differences:\n" << m2 <<"\n"; 
//  }
  error = (m1 - m2).array().abs().maxCoeff();
  gsInfo<<" Distance: "<< error <<"\n";
  if (error > 1e-4)
      ++result;
  
  return result;
}


int run_test( gsBasis<>        & basis,
              gsMatrix<> const & /*pts*/,
              bool const       & vb )
{
  gsMatrix<>   m1; //, m2, m3;
  gsMatrix<index_t> a1;

  gsInfo<<"\n*********** Running for "<< basis ;

  gsInfo<<"* Parametric dimension: "<< basis.dim()         <<"\n" ;
  gsInfo<<"* Support of the whole basis:\n"   << basis.support() <<"\n"  ;
  
  gsMatrix<> points;
  basis.anchors_into(points);
  if (vb)
    gsInfo<<" Evaluate the basis at:\n"   << points <<"\n";
  else
    gsInfo<<" Evaluate the basis at: "   << points.cols() <<" points.\n";


  // For the next line to work, implement:
  // gsBasis::active_into(u,result).
  basis.active_into(points,a1);
  if (vb) 
  {
    gsInfo<<" Result of active basis functions:\n"<< a1 <<"\n";
  }


  // For the next line to work, implement:
  // gsBasis::eval_into(u,result).
  basis.eval_into(points,m1);
  if (vb) 
  {
    gsInfo<<" Result of evaluation:\n"       << m1 <<"\n";
  }
  
  gsInfo<<" Sum of basis functions: :\n"       << m1.colwise().sum()  <<"\n";

  // For the next line to work, implement:
  // gsBasis::deriv_into(u,result).
  basis.deriv_into(points,m1);
  if (vb) 
  {
    gsInfo<<" Result of first derivatives:\n" << m1 <<"\n";
  }

  // For the next line to work, implement:
  // gsBasis::deriv2_into(u,result).
  basis.deriv2_into(points,m1);
  if (vb) 
  {
    gsInfo<<" Result of second derivative:\n"<< m1 <<"\n";
  }

  unsigned single = a1(0,3);

  // For the next line to work, implement:
  // gsBasis::evalSingle_into(u,result).
  basis.evalSingle_into(single, points,m1);
  if (vb) 
  {
      gsInfo<<" Result of evaluation of basis function "<<single<<":\n"       << m1 <<"\n";
  }
  
  // For the next line to work, implement:
  // gsBasis::derivSingle_into(u,result).
  basis.derivSingle_into(single, points, m1);
  if (vb) 
  {
      gsInfo<<" Result of first derivatives of basis function "<<single<<":\n"       << m1 <<"\n";
  }

  // For the next line to work, implement:
  // gsBasis::deriv2Single_into(u,result).

  //basis.deriv2Single_into(single, points, m1);
  //if (vb) 
  //{
  //    gsInfo<<" Result of second derivatives of basis function "<<single<<":\n"       << m1 <<"\n";
  //}

/*
  // For the next line to work, implement:
  // gsBasis::evalAllDersSingle_into(u,result).
  basis.evalAllDersSingle_into(single, points, 2, m1);
  if (vb) 
  {
      gsInfo<<" Result of all derivatives of basis function "
          << single <<":\n"       << m1 <<"\n";
  }
*/
  
    return 0;
}
