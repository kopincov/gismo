
#include <iostream>
#include <set>
#include <map>

#include <gismo.h>



using namespace gismo;


int main(int argc, char *argv[])
{
    unsigned np(1000);
    bool plot = 0;
    //measuring the computational time
    //int clo=clock();
    std::string filename = "thbs_face_3levels.xml";
    memory::unique_ptr< gsTHBSpline<2> > hbs;
    
    gsCmdLine cmd("Hi, give me a file (.xml) with some hierarchical splines.");
    cmd.addPlainString("filename", "File containing hierarchical splines (.xml).", filename);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // The file data
    gsFileData<>  data( filename );
    /*if ( data.has< gsHBSpline<2> >() )
    {
        gsInfo<<"The HB spline functions are not fully functional"<<"\n";
        return 0;
        //hbs = data.getFirst< gsHBSpline<2> >();
    }*/
    if ( data.has< gsTHBSpline<2> >() )
    {
        hbs = data.getFirst< gsTHBSpline<2> >();
    } else {
        gsInfo<<"wrong imput file"<<"\n";
        return 0;
    }
    gsInfo<< "  Got "<< *hbs << "\n";

    gsTHBSplineBasis<2>  HB = hbs->basis();

    HB.printCharMatrix();

  gsInfo<<"\n"<<"Size of the basis "<<HB.size()<<"\n";
  gsMatrix<> anch = HB.anchors();
  gsInfo<<"\n"<<"-----------------------------------------------------------------------";
  gsInfo<<"\n"<<"Greville points check as 2D geometry";
  gsInfo<<"\n"<<"-----------------------------------------------------------------------";
  gsInfo<<"\n"<<"Greville points of THB basis functions (x direction)"<<"\n"<<anch.row(0)<<"\n";
  gsInfo<<"\n"<<"Greville points of THB basis functions (y direction)"<<"\n"<<anch.row(1)<<"\n";
  // Test that the identity function is correctly working in THBS-form
  gsMatrix<> para  = hbs->parameterRange();
  gsVector<> c0 = para.col(0);
  gsVector<> c1 = para.col(1);
  gsMatrix<> pts = uniformPointGrid(c0,c1, 7) ;
  anch.transposeInPlace();
  gsTHBSpline<2> THB_id ( HB, anch );
  gsMatrix<> THB_id_ev;
  THB_id.eval_into( pts, THB_id_ev);
  gsMatrix<> THB_id_der =THB_id.jacobian(pts);
  gsInfo<<"\n"<<"Evaluation points\n"<< pts <<"\n"<< "\nTHB-spline evaluation: check identity function (j-th column => [Px Py])\n"<< THB_id_ev <<"\n";
  //THB_id_der->transposeInPlace();
  gsInfo<<"\n"<<"First derivative: check derivative of identity function (2x2 submatrix for [Px Py] => [dPx/du dPx/dv; dPy/du dPy/dv)\n"<< THB_id_der <<"\n";
  gsMatrix<> THB_id_der2 = THB_id.deriv2(pts);
  //THB_id_der2->transposeInPlace();
  gsInfo<<"\n"<<"Second derivative (i-th column = [dF1/dxx dF1/dyy dF1/dxy; dF2/dxx dF2/dyy dF2/dxy])"<<"\n"<< THB_id_der2;
  gsInfo<<"\n"<<"-----------------------------------------------------------------------";
  gsInfo<<"\n"<<"THB-spline 3D geometry";
  gsInfo<<"\n"<<"-----------------------------------------------------------------------";
  pts = uniformPointGrid(c0,c1, 3);
  gsMatrix<> THB_ev =   hbs->eval( pts);
  gsMatrix<> THB_der = hbs->jacobian(pts);
  //THB_der->transposeInPlace();
  gsMatrix<> THB_der2 = hbs->deriv2(pts);
  gsInfo<<"\n"<<"Evaluation points\n"<< pts<<"\n";
  gsInfo<<"\n"<<"THB-spline evaluation (j-th column => [Px Py Pz])"<<"\n"<< THB_ev<<"\n";
  gsInfo<<"\n"<<"First derivative (3x2 submatrix for [Px Py Pz] => [dPx/du dPx/dv; dPy/du dPy/dv; dPz/du dPz/dv])"<<"\n"<< THB_der<<"\n";
  gsInfo<<"\n"<<"Second derivative (j-th column for [Px Py Pz] => [d2Px/du2 d2Px/dv2 dPx/dudv d2Py/du2 d2Py/dy2 d2Py/dudv dPz/du2 dPz/dv2 dPz/dudv])"<<"\n"<< THB_der2<<"\n";

  gsInfo<<"\n"<<"-----------------------------------------------------------------------";
  gsInfo<<"\n"<<"Parasolid matrix?";
  gsInfo<<"\n"<<"-----------------------------------------------------------------------";
  //
  // P       dP/du     d2P/du2
  // dP/dv   d2P/dudv
  // d2P/dv2
  // ---
  // Px         Py         Pz         dPx/du      dPy/du      dPz/du      d2Px/du2 d2Py/du2 d2Pz/du2
  // dPx/dv     dPy/dv     dPz/dv     d2Px/dudv   d2Py/dudv   d2Pz/dudv
  // d2Px/dv2   d2Py/dv2   d2Pz/dv2
  //
  // P matrix for each point
  gsMatrix<> P(3,9); //P(3,9*pts->cols());
  P.setZero();

  for (int j=0; j < THB_ev.cols(); j++){
      for (int i=0; i < THB_ev.rows(); i++){
          P(0,i)   = THB_ev(i,j);       // P
          P(0,i+3) = THB_der(i,2*j);    // dP/du
          P(0,i+6) = THB_der2(3*i,j);   // d2P/du2
          
          P(1,i)   = THB_der (i,2*j+1);  // dP/dv
          P(1,i+3) = THB_der2(3*i+2,j); // d2P/dudv
          
          P(2,i) = THB_der2(3*i+1,j);   // d2P/dv2
      }
      gsInfo<<"\n"<<"P matrix"<<"\n"<<P<<"\n";
      gsInfo<<"\n"<<"--------------------------";
  }
  gsInfo<<"\n";

  if(plot)
  {
      gsWriteParaview( *hbs , "gsview", np);
      char cmdParaview[100];
      strcpy(cmdParaview,"paraview gsview.pvd\0");
      strcat(cmdParaview," &");
      return system(cmdParaview);
  }
  
  ///////////////////////////////////////////////////////////////////

  //gsMatrix<> eye = gsMatrix<>::Identity(HB.size(), HB.size());
  //gsInfo<<eye<<"\n";
  //gsTHBSpline<2> THB_geometry ( &HB, &eye);
  //gsMatrix<> T_ev_geom_at_anch = THB_geometry.eval(anch);
  //gsInfo<<"\n"<<"Evaluation of the Greville poins with identity matrix "<<"\n"<< T_ev_geom_at_anch<<"\n";
  //delete T_ev_geom;

  return 0;
}

















