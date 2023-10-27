#include <iostream>
#include <set>
#include <map>

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string fn;
    //bool plot;
    //measuring the computational time
    //int clo=clock();
    std::string filename;
    //filename = "thbs_04.xml";
    //filename = "thbs_08.xml";
    //filename = "thbs_15.xml";
    //filename = "gsThbs_MTU_00.xml";
    filename = "thbs_face_3levels.xml";
    memory::unique_ptr< gsTHBSpline<2> > hbs;
    bool passed = true;
    
    gsCmdLine cmd("Hi, give me a file (.xml) with some hierarchical splines.");
    cmd.addPlainString("filename", "File containing hierarchical splines (.xml).", filename);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // The file data
    gsFileData<>  data( filename );      
    if ( data.has< gsTHBSpline<2> >() ){
        hbs = data.getFirst< gsTHBSpline<2> >();
    } else {
        gsInfo<<"THBSpline not found."<<"\n";
        return 0;
    }
    gsInfo << "  Got "<< *hbs << "\n";

    gsTHBSplineBasis<2>  HB = hbs->basis();

    HB.printCharMatrix();
  
  //------------------------------------------------------------------
  // data for parasolid (BSURF format)
  //------------------------------------------------------------------
  /*
  (BSURF_sf :u_degree          2
            :v_degree          2
            :n_u_vertices      4
            :n_v_vertices      5
            :vertex_dim        3
            :is_rational       nil
            :vertex            (vertex-sf cp0)
            :n_u_knots         3
            :u_knot            #(0 1/2 1)
            :u_knot_mult       #(3  1  3)
            :n_v_knots         4
            :v_knot            #(0 1/3 2/3 1)
            :v_knot_mult       #(3  1  1   3)))))*/
  //------------------------------------------------------------------
  // test getBsplinePatches for Parasolid (multiple patches)
  //------------------------------------------------------------------
  gsMatrix<> cp;
  gsVector<index_t> level;
  gsMatrix<index_t> p1, p2;
  gsKnotVector<> u_cknot, v_cknot;
  std::vector<real_t> u_knot, v_knot, uv_box;
  uv_box.resize(4);
  std::vector<unsigned int> u_knot_mult, v_knot_mult;
  unsigned int n_u_knot, n_v_knot, cpstart = 0, cpend;

  // In the parasolid/lisp interface
  // - HB.getBsplinePatches -> is called to get the list of patches
  // - HB.getBsplinePatchGlobal -> is called on each patch to have the information as below

  gsInfo<<"\n"<<"Call getBoxes from gsHDomain"<<"\n";
  HB.tree().getBoxes(p1,p2,level); // geometry from the file
  gsInfo<<"# of boxes from getBoxes: "<<level.size()<<"\n";

  gsInfo<<"\n"<<"Data for Parasolid"<<"\n"<<"\n";
  gsInfo<<"u_degree: "<<HB.degree(0)<<", v_degree: "<<HB.degree(1)<<"\n";
  for (int i = 0; i < level.size(); i++){

      HB.getBsplinePatchGlobal(p1.row(i), p2.row(i), level[i], hbs->coefs(), cp, u_cknot, v_cknot);

      gsInfo<<"----------------------------------"<<"\n";
      gsInfo<<"*** Patch number "<<i<<" ***";
      gsInfo<<"\n"<<"----------------------------------"<<"\n";
      gsInfo<<"bottom-left corner: "<<p1(i,0)<<", "<<p1(i,1)<<"\n";
      gsInfo<<"top-right corner:   "<<p2(i,0)<<", "<<p2(i,1)<<"\n";
      gsInfo<<"level:              "<<level[i]<<"\n";
      gsInfo<<"n_u_vertices:       "<<u_cknot.size()-u_cknot.degree()-1<<"\n";
      gsInfo<<"n_v_vertices:       "<<v_cknot.size()-v_cknot.degree()-1<<"\n";
      //-------------------------------------------------
      // knots (u direction)
      //-------------------------------------------------      
      gsInfo<<"u_cknot:            "<<u_cknot<<"\n";

      u_knot = u_cknot.unique();
      n_u_knot = u_knot.size();
      gsInfo<<"n_u_knots:          "<<n_u_knot<<"\n";
      u_knot_mult.resize(u_knot.size());
      for(unsigned int j = 0; j < u_knot.size(); j++){
          u_knot_mult[j] = u_cknot.u_multiplicityIndex(j);
          //gsInfo<<"u_knot: "<<u_knot[j]<<", u_knot_mult - "<<u_knot_mult[j]<<"\n";
      }
      //-------------------------------------------------
      // knots (v direction)
      //-------------------------------------------------
      gsInfo<<"v_cknot:            "<<v_cknot<<"\n";

      v_knot = v_cknot.unique();
      n_v_knot = v_knot.size();
      gsInfo<<"n_v_knot:           "<<n_v_knot<<"\n";
      v_knot_mult.resize(v_knot.size());
      for(unsigned int j = 0; j < v_knot.size(); j++){
          v_knot_mult[j] = v_cknot.u_multiplicityIndex(j);
          //gsInfo<<"v_knot - "<<v_knot[j]<<", v_knot_mult - "<<v_knot_mult[j]<<"\n";
      }
      //-------------------------------------------------
      // uv_box
      //-------------------------------------------------
      uv_box[0] = HB.tensorLevel(HB.tree().getMaxInsLevel()).component(0).knots().unique()[p1(i,0)]; // u_0
      uv_box[1] = HB.tensorLevel(HB.tree().getMaxInsLevel()).component(1).knots().unique()[p1(i,1)]; // v_0
      uv_box[2] = HB.tensorLevel(HB.tree().getMaxInsLevel()).component(0).knots().unique()[p2(i,0)]; // u_1
      uv_box[3] = HB.tensorLevel(HB.tree().getMaxInsLevel()).component(1).knots().unique()[p2(i,1)]; // v_1

      gsInfo<<"uv_box:             "<<uv_box[0]<<", "<<uv_box[1]<<", "<<uv_box[2]<<", "<<uv_box[3]<<"\n";
      //-------------------------------------------------
      // control points
      //-------------------------------------------------
      cpend = cpstart + (u_cknot.size()-u_cknot.degree()-1)*(v_cknot.size()-v_cknot.degree()-1) - 1;

      gsInfo<<"cpstart             "<<cpstart<<"\n";
      gsInfo<<"cpend               "<<cpend<<"\n";

      cpstart = cpend + 1;
      //------------------------------------------------------------------
      // tensor-product B-spline patch: basis and geometry
      //------------------------------------------------------------------
      gsTensorBSplineBasis<2, real_t> tbasis(u_cknot, v_cknot);
      gsTensorBSpline<2, real_t> tbspline(tbasis, give(cp));

      //------------------------------------------------------------------
      // parametric grid for the given box
      //------------------------------------------------------------------
      gsMatrix<> para = hbs->parameterRange();

      // the correction term is in need in case of discontinuity to
      // avoid the evaluation on the knot value
      para(0,0) = uv_box[0]+0.000001; // u_0
      para(1,0) = uv_box[1]+0.000001; // v_0
      para(0,1) = uv_box[2]-0.000001; // u_1
      para(1,1) = uv_box[3]-0.000001; // v_1

      gsVector<> c0 = para.col(0);
      gsVector<> c1 = para.col(1);
      gsMatrix<> pts = uniformPointGrid(c0,c1, 10);
      gsInfo<<"\n"<<"Parameter values: (u_i,v_i) = i-th column ("<<pts.cols()<<" points)"<<"\n"<<pts<<"\n"<<"\n";
      //------------------------------------------------------------------
      // check THB vs B evaluation
      //------------------------------------------------------------------
      gsMatrix<> THB_ev = hbs->eval(pts);
      gsInfo<<"THB-spline output"<<"\n"<< THB_ev<<"\n"<<"\n";

      gsMatrix<>  Tens_ev = tbspline.eval(pts);
      gsInfo<<"Tensor product B-spline output"<<"\n"<< Tens_ev<<"\n";

      bool compare_eval= true;
      for(int i2 = 0; i2 < THB_ev.cols();i2++){
          for(int j = 0; j < THB_ev.rows();j++){
              if( math::abs( THB_ev(j,i2) - Tens_ev(j,i2) ) > 0.000001)
              {
                  gsInfo<<"Check THB-spline vs. B-spline representation: FAILED."<<"\n";
                  compare_eval = false;
                  passed = false;
              }
          }
      }
      if(compare_eval){
           gsInfo<<"\n"<<"Check THB-spline vs. B-spline representation on "<<pts.cols()<<" points: OK."<<"\n"<<"\n";
      }
      gsInfo<<tbspline;
  }

  gsMultiPatch<> multipatch = HB.getBsplinePatchesToMultiPatch(hbs->coefs());
//  gsFileData<> newdata;
//  newdata << multipatch ;
//  newdata.dump("multipatch_THB");

  gsWriteParaview( multipatch,"multipatch_THB" , 100);

  return passed ? 0 : 1;
}










