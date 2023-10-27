#include <iostream>
#include <set>
#include <map>


#include <gismo.h>
#include <gismo_dev.h>



using namespace gismo;


int main(int argc, char *argv[])
{
    return 0;   // TODO: line 66: thb->getBsplinePatches(b1, b2, level);
    std::string fn;
    //bool plot;
    //measuring the computational time
    //int clo=clock();
    std::string filename;
    filename= "thbs_15.xml"; //default example
    //filename=  "gsThbs_MTU_00.xml";
    //filename=  "thbs_face_3levels.xml";
    gsTHBSpline<2>::uPtr thb;
    
    gsCmdLine cmd("Hi, give me a file (.xml) with some hierarchical splines.");
    cmd.addPlainString("filename", "File containing hierarchical splines (.xml).", filename);
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
        thb = data.getFirst< gsTHBSpline<2> >();
    }
    gsInfo<< "  Got "<< *thb << "\n";

    /*gsTHBSplineBasis<2>  HB = hbs->basis();
    for(int i = 0; i<=HB.maxLevel(); i++){
      gsInfo<<"B-splines of level "<< i<< ":\n";
      for( set<index_t>::iterator  it = HB.m_xmatrix[i].begin(); it != HB.m_xmatrix[i].end(); it++)
          {
              gsInfo << "("<<*it<<"), ";
          }
      gsInfo <<"\n";
    }*/

  //------------------------------------------------------------------
  // ...
  //------------------------------------------------------------------  
    gsMatrix<index_t> b1, b2;
    gsVector<index_t> level;
    std::vector< gsTensorBSpline<2> > bpatches;

    //gsBSplineBasis<> bsu(ku, ku.degree());
    //gsBSplineBasis<> bsv(kv, kv.degree());
    //gsTensorBSplineBasis<2> tbasis(new gsBSplineBasis<>(bsu), new gsBSplineBasis<>(bsv));
    //gsTensorBSpline<2, real_t> tspline(&tbasis, &cp);
    thb->getBsplinePatches(b1, b2, level);  // TODO: this method is not implemented yet.

    gsInfo<<thb->degree(1);

    //gsMatrix<index_t> nvertices;
    //gsCompactKnotVector<> u_cknot, v_cknot;
    //vector<real_t> u_knot, v_knot, uv_box;
    //uv_box.resize(4);
    //vector<index_t> u_knot_mult, v_knot_mult;

    gsInfo<<"\n"<<"Data for Parasolid"<<"\n"<<"\n";
    //gsInfo<<"u_degree: "<<thb.m_deg[0]<<", v_degree: "<<thb.m_deg[1]<<"\n";
    for (int i = 0; i < level.size(); i++){
        gsInfo<<"\n"<<"Patch number "<<i<<"\n";

        gsInfo<<"bottom-left corner: "<<b1(i,0)<<", "<<b1(i,1)<<"\n";
        gsInfo<<"top-right corner:   "<<b2(i,0)<<", "<<b2(i,1)<<"\n";
        gsInfo<<"level:              "<<level[i]<<"\n";
        //gsInfo<<"n_u_vertices:       "<<nvertices(i,0)<<"\n";
        //gsInfo<<"n_v_vertices:       "<<nvertices(i,1)<<"\n";

        //-------------------------------------------------
        // knots (u direction)
        //-------------------------------------------------
        //u_knot = u_cknot.unique();
        //n_u_knot = u_knot.size();
        //gsInfo<<"n_u_knots:          "<<n_u_knot<<"\n";
        //u_knot_mult.resize(u_knot.size());
        //for(unsigned int j = 0; j < u_knot.size(); j++){
        //    u_knot_mult[j] = u_cknot.u_multiplicityIndex(j);
        //    //gsInfo<<"u_knot: "<<u_knot[j]<<", u_knot_mult - "<<u_knot_mult[j]<<"\n";
        //}
        //-------------------------------------------------
        // knots (u direction)
        //-------------------------------------------------

        //-------------------------------------------------
        // uv_box see below?
        //-------------------------------------------------
        //uv_box[0] = HB.m_bases[HB.max_level]->component(0).knots().unique()[b1(i,0)]; // u_0
        //uv_box[1] = HB.m_bases[HB.max_level]->component(1).knots().unique()[b1(i,1)]; // v_0
        //uv_box[2] = HB.m_bases[HB.max_level]->component(0).knots().unique()[b2(i,0)]; // u_1
        //uv_box[3] = HB.m_bases[HB.max_level]->component(1).knots().unique()[b2(i,1)]; // u_0

        //gsInfo<<"uv_box:             "<<uv_box[0]<<", "<<uv_box[1]<<", "<<uv_box[2]<<", "<<uv_box[3]<<"\n";

        //-------------------------------------------------
        // control points
        //-------------------------------------------------


        //------------------------------------------------------------------
        // parametric grid for the given box
        //------------------------------------------------------------------
        /*gsMatrix<> para = hbs->parameterRange();

        para(0,0) = hbs->basis().m_bases[hbs->basis().max_level]->component(0).knots().unique()[b1(i,0)]; // u_0
        para(1,0) = hbs->basis().m_bases[hbs->basis().max_level]->component(1).knots().unique()[b1(i,1)]; // v_0
        para(0,1) = hbs->basis().m_bases[hbs->basis().max_level]->component(0).knots().unique()[b2(i,0)]; // u_1
        para(1,1) = hbs->basis().m_bases[hbs->basis().max_level]->component(1).knots().unique()[b2(i,1)];

        gsVector<> c0 = para.col(0);
        gsVector<> c1 = para.col(1);
        gsMatrix<> * pts = uniformPointGrid(c0,c1, 10);
        gsInfo<<"\n"<<"Parameter values: (u_i,v_i) = i-th column ("<<pts->cols()<<" points)"<<"\n"<<*pts<<"\n"<<"\n";
        //------------------------------------------------------------------
        // check THB vs B evaluation
        //------------------------------------------------------------------
        gsMatrix<> THB_ev = hbs->eval(*pts);
        gsInfo<<"THB-spline output"<<"\n"<<*THB_ev<<"\n"<<"\n";

        gsMatrix<> Tens_ev = tspline.eval(*pts);
        gsInfo<<"Tensor product B-spline output"<<"\n"<<*Tens_ev<<"\n";

        bool compare_eval= true;
        for(int k = 0; k < THB_ev->cols(); k++){
            for(int j = 0; j < THB_ev->rows();j++){
                if( abs( (*THB_ev)(j,k) - (*Tens_ev)(j,k) ) > 0.000001){
                    gsInfo<<"Check THB-spline vs. B-spline representation: FAILED."<<"\n";
                    compare_eval = false;
                }
            }
        }
        if(compare_eval){
             gsInfo<<"\n"<<"Check THB-spline vs. B-spline representation on "<<pts->cols()<<" points: OK."<<"\n";
        }*/

    }

    return 0;
}
