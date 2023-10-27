/// gsFitting_general.cpp
/// Author:Gabor Kiss
/// for testing the class gsFitting

#include <iostream>
#include <time.h>

#include <gismo.h>
#include <gismo_dev.h>


using namespace gismo;


int main(int argc, char *argv[])
{
    bool plot = false;
    //measuring the computational time
    //int clo=clock();
    std::string filename = "face.xml";
    gsGeometry<>::uPtr hbs;
    int np = 500;
    int nd = 60;
    int  err_type = 1;
    
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
    }

    if ( data.has< gsTensorBSpline<2,real_t>  >() )
    {
        hbs = data.getFirst< gsTensorBSpline<2,real_t>  >();
    }
    gsInfo<< "  Got "<< *hbs << "\n";

    gsMatrix<> para  = hbs->parameterRange();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    //the parameter values for the fitting
    gsMatrix<> pts = uniformPointGrid(c0,c1, nd);
    //gsInfo<< "Parameter values used for fitting: "<<"\n"<<*pts<<"\n";
    //the evaluated values of the original surface
    gsMatrix<> hbs_eval = hbs->eval(pts);
    //gsInfo<<"Original surface points: "<<"\n"<< hbs_eval<<"\n";

    //create the initial basis
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"Creating the initial basis for fitting"<<"\n";
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    int deg_x = 3;
    int deg_y = 3;
    int kn = 2;
    gsKnotVector<> T_KV (0, 1, kn , deg_x+1, 1 ) ;
    //kn = (rand()%num_knots_y)+1;
    kn = 5;

    gsKnotVector<> T_KV1 (0, 1,kn,deg_y+1,1) ;
    gsInfo<<"Knot Vector"<<T_KV<<"\n";
    gsInfo<<"Knot Vector"<<T_KV1<<"\n";

    gsTensorBSplineBasis<2> T_tbasis( T_KV, T_KV1 );
    gsInfo<<"basis created"<<"\n";

    //lvl.push_back(4);
    //q(0,0) = 0; q(1,0) = 0; q(0,1) = 1; q(1,1) = 1;//refine everywhere
    //q(0,0) = 0; q(1,0) = 0; q(0,1) = 0.66; q(1,1) = 1;  lvl.push_back(2);//half face to lvl 1
    //q(0,0) = 0; q(1,0) = 0; q(0,1) = 1; q(1,1) = 0.5;  lvl.push_back(3);//upper half of the face

    //q(0,index) = 0.16; q(1,index) = 0.45; q(0,index+1) = 0.83; q(1,index+1) = 0.83;
    //index+=2;
    //lvl.push_back(1); // level 3




    //q(0,index) = 0.167; q(1,index) = 0.16; q(0,index+1) = 0.82; q(1,index+1) = 0.66;
    //lvl.push_back(1);// level 1-globaly
    //index+=2;

    //q(0,index) = 0.22; q(1,index) = 0.16; q(0,index+1) = 0.7; q(1,index+1) = 0.66;
    //lvl.push_back(2);// level 1-globaly
    //index+=2;

    //q(0,index) = 0; q(1,index) = 0.16; q(0,index+1) = 1; q(1,index+1) = 0.57;
    //lvl.push_back(1);// level 1-globaly
    //index+=2;




    int nboxes = 7;
    gsMatrix<> q(2, 2*nboxes);
    std::vector<index_t> lvl;
    int index = 0;
    /////degree 2/////

    q(0,index) = 0; q(1,index) = 0.0; q(0,index+1) = 1; q(1,index+1) = 1;
    lvl.push_back(1);// level 1-globaly
    index+=2;

    q(0,index) = 0; q(1,index) = 0.28; q(0,index+1) = 1; q(1,index+1) = 0.577;
    index+=2;
    lvl.push_back(2); // level 2

   q(0,index) = 0.16; q(1,index) = 0.7; q(0,index+1) = 0.84; q(1,index+1) = 1;
    index+=2;
    lvl.push_back(2);


    q(0,index) = 0.28; q(1,index) = 0.51; q(0,index+1) = 0.72; q(1,index+1) = 0.84;
    index+=2;
    lvl.push_back(3); // level 3

    q(0,index) = 0.68; q(1,index) = 0.248; q(0,index+1) = 1; q(1,index+1) = 0.577;
    index+=2;
    lvl.push_back(2); // level 2

    q(0,index) = 0; q(1,index) = 0; q(0,index+1) = 1; q(1,index+1) = 0.15;
    index+=2;
    lvl.push_back(2); // level 2

    q(0,index) = 0; q(1,index) = 0.5; q(0,index+1) = 1; q(1,index+1) = 1;
    index +=2;
    lvl.push_back(2); // level 2

    //helps to create the refinement
    gsTHBSplineBasis<2>  THB ( T_tbasis ,q,lvl ) ;

    //gsInfo<<THB.m_cvs[0][2]<<"\n";
    //gsInfo<<THB.m_cvs[1][2]<<"\n";
    //gsInfo<<"knot vector"<< THB.m_cvs[0][0]<<"\n";
    //gsInfo<<"knot vector"<< THB.m_cvs[1][0]<<"\n";
    gsInfo<<"Basis has degree "<< deg_x <<" and "<< deg_y <<". The number of levels is "<< THB.maxLevel()+1 <<"\n";
    gsInfo<<"The tree has "<< THB.tree().size() << " nodes.\n" << "\n";


    //create the fitting object
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"Creating the THB fitting object"<<"\n";
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsFitting<>  fitting(pts, hbs_eval, THB);
    gsInfo<<"fit class created"<<"\n";
    gsMatrix<> results(1,6);
    fitting.compute(0.0000001);
    gsGeometry<> * test;
    test = fitting.result();
    gsTHBSpline<2>  * hbs1 = static_cast< gsTHBSpline<2>  *> (test);
    std::vector<real_t> errors;
    fitting.get_Error(errors, err_type);
    real_t error;
    fitting.computeApproxError(error, 0);
    real_t min = 1000000;
    real_t max = -1000000;
    for(unsigned int j =0; j < errors.size();j++){
        if(errors[j]>max){
            max = errors[j];
        }
        if(errors[j]<min){
            min = errors[j];
        }
    }
    results(0,0) = hbs1->basis().maxLevel();
    results(0,1) = hbs1->basis().size();
    results(0,2) = min;
    results(0,3) = max;
    results(0,5) = error;

    real_t num = 0;
    for(unsigned int j = 0; j < errors.size(); j++){
        if(errors[j]< 0.00001){
            num++;
        }
    }
    results(0,4) = (num*100)/errors.size();
    gsFileData<> newdata;
    newdata << *test ;
    plot = true;
    if(plot){
        newdata.dump("gsThbs_loc_face_4");
        gsWriteParaview( *test , "gsThbs_loc_face_4", np);
    }
    gsMatrix<> hbs1_eval =hbs1->eval(pts);
    gsInfo<<"number of points"<<errors.size()<<"\n";
    plot_errors( hbs_eval, hbs1_eval,errors, "gsThbs_loc_face_error_4");
    gsInfo<<"results"<<results<<"\n"<<"\n";
    gsInfo<<results(0,0);
    gsInfo<<" & "<<results(0,1);
    gsInfo.setf(std::ios::scientific);
    for(int j = 2; j < results.cols();j++){
        gsInfo<<" & "<<results(0,j);
    }
    gsInfo<<"\n"<<"Finished"<<"\n";

/*

    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"Creating the Tensor fitting object"<<"\n";
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsMatrix<> results(1,6);
    gsFitting<> * t_fitting = new gsFitting<>(pts, hbs_eval, T_tbasis);

    t_fitting->compute(0.001);
    gsGeometry<> * t_test;
    t_test = t_fitting->result();
    gsInfo<<*t_test<<"\n";

        gsGeometry<> * test;
    test = t_fitting->result();
    std::vector<real_t> errors;
    t_fitting->get_Error(errors, err_type);
    real_t error;
    t_fitting->computeApproxError(error, 0);
    real_t min = 1000000;
    real_t max = -1000000;
    for(unsigned int j =0; j < errors.size();j++){
        if(errors[j]>max){
            max = errors[j];
        }
        if(errors[j]<min){
            min = errors[j];
        }
    }
    results(0,0) = 0;
    results(0,1) = test->basis().size();
    results(0,2) = min;
    results(0,3) = max;
    results(0,5) = error;

    real_t num = 0;
    for(unsigned int j = 0; j < errors.size(); j++){
        if(errors[j]< 0.00001){
            num++;
        }
    }
    results(0,4) = (num*100)/errors.size();
    gsInfo<<results(0,0);
    gsInfo<<" & "<<results(0,1);
    gsInfo.setf(std::ios::scientific);
    for(int j = 2; j < results.cols();j++){
        gsInfo<<" & "<<results(0,j);
    }

    gsFileData<> t_newdata;
    t_newdata << *t_test ;
    t_newdata.dump("t_fitting_1");
    gsWriteParaview( test , "t_fitting_1", np);
    gsInfo<<"\n"<<"Finished"<<"\n";
*/
  return 0;
}
