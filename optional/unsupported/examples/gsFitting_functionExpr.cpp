/// gsFitting_general.cpp
/// Author:Gabor Kiss
/// for testing the class gsFitting

#include <iostream>
#include <algorithm>
#include <time.h>

#include <gismo.h>


;
using namespace gismo;

int main(int argc, char *argv[])
{
    //reading geometry from a file to generate the points and parameter values
    std::string fn;
    bool plot = false;
    //measuring the computational time
    //int clo=clock();
    //std::string filename;
    //gsTHBSpline<2> * hbs;
    index_t np = 100;
    index_t gridpoints = 200;
    index_t  err_type = 1;
    index_t iter = 3;
    index_t function = 7;
    index_t deg_x = 3;
    index_t deg_y = 3;
    real_t lambda = 0.000000001;
    real_t err_threshold = 0.000001;
    real_t tolerance = 0.000001;
    index_t kn = 4;
    index_t kn1 = 4;
    gsMatrix<> para(2,2);
    gsFunctionExpr<> f;
    std::string mfn("file_new_domain_new_insertion");
    index_t extension = 2;
    
    gsCmdLine cmd("Generates random hierarchy configuration and test the partition of unity");
    cmd.addInt("f", "function", "Number of the function", function);
    cmd.addSwitch("plot","Plot result in ParaView format", plot);
    cmd.addInt("p", "points","Number of points for plotting the output", np);
    cmd.addInt("g", "gridpoints", "number of sample points in the grid to be fitted", gridpoints);
    cmd.addInt("e", "error", "error type", err_type);
    cmd.addInt("i", "iter","number of iterations", iter);
    cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
    cmd.addInt("y", "deg_y","degre_in y direction", deg_y);
    cmd.addReal("s", "lambda", "smoothing coefficient", lambda);
    cmd.addReal("t", "threshold", "error threshold", err_threshold);
    cmd.addString("n", "filename", "output filename", mfn);
    cmd.addInt("q", "extension", "extension", extension);
    cmd.addReal("o", "tolerance", "error tolerance", tolerance);
    cmd.addInt("a", "kn_x", "inner knots in x direction", kn);
    cmd.addInt("b", "kn_y", "inner knots in y direction", kn1);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //function = 2;
    
    if (gridpoints < 1)
    {
        std::cerr << "Number of grid points must be positive.\n"; return -1;
    }
   if (np < 1)
    {
        std::cerr << "Number of points must be positive.\n"; return -1;
    }
     if (deg_x < 1)
    {
        std::cerr << "Degree x must be positive.\n"; return -1;
    }
    if (deg_y < 1)
    {
        std::cerr << "Degree y must be positive.\n"; return -1;
    }
    if (kn < 0)
    {
        std::cerr << "Number of inner knots must be non negative.\n"; return -1;
    }
    if (kn1 < 0)
    {
        std::cerr << "Number of inner knots must be non negative.\n"; return -1;
    }
    if (extension < 0)
    {
        std::cerr << "Extension must be non negative.\n"; return -1;
    }


    switch (function) {
    case 1:
        f = gsFunctionExpr<>("(exp(52*sqrt((10*x-2)^2+(10*y-3)^2)))^(-1)",2);//one peak
        para(0,0) = -1;  para(1,0) = -1;  para(0,1) = 1;  para(1,1) = 1;
        gsInfo<<"Source function "<< f <<".\n" << "\n";
        break;
    case 2:
        f = gsFunctionExpr<>("y/2*(cos(4*(x^2+y-1)))^4",2);//waves
        para(0,0) = 0;  para(1,0) = 0;  para(0,1) = 1;  para(1,1) = 1;
        gsInfo<<"Source function "<< f <<".\n" << "\n";
        break;
    case 3:
        f = gsFunctionExpr<>("3/4*exp(-((9*x-2)^2 + (9*y-2)^2)/4)+3/4*exp(-((9*x+1)^2)/49 - (9*y+1)/10)+1/2*exp(-((9*x-7)^2 + (9*y-3)^2)/4)-1/5*exp(-(9*x-4)^2 - (9*y-7)^2)",2);//hils
        para(0,0) = 0;  para(1,0) = 0;  para(0,1) = 1;  para(1,1) = 1;
        gsInfo<<"Source function "<< f <<".\n" << "\n";
        break;
    case 4:
        f = gsFunctionExpr<>("sin((x*x+y*y+2/(5*pi))^(-1))",2);//crater
        para(0,0) = -1;  para(1,0) = -1;  para(0,1) = 1;  para(1,1) = 1;
        gsInfo<<"Source function "<< f <<".\n" << "\n";
        break;
    case 5:
        f = gsFunctionExpr<>("(exp(2*sqrt((x)^2+(y)^2)))^(-1)",2);//cusp
        para(0,0) = -1;  para(1,0) = -1;  para(0,1) = 1;  para(1,1) = 1;
        gsInfo<<"Source function "<< f <<".\n" << "\n";
        break;
    case 6:
        f = gsFunctionExpr<>("(exp(2*sqrt((10*x+3)^2+(10*y-3)^2)))^(-1) + (exp(2*sqrt((10*x-3)^2+(10*y+3)^2)))^(-1)",2);//2 peaks
        para(0,0) = -1;  para(1,0) = -1;  para(0,1) = 1;  para(1,1) = 1;
        gsInfo<<"Source function "<< f <<".\n" << "\n";
        break;
    case 7:
        f = gsFunctionExpr<>("(1.5*exp(sqrt((10*x-3)^2+(10*y-3)^2)))^(-1)+ (1.5*exp(sqrt((10*x+3)^2+(10*y+3)^2)))^(-1) + (1.5*exp(sqrt((10*x)^2+(10*y)^2)))^(-1)",2);//3 peaks
        para(0,0) = -1;  para(1,0) = -1;  para(0,1) = 1;  para(1,1) = 1;
        gsInfo<<"Source function "<< f <<".\n" << "\n";
        break;
    case 8:
        f = gsFunctionExpr<>("(x+y)/2+sqrt(((x-y)/2)^2)",2);//Rvachev functions (R-functions)
        para(0,0) = 0;  para(1,0) = 0;  para(0,1) = 1;  para(1,1) = 1;
        gsInfo<<"Source function "<< f <<".\n" << "\n";
        break;
    default:
        gsInfo<<"Unknown function, please pick one of the functions 1 - 8"<<"\n";
        return 0;
        break;
    }
    /*
       case 8 % 4 peaks
           fun = (exp(1.5*sqrt((10*x-6).^2+(10*y-6).^2))).^(-1) ...
                + (exp(1.5*sqrt((10*x+6).^2+(10*y+6).^2))).^(-1) ...
                + (exp(1.5*sqrt((10*x-2).^2+(10*y-2).^2))).^(-1) ...
                + (exp(1.5*sqrt((10*x+2).^2+(10*y+2).^2))).^(-1);

           FUN = (exp(1.5*sqrt((10*X-6).^2+(10*Y-6).^2))).^(-1) ...
                + (exp(1.5*sqrt((10*X+6).^2+(10*Y+6).^2))).^(-1) ...
                + (exp(1.5*sqrt((10*X-2).^2+(10*Y-2).^2))).^(-1) ...
                + (exp(1.5*sqrt((10*X+2).^2+(10*Y+2).^2))).^(-1);
       case 9 % 5 peaks
           fun = (exp(2*sqrt((10*x-7).^2+(10*y-7).^2))).^(-1) ...
               + (exp(2*sqrt((10*x+7).^2+(10*y+7).^2))).^(-1) ...
               + (exp(2*sqrt((10*x-3).^2+(10*y-3).^2))).^(-1) ...
               + (exp(2*sqrt((10*x+3).^2+(10*y+3).^2))).^(-1) ...
               + (exp(2*sqrt((10*x).^2+(10*y).^2))).^(-1);
   end */

    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> pts = uniformPointGrid(c0,c1, gridpoints);
    //gsInfo<<*pts<<"\n";
    gsMatrix<> f_eval;
    f.eval_into(pts, f_eval);
    gsMatrix<> data_points;
    //create data points from parameter values and evaluated values
    data_points = gsMatrix<>(3,f_eval.cols());
    data_points.row(0) = pts.row(0);
    data_points.row(1) = pts.row(1);
    data_points.row(2) = f_eval.row(0);
    //gsInfo<<*data_points<<"\n";
    //create the initial basis
    /**/gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"Creating the initial basis for fitting"<<"\n";
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";


    gsKnotVector<> T_KV (para(0,0), para(0,1), kn, deg_x+1, 1 ) ;
    gsKnotVector<> T_KV1 (para(1,0), para(1,1),kn1, deg_y+1,1) ;
    gsInfo<<"Knot Vector"<<T_KV<<"\n";
    gsInfo<<"Knot Vector"<<T_KV1<<"\n";

    gsTensorBSplineBasis<2> T_tbasis( T_KV, T_KV1 );

    // write the underlying geometry to a file if desired
    if (0)
    {
        gsWrite( *( gsNurbsCreator<>::BSplineRectangleWithPara(para(0,0), para(1,0), para(0,1), para(1,1)) ),
                 mfn + "_geo.xml" );
    }

    //helps to create the refinement
    gsTHBSplineBasis<2>  THB ( T_tbasis ) ;

    gsInfo<<"Basis has degree "<< deg_x <<" and "<< deg_y <<". The number of levels is "<< THB.maxLevel()+1<<"\n";
    gsInfo<<"The tree has "<< THB.treeSize() << " nodes.\n" << "\n";


    //create the fitting object
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"Creating the THB fitting object"<<"\n";
    gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"error sthreshold "<< err_threshold<<"\n";
    gsInfo<<"lambda "<< lambda<<"\n";
    gsInfo<<"tolerance "<< tolerance<<"\n";
    std::vector<unsigned> ext;
    ext.push_back(extension);
    ext.push_back(extension);
    gsHFitting<2, real_t> ref( pts, data_points , THB, 0.01, ext,lambda);
    gsMatrix<> results(iter+1,7);
    gsStopwatch time;

    for(int i = 0; i <= iter; i++)
    {
        gsInfo<<"iteration "<<i<<" start."<<"\n";

        time.restart();
        ref.nextIteration(tolerance, err_threshold);
        results(i,6) = time.stop();

        gsTHBSpline<2>  * hbs1 =
            static_cast< gsTHBSpline<2>  *> ( ref.result() );

        real_t error;
        ref.computeApproxError(error, 0);
        results(i,0) = hbs1->basis().size();
        results(i,1) = hbs1->basis().maxLevel();
        results(i,2) = ref.minPointError();
        results(i,3) = ref.maxPointError();
        results(i,5) = error;
        const std::vector<real_t> & errors = ref.pointWiseErrors();
        results(i,4) = 100.0 * ref.numPointsBelow(tolerance)/errors.size();

        gsInfo<<"iteration "<<i<<" end."<<"\n";

        plot = true;
        if ( plot )
        {
            std::string mfn1 = mfn;
            std::stringstream ss;
            ss<<mfn1<<i;

            gsWrite(*hbs1, ss.str());

            gsWriteParaview(*hbs1, ss.str(), np, 0, 1);

            //std::string mfn2("gsThbs_global_error_cusp_1");
            //std::stringstream ss2;
            //ss2<<mfn2<<i;
            //gsMatrix<> hbs1_eval;
            //hbs1->eval_into(*pts, hbs1_eval);
            //plot_errors( data_points, hbs1_eval, errors, ss2.str() );
        }
        gsInfo.setf(std::ios::scientific);
        for(int j = 0; j <= i; j++){
	  gsInfo<< j <<" & " << cast<real_t,int>(results(j,0))<<" & "<< cast<real_t,int>(results(j,1))<<" & "<<results(j,2)<<" & "<<results(j,3)<<" & "<<results(j,4)<<" & "<<results(j,5)<<" & "<<results(j,6)<<"\n";
        }
        if(results(i,3) < tolerance)
        {
            break;
        }
    }

    gsInfo<<"\n";
    gsInfo<<"error type "<< err_type<<"\n";
    gsInfo<<"extension "<< ext[0]<<" "<<ext[1]<<"\n";
    gsInfo<<"Basis has degree "<< deg_x <<" and "<< deg_y <<
        ". The number of levels is "<< THB.maxLevel()+1 <<"\n";
    gsInfo<<"The tree has "<< THB.treeSize() << " nodes.\n" << "\n";
    gsInfo<<"error sthreshold "<< err_threshold<<"\n";
    gsInfo<<"lambda "<< lambda<<"\n";
    gsInfo<<"\n"<<"Finished"<<"\n";

    return 0;
}

