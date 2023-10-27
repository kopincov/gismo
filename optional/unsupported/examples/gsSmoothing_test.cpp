/// gsSmoothing_test.cpp
/// Author:Mario Kapl
/// for testing the class gsCurvatureSmoothing
/// takes the point cloud of a curve with the parameter values and the desired knot vector
/// and computes a starting curve with the help of fitting  and smoothes this curve by using
/// Hadenfeld's algorithm and total variation

#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsModeling/gsCurvatureSmoothing.h>
#include <gsUtils/gsUtils.h>

/*
#include <gsCore/gsLinearAlgebra.h>

// curves
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsNurbs.h>
#include <gsBezier/gsBezier.h>

// surfaces or volumes ..
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsBezier/gsTensorBezier.h>

// Read
#include <gsIO/gsCmdLine.h>
#include <gsIO/gsReadFile.h>

// Multipatch container
#include <gsCore/gsMultiPatch.h>

// write to paraview files
#include <gsIO/gsWriteParaview.h>

//closed curve fitting procedure
#include <gsModeling/gsCurveFitting.h>
#include <gsModeling/gsCurvatureSmoothing.h>

// timing
#include <gsUtils/gsStopwatch.h>
*/


using namespace gismo;


int main(int argc, char *argv[])
{
  std::string filename = "fitting_default.xml";
  bool   plot   = false; //for plotting
  index_t    iterationsTV=0; //how many iterations in the total variation part
  real_t lamda = 0.9   ; //the starting lamda in the backtracking line search method in the gradient descent method in the total variation algorithm
  real_t tau = 0.9     ; // the decreasing factor in the backtracking line search method
  index_t smooth_degree = 3; // degree of Hadenfeld
  memory::unique_ptr< gsCurveFitting<> > fitting;

  gsCmdLine cmd("Hi, give me a file (.xml) with some CurveFitting data.");
  cmd.addPlainString("filename", "File containing the input", filename);
  cmd.addInt("i", "iterations",
             "Number of iterations for total variation smoothing", iterationsTV);
  cmd.addReal("l", "lamda",
              "the starting lamda in the backtracking line search method",
              lamda);
  cmd.addReal("t", "tau",
              "the decreasing factor in the backtracking line search method",
              tau);
  cmd.addInt("d", "degree",
             "Smoothing degree for Hadenfeld", smooth_degree);
  cmd.addSwitch("plot", "Plot result in ParaView format", plot);


  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  // The file data
  gsFileData<>  data( filename );
  if ( data.has< gsCurveFitting<> >() )
  {
    fitting = data.getFirst< gsCurveFitting<> >();
    gsInfo<< "  Got a curve fitting problem "<< *fitting << "\n";
  }

  gsStopwatch time;

  // compute the initial curve via fitting
  fitting->setClosedCurve(true); // we need a closed curve
  fitting->compute();
  const gsBSpline<> & curve = fitting->curve();
  gsMatrix<> param_values=fitting->returnParamValues();
  gsMatrix<> points=fitting->returnPoints();

  gsInfo  << "Time of fitting: " << time << "\n";

  gsCurvatureSmoothing<> smoothing(curve, param_values, points);

  // output for Mathematica for the initial curve (after fitting) for the smoothing process
  filename = gsFileManager::getTempPath() + "/output_init.txt"; //default example
  std::ofstream init(filename.c_str() );
  smoothing.write(init);

  real_t error;
  smoothing.computeApproxError(error);
  gsInfo<<"The approximation error of the fitted curve is: "<< error  <<"\n";

  smoothing.computeApproxErrorL2(error);
  gsInfo<<"The L^2-norm approximation error of the fitted curve is: "<< error  <<"\n";

  smoothing.computeApproxErrorLMax(error);
  gsInfo<<"The LMax-norm approximation error of the fitted curve is: "<< error  <<"\n";

  smoothing.computeCurvatureError(error);
  gsInfo<<"The curvature error of the fitted curve is: "<< error  <<"\n";

  gsInfo << "\n";

  gsInfo <<"We are starting to smooth the resulting curve" << "\n";

  gsInfo << "\n";

  gsInfo <<"First: Smoothing by using total variation" << "\n";

  time.restart();

  smoothing.smoothTotalVariation(8,1,lamda,tau,iterationsTV);

  gsInfo  << "Time of smoothing by using total variation: " << time << "\n";

  // output for Mathematica of the smoothed curve after using total variation
  filename = gsFileManager::getTempPath() + "/output_variation.txt"; //default example
  std::ofstream variation(filename.c_str());
  smoothing.write(variation);

  smoothing.computeApproxError(error);
  gsInfo<<"The approximation error of the smoothed curve after using total variation is: "<< error  <<"\n";

  smoothing.computeApproxErrorL2(error);
  gsInfo<<"The L^2-norm approximation error of the smoothed curve after using total variation is: "<< error  <<"\n";

  smoothing.computeApproxErrorLMax(error);
  gsInfo<<"The LMax-norm approximation error of the smoothed curve after using total variation is: "<< error  <<"\n";

  smoothing.computeApproxErrorCoef(error);
  gsInfo<<"The maximal coefficient approximation error of the smoothed curve after using total variation is: "<< error  <<"\n";

  smoothing.computeCurvatureError(error);
  gsInfo<<"The curvature error of the smoothed curve after using total varation is: "<< error  <<"\n";

  gsInfo << "\n";

  gsInfo <<"Second: Smoothing by using Hadenfeld's algorithm" << "\n";

  time.restart();

  //smoothing of the curve by using Hadenfeld's algorithm
  gsVector<index_t> iterated;
  smoothing.smoothHadenfeld(smooth_degree,0.144,400,32000,iterated,false);

  //gsInfo << iterated << "\n";
  gsInfo  << "Time of smoothing by using Hadenfeld's algorithm: " << time << "\n";

  // output for Mathematica of the smoothed curve after using Hadenfeld's algorithm
  filename = gsFileManager::getTempPath() + "/output_haden.txt"; //default example
  std::ofstream haden(filename.c_str());
  smoothing.write(haden);

  const gsBSpline<> & smooth_curve = smoothing.curveSmooth();

  smoothing.computeApproxError(error);
  gsInfo<<"The approximation error of the smoothed curve after using Hadenfeld's algorithm is: "<< error  <<"\n";

  smoothing.computeApproxErrorL2(error);
  gsInfo<<"The L^2-norm approximation error of the smoothed curve after using Hadenfeld's algorithm is: "<< error  <<"\n";

  smoothing.computeApproxErrorLMax(error);
  gsInfo<<"The LMax-norm approximation error of the smoothed curve after using Hadenfeld's algorithm is: "<< error  <<"\n";

  smoothing.computeApproxErrorCoef(error);
  gsInfo<<"The maximal coefficient approximation error of the smoothed curve after using Hadenfeld's algorithm is: "<< error  <<"\n";

  smoothing.computeCurvatureError(error);
  gsInfo<<"The curvature error of the smoothed curve after using Hadenfeld's algorithm is: "<< error  <<"\n";

  // Print the Bspline curve
  gsInfo<< "The resulting curve of the fitting + smoothing problem is a "<< smooth_curve <<"\n";



  // Output a paraview file
  gsWriteParaview( smooth_curve , "bsplinecurve", 1000);


  //Write back to a new file "dump_write.xml"
  gsFileData<> newdata;
  newdata << smooth_curve;
  newdata.dump("dump_write");

  // check if we have given --plot in the command line otherwise no plotting
  if(plot){
    // Call paraview on exit
    char cmdParaview[100];
    strcpy(cmdParaview,"paraview bsplinecurve.vts\0");
    strcat(cmdParaview," &");
    return system(cmdParaview);
  }
  else {
      return 0;
  }
}
