/** @file gsFitting_test.cpp

    @brief File testing the gsCurveFitting class.

    Takes the point cloud of a curve with the parameter values and the
    desired knot vector and computes a curve with the help of fitting
    default type of used curve is a closed curve.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Kapl
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <iostream>


using namespace gismo;

int main(int argc, char *argv[])
{
  std::string filename = "fitting_default.xml"; // for the input file which should be read
  bool plot = false; //for plotting
  memory::unique_ptr< gsCurveFitting<> > fitting;
  gsCmdLine cmd("Hi, give me a file (.xml) with some CurveFitting data.");
  cmd.addPlainString("filename", "File containing input  (.xml).", filename);
  cmd.addSwitch("plot", "Plot result in ParaView format", plot);
  
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
  // The file data
  gsFileData<>  data( filename );
  if ( data.has< gsCurveFitting<> >() )
  {
    fitting = data.getFirst< gsCurveFitting<> >();
    gsInfo<< "  Got a curve fitting problem "<< *fitting << "\n";
  }

  fitting->setClosedCurve(true); // closed or open curve

  //measuring the computational time
  gsStopwatch clo;
  fitting->compute();
  double fittime = clo.stop();

  const gsBSpline<> & curve = fitting->curve();

  //computes the approximation error
  real_t error;
  fitting->computeApproxError(error);
  gsInfo<<"The approximation error of the fitted curve is: "<< error  <<"\n";

  // Print the Bspline curve
  gsInfo<< "The resulting curve of the fitting problem is a "<< curve <<"\n";

  // Output a paraview file
  gsWriteParaview( curve , "bsplinecurve", 1000);

  gsInfo  << "Time of computation in milliseconds: " << fittime * 1000 << "\n";

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
