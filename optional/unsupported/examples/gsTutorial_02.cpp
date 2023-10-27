// Tutorial to make a BSpline curve

#include <gismo.h>

#include <iostream>


using namespace gismo;


int main(int argc, char *argv[])
{

    bool plot = false; // If set to true, paraview file is generated and launched on exit
    
    gsCmdLine cmd("Tutorial 01 shows the use of BSpline curves.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  // This is a pointer to some geometry object
  gsGeometry<>::uPtr geom_1, geom_2;
  
  gsMultiPatch<> mp;
  
  // Make a BSpline curve
  gsKnotVector<> KV(0.0,1.0,1,3);//start,end,interior knots, start/end multiplicites of knots
  gsInfo<< "Knot Vector : =" << KV << "\n"; 
  
  gsMatrix<> coefs(4,3);
  
   // Create a first patch i.e a BSpline Curve
  coefs << 0,0,0,  1,2,3, 2,1,4, 4,4,4 ;
  gsInfo<< "coefs : = \n" << coefs << "\n"; 
  geom_1 = gsGeometry<>::uPtr(new gsBSpline<>( KV, give(coefs)));
  //   <---- after this line the coefs is empty variable! 
  
  gsInfo<< "Knot Vector : =" << KV << "\n";
  
  mp.addPatch(give(geom_1));
   
  // Create a second patch i.e a BSpline Curve
  coefs.resize(4, 3);
  coefs << 0,-1,0,  1, .5, 3, 0,1,4, 1,3,4 ;
  gsInfo<< "coefs : = \n" << coefs << "\n";
  geom_2 = gsGeometry<>::uPtr(new gsBSpline<>( KV, give(coefs)));
  
  mp.addPatch(give(geom_2));
  
   // Print the Bspline curve
  gsInfo<< "I am a "<< mp <<"\n";

   // Output a paraview file
  if (plot) 
      gsWriteParaview( mp , "paraviewout", 100);
  
  if (plot) 
  {
      //Call paraview on exit
      char cmdParaview[100];
      if (mp.nPatches() == 1)
          strcpy(cmdParaview,"paraview paraviewout.vts\0");
      else
          strcpy(cmdParaview,"paraview paraviewout.pvd\0");
      strcat(cmdParaview," &");
      return system(cmdParaview);
  }
  else
  {
      return 0;
  }
}
