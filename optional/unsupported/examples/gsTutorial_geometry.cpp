// Tutorial to create a new geometry from another geometry 
// from an existing geometry  

#include <iostream>

#include <gismo.h>



using namespace gismo;


int main(int argc, char *argv[])
{
    
    bool plot = false; // If set to true, paraview file is generated and launched on exit
    std::string fn ("planar/two_squares.xml");

    int deg = 3;
    int numKnots = 2;
    
    gsCmdLine cmd("Tutorial_Geometry  creates a new geometry from an existing one.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // This is a pointer to some geometry object
    gsGeometry<>::uPtr BS = gsReadFile<>( fn );
    // Print the Bspline curve
    gsInfo<< "I am a "<< *BS <<"\n";
      
    // Make Knot Vector 
    gsKnotVector<> KV(0,1,numKnots,deg+1,deg);//start,end,interior knots, start/end multiplicites of knots, interior mults
    gsTensorBSplineBasis<2> tbsp (new gsBSplineBasis<>(KV), new gsBSplineBasis<>(KV) );
    
    gsGeometry<>::uPtr A = tbsp.interpolateAtAnchors( BS->eval(tbsp.anchors()) );
    // Print the Bspline curve
    gsInfo<< "New geometry: "<< *A <<"\n";
  
    if (plot) 
    {
        // Output a paraview file
        gsWriteParaview( *A , "newgeometry", 100);
        
        // Call paraview on exit
        char cmdParaview[100];
        strcpy(cmdParaview,"paraview newgeometry.vts\0");
        strcat(cmdParaview," &");  
        return system(cmdParaview);
    }
    
    return 0;
}
