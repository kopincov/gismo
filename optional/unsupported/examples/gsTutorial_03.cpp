// Tutorial to create and view a Tensor BSpline  

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
    gsGeometry<> *geom_1, *geom_2;

    gsMultiPatch<> mp;


    // Make a BSpline curve
    gsKnotVector<> KV_1(0,1,1,3);// start,end,interior knots, start/end multiplicites of knots
    gsKnotVector<> KV_2(0,1,1,3); // redundant but ok.
    gsInfo<<"Knot Vector 1"<< KV_1<< "\n";
    gsInfo<<"Knot Vector 2"<< KV_2<< "\n";


    gsMatrix<> coefs(4,3);
    coefs << 0,0,0,  1,2,3, 2,1,4, 4,4,4 ;
    gsInfo<<"coefs = \n "<< coefs << "\n";

    geom_1 = new gsBSpline<>( KV_1, coefs);
    gsInfo<< "geom_1 = " << geom_1 << "\n";

    geom_2 = new gsBSpline<>( KV_2, coefs);
    gsInfo<< "geom_2 = " << geom_2 << "\n";

    // create Tensor product of the BSpline  
    coefs.resize(16,3);
    coefs.setZero();  
    coefs.row(0)<<   1,0,0;
    coefs.row(4)<<   0,1,1;
    coefs.row(8)<<   1,0,1;
    coefs.row(15)<<  0,1,0;
    typename gsTensorBSpline<2>::uPtr T_geom(new gsTensorBSpline<2>( KV_1, KV_2, give(coefs)));
    gsInfo<< "T_geom = " << *T_geom << "\n";


    mp.addPatch(give(T_geom));  // push_back takes only geometry Pointer

    // Interfaces if its domain
    //    mp->addInterface(geom_1, boundary::east, geom_2, boundary::west );

    coefs.resize(16, 3);
    coefs.setZero();
    coefs.row(0)<< -1,0,0;
    coefs.row(4)<< 0, -1, -1;
    coefs.row(8)<< -1,0,-1;
    coefs.row(15)<< 0,-1,0 ;
    gsInfo<<"coeffs = \n "<< coefs << "\n";
    typename gsTensorBSpline<2>::uPtr T_geom_1(new gsTensorBSpline<2>(KV_1, KV_2, give(coefs)));

    mp.addPatch(give(T_geom_1));

    // Print the Tensor Bspline Surface
    gsInfo << "I am a "<< mp <<"\n";

    // Output a paraview file
    if (plot)
        gsWriteParaview( mp , "paraviewout", 100);

    delete geom_1;
    delete geom_2;

    if (plot) 
    {
        // Call paraview on exit
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
