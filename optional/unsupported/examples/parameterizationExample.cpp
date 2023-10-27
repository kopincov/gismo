// Example used for testing parameterization applying BEM and
//predictor corrector segmentation technique
#include <iostream>
#include <math.h>

#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    gsPlanarDomain<>::uPtr Pdomain;
    index_t n_points(20);
    real_t tolerance = 0.0001; //1e-4;
    index_t smpl(60);
    bool plot = false;
    
    std::string fn = "planar/planarDomainPuzzle1.xml";
    
    gsCmdLine cmd("Segmentation in quadrangular patches of a planar domain ");
    cmd.addString("g","geometry","File containing Geometry (.axl, .txt)", fn);
    cmd.addSwitch("plot", "Plot result with ParaView ", plot);
    cmd.addInt("s","samples", "Number of samples", smpl);
    cmd.addInt("p","n_points", "Number of points traced per curve", n_points);
    cmd.addReal("t","tolerance", "Required accuracy", tolerance);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Pdomain =  gsReadFile<>( fn ) ;
    gsFileData<>  filedata(fn);

    if( filedata.has<gsPlanarDomain<> >() )
        Pdomain = filedata.getFirst< gsPlanarDomain<> >();

    else
    {

        if(filedata.has<gsCurve <> >() )
        {
            Pdomain = memory::convert_ptr<gsPlanarDomain<> >(filedata.getAnyFirst< gsCurve<> >());
        }
        else
        {
            gsWarn<< "Did not find any planar domain or geometry in "<< fn<<", quitting.\n";
            return 1;
        }
    }

    GISMO_UNUSED(n_points);
    GISMO_UNUSED(tolerance);

//    gsMatrix<> corners(2,4);
//    corners<< 0, 1, 0, 1,
//              0, 1, 1, 1;
   // gsInfo<<"\n corner matrix : \n"<<corners; v
//corners.conservativeResize(2,4);
//gsInfo<<"\n corner matrix : \n"<<corners;

    //constructing a square
//    char q;
//    gsTemplate<> square(q,1);
// bool cc = false;

 gsMesh<> par;

// To do: put the segment function to a separate source file (is now inside segment2d_example)
// Pdomain->segment( n_points, tolerance, cc, par);


    int exitCommand = 0;
    if(plot)
    {

        gsWriteParaview(*Pdomain, "quadJaka",smpl);

        gsWriteParaview<>( par, "Parameterization");

        //run: paraview
        exitCommand = system("paraview Parameterization.vtp &");
    }

    return exitCommand;
}


