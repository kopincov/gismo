

#include <iostream>
#include <set>
#include <map>

#include <gismo.h>




using namespace gismo;


int main(int argc, char *argv[])
{
    unsigned np(1000);
    bool plot = false;
    std::string filename = "thbs_01.xml";
    gsTHBSpline<2>::uPtr hbs;
    
    gsCmdLine cmd("Hi, give me a file (.xml) with some hierarchical splines.");
    cmd.addPlainString("filename", "File containing hierarchical splines (.xml).", filename);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // The file data
    gsFileData<>  data( filename );
    if ( data.has< gsTHBSpline<2> >() )
    {
        hbs = data.getFirst< gsTHBSpline<2> >();
    } else {
        gsInfo<<"wrong imput file"<<"\n";
        return 0;
    }
    gsInfo<< "  Got "<< *hbs << "\n";

    // output files

    ///////plotting of the surface//////////
    if(plot)
    {
        gsWriteParaview( *hbs , "gsview", np, false, true);//bool mesh, bool cnet
        return system("paraview gsview.pvd &");
    }

    return 0;
}

