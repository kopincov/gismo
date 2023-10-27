
#include <iostream>
#include <vector>

#include <gismo.h>
#include <gismo_dev.h>



using namespace gismo;


int main(int argc, char *argv[])
{   
    // Input options
    bool plot       = false;
    index_t samples(1000);

    // Multipatch object
    gsMultiPatch<> mp;

    // Pointer to a Pde
    memory::unique_ptr< gsPde<> > pde;
    
    std::string fn("surfaces/sphere1.xml");
    std::string fn_pde;
    
    gsCmdLine cmd("Creates a paraview file from a given field on a geometry.");
    cmd.addInt("s","samples", "Number of samples to use for plotting", samples);
    cmd.addString("q","pde","File containing a poisson PDE (.xml)", fn_pde);
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)", fn);
    cmd.addSwitch("plot", "Call ParaView on exit", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    gsReadFile<>(fn, mp);

    if ( ! fn_pde.empty() )
        pde = gsReadFile<>(fn_pde) ;


    if ( mp.coDim() != 1 )
    {
        gsWarn << "Co dimension one is assumed for now.\n";
    }

    if ( pde )
    {
        if ( pde->solutionGiven() )
        {
            gsInfo<<"Outputing Pde solution field to paraview file...\n";
            gsField<> out1( mp, *pde->solution(), false ) ;
            gsWriteParaview<>( out1, "field_output", samples);

            gsInfo<<"Outputing Pde solution gradient field to paraview file...\n";
            gsField<> out2 = gsFieldCreator<>::gradient(mp, *pde->solution());
            gsWriteParaview<>( out2, "field_output_grad", samples);
        }
    }
    else
    {
        gsInfo<<"Outputing normal field to paraview file...\n";
        gsField<> nfield = gsFieldCreator<>::normal(mp);
        //gsField<> nfield = gsFieldCreator<>::jacDet(mp);
        gsWriteParaview<>(nfield, "field_output", samples);
    }
    

    if (plot)
    {
        // Run paraview
        if ( mp.nPatches() == 1 )
            return system("paraview field_output.vts &");
        else
            return system("paraview field_output.pvd &");
    }
    else
        return 0;
}
