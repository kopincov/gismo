#include <iostream>

#include <gismo.h>
#include <gsIO/gsWriteLsdyna.h>



using namespace gismo;

int main(int argc, char *argv[])
{
    std::string input("planar/thbsLSDYNA.xml"); 
    std::string output(""); 
    index_t n(0);
    index_t m(0);
    
    gsCmdLine cmd("Outputs LS-DYNA file.");
    cmd.addPlainString("filename", "Filename containing a geometry.", input);
    cmd.addString("o", "output", "Output file (LSDYNA .k)", output);
    cmd.addInt("n", "dof1", "Grid size in u direction.", n);
    cmd.addInt("m", "dof2", "Grid size in v direction.", m);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "$ -------------------------------------\n"
              << "$ Arguments: \n"
              << "$ input: " << input << "\n"
              << "$ output: " << output << "\n"
              << "$ n: " << n << "\n"
              << "$ m: " << m << "\n"
              << "$ -------------------------------------\n";

    gsGeometry<>::uPtr surf = gsReadFile<>(input);
    if (!surf)
    {
        gsInfo << "No input\n";
        return 0;
    }

    
    gsLsdynaFEM<real_t> writer(*surf);
    gsVector<index_t, 2> gridDimension;
    gridDimension << n, m;

    if (!output.empty())
    {
        std::ofstream file(output.c_str());

        writer.writeKFile(file, "gismo", gridDimension);
        
        file.close();
    }
    else
    {
        writer.writeKFile(gsInfo, "gismo", gridDimension);
    }

    return 0;
}


    

        
        
