#include <iostream>

#include <gismo.h>
#include <gsIO/gsWriteLsdyna.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    std::string input("planar/two_squares.xml");
    std::string output(""); 
    
    gsCmdLine cmd("Outputs LS-DYNA file.");
    cmd.addPlainString("filename", "Filename containing a geometry.", input);
    cmd.addString("o", "output", "Output LS-DYNA (.k) file", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (input.empty())
    {
        gsInfo << "Please provide an input filename or use -h for help.\n";
        return 0;
    }
    
    gsInfo << "$ -------------------------------------\n"
              << "$ Arguments: \n"
              << "$ input: " << input << "\n"
              << "$ output: " << output << "\n"
              << "$ -------------------------------------\n";

    gsGeometry<>::uPtr surf = gsReadFile<>(input);
    if (!surf)
    {
        gsInfo << "No input\n";
        return 0;
    }

    
    gsLsdynaIGA<real_t> writer(*surf);

    if (output != "")
    {
        std::ofstream file(output.c_str());

        writer.writeKFile(file, "gismo");
        
        file.close();
    }
    else
    {
        writer.writeKFile(gsInfo, "gismo");
    }

    return 0;
}
