#include <gismo.h>
#include <gismo_dev.h>


#include <iostream>
#include <string>

using namespace gismo;



int main(int argc, char* argv[])
{
    std::string input("fitting/discontinuous_cylinder.xml");
    std::string output("out");
    
    gsCmdLine cmd("Plotting points.");
    cmd.addString("i", "input", "Input file", input);
    cmd.addString("o", "output", "Output file", output);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsInfo << " \n\nInput arguments: \n\n"
	      << "input: " << input << "\n\n"
	      << "output: " << output << "\n\n"
	      << "--------------------------------------------------\n" << "\n";

    gsFileData<> fd_in(input);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);

    std::string out = output + "_params";
    gsInfo << "Writing: " << out << "\n";
    gsWriteParaviewPoints(uv, out);
    
    out = output + "_points";
    gsInfo << "Writing: " << out << "\n";
    gsWriteParaviewPoints(xyz, out);
    
    return 0;
}
