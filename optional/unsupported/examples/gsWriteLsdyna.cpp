
#include <iostream>

#include <gismo.h>
#include <gsIO/gsWriteLsdyna.h>



using namespace gismo;


int main(int argc, char *argv[])
{
    std::string input("planar/thbsLSDYNA.xml"); 
    std::string output(""); 
    bool solid(false);
    index_t degree(2);
    
    gsCmdLine cmd("Outputs LS-DYNA file.");
    cmd.addPlainString("filename", "Filename containing a geometry.", input);
    cmd.addString("o", "output", "Output LS-DYNA (.k) file", output);
    cmd.addSwitch("solid", "Treat shells as 3D (solid) objects.", solid);
    cmd.addInt("d", "degree", "Degree of through thickness basis.", degree);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // gsInfo << "-------------------------------------\n"
    //           << "Arguments: "
    //           << "input: " << input << "\n\n"
    //           << "output: " << output << "\n\n"
    //           << "-------------------------------------\n";

    gsGeometry<>::uPtr surf = gsReadFile<>(input);
    if (!surf)
    {
        gsInfo << "No input\n";
        return 0;
    }

    gsWriteLsdyna<real_t> writer(*surf);

    if (output == "")
    {
        writer.produceLsdynaInputFile(gsInfo, "out", solid, degree);
    }
    else
    {
        std::ofstream file;
        file.open(output.c_str());
        if (file.is_open())
        {
            writer.produceLsdynaInputFile(file, output, solid, degree);
            file.close();
        }
        else
        {
            gsInfo << "Can not open output file..." << "\n";
        }
    }

    return 0;
}



