// Author: Jaka Speh
//
// Program takes an xml G+SMO file and outputs the .g2 GoTools format.
//
// Currently supported G+SMO files:
// - gsTensorBSpline for dimension 1, 2, 3


#include <iostream>

#include <gismo.h>
#include <gsIO/gsIOUtils.h>


using namespace gismo;


int main(int argc, char* argv[])
{

    std::string inFile = "volumes/cube.xml";
    std::string outFile("goTools");
    gsCmdLine cmd("Hi, give me a xml file and I will write it "
                  "in .g2 format");
    cmd.addPlainString("filename", "G+SMO XML file", inFile);
    cmd.addString("o", "output", "Output file", outFile);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    gsInfo << "------------------------------------------------------------"
                 "\nInput Arguments: \n\n"
                 "Input filename: " << inFile << "\n\n"
                 "Ouput prefix: " << outFile << "\n"
                 "------------------------------------------------------------"
                 "\n\n";

    gsFileData<> data(inFile);

    if (data.has< gsMultiPatch<> >())
    {
        gsMultiPatch<> mp;
        data.getFirst< gsMultiPatch<> >(mp);
        gsWriteGoTools(mp, outFile);
    }
    else if (data.has< gsGeometry<> >())
    {
        gsGeometry<>::uPtr geom = data.getFirst< gsGeometry<> >();
        gsWriteGoTools(*geom, outFile);
    }

    return 0;
}
