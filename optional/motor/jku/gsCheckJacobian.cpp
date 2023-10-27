#ifdef _MSC_VER // to be removed
  #define _USE_MATH_DEFINES
#endif


#include <iostream>

#include <gismo.h>
#include <gsIO/gsIOUtils.h>
//#include <gsOptimizer/gsQualityMeasure.h>
#include "gsQualityMeasure2.h"
#include "gsMotorUtils.h"


using namespace gismo;

int main(int argc, char* argv[])
{

    std::string input(MOTOR_DATA_DIR "/jku/thb_map.xml");
    std::string output("points");

    int jacPts = 10000;

    gsCmdLine cmd("Checking Jacobian determinant sign at many different points");
    cmd.addString("I", "input", "Input prefix", input);
    cmd.addInt("j", "numJacPts",
               "Number of points where we sample Jacobian determinant",
               jacPts);
    cmd.addString("o", "output", "Name of output file", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "------------------------------------------------------------"
                 "\nInput Arguments: \n"
                 "Input: " << input << "\n"
                 "Jacobian Points: " << jacPts << "\n"
                 "Output file: " << output << "\n"
                 "------------------------------------------------------------"
                 "\n\n";

    /*
    if (input == "")
    {
        gsInfo << "provide input geometry" << "\n";
    }
    */

    gsFileData<> data(input);
    gsGeometry<>* geom = NULL;
    if (data.has< gsGeometry<> >())
    {
        geom = data.getFirst< gsGeometry<> >().release();
    }

    if (geom == NULL)
    {
        gsInfo << "Didn't get the input geometry. Aborting..." << "\n";
    }

    checkJacobianDeterminant(*geom, jacPts, true, output, 0);


    delete geom;
    return 0;
}
