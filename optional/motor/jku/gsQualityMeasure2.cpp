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

void saveData(const gsGeometry<>& geom,
              const std::string& output,
              const int number)
{
    gsFileData<> fileData;
    fileData << geom;

    std::string out = output + "_Map_" + util::to_string(number) + ".xml";
    fileData.dump(out);

    gsMesh<> mesh;
    makeMesh(geom.basis(), mesh, 10);
    geom.evaluateMesh(mesh);
    out = output + "_Mesh_" + util::to_string(number);
    gsWriteParaview(mesh, out);

//    gsMatrix<TT> coefs = geom.coefs();
//    coefs.transposeInPlace();
//    out = output + "ControlPoints" + util::to_string(number);
//    gsWriteParaviewPoints(coefs, out);

//    gsMatrix<unsigned> boundary = geom.basis().boundary();
//    gsMatrix<> b(geom.geoDim(), boundary.rows());
//    for (int row = 0; row < boundary.rows(); row++)
//    {
//        b.col(row) = coefs.col((boundary)(row, 0));
//    }

//    out = output + "BoundaryControlPoints" + util::to_string(number);
//    gsWriteParaviewPoints(b, out);

}

int main(int argc, char* argv[])
{

    //std::string input(GISMO_DDATA_DIR "/planar/deformedSquare.xml");
    std::string input(MOTOR_DATA_DIR "/jku/thb_map.xml");
    std::string inputPtsParams(MOTOR_DATA_DIR "/jku/thb_parameters_and_points.xml");
    std::string output("results_optimized");

    int iterations = 0;
    real_t fitting = 0;
    real_t orthogonality = 0;
    real_t skewness = 0;
    real_t eccentricity = 0;
    real_t intersection = 1e-2;
    real_t uniformity = 0;
    real_t area = 0;
    real_t length = 0;
    real_t epsilon = 1e-7;
    int jacPts = 100000;
    bool dumped = false;

    gsCmdLine cmd("Optimization");
    cmd.addString("I", "input", "Input prefix", input);
    cmd.addString("P", "pointsparameters", "Name of file with parameters and points", inputPtsParams);
    cmd.addString("f", "output", "Output prefix", output);
    cmd.addInt("i", "iterations", "Number of iterations", iterations);
    cmd.addReal("F", "fitting",
                "Weight of quality measure: fitting", fitting);
    cmd.addReal("o", "orthogonality",
                "Weight of quality measure: orthogonality", orthogonality);
    cmd.addReal("s", "skewness",
                "Weight of quality measure: skewness", skewness);
    cmd.addReal("e", "eccentricity",
                "Weight of quality measure: eccentricity", eccentricity);
    cmd.addReal("u", "uniformity",
                "Weight of quality measure: uniformity", uniformity);
    cmd.addReal("L", "length",
                "Weight of quality measure: length functional", length);
    cmd.addReal("a", "area",
                "Weight of quality measure: area", area);
    cmd.addReal("n", "intersection",
                "Weight of quality measure: self-intersection", intersection);
    cmd.addReal("p", "epsilon", "Self intersection variable.", epsilon);
    cmd.addInt("j", "numJacPts",
               "Number of points where we sample Jacobian determinant",
               jacPts);
    cmd.addSwitch("dumped", "Use dumped method.", dumped);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "------------------------------------------------------------"
                 "\nInput Arguments: \n"
                 "Input: " << input << "\n"
                 "Ouput prefix: " << output << "\n"
                 "Iterations: " << iterations << "\n"
                 "Fitting: " << fitting << "\n"
                 "Orthogonality: " << orthogonality << "\n"
                 "Skewness: " << skewness << "\n"
                 "Eccentricity: " << eccentricity << "\n"
                 "Uniformity: " << uniformity << "\n"
                 "Length: " << length << "\n"
                 "Area: " << area << "\n"
                 "Intersection: " << intersection << "\n"
                 "Epsilon: " << epsilon << "\n"
                 "Jacobian Points: " << jacPts << "\n"
                 "Dumped: " << dumped << "\n"
                 "------------------------------------------------------------"
                 "\n\n";

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

    gsFileData<> fd_in(inputPtsParams);
    gsMatrix<> uv     = *fd_in.getId<gsMatrix<> >(0);
    gsMatrix<> xy     = *fd_in.getId<gsMatrix<> >(1);

    gsQualityMeasure2<real_t> optimization(*geom, uv, xy, false);

    gsInfo << "Value of functional: "
              << optimization.functional(fitting, orthogonality, skewness,
                                         eccentricity, uniformity,
                                         length, area,
                                         intersection, epsilon)
              << "\n";

    saveData(*geom, output, 0);
    checkJacobianDeterminant(*geom, jacPts, true, output, 0);


    for (int it = 0; it != iterations; it++)
    {
        gsInfo << "iteration: " << it << " / " << iterations - 1 << "\n";

        optimization.optimize(fitting, orthogonality, skewness,
                              eccentricity, uniformity,
                              length, area,
                              intersection, epsilon, dumped);

        gsInfo << "Value of functional: "
                  << optimization.functional(fitting, orthogonality, skewness,
                                             eccentricity, uniformity,
                                             length, area,
                                             intersection, epsilon)
                  << "\n";

        saveData(*geom, output, it + 1);
        checkJacobianDeterminant(*geom, jacPts, true, output, it + 1);

        // Giving different names for files with (possibly) final results
        saveData(*geom, output + "_final", 0);
    }

    delete geom;
    return 0;
}
