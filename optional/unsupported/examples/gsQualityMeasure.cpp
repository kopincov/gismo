#ifdef _MSC_VER // to be removed
  #define _USE_MATH_DEFINES
#endif


#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>
#include <gsIO/gsIOUtils.h>
#include <gsOptimizer/gsQualityMeasure.h>


using namespace gismo;

void saveData(const gsGeometry<>& geom,
              const std::string& output,
              const int number)
{
    gsFileData<> fileData;
    fileData << geom;

    std::string out = output + util::to_string(number) + ".xml";
    fileData.dump(out);

    gsMesh<> mesh;
    makeMesh(geom.basis(), mesh, 10);
    geom.evaluateMesh(mesh);
    out = output + "Mesh" + util::to_string(number);
    gsWriteParaview(mesh, out);

//    gsMatrix<TT> coefs = geom.coefs();
//    coefs.transposeInPlace();
//    out = output + "ControlPoints" + util::to_string(number);
//    gsWriteParaviewPoints(coefs, out);

//    gsMatrix<index_t> boundary = geom.basis().boundary();
//    gsMatrix<> b(geom.geoDim(), boundary.rows());
//    for (int row = 0; row < boundary.rows(); row++)
//    {
//        b.col(row) = coefs.col((boundary)(row, 0));
//    }

//    out = output + "BoundaryControlPoints" + util::to_string(number);
//    gsWriteParaviewPoints(b, out);

}


void checkJacobianDeterminant(const gsGeometry<>& geom,
                              const int points = 100,
                              const bool savePoints = false,
                              const std::string& output = "",
                              const int number = 0)
{
    gsInfo << "Checking Jacobian determinant ..." << "\n";

    gsMatrix<> para  = geom.support();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> pts = uniformPointGrid(c0, c1, points);

    gsMatrix<> plus(pts.rows(), pts.cols());
    plus.setZero();
    gsMatrix<> minus(pts.rows(), pts.cols());
    minus.setZero();
    gsMatrix<> zero(pts.rows(), pts.cols());
    zero.setZero();

    int plusCounter = 0;
    int minusCounter = 0;
    int zeroCounter = 0;

    // calculating the determinant of the Jacobian
    for (int col = 0; col < pts.cols(); col++)
    {
        const real_t determinant = geom.jacobian(pts.col(col)).determinant();
        gsMatrix<> tmp = geom.eval(pts.col(col));

        if (determinant < 0)
        {
            minus.col(minusCounter) = tmp.col(0);
            minusCounter++;
        }
        else if (determinant > 0)
        {
            plus.col(plusCounter) = tmp.col(0);
            plusCounter++;
        }
        else
        {
            zero.col(zeroCounter) = tmp.col(0);
            zeroCounter++;
        }
    }

    gsInfo << "Number of points with \n"
              << "  - positive sign: " << plusCounter << "\n"
              << "  - negative sign: " << minusCounter << "\n"
              << "  - zero sign: " << zeroCounter << "\n" << "\n";

    if (savePoints)
    {
        if (0 < plusCounter)
        {
            std::string out = output + "PositivePoints" + util::to_string(number);
            gsWriteParaviewPoints(plus, out);
        }

        if (0 < minusCounter)
        {
            std::string out = output + "NegativePoints" + util::to_string(number);
            gsWriteParaviewPoints(minus, out);
        }

        if (0 < zeroCounter)
        {
            std::string out = output + "ZeroPoints" + util::to_string(number);
            gsWriteParaviewPoints(zero, out);
        }
    }
}


int main(int argc, char* argv[])
{
    std::string input("planar/deformedSquare.xml");
    std::string output("optimization");
    index_t iterations = 0;
    real_t orthogonality = 0;
    real_t skewness = 0;
    real_t eccentricity = 0;
    real_t intersection = 0;
    real_t uniformity = 0;
    real_t area = 0;
    real_t length = 0;
    real_t epsilon = 1e-7;
    index_t jacPts = 10000;
    bool dumped = false;
    
    gsCmdLine cmd("Optimization");
    cmd.addString("I", "input", "Input", input);
    cmd.addString("f", "output", "Output prefix", output);
    cmd.addInt("i", "iterations", "Number of iterations", iterations);
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
                 "\nInput Arguments: \n\n"
                 "Input: " << input << "\n\n"
                 "Ouput prefix: " << output << "\n\n"
                 "Iterations: " << iterations << "\n\n"
                 "Orthogonality: " << orthogonality << "\n\n"
                 "Skewness: " << skewness << "\n\n"
                 "Eccentricity: " << eccentricity << "\n\n"
                 "Uniformity: " << uniformity << "\n\n"
                 "Length: " << length << "\n\n"
                 "Area: " << area << "\n\n"
                 "Intersection: " << intersection << "\n\n"
                 "Epsilon: " << epsilon << "\n\n"
                 "Jacobian Points: " << jacPts << "\n\n"
                 "Dumped: " << dumped << "\n\n"
                 "------------------------------------------------------------"
                 "\n\n";


    gsFileData<> data(input);
    gsGeometry<>::uPtr geom;
    if (data.has< gsGeometry<> >())
    {
        geom = data.getFirst< gsGeometry<> >();
    }

    if (!geom)
    {
        gsInfo << "Didn't get the input geometry. Aborting..." << "\n";
    }

    gsQualityMeasure<> optimization(*geom);

    gsInfo << "Value of functional: "
              << optimization.functional(orthogonality, skewness,
                                         eccentricity, uniformity,
                                         length, area,
                                         intersection, epsilon)
              << "\n";

    saveData(*geom, output, 0);
    checkJacobianDeterminant(*geom, jacPts, true, output, 0);


    for (int it = 0; it != iterations; it++)
    {
        gsInfo << "iteration: " << it << " / " << iterations - 1 << "\n";

        optimization.optimize(orthogonality, skewness, eccentricity, uniformity,
                              length, area,
                              intersection, epsilon, dumped);

        gsInfo << "Value of functional: "
                  << optimization.functional(orthogonality, skewness,
                                             eccentricity, uniformity, length,
                                             area, intersection, epsilon)
                  << "\n";

        saveData(*geom, output, it + 1);
        checkJacobianDeterminant(*geom, jacPts, true, output, it + 1);
    }

    return 0;
}
