/** @file gsScrewMachineGeometry.cpp

    @brief Reads the data for the compressible fluid flow solver
    (coefficients, domain, boundary conditions) from an XML file and
    solves the (initial) boundary value problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Hinz, M. Moeller
*/

#include <gismo.h>

using namespace gismo;

template<typename T>
class gsScrewMachine
{
public:

    /// \brief Constructor (default)
    gsScrewMachine()
    : fn(""),
      fn_points(""),
      fn_segments("")
    {}

    /// \brief Constructor (from command line arguments)
    gsScrewMachine(int argc, char *argv[]) : fn(""), fn_points(""), fn_segments("")
    {
        // Define list of command line arguments.
        gsCmdLine cmd("gsScrewMachine mesh generator.");

        cmd.addPlainString("filename",
                           "File containing the complete configuration (.xml)",
                           fn);

        cmd.addString("p", "fn_points",
                      "File containing the point cloud of the profile (.xml)",
                      fn_points);

        cmd.addString("s", "fn_segments",
                      "File containing the segmentation of the profile (.xml)",
                      fn_segments);

        // Get command line arguments
        cmd.getValues(argc,argv);

        if (fn.empty())
        {
            gsInfo << cmd.getMessage();
            gsInfo
                << "\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
            exit(EXIT_FAILURE);
        }

        // Open xmlData
        gsFileData<T> xmlData(fn);

        // Read point cloud of the profile
        gsDebug << "Reading point cloud from file...\n";

        if (!fn_points.empty())
        {
            gsFileData<T> xmlData(fn_points);
            rotor  = *xmlData.template getId<gsMatrix<T> >(0);
            casing = *xmlData.template getId<gsMatrix<T> >(1);
        }
        else
        {
            rotor  = *xmlData.template getId<gsMatrix<T> >(0);
            casing = *xmlData.template getId<gsMatrix<T> >(1);
        }

        rotor = rotor.transpose();
        casing = casing.transpose();

        gsInfo << rotor.rows() << "," << rotor.cols();

        // Read segmentation points of the profile
        gsDebug << "Reading segmentation points from file...\n";

        if (!fn_segments.empty())
        {
            gsFileData<T> xmlData(fn_segments);
            segments = *xmlData.template getId<gsMatrix<T> >(0);
        }
        else
        {
            segments = *xmlData.template getId<gsMatrix<T> >(2);
        }

        for (index_t row = 1; row < rotor.rows(); row++)
            {
                if ((rotor(row-1, 0) == rotor(row, 0)) &&
                    (rotor(row-1, 1) == rotor(row, 1)))
                    gsInfo << rotor(row,0) << "," << rotor(row,1) << std::endl;
            }

    }

private:
    /// Attributes
    std::string fn;
    std::string fn_points;
    std::string fn_segments;

    /// Rotor + casing (Point clouds)
    gsMatrix<T> rotor;
    gsMatrix<T> casing;
    gsMatrix<T> segments;

    /*
    //this function adds the first point of the matrix to its end
    gsMatrix<>  TuDoScrewMachine::closeContour(gsMatrix<> unclosedContour)
    {

        gsMatrix<> matrixClosed(unclosedContour.rows() + 1, unclosedContour.cols() + 1);
        for (int c = 0; c < unclosedContour.rows(); c++)
        {
            matrixClosed(c, 0) = unclosedContour(c, 0);
            matrixClosed(c, 1) = unclosedContour(c, 1);
            matrixClosed(c, 2) = 0;
        }
        matrixClosed(unclosedContour.rows(), 0) = unclosedContour(0, 0);
        matrixClosed(unclosedContour.rows(), 1) = unclosedContour(0, 1);
        matrixClosed(unclosedContour.rows(), 2) = 0;

        return matrixClosed;
        }*/

    /*
    // Rotate given contour around axialDistance with angle
    gsMatrix<T> rotateContour(gsMatrix<T> contour, T angle, T axialDistance)
    {

        gsMatrix<T> coefs(contour.rows(), 2);
        for (index_t row = 0; row < contour.rows(); row++)
        {
            coefs(row, 0) = (cos(angle / 180.0 * EIGEN_PI)*(contour(row, 0) - axialDistance) - sin(angle / 180.0 * EIGEN_PI)*contour(row, 1)) + axialDistance;
            coefs(row, 1) = sin(angle / 180.0 * EIGEN_PI)*(contour(row, 0) - axialDistance) + cos(angle / 180.0 * EIGEN_PI)*contour(row, 1);

        }
        //crossSections.push_back(matrix);
        return coefs;
        }*/

};

int main(int argc, char *argv[])
{
    try { gsScrewMachine<real_t> geo(argc, argv); } catch (int& e) { return e; }
}
