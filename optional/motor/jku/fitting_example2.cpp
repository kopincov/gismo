///** @file fitting_example2.cpp

//    @brief Demonstrates fitting of data samples

//    This file is part of the G+Smo library.

//    This Source Code Form is subject to the terms of the Mozilla Public
//    License, v. 2.0. If a copy of the MPL was not distributed with this
//    file, You can obtain one at http://mozilla.org/MPL/2.0/.

//    Author(s): ...
//*/


//#include <iostream>
//#include <gismo.h>
//#include <stdio.h>
//#include <stdlib.h>

//using namespace gismo;

//// forward declaration of some utility functions
//void print(const gsBSplineBasis<>& bsb, const std::string& name);
//gsMatrix<> points(const real_t x_1, const real_t y_1,
//                  const real_t x_2, const real_t y_2,
//                  const int numPoints);

int main()
{
    return 0;
}

//int main(int argc, char *argv[])
//{

//    // Options with default values
//    bool save     = false;
//    int numURef   = 0;
//    int iter      = 2;
//    index_t k     = 3;
//    index_t deg    = 3;
//    real_t lambda = 1e-07;
//    real_t threshold = 1e-02;
//    real_t tolerance   = 1e-02;
//    int extension = 2;
//    real_t refPercent  = 0.1;
//    index_t numSamples = 19; // Number of data samples

//    index_t numSeg = 9;  // Number of segments
//    std::string it_num;

//    // Reading options from the command line
//    gsCmdLine cmd("Example of data fitting using B-splines");
//    cmd.addSwitch("save", "Save result in XML format", save);
//    cmd.addInt("i", "iter", "number of iterations", iter);
//    cmd.addInt("k", "knots", "this is the number of interior knots", k);
//    cmd.addInt("d", "degree", "this is the degree", deg);
//    cmd.addReal("s", "lambda", "smoothing coefficient", lambda);
//    cmd.addReal("t", "threshold", "error threshold (special valule -1)", threshold);
//    cmd.addReal("p", "refPercent", "percentage of points to refine in each iteration", refPercent);
//    cmd.addInt("q", "extension", "extension size", extension);
//    cmd.addInt("r", "urefine", "initial uniform refinement steps", numURef);
//    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tolerance);
//    cmd.addInt("n", "samples", "This is the number of samples", numSamples);

//    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

//    if (deg < 1)
//    { gsInfo << "Degree x must be positive.\n";  return 1;}

//    if ( tolerance < 0 )
//    {
//        gsInfo << "Error tolerance cannot be negative, setting it to default value.\n";
//        tolerance = 1e-02;
//    }

//    if (threshold > 0 && threshold > tolerance )
//    {
//        gsInfo << "Refinement threshold is over tolerance, setting it the same as tolerance.\n";
//        threshold = tolerance;
//    }

//    if ((numSamples-1) % numSeg != 0)
//    {
//        numSamples -= ((numSamples-1) % numSeg);
//        gsWarn << "Number of samples is changed (decreased): " << numSamples << "\n";
//    }

//    gsMatrix<real_t> m_points(3, numSamples);
//    m_points.setZero();

//    gsMatrix<real_t> m_param_values(1, numSamples);
//    m_param_values.setZero();

//    gsInfo << "Generating data points ...\n";
//    real_t dist  = 1.0*numSeg/(numSamples-1);
//    index_t numPointsSeg = (numSamples-1)/numSeg; // The last point is not included

//    m_points.block(0, 0*numPointsSeg, 3, numPointsSeg) = points(0.0, 0.0, 1.0-dist, 0.0, numPointsSeg);
//    m_points.block(0, 1*numPointsSeg, 3, numPointsSeg) = points(1.0, 0.0, 1.0, 1.0-dist, numPointsSeg);
//    m_points.block(0, 2*numPointsSeg, 3, numPointsSeg) = points(1.0, 1.0, 2.0-dist, 1.0, numPointsSeg);
//    m_points.block(0, 3*numPointsSeg, 3, numPointsSeg) = points(2.0, 1.0, 2.0, 0.0+dist, numPointsSeg);
//    m_points.block(0, 4*numPointsSeg, 3, numPointsSeg) = points(2.0, 0.0, 3.0-dist, 0.0, numPointsSeg);
//    m_points.block(0, 5*numPointsSeg, 3, numPointsSeg) = points(3.0, 0.0, 3.0, 1.0-dist, numPointsSeg);
//    m_points.block(0, 6*numPointsSeg, 3, numPointsSeg) = points(3.0, 1.0, 4.0-dist, 1.0, numPointsSeg);
//    m_points.block(0, 7*numPointsSeg, 3, numPointsSeg) = points(4.0, 1.0, 4.0, 0.0+dist, numPointsSeg);
//    m_points.block(0, 8*numPointsSeg, 3, numPointsSeg+1) = points(4.0, 0.0, 5.0, 0.0,      numPointsSeg+1);

///*  Silly way ****************
//    gsInfo << "Generating data points ...\n";

//    const real_t eps = 1e-14;
//    real_t dist  = 1.0*numSeg/(s-1);
//    real_t current_x = 0.0;
//    real_t current_y = 0.0;
//    index_t i = 0;

//    // Segment 1
//    do
//    {
//        m_points(0, i) = current_x;
//        m_points(1, i) = current_y;
//        gsInfo << m_points(0, i) << " " << m_points(1, i) << "\n";
//        current_x += dist;
//        i++;
//    } while (current_x < 1.0 - eps);
//    gsInfo << "Done with Segment 1\n";

//    // Segment 2
//    do
//    {
//        m_points(0, i) = current_x;
//        m_points(1, i) = current_y;
//        gsInfo << m_points(0, i) << " " << m_points(1, i) << "\n";
//        current_y += dist;
//        i++;
//    } while (current_y < 1.0 - eps);
//    gsInfo << "Done with Segment 2\n";

//    // Segment 3
//    do
//    {
//        m_points(0, i) = current_x;
//        m_points(1, i) = current_y;
//        gsInfo << m_points(0, i) << " " << m_points(1, i) << "\n";
//        current_x += dist;
//        i++;
//    } while (current_x < 2.0 - eps);
//    gsInfo << "Done with Segment 3\n";

//    // Segment 4
//    do
//    {
//        m_points(0, i) = current_x;
//        m_points(1, i) = current_y;
//        gsInfo << m_points(0, i) << " " << m_points(1, i) << "\n";
//        current_y -= dist;
//        i++;
//    } while (current_y > 0.0 + eps);
//    gsInfo << "Done with Segment 4\n";

//    // Segment 5
//    do
//    {
//        m_points(0, i) = current_x;
//        m_points(1, i) = current_y;
//        gsInfo << m_points(0, i) << " " << m_points(1, i) << "\n";
//        current_x += dist;
//        i++;
//    } while (current_x < 3.0 - eps);
//    gsInfo << "Done with Segment 5\n";

//    // Segment 6
//    do
//    {
//        m_points(0, i) = current_x;
//        m_points(1, i) = current_y;
//        gsInfo << m_points(0, i) << " " << m_points(1, i) << "\n";
//        current_y += dist;
//        i++;
//    } while (current_y < 1.0 - eps);
//    gsInfo << "Done with Segment 6\n";

//    // Segment 7
//    do
//    {
//        m_points(0, i) = current_x;
//        m_points(1, i) = current_y;
//        gsInfo << m_points(0, i) << " " << m_points(1, i) << "\n";
//        current_x += dist;
//        i++;
//    } while (current_x < 4.0 - eps);
//    gsInfo << "Done with Segment 7\n";

//    // Segment 8
//    do
//    {
//        m_points(0, i) = current_x;
//        m_points(1, i) = current_y;
//        gsInfo << m_points(0, i) << " " << m_points(1, i) << "\n";
//        current_y -= dist;
//        i++;
//    } while (current_y > 0.0 + eps);
//    gsInfo << "Done with Segment 8\n";

//    // Segment 9
//    do
//    {
//        m_points(0, i) = current_x;
//        m_points(1, i) = current_y;
//        gsInfo << m_points(0, i) << " " << m_points(1, i) << "\n";
//        current_x += dist;
//        i++;
//    } while (current_x <= 5.0 + eps);
//    gsInfo << "Done with Segment 9\n";
//*/

////    gsInfo << "Points (x):" << "\n" << m_points.row(0) << "\n";
////    gsInfo << "Points (y):" << "\n" << m_points.row(1) << "\n";

//    gsInfo << "Generating parameters ...\n";
//    for (int i = 0; i < numSamples; i++)
//    {
//        m_param_values(0, i) = 1.0*i/(numSamples-1);
//    }
////    gsInfo << "Param. values:" << "\n" << m_param_values << "\n";

//    // This is for outputing an XML file, if requested
//    gsFileData<> fd;

//    gsInfo << "Generating the basis ...\n";
//    gsKnotVector<> kv (0, 1, k, deg+1 );

//    gsBSplineBasis<> T_tbasis( kv );
//    T_tbasis.uniformRefine( (1<<numURef)-1 );
//    //print(bsb, "bsb");

//    // Create Initial hierarchical basis

//    gsTHBSplineBasis<1, real_t> HB( T_tbasis );
//    //print(hbsb, "hbsb");

//    // Specify extension size
//    std::vector<int> ext;
//    ext.push_back(extension);

//    gsHFitting<1, real_t> ref (m_param_values, m_points, HB, refPercent, ext, lambda);
//    //ref.compute(lambda);

//    const std::vector<real_t> & errors = ref.pointWiseErrors();

//    // Print settings summary
//    gsInfo << "Fitting "<< m_points.cols() <<" samples.\n";
//    gsInfo << "----------------\n";
//    gsInfo << "Cell extension     : " << ext[0] << " " << ext[1] << ".\n";
//    if ( threshold >= 0.0 )
//        gsInfo << "Ref. threshold     : " << threshold << ".\n";
//    else
//        gsInfo << "Cell refinement    : " << 100*refPercent << "%%.\n";
//    gsInfo << "Error tolerance    : " << tolerance << ".\n";
//    gsInfo << "Smoothing parameter: " << lambda << ".\n";

//    gsStopwatch time;

//    for(int i = 0; i <= iter; i++)
//    {
//        gsInfo << "----------------\n";
//        gsInfo << "Iteration " << i << ".." << "\n";

//        time.restart();
//        ref.nextIteration(tolerance, threshold);
//        const double clock = time.stop();
//        gsInfo << "Fitting time: " << clock << "\n";
////        ref.computeErrors();

//        gsInfo << "Fitted with " << ref.result()->basis() << "\n";
//        gsInfo << "Min distance : " << ref.minPointError() << " / ";
//        gsInfo << "Max distance : " << ref.maxPointError() << "\n";
//        gsInfo << "Points below tolerance: " << 100.0 * ref.numPointsBelow(tolerance)/errors.size() << "%.\n";

//        if ( save )
//        {
//            std::stringstream ss;
//            ss << "fitting_out2_" << i;
//            gsWriteParaview(*ref.result(), ss.str());
//            gsInfo << "Writing fitting_out2_" << i << ".xml\n";
//        }

//        if ( ref.maxPointError() < tolerance )
//        {
//            gsInfo << "Error tolerance achieved after " << i << " iterations.\n";
//            break;
//        }

//    }
///*
//    gsInfo << "----------------\nFinished.\n";

//    if ( save )
//    {
//        gsInfo << "Writing fitting_data2.vtp\n";
//        gsWriteParaviewPoints(m_points, "fitting_data2");
//    }
//*/
//    return 0;
//}

//void print(const gsBSplineBasis<>& bsb,
//           const std::string& name)
//{
//    gsInfo << name << ": \n";
//    bsb.print(gsInfo);
//    gsInfo << "\n\n";
//}

//gsMatrix<> points(const real_t x_1, const real_t y_1,
//                  const real_t x_2, const real_t y_2,
//                  const int numPoints)
//{
//    //GISMO_ENSURE((x_1 <= x_2) && (y_1 <= y_2), "Incorrect coordinates of points");
//    gsMatrix<> m(3, numPoints);
//    m.setZero();
//    // Distance between two points:
//    real_t dist = math::sqrt((x_2-x_1)*(x_2-x_1) + (y_2-y_1)*(y_2-y_1))/(numPoints-1);
//    if (!(x_1 <= x_2) || !(y_1 <= y_2))
//    {
//        dist = -dist; // Important step
//    }
//    real_t coord_x = x_1;
//    real_t coord_y = y_1;

//    if (y_1 == y_2)
//    {
//        for (int i = 0; i < numPoints; i++)
//        {
//            m(0, i) = coord_x;
//            m(1, i) = coord_y;
//            coord_x += dist;
//        }
//    }
//    else if (x_1 == x_2)
//    {
//        for (int i = 0; i < numPoints; i++)
//        {
//            m(0, i) = coord_x;
//            m(1, i) = coord_y;
//            coord_y += dist;
//        }
//    }
//    else
//    {
//        GISMO_ERROR("Points should lie in horizontal or vertical line");
//    }
//    return m;
//}

