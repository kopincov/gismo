/* Tests functions for creating a periodic B-spline basis.
 *
 * Author: Dominik Mokris, dominik.mokris@jku.at
 *
 * First, some values of the periodic basis are automatically compared with the corresponding yielded by the non-periodic counterpart.
 * Second, several files are prepared to visually check and compare in Paraview.
*/

#include <iostream>
#include <string>
#include <vector>

#include <gismo.h>
#include <gismo_dev.h>


using namespace gismo;


// Error message, when periodic and non-periodic bases give different results where they should agree.
void printError( gsMatrix<>& result, gsMatrix<>& presult, int& test_result, std::string message = "" )
{
    gsWarn << "Periodic and nonperiodic basis yield different " << message << ".\n";
    gsInfo << "Normal basis result:" << "\n" << result << "\n";
    gsInfo << "Periodic basis result:" << "\n" << presult << "\n";

    test_result ++;
}

// Checks that the number of basis functions and values of the basis functions, first and second derivatives are the same
//  in some points.
void automatedTesting( int& test_result )
{
    // We prepare a basis and a second one that will be converted to periodic.
    gsBSplineBasis<>::uPtr basis;
    gsBSplineBasis<>::uPtr pbasis;
    std::string filename = "hat_closedcurve.xml";
    gsFileData<> data(filename);
    basis = data.getAnyFirst< gsBSplineBasis<> >();
    pbasis = data.getAnyFirst< gsBSplineBasis<> >();

    int bsize = basis -> size();
    // Now convert to periodic and check whether the sizes agree.
    pbasis -> setPeriodic(true);

    int pm_periodic = pbasis -> numCrossingFunctions();
    int psize = pbasis->size();

    // If the sizes disagree, the test has failed.
    if( bsize!= 84 || bsize != psize + pm_periodic || pm_periodic == 0 )
    {
        gsWarn << "Wrong number of basis functions.\n";
        test_result ++;
    }
    gsInfo << "basis size: " << bsize << "\n";
    gsInfo << "pbasis size: " << psize << "\n";
    gsInfo << "pm_periodic: " << pm_periodic << "\n";

    gsInfo << "pbasis twin of "<< 1 <<": " << pbasis->twin(1) << "\n";
    gsInfo << "pbasis twin of "<< pbasis->twin(1) <<": " <<
        pbasis->twin(pbasis->twin(1)) << "\n";


    // Now we compare evaluation procedures
    gsMatrix<> eval_pts(1,4);
    eval_pts(0,0) = -0.04;
    eval_pts(0,1) = 0;
    eval_pts(0,2) = 0.04;
    eval_pts(0,3) = 0.4;

    gsMatrix<> result(2,4);
    gsMatrix<> presult(2,4);

    // Compare evaluations.
    basis -> eval_into(eval_pts, result);
    pbasis -> eval_into(eval_pts, presult);

    if( result != presult )
        printError( result, presult, test_result, "values" );

    // Check the derivatives as well.
    basis -> deriv_into(eval_pts, result);
    pbasis -> deriv_into(eval_pts, presult);

    if( result != presult )
        printError( result, presult, test_result, "derivatives");

    // And also the second derivatives.
    basis -> deriv2_into(eval_pts, result);
    pbasis -> deriv2_into(eval_pts, presult);

    if( result != presult )
        printError( result, presult, test_result, "second derivatives");

}

// Creates C0 cubic curve, refines and samples several points so that user can check in Paraview that the periodic B-spline looks as expected.
void visualTestingC0()
{
    gsKnotVector<> kv(0,1,2,4);//start,end,interior knots, start/end multiplicites of knots1
    gsMatrix<> coefs(6,3);
    coefs << 0,0,0,  0,2,0, 2,2,0, 4,4,0, 4,2,0, -3,1,0 ;

    gsBSpline<> g1( kv, give(coefs), true);

    gsInfo << "Knot vector: " << g1.knots() << "\n";

    gsMatrix<> points = g1.sample(150);
    gsMatrix<> X = points.row(0);
    gsMatrix<> Y = points.row(1);


    // Compare the curve and sampled points on it.
    //gsWriteParaviewPoints( X, Y, "pointsonbsplinecurve0" );
    gsWriteParaview( g1 , "bsplinecurve0", 100, false, true);

    // Display also the control points.
    //g1.controlPointsForParaview( X, Y);
    //gsWriteParaviewPoints( X, Y, "controlPoints0", );
    //gsInfo << "X: " << X << "\n";
    //gsInfo << "Y: " << Y << "\n";

    // Test, whether adding several knots at a time visually changes the curve.
    // Commented values are for testing assertions.
    gsInfo << "knot vector before refinement:\n" << g1.knots() << "\n";
    gsKnotVector<> newKnots;
    //newKnots.insert( -1 );
    //newKnots.insert( 0 );
    newKnots.insert(-0.1 );
    newKnots.insert( 0.1 );
    newKnots.insert( 0.2 );
    newKnots.insert( 0.3 );
    //newKnots.insert( 1 );
    newKnots.insert( 1.05 );
    //newKnots.insert( 5 );
    g1.insertKnots(newKnots.begin(), newKnots.end() );


    // We can try inserting one value alone as well.
    g1.insertKnot(1.04);
    gsInfo << "Knot vector after refinement: " << g1.knots() << "\n";

    // For displaying the curve and its control points.
    gsWriteParaview( g1, "bsplinecurve1", 100, false, true );
    //g1.controlPointsForParaview( X, Y );
    //gsWriteParaviewPoints( X, Y, "controlPoints1" );

    //Call paraview on exit
    char cmd1[100];
    strcpy(cmd1,"paraview bsplinecurve0.vts\0");
    strcat(cmd1," &");
    gsInfo << "paraview returned: " << system(cmd1) << "\n";
}

// C3 quartic curve is twice refined and user then can check in Paraview that the shape hasn't changed."
void visualTestingC3()
{
    std::string filename = "hat_closedcurve.xml";
    gsFileData<> data(filename);
    gsBSpline<> curve= *data.getAnyFirst< gsBSpline<> >();
    curve.setPeriodic();

    gsWriteParaview( curve, "hatOriginal", 120 );
    curve.insertKnot( 0.71 );
    gsWriteParaview( curve, "hatRefined", 120 );

    // Refinement now with the standard vector to check that the iterator stuff works in more than one situation only.
    std::vector<real_t> newKnots(5,0);
    newKnots[0] = -0.01;
    newKnots[1] =  0.1 ;
    newKnots[2] =  0.2 ;
    newKnots[3] =  0.95;
    newKnots[4] =  1.04;

    curve.insertKnots( newKnots.begin(), newKnots.end() );
    gsWriteParaview( curve, "hatDoubleRefined", 120, false, true );

    // Compare hatOriginal.vts, hatRefined.vts and hatDoubleRefined.vts in Paraview.
    // If everything works well, they should look the same.
}

int main(int argc, char *argv[])
{
    bool plot = false; // If set to true, paraview file is generated and launched on exit
    gsCmdLine cmd("gsPeriodicBSpline_test tests evaluation, refinement and plotting of periodic (closed) BSpline curves.");
    cmd.addSwitch( "plot", "Perform additional tests and launch ParaView with some of the results.", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    int test_result = 0; // Number of failed tests in this file.

    automatedTesting( test_result );

    if( plot )
    {
        visualTestingC0();
        visualTestingC3();
    }

    if( test_result > 0 )
    {
        gsWarn << test_result << " failed tests. Please, try and correct the code.\n";
        return 1;
    }
    else
        return 0;
}
