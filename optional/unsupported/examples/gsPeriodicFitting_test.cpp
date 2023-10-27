/// gsPeriodicFitting_test.cpp
/// Author:Dominik Mokris, dominik.mokris@jku.at, based on gsFitting_test.cpp by Mario Kapl
/// Tests curve fitting with a periodic curve.
/// Takes the point cloud of a curve with the parameter values and the desired knot vector
/// and computes a curve approximation curve with the help of least squares fitting.

#include <iostream>

#include <string>

#include <time.h> // Note: use gsUtils/gsStopWatch.h instead (example eg. gsSmoothing_test.cpp)

#include <gismo.h>
#include <gismo_dev.h>


using namespace gismo;

memory::unique_ptr< gsCurveFitting<> > readFittingFromFile( std::string nameOfFileWithFitting, bool isApproximationCurvePeriodic )
{
    std::string filename = nameOfFileWithFitting + ".xml";
    memory::unique_ptr< gsCurveFitting<> > fitting;
    gsFileData<>  data( filename );
    if ( data.has< gsCurveFitting<> >() )
        fitting = data.getFirst< gsCurveFitting<> >();
    else
        gsWarn << "File " << filename << " is missing or corrupted.\n";

    fitting->setClosedCurve( isApproximationCurvePeriodic );
    return fitting;
}

class fittingTest
{
public:
    fittingTest( std::string nameOfFileWithFitting, std::string testId, real_t tolerance, bool isApproximatedCurvePeriodic, bool plot, bool print )
    {
        m_fitting = readFittingFromFile( nameOfFileWithFitting, isApproximatedCurvePeriodic );
        m_testId = testId;
        m_tolerance = tolerance;
        m_plot = plot;
        m_print = print;
    }

    bool result();

private:
    fittingTest(){} // Private to prevent someone calling it without propper initialisation.
    void printResults();
    void performFittingAndMeasureElapsedTime();
    void postprocessResults();

    // Data members
    memory::unique_ptr< gsCurveFitting<> > m_fitting;
    std::string m_testId;
    real_t m_tolerance;
    real_t m_approxError;
    double m_elapsedTime;
    bool m_plot;
    bool m_print;
    bool m_passed;
    gsBSpline<> m_approxCurve;
};

bool fittingTest::result()
{
    performFittingAndMeasureElapsedTime();
    postprocessResults();
    return m_passed;
}

void fittingTest::printResults()
{
    gsInfo << "Got a curve fitting problem "<< *m_fitting << "\n";
    if( !m_passed )
        gsInfo << "FAILED." << "\n";
    else
        gsInfo << "SUCCESS." << "\n";
    gsInfo << "The approximation error of the fitted curve is: "<< m_approxError  << "\n";
    gsInfo << "The tolerance was: " << m_tolerance << "\n";
    gsInfo << "The resulting curve of the fitting problem is a "<< m_approxCurve << "\n";
    gsInfo << "Time of computation in milliseconds: " << m_elapsedTime << "\n";
}

void fittingTest::performFittingAndMeasureElapsedTime()
{
    int initialTime=clock(); //measuring the computational time
    m_fitting->compute_periodic();
    m_elapsedTime = (clock() - initialTime)/1000;
}

void fittingTest::postprocessResults()
{
    m_approxCurve = m_fitting->curve();
    m_fitting->computeApproxError(m_approxError);
    m_passed = (m_approxError <= m_tolerance );

    if( m_plot )
    {
        gsWriteParaview( m_approxCurve , "bsplinecurve" + m_testId, 1000);
    }
    if( m_print || !m_passed )
        printResults();
}

void readCommandLineArguments( int argc, char *argv[], bool& plot, bool& print )
{
    plot = false;
    print = false;
    gsCmdLine cmd("Hi, I perform curve fitting by least squares. Do you want to see the graph of the result? And the approximation error etc.?");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("print", "Describe the results on the screen", print);
    cmd.getValues(argc,argv);
}

int main(int argc, char *argv[])
{    
    bool plot = false;
    bool print = false;
    try { readCommandLineArguments(argc, argv, plot, print); } catch (int rv) { return rv; }

    real_t tol = std::pow(10.0, - REAL_DIG * 0.74);
    gsDebugVar(tol);
    gsDebugVar(math::sqrt(tol));

    fittingTest test0( "fitting_default", "0", tol, true, plot, print );
    fittingTest test1( "fitting1", "1", math::sqrt(tol), true, plot, print );
    fittingTest test2( "fitting2", "2", tol, true, plot, print );
    fittingTest test3( "fitting3", "3", tol, true, plot, print );

    bool passedEverything =
            test0.result() &&
            test1.result() &&
            test2.result() &&
            test3.result();

    return passedEverything ? 0 : 1;
}
