
#include <gsBoxSplines/gsBoxSplineBasis.h>
#include <gsIO/gsWriteParaview.h>
#include <gsUtils/gsStopwatch.h>
#include <gsIO/gsCmdLine.h>

#include <iostream>


using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false;
    index_t samples = 1000;
    gsCmdLine cmd("Testing gsBoxSplineBasis.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addInt("s","samples", "Number of samples to use for viewing", samples);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    

    gsInfo<< "\n------------ Box Spline Test ------------\n" ;
    gsStopwatch time;
    double timer = time.stop();
    
    // Direction Matrix
    gsMatrix<> C(2,2);
    C << 1,0,  0,1 ;
    
    // Multiplicity vector
    gsVector<int> M(2);
    M << 3,3 ;
    
    // Some evaluation points
//    gsMatrix<> u(2,2);
//    u << 1.5, 1.0,
//        1.0, 1.0 ;
/*    gsMatrix<> u(2,2);
        u << 1.0, 1.0,
            1.0, 1.0 ; */
    gsMatrix<> u(2,12);
    for(int i = 0; i < 12; i++)
    {
        u(0,i) = 3+((0)/12);
        u(1,i) = 2+((i)/12);
    }

    // Tensor-product box-spline
    gsBoxSplineBasis<2,real_t> BS1(C,M);
    gsInfo << "\n-----------------------------------------\n"
         << BS1 << "\n"
         << "eval: "<<  BS1.eval(u) <<"\n"
         << "-----------------------------------------\n";
    
    // Courant
    gsBoxSplineBasis<2,real_t> BS2(1,1,1);
    gsInfo << "\n--------- Courant element ---------------\n"
         << BS2 << "\n"
         << "eval: "<<  BS2.eval(u) << "\n"
         << "-----------------------------------------\n";
  
    gsBoxSplineBasis<2,real_t> BS3(2,1,1);
    gsInfo << "\n-----------------------------------------\n"
         << BS3 << "\n"
         << "eval: "<<  BS3.eval(u) <<"\n"
         << "-----------------------------------------\n";
    
    gsBoxSplineBasis<2,real_t> BS4(2,2,1);
    gsInfo << "\n-----------------------------------------\n"
         << BS4 <<"\n"
         << "eval: "<<  BS4.eval(u) <<"\n"
         << "-----------------------------------------\n";

    // Loop
    gsBoxSplineBasis<2,real_t> BS5(2,2,2);
    gsInfo << "\n----------- Loop element ----------------\n"
         << BS5 <<"\n"
         << "eval: "<< BS5.eval(u) <<"\n"
         << "-----------------------------------------\n";

    real_t testsum = 0;
    gsMatrix<> result = BS5.evalresult(u);
    for(int i = 0; i < result.cols(); i++)
    {
        gsInfo << result(0,i) << "\n";
        testsum = testsum + result(0,i);
    }
    gsInfo <<"Testsum: " << testsum << "\n";
    // test by Dominik
    gsBoxSplineBasis<2,real_t> BSDerv(2,2,2);
    gsInfo << "\n----------- Derivation ----------------\n"
         << BSDerv <<"\n"
         << "eval: "<< BSDerv.deriv(u) <<"\n"
         << "-----------------------------------------\n";
    // test end

    // ZP
    gsBoxSplineBasis<2,real_t> BS6(1,1,1,1);
    gsInfo << "\n------- Zwart-Powell element ------------\n"
         << BS6 <<"\n"
         << "eval: "<< BS6.eval(u) <<"\n"
         << "-----------------------------------------\n";

    // ZP extended
    gsBoxSplineBasis<2,real_t> BS61(2,2,1,1);
    gsInfo << "\n-----------------------------------------\n"
         << BS61 <<"\n"
         << "eval: "<< BS61.eval(u) <<"\n"
         << "-----------------------------------------\n";

    // test - full matrix given
    gsMatrix<> D(2,4);
    D << 1,1,0,0, 0,0,1,1 ; // rowwise
    gsBoxSplineBasis<2,real_t> BS7(D);
    gsInfo << "\n-----------------------------------------\n"
         << BS7 <<"\n"
         << "eval: "<< BS7.eval(u) <<"\n"
         << "-----------------------------------------\n";


    // Some evaluation points for 3D
    gsMatrix<> u1(3,3);
    u1 << 1.5, 1.0, 2.0,
         1.0, 1.0, 2.0,
         1.0, 1.0, 2.0 ;

    // test - full matrix given 3D!!! (no plotting!)
    gsMatrix<> D1(3,6);
    D1 << 1,1,0,0,0,0, 0,0,1,1,0,0, 0,0,0,0,1,1 ; // rowwise
    gsBoxSplineBasis<3,real_t> BS8(D1);
    gsInfo << "\n-----------------------------------------\n"
         << BS8 <<"\n"
         << "eval: "<< BS8.eval(u1) <<"\n"
         << "-----------------------------------------\n";

    timer = time.stop() - timer;
    gsInfo << "time: "<< timer <<" sec \n" ;
    
    int exitCommand(0);
    if (plot){
        // Write approximate and exact solution to paraview files

        gsInfo<<"Plotting in Paraview...\n";
        double timer2 = time.stop();
        gsWriteParaview<>( BS5, "boxspline_paraview", samples);
        timer2 = time.stop() - timer2;
        gsInfo << "time to plot: "<< timer2 <<" sec \n" ;

        // Run paraview on exit
        exitCommand = system("paraview boxspline_paraview.pvd &");
    }
    
    return exitCommand;
}
