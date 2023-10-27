// Optimize a Parameterization with respect to the numerical stability
// Elisabeth Pilgerstorfer

#include <iostream>

#include <gismo.h>

#include <gsOptimizer/gsOptParamConditionNumber.h>



using namespace gismo;


int main(int argc, char *argv[])
{
    bool plot(false); // If set to true, paraview file is generated and launched on exit
    GISMO_UNUSED(plot);

    std::string fn("");
    
    int n_iter = 40;
    
    gsCmdLine cmd("Tutorial 01 shows the use of BSpline curves.");
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)",fn);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
     
    // small example
    std::vector<real_t> knots(12);
    knots[0]=0;
    knots[1]=0;
    knots[2]=0;
    knots[3]=0.5;
    knots[4]=0.75;
    knots[5]=0.875;
    knots[6]=0.9375;
    knots[7]=0.96875;
    knots[8]=0.984375;
    knots[9]=1;
    knots[10]=1;
    knots[11]=1;

    gsKnotVector<> KKV(give(knots),2);
    gsMatrix<> coefs(9,1);
    coefs << 0,0.25,0.625,0.8125,0.90625,0.953125,0.976563, 0.992188,1;
    
    gsInfo<<"The knots are :"<< KKV <<"\n";
    gsInfo<<"The coefs are :\n"<< coefs <<"\n";
    
    gsGeometry<> * bsp = new gsBSpline<>( KKV, coefs ); // =  gsReadfile<>(...)
    gsOptParamConditionNumber<real_t> optimizer(*bsp, n_iter);
    
    optimizer.run();

    gsInfo << "thse are the final coefs:\n" << optimizer.result().coefs() << "\n";    
    gsInfo << "this is the final bound " << optimizer.lastComputedBound() << "\n";    



    delete bsp;

    return 0;
}

