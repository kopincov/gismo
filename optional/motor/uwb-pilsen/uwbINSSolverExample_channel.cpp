/** @file uwbINSSolversExample.cpp

Author(s): E. Turnerova
*/

//#pragma once

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>

// solvers
#include "uwbINSSolverSteady.h"
#include "uwbINSSolverUnsteady.h"

//#include "uwbRANSExamplesSetting.h"
#include "uwbGeometryCreators.h"
#include "uwbDataFileReader.h"


using namespace gismo;

template<class T> gsMultiPatch<T> BSplineGeometry(T const & ax = 1, T const & bx = 2, T const & ay = 1, T const & by = 2, int const & numPatchDiv = 0);
template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, T Re = 1, bool Hughes = false, bool classicalChannel = false);
template<class T> void refineBasis_NS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal);

int main(int argc, char *argv[])
{
    // ========================================= Settings =========================================
    //---------------------------------------------------------------------------------------------
    gsInfo << "Reading file settings. \n";
    gsVector<int> inInt(11);
    gsVector<real_t> inRealT(3);
    gsVector<bool> inBool(9);

    //readInitialSetting(MOTOR_DATA_DIR "uwb-pilsen/initialSettingsINSchannel.txt", inInt, inRealT, inString, inBool);
    readInitialSetting("initialSettingsINSchannel.txt", inInt, inRealT, inBool);

    // ========================================= Settings =========================================
    bool plot = inBool(0);
    bool INSafc = inBool(1);
    bool INSafcHO = inBool(2);
    bool SUPG = inBool(3);
    bool SRBAV = inBool(4);
    bool dirElemLength = inBool(5);
    bool TCSD = inBool(6);
    bool SRBAVscaleRes = inBool(7);
    bool loadIC = inBool(8);

    int plot_pts = inInt(0);
    int numIter = inInt(1); // max number of iterations
    int numRefine = inInt(2);
    int numRefineLocal = inInt(3);
    int numThreads = inInt(4); // number of threads for assembly
    int tauType = inInt(5);
    int SRBAVtype = inInt(6);
    int SRBAVresidualType = inInt(7);
    int hDirType = inInt(8);
    int deg = inInt(9);
    int numElevate = inInt(10);

    //real_t viscosity = inRealT(0);
    real_t tol = inRealT(0); // stopping tolerance
    real_t SRBAValpha = inRealT(1);
    real_t Re = inRealT(2);
    //real_t SRBAVscaleFactorRes = inRealT(6);

    real_t viscosity = 1. / Re;//2.6 * 2. / Re;

    bool Hughes = true;
    bool classicalChannel = true;

    //==================================================

    //command line
    gsCmdLine cmd("Solves Navier-Stokes problem with an isogeometric discretization.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addInt("r", "uniformRefine",
        "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("l", "localRefine",
        "Number of local h-refinement steps to perform before solving", numRefineLocal);
    cmd.addInt("s", "plotSamples",
        "Number of sample points to use for plotting", plot_pts);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (numRefine < 0)
    {
        gsInfo << "Number of refinements must be non-negative, quitting.\n";
        return -1;
    }

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    //opt.dirStrategy = dirichlet::none;//elimination;
    opt.intStrategy = iFace::glue;

    gsInfo << "Solving the steady channel example.\n";

    // ========================================= Define problem =========================================

    gsBoundaryConditions<> bcInfo;
    defineBCs_NS(bcInfo, Re, Hughes, classicalChannel);

    //lambda = -0.963740544195767032178
    //3.141592653589793238463
    //nu=1/Re = 0.025// 0.13

    std::string f1, f2;
    if (classicalChannel)
    {
        f1 = "0";
        f2 = "0";
    }
    else
    {
        if (Hughes)
            f1 = "(0.025*((-0.963740544195767032178)*(-0.963740544195767032178) - 4*3.141592653589793238463*3.141592653589793238463) + 0.963740544195767032178)*exp(-0.963740544195767032178*x)*cos(2*3.141592653589793238463*y)";
        else
            f1 = "(0.025*((-0.963740544195767032178)*(-0.963740544195767032178) - 4*3.141592653589793238463*3.141592653589793238463) + 0.963740544195767032178)*exp(-0.963740544195767032178*x)*cos(2*3.141592653589793238463*y) - 2*0.963740544195767032178*exp(-0.963740544195767032178*2*x)";

        f2 = "(0.025*(2*3.141592653589793238463*(-0.963740544195767032178) - (-0.963740544195767032178)*(-0.963740544195767032178)*(-0.963740544195767032178)/(2*3.141592653589793238463)) + (-0.963740544195767032178)*(-0.963740544195767032178)/(2*3.141592653589793238463))*exp(-0.963740544195767032178*x)*sin(2*3.141592653589793238463*y)";
    }


    //gsFunctionExpr<> f("0", "0", 2);
    gsFunctionExpr<> f(f1, f2, 2);

    // ========================================= Define geometry =========================================

    gsMultiPatch<> patches;

    real_t ax, bx, ay, by;

    if (Hughes)
    {
        ax = 0;
        bx = 1.0;
        ay = -0.5;
        by = 0.5;
    }
    else
    {
        ax = -0.5;
        bx = 1.;
        ay = -0.5;
        by = 1.5;
    }

    int numPatchDiv = 0;

    patches = BSplineGeometry<real_t>(ax, bx, ay, by, numPatchDiv);

    gsInfo << patches << "\n";

    // ========================================= Define basis =========================================

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> tbasis(patches);

    tbasis.setDegree(deg);

    //gsInfo << "tbasis.basis(0) = \n" << tbasis.basis(0) << "\n";
    gsInfo << "\n-----------------------------------------\n";
    for (unsigned p = 0; p < patches.nPatches(); p++)
        gsInfo << "tbasis.basis(p) = \n" << tbasis.basis(p) << "\n";
    gsInfo << "-----------------------------------------\n\n";

    refineBasis_NS(tbasis, numRefine, numRefineLocal);

    if (numElevate > 0)
    {
        tbasis.degreeElevate(numElevate);
        gsInfo << "\nBasis elevated.\n\n";
    }
    gsInfo << "tbasis.basis(0) = \n" << tbasis.basis(0) << "\n";

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    gsInfo << "\nVelocity basis degree: " << discreteBases[0].degree() << "\n";
    gsInfo << "Pressure basis degree: " << discreteBases[1].degree() << "\n\n";

    // ========================================= Compute Steady NS =========================================
    //uwbINSPde<real_t> NSpdeSteady(patches, bcInfo, f, viscositySteady);
    //uwbINSSolverParams<real_t> paramsSteady(NSpdeSteady, discreteBases, opt);
    //uwbINSSolverSteady<real_t> navStokesSteady(paramsSteady); // steady coupled solver

    //if (!loadIC && solveSteadyNS)
    //    computeSteadyNS(navStokesSteady, numIterNSSteady, plot_pts);

    //==================================== compute aspect ratio  ==================================================
    //real_t maxAspectRatio = navStokesSteady.computeAspectRatio();
    //gsInfo << "maxAspectRatio = " << maxAspectRatio << "\n";

    // ========================================= Define solver =========================================

    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    //params.settings().set(constantsINS::timeStep, timeStep);
    //params.settings().set(constantsINS::unst_innerIt, maxRANSPicardIt);

    // --------- other options ---------
    params.setNumThreads(numThreads); // set the number of threads for assembly
 //   params.settings().set(constantsINS::unst_innerIt, 5); // max number of inner Picard iter. for unstready
 //   params.settings().set(constantsINS::unst_innerTol, 1e-4); // stopping tolerance for inner Picard iter. in unstready
 // ---------------------------------

    if (INSafc)
        params.settings().set(constantsINS::AFC, INSafc);
    if (INSafcHO)
        params.settings().set(constantsINS::AFC_HO, INSafcHO);
    if (SUPG)
    {
        params.settings().set(constantsINS::SUPG, SUPG);
        params.settings().set(constantsINS::tauStabTypeSUPG, tauType);
    }
    else if (SRBAV)
    {
        params.settings().set(constantsINS::SRBAV, SRBAV);
        params.settings().set(constantsINS::SRBAVtype, SRBAVtype);
        params.settings().set(constantsINS::tauStabTypeSRBAV, tauType);
        params.settings().set(constantsINS::SRBAVresidualType, SRBAVresidualType);
        params.settings().set(constantsINS::SRBAValpha, SRBAValpha);
        //real_t L_step = b/2.;
        real_t SRBAVscaleFactorRes = 1.;
        //if (SRBAVscaleRes)
        //    SRBAVscaleFactorRes = L_step / math::pow(uMax, 2);
        params.settings().set(constantsINS::srbavScaleFactorRes, SRBAVscaleFactorRes);
    }
    else if (TCSD)
    {
        params.settings().set(constantsINS::TCSD, TCSD);
        params.settings().set(constantsINS::tauStabTypeSUPG, tauType);
    }

    params.settings().set(constantsINS::dirElemLength, dirElemLength);
    params.settings().set(constantsINS::hDirType, hDirType);

    //params.settings().set(constantsINS::CrankNicholson, CrankNicholson);

    // solver
    uwbINSSolverSteady<real_t> navStokes(params);

    // ========================================= Solving =========================================
    gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

    gsMatrix<real_t> NSSolution;
    //------------- read from file ---------------------
    if (loadIC)
    {
        gsFileData<> fdRead("ICsolution.xml");
        NSSolution = *(fdRead.getFirst< gsMatrix<real_t> >());
        gsInfo << "Solution loaded.\n";
    }
    //------------------------------------------------

    gsInfo << "initialization...\n";
    //navStokes.initialize();

    if (loadIC)
    {
        navStokes.setSolution(NSSolution);
        gsInfo << "Solution setted.\n";
    }
    /*else
    {
        navStokes.setStokesInitialCondition();
    }*/

    navStokes.initialize();

    navStokes.solve(numIter, tol); // solution change norm tol = 10^(-5)

    real_t Tassembly = navStokes.getAssemblyTime();
    real_t Tsolve = navStokes.getSolverSetupTime();
    real_t Tsetupsolve = navStokes.getSolveTime();
    gsInfo << "Assembly time:" << Tassembly << "\n";
    gsInfo << "Solve time:" << Tsolve << "\n";
    gsInfo << "Solver setup time:" << Tsetupsolve << "\n";

    gsField<> velocity = navStokes.constructSolution(0);
    gsField<> pressure = navStokes.constructSolution(1);
    //navStokes.getAssembler()->computeError(velocity, pressure, Hughes);

    gsMatrix<> sol = navStokes.getSolution();
    //navStokes.getAssembler()->plotError("error", sol, plot_pts, Hughes);

    gsFileData<real_t> fd;
    fd << navStokes.getSolution();
    fd.save("channel_NSsolution.xml");

    // Optionally plot solution in paraview
    if (plot)
    {
        //gsField<> velocity = navStokes.constructSolution(0);
        //gsField<> pressure = navStokes.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, "channel_velocity", plot_pts);
        gsWriteParaview<>(pressure, "channel_pressure", plot_pts);
        gsWriteParaview<>(pressure, "channel_pressure_mesh", plot_pts, true);

        //gsMatrix<> sol = navStokes.getSolution();
        navStokes.getAssembler()->plotResiduum("reziduum", sol, plot_pts);

        /*std::vector<real_t> energy = navStokes.getEnergy();
        gsMatrix<> energyGs;
        int numTimeSteps = energy.size();
        energyGs.setZero(numTimeSteps, 1);
        for (int i = 0; i < numTimeSteps; i++)
            energyGs(i, 0) = energy[i];
        gsFileData<real_t> en;
        en << energyGs;
        en.save("step_energy.xml");*/
    }


    return 0;
}

template<class T> gsMultiPatch<T> BSplineGeometry(T const & ax, T const & bx, T const & ay, T const & by, int const & numPatchDiv)
{
    gsMultiPatch<T> mp;

    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(ax, ay, bx, by));

    //for (int i = 0; i < numPatchDiv; i++)
    //{

    //}

    //mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, 0.0, a, b / 2));
    //mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, b / 2, a, b));
    //mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-a_in, b / 2, 0.0, b));

    //mp.addInterface(0, boundary::north, 1, boundary::south);
    //mp.addInterface(1, boundary::west, 2, boundary::east);

    mp.addAutoBoundaries();

    return mp;
}

template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, T Re, bool Hughes, bool classicalChannel)
{
    if (classicalChannel)
    {
        std::string xVel_in = "1 - 4*y*y";

        gsFunctionExpr<T> Uin(xVel_in, "0", 2); // inlet velocity
        gsFunctionExpr<T> Uwall("0", "0", 2); // wall velocity

        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);
        bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
    }
    else
    {
        /*const T PI = 3.141592653589793238463;
        T PIsquared = PI*PI;
        T ReSquared = Re * Re;
        std::string PIStr = std::to_string(PI);
        std::string PIsquaredStr = std::to_string(PIsquared);
        std::string ReStr = std::to_string(Re);
        std::string ReSquaredStr = std::to_string(ReSquared);
        std::string lambda = ReStr + "/2 - sqrt(" + ReSquaredStr + "/4 + 4*" + PIsquaredStr + ")";
        */
        //std::string xVel = "1 - exp(" + lambda + "*x)*cos(2*" + PIStr + "*y)";
        //std::string yVel = lambda + "/(2*" + PIStr + ")*exp(" + lambda + "*x)*sin(2*" + PIStr + "*y)";
        //std::string xVel_in = "1 - exp(" + lambda + "*" + std::to_string(-0.5) + ")*cos(2*" + PIStr + "*y)";
        //std::string yVel_in = lambda + "/(2*" + PIStr + ")*exp(" + lambda + "*" + std::to_string(-0.5) + ")*sin(2*" + PIStr + "*y)";
        //std::string xVel_low = "1 - exp(" + lambda + "*x)*cos(2*" + PIStr + "*" + std::to_string(-0.5) + ")";
        //std::string xVel_up = "1 - exp(" + lambda + "*x)*cos(2*" + PIStr + "*" + std::to_string(1.5) + ")";

        //lambda = -0.963740544195767032178
        //3.141592653589793238463
        //*
        //std::string xVel = "1 - exp(-0.963740544195767032178*x)*cos(2*3.141592653589793238463*y)";
        //std::string yVel = "(-0.963740544195767032178/(2*3.141592653589793238463))*exp(-0.963740544195767032178*x)*sin(2*3.141592653589793238463*y)";

        std::string xVel_in, yVel_in, xVel_out, yVel_out, xVel_low, xVel_up;

        if(Hughes)
        {
            xVel_in = "1 - cos(2*3.141592653589793238463*y)";
            yVel_in = "(-0.963740544195767032178/(2*3.141592653589793238463))*sin(2*3.141592653589793238463*y)";
            xVel_out = "1 - exp(-0.963740544195767032178*1.)*cos(2*3.141592653589793238463*y)";
            yVel_out = "(-0.963740544195767032178/(2*3.141592653589793238463))*exp(-0.963740544195767032178*1.)*sin(2*3.141592653589793238463*y)";
            xVel_low = "1 - exp(-0.963740544195767032178*x)*cos(2*3.141592653589793238463*(-0.5))";
            xVel_up = "1 - exp(-0.963740544195767032178*x)*cos(2*3.141592653589793238463*0.5)";
        }
        else
        {
            xVel_in = "1 - exp(-0.963740544195767032178*(-0.5))*cos(2*3.141592653589793238463*y)";
            yVel_in = "(-0.963740544195767032178/(2*3.141592653589793238463))*exp(-0.963740544195767032178**(-0.5))*sin(2*3.141592653589793238463*y)";
            xVel_out = "1 - exp(-0.963740544195767032178*1.)*cos(2*3.141592653589793238463*y)";
            yVel_out = "(-0.963740544195767032178/(2*3.141592653589793238463))*exp(-0.963740544195767032178*1.)*sin(2*3.141592653589793238463*y)";
            xVel_low = "1 - exp(-0.963740544195767032178*x)*cos(2*3.141592653589793238463*(-0.5))";
            xVel_up = "1 - exp(-0.963740544195767032178*x)*cos(2*3.141592653589793238463*1.5)";
        }

        //std::string xVel_in = "1 - exp(-((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*0.5)*cos(2*3.141592653589793238463*y)";
        //std::string yVel_in = "(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))/(2*3.141592653589793238463))*exp(-((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*0.5)*sin(2*3.141592653589793238463*y)";
        //std::string xVel_low = "1 - exp(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*x)*cos(-2*3.141592653589793238463*0.5)";
        //std::string xVel_up = "1 - exp(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*x)*cos(-2*3.141592653589793238463*1.5)";
        //std::string xVel_out = "1 - exp(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*1.)*cos(2*3.141592653589793238463*y)";
        //std::string yVel_out = "(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))/(2*3.141592653589793238463))*exp(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*1.)*sin(2*3.141592653589793238463*y)";


        gsFunctionExpr<T> Uin(xVel_in, yVel_in, 2); // inlet velocity
        gsFunctionExpr<T> UwallLow(xVel_low, "0", 2); // wall velocity
        gsFunctionExpr<T> UwallUp(xVel_up, "0", 2); // wall velocity
        gsFunctionExpr<T> Uout(xVel_out, yVel_out, 2); // inlet velocity

        //bool pressureBC = false;

        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);
        bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, UwallLow, 0);
        bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, UwallUp, 0);
        //bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uout, 0);


        //lambda=((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))
        //(2*3.141592653589793238463)
        //std::string P_out = "(0.13*((2*3.141592653589793238463)*(2*3.141592653589793238463)/((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463))) - ((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463))))+1.)*exp(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*1.)*cos(2*3.141592653589793238463*y) - 0.5*exp(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*2.)*(sin(2*3.141592653589793238463*y)*sin(2*3.141592653589793238463*y) + cos(2*3.141592653589793238463*y)*cos(2*3.141592653589793238463*y))";
        //std::string P_out = "(0.13*(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463))) - ((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))/((2*3.141592653589793238463)*(2*3.141592653589793238463)))+((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))/((2*3.141592653589793238463)*(2*3.141592653589793238463)))*exp(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*1.)*cos(2*3.141592653589793238463*y) - 0.5*exp(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*1.)";


                //"(0.13*(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463))) - ((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))/((2*3.141592653589793238463)*(2*3.141592653589793238463))) + ((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))/((2*3.141592653589793238463)*(2*3.141592653589793238463)).)*exp(((40/2)-sqrt((40/2)*(40/2) + (2*3.141592653589793238463)*(2*3.141592653589793238463)))*1.)*cos(2*3.141592653589793238463*y)";


        //lambda = -0.963740544195767032178
        //3.141592653589793238463
        //nu=1/Re = 0.025// 0.13
        //*
        std::string P_out = "0.5*(1. - exp(2*(-0.963740544195767032178)))";
        //std::string P_out = "0.5*exp(2*(-0.963740544195767032178))";
        //std::string P_out = "(0.13*(-0.963740544195767032178) - 0.13*(-0.963740544195767032178)*(-0.963740544195767032178)*(-0.963740544195767032178)/(4*3.141592653589793238463*3.141592653589793238463) + (-0.963740544195767032178)*(-0.963740544195767032178)/(4*3.141592653589793238463*3.141592653589793238463))*exp(-0.963740544195767032178)*cos(2*3.141592653589793238463*y) + (0.13*(-0.963740544195767032178) - 0.13*(-0.963740544195767032178)*(-0.963740544195767032178)*(-0.963740544195767032178)/(4*3.141592653589793238463*3.141592653589793238463) + (-0.963740544195767032178)*(-0.963740544195767032178)/(4*3.141592653589793238463*3.141592653589793238463))*exp(-0.963740544195767032178)";

        gsFunctionExpr<> P(P_out, 2);
        bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, P, 1);

        /*if (pressureBC)
        {
            gsFunctionExpr<> P("0", 2);
            gsFunctionExpr<> Pin("1.", 2);
            bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, P, 1);
            bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, P, 1);
            bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Pin, 1);
        }
        else
            bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uin, 0);*/
    }
}

template<class T> void refineBasis_NS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal)
{
    //gsMatrix<> box(2, 2);
    //box << 0, 1, 0, 0;
    //basis.refine(0, box);

    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();
}

