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

template<class T> gsMultiPatch<T> BSplineStep2D(T const & a = 1, T const & b = 2, T const & a_in = 1);
template<class T> gsTensorBSpline<2, T> BSplineRect(int deg = 1, const T llx = 0, const T lly = 0, const T a = 1, const T b = 1);
template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, T uMax = 1);
template<class T> void refineBasis_NS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal, bool doubleChannel = false);
template<class T> void computeSteadyNS(uwbINSSolverSteady<T>& navStokesSteady, int numIterNSSteady, int plot_pts = 10000);

int main(int argc, char *argv[])
{
    // ========================================= Settings ========================================= 
    //---------------------------------------------------------------------------------------------
    gsInfo << "Reading file settings. \n";
    gsVector<int> inInt(14);
    gsVector<real_t> inRealT(6);
    gsVector<bool> inBool(12);

    //readInitialSetting(MOTOR_DATA_DIR "uwb-pilsen/initialSettingsINSstep_AFC.txt", inInt, inRealT, inString, inBool);
    readInitialSetting("initialSettingsINSstep_AFC.txt", inInt, inRealT, inBool);

    // ========================================= Settings =========================================
    bool plot = inBool(0);
    bool INSafc = inBool(1);
    bool INSafcHO = inBool(2);
    bool solveSteadyNS = inBool(3);
    bool animate = inBool(4);
    bool SUPG = inBool(5);
    bool SRBAV = inBool(6);
    bool dirElemLength = inBool(7);
    bool TCSD = inBool(8);
    bool SRBAVscaleRes = inBool(9);
    bool loadIC = inBool(10);
    bool CrankNicholson = inBool(11);

    int plot_pts = inInt(0);
    int numIterNSSteady = inInt(1); // max number of time steps
    int numIter = inInt(2); // max number of time steps
    int maxRANSPicardIt = inInt(3);
    int numRefine = inInt(4);
    int numRefineLocal = inInt(5);
    int numThreads = inInt(6); // number of threads for assembly
    int animateStep = inInt(7);
    int tauType = inInt(8);
    int SRBAVtype = inInt(9);
    int SRBAVresidualType = inInt(10);
    int hDirType = inInt(11);
    int deg = inInt(12);
    int numElevate = inInt(13);

    real_t viscositySteady = inRealT(0);
    real_t viscosity = inRealT(1);
    real_t timeStep = inRealT(2); // time step for unsteady computation
    real_t tol = inRealT(3); // stopping tolerance
    real_t uMax = inRealT(4); // inlet velocity maximum
    real_t SRBAValpha = inRealT(5);
    //real_t SRBAVscaleFactorRes = inRealT(6);


    bool doubleChannel = false;

    //----------------------------------------------------

    /*bool plot = true;
    int plot_pts = 10000;

    real_t viscositySteady = 0.1;
    int numIterNSSteady = 10; // max number of time steps

    real_t viscosity = 0.0001;
    real_t timeStep = 0.01; // time step for unsteady computation
    real_t tol = 1e-5; // stopping tolerance
    int numIter = 500; // max number of time steps
    int maxRANSPicardIt = 5;

    bool INSafc = true;
    bool INSafcHO = false;

    int numRefine = 2;
    int numRefineLocal = 2;
    bool solveSteadyNS = true;
    real_t uMax = 2; // inlet velocity maximum
    int numThreads = 1; // number of threads for assembly

    bool animate = false;
    int animateStep = 250;*/
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

    gsInfo << "Solving the backward step RANS example.\n";

    // ========================================= Define problem ========================================= 

    gsBoundaryConditions<> bcInfo;
    defineBCs_NS(bcInfo, uMax);

    gsFunctionExpr<> f("0", "0", 2);

    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;

    real_t a = 8;
    if (doubleChannel)
        a = 16;
    real_t b = 2;
    real_t a_in = 1;

    patches = BSplineStep2D<real_t>(a, b, a_in);

    gsInfo << patches << "\n";

    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> tbasis(patches);

    tbasis.setDegree(deg);

    gsInfo << "tbasis.basis(0) = \n" << tbasis.basis(0) << "\n";
    
    refineBasis_NS(tbasis, numRefine, numRefineLocal, doubleChannel);

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
    uwbINSPde<real_t> NSpdeSteady(patches, bcInfo, f, viscositySteady);
    uwbINSSolverParams<real_t> paramsSteady(NSpdeSteady, discreteBases, opt);
    uwbINSSolverSteady<real_t> navStokesSteady(paramsSteady); // steady coupled solver

    if (!loadIC && solveSteadyNS)
        computeSteadyNS(navStokesSteady, numIterNSSteady, plot_pts);

    //==================================== compute aspect ratio  ==================================================
    real_t maxAspectRatio = navStokesSteady.computeAspectRatio();
    gsInfo << "maxAspectRatio = " << maxAspectRatio << "\n";

    // ========================================= Define solver ========================================= 

    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    params.settings().set(constantsINS::timeStep, timeStep);
    params.settings().set(constantsINS::unst_innerIt, maxRANSPicardIt);

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
        real_t L_step = b/2.;
        real_t SRBAVscaleFactorRes = 1.;
        if (SRBAVscaleRes)
            SRBAVscaleFactorRes = L_step / math::pow(uMax, 2);
        params.settings().set(constantsINS::srbavScaleFactorRes, SRBAVscaleFactorRes);
    }
    else if (TCSD)
    {
        params.settings().set(constantsINS::TCSD, TCSD);
        params.settings().set(constantsINS::tauStabTypeSUPG, tauType);
    }

    params.settings().set(constantsINS::dirElemLength, dirElemLength);
    params.settings().set(constantsINS::hDirType, hDirType);

    // solver
    uwbINSSolverUnsteady<real_t> navStokes(params);

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
        navStokes.setInitialCondition(NSSolution);
        gsInfo << "Solution setted.\n";
    }
    else if (solveSteadyNS)
    {
        gsMatrix<> steadySol = navStokesSteady.getSolution();
        if (INSafc || INSafcHO)
            steadySol = navStokes.getSolution_full(steadySol);
        navStokes.setInitialCondition(steadySol);
        gsInfo << "Steady solution setted.\n";
    }
    /*else
    {
        navStokes.setStokesInitialCondition();
    }*/

    navStokes.initialize();

    if (animate)
        navStokes.solveWithAnimation(numIter, animateStep, tol, plot_pts);
    else
        navStokes.solve(numIter, tol); // solution change norm tol = 10^(-5)

    real_t Tassembly = navStokes.getAssemblyTime();
    real_t Tsolve = navStokes.getSolverSetupTime();
    real_t Tsetupsolve = navStokes.getSolveTime();
    gsInfo << "Assembly time:" << Tassembly << "\n";
    gsInfo << "Solve time:" << Tsolve << "\n";
    gsInfo << "Solver setup time:" << Tsetupsolve << "\n";

    gsFileData<real_t> fd;
    fd << navStokes.getSolution();
    fd.save("step_NSsolution.xml");

    // Optionally plot solution in paraview
    if (plot)
    {
        gsField<> velocity = navStokes.constructSolution(0);
        gsField<> pressure = navStokes.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, "step_velocity", plot_pts);
        gsWriteParaview<>(pressure, "step_pressure", plot_pts);

        gsMatrix<> sol = navStokes.getSolution();
        gsMatrix<> oldSol = navStokes.getAssembler()->getSolution();
        navStokes.getAssembler()->plotResiduum("reziduum", sol, oldSol, plot_pts);
    }


    return 0; 
}

template<class T> gsMultiPatch<T> BSplineStep2D(T const & a, T const & b, T const & a_in)
{
    gsMultiPatch<T> mp;

    //mp.addPatch(BSplineRect(2, 0.0, 0.0, a, b / 2));
    //mp.addPatch(BSplineRect(2, 0.0, b / 2, a, b));
    //mp.addPatch(BSplineRect(2, -a_in, b / 2, 0.0, b));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, 0.0, a, b / 2));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, b / 2, a, b));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-a_in, b / 2, 0.0, b));

    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(1, boundary::west, 2, boundary::east);
    mp.addAutoBoundaries();

    return mp;
}

template<class T>
gsTensorBSpline<2, T> BSplineRect(int deg, const T llx, const T lly, const T a, const T b) // llx - lower left x, lly - lower left y
{
    gsKnotVector<T> kv(0, 1, 0, deg + 1); // first, last, inter, mult_end

    int n = kv.size() - deg - 1;
    gsMatrix<T> coef(n*n, 2);

    switch (deg)
    {
    case 1:
    {
        coef << llx + 0, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + b,
            llx + a, lly + b;
        break;
    }
    case 2:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 2), lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 2),
            llx + (a / 2), lly + (b / 2),
            llx + a, lly + (b / 2),
            llx + 0, lly + b,
            llx + (a / 2), lly + b,
            llx + a, lly + b;
        break;
    }
    case 3:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 3), lly + 0,
            llx + (2. / 3) * a, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 3),
            llx + (a / 3), lly + (b / 3),
            llx + (2. / 3) * a, lly + (b / 3),
            llx + a, lly + (b / 3),
            llx + 0, lly + (2. / 3) * b,
            llx + (a / 3), lly + (2. / 3) * b,
            llx + (2. / 3) * a, lly + (2. / 3) * b,
            llx + a, lly + (2. / 3) * b,
            llx + 0, lly + b,
            llx + (a / 3), lly + b,
            llx + (2. / 3) * a, lly + b,
            llx + a, lly + b;
        break;
    }
    case 4:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 4), lly + 0,
            llx + (2. / 4) * a, lly + 0,
            llx + (3. / 4) * a, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 4),
            llx + (a / 4), lly + (b / 4),
            llx + (2. / 4) * a, lly + (b / 4),
            llx + (3. / 4) * a, lly + (b / 4),
            llx + a, lly + (b / 4),
            llx + 0, lly + (2. / 4) * b,
            llx + (a / 4), lly + (2. / 4) * b,
            llx + (2. / 4) * a, lly + (2. / 4) * b,
            llx + (3. / 4) * a, lly + (2. / 4) * b,
            llx + a, lly + (2. / 4) * b,
            llx + 0, lly + (3. / 4) * b,
            llx + (a / 4), lly + (3. / 4) * b,
            llx + (2. / 4) * a, lly + (3. / 4) * b,
            llx + (3. / 4) * a, lly + (3. / 4) * b,
            llx + a, lly + (3. / 4) * b,
            llx + 0, lly + b,
            llx + (a / 4), lly + b,
            llx + (2. / 4) * a, lly + b,
            llx + (3. / 4) * a, lly + b,
            llx + a, lly + b;
        break;
    }
    case 5:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 5), lly + 0,
            llx + (2. / 5) * a, lly + 0,
            llx + (3. / 5) * a, lly + 0,
            llx + (4. / 5) * a, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 5),
            llx + (a / 5), lly + (b / 5),
            llx + (2. / 5) * a, lly + (b / 5),
            llx + (3. / 5) * a, lly + (b / 5),
            llx + (4. / 5) * a, lly + (b / 5),
            llx + a, lly + (b / 5),
            llx + 0, lly + (2. / 5) * b,
            llx + (a / 5), lly + (2. / 5) * b,
            llx + (2. / 5) * a, lly + (2. / 5) * b,
            llx + (3. / 5) * a, lly + (2. / 5) * b,
            llx + (4. / 5) * a, lly + (2. / 5) * b,
            llx + a, lly + (2. / 5) * b,
            llx + 0, lly + (3. / 5) * b,
            llx + (a / 5), lly + (3. / 5) * b,
            llx + (2. / 5) * a, lly + (3. / 5) * b,
            llx + (3. / 5) * a, lly + (3. / 5) * b,
            llx + (4. / 5) * a, lly + (3. / 5) * b,
            llx + a, lly + (3. / 5) * b,
            llx + 0, lly + (4. / 5) * b,
            llx + (a / 5), lly + (4. / 5) * b,
            llx + (2. / 5) * a, lly + (4. / 5) * b,
            llx + (3. / 5) * a, lly + (4. / 5) * b,
            llx + (4. / 5) * a, lly + (4. / 5) * b,
            llx + a, lly + (4. / 5) * b,
            llx + 0, lly + b,
            llx + (a / 5), lly + b,
            llx + (2. / 5) * a, lly + b,
            llx + (3. / 5) * a, lly + b,
            llx + (4. / 5) * a, lly + b,
            llx + a, lly + b;
        break;
    }
    case 6:
    {
        coef << llx + 0, lly + 0,
            llx + (a / 6), lly + 0,
            llx + (2. / 6) * a, lly + 0,
            llx + (3. / 6) * a, lly + 0,
            llx + (4. / 6) * a, lly + 0,
            llx + (5. / 6) * a, lly + 0,
            llx + a, lly + 0,
            llx + 0, lly + (b / 6),
            llx + (a / 6), lly + (b / 6),
            llx + (2. / 6) * a, lly + (b / 6),
            llx + (3. / 6) * a, lly + (b / 6),
            llx + (4. / 6) * a, lly + (b / 6),
            llx + (5. / 6) * a, lly + (b / 6),
            llx + a, lly + (b / 6),
            llx + 0, lly + (2. / 6) * b,
            llx + (a / 6), lly + (2. / 6) * b,
            llx + (2. / 6) * a, lly + (2. / 6) * b,
            llx + (3. / 6) * a, lly + (2. / 6) * b,
            llx + (4. / 6) * a, lly + (2. / 6) * b,
            llx + (5. / 6) * a, lly + (2. / 6) * b,
            llx + a, lly + (2. / 6) * b,
            llx + 0, lly + (3. / 6) * b,
            llx + (a / 6), lly + (3. / 6) * b,
            llx + (2. / 6) * a, lly + (3. / 6) * b,
            llx + (3. / 6) * a, lly + (3. / 6) * b,
            llx + (4. / 6) * a, lly + (3. / 6) * b,
            llx + (5. / 6) * a, lly + (3. / 6) * b,
            llx + a, lly + (3. / 6) * b,
            llx + 0, lly + (4. / 6) * b,
            llx + (a / 6), lly + (4. / 6) * b,
            llx + (2. / 6) * a, lly + (4. / 6) * b,
            llx + (3. / 6) * a, lly + (4. / 6) * b,
            llx + (4. / 6) * a, lly + (4. / 6) * b,
            llx + (5. / 6) * a, lly + (4. / 6) * b,
            llx + a, lly + (4. / 6) * b,
            llx + 0, lly + (5. / 6) * b,
            llx + (a / 6), lly + (5. / 6) * b,
            llx + (2. / 6) * a, lly + (5. / 6) * b,
            llx + (3. / 6) * a, lly + (5. / 6) * b,
            llx + (4. / 6) * a, lly + (5. / 6) * b,
            llx + (5. / 6) * a, lly + (5. / 6) * b,
            llx + a, lly + (5. / 6) * b,
            llx + 0, lly + b,
            llx + (a / 6), lly + b,
            llx + (2. / 6) * a, lly + b,
            llx + (3. / 6) * a, lly + b,
            llx + (4. / 6) * a, lly + b,
            llx + (5. / 6) * a, lly + b,
            llx + a, lly + b;
        break;
    }
    default:
        GISMO_ERROR("Degree not implemented.");
        break;
    }

    return gsTensorBSpline<2, T>(kv, kv, give(coef));
}

template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, T uMax)
{
    std::string xVel = util::to_string(uMax) + " * (-4*(y-1.5)^2 + 1)";
    gsFunctionExpr<T> Uin(xVel, "0", 2); // inlet velocity
    gsFunctionExpr<T> Uwall("0", "0", 2); // wall velocity

    bool pressureBC = false;

    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Uwall, 0);
    //bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uin, 0);

    if (pressureBC)
    {
        gsFunctionExpr<> P("0", 2);
        gsFunctionExpr<> Pin("1.", 2);
        bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, P, 1);
        bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, P, 1);
        bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Pin, 1);
    }
    else
        bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uin, 0);
}

template<class T>
void computeSteadyNS(uwbINSSolverSteady<T>& navStokesSteady, int numIterNSSteady, int plot_pts)
{
    gsInfo << "Solving Steady case: \n";
    gsInfo << "numDofs: " << navStokesSteady.numDofs() << "\n";

    navStokesSteady.initialize(); // steady solver
    navStokesSteady.solve(numIterNSSteady, 1e-5);

    gsField<> velocitySteady = navStokesSteady.constructSolution(0);
    gsField<> pressureSteady = navStokesSteady.constructSolution(1);
    gsWriteParaview<>(velocitySteady, "step_SteadyVelocity", plot_pts, true);
    gsWriteParaview<>(pressureSteady, "step_SteadyPressure", plot_pts);

    gsFileData<real_t> fd;
    fd << navStokesSteady.getSolution();
    fd.save("step_NSsteadySolution.xml");
}

template<class T> void refineBasis_NS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal, bool doubleChannel)
{
    gsMatrix<> box(2, 2);
    box << 0, 1, 0, 0;

    int refP = 0;
    if (doubleChannel)
        refP = 1;
    for (int i = 0; i < 3+refP; i++)
    {
        basis.refine(0, box);
        basis.refine(1, box);
    }

    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    // refinement near wal
    int numRefineLocalWal = 1;

    real_t parArea = 0.1;
    real_t parArea2 = 0.02;

    for (int j = 0; j < numRefineLocal; j++)
    {
        box << 0, 0, 0, parArea;
        for (int i = 0; i < numRefineLocalWal; i++)
        {
            basis.refine(0, box);
            basis.refine(1, box);
            basis.refine(2, box);
        }

        box << 0, 0, 1 - parArea, 1;
        for (int i = 0; i < numRefineLocalWal; i++)
        {
            basis.refine(0, box);
            basis.refine(1, box);
            basis.refine(2, box);
        }

        box << 0, parArea2, 0, 0;
        basis.refine(0, box);
        basis.refine(1, box);

        box << 1-parArea, 1, 0, 0;
        basis.refine(2, box);

        parArea = parArea / 2;
        parArea2 = parArea2 / 2;
    }
}

