/** @file uwbINSSolverCircle.cpp

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
template<class TT> gsMultiPatch<TT>  BSplineCircle2D4Patches(int degree = 2, bool plotMeshes = false, TT const & length = 2.2, TT const & width = 0.4, TT const & widthExpand = 0.01, TT const & radius = 0.05, TT const & centreX = 0.0, TT const & centreY = 0.0, bool const & prepatch = true, TT const & prepatchWidth = 0.2);
template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, T uMax = 1, bool prepatch = false);
template<class T> void refineBasis_NS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal, bool prepatch = false);//, bool doubleChannel = false);
template<class T> void computeSteadyNS(uwbINSSolverSteady<T>& navStokesSteady, int numIterNSSteady, int plot_pts = 10000);

int main(int argc, char *argv[])
{
    // ========================================= Settings ========================================= 
    //---------------------------------------------------------------------------------------------
    gsInfo << "Reading file settings. \n";
    gsVector<int> inInt(14);
    gsVector<real_t> inRealT(12);
    gsVector<bool> inBool(15);

    //readInitialSetting(MOTOR_DATA_DIR "uwb-pilsen/initialSettingsINScircle_AFC.txt", inInt, inRealT, inString, inBool);
    readInitialSetting("initialSettingsINScircle_AFC.txt", inInt, inRealT, inBool);

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
    bool prepatch = inBool(11);
    bool plotMeshes = inBool(12);
    bool reduceContinuity = inBool(13);
    bool CrankNicholson = inBool(14);

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

    real_t reynoldsSteady = inRealT(0);
    real_t reynolds = inRealT(1);
    real_t timeStep = inRealT(2); // time step for unsteady computation
    real_t tol = inRealT(3); // stopping tolerance
    real_t uMax = inRealT(4); // inlet velocity maximum
    real_t SRBAValpha = inRealT(5);
    real_t radius = inRealT(6);
    real_t centreX = inRealT(7);
    real_t centreY = inRealT(8);
    real_t width = inRealT(9);
    real_t widthExpand = inRealT(10);
    real_t length = inRealT(11);
    //real_t SRBAVscaleFactorRes = inRealT(6);

    real_t viscositySteady = uMax*2*radius/reynoldsSteady;
    real_t viscosity = uMax*2*radius/reynolds;


    //bool doubleChannel = false;
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

    gsInfo << "Solving flow past cylinder NS example.\n";

    // ========================================= Define problem ========================================= 

    gsBoundaryConditions<> bcInfo;
    defineBCs_NS(bcInfo, uMax, prepatch);

    gsFunctionExpr<> f("0", "0", 2);

    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;

    //real_t radius = 0.05;
    //real_t centreX = 0.0;
    //real_t centreY = 0.0;
    //real_t width = 0.4;
    //real_t widthExpand = 0.01;
    //real_t length = 2.2;
    real_t prepatchWidth = width;

    patches = BSplineCircle2D4Patches<real_t>(deg, plotMeshes, length, width, widthExpand, radius, centreX, centreY, prepatch, prepatchWidth);

    gsInfo << patches << "\n";

    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> tbasis(patches);

    //tbasis.setDegree(deg);

    gsInfo << "tbasis.basis(0) = \n" << tbasis.basis(0) << "\n";
    
    refineBasis_NS(tbasis, numRefine, numRefineLocal, prepatch);//, doubleChannel);

    if (numElevate > 0)
    {
        tbasis.degreeElevate(numElevate);
        gsInfo << "\nBasis elevated.\n\n";
    }
    //gsInfo << "tbasis.basis(0) = \n" << tbasis.basis(0) << "\n";
    
    gsInfo << "-----------------------------------------------\n";
    gsInfo << "tbasis.basis(0) = \n" << tbasis.basis(0) << "\n";

    if (reduceContinuity)
    {
        //tbasis.reduceContinuity(1);
        //gsInfo << "Continuity reduced by 1. \n\ntbasis.basis(0) = \n" << tbasis.basis(0) << "\n";
        //gsInfo << "Continuity reduced by 1 in component 1. \n\n";
        for (unsigned p = 0; p < patches.nPatches(); p++)
        {
            for (int comp = 0; comp < 2; comp++)
            {
                int deg = tbasis.basis(p).component(comp).maxDegree();
                tbasis.basis(p).component(comp).reduceContinuity(deg-1);
            }
            gsInfo << "tbasis.basis(p) = \n" << tbasis.basis(p) << "\n";
        }
    }
    gsInfo << "-----------------------------------------------\n";
    
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    gsInfo << "\n==================================================\n";
    gsInfo << "Velocity basis max degree in component 0: " << discreteBases[0].maxDegree(0) << "\n";
    gsInfo << "Velocity basis max degree in component 1: " << discreteBases[0].maxDegree(1) << "\n";
    gsInfo << "Velocity basis min degree in component 0: " << discreteBases[0].minDegree(0) << "\n";
    gsInfo << "Velocity basis min degree in component 1: " << discreteBases[0].minDegree(1) << "\n";
    gsInfo << "Pressure basis max degree in component 0: " << discreteBases[1].maxDegree(0) << "\n";
    gsInfo << "Pressure basis max degree in component 1: " << discreteBases[1].maxDegree(1) << "\n";
    gsInfo << "Pressure basis min degree in component 0: " << discreteBases[1].minDegree(0) << "\n";
    gsInfo << "Pressure basis min degree in component 1: " << discreteBases[1].minDegree(1) << "\n";
    gsInfo << "==================================================\n\n";


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
        real_t diam = 2.*radius;
        real_t SRBAVscaleFactorRes = 1.;
        if (SRBAVscaleRes)
            SRBAVscaleFactorRes = diam / math::pow(uMax, 2);
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

    //--------------------------------------------------------------------------------------
    //sides at which is computed the lift
    /*std::vector<boxSide> circleSides;
    circleSides.push_back(boundary::east);
    circleSides.push_back(boundary::east);
    circleSides.push_back(boundary::east);
    circleSides.push_back(boundary::east);
    gsVector<int> circlePatchNumbers(4);
    circlePatchNumbers << 0, 1, 2, 3;

    gsVector<real_t> normalStream(2);
    normalStream << 0.0, 1.0;
    normalStream = normalStream / normalStream.norm();
    navStokes.setLiftParams(circlePatchNumbers, circleSides, uMax, 2*radius, normalStream);*/
    //--------------------------------------------------------------------------------------

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
    fd.save("circle_NSsolution.xml");

    // Optionally plot solution in paraview
    if (plot)
    {
        gsField<> velocity = navStokes.constructSolution(0);
        gsField<> pressure = navStokes.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, "circle_velocity", plot_pts);
        if (plotMeshes)
            gsWriteParaview<>(velocity, "circle_velocity_mesh", plot_pts, true);
        gsWriteParaview<>(pressure, "circle_pressure", plot_pts);

        gsMatrix<> sol = navStokes.getSolution();
        gsMatrix<> oldSol = navStokes.getAssembler()->getSolution();
        navStokes.getAssembler()->plotResiduum("reziduum", sol, oldSol, plot_pts);
    }

    //sides at which is computed the lift
    /*std::vector<boxSide> circleSides;
    circleSides.push_back(boundary::east);
    circleSides.push_back(boundary::east);
    circleSides.push_back(boundary::east);
    circleSides.push_back(boundary::east);
    gsVector<int> circlePatchNumbers(4);
    circlePatchNumbers << 0, 1, 2, 3;

    gsVector<real_t> normalStream(2);
    normalStream << 0.0, 1.0;
    normalStream = normalStream / normalStream.norm();
    navStokes.setLiftParams(circlePatchNumbers, circleSides, uMax, 2*radius, normalStream);*/

    //gsField<> pressureField = navStokes.constructSolution(1);
    //real_t cl = navStokes.getAssembler()->computeLift(pressureField);
    //gsInfo << "\nlift coefficient = " << cl << "\n";

    //navStokes.closeFile();

    return 0; 
}

//4 patches around circle
template<class TT>
gsMultiPatch<TT>  BSplineCircle2D4Patches(int degree, bool plotMeshes, const TT & length, TT const & width, TT const & widthExpand, TT const & radius, TT const & centreX, TT const & centreY, bool const &prepatch, TT const & prepatchWidth)
{
    const real_t PI = 3.141592653589793238463;
    //unsigned degree = 3;
    real_t alpha = PI/2.0;
    real_t const K = (4.0/3.0)*(tan(alpha/4.0)); //constant for circle
    //bool plotMeshes = false;

    gsMultiPatch<TT> mp;// = new gsMultiPatch<TT>;

    //initial data for patches
    gsKnotVector<TT> kv_cub(0, 1, 0, 4);
    gsKnotVector<TT> kv_lin(0, 1, 0, 2);
    gsTensorBSplineBasis<2, TT> basis(kv_lin, kv_cub);

    //--------------------------------patch 0-------------------------------------------
    gsMatrix<TT> coef_patch0(8, 2);
    coef_patch0 <<  - width/2, - width/2,
                    - radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY,
                     - width/2, - (width + widthExpand)/4,
                    - radius*cos(alpha/2) + centreX +  K * ( - radius * cos(alpha/2)), - radius*sin(alpha/2) + centreY +   K * (radius * sin(alpha/2)),
                      - width/2, (width + widthExpand)/4,
                    - radius*cos(alpha/2) + centreX +  K * ( - radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +   K * (- radius *  sin(alpha/2)),
                    - width/2,  width/2 + widthExpand,
                    - radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY;
    gsTensorBSpline<2, TT> patch0(basis, coef_patch0);
    patch0.degreeElevate(degree - 1, 0);


    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1(8, 2);
    coef_patch1 <<  - width/2, width/2 + widthExpand,
            - radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY,
                    - width/4, width/2 + widthExpand,
            - radius*cos(alpha/2) + centreX + K * ( radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +  K * ( radius * sin(alpha/2)),
                     width/4, width/2 + widthExpand,
                        radius*cos(alpha/2) + centreX + K * (- radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +  K * ( radius * sin(alpha/2)),
                     + width/2, width/2 + widthExpand,

                                radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY;
    gsTensorBSpline<2, TT> patch1(basis, coef_patch1);
    patch1.degreeElevate(degree - 1, 0);

    //--------------------------------patch 2-------------------------------------------
    gsMatrix<TT> coef_patch2(8, 2);
    coef_patch2 <<  - width/2,- width/2,
                   - radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY,
             - width/4, - width/2,
            - radius*cos(alpha/2) + centreX + K * ( radius * cos(alpha/2)),  - radius*sin(alpha/2) + centreY +  K * ( - radius * sin(alpha/2)),
             width/4, - width/2,
            radius*cos(alpha/2) + centreX + K * (- radius * cos(alpha/2)),  - radius*sin(alpha/2) + centreY +  K * ( - radius * sin(alpha/2)),
              + width/2,- width/2,
            radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY;

    gsTensorBSpline<2, TT> patch2(basis, coef_patch2);
    patch2.degreeElevate(degree - 1, 0);

    //--------------------------------patch 3-------------------------------------------
    gsMatrix<TT> coef_patch3(8, 2);
    coef_patch3 << width/2, - width/2,
                   radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY,
                      width/2, -(width + widthExpand)/4,
                    radius*cos(alpha/2) + centreX +  K * ( radius * cos(alpha/2)), - radius*sin(alpha/2) + centreY +   K * (radius*sin(alpha/2)),
                   width/2,  (width + widthExpand)/4,
                    radius*cos(alpha/2) + centreX +  K * ( radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +   K * (- radius * sin(alpha/2)),
                    width/2,  (width/2 + widthExpand),
                     radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY;
    gsTensorBSpline<2, TT> patch3(basis, coef_patch3);
    patch3.degreeElevate(degree - 1, 0);

    //--------------------------------patch 4-------------------------------------------
    gsMatrix<TT> coef_patch4(8, 2);
    coef_patch4 << length - width/2, - width/2,
                   width/2, - width/2,
                      length - width/2, - (width + widthExpand)/4,
                    width/2, - (width + widthExpand)/4,
                   length - width/2,  (width + widthExpand)/4,
                    width/2, + (width + widthExpand)/4,
                    length - width/2,  (width/2 + widthExpand),
                   width/2,  (width/2 + widthExpand);
    gsTensorBSpline<2, TT> patch4(basis, coef_patch4);
    patch4.degreeElevate(degree - 1, 0);

    //refining last patch according to length/width
    real_t insertKnotNum = trunc(length/width)-1;
    if(insertKnotNum>1.0){
        for(real_t i = 1/(2*insertKnotNum); i < 1.0; i+=1/(2*insertKnotNum)){
        patch4.insertKnot(i,0,1);
        }
    }


    mp.addPatch(patch0);
    mp.addPatch(patch1);
    mp.addPatch(patch2);
    mp.addPatch(patch3);
    mp.addPatch(patch4);

    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(0, boundary::south, 2, boundary::south);
    mp.addInterface(3, boundary::north, 1, boundary::north);
    mp.addInterface(3, boundary::south, 2, boundary::north);
    mp.addInterface(4, boundary::east, 3, boundary::west);

    mp.addInterface(2, boundary::west, 1, boundary::west);
    mp.addInterface(4, boundary::north, 4, boundary::south);


    //add start patch
    if(prepatch){
        //--------------------------------patch 5-------------------------------------------
        gsMatrix<TT> coef_patch5(8, 2);
        coef_patch5 << -prepatchWidth - width/2, - width/2,
                - width/2, - width/2,
                 - prepatchWidth - width/2, - (width + widthExpand)/4,
                - width/2 , - (width + widthExpand)/4,
                  - prepatchWidth - width/2, (width + widthExpand)/4,
                - width/2 , (width + widthExpand)/4,
                - prepatchWidth - width/2,  width/2 + widthExpand,
                - width/2,  width/2 + widthExpand;
        gsTensorBSpline<2, TT> patch5(basis, coef_patch5);
        patch5.degreeElevate(degree - 1, 0);

        real_t insertKnotNumPrepatch = trunc(prepatchWidth/width);

        if(insertKnotNumPrepatch > 0.00001){
            for(real_t i = 1/(2*insertKnotNumPrepatch); i < 1.0; i+=1/(2*insertKnotNumPrepatch)){

            patch5.insertKnot(i,0,1);
            }
        }


        mp.addPatch(patch5);
        mp.addInterface(0, boundary::west, 5, boundary::east);
        mp.addInterface(5, boundary::north, 5, boundary::south);

        if (plotMeshes){
            gsMesh<> mesh5;
            patch5.controlNet(mesh5);
            gsWriteParaview(mesh5,"cpCirclePatch5");
        }
    }




    mp.addAutoBoundaries();

    if (plotMeshes) {
        gsMesh<> mesh0;
        patch0.controlNet(mesh0);
        gsWriteParaview(mesh0,"cpCirclePatch0");

        gsMesh<> mesh1;
        patch1.controlNet(mesh1);
        gsWriteParaview(mesh1,"cpCirclePatch1");

        gsMesh<> mesh2;
        patch2.controlNet(mesh2);
        gsWriteParaview(mesh2,"cpCirclePatch2");

        gsMesh<> mesh3;
        patch3.controlNet(mesh3);
        gsWriteParaview(mesh3,"cpCirclePatch3");

        gsMesh<> mesh4;
        patch4.controlNet(mesh4);
        gsWriteParaview(mesh4,"cpCirclePatch4");


   }

    return mp;
}



//----------------
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
//----------------

template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, T uMax, bool prepatch)
{
    std::string xVel = util::to_string(uMax);// + " * (-4*(y-1.5)^2 + 1)";
    gsFunctionExpr<T> Uin(xVel, "0", 2); // inlet velocity
    gsFunctionExpr<T> Uwall("0", "0", 2); // wall velocity

    //----boundary conditions for 4 patches around circle----
    if(prepatch)
    {
        bcInfo.addCondition(5, boundary::west, condition_type::dirichlet, Uin, 0);
        //bcInfo.addCondition(5, boundary::north, condition_type::dirichlet, Uwall, 0);
        //bcInfo.addCondition(5, boundary::south, condition_type::dirichlet, Uwall, 0);
    }
    else
        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);

    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);
    //bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, Uwall, 0);
    //bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(3, boundary::east, condition_type::dirichlet, Uwall, 0);
    //bcInfo.addCondition(4, boundary::north, condition_type::dirichlet, Uwall, 0);
    //bcInfo.addCondition(4, boundary::south, condition_type::dirichlet, Uwall, 0);

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
    gsWriteParaview<>(velocitySteady, "circle_SteadyVelocity", plot_pts, true);
    gsWriteParaview<>(pressureSteady, "circle_SteadyPressure", plot_pts);

    gsFileData<real_t> fd;
    fd << navStokesSteady.getSolution();
    fd.save("circle_NSsteadySolution.xml");
}

//finerMesh

template<class T> void refineBasis_NS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal, bool prepatch)//, bool doubleChannel)
{
    gsMatrix<> box(2, 2);
    box << 0, 1, 0, 0;

    //int refP = 0;
    //if (doubleChannel)
    //    refP = 1;
    for (int i = 0; i < 2; i++)
        basis.refine(4, box);

    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    // refinement near wal
    int numRefineLocalWal = 3;

    real_t parArea = 0.0625;//125;
    real_t parArea2 = 0.15;

    for (int j = 0; j < numRefineLocalWal; j++)
    {
        box << 1-parArea, 1, 0, 0;
            basis.refine(0, box);
            basis.refine(1, box);
            basis.refine(2, box);
            basis.refine(3, box);

        parArea = parArea / 2;
    }
    //for (int j = 0; j < numRefineLocal-1; j++)
    //{
        //box << 0, 0, 0.5 - parArea2, 0.5+parArea2;
        box << 0, 0, 0, 1;
        //for (int i = 0; i < numRefineLocalWal; i++)
        //{
            basis.refine(4, box);
            basis.refine(3, box);
        //}
        //parArea2 = parArea2 / 2;
    //}

     box << 0, 0, 0.5 - parArea2, 0.5+parArea2;
     basis.refine(4, box);
     basis.refine(3, box);

     box << 0, 0, 0.8, 1;
     basis.refine(1, box);
     basis.refine(2, box);
}

/*//coarser mesh
template<class T> void refineBasis_NS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal, bool prepatch)//, bool doubleChannel)
{
    gsMatrix<> box(2, 2);
    box << 0, 1, 0, 0;

    //int refP = 0;
    //if (doubleChannel)
    //    refP = 1;
    for (int i = 0; i < 1; i++)
        basis.refine(4, box);

    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    // refinement near wal
    int numRefineLocalWal = 1;

    real_t parArea = 0.0625;//125;
    real_t parArea2 = 0.1;

    for (int j = 0; j < numRefineLocal; j++)
    {
        box << 1-parArea, 1, 0, 0;
        for (int i = 0; i < numRefineLocalWal; i++)
        {
            basis.refine(0, box);
            basis.refine(1, box);
            basis.refine(2, box);
            basis.refine(3, box);
        }

        parArea = parArea / 2;
    }
    for (int j = 0; j < numRefineLocal-1; j++)
    {
        box << 0, 0, 0.5 - parArea2, 0.5+parArea2;
        for (int i = 0; i < numRefineLocalWal; i++)
        {
            basis.refine(4, box);
            basis.refine(3, box);
        }
        parArea2 = parArea2 / 2;
    }
}*/
