/** @file uwbRANSTestStep2D4P.cpp

Author(s): E. Turnerova
*/
#ifdef _MSC_VER // to be removed
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>



// solvers
#include "uwbINSSolverSteady.h"
#include "uwbRANSSolver.h"
#include "uwbRANSSolverIterative.h"
#include "uwbTMSolverKOmegaLinSteady.h"
#include "uwbTMSolverKOmega.h"

//problem settings
#include "uwbRANSExamplesSetting.h"
#include "uwbGeometryCreators.h"
#include "uwbDataFileReader.h"

#include <math.h>
#include <string.h>
#include <sstream>


using namespace gismo;

int main(int argc, char *argv[])
{
    //========================================== load data from file =======================================================
    gsInfo << "Reading file settings. \n";
    gsVector<int> inInt(27);
    gsVector<real_t> inRealT(12);
    std::string inString;
    gsVector<bool> inBool(22);

    //readInitialSetting(MOTOR_DATA_DIR "uwb-pilsen/initialSettingsRANS2D4Pstep.txt", inInt, inRealT, inString, inBool);
    readInitialSetting("initialSettingsRANS2D4Pstep.txt", inInt, inRealT, inString, inBool);

    // ========================================= Settings =========================================
    //====================================================
    bool plot = inBool(0);
    bool plotMeshes = inBool(1);

    bool SUPG = inBool(2);
    bool TCSD = inBool(3);
    bool RANSad = inBool(4);
    //...................
    bool RANSsrbav = inBool(5);
    //...................

    bool TMsupg = inBool(6);
    bool TMad = inBool(7);
    //...................
    bool TMsrbav = inBool(8);
    //...................
    //bool TMafc = inBool(9);
    //bool TMafcHO = inBool(10);

    bool animate = inBool(11);

    bool limitTMProduction = inBool(12);

    bool loadIC = inBool(13);
    bool NSsteadyInitial = inBool(14);

    bool computeTMfirst = inBool(15);

    bool bOFbc = inBool(16);

    bool reSaveSolution = inBool(17);

    bool dirElemLengthRANS = inBool(18);
    bool dirElemLengthTM = inBool(19);

    //bool RANSafc = inBool(20);
    //bool RANSafcHO = inBool(21);

    int plot_pts = inInt(0);

    int numIterSteadyNS = inInt(1);
    int numIterKOmegaSteady = inInt(2);
    int numIterRANS = inInt(3);
    int minNumIterRANS = inInt(4);
    int maxRANSPicardIt = inInt(5);
    int maxTMPicardFirstIt = inInt(6);
    int maxTMPicardIt = inInt(7);

    int tauSUPGType = inInt(8);
    int RANSsrbavType = inInt(9);
    int RANSsrbavResidualType = inInt(10);
    int tauSRBAV = inInt(11);
    int tauStabTypeAD = inInt(12);

    int tauTMStabTypeSUPG = inInt(13);

    int TMsrbavType = inInt(14);
    int TMsrbavResidualType = inInt(15);
    int TMtauSRBAV = inInt(16);

    int animateStep = inInt(17);

    int refineType = inInt(18);//2;
    int numRefine = inInt(19);//3;
    int numRefineLocal = inInt(20);//1;

    int deg = inInt(21);

    int plot_pts_atWalls = inInt(22);

    int numFinalYplusRef = inInt(23);
    int numFinalSingularPointRef = inInt(24);

    int hDirTypeRANS = inInt(25);
    int hDirTypeTM = inInt(26);

    real_t timeStep = inRealT(0);
    real_t tolRelNorm = inRealT(1);

    real_t viscosity = inRealT(2);
    real_t viscositySteady = inRealT(3);

    real_t turbIntensity = inRealT(4);
    real_t viscosityRatio = inRealT(5);

    real_t kIn = inRealT(6);
    real_t oIn = inRealT(7);

    real_t productionXPoint = inRealT(8);

    real_t uMax = inRealT(9);

    real_t RANSsrbavAlpha = inRealT(10);
    real_t TMsrbavAlpha = inRealT(11);

    std::string tmEvaluator = inString;

    //-----------
    bool printInfo = true;
    if (reSaveSolution)
    {
        plot = true;
        loadIC = true;
        numIterRANS = 0;
        plotMeshes = false;
        animate = false;
        printInfo = false;
    }

    // ========================================= Geometry =========================================
    gsInfo << "Reading file geometry. \n";

    gsVector<int> inGeoInt(7);
    gsMatrix<int> refineUniformSettings(0,0);
    gsMatrix<int> refineLocalSettings(0,0);
    gsVector<real_t> geomParams(0);
    gsVector<bool> inGeoBool(2);
    std::vector<real_t> kvfit_knots;

    //readInitialGeometry(MOTOR_DATA_DIR "uwb-pilsen/initialGeometryStep.txt", inGeoInt, refineUniformSettings, refineLocalSettings, geomParams, inGeoBool, kvfit_knots);
    readInitialGeometry("initialGeometryStep.txt", inGeoInt, refineUniformSettings, refineLocalSettings, geomParams, inGeoBool, kvfit_knots);

    if (!reSaveSolution)
    {
        gsInfo << "inGeoInt = \n" << inGeoInt << "\n";
        gsInfo << "inGeoBool = \n" << inGeoBool << "\n";
    }

    int uniformRefine = inGeoInt(3);
    //========================================================================================================================

    //------------direct solver----------------------
    typedef uwbRANSSolver<real_t> SolverType;

    //------------linear solvers---------------------
    bool iterative = false;
    //typedef uwbRANSSolverIterative<real_t, uwbGMResRight<real_t> > SolverType;

    // preconditioner type options:
    // LSC_AdiagEqual, LSC_Adiag, LSC_Awhole, LSC_Amod
    // AL_Awhole, AL_Amod
    // SIMPLE_AdiagEqual, SIMPLE_Adiag, SIMPLE_Awhole, SIMPLE_Amod
    // SIMPLER_AdiagEqual, SIMPLER_Adiag, SIMPLER_Awhole, SIMPLER_Amod
    // MSIMPLER_AdiagEqual, MSIMPLER_Adiag, MSIMPLER_Awhole, MSIMPLER_Amod
    //std::string precType = "LSC_AdiagEqual";
    std::string precType = "MSIMPLER_AdiagEqual";

    real_t gamma = 2.5; // parameter for AL precond.

    int linMaxIt = 100;
    real_t linTol = 1e-7;

    int numThreads = 1; // number of threads for assembly

    if (!reSaveSolution)
        gsInfo << "Solving turbulent flow for backward facing step.\n";
    else
        gsInfo << "\nResave solution from loaded data.\n\n";

    //==================================================================================================
    //==================================================================================================
    //==================================================================================================
    //==================================================================================================
    // ========================================= Define problem =========================================
    gsInfo << "===================================================\n";
    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.intStrategy = iFace::glue;
    opt.dirValues = dirichlet::interpolation;
    //opt.dirValues = dirichlet::l2Projection; //default

    real_t stepHeight = 0.0127;
    real_t preStepHeight = 8 * stepHeight;
    real_t preStepLength = 1.397;//4 * stepHeight;
    real_t prePatchLength = 0.254;//0.401;
    real_t length = 50 * stepHeight;

    uwbRANSBackwardStep2D4PExample<real_t> problemSettings(viscosity, stepHeight, preStepHeight,
                                                           preStepLength, prePatchLength, length,
                                                           deg, uMax, plot_pts);

    //T viscosity, T stepHeight, T preStepHeight, T preStepLength, T prePatchLength, T length, T uMax, int plot_pts
    //========================================= Define geometry =========================================
    gsBoundaryConditions<> bcInfo;
    problemSettings.defineBCs_NS(bcInfo);
    gsFunctionExpr<> f("0", "0", 2);

    gsMultiPatch<real_t> patches = problemSettings.makeMultiPatch();
    if (!reSaveSolution)
        gsInfo << patches << "\n";

    // ========================================= Define basis and refine =========================================
    gsMultiBasis<> tbasis(patches); // basis for RANS equations
    //real_t wallRefineKnot = 0.15;
    switch (refineType) {
      case 1:
        numRefine = 2;
        numRefineLocal = 6;

        problemSettings.refineBasis(tbasis, numRefine, numRefineLocal);

        break;

    case 2:
      problemSettings.refineBasis(tbasis, numRefine, numRefineLocal, numFinalYplusRef, numFinalSingularPointRef);
      //for (int i = 0; i < numRefine; ++i)
      //    tbasis.uniformRefine();

      break;

    case 3:
      for(int i = 0; i < uniformRefine; i++)
          tbasis.uniformRefine();

      //problemSettings.refineBasisType2(tbasis, numRefine);

      for(int i = 0; i < refineUniformSettings.rows(); i++)
          problemSettings.refineBasisUniformZones(tbasis, refineUniformSettings(i,0), refineUniformSettings(i,1), refineUniformSettings(i,2), refineUniformSettings(i,3), refineUniformSettings(i,4));

      for(int i = 0; i < refineLocalSettings.rows(); i++)
          problemSettings.refineBasisLocalZones(tbasis, refineLocalSettings(i,0), refineLocalSettings(i,1), refineLocalSettings(i,2), refineLocalSettings(i,3), refineLocalSettings(i,4));

      break;
    default: //default uniform refinement

      //numRefine = 3;
      //uniform refinement:
      for (int i = 0; i < numRefine; ++i)
          tbasis.uniformRefine();
      if (refineType != 2)
          break;
  }

    int numElement = 0;
    for (int p = 0; p < 3; p++)
        numElement += tbasis.basis(p).numElements();
    gsInfo << "\nnumElements = " << numElement << "\n\n";

    /*const gsTensorBSplineBasis<2, real_t>*  basis0 = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&tbasis.piece(0));
    //const gsKnotVector<>& knotVectorY = basis->knots(1);
    gsInfo << "basis0->knots(1) = \n" << basis0->knots(1) << "\n";
    const gsTensorBSplineBasis<2, real_t>*  basis2 = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&tbasis.piece(2));
    //const gsKnotVector<>& knotVectorY = basis->knots(1);
    gsInfo << "basis2->knots(1) = \n" << basis2->knots(1) << "\n";
    getchar();*/

    //-------------------------------------------------------------------------------------------------
    //==================================================================================================

    std::ofstream file;
    file.open("Step2DOutputInfo.txt");
    if (printInfo)
    {
        file << "============parameters===================" << "\n";
        file << "timeStep: " << timeStep << "\n";
        file << "viscosity: " << viscosity << "\n";
        file << "viscositySteady: " << viscositySteady << "\n";
        file << "turbIntensity: " << turbIntensity << "\n";
        file << "numIterSteadyNS : " << numIterSteadyNS << "\n";
        file << "numIterKOmegaSteady : " << numIterKOmegaSteady << "\n";
        file << "minNumIterRANS : " << minNumIterRANS << "\n";
        file << "maxRANSPicardIt : " << maxRANSPicardIt << "\n";
        file << "tolRelNorm : " << tolRelNorm << "\n";
        file << "tmEvaluator = " << tmEvaluator << "\n";
        file << "maxTMPicardFirstIt : " << maxTMPicardFirstIt << "\n";
        file << "iterative : " << iterative << "\n";
        file << "lin. max iter. : " << linMaxIt << "\n";
        file << "lin. tol : " << linTol << "\n";
        file << "numRefine : " << numRefine << "\n";
        file << "numRefineLocal : " << numRefineLocal << "\n";
        file << "=======================================\n";
    }

    int num_meshes;
    num_meshes = patches.nPatches();

    if (plotMeshes)
    {
        for (int index_of_mesh = 0; index_of_mesh < num_meshes; index_of_mesh++)
        {
            gsMesh<> mesh;
            std::ostringstream strs_mesh;
            strs_mesh << index_of_mesh;
            makeMesh(tbasis.at(index_of_mesh), mesh, 10);
            patches.patch(index_of_mesh).evaluateMesh(mesh);
            gsWriteParaview(mesh, "step_meshPatch" + strs_mesh.str());
        }
        gsInfo << "\nMesh plotted.\n\n";
    }

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    if (!reSaveSolution)
        gsInfo << "Velocity basis degree: " << discreteBases[0].degree() << "\n\n";
    if (printInfo)
        file << "Velocity basis degree: " << discreteBases[0].degree() << "\n";

    // ========================================= Compute Steady NS =========================================
    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscositySteady);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    params.setNumThreads(numThreads); // set the number of threads for assembly

    /*if (SUPG)
    {
        params.settings().set(constantsINS::SUPG, SUPG);
        params.settings().set(constantsINS::tauStabType, tauSUPGType);
    }*/

    uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver

    if (!reSaveSolution)
        gsInfo << "numDofs = " << navStokes.numDofs() << "\n";

    if (!loadIC && NSsteadyInitial)
        problemSettings.computeSteadyNS(navStokes, numIterSteadyNS);

    //gsMatrix<> sol = navStokes.getSolution();
    //navStokes.getAssembler()->plotResiduum("residuumSteadyNS", sol, plot_pts);

    //================================= wall distance estimation ==================================================
    real_t wallDistance = problemSettings.computeWallDistance(navStokes);

    if (!reSaveSolution)
    {
       gsInfo << "\nminimum wallDistance = " << wallDistance << "\n";
       if (printInfo)
           file << "\nminimum wallDistance = " << wallDistance << "\n";

    //==================================== compute aspect ratio  ==================================================
        real_t maxAspectRatio = navStokes.computeAspectRatio();
        gsInfo << "maxAspectRatio = " << maxAspectRatio << "\n";
        real_t minAspectRatio = navStokes.computeAspectRatio(true);
        gsInfo << "minAspectRatio = " << minAspectRatio << "\n";
        if (printInfo)
        {
            file << "maxAspectRatio = " << maxAspectRatio << "\n";
            file << "minAspectRatio = " << minAspectRatio << "\n";
        }

    }
    // ========================================= Define turbulence solver =========================================
    gsBoundaryConditions<> bcInfoTurb;

    real_t uFreeStream = problemSettings.getFreeStream();
    real_t kInConst, oInConst;
    if (bOFbc)
    {
        kInConst = kIn;//0.00109;
        oInConst = oIn;//181728;
    }
    else
    {
        kInConst = 1.5 * math::pow(uFreeStream * turbIntensity, 2); // (3/2)*(UI)^2
        oInConst = kInConst / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at inlet
    }

    real_t beta = 0.0708;
    real_t kWall, oWall;
    /*if (bOFbc)
    {
        kWall = kInConst;
        oWall = oInConst;
    }
    else
    {*/
        kWall = 0.;
        oWall = 6 * viscosity / (beta * math::pow(wallDistance, 2));
        oWall = oWall / 10.;
    //}

    if (!reSaveSolution)
    {
        gsInfo << "\nkInConst = " << kInConst << "\n";
        gsInfo << "oInConst = " << oInConst << "\n";
        gsInfo << "kWall = " << kWall << "\n";
        gsInfo << "oWall = " << oWall << "\n\n";
    }
    if (printInfo)
    {
        file << "\nkInConst = " << kInConst << "\n";
        file << "oInConst = " << oInConst << "\n";
        file << "kWall = " << kWall << "\n";
        file << "oWall = " << oWall << "\n\n";
    }

    problemSettings.defineBCs_TM(bcInfoTurb, kInConst, kWall, oInConst, oWall, bOFbc);

    gsDofMapper koMapper;
    discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfoTurb, koMapper, 0);

    uwbINSPde<real_t> koPdeSteady(patches, bcInfoTurb, f, viscositySteady);
    uwbINSSolverParams<real_t> koParamsSteady(koPdeSteady, discreteBases, opt);

    koParamsSteady.setNumThreads(numThreads); // set the number of threads for assembly

    //if (TMsupg)
    //    koParamsSteady.settings().set(constantsINS::TMsupg, TMsupg);

    gsMatrix<> koInitial(koMapper.freeSize(), 2);
    koInitial.col(0).setConstant(kInConst);//(kWall);
    koInitial.col(1).setConstant(oInConst);//(oWall);

    uwbTMSolverKOmegaLinSteady<real_t> turbSolver(koParamsSteady);
    turbSolver.setInitialCondition(koInitial);
    // ========================================= Solving k-omega =========================================

    if (!loadIC && NSsteadyInitial)
        problemSettings.computeSteadyKOmegaTM(navStokes, turbSolver, numIterKOmegaSteady, koMapper.freeSize());

    // ========================================= Solving RANS =========================================
    if (!reSaveSolution)
    {
        gsInfo << "\nSolving RANS...\n";

        if (iterative)
        {
            gsInfo << "SOLVING with iterative solver, preconditioner type: " << precType << "\n\n";
            file << "SOLVING with iterative solver, preconditioner type: " << precType << "\n\n";

            if (precType == "AL_Awhole" || precType == "AL_Amod")
                file << "gamma: " << gamma << "\n\n";
        }
        else
        {
            gsInfo << "SOLVING with direct solver.\n\n";
            file << "SOLVING with direct solver.\n\n";
        }
    }

    uwbINSPde<real_t> RANSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> RANSparams(RANSpde, discreteBases, opt);

    RANSparams.setNumThreads(numThreads); // set the number of threads for assembly

    RANSparams.settings().set(constantsINS::timeStep, timeStep);
    RANSparams.settings().set(constantsINS::unst_innerIt, maxRANSPicardIt);
    //RANSparams.settings().set(constantsINS::unst_innerTol, picardTol);

    RANSparams.settings().set(constantsINS::dirElemLength, dirElemLengthRANS);
    RANSparams.settings().set(constantsINS::hDirType, hDirTypeRANS);

    if (SUPG)
    {
        RANSparams.settings().set(constantsINS::SUPG, SUPG);
        RANSparams.settings().set(constantsINS::tauStabTypeSUPG, tauSUPGType);
    }
    if (TCSD)
    {
        RANSparams.settings().set(constantsINS::TCSD, TCSD);
        RANSparams.settings().set(constantsINS::tauStabTypeSUPG, tauSUPGType);
    }

    if (RANSad)
    {
        RANSparams.settings().set(constantsINS::RANSad, RANSad);
        RANSparams.settings().set(constantsINS::tauStabTypeAD, tauStabTypeAD);
    }
    if (RANSsrbav)
    {
        RANSparams.settings().set(constantsINS::SRBAV, RANSsrbav);
        RANSparams.settings().set(constantsINS::SRBAVtype, RANSsrbavType);
        RANSparams.settings().set(constantsINS::tauStabTypeSRBAV, tauSRBAV);
        RANSparams.settings().set(constantsINS::SRBAVresidualType, RANSsrbavResidualType);
        RANSparams.settings().set(constantsINS::SRBAValpha, RANSsrbavAlpha);
    }
    //if (RANSafc)
    //    RANSparams.settings().set(constantsINS::AFC, RANSafc);
    //else if (RANSafcHO)
    //    RANSparams.settings().set(constantsINS::AFC_HO, RANSafcHO);

    //if (TMafc || TMafcHO)
    //    discreteBases[1].getMapper(dirichlet::none, opt.intStrategy, bcInfoTurb, koMapper, 0);
    //else
    //    discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfoTurb, koMapper, 0);

    uwbINSPde<real_t> koPde(patches, bcInfoTurb, f, viscosity);
    uwbINSSolverParams<real_t> koParams(koPde, discreteBases, opt);

    koParams.setNumThreads(numThreads); // set the number of threads for assembly

    koParams.settings().set(constantsINS::timeStep, timeStep);
    koParams.settings().set(constantsINS::turb_innerFirstIt, maxTMPicardFirstIt);
    koParams.settings().set(constantsINS::turb_innerIt, maxTMPicardIt);
    koParams.settings().setTurbulenceEvaluator(tmEvaluator);

    koParams.settings().set(constantsINS::dirElemLength, dirElemLengthTM);
    koParams.settings().set(constantsINS::hDirType, hDirTypeTM);

    if (TMsupg)
    {
        koParams.settings().set(constantsINS::TMsupg, TMsupg);
        koParams.settings().set(constantsINS::tauStabTypeSUPG, tauTMStabTypeSUPG); //set formula for stabilization parameter tau
    }
    //else if (TMafc)
    //    koParams.settings().set(constantsINS::TMafc, TMafc);
    //else if (TMafcHO)
    //    koParams.settings().set(constantsINS::TMafcHO, TMafcHO);
    else if (TMad)
    {
        koParams.settings().set(constantsINS::TMad, TMad);
        koParams.settings().set(constantsINS::tauStabTypeAD, tauTMStabTypeSUPG);
    }
    else if (TMsrbav)
    {
        koParams.settings().set(constantsINS::SRBAV, TMsrbav);
        koParams.settings().set(constantsINS::SRBAVtype, TMsrbavType);
        koParams.settings().set(constantsINS::tauStabTypeSRBAV, TMtauSRBAV);
        koParams.settings().set(constantsINS::SRBAVresidualType, TMsrbavResidualType);
        koParams.settings().set(constantsINS::SRBAValpha, TMsrbavAlpha);
    }

    if (limitTMProduction)
    {
        koParams.settings().set(constantsINS::limitTMProduction, limitTMProduction);
        koParams.settings().set(constantsINS::productionXPoint, productionXPoint);
    }

    gsMatrix<real_t> RANSSolution, TMSolution;
    //------------- read from file ---------------------
    if (loadIC)
    {
        gsFileData<> fdRead("RANS_IC_solution.xml");
        RANSSolution = *(fdRead.getFirst< gsMatrix<real_t> >());

        gsFileData<> fdReadTM("TM_IC_solution.xml");
        TMSolution = *(fdReadTM.getFirst< gsMatrix<real_t> >());
    }
    //------------------------------------------------

    uwbTMSolverKOmega<real_t> turbSolver_unsteady(koParams);

    gsMatrix<> koInitialUnsteady(koMapper.freeSize(), 2);
    koInitialUnsteady.col(0).setConstant(kInConst);
    koInitialUnsteady.col(1).setConstant(oInConst);
    if (loadIC)
        turbSolver_unsteady.setInitialCondition(TMSolution);
    else if (NSsteadyInitial)
        turbSolver_unsteady.setInitialCondition(turbSolver.getSolution());
    else
        turbSolver_unsteady.setInitialCondition(koInitialUnsteady);//(turbSolver.getSolution());//

    if (tmEvaluator != "koWilcoxLRN")
        problemSettings.solvePoissonEquation(turbSolver_unsteady);

    if (iterative)
    {
        RANSparams.settings().setPrecondType(precType);
        RANSparams.getPrecOptions().setReal("gamma", gamma);
        RANSparams.settings().set(constantsINS::iter_maxIt, linMaxIt);
        RANSparams.settings().set(constantsINS::iter_tol, linTol);
    }
    RANSparams.settings().setTurbulenceSolver(&turbSolver_unsteady);

    //--------
    //uwbRANSSolver<real_t> ransSolver(RANSparams);
    //uwbRANSSolverIterative<real_t, uwbGMResRight<real_t> > ransSolver(RANSparams);
    SolverType ransSolver(RANSparams);
    //--------

    if (!reSaveSolution)
        gsInfo << "initialization...\n";
    ransSolver.setComputationSequence(computeTMfirst);

    if (loadIC)
    {
        //if (RANSafc || RANSafcHO)
        //    RANSSolution = ransSolver.getSolution_full(RANSSolution);
        ransSolver.setInitialCondition(RANSSolution);
    }
    else if (NSsteadyInitial)
    {
        gsMatrix<> steadySol = navStokes.getSolution();
        //if (RANSafc || RANSafcHO)
        //   steadySol = ransSolver.getSolution_full(steadySol);
        ransSolver.setInitialCondition(steadySol);
    }
    else
    {
        gsDofMapper uMapper, pMapper;
        /*if (RANSafc || RANSafcHO)
        {
            discreteBases[0].getMapper(dirichlet::none, opt.intStrategy, bcInfo, uMapper, 0);
            discreteBases[1].getMapper(dirichlet::none, opt.intStrategy, bcInfo, pMapper, 1);
        }
        else
        {
            discreteBases[0].getMapper(opt.dirStrategy, opt.intStrategy, bcInfo, uMapper, 0);
            discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfo, pMapper, 1);
        }*/
        discreteBases[0].getMapper(opt.dirStrategy, opt.intStrategy, bcInfo, uMapper, 0);
        discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfo, pMapper, 1);
        int uDofs = uMapper.freeSize();
        int pDofs = pMapper.freeSize();
        gsVector<> constV1, constV2;
        constV1.setConstant(uDofs, 1e-8);
        constV2.setConstant(uDofs + pDofs, 0.);
        gsMatrix<> upInitial(2*uDofs + pDofs, 1);
        upInitial.middleRows(0, uDofs) = constV1;
        upInitial.middleRows(uDofs, uDofs + pDofs) = constV2;

        ransSolver.setInitialCondition(upInitial);
    }

    ransSolver.initialize();
    /*if (loadIC)
        ransSolver.setInitialCondition(RANSSolution);
    else
        ransSolver.setStokesInitialCondition();//.setInitialCondition(navStokes.getSolution());*/

    if (!reSaveSolution)
        gsInfo << "numDofs: " << ransSolver.numDofs() << "\n";
    if (printInfo)
        file << "numDofs: " << ransSolver.numDofs() << "\n";

    //ransSolver.getAssembler()->plotRANSresidual("rezidual", ransSolver.getSolution(), plot_pts);

    //=== solving with solution change relative norm as stopping criterion ============
    if (!reSaveSolution)
    {

        if (animate)
            ransSolver.solveWithAnimation(numIterRANS, animateStep, tolRelNorm, plot_pts, true);
        else
            ransSolver.solve(numIterRANS, tolRelNorm); // solution change norm tol = 10^(-5)

        real_t Tassembly = ransSolver.getAssemblyTime();
        real_t Tsetupsolve = ransSolver.getSolverSetupTime();
        real_t Tsolve = ransSolver.getSolveTime();
        real_t Tturbassembly = ransSolver.getTurbAssemblyTime();
        real_t Tturbsetupsolve = ransSolver.getTurbSolverSetupTime();
        real_t Tturbsolve = ransSolver.getTurbSolveTime();

        gsInfo << "Assembly time:" << Tassembly << "\n";
        gsInfo << "Solve time:" << Tsolve << "\n";
        gsInfo << "Solver setup time:" << Tsetupsolve << "\n";
        gsInfo << "Turbulent assembly time:" << Tturbassembly << "\n";
        gsInfo << "Turbulent solver setup time:" << Tturbsetupsolve << "\n";
        gsInfo << "Turbulent solve time:" << Tturbsolve << "\n";

        if (printInfo)
        {
            file << "Assembly time:" << Tassembly << "\n";
            file << "Solve time:" << Tsolve << "\n";
            file << "Solver setup time:" << Tsetupsolve << "\n";
            file << "Turbulent assembly time:" << Tturbassembly << "\n";
            file << "Turbulent solver setup time:" << Tturbsetupsolve << "\n";
            file << "Turbulent solve time:" << Tturbsolve << "\n";
            file.close();
        }

        //--- save solution into file ---
        gsFileData<> fd, fd_TM;
        fd << ransSolver.getSolution();
        fd_TM << turbSolver_unsteady.getSolution();
        fd.save(problemSettings.getGeometryName() + "RANS_solution.xml");
        fd_TM.save(problemSettings.getGeometryName() + "TM_solution.xml");

    }
    //=====================================================================================================

    // Optionally plot solution in paraview
    if (plot)
    {
        gsInfo << "Plotting in Paraview...\n";

        gsField<> velocity = ransSolver.constructSolution(0);
        gsField<> pressure = ransSolver.constructSolution(1);
        gsField<> kOmega = turbSolver_unsteady.constructSolution();

        problemSettings.plotSolutionField(velocity, "velocity");
        problemSettings.plotSolutionField(pressure, "pressure");
        problemSettings.plotSolutionField(kOmega, "komega");
        ransSolver.plotTurbulentViscosity(problemSettings.getGeometryName() + "turb_viscosity", plot_pts);

        const gsGeometry<> * geoResult = & velocity.igaFunction();
        gsWriteParaview(* geoResult,"geoResult",10000,true,true);

        gsVector<> referencePoint(2);
        referencePoint << -0.0508, 0.0508;//-preStepLength, preStepHeight/2.;
        //referencePoint << 0, 0.5;//-preStepLength, preStepHeight/2.;
        real_t density = 1000.;
        ransSolver.plotPressureCoefficient("pressureCoefficient_Cp", 2, referencePoint, uMax, density, plot_pts_atWalls);
        ransSolver.plot2DVorticity("vorticity", plot_pts);

        gsMatrix<> ransSol = ransSolver.getSolution();
        gsMatrix<> ransOldSol = ransSolver.getAssembler()->getSolution();
        gsMatrix<> tmSol = turbSolver_unsteady.getSolution();
        ransSolver.plotResiduum("RANSresiduum", ransSol, ransOldSol, tmSol, plot_pts);

        //ransSolver.plotShearStress("shearStress", ransSol, plot_pts);
        //ransSolver.saveCfAtWall("CfAtWall", ransSol, 0, referencePoint, plot_pts_atWalls);

        std::vector<boxSide> sides;
        sides.push_back(boundary::south);
        sides.push_back(boundary::south);
        gsVector<int> patchNumbers(2);
        patchNumbers << 0, 2;
        //ransSolver.plot2DSkinFrictionCoefficient("skinFrictionCoefficient_Cf", patchNumbers, sides, uMax, density, plot_pts_atWalls);//math::round(math::sqrt(plot_pts)));
    }

    return 0;
}
