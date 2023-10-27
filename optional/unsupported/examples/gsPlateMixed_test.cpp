/**  gsPlateMixed_test.cpp
    This file is a test of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): K. Rafetseder
*/


#include <gismo.h>
#include <gsPlateMixed/gsPlateMixedAssembler.h>
#include <gsPlateMixed/gsNormL2M.h>
#include <gsPde/gsShellMixedPde.h>
#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsSeminormH1.h>


using namespace gismo;
using std::cout;
using std::endl;


int main (int argc, char** args)
{
    int result = 0;

    /// Options ///

    int testcase_geometry = 1;
    // testcase_geometry: 1 rectangle [-1,1]*[-1,1] (1 patch)
    //                    2 custom geometry (1 patch)
    //                    3 rectangle [-1,1]*[-1,1] (4 patches)
    //                    4 L-shape p=1 (2 patches)
    int testcase_bc = 5;
    // testcase_bc: 1 whole boundary clamped
    //              2 whole boundary simplysupp
    //              3 two opposite edges simplysupp others free
    //              4 two opposite edges simplysupp others clamped
    //              5 two opposite edges simplysupp others free and clamped

    //              6 constant load f = 1

    //              7 disc whole boundary simplysupp, f=1

    int nonHomogeneous_bc = false;


    bool plot = false;
    bool analytSol = true;
    bool convRates_p_phi = true;
    bool relative_error = false;
    bool computeEValsA11 = false;

    bool outputEvals = false;
    bool outputMatrix = false;
    bool outputBlock = false;
    
    int numIncreaseDegree = 0;

    int numHRefine_start = 3;// 3
    int numHRefine_end  = 5; // 7
    int fineHRefine_add = 1; // 2
    int sizeError = numHRefine_end - numHRefine_start+1;
    bool keepC0 = false;

    gsOptionList options = gsPlateMixedAssembler<>::defaultOptions();
    //options.setInt("DirichletValues"  , dirichlet::homogeneous);
    options.setInt("DirichletValues"  , dirichlet::interpolation);
    //options.setInt("DirichletValues"  , dirichlet::l2Projection);
    options.setInt("DirichletStrategy", dirichlet::elimination );
    options.setInt("InterfaceStrategy", iFace    ::conforming  );
    options.setReal("quA", 1.0);
    options.setInt ("quB", 1  );
    options.setReal("bdA", 2.0);
    options.setInt ("bdB", 1  );


    options.setInt("phiBcMethodSs", gsPlateMixedAssembler<>::nitsche);
    options.setInt("phiBcMethodF", gsPlateMixedAssembler<>::nitsche);
    options.setSwitch("wHomogenPsiq", true);
    options.setSwitch("wPenaltyPsiq", true);
    options.setSwitch("solveWholeSystem",  false);
    options.setSwitch("pBoundary0", false);
    options.setSwitch("lambdaTMean0", false);



    /*
    options.setInt("phiBcMethodSs", gsPlateMixedAssembler<>::nitscheDerivative);
    options.setInt("phiBcMethodF", gsPlateMixedAssembler<>::nitscheDerivative);
    options.setSwitch("wHomogenPsiq", true);
    options.setSwitch("wPenaltyPsiq", false);
    options.setSwitch("solveWholeSystem",  false);
    options.setSwitch("pBoundary0", false);
*/


/*
    options.setInt("phiBcMethodSs", gsPlateMixedAssembler<>::lagrange);
    options.setInt("phiBcMethodF", gsPlateMixedAssembler<>::lagrange);
    options.setSwitch("wHomogenPsiq", false);
    options.setSwitch("wPenaltyPsiq", false);
    options.setSwitch("solveWholeSystem", false);
    options.setSwitch("pBoundary0", false);
    options.setSwitch("lambdaTMean0", false);
*/




    // Options for fine solution
    gsOptionList optionsFine(options);

/*
    optionsFine.setInt("phiBcMethodSs", gsPlateMixedAssembler<>::nitsche);
    optionsFine.setInt("phiBcMethodF", gsPlateMixedAssembler<>::nitsche);
    optionsFine.setSwitch("wHomogenPsiq", true);
    optionsFine.setSwitch("wPenaltyPsiq", true);
    optionsFine.setSwitch("solveWholeSystem", false);
    optionsFine.setSwitch("pBoundary0", false);
*/

    /*
    optionsFine.setInt("phiBcMethodSs", gsPlateMixedAssembler<>::lagrange);
    optionsFine.setInt("phiBcMethodF", gsPlateMixedAssembler<>::lagrange);
    optionsFine.setSwitch("wHomogenPsiq", false);
    optionsFine.setSwitch("wPenaltyPsiq", false);
    optionsFine.setSwitch("solveWholeSystem", false);
    optionsFine.setSwitch("pBoundary0", false);
*/



    //gsInfo << "Assembler "<< options;


    /// Variables ///

    gsMatrix<real_t> data;
    gsMatrix<real_t> data_timing;
    data.setZero(sizeError, 9);
    data_timing.setZero(sizeError, 9);

    int indRefine = 0;

    gsMultiPatch<> sol_w, fineSol_w, sol_p, fineSol_p, sol_phi, fineSol_phi, sol_M, fineSol_M;
    gsMatrix<> solVector_p, solVector_phi, solVector_w, solVector_lambdaN, solVector_lambdaT, solVector_lambdaMean, solVector_lambda, solVector_phi_lambda, solVector;

    gsStopwatch watch;


    /// Right-hand side, analytical solution ///

    gsFunctionExpr<> f,analytSol_w, analytSol_M;
    gsFunctionExpr<> dataClamped; // grad w
    gsFunctionExpr<> dataSimplySupp; // grad w, M
    gsFunctionExpr<> dataFree; // M, grad M

    //gsFunctionExp<>

    real_t thickness = 0.05;
    real_t E = 0;
    real_t nu = 0;
    //real_t nu = 3./10;


    if(testcase_bc == 1)
    {


        f=gsFunctionExpr<>("16*Pi^4*(-16*cos(4*Pi*y)+cos(2*Pi*x)*(-1+25*cos(4*Pi*y)))",2);
        analytSol_w=gsFunctionExpr<>("(1-cos(2*pi*x))*(1-cos(4*pi*y))",2);
        //M = [M_00, M_01, M_10, M_11]
        analytSol_M=gsFunctionExpr<>("-8*Pi^2*Cos(2*Pi*x)*Sin(2*Pi*y)^2","-8*Pi^2*Sin(2*Pi*x)*Sin(4*Pi*y)","-8*Pi^2*Sin(2*Pi*x)*Sin(4*Pi*y)","-32*Pi^2*Cos(4*Pi*y)*Sin(Pi*x)^2",2);
        dataClamped=gsFunctionExpr<>("0","0",2);


        // higher oszillation
        /*
        f=gsFunctionExpr<>("16*Pi^4*(-256*Cos(8*Pi*y) + Cos(6*Pi*x)*(-81 + 625*Cos(8*Pi*y)))",2);
        analytSol_w=gsFunctionExpr<>("(1-cos(6*pi*x))*(1-cos(8*pi*y))",2);
        analytSol_M=gsFunctionExpr<>("-72*Pi^2*Cos(6*Pi*x)*Sin(4*Pi*y)*Sin(4*Pi*y)","-48*Pi^2*Sin(6*Pi*x)*Sin(8*Pi*y)","-48*Pi^2*Sin(6*Pi*x)*Sin(8*Pi*y)","64*Pi^2*(-1 + Cos(6*Pi*x))*Cos(8*Pi*y)",2);
        */

        // with c(u,v) = u v
        /*
        f = gsFunctionExpr<>("1 - (1 + 256*Pi^4)*Cos(4*Pi*y) + Cos(2*Pi*x)*(-1 - 16*Pi^4 + (1 + 400*Pi^4)*Cos(4*Pi*y))",2);
        analytSol_w=gsFunctionExpr<>("(1-cos(2*pi*x))*(1-cos(4*pi*y))",2);
        analytSol_M=gsFunctionExpr<>("-8*Pi^2*Cos(2*Pi*x)*Sin(2*Pi*y)^2","-8*Pi^2*Sin(2*Pi*x)*Sin(4*Pi*y)","-8*Pi^2*Sin(2*Pi*x)*Sin(4*Pi*y)","-32*Pi^2*Cos(4*Pi*y)*Sin(Pi*x)^2",2);
        */

    }
    else if (testcase_bc == 2)
    {
        /*
        f=gsFunctionExpr<>("400*Pi^4*sin(2*Pi*x)*sin(4*Pi*y)",2);
        analytSol_w=gsFunctionExpr<>("sin(2*Pi*x)*sin(4*Pi*y)",2);
        analytSol_M=gsFunctionExpr<>("4*Pi^2*Sin(2*Pi*x)*Sin(4*Pi*y)","-8*Pi^2*Cos(2*Pi*x)*Cos(4*Pi*y)","-8*Pi^2*Cos(2*Pi*x)*Cos(4*Pi*y)","16*Pi^2*Sin(2*Pi*x)*Sin(4*Pi*y)",2);
        dataSimplySupp=gsFunctionExpr<>("0","0","0","0","0","0","0","0","0",2);
        */
    }
    else if (testcase_bc == 3)
    {
        // E W ss, N S free
        f=gsFunctionExpr<>("16*Pi^4*Sin(2*pi*x)",2);
        analytSol_w=gsFunctionExpr<>("sin(2*Pi*x)",2);
        analytSol_M=gsFunctionExpr<>("0","0","0","0",2);
    }
    else if (testcase_bc == 4)
    {
        // E W ss, N S clamped
        /*
        f=gsFunctionExpr<>("25*Pi^3*Cos(Pi*y)*(-Cosh(Pi) + Pi*1./Sinh(Pi))*Sin(Pi*x)*Sin(Pi*y)",2);
        analytSol_w=gsFunctionExpr<>("Sin(Pi*x)*(y*Cosh(Pi*y) + ((-(Cosh(Pi)/Pi) + 1./Sinh(Pi))*Sin(2*Pi*y))/2. - Cosh(Pi)/Sinh(Pi)*Sinh(Pi*y))",2);
        analytSol_M = gsFunctionExpr<>("(Pi*Sin(Pi*x)*(2*Pi*y*Cosh(Pi*y) + (-Cosh(Pi) + Pi*1./Sinh(Pi))*Sin(2*Pi*y) - 2*Pi*Cosh(Pi)/Sinh(Pi)*Sinh(Pi*y)))/2.", "-(Pi*Cos(Pi*x)*(Cosh(Pi*y)*(1 - Pi*Cosh(Pi)/Sinh(Pi)) + Cos(2*Pi*y)*(-Cosh(Pi) + Pi*1./Sinh(Pi)) + Pi*y*Sinh(Pi*y)))", "-(Pi*Cos(Pi*x)*(Cosh(Pi*y)*(1 - Pi*Cosh(Pi)/Sinh(Pi)) + Cos(2*Pi*y)*(-Cosh(Pi) + Pi*1./Sinh(Pi)) + Pi*y*Sinh(Pi*y)))", "-(Pi*Sin(Pi*x)*(Pi*y*Cosh(Pi*y) + 2*(Cosh(Pi) - Pi*1./Sinh(Pi))*Sin(2*Pi*y) + (2 - Pi*Cosh(Pi)/Sinh(Pi))*Sinh(Pi*y)))",2);
        */

    }
    else if (testcase_bc == 5)
    {
        // E W ss, N free, S clamped
        /*
        f=gsFunctionExpr<>("(25*Pi^4*(5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(Pi*x)*Sin(2*Pi*y))/(58*Pi*Cosh(Pi) - 18*Pi*Cosh(3*Pi) + 8*(2 + 5*Pi^2)*Sinh(Pi) - 24*Sinh(3*Pi))",2);
        analytSol_w=gsFunctionExpr<>("(Sin(Pi*x)*((5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(2*Pi*y) + 2*Cosh(Pi*y)*(Pi*(29 + 17*y)*Cosh(Pi) + 3*Pi*(-3 + y)*Cosh(3*Pi) + 4*(2 + Pi^2*(5 + 7*y))*Sinh(Pi) - 12*Sinh(3*Pi)) - 2*(4*(4 + Pi^2*(7 + 5*y))*Cosh(Pi) + 12*Cosh(3*Pi) + Pi*(31 + 19*y)*Sinh(Pi) + 3*Pi*(1 - 3*y)*Sinh(3*Pi))*Sinh(Pi*y)))/(58*Pi*Cosh(Pi) - 18*Pi*Cosh(3*Pi) + 8*(2 + 5*Pi^2)*Sinh(Pi) - 24*Sinh(3*Pi))",2);
        analytSol_M=gsFunctionExpr<>("(Pi^2*Sin(Pi*x)*((5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(2*Pi*y) + 2*Cosh(Pi*y)*(Pi*(29 + 17*y)*Cosh(Pi) + 3*Pi*(-3 + y)*Cosh(3*Pi) + 4*(2 + Pi^2*(5 + 7*y))*Sinh(Pi) - 12*Sinh(3*Pi)) - 2*(4*(4 + Pi^2*(7 + 5*y))*Cosh(Pi) + 12*Cosh(3*Pi) + Pi*(31 + 19*y)*Sinh(Pi) + 3*Pi*(1 - 3*y)*Sinh(3*Pi))*Sinh(Pi*y)))/(58*Pi*Cosh(Pi) - 18*Pi*Cosh(3*Pi) + 8*(2 + 5*Pi^2)*Sinh(Pi) - 24*Sinh(3*Pi))","(Pi^2*Cos(Pi*x)*(-(Cos(2*Pi*y)*(5 + 8*Pi^2 + 3*Cosh(4*Pi))) + Cosh(Pi*y)*((-1 + 4*Pi^2*(7 + 5*y))*Cosh(Pi) + 9*Cosh(3*Pi) + Pi*(3 + 19*y)*Sinh(Pi) + 3*Pi*(1 - 3*y)*Sinh(3*Pi)) - (Pi*(9 + 17*y)*Cosh(Pi) + 3*Pi*(-3 + y)*Cosh(3*Pi) + (-11 + 4*Pi^2*(5 + 7*y))*Sinh(Pi) - 3*Sinh(3*Pi))*Sinh(Pi*y)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi))","(Pi^2*Cos(Pi*x)*(-(Cos(2*Pi*y)*(5 + 8*Pi^2 + 3*Cosh(4*Pi))) + Cosh(Pi*y)*((-1 + 4*Pi^2*(7 + 5*y))*Cosh(Pi) + 9*Cosh(3*Pi) + Pi*(3 + 19*y)*Sinh(Pi) + 3*Pi*(1 - 3*y)*Sinh(3*Pi)) - (Pi*(9 + 17*y)*Cosh(Pi) + 3*Pi*(-3 + y)*Cosh(3*Pi) + (-11 + 4*Pi^2*(5 + 7*y))*Sinh(Pi) - 3*Sinh(3*Pi))*Sinh(Pi*y)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi))","-(Sin(Pi*x)*(Pi^2*Cosh(Pi*y)*(1 + (Pi*y*(17*Cosh(Pi) + 3*Cosh(3*Pi) + 28*Pi*Sinh(Pi)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi))) - (4*Pi^2*(5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(2*Pi*y))/(58*Pi*Cosh(Pi) - 18*Pi*Cosh(3*Pi) + 8*(2 + 5*Pi^2)*Sinh(Pi) - 24*Sinh(3*Pi)) + (2*Pi^2*Cosh(Pi*y)*(-20*Pi*Cosh(Pi) - 19*Sinh(Pi) + 9*Sinh(3*Pi)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi)) + (2*Pi^2*(17*Cosh(Pi) + 3*Cosh(3*Pi) + 28*Pi*Sinh(Pi))*Sinh(Pi*y))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi)) + Pi^2*((Pi*y*(-20*Pi*Cosh(Pi) - 19*Sinh(Pi) + 9*Sinh(3*Pi)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi)) - (4*(4 + 7*Pi^2)*Cosh(Pi) + 12*Cosh(3*Pi) + 31*Pi*Sinh(Pi) + 3*Pi*Sinh(3*Pi))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi)))*Sinh(Pi*y)))",2);
        */

        // N S ss, E free, W clamped
        /*
        f=gsFunctionExpr<>("(25*Pi^4*(5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(2*Pi*x)*Sin(Pi*y))/(58*Pi*Cosh(Pi) - 18*Pi*Cosh(3*Pi) + 8*(2 + 5*Pi^2)*Sinh(Pi) - 24*Sinh(3*Pi))",2);
        analytSol_w=gsFunctionExpr<>("(Sin(Pi*y)*((5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(2*Pi*x) + 2*Cosh(Pi*x)*(Pi*(29 + 17*x)*Cosh(Pi) + 3*Pi*(-3 + x)*Cosh(3*Pi) + 4*(2 + Pi^2*(5 + 7*x))*Sinh(Pi) - 12*Sinh(3*Pi)) - 2*(4*(4 + Pi^2*(7 + 5*x))*Cosh(Pi) + 12*Cosh(3*Pi) + Pi*(31 + 19*x)*Sinh(Pi) + 3*Pi*(1 - 3*x)*Sinh(3*Pi))*Sinh(Pi*x)))/(58*Pi*Cosh(Pi) - 18*Pi*Cosh(3*Pi) + 8*(2 + 5*Pi^2)*Sinh(Pi) - 24*Sinh(3*Pi))",2);
        analytSol_M=gsFunctionExpr<>("-(Sin(Pi*y)*(Pi^2*Cosh(Pi*x)*(1 + (Pi*x*(17*Cosh(Pi) + 3*Cosh(3*Pi) + 28*Pi*Sinh(Pi)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi))) - (4*Pi^2*(5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(2*Pi*x))/(58*Pi*Cosh(Pi) - 18*Pi*Cosh(3*Pi) + 8*(2 + 5*Pi^2)*Sinh(Pi) - 24*Sinh(3*Pi)) + (2*Pi^2*Cosh(Pi*x)*(-20*Pi*Cosh(Pi) - 19*Sinh(Pi) + 9*Sinh(3*Pi)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi)) + (2*Pi^2*(17*Cosh(Pi) + 3*Cosh(3*Pi) + 28*Pi*Sinh(Pi))*Sinh(Pi*x))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi)) + Pi^2*((Pi*x*(-20*Pi*Cosh(Pi) - 19*Sinh(Pi) + 9*Sinh(3*Pi)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi)) - (4*(4 + 7*Pi^2)*Cosh(Pi) + 12*Cosh(3*Pi) + 31*Pi*Sinh(Pi) + 3*Pi*Sinh(3*Pi))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi)))*Sinh(Pi*x)))","(Pi^2*Cos(Pi*y)*(-(Cos(2*Pi*x)*(5 + 8*Pi^2 + 3*Cosh(4*Pi))) + Cosh(Pi*x)*((-1 + 4*Pi^2*(7 + 5*x))*Cosh(Pi) + 9*Cosh(3*Pi) + Pi*(3 + 19*x)*Sinh(Pi) + 3*Pi*(1 - 3*x)*Sinh(3*Pi)) - (Pi*(9 + 17*x)*Cosh(Pi) + 3*Pi*(-3 + x)*Cosh(3*Pi) + (-11 + 4*Pi^2*(5 + 7*x))*Sinh(Pi) - 3*Sinh(3*Pi))*Sinh(Pi*x)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi))","(Pi^2*Cos(Pi*y)*(-(Cos(2*Pi*x)*(5 + 8*Pi^2 + 3*Cosh(4*Pi))) + Cosh(Pi*x)*((-1 + 4*Pi^2*(7 + 5*x))*Cosh(Pi) + 9*Cosh(3*Pi) + Pi*(3 + 19*x)*Sinh(Pi) + 3*Pi*(1 - 3*x)*Sinh(3*Pi)) - (Pi*(9 + 17*x)*Cosh(Pi) + 3*Pi*(-3 + x)*Cosh(3*Pi) + (-11 + 4*Pi^2*(5 + 7*x))*Sinh(Pi) - 3*Sinh(3*Pi))*Sinh(Pi*x)))/(29*Pi*Cosh(Pi) - 9*Pi*Cosh(3*Pi) + 4*(2 + 5*Pi^2)*Sinh(Pi) - 12*Sinh(3*Pi))","(Pi^2*Sin(Pi*y)*((5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(2*Pi*x) + 2*Cosh(Pi*x)*(Pi*(29 + 17*x)*Cosh(Pi) + 3*Pi*(-3 + x)*Cosh(3*Pi) + 4*(2 + Pi^2*(5 + 7*x))*Sinh(Pi) - 12*Sinh(3*Pi)) - 2*(4*(4 + Pi^2*(7 + 5*x))*Cosh(Pi) + 12*Cosh(3*Pi) + Pi*(31 + 19*x)*Sinh(Pi) + 3*Pi*(1 - 3*x)*Sinh(3*Pi))*Sinh(Pi*x)))/(58*Pi*Cosh(Pi) - 18*Pi*Cosh(3*Pi) + 8*(2 + 5*Pi^2)*Sinh(Pi) - 24*Sinh(3*Pi))",2);
        */

        // N S ss, E free, W clamped (less oszillation)

        f=gsFunctionExpr<>("4*Pi^4*Sin(Pi*x)*Sin(Pi*y)",2);
        analytSol_w=gsFunctionExpr<>("(Sin(Pi*y)*((5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(Pi*x) - 2*Cosh(Pi*x)*(Sinh(Pi) + Pi*((7 + 4*x)*Cosh(Pi) - 3*Cosh(3*Pi) + 4*Pi*(1 + 2*x)*Sinh(Pi)) - 3*Sinh(3*Pi)) + 2*((5 + 4*Pi^2*(2 + x))*Cosh(Pi) + 3*Cosh(3*Pi) + Pi*(8 + 5*x)*Sinh(Pi) - 3*Pi*x*Sinh(3*Pi))*Sinh(Pi*x)))/(5 + 8*Pi^2 + 3*Cosh(4*Pi))",2);
        analytSol_M=gsFunctionExpr<>("(Pi^2*Sin(Pi*y)*((5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(Pi*x) + 4*Cosh(Pi*x)*(Pi*Cosh(Pi)*(1 + 2*x - 3*Cosh(2*Pi)) + 2*Pi^2*(1 + 2*x)*Sinh(Pi) + 6*Sinh(Pi)^3) - 4*(2*Pi^2*(2 + x)*Cosh(Pi) + Pi*(-4 + x - 3*x*Cosh(2*Pi))*Sinh(Pi) + 6*Cosh(Pi)*Sinh(Pi)^2)*Sinh(Pi*x)))/(5 + 8*Pi^2 + 3*Cosh(4*Pi))", "(Pi^2*Cos(Pi*y)*(-(Cos(Pi*x)*(5 + 8*Pi^2 + 3*Cosh(4*Pi))) - 2*Cosh(Pi*x)*(Cosh(Pi) + 3*Cosh(3*Pi) + Pi*(4*Pi*(2 + x)*Cosh(Pi) + 5*x*Sinh(Pi) - 3*x*Sinh(3*Pi))) + 8*(-Sinh(Pi) + Pi^2*(1 + 2*x)*Sinh(Pi) + Pi*Cosh(Pi)*(x - 3*Sinh(Pi)^2))*Sinh(Pi*x)))/(5 + 8*Pi^2 + 3*Cosh(4*Pi))", "(Pi^2*Cos(Pi*y)*(-(Cos(Pi*x)*(5 + 8*Pi^2 + 3*Cosh(4*Pi))) - 2*Cosh(Pi*x)*(Cosh(Pi) + 3*Cosh(3*Pi) + Pi*(4*Pi*(2 + x)*Cosh(Pi) + 5*x*Sinh(Pi) - 3*x*Sinh(3*Pi))) + 8*(-Sinh(Pi) + Pi^2*(1 + 2*x)*Sinh(Pi) + Pi*Cosh(Pi)*(x - 3*Sinh(Pi)^2))*Sinh(Pi*x)))/(5 + 8*Pi^2 + 3*Cosh(4*Pi))", "(Pi^2*Sin(Pi*y)*((5 + 8*Pi^2 + 3*Cosh(4*Pi))*Sin(Pi*x) - 2*Cosh(Pi*x)*(Sinh(Pi) + Pi*((7 + 4*x)*Cosh(Pi) - 3*Cosh(3*Pi) + 4*Pi*(1 + 2*x)*Sinh(Pi)) - 3*Sinh(3*Pi)) + 2*((5 + 4*Pi^2*(2 + x))*Cosh(Pi) + 3*Cosh(3*Pi) + Pi*(8 + 5*x)*Sinh(Pi) - 3*Pi*x*Sinh(3*Pi))*Sinh(Pi*x)))/(5 + 8*Pi^2 + 3*Cosh(4*Pi))", 2);
        dataClamped=gsFunctionExpr<>("0","0",2);
        dataSimplySupp=gsFunctionExpr<>("0","0","0","0","0","0","0","0","0",2);
        dataFree=gsFunctionExpr<>("0","0","0","0","0","0","0","0","0",2);
    }

    else if (testcase_bc == 6)
    {
        //f=gsFunctionExpr<>("sin(Pi*x)*cos(2*Pi*y)",2);
        f=gsFunctionExpr<>("1",2);
        analytSol_w=gsFunctionExpr<>("0",2);
        analytSol_M=gsFunctionExpr<>("0","0","0","0",2);
        dataClamped=gsFunctionExpr<>("0","0",2);
        dataSimplySupp=gsFunctionExpr<>("0","0","0","0","0","0","0","0","0",2);
        dataFree=gsFunctionExpr<>("0","0","0","0","0","0","0","0","0",2);

    }

    else if (testcase_bc == 7)
    {
        //f=gsFunctionExpr<>("sin(Pi*x)*cos(2*Pi*y)",2);
        f=gsFunctionExpr<>("1",2);

        // nu=0
        //analytSol_w=gsFunctionExpr<>("5/64-3/32*(x^2+y^2)+1/64*(x^2+y^2)^2",2);
        //analytSol_M=gsFunctionExpr<>("(3 - 3*x^2 - y^2)/16","-(x*y)/8","-(x*y)/8","(3 - x^2 - 3*y^2)/16",2);

        // nu=1/3
        //analytSol_w=gsFunctionExpr<>("1/16-5/64*(x^2+y^2)+1/64*(x^2+y^2)^2",2);
        //analytSol_M=gsFunctionExpr<>("(5 - 5*x^2 - 3*y^2)/24","-(x*y)/12","-(x*y)/12","(5 - 3*x^2 - 5*y^2)/24",2);

        // nu = 0.3
        analytSol_w=gsFunctionExpr<>("53/832-33/416*(x^2+y^2)+1/64*(x^2+y^2)^2",2);
        analytSol_M=gsFunctionExpr<>("(33 - 33*x^2 - 19*y^2)/160","-(7*x*y)/80","-(7*x*y)/80","(33 - 19*x^2 - 33*y^2)/160",2);

    }

    if(nonHomogeneous_bc)
    {

/*
        analytSol_w=gsFunctionExpr<>("Cos(Pi*x)*Cos(2*Pi*y)",2);
        f=gsFunctionExpr<>("25*Pi^4*Cos(Pi*x)*Cos(2*Pi*y)",2);
        //M = [M_00, M_01, M_10, M_11]
        analytSol_M=gsFunctionExpr<>("Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","-2*Pi^2*Sin(Pi*x)*Sin(2*Pi*y)","-2*Pi^2*Sin(Pi*x)*Sin(2*Pi*y)","4*Pi^2*Cos(Pi*x)*Cos(2*Pi*y)",2);
        dataClamped=gsFunctionExpr<>("-(Pi*Cos(2*Pi*y)*Sin(Pi*x))","-2*Pi*Cos(Pi*x)*Sin(2*Pi*y)",2);
        dataSimplySupp=gsFunctionExpr<>("-(Pi*Cos(2*Pi*y)*Sin(Pi*x))","-2*Pi*Cos(Pi*x)*Sin(2*Pi*y)","Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","4*Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","-2*Pi^2*Sin(Pi*x)*Sin(2*Pi*y)","0","0","0","0",2);
        dataFree=gsFunctionExpr<>("Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","4*Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","-2*Pi^2*Sin(Pi*x)*Sin(2*Pi*y)","-(Pi^3*Cos(2*Pi*y)*Sin(Pi*x))","-4*Pi^3*Cos(2*Pi*y)*Sin(Pi*x)","-2*Pi^3*Cos(Pi*x)*Sin(2*Pi*y)","-2*Pi^3*Cos(Pi*x)*Sin(2*Pi*y)","-8*Pi^3*Cos(Pi*x)*Sin(2*Pi*y)","-4*Pi^3*Cos(2*Pi*y)*Sin(Pi*x)",2);
*/


        analytSol_w=gsFunctionExpr<>("Cos(Pi*x)*Cos(2*Pi*y)+2*x^4*y^4",2);
        f=gsFunctionExpr<>("48*(x^4 + 12*x^2*y^2 + y^4) + 25*Pi^4*Cos(Pi*x)*Cos(2*Pi*y)",2);
        //M = [M_00, M_01, M_10, M_11]
        analytSol_M=gsFunctionExpr<>("-24*x^2*y^4 + Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","-32*x^3*y^3 - 2*Pi^2*Sin(Pi*x)*Sin(2*Pi*y)","-32*x^3*y^3 - 2*Pi^2*Sin(Pi*x)*Sin(2*Pi*y)","-24*x^4*y^2 + 4*Pi^2*Cos(Pi*x)*Cos(2*Pi*y)",2);
        dataClamped=gsFunctionExpr<>("8*x^3*y^4 - Pi*Cos(2*Pi*y)*Sin(Pi*x)","8*x^4*y^3 - 2*Pi*Cos(Pi*x)*Sin(2*Pi*y)",2);
        dataSimplySupp=gsFunctionExpr<>("8*x^3*y^4 - Pi*Cos(2*Pi*y)*Sin(Pi*x)","8*x^4*y^3 - 2*Pi*Cos(Pi*x)*Sin(2*Pi*y)","-24*x^2*y^4 + Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","-24*x^4*y^2 + 4*Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","-32*x^3*y^3 - 2*Pi^2*Sin(Pi*x)*Sin(2*Pi*y)","0","0","0","0",2);
        dataFree=gsFunctionExpr<>("-24*x^2*y^4 + Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","-24*x^4*y^2 + 4*Pi^2*Cos(Pi*x)*Cos(2*Pi*y)","-32*x^3*y^3 - 2*Pi^2*Sin(Pi*x)*Sin(2*Pi*y)","-48*x*y^4 - Pi^3*Cos(2*Pi*y)*Sin(Pi*x)","-96*x^3*y^2 - 4*Pi^3*Cos(2*Pi*y)*Sin(Pi*x)","-96*x^2*y^3 - 2*Pi^3*Cos(Pi*x)*Sin(2*Pi*y)","-96*x^2*y^3 - 2*Pi^3*Cos(Pi*x)*Sin(2*Pi*y)","-48*x^4*y - 8*Pi^3*Cos(Pi*x)*Sin(2*Pi*y)","-96*x^3*y^2 - 4*Pi^3*Cos(2*Pi*y)*Sin(Pi*x)",2);

    }
    else
    {
        dataClamped=gsFunctionExpr<>("0","0",2);
        dataSimplySupp=gsFunctionExpr<>("0","0","0","0","0","0","0","0","0",2);
        dataFree=gsFunctionExpr<>("0","0","0","0","0","0","0","0","0",2);

    }

    
    /// Setup geometry ///

    gsMultiPatch<> mpptr;
    gsMultiPatch<> mp;

    if (testcase_geometry == 1)
    {
        mpptr = gsNurbsCreator<>::BSplineSquareGrid(1, 1, 2, -1, -1);

        //gsReadFile<>("square_shifted.xml", mpptr);
    }
    else if(testcase_geometry == 2)
    {

        // BSplineFatQuarterAnnulus
        /*
        //mpptr = gsNurbsCreator<>::BSplineFatQuarterAnnulus();
        mpptr = gsMultiPatch<>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus());

        // s.t. p=2 in both directions
        mpptr.patch(0).degreeElevate(1,0);
        */

        // L-shape
        //mpptr = gsMultiPatch<>(*gsNurbsCreator<>::BSplineLShape_p1(1));

        // Disc (with radius 1)
        mpptr = gsMultiPatch<>(*gsNurbsCreator<>::NurbsDisk(0.5));

    }
    // unit square (4 patches) p=1
    else if (testcase_geometry == 3)
    {
        mpptr = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 1, -1, -1);
    }

    // L-shape (2 patches) p=1
    else if (testcase_geometry == 4)
    {
        gsReadFile<>("planar/lshape2d_2patches.xml", mpptr);
    }

    mp=mpptr;

    cout<<"geometry "<< mp<<endl;
    cout<<mp[0].basis()<<endl;
    //cout<<mp[0].coefs()<<endl;

    // geometry: elevate degree (smoothness is increased, if no interior knots are present)
    //mp.degreeElevate(numIncreaseDegree);
    //mpptr.patch(0).degreeReduce(1);

    // Shift control points
    /*
    //mp.degreeElevate();
    //mp.uniformRefine();
    gsMatrix<> & coefs = mpptr->patch(0).coefs();
    cout<<coefs<<endl<<endl;
    coefs(0,1)=-3;
    //coefs(0,0)=-3;
    cout<<coefs<<endl;


    // print
    cout<<endl;
    mpptr->piece(0).print(cout);
    cout<<endl;
*/

    
    /// Setup boundary conditions ///

    enum bcType_plate
    {
        clamped = 0,
        simplysupp   = 1,
        free     = 2
    };
    
    gsBoundaryConditions<> BCs;
    int BCs_patch[5];

    std::vector<std::vector<gsShellMixedPde<>::corner_type> > corners(mp.size());
    std::vector<std::vector<gsMatrix<> > > cornerCouplingCoefs(mp.size());
    std::vector<std::vector<int >  > indFreeComp(mp.size());
    corners[0].resize(5);
    cornerCouplingCoefs[0].resize(5);

    corners[0][1]= gsShellMixedPde<>::other;
    corners[0][2]= gsShellMixedPde<>::other;
    corners[0][3]= gsShellMixedPde<>::other;
    corners[0][4]= gsShellMixedPde<>::other;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();


    if(testcase_bc == 1)
    {
        BCs_patch[boundary::north]=clamped;
        BCs_patch[boundary::east]=clamped;
        BCs_patch[boundary::south]=clamped;
        BCs_patch[boundary::west]=clamped;

    }
    else if(testcase_bc == 2)
    {
        BCs_patch[boundary::north]=simplysupp;
        BCs_patch[boundary::east]=simplysupp;
        BCs_patch[boundary::south]=simplysupp;
        BCs_patch[boundary::west]=simplysupp;
    }

    else if(testcase_bc == 3)
    {
        // fsfs (N E S W)

        BCs_patch[boundary::north]=free;
        BCs_patch[boundary::east]=simplysupp;
        BCs_patch[boundary::south]=free;
        BCs_patch[boundary::west]=simplysupp;


        corners[0][1]= gsShellMixedPde<>::sf;
        corners[0][2]= gsShellMixedPde<>::sf;
        corners[0][3]= gsShellMixedPde<>::sf;
        corners[0][4]= gsShellMixedPde<>::sf;

        /*
        cornerCouplingCoefs[0][1].resize(1,4);
        cornerCouplingCoefs[0][1]<< 0, 1, -1, 0;
        cornerCouplingCoefs[0][2].resize(1,4);
        cornerCouplingCoefs[0][2]<< 1, 0, 0, 1;
        cornerCouplingCoefs[0][3].resize(1,4);
        cornerCouplingCoefs[0][3]<< 1, 0, 0, 1;
        cornerCouplingCoefs[0][4].resize(1,4);
        cornerCouplingCoefs[0][4]<< 0, 1, -1, 0;
        */

        // coef(0,1) = -3 (free not prallel)

        cornerCouplingCoefs[0][1].resize(1,4);
        cornerCouplingCoefs[0][1]<< 0, 1, -2./sqrt(2), 1;
        cornerCouplingCoefs[0][2].resize(1,4);
        cornerCouplingCoefs[0][2]<< 1, 0, 1, 2./sqrt(2);
        cornerCouplingCoefs[0][3].resize(1,4);
        cornerCouplingCoefs[0][3]<< 1, 0, 0, 1;
        cornerCouplingCoefs[0][4].resize(1,4);
        cornerCouplingCoefs[0][4]<< 0, 1, -1, 0;



        // coef(0,0) = -3 (ss not prallel)
        /*
        cornerCouplingCoefs[0][1].resize(1,4);
        cornerCouplingCoefs[0][1]<< 0, 1, -2./sqrt(2), 1;
        cornerCouplingCoefs[0][2].resize(1,4);
        cornerCouplingCoefs[0][2]<< 1, 0, 0, 1;
        cornerCouplingCoefs[0][3].resize(1,4);
        cornerCouplingCoefs[0][3]<< 1, 0, 1, 2./sqrt(2);
        cornerCouplingCoefs[0][4].resize(1,4);
        cornerCouplingCoefs[0][4]<< 0, 1, -1, 0;
*/


    }
    else if(testcase_bc == 4)
    {

        BCs_patch[boundary::north]=clamped;
        BCs_patch[boundary::east]=simplysupp;
        BCs_patch[boundary::south]=clamped;
        BCs_patch[boundary::west]=simplysupp;
    }
    else if(testcase_bc == 5)
    {


        // sfsc (N E S W)
        BCs_patch[boundary::north]=simplysupp;
        BCs_patch[boundary::east]=free;
        BCs_patch[boundary::south]=simplysupp;
        BCs_patch[boundary::west]=clamped;

        // set corner coupling for lagrange multipliers
        corners[0][1]= gsShellMixedPde<>::other;
        corners[0][2]= gsShellMixedPde<>::sf;
        cornerCouplingCoefs[0][2].resize(1,4);
        cornerCouplingCoefs[0][2]<< 0, 1, -1, 0;
        corners[0][3]= gsShellMixedPde<>::other;
        corners[0][4]= gsShellMixedPde<>::sf;
        cornerCouplingCoefs[0][4].resize(1,4);
        cornerCouplingCoefs[0][4]<< 1, 0, 0, 1;



        // fscs (N E S W)
        /*
        BCs_patch[boundary::north]=free;
        BCs_patch[boundary::east]=simplysupp;
        BCs_patch[boundary::south]=clamped;
        BCs_patch[boundary::west]=simplysupp;

        corners[0][1]= gsShellMixedPde<>::other;
        corners[0][2]= gsShellMixedPde<>::other;
        corners[0][3]= gsShellMixedPde<>::sf;
        cornerCouplingCoefs[0][3].resize(1,4);
        cornerCouplingCoefs[0][3]<< 1, 0, 0, 1;
        corners[0][4]= gsShellMixedPde<>::sf;
        cornerCouplingCoefs[0][4].resize(1,4);
        cornerCouplingCoefs[0][4]<< 0, 1, -1, 0;
        */

    }
    else if(testcase_bc == 6)
    {

        // cfcc (N E S W)
/*
        BCs_patch[boundary::north]=clamped;
        BCs_patch[boundary::east]= free;
        BCs_patch[boundary::south]= clamped;
        BCs_patch[boundary::west]=clamped;
*/

        // cscc (N E S W)
        /*
        BCs_patch[boundary::north]=clamped;
        BCs_patch[boundary::east]= simplysupp;
        BCs_patch[boundary::south]=clamped;
        BCs_patch[boundary::west]=clamped;
*/
        // cfcs (N E S W)

        BCs_patch[boundary::north]= clamped;
        BCs_patch[boundary::east]= free;
        BCs_patch[boundary::south]= clamped;
        BCs_patch[boundary::west]= clamped;



        // cffc (N E S W)
        // southeast corner lambdaT eliminieren
        /*
        BCs_patch[boundary::north]=clamped;
        BCs_patch[boundary::east]= free;
        BCs_patch[boundary::south]=free;
        BCs_patch[boundary::west]=clamped;

        corners[0][1]= gsShellMixedPde<>::other;
        corners[0][2]= gsShellMixedPde<>::ff;
        cornerCouplingCoefs[0][2].resize(2,4);
        cornerCouplingCoefs[0][2]<< 1, 0, 0, 1, 0, 1, -1, 0;
        corners[0][3]= gsShellMixedPde<>::other;
        corners[0][4]= gsShellMixedPde<>::other;
        */



/*
        // fcfs (N E S W)
        BCs_patch[boundary::north]=free;
        BCs_patch[boundary::east]=clamped;
        BCs_patch[boundary::south]=free;
        BCs_patch[boundary::west]=simplysupp;

        // set corner coupling for lagrange multipliers
        corners[0][1]= gsShellMixedPde<>::sf;
        corners[0][2]= gsShellMixedPde<>::none;
        corners[0][3]= gsShellMixedPde<>::sf;
        corners[0][4]= gsShellMixedPde<>::none;

        cornerCouplingCoefs[0][1].resize(1,4);
        cornerCouplingCoefs[0][1]<< 0, 1, -1, 0;
        cornerCouplingCoefs[0][2].resize(1,4);
        cornerCouplingCoefs[0][2]<< 1, 0, 0, 1;
        cornerCouplingCoefs[0][3].resize(1,4);
        cornerCouplingCoefs[0][3]<< 1, 0, 0, 1;
        cornerCouplingCoefs[0][4].resize(1,4);
        cornerCouplingCoefs[0][4]<< 0, 1, -1, 0;
*/

    }
    else if(testcase_bc == 7)
    {
        BCs_patch[boundary::north]=simplysupp;
        BCs_patch[boundary::east]= simplysupp;
        BCs_patch[boundary::south]=simplysupp;
        BCs_patch[boundary::west]=simplysupp;


        corners[0][1]= gsShellMixedPde<>::none;
        corners[0][2]= gsShellMixedPde<>::none;
        corners[0][3]= gsShellMixedPde<>::none;
        corners[0][4]= gsShellMixedPde<>::none;

    }

    if(testcase_geometry == 1 || testcase_geometry == 2)
    {
        /*
        // p
        if(!options.getSwitch("pBoundary0"))
        {
            if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
                BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0 );
            if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
                BCs.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0 );
            if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
                BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0 );
            if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
                BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 );
        }
        else
        {
            BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0 );
            BCs.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0 );
            BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0 );
            BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 );
        }
        */

        // w
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(0,boundary::north, condition_type::dirichlet, analytSol_w, 3 );
        if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(0,boundary::east, condition_type::dirichlet, analytSol_w, 3 );
        if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::dirichlet, analytSol_w, 3 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::dirichlet, analytSol_w, 3 );


        if(BCs_patch[boundary::north]==clamped)
            BCs.add(0,boundary::north, "Clamped", dataClamped, 3 );
        if(BCs_patch[boundary::east]==clamped)
            BCs.add(0,boundary::east, "Clamped", dataClamped, 3 );
        if(BCs_patch[boundary::south]==clamped)
            BCs.add(0,boundary::south, "Clamped", dataClamped, 3 );
        if(BCs_patch[boundary::west]==clamped)
            BCs.add(0,boundary::west, "Clamped", dataClamped, 3 );



        // phi
        if(BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(0,boundary::north, condition_type::neumann, dataSimplySupp, 1);
        if(BCs_patch[boundary::north]==free)
            BCs.addCondition(0,boundary::north, condition_type::robin, dataFree, 1 );

        if(BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(0,boundary::east, condition_type::neumann, dataSimplySupp, 1);
        if(BCs_patch[boundary::east]==free)
            BCs.addCondition(0,boundary::east, condition_type::robin, dataFree, 1 );

        if(BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::neumann, dataSimplySupp, 1);
        if(BCs_patch[boundary::south]==free)
            BCs.addCondition(0,boundary::south, condition_type::robin, dataFree, 1 );

        if(BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::neumann, dataSimplySupp, 1);
        if(BCs_patch[boundary::west]==free)
            BCs.addCondition(0,boundary::west, condition_type::robin, dataFree, 1 );


        // lagrange multiplier lambdaN
        if(BCs_patch[boundary::north]==clamped)
            BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 4);
        if(BCs_patch[boundary::east]==clamped)
            BCs.addCondition(0,boundary::east, condition_type::dirichlet, 0, 4);
        if(BCs_patch[boundary::south]==clamped)
            BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 4);
        if(BCs_patch[boundary::west]==clamped)
            BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 4);

        /*
        BCs.addCornerValue(boundary::northeast, 0, 0, 4);
        BCs.addCornerValue(boundary::northwest, 0, 0, 4);
        BCs.addCornerValue(boundary::southeast, 0, 0, 4);
        BCs.addCornerValue(boundary::southwest, 0, 0, 4);
        */

    }

    else if(testcase_geometry == 3)
    {
        // patch 0
        // p
        if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 );

        // w
        if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 3 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 3 );

        // phi
        if(BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::south]==free)
            BCs.addCondition(0,boundary::south, condition_type::robin, 0, 1 );

        if(BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::west]==free)
            BCs.addCondition(0,boundary::west, condition_type::robin, 0, 1 );

        // patch 1
        // p
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0 );

        // w
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(1,boundary::north, condition_type::dirichlet, 0, 3 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(1,boundary::west, condition_type::dirichlet, 0, 3 );

        // phi
        if(BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(1,boundary::north, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::north]==free)
            BCs.addCondition(1,boundary::north, condition_type::robin, 0, 1 );

        if(BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(1,boundary::west, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::west]==free)
            BCs.addCondition(1,boundary::west, condition_type::robin, 0, 1 );

        // patch 2
        // p
        if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0 );

        // w
        if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(2,boundary::east, condition_type::dirichlet, 0, 3 );
        if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(2,boundary::south, condition_type::dirichlet, 0, 3 );

        // phi
        if(BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(2,boundary::east, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::east]==free)
            BCs.addCondition(2,boundary::east, condition_type::robin, 0, 1 );

        if(BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(2,boundary::south, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::south]==free)
            BCs.addCondition(2,boundary::south, condition_type::robin, 0, 1 );

        // patch 3
        // p
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0 );

        // w
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(3,boundary::north, condition_type::dirichlet, 0, 3 );
        if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(3,boundary::east, condition_type::dirichlet, 0, 3 );

        // phi
        if(BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(3,boundary::north, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::north]==free)
            BCs.addCondition(3,boundary::north, condition_type::robin, 0, 1 );

        if(BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(3,boundary::east, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::east]==free)
            BCs.addCondition(3,boundary::east, condition_type::robin, 0, 1 );
    }

    if(testcase_geometry == 4)
    {
        // patch 0
        // p
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 );

        // w
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 3 );
        if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 3 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 3 );

        // phi
        if(BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(0,boundary::north, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::north]==free)
            BCs.addCondition(0,boundary::north, condition_type::robin, 0, 1 );

        if(BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::south]==free)
            BCs.addCondition(0,boundary::south, condition_type::robin, 0, 1 );

        if(BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::west]==free)
            BCs.addCondition(0,boundary::west, condition_type::robin, 0, 1 );

        // patch 1
        // p
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0 );

        // w
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(1,boundary::north, condition_type::dirichlet, 0, 3 );
        if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(1,boundary::east, condition_type::dirichlet, 0, 3 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(1,boundary::west, condition_type::dirichlet, 0, 3 );

        // phi
        if(BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(1,boundary::north, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::north]==free)
            BCs.addCondition(1,boundary::north, condition_type::robin, 0, 1 );

        if(BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(1,boundary::east, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::east]==free)
            BCs.addCondition(1,boundary::east, condition_type::robin, 0, 1 );

        if(BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(1,boundary::west, condition_type::neumann, 0, 1);
        if(BCs_patch[boundary::west]==free)
            BCs.addCondition(1,boundary::west, condition_type::robin, 0, 1 );



    }


    gsBoundaryConditions<> BCsFine(BCs);

    // p BCs
    if(!options.getSwitch("pBoundary0"))
    {
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 );
    }
    else
    {
        BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0 );
        BCs.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0 );
        BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0 );
        BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 );
    }

    // p BCsFine
    if(!optionsFine.getSwitch("pBoundary0"))
    {
        if(BCs_patch[boundary::north]==clamped|| BCs_patch[boundary::north]==simplysupp)
            BCsFine.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::east]==clamped|| BCs_patch[boundary::east]==simplysupp)
            BCsFine.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::south]==clamped|| BCs_patch[boundary::south]==simplysupp)
            BCsFine.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0 );
        if(BCs_patch[boundary::west]==clamped|| BCs_patch[boundary::west]==simplysupp)
            BCsFine.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 );
    }
    else
    {
        BCsFine.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0 );
        BCsFine.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0 );
        BCsFine.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0 );
        BCsFine.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 );
    }

    // lagrange multiplier lambdaT BCs
    if(!options.getSwitch("pBoundary0"))
    {

        if(BCs_patch[boundary::north]==clamped || BCs_patch[boundary::north]==simplysupp)
            BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 5);
        if(BCs_patch[boundary::east]==clamped || BCs_patch[boundary::east]==simplysupp)
            BCs.addCondition(0,boundary::east, condition_type::dirichlet, 0, 5);
        if(BCs_patch[boundary::south]==clamped || BCs_patch[boundary::south]==simplysupp)
            BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 5);
        if(BCs_patch[boundary::west]==clamped || BCs_patch[boundary::west]==simplysupp)
            BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 5);


        BCs.addCornerValue(boundary::northeast, 0, 0, 5);
        BCs.addCornerValue(boundary::northwest, 0, 0, 5);
        BCs.addCornerValue(boundary::southeast, 0, 0, 5);
        BCs.addCornerValue(boundary::southwest, 0, 0, 5);

    }
    else
    {
        BCs.addCondition(0,boundary::north, condition_type::dirichlet, 0, 5);
        BCs.addCondition(0,boundary::east, condition_type::dirichlet, 0, 5);
        BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 5);
        BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 5);
    }



    // lagrange multiplier lambdaT BCsFine
    if(!optionsFine.getSwitch("pBoundary0"))
    {
        if(BCs_patch[boundary::north]==clamped || BCs_patch[boundary::north]==simplysupp)
            BCsFine.addCondition(0,boundary::north, condition_type::dirichlet, 0, 5);
        if(BCs_patch[boundary::east]==clamped || BCs_patch[boundary::east]==simplysupp)
            BCsFine.addCondition(0,boundary::east, condition_type::dirichlet, 0, 5);
        if(BCs_patch[boundary::south]==clamped || BCs_patch[boundary::south]==simplysupp)
            BCsFine.addCondition(0,boundary::south, condition_type::dirichlet, 0, 5);
        if(BCs_patch[boundary::west]==clamped || BCs_patch[boundary::west]==simplysupp)
            BCsFine.addCondition(0,boundary::west, condition_type::dirichlet, 0, 5);


        BCsFine.addCornerValue(boundary::northeast, 0, 0, 5);
        BCsFine.addCornerValue(boundary::northwest, 0, 0, 5);
        BCsFine.addCornerValue(boundary::southeast, 0, 0, 5);
        BCsFine.addCornerValue(boundary::southwest, 0, 0, 5);

    }
    else
    {
        BCsFine.addCondition(0,boundary::north, condition_type::dirichlet, 0, 5);
        BCsFine.addCondition(0,boundary::east, condition_type::dirichlet, 0, 5);
        BCsFine.addCondition(0,boundary::south, condition_type::dirichlet, 0, 5);
        BCsFine.addCondition(0,boundary::west, condition_type::dirichlet, 0, 5);
    }

    // phi corner values BCs
    //BCs.addCornerValue(boundary::northeast, 0, 0, 1);
    //BCs.addCornerValue(boundary::northeast, 0, 0, 2);

    //BCs.addCornerValue(boundary::southwest, 0, 0, 1); // cscs: southeast gives nan
    //BCs.addCornerValue(boundary::southwest, 0, 0, 2);

    // phi corner values BCsFine
    //BCsFine.addCornerValue(boundary::northeast, 0, 0, 1);
    //BCsFine.addCornerValue(boundary::northeast, 0, 0, 2);

    //BCsFine.addCornerValue(boundary::southwest, 0, 0, 1);
    //BCsFine.addCornerValue(boundary::southwest, 0, 0, 2);


    /// Order boundary ///

    std::vector<patchSide>& boundaries = mpptr.boundaries();
    std::vector<patchSide> boundaries_copy(boundaries);


    if(testcase_geometry == 1)
    {
        // gsNurbsCreator<>::BSplineSquareGrid
        // starting with south (counter-clockwise)


        boundaries[0] = boundaries_copy[1];
        boundaries[1] = boundaries_copy[2];
        boundaries[2] = boundaries_copy[0];
        boundaries[3] = boundaries_copy[3];



        // gsNurbsCreator<>::BSplineSquareGrid
        // starting with east (counter-clockwise)
/*
        boundaries[0] = boundaries_copy[2];
        boundaries[1] = boundaries_copy[0];
        boundaries[2] = boundaries_copy[3];
        boundaries[3] = boundaries_copy[1];
*/


        // square_shifted.xml
        // starting with south (counter-clockwise)
        /*
        boundaries[0] = boundaries_copy[2];
        boundaries[1] = boundaries_copy[1];
        boundaries[2] = boundaries_copy[3];
        boundaries[3] = boundaries_copy[0];
*/

    }

    else if(testcase_geometry == 2)
    {
        // quarter annulus
        // starting with south (counter-clockwise)
        boundaries[0] = boundaries_copy[2];
        boundaries[1] = boundaries_copy[1];
        boundaries[2] = boundaries_copy[3];
        boundaries[3] = boundaries_copy[0];

    }



    /// PlateMixedPde ///

    gsShellMixedPde<> pde(mpptr, BCs, f, pLoads, E, nu, thickness, corners, cornerCouplingCoefs, indFreeComp);
    gsShellMixedPde<> pdeFine(mpptr, BCsFine, f, pLoads, E, nu, thickness, corners, cornerCouplingCoefs, indFreeComp);



    /// Refinement h and p ///

    // Copy basis from the geometry
    gsMultiBasis<> basis(mpptr);

    // For NURBS geometry use associated B-spline basis

    gsBasis<>& nBasis = mpptr.basis(0);
    gsTensorNurbsBasis<2, real_t>* nBasis_cast = dynamic_cast< gsTensorNurbsBasis<2, real_t>* >(&nBasis);
    if(nBasis_cast != NULL)
    {
        cout<< "NURBS geometry given, corresponding BSpline basis used for compuatitions"<<endl;
        basis = nBasis_cast->source();
    }


    // elevate degree
    //basis.degreeElevate(numElevateDegree);

    if(testcase_geometry == 2) // Disc
        basis.degreeReduce(1);


    // increase degree
    basis.degreeIncrease(numIncreaseDegree);

    gsMultiBasis<> basis_degMinus1(basis);
    basis_degMinus1.degreeReduce(1);


    // h-refine each basis
    for (int i = 0; i < numHRefine_start; ++i)
    {
        if(keepC0)
            basis.uniformRefine(1, basis.maxDegree(0)); // multipatch: basis of patch 0 is considered
        else
        {
            basis.uniformRefine();
            basis_degMinus1.uniformRefine();
        }

    }

    /// Compute fine solution ///

    if(!analytSol || convRates_p_phi)
    {
        cout<<"numHRefine "<<numHRefine_end+fineHRefine_add<<endl;

        gsMultiBasis<> basisFine( basis);
        gsMultiBasis<> basisFine_degMinus1( basis_degMinus1);

        for (int i = numHRefine_start; i < numHRefine_end+fineHRefine_add; ++i)
        {
            if(keepC0)
                basisFine.uniformRefine(1, basisFine.maxDegree(0)); // multipatch: basis of patch 0 is considered
            else
            {
                basisFine.uniformRefine();
                basisFine_degMinus1.uniformRefine();
            }
        }

        /// Assemble and Solve ///

        // construct bases

        std::vector< gsMultiBasis<> >  bases;

        bases.push_back(basisFine); // basis for p, w
        bases.push_back(basisFine_degMinus1); // basis for lambda

        // basis for phi
        /*
        gsMultiBasis<> basisFine_copy( basisFine);
        basisFine_copy.degreeIncrease(1);
        bases.push_back( basisFine_copy);
        */


        gsPlateMixedAssembler<real_t> assembler(pdeFine, bases, optionsFine);

        // Generate system matrix and load vector
        watch.restart();
        assembler.assemble();
        cout<<"Assembling time: "<<watch.stop()<<endl;

        gsSparseSystem<> spSystem = assembler.system();

        if(! options.getSwitch("solveWholeSystem"))
        {

            gsVector<index_t> blocks(7);
            blocks<< 0,1,1,2,3,4,5;
            gsSparseMatrix<>::BlockView matrixBlockView = spSystem.blockView(6, 6, blocks, blocks);
            gsMatrix<>::BlockView rhsBlockView = spSystem.blockViewRhs(6, blocks);

            gsSparseMatrix<> A00 = matrixBlockView(0,0);
            gsSparseMatrix<> A10 = matrixBlockView(1,0);
            gsSparseMatrix<> A01 = matrixBlockView(0,1);
            gsSparseMatrix<> A11 = matrixBlockView(1,1);
            gsSparseMatrix<> B = matrixBlockView(2,0);
            gsMatrix<> fq = rhsBlockView(0,0);
            gsMatrix<> fpsi = rhsBlockView(1,0);
            gsMatrix<> fv = rhsBlockView(2,0);

            gsSparseMatrix<> Ln0 = matrixBlockView(3,0);
            gsSparseMatrix<> Ln1 = matrixBlockView(3,1);
            gsSparseMatrix<> Lt0 = matrixBlockView(4,0);
            gsSparseMatrix<> Lt1 = matrixBlockView(4,1);
            gsMatrix<> fmuN = rhsBlockView(3,0);
            gsMatrix<> fmuT = rhsBlockView(4,0);

            gsSparseMatrix<> LnMean = matrixBlockView(5,3);
            gsSparseMatrix<> LtMean = matrixBlockView(5,4);


            // add Ln0, Lt0 to A10 and Ln1, Lt1 to A11
            index_t A11_rows= A11.rows();
            index_t A11_cols= A11.cols();
            index_t A10_rows= A10.rows();
            index_t A10_cols= A10.cols();
            A11.conservativeResize(A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+spSystem.colMapper(6).freeSize(), A11_cols+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+spSystem.colMapper(6).freeSize());
            A10.conservativeResize(A10_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+spSystem.colMapper(6).freeSize(), A10_cols);
            fpsi.conservativeResize(A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+spSystem.colMapper(6).freeSize(),1);
            fpsi.middleRows(A11_rows, spSystem.colMapper(4).freeSize()) = fmuN;
            fpsi.middleRows(A11_rows + spSystem.colMapper(4).freeSize(), spSystem.colMapper(5).freeSize()) = fmuT;
            fpsi.bottomRows(spSystem.colMapper(6).freeSize()).setZero();

            // reserve memory
            gsVector<int> nonZerosPerCol(A11.outerSize());
            for(index_t i=0; i<A11.outerSize();i++)
                nonZerosPerCol(i) = cast<double,int>((A11.innerVector(i).nonZeros()+Ln1.rows()+Lt1.rows())*1.333); //Todo: improve performace
            A11.reserve(nonZerosPerCol);

            nonZerosPerCol.resize(A10.outerSize());
            for(index_t i=0; i<A10.outerSize();i++)
                nonZerosPerCol(i) = cast<double,int>((A10.innerVector(i).nonZeros()+Ln0.rows()+Lt0.rows())*1.333); //Todo: improve performace
            A10.reserve(nonZerosPerCol);

            for (int k=0; k<Ln1.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(Ln1,k); it; ++it)
                {
                    A11.insert(A11_rows+it.row(),it.col()) = it.value();
                    A11.insert(it.col(),A11_cols+it.row()) = it.value();
                }

            for (int k=0; k<Ln0.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(Ln0,k); it; ++it)
                {
                    A10.insert(A10_rows+it.row(),it.col()) = it.value();
                }


            for (int k=0; k<Lt1.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(Lt1,k); it; ++it)
                {
                    A11.insert(A11_rows+spSystem.colMapper(4).freeSize()+it.row(),it.col()) = it.value();
                    A11.insert(it.col(),A11_cols+spSystem.colMapper(4).freeSize()+it.row()) = it.value();
                }

            for (int k=0; k<Lt0.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(Lt0,k); it; ++it)
                {
                    A10.insert(A10_rows+spSystem.colMapper(4).freeSize()+it.row(),it.col()) = it.value();
                }

            for (int k=0; k<LnMean.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(LnMean,k); it; ++it)
                {
                    A11.insert(A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+it.row(),A11_cols+it.col()) = it.value();
                    A11.insert(A11_cols+it.col(),A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+it.row()) = it.value();
                }

            for (int k=0; k<LtMean.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(LtMean,k); it; ++it)
                {
                    A11.insert(A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+it.row(),A11_cols+spSystem.colMapper(4).freeSize()+it.col()) = it.value();
                    A11.insert(A11_cols+spSystem.colMapper(4).freeSize()+it.col(),A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+it.row()) = it.value();
                }


            // Initialize direct solver

#if (defined(GISMO_WITH_PARDISO))
            gsSparseSolver<>::PardisoLLT solverB;
            //gsSparseSolver<>::PardisoLDLT solverA11;
            gsSparseSolver<>::PardisoLU solverA11;
#else
            gsSparseSolver<>::LU solverB;
            //gsSparseSolver<>::LU solverBaddN;
            gsSparseSolver<>::LU solverA11;
#endif

            // Solve

            watch.restart();
            // p problem
            solverB.compute(-B); // -B is positive definite matrix -> solver
            solVector_p = solverB.solve(-fv); // -fv because solverB = -B (pos. def.)

            // phi problem
            fpsi = fpsi - A10*solVector_p;

            solverA11.compute(A11);
            solVector_phi_lambda = solverA11.solve(fpsi);
            solVector_phi = solVector_phi_lambda.block(0,0,spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize(),1);
            solVector_lambdaN = solVector_phi_lambda.block(spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize(),0,spSystem.colMapper(4).freeSize(),1);
            solVector_lambdaT = solVector_phi_lambda.block(spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize()+ spSystem.colMapper(4).freeSize(),0, spSystem.colMapper(5).freeSize(),1);
            solVector_lambda = solVector_phi_lambda.block(spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize(),0,spSystem.colMapper(4).freeSize() + spSystem.colMapper(5).freeSize(),1);
            solVector_lambdaMean = solVector_phi_lambda.block(spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize()+spSystem.colMapper(4).freeSize() + spSystem.colMapper(5).freeSize(),0,spSystem.colMapper(6).freeSize(),1);

            // w problem
            fq = fq - A00*solVector_p - A01*solVector_phi - Ln0.transpose()*solVector_lambdaN - Lt0.transpose()*solVector_lambdaT;
            solVector_w = solverB.solve(-fq); // -fq because solverB = -B (pos. def.)

            cout<<"Solving time: "<<watch.stop()<<endl;
        }
        else
        {
#if (defined(GISMO_WITH_PARDISO))
            gsSparseSolver<>::PardisoLU solverS;
#else
            gsSparseSolver<>::LU solverS;
#endif

            solverS.compute(assembler.matrix());
            solVector = solverS.solve(-assembler.rhs());

        }

        /// Construct solution ///
        if(! options.getSwitch("solveWholeSystem"))
        {
            index_t numDofs= spSystem.rows();
            solVector.resize(numDofs,1);
            solVector << solVector_p, solVector_phi, solVector_w, solVector_lambdaN, solVector_lambdaT, solVector_lambdaMean;
        }

        gsVector<index_t> unk_p(1);
        unk_p << 0;
        assembler.constructSolution(solVector, fineSol_p, unk_p);
        gsVector<index_t> unk_phi(2);
        unk_phi << 1,2;
        assembler.constructSolution(solVector, fineSol_phi, unk_phi);
        gsVector<index_t> unk_w(1);
        unk_w << 3;
        assembler.constructSolution(solVector, fineSol_w, unk_w);

        gsVector<index_t> unk_M(3);
        unk_M << 0,1,2;
        assembler.constructSolution(solVector, fineSol_M, unk_M);
    }

    cout<<endl;


    /// Compute solution ///

    for(index_t numHRefine=numHRefine_start; numHRefine<=numHRefine_end; ++numHRefine)
    {

        cout<<"numHRefine "<<numHRefine<<endl;
        cout<<"degree basis "<<basis[0].degree(0)<<" "<<basis[0].degree(1)<<endl;

        for(size_t i = 0; i < basis.nBases(); ++i)
            gsInfo<<"Approximation basis (patch "<<i<<"): " << basis.piece(i) << "\n";
        for(size_t i = 0; i < basis.nBases(); ++i)
            gsInfo<<"Approximation basis_degMinus1 (patch "<<i<<"): " << basis_degMinus1.piece(i) << "\n";

        /// Assemble and Solve ///

        // construct bases

        std::vector< gsMultiBasis<> >  bases;

        bases.push_back( basis);// basis for p, w
        bases.push_back( basis_degMinus1);// basis for lambda

        // basis for phi
        /*
        gsMultiBasis<> basis_copy( basis);
        basis_copy.degreeIncrease(1);
        bases.push_back( basis_copy);
        */


        gsPlateMixedAssembler<real_t> assembler(pde, bases, options);

        // Assemble
        watch.restart();
        assembler.assemble();
        cout<<"Assembling time: "<<watch.stop()<<endl;
        data_timing(indRefine, 0) = watch.stop();

        // Build decoupled system
        gsSparseSystem<> spSystem = assembler.system();

        if(! options.getSwitch("solveWholeSystem"))
        {
            watch.restart();
            cout<<"Dofs "<<spSystem.cols()<<" = "<<spSystem.colMapper(0).freeSize()<<" + "<<spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize()<<" + "<<spSystem.colMapper(0).freeSize()<<endl;
            gsVector<index_t> blocks(7);
            blocks<< 0,1,1,2,3,4,5;
            gsSparseMatrix<>::BlockView matrixBlockView = spSystem.blockView(6, 6, blocks, blocks);
            gsMatrix<>::BlockView rhsBlockView = spSystem.blockViewRhs(6, blocks);

            gsSparseMatrix<> A00 = matrixBlockView(0,0);
            gsSparseMatrix<> A10 = matrixBlockView(1,0);
            gsSparseMatrix<> A01 = matrixBlockView(0,1);
            gsSparseMatrix<> A11 = matrixBlockView(1,1);
            gsSparseMatrix<> B = matrixBlockView(2,0);
            gsMatrix<> fq = rhsBlockView(0,0);
            gsMatrix<> fpsi = rhsBlockView(1,0);
            gsMatrix<> fv = rhsBlockView(2,0);

            gsSparseMatrix<> Ln0 = matrixBlockView(3,0);
            gsSparseMatrix<> Ln1 = matrixBlockView(3,1);
            gsSparseMatrix<> Lt0 = matrixBlockView(4,0);
            gsSparseMatrix<> Lt1 = matrixBlockView(4,1);
            gsMatrix<> fmuN = rhsBlockView(3,0);
            gsMatrix<> fmuT = rhsBlockView(4,0);

            gsSparseMatrix<> LnMean = matrixBlockView(5,3);
            gsSparseMatrix<> LtMean = matrixBlockView(5,4);


            /*
            cout<<"fpsi \n"<<fpsi<<endl;
            cout<<"fv \n"<<fv<<endl;
            cout<<"fmuN \n"<<fmuN<<endl;
            cout<<"fmuT \n"<<fmuT<<endl;
            */




            // add Ln0, Lt0 to A10 and Ln1, Lt1 to A11
            index_t A11_rows= A11.rows();
            index_t A11_cols= A11.cols();
            index_t A10_rows= A10.rows();
            index_t A10_cols= A10.cols();
            A11.conservativeResize(A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+spSystem.colMapper(6).freeSize(), A11_cols+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+spSystem.colMapper(6).freeSize());
            A10.conservativeResize(A10_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+spSystem.colMapper(6).freeSize(), A10_cols);
            fpsi.conservativeResize(A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+spSystem.colMapper(6).freeSize(),1);
            fpsi.middleRows(A11_rows, spSystem.colMapper(4).freeSize()) = fmuN;
            fpsi.middleRows(A11_rows + spSystem.colMapper(4).freeSize(), spSystem.colMapper(5).freeSize()) = fmuT;
            fpsi.bottomRows(spSystem.colMapper(6).freeSize()).setZero();

            //cout<<"fpsi \n"<<fpsi<<endl;


            // reserve memory
            gsVector<int> nonZerosPerCol(A11.outerSize());
            for(index_t i=0; i<A11.outerSize();i++)
                nonZerosPerCol(i) = cast<double,int>((A11.innerVector(i).nonZeros()+Ln1.rows()+Lt1.rows())*1.333); //Todo: improve performace
            A11.reserve(nonZerosPerCol);

            nonZerosPerCol.resize(A10.outerSize());
            for(index_t i=0; i<A10.outerSize();i++)
                nonZerosPerCol(i) = cast<double,int>((A10.innerVector(i).nonZeros()+Ln0.rows()+Lt0.rows())*1.333); //Todo: improve performace
            A10.reserve(nonZerosPerCol);


            for (int k=0; k<Ln1.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(Ln1,k); it; ++it)
                {
                    A11.insert(A11_rows+it.row(),it.col()) = it.value();
                    A11.insert(it.col(),A11_cols+it.row()) = it.value();
                }

            for (int k=0; k<Ln0.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(Ln0,k); it; ++it)
                {
                    A10.insert(A10_rows+it.row(),it.col()) = it.value();
                }

            for (int k=0; k<Lt1.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(Lt1,k); it; ++it)
                {
                    A11.insert(A11_rows+spSystem.colMapper(4).freeSize()+it.row(),it.col()) = it.value();
                    A11.insert(it.col(),A11_cols+spSystem.colMapper(4).freeSize()+it.row()) = it.value();
                }

            for (int k=0; k<Lt0.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(Lt0,k); it; ++it)
                {
                    A10.insert(A10_rows+spSystem.colMapper(4).freeSize()+it.row(),it.col()) = it.value();
                }

            for (int k=0; k<LnMean.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(LnMean,k); it; ++it)
                {
                    A11.insert(A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+it.row(),A11_cols+it.col()) = it.value();
                    A11.insert(A11_cols+it.col(),A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+it.row()) = it.value();
                }

            for (int k=0; k<LtMean.outerSize(); ++k)
                for (gsSparseMatrix<>::InnerIterator it(LtMean,k); it; ++it)
                {
                    A11.insert(A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+it.row(),A11_cols+spSystem.colMapper(4).freeSize()+it.col()) = it.value();
                    A11.insert(A11_cols+spSystem.colMapper(4).freeSize()+it.col(),A11_rows+spSystem.colMapper(4).freeSize()+spSystem.colMapper(5).freeSize()+it.row()) = it.value();
                }


            cout<<"System building time: "<<watch.stop()<<endl;

            /*
                // Compute singular values C

                gsSparseMatrix<> CC_tr =C_boundary*C_boundary.transpose();
                gsMatrix<> CC_trd = CC_tr;

                //cout<<"CC_trd"<<endl<<CdCd_tr<<endl;

                cout<<"singvals Cd_boundary"<<endl;
                //cout<<CdCd_tr.eigenvalues()<<endl;

                gsMatrix<> singval = CC_trd.selfadjointView<Lower>().eigenvalues();

                if(singval.rows()<10)
                {
                    cout<<singval<<endl;
                }
                else
                {
                    for(index_t i=0; i<10; i++)
                        cout<<sqrt(singval(i,0))<<endl;
                }
                */



            // Compute eigenvalues
            if(computeEValsA11)
            {
                cout<<"eigs A11: \n"<<A11.toDense().eigenvalues()<<"\n";
            }


            // Initialize direct solver

#if (defined(GISMO_WITH_PARDISO))
            gsSparseSolver<>::PardisoLLT solverB;
            //gsSparseSolver<>::PardisoLDLT solverA11;
            gsSparseSolver<>::PardisoLU solverA11;
#else

            gsSparseSolver<>::LU solverB;
            //gsSparseSolver<>::LU solverBaddN;
            gsSparseSolver<>::LU solverA11;
#endif

            // Solve
            watch.restart();

            // p problem
            //gsInfo<<"B: "<<B.toDense()<<"\n";
            solverB.compute(-B); // -B is positive definite matrix -> solver
            solVector_p = solverB.solve(-fv); // -fv because solverB = -B (pos. def.)

            // phi problem
            fpsi = fpsi - A10*solVector_p;

            solverA11.compute(A11);
            cout<<"LU factorization succeded "<<solverA11.succeed()<<endl;
            solVector_phi_lambda = solverA11.solve(fpsi);
            solVector_phi = solVector_phi_lambda.block(0,0,spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize(),1);
            solVector_lambdaN = solVector_phi_lambda.block(spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize(),0,spSystem.colMapper(4).freeSize(),1);
            solVector_lambdaT = solVector_phi_lambda.block(spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize()+ spSystem.colMapper(4).freeSize(),0, spSystem.colMapper(5).freeSize(),1);
            solVector_lambda = solVector_phi_lambda.block(spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize(),0,spSystem.colMapper(4).freeSize() + spSystem.colMapper(5).freeSize(),1);
            solVector_lambdaMean = solVector_phi_lambda.block(spSystem.colMapper(1).freeSize() + spSystem.colMapper(2).freeSize()+spSystem.colMapper(4).freeSize() + spSystem.colMapper(5).freeSize(),0,spSystem.colMapper(6).freeSize(),1);

            // w problem
            fq = fq - A00*solVector_p - A01*solVector_phi - Ln0.transpose()*solVector_lambdaN - Lt0.transpose()*solVector_lambdaT;
            solVector_w = solverB.solve(-fq); // -fq because solverB = -B (pos. def.)

            cout<<"Solving time: "<<watch.stop()<<endl;
        }
        else
        {
#if (defined(GISMO_WITH_PARDISO))
            gsSparseSolver<>::PardisoLU solverS;
#else
            gsSparseSolver<>::LU solverS;
#endif

            solverS.compute(assembler.matrix());
            solVector = solverS.solve(-assembler.rhs());
        }

        if(outputMatrix)
        {
            cout<<assembler.matrix().toDense()<<endl;
        }

        if(outputEvals)
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs(assembler.matrix().toDense());
            cout<<"eigs: \n"<<eigs.eigenvalues()<<"\n";
        }

        if(outputBlock)
        {

            //gsSparseSystem<> spSystem = assembler.system();
            gsVector<index_t> blocks_L(7);
            blocks_L<< 0,0,0,1,2,2,3;
            gsSparseMatrix<>::BlockView matrixBlockView_L = spSystem.blockView(4, 4, blocks_L, blocks_L);
            gsSparseMatrix<> L = matrixBlockView_L(2,0);

            //gsSparseSystem<> spSystem = assembler.system();
            gsVector<index_t> blocks_LnLt(7);
            blocks_LnLt<< 0,0,0,1,2,3,4;
            gsSparseMatrix<>::BlockView matrixBlockView_LnLt = spSystem.blockView(5, 5, blocks_LnLt, blocks_LnLt);
            gsSparseMatrix<> Ln = matrixBlockView_LnLt(2,0);
            gsSparseMatrix<> Lt = matrixBlockView_LnLt(3,0);

            /*
            cout<<"Ln \n"<<Ln.toDense()<<endl;
            cout<<"Lt \n"<<Lt.toDense()<<endl;
            */

            gsVector<index_t> blocks_L0L1(7);
            blocks_L0L1<< 0,1,1,2,3,3,4;
            gsSparseMatrix<>::BlockView matrixBlockView_L0L1 = spSystem.blockView(5, 5, blocks_L0L1, blocks_L0L1);
            gsSparseMatrix<> L0 = matrixBlockView_L0L1(3,0);
            gsSparseMatrix<> L1 = matrixBlockView_L0L1(3,1);

            //gsSparseSystem<> spSystem = assembler.system();
            gsVector<index_t> blocks(7);
            blocks<< 0,1,1,2,3,4,5;
            gsSparseMatrix<>::BlockView matrixBlockView = spSystem.blockView(6, 6, blocks, blocks);
            gsSparseMatrix<> Ln0 = matrixBlockView(3,0);
            gsSparseMatrix<> Ln1 = matrixBlockView(3,1);
            gsSparseMatrix<> Lt0 = matrixBlockView(4,0);
            gsSparseMatrix<> Lt1 = matrixBlockView(4,1);

            /*
            cout<<"Ln0 \n"<<Ln0.toDense()<<endl;
            cout<<"Ln1 \n"<<Ln1.toDense()<<endl;
            cout<<"Lt0 \n"<<Lt0.toDense()<<endl;
            cout<<"Lt1 \n"<<Lt1.toDense()<<endl;
            */

            Eigen::JacobiSVD<Eigen::MatrixXd> svd_L(L.toDense());
            cout<<"L:"<< " rows "<<L.rows()<<" cols "<<L.cols()<<" rank "<<svd_L.rank()<<"\n";

            if(L0.cols()!=0)
            {
                Eigen::JacobiSVD<Eigen::MatrixXd> svd_L0(L0.toDense());
                cout<<"L0:"<< " rows "<<L0.rows()<<" cols "<<L0.cols()<<" rank "<<svd_L0.rank()<<"\n";
            }
            else
            {
                cout<<"L0:"<< " rows "<<L0.rows()<<" cols "<<L0.cols()<<" rank "<<0<<"\n";
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> svd_L1(L1.toDense());
            cout<<"L1:"<< " rows "<<L1.rows()<<" cols "<<L1.cols()<<" rank "<<svd_L1.rank()<<"\n";

            Eigen::JacobiSVD<Eigen::MatrixXd> svd_Ln(Ln.toDense());
            cout<<"Ln:"<< " rows "<<Ln.rows()<<" cols "<<Ln.cols()<<" rank "<<svd_Ln.rank()<<"\n";
            if(Lt.rows()!=0)
            {
                Eigen::JacobiSVD<Eigen::MatrixXd> svd_Lt(Lt.toDense());
                cout<<"Lt:"<< " rows "<<Lt.rows()<<" cols "<<Lt.cols()<<" rank "<<svd_Lt.rank()<<"\n";
            }
            else
            {
                cout<<"Lt:"<< " rows "<<Lt.rows()<<" cols "<<Lt.cols()<<" rank "<<0<<"\n";
            }

            if(L0.cols()!=0)
            {
                Eigen::JacobiSVD<Eigen::MatrixXd> svd_Ln0(Ln0.toDense());
                cout<<"Ln0:"<< " rows "<<Ln0.rows()<<" cols "<<Ln0.cols()<<" rank "<<svd_Ln0.rank()<<"\n";
            }
            else
            {
                cout<<"Ln0:"<< " rows "<<Ln0.rows()<<" cols "<<Ln0.cols()<<" rank "<<0<<"\n";
            }

            Eigen::JacobiSVD<Eigen::MatrixXd> svd_Ln1(Ln1.toDense());
            cout<<"Ln1:"<< " rows "<<Ln1.rows()<<" cols "<<Ln1.cols()<<" rank "<<svd_Ln1.rank()<<"\n";

            if(Lt.rows()!=0)
            {
                Eigen::JacobiSVD<Eigen::MatrixXd> svd_Lt0(Lt0.toDense());
                cout<<"Lt0:"<< " rows "<<Lt0.rows()<<" cols "<<Lt0.cols()<<" rank "<<svd_Lt0.rank()<<"\n";
                Eigen::JacobiSVD<Eigen::MatrixXd> svd_Lt1(Lt1.toDense());
                cout<<"Lt1:"<< " rows "<<Lt1.rows()<<" cols "<<Lt1.cols()<<" rank "<<svd_Lt1.rank()<<"\n";
            }
            else
            {
                cout<<"Lt0:"<< " rows "<<Lt0.rows()<<" cols "<<Lt0.cols()<<" rank "<<0<<"\n";
                cout<<"Lt1:"<< " rows "<<Lt1.rows()<<" cols "<<Lt1.cols()<<" rank "<<0<<"\n";
            }

        }

        /// Construct solution ///

        // Construct solution (as gsMultiPatch)
        if(! options.getSwitch("solveWholeSystem"))
        {
            index_t numDofs= spSystem.rows();
            solVector.resize(numDofs,1);
            solVector << solVector_p, solVector_phi, solVector_w, solVector_lambdaN, solVector_lambdaT, solVector_lambdaMean;
        }

        gsVector<index_t> unk_p(1);
        unk_p << 0;
        assembler.constructSolution(solVector, sol_p, unk_p);
        gsVector<index_t> unk_phi(2);
        unk_phi << 1,2;
        assembler.constructSolution(solVector, sol_phi, unk_phi);
        gsVector<index_t> unk_w(1);
        unk_w << 3;
        assembler.constructSolution(solVector, sol_w, unk_w);

        gsVector<index_t> unk_M(3);
        unk_M << 0,1,2;
        assembler.constructSolution(solVector, sol_M, unk_M);

        // Construct solution fields
        gsField<> solField_w(mpptr, sol_w);
        gsField<> solField_p(mpptr, sol_p);
        gsField<> solField_phi(mpptr, sol_phi);
        gsField<> solField_M(mpptr, sol_M);

        if(analytSol && !convRates_p_phi)
        {
            fineSol_w = sol_w;
            fineSol_p = sol_p;
            fineSol_phi = sol_phi;
        }



        /// Compute point values of solution ///
        /*

        gsMatrix<> point_param(2,1);
        point_param<<0.5, 0.5;
        //point_param<<0, 0;
        gsMatrix<> point_physical(solField_w.point(point_param,0));
        gsMatrix<> val(solField_w.value(point_param,0));
        cout<<"point_physical"<<endl;
        cout<<point_physical<<endl;
        cout<<"w"<<endl;
        cout<<val<<endl;
        */



        /// Compute discretization errors ///

        data(indRefine, 0) = numHRefine;
        //data(indRefine, 1) = assembler.numDofs();
        data(indRefine, 1) = spSystem.colMapper(0).freeSize();

        if(analytSol) {
            gsNormL2<real_t> L2Norm(solField_w, analytSol_w);
            gsSeminormH1<real_t> H1Seminorm(solField_w, analytSol_w);
            gsNormL2M<real_t> L2NormM(solField_M, analytSol_M);

            real_t errorL2 = L2Norm.compute();
            real_t errorH1 = H1Seminorm.compute();

            data(indRefine, 2) = errorL2;
            data(indRefine, 3) = sqrt(errorL2 * errorL2 + errorH1 * errorH1);
            data(indRefine, 4) = L2NormM.compute();
        }
        else
        {
            for (size_t pn=0; pn < mpptr.nPatches(); ++pn )
            {
                gsMultiPatch<> mpPatch(mpptr.patch(pn));
                gsMultiPatch<> fineSolPatch_w(fineSol_w[pn]);
                gsField<> fineFieldPatch_w(mpPatch, fineSolPatch_w);

                gsNormL2<real_t> L2Norm(fineFieldPatch_w, sol_w[pn], true);
                gsSeminormH1<real_t> H1Seminorm(fineFieldPatch_w, sol_w[pn], true);

                gsMultiPatch<> fineSolPatch_M(fineSol_M[pn]);
                gsField<> fineFieldPatch_M(mpPatch, fineSolPatch_M);

                gsNormL2M<real_t> L2NormM(fineFieldPatch_M, sol_M[pn], true);

                real_t errorL2 = L2Norm.compute();
                real_t errorH1 = H1Seminorm.compute();
                real_t errorL2M = L2NormM.compute();

                data(indRefine, 2) += errorL2*errorL2;
                data(indRefine, 3) += errorL2*errorL2 + errorH1*errorH1;
                data(indRefine, 4) += errorL2M*errorL2M;
            }

            data(indRefine, 2) = sqrt(data(indRefine, 2));
            data(indRefine, 3) = sqrt(data(indRefine, 3));
            data(indRefine, 4) = sqrt(data(indRefine, 4));


        }

        if(convRates_p_phi)
        {
            for (size_t pn=0; pn < mpptr.nPatches(); ++pn )
            {
                gsMultiPatch<> mpPatch(mpptr.patch(pn));
                gsMultiPatch<> fineSolPatch_p(fineSol_p[pn]);
                gsField<> fineFieldPatch_p(mpPatch, fineSolPatch_p);

                gsMultiPatch<> fineSolPatch_phi(fineSol_phi[pn]);
                gsField<> fineFieldPatch_phi(mpPatch, fineSolPatch_phi);

                gsNormL2<real_t> L2Norm_p(fineFieldPatch_p, sol_p[pn], true);
                gsSeminormH1<real_t> H1Seminorm_p(fineFieldPatch_p, sol_p[pn], true);

                gsNormL2<real_t> L2Norm_phi(fineFieldPatch_phi, sol_phi[pn], true);
                gsSeminormH1<real_t> H1Seminorm_phi(fineFieldPatch_phi, sol_phi[pn], true);

                real_t errorL2_p = L2Norm_p.compute();
                real_t errorH1_p = H1Seminorm_p.compute();

                real_t errorL2_phi = L2Norm_phi.compute();
                real_t errorH1_phi = H1Seminorm_phi.compute();

                data(indRefine, 5) += errorL2_p*errorL2_p;
                data(indRefine, 6) += errorL2_p*errorL2_p + errorH1_p*errorH1_p;
                data(indRefine, 7) += errorL2_phi*errorL2_phi;
                data(indRefine, 8) += errorL2_phi*errorL2_phi + errorH1_phi*errorH1_phi;
            }

            data(indRefine, 5) = sqrt(data(indRefine, 5));
            data(indRefine, 6) = sqrt(data(indRefine, 6));
            data(indRefine, 7) = sqrt(data(indRefine, 7));
            data(indRefine, 8) = sqrt(data(indRefine, 8));
        }


        /// Plot ///

        //gsField<> fineField_w(mpptr, fineSol_w);
        //gsField<> fineField_p(mpptr, fineSol_p);
        //gsField<> fineField_phi(mpptr, fineSol_phi);

        gsField<> errorField_wExact = gsFieldCreator<>::absError(solField_w, analytSol_w);
        gsField<> errorField_phiFine = gsFieldCreator<>::absError(solField_phi, fineSol_phi[0], true);
        gsField<> errorField_wFine = gsFieldCreator<>::absError(solField_w, fineSol_w[0], true);


        gsField<> analytField_w( mpptr,analytSol_w, false );

        // Plot solution in paraview
        if (plot && numHRefine == numHRefine_end)
        {
            // Write solution to paraview file
            std::cout<<"Plotting in Paraview...\n";
            gsWriteParaview<>( mpptr, "geometry", 1000, true, true);
            //gsWriteParaview<>( basis.front(), "basis", 1000);
            gsWriteParaview<>( solField_w, "sol_w", 1000);
            gsWriteParaview<>( solField_p, "sol_p", 1000);
            gsWriteParaview<>( solField_phi, "sol_phi", 1000);

            gsWriteParaview<>( errorField_wExact, "error_wExact", 1000);
            gsWriteParaview<>( errorField_phiFine, "error_phiFine", 1000);
            gsWriteParaview<>( errorField_wFine, "error_wFine", 1000);

            gsWriteParaview<>( analytField_w, "analytSol_w", 1000);

            // Run paraview on exit
            result = system("paraview sol_w.pvd &");
        }


        ++indRefine;
        if(keepC0)
            basis.uniformRefine(1, basis.maxDegree(0)); // multipatch: basis of patch 0 is considered
        else
        {
            basis.uniformRefine();
            basis_degMinus1.uniformRefine();
        }

    } //end refinement loop


    /// Compute norms for relative error ///

    if(relative_error)
    {
        real_t ValL2Norm_w = 0;
        real_t ValH1Norm_w = 0;
        real_t ValL2Norm_M = 0;

        real_t ValL2Norm_p = 0;
        real_t ValH1Norm_p = 0;
        real_t ValL2Norm_phi = 0;
        real_t ValH1Norm_phi = 0;

        // w, M
        if(analytSol)
        {

            const int sz = basis[0].size();
            index_t dim = 1;
            gsMatrix<> coeffs;
            coeffs.setZero(sz,dim);
            gsMultiPatch<> sol_zero;
            sol_zero.addPatch(basis[0].makeGeometry( give(coeffs) ));
            gsField<> solField_zero(mpptr, sol_zero);

            dim = 4;
            gsMatrix<> coeffs_4;
            coeffs_4.setZero(sz,dim);
            gsMultiPatch<> sol_zero_4;
            sol_zero_4.addPatch(basis[0].makeGeometry( give(coeffs_4) ));
            gsField<> solField_zero_4(mpptr, sol_zero_4);

            gsNormL2<real_t> L2Norm(solField_zero, analytSol_w);
            gsSeminormH1<real_t> H1Seminorm(solField_zero, analytSol_w);
            gsNormL2<real_t> L2NormM(solField_zero_4, analytSol_M);

            real_t errorL2 = L2Norm.compute();
            real_t errorH1 = H1Seminorm.compute();

            ValL2Norm_w = errorL2;
            ValH1Norm_w = sqrt(errorL2*errorL2 + errorH1*errorH1);
            ValL2Norm_M = L2NormM.compute();
        }

        else
        {
            for (size_t pn=0; pn < mpptr.nPatches(); ++pn )
            {
                gsMultiPatch<> mpPatch(mpptr.patch(pn));
                gsMultiPatch<> fineSolPatch_w(fineSol_w[pn]);
                gsField<> fineFieldPatch_w(mpPatch, fineSolPatch_w);

                gsNormL2<real_t> L2Norm(fineFieldPatch_w);
                gsSeminormH1<real_t> H1Seminorm(fineFieldPatch_w);

                const gsConstantFunction<> zeroFunction(gsVector<>::Zero(4), 2);
                gsMultiPatch<> fineSolPatch_M(fineSol_M[pn]);
                gsField<> fineFieldPatch_M(mpPatch, fineSolPatch_M);

                gsNormL2M<real_t> L2NormM(fineFieldPatch_M, zeroFunction);

                real_t errorL2 = L2Norm.compute();
                real_t errorH1 = H1Seminorm.compute();
                real_t errorL2M = L2NormM.compute();

                ValL2Norm_w += errorL2*errorL2;
                ValH1Norm_w += errorL2*errorL2 + errorH1*errorH1;
                ValL2Norm_M += errorL2M*errorL2M;
            }

            ValL2Norm_w = sqrt(ValL2Norm_w);
            ValH1Norm_w = sqrt(ValH1Norm_w);
            ValL2Norm_M = sqrt(ValL2Norm_M);
        }


        // p, phi
        if(convRates_p_phi)
        {

            for (size_t pn=0; pn < mpptr.nPatches(); ++pn )
            {
                gsMultiPatch<> mpPatch(mpptr.patch(pn));
                gsMultiPatch<> fineSolPatch_p(fineSol_p[pn]);
                gsField<> fineFieldPatch_p(mpPatch, fineSolPatch_p);

                gsMultiPatch<> fineSolPatch_phi(fineSol_phi[pn]);
                gsField<> fineFieldPatch_phi(mpPatch, fineSolPatch_phi);

                gsNormL2<real_t> L2Norm_p(fineFieldPatch_p);
                gsSeminormH1<real_t> H1Seminorm_p(fineFieldPatch_p);

                gsNormL2<real_t> L2Norm_phi(fineFieldPatch_phi);
                gsSeminormH1<real_t> H1Seminorm_phi(fineFieldPatch_phi);

                real_t errorL2_p = L2Norm_p.compute();
                real_t errorH1_p = H1Seminorm_p.compute();
                real_t errorL2_phi = L2Norm_phi.compute();
                real_t errorH1_phi = H1Seminorm_phi.compute();

                ValL2Norm_p += errorL2_p*errorL2_p;
                ValH1Norm_p += errorL2_p*errorL2_p + errorH1_p*errorH1_p;
                ValL2Norm_phi += errorL2_phi*errorL2_phi;
                ValH1Norm_phi += errorL2_phi*errorL2_phi + errorH1_phi*errorH1_phi;
            }

            ValL2Norm_p = sqrt(ValL2Norm_p);
            ValH1Norm_p = sqrt(ValH1Norm_p);
            ValL2Norm_phi = sqrt(ValL2Norm_phi);
            ValH1Norm_phi = sqrt(ValH1Norm_phi);
        }

        cout<<"ValL2Norm_w "<<ValL2Norm_w<<endl;
        cout<<"ValH1Norm_w "<<ValH1Norm_w<<endl;
        cout<<"ValL2Norm_M "<<ValL2Norm_M<<endl;
        cout<<"ValL2Norm_p "<<ValL2Norm_p<<endl;
        cout<<"ValH1Norm_p "<<ValH1Norm_p<<endl;
        cout<<"ValL2Norm_phi "<<ValL2Norm_phi<<endl;
        cout<<"ValH1Norm_phi "<<ValH1Norm_phi<<endl;


        for (index_t i = 0; i<sizeError; ++i)
        {
            data(i,2)/=ValL2Norm_w;
            data(i,3)/=ValH1Norm_w;
            data(i,4)/=ValL2Norm_M;
            data(i,5)/=ValL2Norm_p;
            data(i,6)/=ValH1Norm_p;
            data(i,7)/=ValL2Norm_phi;
            data(i,8)/=ValH1Norm_phi;

        }

    }



    /// Generate table with discretization errors  ///

    cout<<std::setw (7)<<"L"<<" | " <<std::setw (7)<<"dof"<<" | "<<std::setw (13)<<"L2 error w" <<std::setw (13)<<" "<<" | "<<std::setw (13)<<"H1 error w"<<std::setw (13)<<" "<<" | "<<std::setw (13)<<"L2 error M"<<std::setw (13)<<" "<<" | "<<std::setw (13)<<"L2 error p"<<std::setw (13)<<" "<<" | "<<std::setw (13)<<"H1 error p"<<std::setw (13)<<" "<<" | "<<std::setw (13)<<"L2 error phi"<<std::setw (13)<<" "<<" | "<<std::setw (13)<<"H1 error phi"<<std::setw (13)<<" "<<" | "<<endl;
    for (index_t i = 0; i<sizeError; ++i)
    {
        if(i > 0)
        {
            cout<<std::setw (7)<< data(i,0)<<" | "<<std::setw (7)<< data(i,1)<<" | "<<std::setw (13)<<std::scientific<< data(i,2);
            cout.unsetf(std::ios_base::floatfield);
            cout<<std::setw (13)<< math::log(data(i-1,2)/data(i,2))/std::log(2.0) <<" | "<<std::setw (13)<<std::scientific<< data(i,3);
            cout.unsetf(std::ios_base::floatfield);
            cout<<std::setw (13)<<math::log(data(i-1,3)/data(i,3))/std::log(2.0)<<" | "<<std::setw (13)<<std::scientific<< data(i,4);
            cout.unsetf(std::ios_base::floatfield);
            cout<<std::setw (13)<<math::log(data(i-1,4)/data(i,4))/std::log(2.0)<<" | "<<std::setw (13)<<std::scientific<< data(i,5);
            cout.unsetf(std::ios_base::floatfield);
            cout<<std::setw (13)<<math::log(data(i-1,5)/data(i,5))/std::log(2.0)<<" | "<<std::setw (13)<<std::scientific<< data(i,6);
            cout.unsetf(std::ios_base::floatfield);
            cout<<std::setw (13)<<math::log(data(i-1,6)/data(i,6))/std::log(2.0)<<" | "<<std::setw (13)<<std::scientific<< data(i,7);
            cout.unsetf(std::ios_base::floatfield);
            cout<<std::setw (13)<<math::log(data(i-1,7)/data(i,7))/std::log(2.0)<<" | "<<std::setw (13)<<std::scientific<< data(i,8);
            cout.unsetf(std::ios_base::floatfield);
            cout<<std::setw (13)<<math::log(data(i-1,8)/data(i,8))/std::log(2.0)<<" | "<<endl;
        }
        else
        {
            cout<<std::setw (7)<< data(i,0)<<" | "<<std::setw (7)<< data(i,1)<<" | "<<std::setw (13)<<std::scientific<< data(i,2)<<std::setw (13)<<""<<" | "<<std::setw (13)<<std::scientific << data(i,3)<<std::setw (13)<<""<<" | "<<std::setw (13)<<std::scientific << data(i,4)<<std::setw (13)<<""<<" | "<<std::setw (13)<<std::scientific << data(i,5)<<std::setw (13)<<""<<" | "<<std::setw (13)<<std::scientific << data(i,6)<<std::setw (13)<<""<<" | "<<std::setw (13)<<std::scientific << data(i,7)<<std::setw (13)<<""<<" | "<<std::setw (13)<<std::scientific << data(i,8)<<std::setw (13)<<""<<" | "<<endl;
            cout.unsetf(std::ios_base::floatfield);
        }
    }

    cout<<endl<<endl;
    cout<<std::setw (7)<<"L"<<" | " <<std::setw (7)<<"dof"<<" | "<<std::setw (13)<<"assemb. time" <<std::setw (13)<<" "<<" | "<<endl;
    for (index_t i = 0; i<sizeError; ++i)
    {
        if(i > 0)
        {
            cout<<std::setw (7)<< data(i,0)<<" | "<<std::setw (7)<< data(i,1)<<" | "<<std::setw (13)<< data_timing(i,0);
            cout<<std::setw (13)<< math::log(data_timing(i-1,0)/data_timing(i,0))/std::log(2.0) <<" | "<<endl;
        }
        else
        {
            cout<<std::setw (7)<< data(i,0)<<" | "<<std::setw (7)<< data(i,1)<<" | "<<std::setw (13)<< data_timing(i,0)<<std::setw (13)<<""<<" | "<<endl;
        }
    }

    return result;
}




