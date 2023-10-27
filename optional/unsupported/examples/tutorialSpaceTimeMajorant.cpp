/** @file tutorialSpaceTimeMajorant.cpp

    @brief Example for testing the gsSpaceTimeSolver (globally stabilised space-time scheme)
    with adaptive refinement based on THB-splines using functional error estimates for the error
    control and following from them indicators as a refinement criteria.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/
#include <iostream>
#include <gismo.h>
#include <gismo_dev.h>

#include <gsAssembler/gsAdaptiveRefUtils.h>                     // Util with implemetation of the marking procedures
#include <gsAssembler/gsSeminormH1.h>                           // Class implementting H1-seminorm
#include <gsAssembler/gsNormL2.h>                               // Class implementting L2-norm

#include <gsErrorEstimates/gsTestSpaceTimeMajorant.h>           // Class that tests majorant

// Classes targeted to the globally stabilised scheme and the correcponding norms
#include <gsErrorEstimates/gsSpaceTimeAssembler.h>
#include <gsErrorEstimates/gsSpaceTimeNorm.h>                   // || grad_x w ||^2_Q + h || w_t ||^2_Q
#include <gsErrorEstimates/gsSpaceTimeSpaceGradNorm.h>          // || grad_x w ||^2_Q
#include <gsErrorEstimates/gsSpaceTimeSliceNorm.h>              // || w ||^2_{\Sigma_T} + h || grad_x w ||^2_{\Sigma_T}
#include <gsErrorEstimates/gsSpaceTimeSpaceGradSliceNorm.h>     // || grad_x w ||^2_{\Sigma_T}
#include <gsErrorEstimates/gsSpaceTimeSigmaTNorm.h>             // || w ||^2_{\Sigma_T}

#include <gsErrorEstimates/gsNormFields.h>

// Classes related to the error identity reconstruction
#include <gsErrorEstimates/gsErrEstSpaceTimeIdentity.h>         // ||f + Delta_x v - v_t||^2_Q
#include <gsErrorEstimates/gsSpaceTimeSolOperNorm.h>            // || Delta_x w ||^2_Q + || w_t ||^2_Q
#include <gsErrorEstimates/gsSpaceTimeDeltaxNorm.h>             // || Delta_x w ||^2_Q
#include <gsErrorEstimates/gsSpaceTimeDtNorm.h>                 // || w_t ||^2_Q

#include <gsErrorEstimates/gsErrEstSpaceTimeResidual.h>

using namespace gismo;

// S. Matculevich
//
// This is a test example for a illustrating the adaptive refinement procedure implemented for the gsSpaceTimeAssembler
// Flags, parameters, geometry and prescribed exact solution are specified within the main() function

void gsParseCommandLine(int argc, char **argv, bool plot)
{
    //! [Parse command line]
    gsCmdLine cmd("Tutorial on solving a heat eqaution with guaranteed error control using the funcional error estimate.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the gmresSolver", plot);
    cmd.getValues(argc,argv);
}

int main(int argc, char *argv[])
{

    // Define stopwatch to measure the performance of the routines
    gsCPUStopwatch clock, clock_total;
    clock_total.restart();

    //! [Initialize Test Parameters]
    // -------------------------------------------------------------------------------------------------------------- //
    // Define constants and preliminaries
    const int NUM_PATCHES = 1;     // All the geometries are single-patch geometries

    real_t rho = 1.0 / 16.0;    // Parameter for the hierarchical linear solvers
    real_t TOL = 1e-3;          // Relative tolerance for linear solvers

    // Define parameters of I/O
    bool plotToParaview            = false; // Flag indicating whether objects must be plotted in ParaView
    bool saveToFile                = false;  // Flag indicating whether the results should be saved to the files
    bool isAdaptive                = false;  // Flag if the adaptive refinement mode is on
    bool withMajorant              = true;  // Flag if the error control is driven by the majorant
    bool withMajorantOptimization  = true;  // Flag if the majorant should be minimized
    bool withMajorantEquilibration = false; // Flag if the adaptivity should be driven by equilibrated majorant
    bool withMajorantII            = true;  // Flag is the 2nd majorant must be reconstructed

    try { gsParseCommandLine(argc, argv, plotToParaview); } catch (int rv) { return rv; }

    // Define test-case parameters (number and dimension)
    const unsigned
    // 2d examples:
    // exampleNumber(2), d(2);        // 2 example: 2d unit square, u = (1 - x)*x*x*(1 - t)*t
    // exampleNumber(3), d(2);        // 3 example: 2d unit square, u = sin(pi*x)*sin(pi*t)
    // exampleNumber(4), d(2);        // 4 example: 2d unit square, u = sin(6.0*pi*x)*sin(3.0*pi*t)
    // exampleNumber(5), d(2);        // 5 example: 2d unit square, u = cos(x)*exp(t)
    // exampleNumber(6), d(2);        // 6 example: 2d unit square, u = (x^2 - x) * (y^2 - y) * exp(-100 * ((x - 0.8)^2 + (y - 0.05)^2))
    // exampleNumber(7), d(2);        // 7 example: 2d rectangle $(0, 2) x (0, 1)$, u = (2 - x)*x*x*(1 - y)*y
    // exampleNumber(8), d(2);        // 8 example: 2d rectangle $(0, 2) x (0, 1)$, u = (x^2 - 2*x) * (y^2 - y) * exp(-100 * ((x - 1.4)^2 + (y - 0.95)^2))
    // exampleNumber(9), d(2);        // 9 example: 2d rectangle (0, 1) x (0, 2), u = (x^2 + t^2)^(1.0/3.0) * sin(2.0/3.0*atan2(t,x) + pi)
    // exampleNumber(10), d(2);       // 10 example: 2d rectangle (0, 1) x (0, 2), u = cos(x)*exp(t)
    // exampleNumber(11), d(2);       // 11 example: 2d unit square, u = sin(1 / (1/pi/10 + (x^2 + y^2)^(1.0/2.0)))
    // exampleNumber(12), d(2);       // 12 example: 2d unit square, u =
    // exampleNumber(16), d(2);       // 22 example: 2d unit square, u = (x^2 - x)*(y^2 - y)*exp(-100*((x - 0.25)^2 + (y - 0.25)^2))
    // exampleNumber(33), d(2);       // 33 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)
    // exampleNumber(34), d(2);       // 34 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)^0.5
    // exampleNumber(35), d(2);       // 35 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)^1.5
    // exampleNumber(37), d(2);       // 37 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)^(3/4)
    // exampleNumber(38), d(2);       // 38 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)^(2/3)
    // exampleNumber(39), d(2);       // 39 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)^(4/5)
    // exampleNumber(40), d(2);       // 40 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)^(9/10)
    exampleNumber(41), d(2);          // 41 example: 2d [0, 1] x (0, 1), u = sin(pi*x)*abs(1 - y)^(3/4)
    // exampleNumber(42), d(2);       // 42 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y) with C^0 continuaty
    // exampleNumber(43), d(2);       // 43 example: 2d [0, 1] x (0, 1), u = sin(pi*x)*abs(1 - y)
    // exampleNumber(44), d(2);       // 44 example: 2d [0, 1] x (0, 1), u = x^(5/2)(1-x)*t^(3/4)
    // exampleNumber(45), d(2);       // 45 example: 2d [0, 1] x (0, 1), u = x^(5/2)(1-x)*t^(3/4)

    // 3d examples:
    // exampleNumber(20), d(3);     // 20 example: unit cube, u = sin(pi*x)*sin(pi*y)*sin(pi*t),        non-homogeneous BC
    // exampleNumber(25), d(3);     // 25 example: unit cube, u = (1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2,  homogeneous BC
    // exampleNumber(21), d(3);     // 21 example: unit cube, u = cos(x)*exp(y)*sin(pi*t),              non-homogeneous BC
    // exampleNumber(23), d(3);     // 23 example: 2d+1 quater annulus + [0, 1] in time, u = (1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2
    // exampleNumber(24), d(3);     // 24 example: 2d+1 quater annulus + [0, 1] in time, u = sin(pi*x)*sin(pi*y)*sin(pi*z)
    // exampleNumber(26), d(3);     // 26 example: 2d+1 quater annulus + [0, 1] in time, u = (1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2
    // exampleNumber(27), d(3);     // 27 example: square x [0, 2], u = (1 - x)*x^2*(1 - y)*y^2*(2 - z)*z^2,  homogeneous BC
    // exampleNumber(28), d(3);     // 28 example: [0, 2] x [0, 1] x [0, 1], u = (2 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2,  homogeneous BC
    // exampleNumber(29), d(3);     // 29 example: [0, 2] x [0, 3] x [0, 1], u = (2 - x)*x^2*(3 - y)*y^2*(1 - z)*z^2,  homogeneous BC
    // exampleNumber(30), d(3);     // 30 example: 2d+1 quater annulus + [0, 1] in time, u = (x^2 + y^2 - 1)*(x^2 + y^2 - 4)*(1 - z)*z^2,  homogeneous BC
    // exampleNumber(31), d(3);     // 31 example: G-domain + [0, 1] in time, u = (x^2 + y^2 - 1)*(x^2 + y^2 - 4)*exp(-100 * ((z - 0.8)^2)),  homogeneous BC
    // exampleNumber(15), d(3);     // 15 example: 3d l-shape x (0, 2), u = if( y > 0, (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (x^2+y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    // exampleNumber(31), d(3);     // 31 example: 3d unit cube, (x^2 - x)*(y^2 - y)*(z^2 - z)*exp(-100*((x - 0.25)^2 + (y - 0.25)^2 + (z - 0.25)^2))
    // exampleNumber(32), d(3);     // 32 example: 3d, unit square + [0, 2] in time, u = if( y > 0, (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    // exampleNumber(36), d(3);     // 36 example: 3d, unit square + [0, 2] in time, u = if( y > 0, (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )

    // Initialise the class that tests SpaceTimeMajorant
    gsTestSpaceTimeMajorant<d> testSpaceTime(exampleNumber,
                                             isAdaptive,
                                             withMajorant,
                                             withMajorantOptimization,
                                             withMajorantEquilibration,
                                             withMajorantII);

    // Init the degree of basis S^{p, p}, where p = vDegree
    int vDegree(2), m(1), l(1);
    // Init the degree of the basis for the flux: y \in S^{p + k, p + k} + S^{p + k, p + k}, p = vDegree
    // yDegree must be >= 2 to ensure that S^{p + k, p + k} + S^{p + k, p + k} \subset H(\Omega, div)
    int yDegree(vDegree + m);
    int wDegree(vDegree + l);

    // Setting up the refinement strategy:
    // Number of initial uniform and total unif./adapt. refinement steps
    unsigned int numInitUniformRefV(1), numInitUniformRefY(1), numInitUniformRefW(1), totalRef(5);

    // Initialising the marking strategy
    MarkingStrategy adaptRefCrit(GARU); // with alternatives GARU, PUCA, and BULK
    real_t markingParamTheta(0.4);  // parameter theta in marking strategy

    // Initializing the stratery to refine the mesh for y and w
    unsigned int yBasisRefDelay(0); // parameter for the delay in updating y_h basis for the refinement
    unsigned int wBasisRefDelay(0); // parameter for the delay in updating w-h basis for the refinement

    // Create the folder with results of the current configuration (test)
    testSpaceTime.gsCreateResultsFolder(saveToFile, vDegree, yDegree, wDegree, yBasisRefDelay, wBasisRefDelay,
                                        totalRef, adaptRefCrit, markingParamTheta);
    // -------------------------------------------------------------------------------------------------------------- //
    //! [Initialize Test Parameters]

    //! [Define Problem Data]
    // -------------------------------------------------------------------------------------------------------------- //
    gsFunctionExpr<> uDFunc, fFunc, uFunc;
    gsBoundaryConditions<> bcInfo;
    testSpaceTime.gsInitializeProblemData(uDFunc, fFunc, uFunc, bcInfo);
    testSpaceTime.gsLogProblemData();

    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Problem Data]

    //! [Define Basis]
    // -------------------------------------------------------------------------------------------------------------- //
    gsMultiBasis<> basisV, basisY, basisW;
    testSpaceTime.gsGetInitialBasis(vDegree, yDegree, wDegree,
                                   basisV, basisY, basisW,
                                   numInitUniformRefV,
                                   numInitUniformRefY,
                                   numInitUniformRefW);
    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Basis]

    //! [Define Auxiliary Structures for Storing of the Results]
    // -------------------------------------------------------------------------------------------------------------- //
    // Initialize arrays to store the error and estimators on each refinement step
    gsVector<real_t> eL2Vector(totalRef),           // || e_i ||_Q        := || u - u_i ||_Q
            eH1Vector(totalRef),                    // || grad e_i ||_Q   := || grad(u - u_i) ||_Q
            eSpaceTimeVector(totalRef),             // ||| e_i |||_h      := (|| grad_x e_i ||^2_Q + detla_h * || (u - u_i)_t ||^2_Q)^1/2
            eFullSpaceTimeVector(totalRef),         // ||| e_i |||_h      := (|| grad_x e_i ||^2_Q + detla_h * || (e_i)_t ||^2_Q + || e_i ||^2_T + detla_h * || grad_x e_i ||^2_T)^1/2
            eSpaceTimeSpaceGradVector(totalRef),    // || grad_x e_i ||_Q
            eFullSpaceTimeSpaceGradVector(totalRef),// (|| grad_x e_i ||^2_Q + || e_i ||^2_T)^2
            eSpaceTimeSolOperVector(totalRef),      // (|| Delta_x e_i ||^2_Q + || (e_i)_t||^2_Q )^1/2
            eSpaceTimeDeltaxVector(totalRef),       // || Delta_x e_i ||_Q
            eSpaceTimeDtVector(totalRef),           // || (e_i)_t||^2_Q
            eFullSpaceTimeSolOperVector(totalRef),  // || e_i ||_{L, Q}   := (|| Delta_x e_i ||^2_Q + || (e_i)_t||^2_Q + || grad_x e_i ||^2_T)^1/2
            eIdentVector(totalRef),                 // || Delta_x u_i + f - (u_i)_t ||_Q
            eFullIdentVector(totalRef),             // (|| Delta_x u_i + f - (u_i)_t ||_Q + || grad_x e_i ||_0)^1/2
            relErrorVector(totalRef-1),             // delta_i = || u_i - u_{i-1} ||, i > 0
            relError0Vector(totalRef-1),            // eps0_i  = || u_i - u_0 ||, i > 0
            thetaVector(totalRef-1),                // theta_i  = || u - u_i || / || u - u_{i-1}||, i > 0
            stopcritVector(totalRef-1),             // stop_crit_i = (1 - theta_i) * TOL * eps0_i * rho, i > 0
            majVector(totalRef),
            mdVector(totalRef),
            meqVector(totalRef),
            majhVector(totalRef),
            minVector(totalRef),
            majDeltaHVector(totalRef),
            etaVector(totalRef),
            hmaxVector(totalRef),
            hminVector(totalRef);

    // Other componenst of the majorant
    gsVector<real_t> gradxeTVector(totalRef), gradxe0Vector(totalRef),
                     eTVector(totalRef), e0Vector(totalRef), ehSigmaTVector(totalRef),
                     majIIVector(totalRef), majIIGapVector(totalRef),
                     eWL2Vector(totalRef), eWH1Vector(totalRef);

    // Initialize arrays to DOFs for v and y on each refinement step
    gsVector<index_t> vDOFs(totalRef), yDOFs(totalRef), wDOFs(totalRef);
    // Initialize vectors with assembling and computation times on each refinement step
    gsVector<double> timeAsmbV(totalRef),
            timeAsmbDivDivY(totalRef),
            timeAsmbMMY(totalRef),
            timeAsmbY(totalRef),
            timeAsmbW(totalRef),
            timeAsmbH1Error(totalRef),
            timeAsmbL2Error(totalRef),
            timeAsmbSpaceTimeError(totalRef),
            timeAsmbSpaceTimeGradSpaceError(totalRef),
            timeAsmbSpaceTimeSolOperError(totalRef),
            timeAsmbSpaceTimeDeltaxError(totalRef),
            timeAsmbSpaceTimeDtError(totalRef),
            timeAsmbSpaceTimeErrorIdentity(totalRef),
            timeAsmbMaj(totalRef),
            timeAsmbMajII(totalRef),
            timeAsmbDeltaHMajorant(totalRef),
            timeAsmbMinorant(totalRef),
            timeAsmbEtaIndicator(totalRef),
            timeAsmbEquilY(totalRef);
    timeAsmbMaj.setZero(totalRef);
    timeAsmbDeltaHMajorant.setZero(totalRef);
    timeAsmbMajII.setZero(totalRef);

    int numOfSolvers = 2;
    // matrix are used since we compare the performance of different solvers
    // [0] is the direct solver
    // [1] is the iterative solver
    gsMatrix<double> timeSolvV(totalRef, numOfSolvers),
            timeSolvY(totalRef, numOfSolvers),
            timeSolvW(totalRef, numOfSolvers),
            timeSolvEquilY(totalRef, numOfSolvers);
    timeSolvY.setZero(totalRef, 2);
    timeSolvV.setZero(totalRef, 2);
    timeSolvW.setZero(totalRef, 2);

    // Initialize auxiliary vectors of all the basis' for y to store the history along the refinement
    std::vector< gsMultiBasis<> > basisYVector;
    std::vector< gsMultiBasis<> > basisWVector;

    // Initialize auxiliary vectors of with reconstructed solutions fields and multipatches
    std::vector< gsField<> > solutionFieldVector;
    std::vector< gsMultiPatch<> > solutionMPVector;

    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Auxiliary Structures to Store the Results]
    gsPoissonPde<> heatPde(testSpaceTime.patches, bcInfo, fFunc);

    real_t theta(1);
    gsFunctionExpr<real_t> thetaFunc(std::to_string(theta), testSpaceTime.dim);

    gsSpaceTimeAssembler<real_t> spaceTimeAssemblerV;
    spaceTimeAssemblerV = gsSpaceTimeAssembler<real_t>(testSpaceTime.patches, basisV, bcInfo, *heatPde.rhs(), thetaFunc);
    spaceTimeAssemblerV.options().setInt("DirichletValues", dirichlet::l2Projection);
    spaceTimeAssemblerV.options().setInt("InterfaceStrategy", iFace::glue);

    gsSpaceTimeAssembler<real_t> spaceTimeAssemblerW;
    spaceTimeAssemblerW = gsSpaceTimeAssembler<real_t>(testSpaceTime.patches, basisW, bcInfo, *heatPde.rhs(), thetaFunc);
    spaceTimeAssemblerW.options().setInt("DirichletValues", dirichlet::l2Projection);
    spaceTimeAssemblerW.options().setInt("InterfaceStrategy", iFace::glue);

    gsMultiPatch<> mpV, mpY, mpW;
    gsMatrix<> vVector(1, 1), yVector(1, 1), wVector(1, 1);
    //gsMatrix<> vRefVector(1, 1);
    gsField<> w;
    gsField<> v;

    //! [Refinement Iterations]
    // -------------------------------------------------------------------------------------------------------------- //
    for( unsigned int refCount = 0; refCount < totalRef ; refCount++ )
    {
        // Logging the current state of the refinement of basises
        testSpaceTime.gsLogRefinementBasisInfo(refCount, NUM_PATCHES, totalRef,
                                                   spaceTimeAssemblerV.multiBasis(), basisY,
                                                   spaceTimeAssemblerW.multiBasis());
        // Get the max and min h of the mesh (must be the same for the uniform refinement)
        hmaxVector[refCount] = spaceTimeAssemblerV.multiBasis(0).basis(0).getMaxCellLength();
        hminVector[refCount] = spaceTimeAssemblerV.multiBasis(0).basis(0).getMinCellLength();

        // Update theta function for the adaitive refinement case
        if (testSpaceTime.isAdaptive) {
            theta = hminVector[refCount];
            gsFunctionExpr <real_t> updatedThetaFunc(std::to_string(theta), testSpaceTime.dim);
            spaceTimeAssemblerV.thetaUpdate(fFunc, updatedThetaFunc);
            spaceTimeAssemblerW.thetaUpdate(fFunc, updatedThetaFunc);
            thetaFunc = updatedThetaFunc;
        }
        /*
        else{ // TODO: check if this update is needed
            theta = 1;
            gsFunctionExpr <real_t> updatedThetaFunc(std::to_string(theta), testSpaceTime.dim);
            spaceTimeAssemblerV.thetaUpdate(fFunc, updatedThetaFunc);
            spaceTimeAssemblerW.thetaUpdate(fFunc, updatedThetaFunc);
            thetaFunc = updatedThetaFunc;
        }
        */
        // Reconstruct degrees of freedom, gsMultiPatch, and gsField representation of the approximation v
        testSpaceTime.gsRecontructV(refCount, spaceTimeAssemblerV, bcInfo, vVector, mpV, v, vDOFs,
                                    stopcritVector, timeAsmbV, timeSolvV);
        // Reconstruct degrees of freedom, gsMultiPatch, and gsField representation of the approximation w
        testSpaceTime.gsRecontructW(refCount, spaceTimeAssemblerW, bcInfo, wVector, mpW, w, wDOFs,
                                    stopcritVector, timeAsmbW, timeSolvW);
        // Store the last gsField- and gsMultipatch-representationa of v
        solutionFieldVector.push_back(v);
        solutionMPVector.push_back(mpV);

        //! [Error, Majorant (Optimal Flux), and Residual Estimate Computation]
        // ---------------------------------------------------------------------------------------------------------- //
        // ---------------------------------------------------------------------------------------------------------- //
        int elemNumber = spaceTimeAssemblerV.multiBasis(0).basis(0).numElements();
        std::vector<real_t> mdDistr, mIIdDistr, eH1Distr, eL2Distr, eH2Distr,
                eSpaceTimeDistr, eSpaceTimeSpaceGradDistr,
                eSpaceTimeSolOperDistr, eSpaceTimeDeltaxDistr, eSpaceTimeDtDistr,
                etaDistr, minDistr, eIdentDistr;
        mdDistr.resize(elemNumber);
        eH1Distr.resize(elemNumber);
        etaDistr.resize(elemNumber);
        minDistr.resize(elemNumber);
        eIdentDistr.resize(elemNumber);
        mIIdDistr.resize(elemNumber);

        // Deep copies of the FunctionExpression to calculate norms and residuals with OpenMP
        std::vector<gsFunctionExpr<> *> uForOMP;
        std::vector<gsFunctionExpr<> *> fForOMP;
        index_t numOfU(15);
        index_t numOfF(1);
        index_t numOfWatches(21);
        std::vector<gsCPUStopwatch> watches(numOfWatches);
        for (index_t i = 0; i < numOfU; i++)    uForOMP.emplace_back(new gsFunctionExpr<>(uFunc));
        for (index_t i = 0; i < numOfF; i++)    fForOMP.emplace_back(new gsFunctionExpr<>(fFunc));

        // add top and bottom m_sides into the collections
        std::vector<patchSide> topSides, bottomSides;
        topSides.push_back(spaceTimeAssemblerV.patches().boundaries()[d == 2 ? 0 : d*(d-1) - 1]); // spaceTimeAssemblerV.patches().boundaries()[5]:back (top of the cylinder)
        bottomSides.push_back(spaceTimeAssemblerV.patches().boundaries()[d == 2 ? 1 : d*(d-1) - 2]); // spaceTimeAssemblerV.patches().boundaries()[4]:front (bottom of the cylinder)

        #pragma omp parallel sections
        {
            #pragma omp section
            {
                gsSpaceTimeSolOperNorm<real_t> eSpaceTimeSolOper(v, * uForOMP[1]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeSolOper, eSpaceTimeSolOperDistr, elemNumber,
                                                          eSpaceTimeSolOperVector, timeAsmbSpaceTimeSolOperError, refCount);
                gsInfo << "t_{e/w} (||| e |||_L)        = " <<  timeAsmbSpaceTimeSolOperError[refCount] << " sec.\n";
            }
            #pragma omp section
            {
                gsSpaceTimeDeltaxNorm<real_t> eSpaceTimeDeltax(v, * uForOMP[2]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeDeltax, eSpaceTimeDeltaxDistr, elemNumber,
                                                      eSpaceTimeDeltaxVector, timeAsmbSpaceTimeDeltaxError, refCount);
                gsInfo << "t_{e/w} (|| Delta_x e ||)    = " <<  timeAsmbSpaceTimeDeltaxError[refCount] << " sec.\n";
            }
            #pragma omp section
            {
                gsSpaceTimeDtNorm<real_t> eSpaceTimeDt(v, * uForOMP[3]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeDt, eSpaceTimeDtDistr, elemNumber,
                                                      eSpaceTimeDtVector, timeAsmbSpaceTimeDtError, refCount);
                gsInfo << "t_{e/w} (|| Dt e ||)         = " <<  timeAsmbSpaceTimeDtError[refCount] << " sec.\n";
            }
            #pragma omp section
            {
                gsErrEstSpaceTimeIdentity<real_t> eIdentity(v, * fForOMP[0]);
                testSpaceTime.gsCalculateDistribution(eIdentity, eIdentDistr, elemNumber,
                                                          eIdentVector, timeAsmbSpaceTimeErrorIdentity, refCount);
                gsInfo << "t_{e/w} ( Id )               = " <<  timeAsmbSpaceTimeErrorIdentity[refCount] << " sec.\n";
            }
            #pragma omp section
            {
                gsSeminormH1<real_t> eH1Seminorm(v, * uForOMP[5]);
                testSpaceTime.gsCalculateDistribution(eH1Seminorm, eH1Distr, elemNumber,
                                                          eH1Vector, timeAsmbH1Error, refCount);
                gsInfo << "t_{e/w} (| e |_H1)           = " <<  timeAsmbH1Error[refCount] << " sec.\n";
            }
            #pragma omp section
            {
                gsNormL2<real_t> eL2Norm(v, * uForOMP[6]);
                testSpaceTime.gsCalculateDistribution(eL2Norm, eL2Distr, elemNumber,
                                                          eL2Vector, timeAsmbL2Error, refCount);
                gsInfo << "t_{e/w} (|| e ||_L2)         = " <<  timeAsmbL2Error[refCount] << " sec.\n";
            }
            #pragma omp section
            {
                gsSpaceTimeNorm<real_t> eSpaceTimeNorm(v, * uForOMP[7], thetaFunc);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeNorm, eSpaceTimeDistr, elemNumber,
                                                          eSpaceTimeVector, timeAsmbSpaceTimeError, refCount);
                gsInfo << "t_{e/w} (|| e ||_{s, h, Q})  = " <<  timeAsmbSpaceTimeError[refCount] << " sec.\n";
            }
            #pragma omp section
            {
                gsSpaceTimeSpaceGradNorm<real_t> eSpaceTimeSpaceGrad(v, * uForOMP[8]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeSpaceGrad, eSpaceTimeSpaceGradDistr, elemNumber,
                                                      eSpaceTimeSpaceGradVector, timeAsmbSpaceTimeGradSpaceError,
                                                      refCount);
                gsInfo << "t_{e/w} (|| grad_x e ||)     = " <<  timeAsmbSpaceTimeGradSpaceError[refCount] << " sec.\n";
            }
            #pragma omp section
            {
                watches[0].restart();
                gsSpaceTimeSpaceGradSliceNorm<real_t> eSpaceTimeSpaceGradTopSliceNorm(solutionFieldVector[refCount], * uForOMP[10], spaceTimeAssemblerV.patches().interfaces(), topSides);
                gradxeTVector[refCount] = eSpaceTimeSpaceGradTopSliceNorm.compute();
                gsInfo << "t_{e/w} (|| grad_x e ||_T)   = " << watches[0].stop() << " sec.\t || grad_x e ||_T   = "<< gradxeTVector[refCount] << "\n";
            }
            #pragma omp section
            {
                watches[1].restart();
                gsSpaceTimeSpaceGradSliceNorm<real_t> eSpaceTimeSpaceGradBottomSliceNorm(solutionFieldVector[refCount], * uForOMP[11], spaceTimeAssemblerV.patches().interfaces(), bottomSides);
                gradxe0Vector[refCount] = eSpaceTimeSpaceGradBottomSliceNorm.compute();
                gsInfo << "t_{e/w} (|| grad_x e ||_0)   = " << watches[1].stop() << " sec.\t || grad_x e ||_0   = "<< gradxe0Vector[refCount] << "\n";
            }
            #pragma omp section
            {
                watches[2].restart();
                gsSpaceTimeSliceNorm<real_t> eSpaceTimeTopSliceNorm(solutionFieldVector[refCount], * uForOMP[0], spaceTimeAssemblerV.patches().interfaces(), topSides);
                eTVector[refCount] = eSpaceTimeTopSliceNorm.compute();
                gsInfo << "t_{e/w} (|| e ||_T)          = " << watches[2].stop() << " sec.\t || e ||_T   = "<< eTVector[refCount] << "\n";
            }
            #pragma omp section
            {
                watches[3].restart();
                gsSpaceTimeSliceNorm<real_t> eSpaceTimeBottomSliceNorm(solutionFieldVector[refCount], *uForOMP[12], spaceTimeAssemblerV.patches().interfaces(), bottomSides);
                e0Vector[refCount] = eSpaceTimeBottomSliceNorm.compute();
                gsInfo << "t_{e/w} (|| e ||_0)          = " << watches[3].stop() << " sec.\t || e ||_0   = "<< e0Vector[refCount] << "\n";
            }
            #pragma omp section
            {
                watches[20].restart();
                gsSpaceTimeSigmaTNorm<real_t> eSigmaTNorm(v, * uForOMP[14], spaceTimeAssemblerV.patches().interfaces(), bottomSides, theta);
                ehSigmaTVector[refCount] = eSigmaTNorm.compute();
                gsInfo << "t_{e/w} (|| e ||_{s,h, T})   = " << watches[20].stop() << " sec.\t || e ||_{s,h, T}   = "<< ehSigmaTVector[refCount] << "\n";

            }
            #pragma omp section
            {
                watches[4].restart();
                gsSeminormH1<real_t> eWH1Seminorm(w, * uForOMP[12]);
                testSpaceTime.gsCalculateDistribution(eWH1Seminorm, eH1Distr, elemNumber, eWH1Vector, timeAsmbH1Error,
                                                      refCount);
                gsInfo << "t_{e/w} (| u - w |_H1)      = " << watches[4].stop() << " sec.\t || grad e_w ||_Q   = "<< eWH1Vector[refCount] << "\n";
            }
            #pragma omp section
            {
                watches[5].restart();
                gsNormL2<real_t> eWL2Norm(w, * uForOMP[13]);
                testSpaceTime.gsCalculateDistribution(eWL2Norm, eL2Distr, elemNumber, eWL2Vector, timeAsmbL2Error, refCount);
                gsInfo << "t_{e/w} (|| u - w ||_L2)    = " << watches[5].stop() << " sec.\t" << "|| e_w ||_Q   = " << eWL2Vector[refCount] << "\n";

            }
            #pragma omp section
            {
                // Compute the residual between two successive iterations
                if (refCount > 0) {
                    gsField<> u_cur = solutionFieldVector[refCount];
                    gsMultiPatch<> u_cur_MP = solutionMPVector[refCount];
                    gsMultiPatch<> u_prev_MP = solutionMPVector[refCount-1];

                    gsMatrix<> u_curIntPoints = u_cur.igaFunction(0).basis().anchors();

                    // truncate the fixed Dirichlet BC node and leave only free coefficients
                    gsMatrix<> u_prevIntValsFree;
                    testSpaceTime.interpolateToRefVectorWithDirichletBC(spaceTimeAssemblerV, u_prev_MP, u_curIntPoints,
                                                                        u_prevIntValsFree);
                    // reconstruct a multi-patch function from the obtained free coeffs.
                    gsMultiPatch<> u_prevIntMP;
                    spaceTimeAssemblerV.constructSolution(u_prevIntValsFree, u_prevIntMP);

                    //const gsGeometry<> &geom = u_cur_MP.piece(0);
                    //const gsMultiPatch<> mp(geom);
                    const gsMultiPatch<> mp(u_cur_MP.piece(0));

                    // reconstruct a field multi-patch function from the obtained free coeffs.
                    gsField<> u_prevInt(mp, u_prevIntMP);

                    watches[6].restart();
                    gsNormFields<real_t> ucur_minus_uprev(u_cur, u_prevInt);
                    relErrorVector[refCount - 1] = ucur_minus_uprev.compute(false, elemNumber);
                    gsInfo << "t_{e/w} (|| u_i - u_{i-1} ||_L2)    = " << watches[6].stop() << " sec.\n";
                }
            }
            #pragma omp section
            {
                // Compute the residual between the current iteration and initial u_0
                if (refCount > 0) {
                    gsField<> u_cur = solutionFieldVector[refCount];

                    gsMultiPatch<> u_0_MP = solutionMPVector[0];
                    gsMultiPatch<> u_cur_MP = solutionMPVector[refCount];

                    gsMatrix<> u_curIntPoints = u_cur.igaFunction(0).basis().anchors();

                    // truncate the fixed Dirichlet BC node and leave only free coefficients
                    gsMatrix<> u_0IntValsFree;
                    testSpaceTime.interpolateToRefVectorWithDirichletBC(spaceTimeAssemblerV, u_0_MP, u_curIntPoints,
                                                                            u_0IntValsFree);
                    // reconstruct a multi-patch function from the obtained free coeffs.
                    gsMultiPatch<> u_0IntMP;
                    spaceTimeAssemblerV.constructSolution(u_0IntValsFree, u_0IntMP);

                    //const gsGeometry<> & geom = u_cur_MP.piece(0);
                    //const gsMultiPatch<> mp(geom);
                    const gsMultiPatch<> mp(u_cur_MP.piece(0));

                    // reconstruct a field multi-patch function from the obtained free coeffs.
                    gsField<> u_0Int(mp, u_0IntMP);

                    watches[7].restart();
                    gsNormFields<real_t> ucur_minus_u0(u_cur, u_0Int);
                    relError0Vector[refCount-1] = ucur_minus_u0.compute(false, elemNumber);
                    gsInfo << "t_{e/w} (|| u_i - u_0 ||_L2)        = " << watches[7].stop() << " sec.\n";
                }
            }
        }

        // || e ||_{s, h}
        eFullSpaceTimeVector[refCount]          = math::sqrt(math::pow(eSpaceTimeVector[refCount], 2) +  math::pow(eTVector[refCount], 2) + theta * hmaxVector[refCount] * math::pow(gradxeTVector[refCount], 2));
        // || e ||_{L}
        eFullSpaceTimeSolOperVector[refCount]   = math::sqrt(math::pow(eSpaceTimeSolOperVector[refCount], 2) + math::pow(gradxeTVector[refCount], 2));
        // Id
        eFullIdentVector[refCount]              = math::sqrt(math::pow(eIdentVector[refCount], 2) + math::pow(gradxe0Vector[refCount], 2));
        // [e]
        eFullSpaceTimeSpaceGradVector[refCount] = math::sqrt(math::pow(eSpaceTimeSpaceGradVector[refCount], 2) + math::pow(eTVector[refCount], 2));

        // Logging the summary on the various errors
        testSpaceTime.gsLogRefinementIterationErrorReport(refCount, hmaxVector, hminVector,
                                                          eL2Vector, eH1Vector,
                                                          eSpaceTimeSpaceGradVector, eFullSpaceTimeSpaceGradVector,
                                                          eSpaceTimeVector, eFullSpaceTimeVector,
                                                          eSpaceTimeSolOperVector, eFullSpaceTimeSolOperVector,
                                                          eIdentVector, eFullIdentVector,
                                                          gradxe0Vector, gradxeTVector, e0Vector, eTVector, ehSigmaTVector);

        testSpaceTime.gsRecontructMajorantBasedOnOptimalFlux(refCount, basisY, yDegree,
                                                                yVector, mpY, yDOFs,
                                                                mpV, v, mpW, w,
                                                                stopcritVector,
                                                                fFunc, uFunc,
                                                                hmaxVector, hminVector,
                                                                timeAsmbDivDivY, timeAsmbMMY, timeAsmbY,
                                                                timeSolvY,
                                                                timeAsmbMaj, timeAsmbDeltaHMajorant, timeAsmbMajII,
                                                                majVector, mdVector, meqVector, majhVector, majIIVector, majIIGapVector,
                                                                mdDistr, mIIdDistr,
                                                                e0Vector,
                                                                elemNumber,
                                                                spaceTimeAssemblerV, topSides, bottomSides, theta);

        if (refCount > 0) {
            thetaVector[refCount-1]    = eH1Vector[refCount] / eH1Vector[refCount-1];
            stopcritVector[refCount-1] = (1 - thetaVector[refCount-1]) * TOL * rho * relError0Vector[refCount-1];
        }
        // Log computational costs results
        testSpaceTime.gsLogAssemblingSolvingTime(refCount, timeAsmbDivDivY, timeAsmbMMY, timeAsmbY, timeAsmbW,
                                                 timeSolvV, timeSolvY, timeSolvW,
                                                 timeAsmbMaj, timeAsmbMajII, timeAsmbSpaceTimeError);
        //! [Error and Residual Estimate Computation]

        // Log and plotToParaview the results
        testSpaceTime.gsLogRefinementIterationInfo(refCount, vDOFs, yDOFs, wDOFs,
                                                       eH1Vector, eL2Vector,
                                                       eFullSpaceTimeVector, eSpaceTimeSpaceGradVector, eFullSpaceTimeSpaceGradVector, eFullSpaceTimeSolOperVector,
                                                       relErrorVector, relError0Vector, thetaVector, stopcritVector,
                                                       majVector, majhVector, majIIVector, majIIGapVector, minVector,
                                                       etaVector, eFullIdentVector);

        if (refCount <= totalRef - 1) {
            testSpaceTime.gsSaveToFileRefinementIterationInfo(saveToFile,
                                                              v, spaceTimeAssemblerV.multiBasis(),
                                                              eSpaceTimeSpaceGradDistr, mdDistr,
                                                              eSpaceTimeSolOperDistr, eIdentDistr,
                                                              e0Vector, eTVector, gradxe0Vector, gradxeTVector,
                                                              refCount,
                                                              refCount, totalRef);
        }
        //! [Refine]
        if (refCount < totalRef - 1) {
            testSpaceTime.gsExecuteRefinement(spaceTimeAssemblerV,
                                             basisY, spaceTimeAssemblerW,
                                             basisYVector, basisWVector, mdDistr, //mdDistr, mdDistr, //mIIdDistr, // eIdentDistr, //eSpaceTimeSolOperDistr, //eH1Distr, //mdDistr, //eSpaceTimeDistr, //eSpaceTimeSolOperDistr, //eH1Distr, //eIdentDistr, //mdDistr,
                                             adaptRefCrit,
                                             markingParamTheta,
                                             refCount,
                                             yBasisRefDelay,
                                             wBasisRefDelay);
            spaceTimeAssemblerV.refresh();
            spaceTimeAssemblerW.refresh();

            // get new interpolation points for v, y, w from new reconstructed basises
            gsMatrix<> vInterpPoints = spaceTimeAssemblerV.multiBasis().basis(0).anchors();
            gsMatrix<> yInterpPoints = basisY.basis(0).anchors();
            gsMatrix<> wInterpPoints = spaceTimeAssemblerW.multiBasis().basis(0).anchors();
            gsMatrix<> vRefVector(1, 1), yRefVector(1, 1), wRefVector(1, 1);

            // evaluate new values of v, y, w based on the obtained interpolation points
            testSpaceTime.interpolateToRefVectorWithDirichletBC(spaceTimeAssemblerV, mpV, vInterpPoints, vRefVector);
            testSpaceTime.gsSetVRefVector(vRefVector);

            yRefVector = (mpY.patch(0).eval(yInterpPoints)); // returns the [1 x N] vector that needed to be transposed
            yRefVector.transposeInPlace();
            testSpaceTime.gsSetYRefVector(yRefVector);

            testSpaceTime.interpolateToRefVectorWithDirichletBC(spaceTimeAssemblerW, mpW, wInterpPoints, wRefVector);
            testSpaceTime.gsSetWRefVector(wRefVector);
        }
    }

    // -------------------------------------------------------------------------------------------------------------- //
    //! [Refinement Iterations]

    testSpaceTime.gsLogTestResults(vDegree, yDegree, wDegree,
                                  m, l,
                                  yBasisRefDelay, wBasisRefDelay,
                                  markingParamTheta, numInitUniformRefV, totalRef,
                                  vDOFs, yDOFs, wDOFs,
                                  timeAsmbV, timeAsmbDivDivY, timeAsmbMMY, timeAsmbY, timeAsmbW,
                                  timeSolvV, timeSolvY, timeSolvW,
                                  timeAsmbH1Error, timeAsmbSpaceTimeSolOperError, timeAsmbSpaceTimeDeltaxError, timeAsmbSpaceTimeDtError,
                                  timeAsmbMaj, timeAsmbDeltaHMajorant, timeAsmbMajII, timeAsmbEtaIndicator, timeAsmbSpaceTimeErrorIdentity,
                                  eL2Vector, eH1Vector,
                                  eFullSpaceTimeVector, eSpaceTimeSpaceGradVector, eFullSpaceTimeSpaceGradVector, eFullSpaceTimeSolOperVector, eSpaceTimeDeltaxVector, eSpaceTimeDtVector,
                                  relErrorVector, relError0Vector,
                                  majVector, mdVector, meqVector, majhVector, majIIVector, majIIGapVector, minVector, etaVector, eIdentVector);


    testSpaceTime.gsSaveToFileTestResults(saveToFile, vDOFs, yDOFs, wDOFs,
                                              eL2Vector, eH1Vector, eSpaceTimeSpaceGradVector, eSpaceTimeVector, eSpaceTimeSolOperVector,
                                              majVector, majhVector, majIIVector, majIIGapVector, eIdentVector,
                                              totalRef);

    gsInfo << "\nTotal execution time : " << clock_total.stop() << "\n";
    return 0;
}
