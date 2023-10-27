/** @file gsCoDiPacktest

    @brief Test of AD in reverse mode with use of CoDiPack with memory
    efficient solution of adjoint linear system
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Jaeschke
*/

#include <iostream>

#undef GISMO_BUILD_LIB
#include <gsCore/gsConfig.h>
#include <gsCore/gsDevConfig.h>

#ifdef GISMO_WITH_ADDSL
#include <gsAdDSL/gsAdDSL.h>
#endif

#include <gismo.h>

#include <iomanip>
#include <gsAssembler/gsCDRAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>

using namespace gismo;

#ifdef GISMO_WITH_ADDSL
// Implementation of the system solver for adDSL
void solveSystem_b(codi::BaseTape                 *tape,
                   codi::DataStore                *data,
                   codi::AdjointInterface<double> *interface) {

    gsSparseMatrix<adDSL::DSLActiveReal>* matrix;
    gsMatrix<adDSL::DSLActiveReal>* rhs;
    gsMatrix<adDSL::DSLActiveReal>* sol;

    data->getData(matrix);
    data->getData(rhs);
    data->getData(sol);

    gsSparseSolver<adDSL::DSLActiveReal>::LU solver;

    gsMatrix<adDSL::DSLActiveReal> solAdj(*sol);
    gsMatrix<adDSL::DSLActiveReal> rhsAdj(*sol);
    for(index_t i = 0; i < sol->size(); ++i) {
        solAdj[i] = adDSL::DSLActiveReal::getTape().getAdjoint((*sol)[i].gradientData());
    }

    solver.compute(matrix->transpose());
    rhsAdj = solver.solve(solAdj);

    for(index_t i = 0; i < sol->size(); ++i) {
        auto index = (*rhs)[i].gradientData();
        adDSL::DSLActiveReal::getTape().getAdjoint(index) += rhsAdj[i].getValue();
    }
    for (int e=0; e<matrix->outerSize(); ++e) {
        for (gsSparseMatrix<adDSL::DSLActiveReal>::InnerIterator it(*matrix,e); it; ++it) {
            int k = it.row();
            int l = it.col();
            adDSL::DSLActiveReal& temp1 = matrix->at(k,l);
            adDSL::DSLActiveReal::getTape().getAdjoint(temp1.gradientData())
                += -rhsAdj[l].getValue() * (*sol)[k].getValue();
        }
    }

}
#endif

// Implementation of the system solver for CoDiPack
void solveSystem_b(codi::RealReverse::TapeType    *tape,
                   codi::DataStore                *data,
                   codi::AdjointInterface<double> *interface) {

    gsSparseMatrix<codi::RealReverse>* matrix;
    gsMatrix<codi::RealReverse>* rhs;
    gsMatrix<codi::RealReverse>* sol;

    data->getData(matrix);
    data->getData(rhs);
    data->getData(sol);

    gsSparseSolver<codi::RealReverse>::LU solver;

    gsMatrix<codi::RealReverse> solAdj(*sol);
    gsMatrix<codi::RealReverse> rhsAdj(*sol);
    for(index_t i = 0; i < sol->size(); ++i) {
        solAdj[i] = (*sol)[i].getGradient();
    }

    solver.compute(matrix->transpose());
    rhsAdj = solver.solve(solAdj);

    for(index_t i = 0; i < sol->size(); ++i) {
        auto index = (*rhs)[i].getGradientData();
        tape->gradient(index) += rhsAdj[i].getValue();
    }
    for (int e=0; e<matrix->outerSize(); ++e) {
        for (gsSparseMatrix<codi::RealReverse>::InnerIterator it(*matrix,e); it; ++it) {
            int k = it.row();
            int l = it.col();
            codi::RealReverse& temp1 = matrix->at(k,l);
            tape->gradient(temp1.getGradientData()) += -rhsAdj[l].getValue() * (*sol)[k].getValue();
        }
    }
}

#ifdef GISMO_WITH_ADDSL
// Implementation of the system deleter for adDSL
void solveSystem_delete(codi::BaseTape* tape, codi::DataStore* data) {}
#endif

// Implementation of the system deleter for CoDiPack
void solveSystem_delete(codi::RealReverse::TapeType* tape, codi::DataStore* data) {}

#ifdef GISMO_WITH_ADDSL
// Implementation of the system solver for adDSL
void solveSystem(gsSparseSolver<adDSL::DSLActiveReal>::LU   &solver,
                 const gsSparseMatrix<adDSL::DSLActiveReal> &matrix,
                 const gsMatrix<adDSL::DSLActiveReal>       &rhs,
                 gsMatrix<adDSL::DSLActiveReal>             &sol) {

    codi::BaseTape& tape = codi::DSLExpression::getGlobalTape();
    tape.setPassive();

    codi::DataStore* dataHandler = new codi::DataStore();
    dataHandler->addData(&matrix);
    dataHandler->addData(&rhs);
    solver.compute( matrix );
    sol = solver.solve( rhs );

    dataHandler->addData(&sol);

    tape.pushExternalFunction(&solveSystem_b, dataHandler, &solveSystem_delete);
    tape.setActive();
    for(index_t i = 0; i < sol.size(); ++i) {
        tape.registerInput(sol[i]);
    }
}
#endif

// Implementation of the system solver for adDSL
void solveSystem(gsSparseSolver<codi::RealReverse>::LU   &solver,
                 const gsSparseMatrix<codi::RealReverse> &matrix,
                 const gsMatrix<codi::RealReverse>       &rhs,
                 gsMatrix<codi::RealReverse>             &sol) {

    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setPassive();

    codi::DataStore* dataHandler = new codi::DataStore();
    dataHandler->addData(&matrix);
    dataHandler->addData(&rhs);
    solver.compute( matrix );
    sol = solver.solve( rhs );

    dataHandler->addData(&sol);

    tape.pushExternalFunction(&solveSystem_b, dataHandler, &solveSystem_delete);
    tape.setActive();
    for(index_t i = 0; i < sol.size(); ++i) {
        tape.registerInput(sol[i]);
    }
}

int main(int argc, char* argv[])
{
    // Input options
    int numElevate  = 0;
    int numHref     = 0;
    int basisDegree = 0;
    bool EffSolv    = false;

    gsCmdLine cmd("Testing compressible Euler problem.");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("p","degree",
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    cmd.addInt("e","degreeElevation",
               "Number of degree elevation steps to perform on the Geometry's basis before solving",
               numElevate);
    cmd.addSwitch("effSolv", "Solve the system efficiently", EffSolv);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

#ifdef GISMO_WITH_ADDSL
#define RealReverse_t adDSL::DSLActiveReal
    codi::BaseTape& tape = codi::DSLExpression::getGlobalTape();
#else
#define RealReverse_t codi::RealReverse
    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
#endif

    RealReverse_t a = 3.0;
    RealReverse_t b = 2.0;

    tape.setActive();
    tape.registerInput(a);
    tape.registerInput(b);

    gsMultiPatch<RealReverse_t> mp( *gsNurbsCreator<RealReverse_t>::BSplineRectangle(0.0,0.0,a,b) );

    gsFunctionExpr<RealReverse_t>  f("0.0", 2);
    gsFunctionExpr<RealReverse_t>  g("1.0", 2);
    gsFunctionExpr<RealReverse_t>  coeff_A("1.0","0","0","1.0", 2);
    gsFunctionExpr<RealReverse_t>  coeff_b("0.2","0.4", 2);
    gsFunctionExpr<RealReverse_t>  coeff_c("0", 2);

    gsBoundaryConditions<RealReverse_t> BCs;
    for (gsMultiPatch<RealReverse_t>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        BCs.addCondition( *bit, condition_type::dirichlet, &g);
    }

    gsMultiBasis<RealReverse_t> bases(mp);
    if (basisDegree)
        bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate);
    for (int i=0; i<numHref ; i++)
        bases.uniformRefine();

    gsSparseSolver<RealReverse_t>::LU solver;
    gsCDRAssembler<RealReverse_t> galerkin(mp, bases, BCs, f, coeff_A, coeff_b, coeff_c,
                                           dirichlet::elimination, iFace::glue, stabilizerCDR::none);
    galerkin.assemble();
    gsMatrix<RealReverse_t> solVector;

    std::cout << "\n\nTape statistics before solving the system:\n\n";
    tape.printStatistics();

    if (EffSolv)
    {
        // Efficient way of solving the system
        solveSystem(solver, galerkin.matrix(), galerkin.rhs(), solVector);
    }
    else
    {
        // Inefficient way of solving the system
        solver.compute(galerkin.matrix());
        solVector = solver.solve(galerkin.rhs());
    }

    std::cout << "\n\nTape statistics after solving the system:\n\n";
    tape.printStatistics();

    gsField<RealReverse_t> sol = galerkin.constructSolution(solVector);

    gsConstantFunction<RealReverse_t> zero (0.0,2);
    RealReverse_t result = sol.distanceL2(zero,true);

    result = result * result;
#ifndef GISMO_WITH_ADDSL
    tape.registerOutput(result);
#endif
    tape.setPassive();

    std::cout << "\n\nTape statistics at the end:\n\n";
    tape.printStatistics();

#ifdef GISMO_WITH_ADDSL
    adDSL::DSLActiveReal::getTape().getAdjoint(result.index) = 1.0;
#else
    result.setGradient(1.0);
#endif

    tape.evaluate();
    std::cout << "\n\nThis code computes the derivatives of area of a 3x2 rectangle with respect to lengths of its sides:\n\n";

#ifdef GISMO_WITH_ADDSL
    std::cout << adDSL::DSLActiveReal::getTape().getAdjoint(a.index) << std::endl
              << adDSL::DSLActiveReal::getTape().getAdjoint(b.index) << std::endl;
#else
    std::cout << a.getGradient() << std::endl
              << b.getGradient() << std::endl;
#endif

    return 0;
}
