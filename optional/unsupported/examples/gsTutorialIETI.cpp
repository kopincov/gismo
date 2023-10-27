/**  gsTutorialIETI.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):     C. Hofer
    Created on:  2014-12-03

    Simple Example file how to use IETI

*/

#include <iostream>
#include <gismo.h>
#include <gismo_dev.h>
#include <gsIETI/gsIETIUtils.h>
#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIScaledDirichlet.h>

#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>

#include <gsSolver/gsConjugateGradient.h>

using namespace gismo;

//Important Note: The boundary condition has to be called with
// bcInfo.addCondition(***, *** , *** , g);
// in order that the parallelization works correctly (also with OpenMP)! So no &g !!!
// Otherwise, you have a data race!

int main (int argc, char** args)
{
    int result = 0;
    bool plot = false;
    bool nitsche = false;
    bool timings = false;
    index_t jump = 0;

    dirichlet::strategy dstrat = dirichlet::nitsche;

    //options for IETI, see also gsIETIUtils
    std::string m = "C"; //Possible options A,B,C,D -- (B==C in 2D), (A does not work good in 3D)
    std::string scaling = "coeff"; // Possible options none,mult,coeff,stiff,stiffM

    // Number for h-refinement of the computational (trail/test) basis.
    index_t numRefine  = 1;
    // Number for p-refinement of the computational (trail/test) basis.
    index_t numElevate = 0;

    gsCmdLine cmd("Hi, I will test IETIdG");
    cmd.addSwitch("n","nitsche","Use Nitsche if chosen; otherwise elimination",nitsche);
    cmd.addSwitch("p","plot","Plot solution", plot);
    cmd.addInt("r","refine","Number of refinements",numRefine);
    cmd.addInt("e","elevate","maximal elevation",numElevate);
    cmd.addInt("j","jumps","jumping coeff from 10^-j to 10^j",jump);
    cmd.addString("s", "IETI.Strategy", "Choosen strategy",m);
    cmd.addString("",  "IETI.Scaling", "chosen scaling",scaling);
    cmd.addSwitch("",  "IETI.ExtraTiming", "enable extraTimings", timings);

    try { cmd.getValues(argc,args); } catch (int rv) { return rv; }

    gsOptionList option = cmd.getOptionList();
    std::string jumpStr = util::to_string(jump);

    //Store the handling of BC
    if(nitsche)
        dstrat = dirichlet::nitsche;
    else
        dstrat = dirichlet::elimination;

    //Setup the rhs f and the dirichlet BC g
    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)",2);
    gsFunctionExpr<> g("sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)+x+y",2);

    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Dirichlet condition"<< g <<"\n" << "\n";

    //Create the Multipatch domain
    gsFileData<> fileData("yeti_mp2.xml");
    gsMultiPatch<>::uPtr patches;
    if (fileData.has< gsMultiPatch<> >())
        patches = fileData.getFirst< gsMultiPatch<> >();
    else
        return 1;

    //Create the boundary conditions:
    //Note if you use the YETI-Footprint and you select the whole boundary as dirichlet boundary (an eliminate it), then you do not
    //have any interior patch vertices. IETI still works in that case, but it does not really represent the common case.
    gsBoundaryConditions<> bcInfo;
    for(gsMultiPatch<real_t>::const_biterator it=patches->bBegin();it!=patches->bEnd();it++)
        if((*it).patch == 1)
            bcInfo.addCondition(*it, condition_type::dirichlet, g);
   gsInfo<<bcInfo<<"\n";

    //Specify the diffusion coefficient on each patch, here it is piecewise constant with value 10^j, 10^{-j}, j is an integer input.
    gsPiecewiseFunction<real_t>  alpha;
    gsFunctionExpr<> a1("1.e-"+jumpStr,2);
    gsFunctionExpr<> a2("1.e+"+jumpStr,2);
    for(size_t np = 0; np<patches->nPatches();np++)
        if(np%2 == 0)
            alpha.addPiece(a1);
        else
            alpha.addPiece(a2);
    gsInfo<<"Using jumping coeff: 10^"+jumpStr +" - 10^-"+jumpStr<<"\n";


    // Refinement:
    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( *patches );

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int max_tmp = refine_bases.maxDegree(0);
        for (int j = 1; j < patches->parDim(); ++j )
            if ( max_tmp < refine_bases.maxDegree(j) )
                max_tmp = refine_bases.maxDegree(j);

        // Elevate all degrees uniformly
        max_tmp += numElevate;

        gsInfo<<"Degree set to: "<<max_tmp<<"\n";
        refine_bases.setDegree(max_tmp);
    }

    // h-refine each basis
    for (int i = 0; i < numRefine; ++i)
        refine_bases.uniformRefine();


    //--------------------------------------
    //Setup the PDE and Assembler
    //The PDE stores the patches, BC, rhs and the other parameters of the PDE, like diffusion coefficient.

    //The PDE needs to overload two functions
    // virtual gsPde<T>* restrictToPatch(unsigned np) const
    // virtual T getCoeffForIETI(unsigned np) const
    gsPoissonHeterogeneousPde<real_t> pde(*patches, bcInfo,f,alpha);

    //The assembler must be initializeable via void initialize(..), hence, needs to implement
    // virtual void refresh();
    gsPoissonHeterogeneousAssembler<real_t> assembler(pde,refine_bases,dstrat);

    //Set the options for IETI
    //if you solve a non symmetric problem, then you have to specify the options
    // "NonSymmetric" as true, e.g. defaultOpt.setSwitch("NonSymmetric",true);
    gsOptionList defaultOpt= gsIETIAssembler<real_t>::defaultOptions();
    defaultOpt.update(option.getGroup("IETI"),gsOptionList::addIfUnknown);

    //Setup the IETI solver
    gsStopwatch time;
    gsIETIAssembler<real_t> ass(assembler); //IETI needs as input just an fully working assembler
    ass.setOptions(defaultOpt);
    ass.init();
    gsInfo << "Time for IETI Initialization: "<< time<<"\n"<<std::flush;

    //Assemble the local Matrices and the local LU factorizations
    time.restart();
    ass.assemble();
    time.stop();
    gsInfo << "Time for assembling: "<< time<<"\n"<<std::flush;

    //IETI solves for the lagrange multipliers, enforcing continuity across the interface

    //Setup the linearoperator for the iterative solver
    gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(ass);
    solv->init();
    //Setup the preconditioner
    gsScaledDirichletPrecond<real_t>::Ptr precDir = gsScaledDirichletPrecond<real_t>::make(ass);
    //Setup the CG
    gsConjugateGradient<> PCG(solv,precDir);
    PCG.setMaxIterations(100);
    PCG.setCalcEigenvalues(true);

    //Intial guess
    gsMatrix<> solVector(ass.systemSize(),ass.numberRhs());
    solVector.setZero();

    //Solve the IETI system, we obtain the lagrange multipliers
    time.restart();
    PCG.solve(solv->getRhs(),solVector);
    time.stop();
    gsInfo << "Time for solving: "<< time<<"\n"<<std::flush;

    //Information about the eigenvalues of the preconditioned system
    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo<<"eigenvalues: "<<"\n"<<eigs.minCoeff()<<" -- "<<eigs.maxCoeff()<<"\n";
    gsInfo<<"Number of iterations: "<<PCG.iterations()<<"\n";

    //Backsubstituting the lagrange multipliers in order to obtain the solution from Ku=f
    gsMatrix<> solution;
    time.restart();
    solv->calculateSolution(solVector,solution);
    time.stop();
    gsInfo << "Time for post processing: "<< time<<"\n"<<std::flush;

    //-----------------------------------------------------------------
    //Comparism with original system
    assembler.assemble();
    gsSparseMatrix<> K = assembler.fullMatrix();
    gsSparseSolver<>::LU solver;
    solver.compute(K);
    gsMatrix<> solutionDirectsolver = solver.solve( assembler.rhs() );

    //Error between u_IETI and u_Direct
    gsInfo<<"Error ||u_IETI - u_Direct||= "<<(solution - solutionDirectsolver).norm()<<"\n\n";
    //--------------------------------------------------------------------------------------
    //If specified, plot the domain
    if (plot)
    {
        gsField<> solIETI = assembler.constructSolution(solution);
        gsWriteParaview<>( solIETI, "solutionIETI", 1000);

        // Run paraview
        result = system("paraview IETI.pvd  &");
    }

    //return  1 for failures and 0 for success
    return result;
}
