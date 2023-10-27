/** @file stokesLidCavityFlow.cpp

    @brief Example solving the lid cavity flow benchmark

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>
#include <ctime>
#include <cstdlib>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsAssembler/gsStokesAssembler2.h>
#include <gsAssembler/gsStokesAssemblerNew.h>
#include <gsSolver/gsStokesIterativeSolver.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine = 0; // defaults to 1
    bool plot = false;
    index_t plot_pts = 1000;   // defaults to 1000
    
    gsCmdLine cmd("Solves Stationary Stoke's problem with an isogeometric discretization.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addInt("r","uniformRefine", 
               "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("s","plotSamples", 
               "Number of sample points to use for plotting", plot_pts);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (numRefine<0)
    { 
        gsInfo<<"Number of refinements must be non-negative, quitting.\n"; 
        return -1;
    }

    gsInfo<<"Solving the lid-driven cavity flow example, on a square.\n";

    // Source function
    //2D_
    gsFunctionExpr<> f("0","0",2) ;
    //gsFunctionExpr<> f("0","0","0",3) ;

    // Boundary condition
    //2D_
    gsFunctionExpr<> U("if(y==0, 1, 0)", "0",2) ;
    gsFunctionExpr<> U0("if(y==0, 1, 0)",2) ;
    gsFunctionExpr<> U1( "0",2) ;
    //gsFunctionExpr<> U("if(y==0, 1, 0)","0", "0",3);
    gsInfo<<"Boundary data "<< U <<".\n\n";

    // Define Geometry
    //2D_
    gsGeometry<>::uPtr geo( gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(1)) );
    //gsGeometry<>::uPtr geo = gsNurbsCreator<>::BSplineCube(1);

    // Dirichlet BCs
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &U );
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &U );
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &U );
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &U );
    //bcInfo.addCondition( boundary::front, condition_type::dirichlet, &U );
    //bcInfo.addCondition( boundary::back, condition_type::dirichlet, &U );

    // Dirichlet BCs
    gsBoundaryConditions<> bcInfoNew;
    bcInfoNew.addCondition( boundary::west,  condition_type::dirichlet, &U0,0 );
    bcInfoNew.addCondition( boundary::east,  condition_type::dirichlet, &U0,0 );
    bcInfoNew.addCondition( boundary::north, condition_type::dirichlet, &U0,0 );
    bcInfoNew.addCondition( boundary::south, condition_type::dirichlet, &U0,0 );
    bcInfoNew.addCondition( boundary::west,  condition_type::dirichlet, &U1,1 );
    bcInfoNew.addCondition( boundary::east,  condition_type::dirichlet, &U1,1 );
    bcInfoNew.addCondition( boundary::north, condition_type::dirichlet, &U1,1 );
    bcInfoNew.addCondition( boundary::south, condition_type::dirichlet, &U1,1 );
    //bcInfo.addCondition( boundary::front, condition_type::dirichlet, &U );
    //bcInfo.addCondition( boundary::back, condition_type::dirichlet, &U );

    // Define discretization space by refining the basis of the geometry
    gsBasis<>::Ptr tbasis = geo->basis().clone();
    tbasis->degreeElevate(1);
    for (int i = 0; i < numRefine; ++i)
      tbasis->uniformRefine();
    tbasis->degreeElevate(1);

    gsInfo<<"Discretization Space:" << *tbasis << "\n";

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back( gsMultiBasis<>(*tbasis) );//Basis for velocity
    discreteBases.push_back( gsMultiBasis<>(*tbasis) );//Basis for pressure
    //p-refine the velocity space (Taylor-Hood type)
    discreteBases[0].degreeElevate(1);

    gsStokesPde<real_t> stokesPDE(gsMultiPatch<>(*geo),bcInfoNew,&f, NULL,0.001);


    gsStokesAssembler2<real_t> stokes(gsMultiPatch<>(*geo), discreteBases,  bcInfo, f, dirichlet::elimination);
    stokes.setViscosity(0.001);

    gsStokesAssemblerNew<real_t> stokesNew(stokesPDE,discreteBases,dirichlet::elimination);
    stokesNew.assemble();
    gsInfo << "System size: " << stokes.numDofs() << "\n";

    // Compute solution field
    gsStopwatch time;
    stokes.assemble();
    const double assmTime = time.stop();
    gsInfo << "Assembling time: " <<  assmTime << " s" << "\n";

    gsInfo<<"Difference of matrices and rhs:\n matrix: "<<(stokes.matrix()-stokesNew.matrix()).norm()
         <<"\n rhs: "<< (stokes.rhs()-stokesNew.rhs()).norm()<<"\n";

    gsInfo<<"Stokes ddof:\n"<<stokes.dirValues()<< "\n StokesNew ddof:\n"<<stokesNew.dirValues()<<"\n";
    gsInfo<<"Stokes rhs:\n"<<stokes.rhs()<< "\n StokesNew rhs:\n"<<stokesNew.rhs()<<"\n";

    // Note: This discretization produces a singular matrix, needs to be
    // solved with an appropriate method that can handle that
    gsStokesIterativeSolver<real_t> StokesSolver(stokes.matrix(), stokes.rhs());
    time.restart();
    StokesSolver.solveWithDirectSolver();
    gsMatrix<> solVector = StokesSolver.systemSolution();
    const double solveTime = time.stop();
    gsField<> velocity  = stokes.constructSolution(solVector, 0);
    gsField<> pressure  = stokes.constructSolution(solVector, 1);
    gsInfo << "Solver time:     " <<  solveTime << " s" << "\n";
    //gsInfo <<"matrix:\n"<<  stokes.matrix().toDense() <<"\n";
    //gsInfo <<"rhs: \n"<< stokes.rhs().transpose() <<"\n";
    //gsInfo <<"solution: \n"<< solVector.transpose() <<"\n";
    
    gsInfo << "Velocity: "<<  velocity.function()  <<"\n";
    gsInfo << "Pressure: "<<  pressure.function()  <<"\n";
    //gsInfo << "Coefs: \n"<< static_cast< gsGeometry<> * >( & pressure->function() )->coefs().transpose() <<"\n";


    // Optionally plot solution in paraview
    int result = 0;
    if (plot)
    {
        // Write solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( velocity, "lidCavity_velocity", plot_pts);
        gsWriteParaview<>( pressure, "lidCavity_pressure", plot_pts);

        // Run paraview
        result = system("paraview lidCavity_velocity.pvd &");
    }

    gsStokesIterativeSolver<real_t> StokesSolverNew(stokes.matrix(), stokes.rhs());
    time.restart();
    StokesSolverNew.solveWithDirectSolver();
    solVector = StokesSolverNew.systemSolution();
    const double solveTime2 = time.stop();
    gsMultiPatch<real_t> vel;
   // velocity  = stokesNew.constructSolution(solVector, 0);
    gsVector<index_t> dims = gsVector<index_t>::LinSpaced(f.targetDim(),0,f.targetDim()-1);
    stokesNew.constructSolution(solVector,vel, dims);
    gsField<real_t> velField(stokesPDE.domain(),vel);

    pressure  = stokesNew.constructSolution(solVector, 1);
    gsInfo << "Solver time:     " <<  solveTime2 << " s" << "\n";
    //gsInfo <<"matrix:\n"<<  stokes.matrix().toDense() <<"\n";
    //gsInfo <<"rhs: \n"<< stokes.rhs().transpose() <<"\n";
    //gsInfo <<"solution: \n"<< solVector.transpose() <<"\n";

    gsInfo << "Velocity: "<<  velocity.function()  <<"\n";
    gsInfo << "Pressure: "<<  pressure.function()  <<"\n";
    //gsInfo << "Coefs: \n"<< static_cast< gsGeometry<> * >( & pressure->function() )->coefs().transpose() <<"\n";


    // Optionally plot solution in paraview
    if (plot)
    {
        // Write solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( velField, "lidCavity_velocity2", plot_pts);
        gsWriteParaview<>( pressure, "lidCavity_pressure2", plot_pts);

        // Run paraview
        result = system("paraview lidCavity_velocity2.pvd &");
    }

    return result;
}
