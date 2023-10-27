/** @file uwbScrewMachineTest.cpp

    Author(s): H. Hornikova
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>

#include "../uwb-pilsen/uwbINSSolverSteady.h"

using namespace gismo;

int main(int argc, char *argv[])
{

        // ========================================= Settings =========================================

        bool plot = true;
        int numRefine = 0;
        int plot_pts = 50000;
        std::string filename;

        real_t viscosity = 0.1;

        //command line
        gsCmdLine cmd("Solves Navier-Stokes problem with an isogeometric discretization.");
        cmd.addSwitch("plot", "Plot result in ParaView format", plot);
        cmd.addInt("r", "uniformRefine",
                   "Number of Uniform h-refinement steps to perform before solving", numRefine);
        cmd.addInt("s", "plotSamples",
                   "Number of sample points to use for plotting", plot_pts);
        cmd.addPlainString("filename",
                           "File containing the complete problem configuration (.xml)",
                           filename);

        try { cmd.getValues(argc,argv); } catch (int& e) { return e; }

        if (numRefine<0)
        {
                gsInfo << "Number of refinements must be non-negative, quitting.\n";
                return -1;
        }

        gsAssemblerOptions opt;
        opt.dirStrategy = dirichlet::elimination;

        gsInfo << "Solving flow in the TUD geometry.\n";

        // ========================================= Define problem =========================================

        gsBoundaryConditions<> bcInfo;
        gsFunctionExpr<> *f; // external force
        gsFunctionExpr<> *Uin, *Uwall; // relative velocity BCs for inlet/rotating part/stationary part

        f = new gsFunctionExpr<>("0", "0", 2);

        // Boundary conditions
        //Uin = new gsFunctionExpr<>("0.8948", "0.4464", 2); // unit vector normal to boundary
        Uin = new gsFunctionExpr<>("0.2*0.8948", "0.2*0.4464", 2);
        Uwall = new gsFunctionExpr<>("0", "0", 2);

        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);
        bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);

        // ========================================= Define geometry =========================================

        gsFileData<> fileData(filename);

        gsMultiPatch<>* patches = new gsMultiPatch<>;
        if (fileData.has< gsTensorBSpline<2, real_t> >())
        {
                patches->addPatch(*fileData.getFirst< gsTensorBSpline<2, real_t> >());
        }
        else
        {
                gsWarn << "Input file doesn't have a gsTensorBSpline<2, real_t> inside.\n";
                return -1;
        }

        patches->computeTopology();

        gsInfo << *patches << "\n";


        // ========================================= Define basis =========================================

        gsMultiBasis<> tbasis(*patches);
        for (int i = 0; i < numRefine; ++i)
                tbasis.uniformRefine();

        std::vector< gsMultiBasis<> >  discreteBases;
        discreteBases.push_back(tbasis);//Basis for velocity
        discreteBases.push_back(tbasis);//Basis for pressure
        discreteBases[0].degreeElevate(1); //elevate the velocity space

        // ========================================= Define solver =========================================

        uwbINSPde<real_t> NSpde(*patches, bcInfo, *f, viscosity);
        uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

        uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver

        // ========================================= Solving =========================================

        gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

        gsStopwatch time;

        gsInfo << "initialization...\n";

        time.restart();
        navStokes.initialize();
        real_t Tassembly = time.stop();

        gsInfo << "Assembly time: " << Tassembly << "\n";

        time.restart();
        navStokes.solve(30, 1e-5);
        real_t Tsolve = time.stop();

        gsInfo << "Solve time: " << Tsolve << "\n";


        //// Optionally plot solution in paraview
        if (plot)
        {

                gsField<> velocity = navStokes.constructSolution(0);
                gsField<> pressure = navStokes.constructSolution(1);

                // Write solution to paraview files
                gsInfo << "Plotting in Paraview...\n";
                gsWriteParaview<>(velocity, "tud_velocity", plot_pts);
                gsWriteParaview<>(pressure, "tud_pressure", plot_pts);
        }

        delete f;
        delete Uin;
        delete Uwall;
        delete patches;

        return 0;
}
