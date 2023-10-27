/** @file gsFluxPdeTest.cpp

    @brief Test for the correctness of the Flux PDE and is assembly

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C.Hofer
*/

#include <gismo.h>
#include <gsAssembler/gsFluxAssembler.h>
#include <gsPde/gsFluxPde.h>
#include <gsPde/gsNewtonIterator.h>
#include <gsIETI/gsNewtonIETI.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    int testcase = 0;
    real_t eps = 0.1;
    real_t p= 1.1;

    std::string m = "B";
    std::string scaling = "stiff";

    gsCmdLine cmd("Solving the Poisson problem via gsFluxPde<T>");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addString("s", "IETI.Strategy", "Choosen strategy",m);
    cmd.addString("",  "IETI.Scaling", "chosen scaling",scaling);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsOptionList option = cmd.getOptionList();
    //! [Parse command line]

    //! [Function data]
    // Define source function
   // gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",2);
    // For homogeneous term, we can use this (last argument is the dimension of the domain)
    //gsConstantFunction<> f(0.0, 0.0, 2);

    std::ostringstream strs;
    strs << eps;
    std::string eps_s = strs.str();
    strs.str("");
    strs.clear();
    strs<< p;
    std::string p_s = strs.str();

    // Define exact solution (optional)
   // gsFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)+pi/10",2);

    gsFunctionExpr<> g("sin(x+2)",2);
    gsFunctionExpr<> f("u:="+eps_s+";v:="+p_s+";w:=(u^2+cos(x+2)^2);-(v-2)*w^((v-4)/2)*cos(x+2)^2*sin(x+2)- w^((v-2)/2)*sin(x+2)",2);

  //  gsFunctionExpr<> g("exp(x)",2);
  //  gsFunctionExpr<> f("u:="+eps_s+";v:="+p_s+";w:=+(u^2+exp(2*x));(v-2)*w^((v-4)/2)*exp(3*x)+ w^((v-2)/2)*exp(x)",2);

    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n"<<std::flush;
    gsInfo<<"Exact solution "<< g <<"\n\n"<<std::flush;
    //! [Function data]

    //! [Geometry data]
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;
    // Create 4 (2 x 2) patches of squares:
    //
    // Square/patch 0 is in lower left  corner
    // Square/patch 1 is in upper left  corner
    // Square/patch 2 is in lower right corner
    // Square/patch 3 is in upper right corner
    //
    // The last argument scale the squares such that we
    // get the unit square as domain.
    switch(testcase)
    {
    case 0:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);
        break;
    case 1:
        gsFileData<> fileData("yeti_mp2.xml");
        if (fileData.getFirst(patches))
        {
        }
        break;
    }

    // patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);
    gsInfo << "The domain is a "<< patches <<"\n";

    gsBoundaryConditions<> bcInfo;

    for (gsMultiPatch<>::const_biterator
             bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
    }

    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );

    // Number for h-refinement of the computational (trial/test) basis.
    int numRefine  =1;

    // Number for p-refinement of the computational (trial/test) basis.
    int numElevate = 2;

    // h-refine each basis (4, one for each patch)
    for (int i = 0; i < numRefine; ++i)
      refine_bases.uniformRefine();

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int max_tmp = refine_bases.minCwiseDegree();

        // Elevate all degrees uniformly
        max_tmp += numElevate;
        refine_bases.setDegree(max_tmp);
        gsInfo<< "degree set to: "<<max_tmp<<"\n";
    }

    //! [Assemble]
    gsFunctionExpr<real_t> a("1",2);
    gsPiecewiseFunction<real_t> alpha;
    for(size_t np=0; np< patches.nPatches();np++)
        alpha.addPiece(a);

  //  gsPoissonFlux<real_t> PFlux(2,alpha);
    gspLaplaceFlux<real_t> pLFlux(2,eps,p);
    gsFluxPde<real_t> fluxPde(patches,bcInfo,&pLFlux,f);
//    gsFluxPde<real_t> fluxPdeP(*patches,bcInfo,&PFlux,f);
    gsFluxAssembler<real_t> fluxAssembler (fluxPde,refine_bases, dirichlet::elimination,iFace::glue);
  //  gsFluxAssembler<real_t> fluxAssemblerP (fluxPdeP,refine_bases, dirichlet::elimination,iFace::glue);

//    gsNewtonIterator<real_t> NIP(fluxAssemblerP);

//    NIP.solve();
    gsNewtonIterator<real_t> NI2(fluxAssembler);
    NI2.setTolerance(1.e-8);
    NI2.solve();

    gsInfo<<"--------------------------------------------\n\n";

  // gsNewtonIterator<real_t> NI(fluxAssembler);

    gsNewtonIETI<real_t> NI(fluxAssembler,option.getGroup("IETI"));
    NI.setTolerance(1.e-8);
    NI.setPCGTol(1.e-9);
   NI.solve();


    gsMultiPatch<real_t> solMP = NI.solution();

/*
    gsPoissonPde<real_t> ppde(*patches,bcInfo,f);
    gsPoissonAssembler<real_t>pass(ppde,refine_bases,dirichlet::elimination,iFace::glue);
    pass.assemble();

    gsInfo<<"Stiff:\n"<<pass.matrix().toDense()<<"\n";
    gsInfo<<"Stiff Linearized:\n"<<fluxAssembler.matrix().toDense()<<"\n";

    gsInfo<<"Rhs Poisson:\n"<<pass.rhs().transpose()<<"\n";
*/

    gsField<real_t> sol(patches,solMP);

    if (plot)
    {
        //! [Plot in Paraview]
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>(sol, "Flux2d", 1000);
        const gsField<> exact( fluxAssembler.patches(), g, false );
        gsWriteParaview<>( exact, "Flux2d_ex", 1000);

        // Run paraview
        return system("paraview Flux2d_ex.pvd &");
        //! [Plot in Paraview]
    }
    else
    {
        gsInfo<<"Quitting.. No output created, re-run with --plot to get a ParaView "
                "file containing Plotting image data.\n";
        return 0;
    }

}// end main

