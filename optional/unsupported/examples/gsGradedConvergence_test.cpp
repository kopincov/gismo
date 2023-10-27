/** @file gsGradedConvergence_test.cpp

    @brief Testing isogeometric graded meshes

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, S. Moore
*/

#include <gismo.h>
#include <gismo_dev.h>


using namespace gismo;

template <class T>
gsBasis<T> *makeGradedBasis(short_t d, real_t u0, real_t u1,
                           unsigned interior, int degree, 
                           unsigned mult_interior = 1,
                            real_t grading=1.0,
                            bool twopointfive = false)
  {
      gsKnotVector<T> KV(degree);
      KV.initGraded(u0, u1, interior,degree, grading, mult_interior);

      switch (d)
      {
      case 1:
          //return new gsTensorBSplineBasis<1>(KV);
          return new gsBSplineBasis<T>(KV);
      case 2:
          return new gsTensorBSplineBasis<2,T>(KV,KV);
      case 3:
          if ( twopointfive )
          {
              gsKnotVector<T> KVz;
              KVz.initUniform(u0, u1, interior, degree+1, mult_interior);
              return new gsTensorBSplineBasis<3,T>(KV,KV,KVz);
          }
          else
              return new gsTensorBSplineBasis<3,T>(KV,KV,KV);
      default:
          gsWarn<<"problem.\n";
          return NULL;
      }
  }

template <class T>
gsBasis<T> *makeQuasiGradedBasis(short_t d, real_t u0, real_t u1,
                                unsigned interior, int degree, 
                                unsigned mult_interior = 1,
                                real_t grading=1.0,
                                bool twopointfive = false,
                                T border = 1.0
    )
  {
      gsKnotVector<T> KV(degree);//
      KV.insert(u0, degree+1);
      
      gsDebugVar(interior);

      unsigned interiorc = interior;
      interior   = math::floor( interior / 3.0 );
      interiorc -= interior;

      gsDebugVar(interior);
      gsDebugVar(interiorc);

      const T h2 = (border-u0) / (interiorc+1);

      gsDebugVar(h2);

      KV.insert( border, mult_interior );

      for ( unsigned i=1; i<=interiorc; i++ )
      {
          KV.insert(math::pow(i*h2, border/grading), mult_interior );
          gsInfo<<"Insert (h2): "<< math::pow(i*h2, 1.0/grading) <<"\n";
      }

      const T h1 = ( (border-u0) - math::pow( (interiorc)*h2, border/grading) ) / 2.0;
      gsDebugVar(h1);
          
      int i = 1;
      while ( border + i*h1 < u1 )
          {
              KV.insert( border + i*h1, mult_interior );
              gsInfo<<"Insert (h1): "<< border + i*h1 <<"\n";
              i++;
          }

      KV.insert(u1, degree+1);
      gsInfo<<"KV = "<< KV.detail() <<"\n";

      switch (d)
      {
      case 1:
          //return new gsTensorBSplineBasis<1>(KV);
          return new gsBSplineBasis<T>(KV);
      case 2:
          return new gsTensorBSplineBasis<2,T>(KV,KV);
      case 3:
          if ( twopointfive )
          {
              gsKnotVector<T> KVz;
              KVz.initUniform(u0, u1, interior, degree+1, mult_interior);
              return new gsTensorBSplineBasis<3,T>(KV,KV,KVz);
          }
          else
              return new gsTensorBSplineBasis<3,T>(KV,KV,KV);
      default:
          gsWarn<<"problem.\n";
          return NULL;
      }
  }

  //----------------------------------  

    real_t penalty(const gsBasis<> & B)
    {
        const int deg = B.maxDegree();
        return (deg + B.dim()) * (deg + 1) * (2.0);
    }
  
int main(int argc, char *argv[])
{   
    // Input options
    index_t numHref      = 2;
    index_t basisDegree  = 1;
    index_t basisCont    = -2;
    real_t grading   = 1.0;
    real_t scaling   = 1.0;
    bool twopointfive = false;
    bool plot        = false;

    // Multipatch object
    gsMultiPatch<> mp;

    // Pde
    memory::unique_ptr< gsPoissonPde<> > pde;

    int result = 0;
    std::string fn("planar/lshape2d_2patches_new.xml");
    std::string fn_pde("");
    
    gsCmdLine cmd("Testing a Graded Multipatch L-Shape problem.");
    cmd.addInt("r","hRefine", 
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("p","degree", 
               "Degree of the basis functions to use for solving",
               basisDegree);
    cmd.addInt("c","continuity", 
               "Continuity of the basis functions to use for solving", basisCont);
    cmd.addString("g","geometry",
                   "File containing Geometry (.xml, .axl, .txt)", fn);
    cmd.addString("q","pde","File containing a poisson PDE (.xml)", fn_pde);
    cmd.addSwitch("anisotropic", "Make anisotropic mesh", twopointfive);
    cmd.addReal("s","scale", "Scaling parameter", scaling);
    cmd.addReal("m","grade", "Grading parameter", grading);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (basisCont == -2 )
        basisCont = std::max(basisDegree-1, static_cast<index_t>(0));

    gsReadFile<>(fn, mp);
    gsDebugVar(mp);

    // Read PDE (and known exact solution)
    if ( fn_pde.empty() )
    switch ( mp.geoDim() )
    {
        case 1:
            fn_pde = "pde/poisson1d_sing.xml";
            // {
            // gsFunctionExpr<> * ff = 
            //     dynamic_cast<gsFunctionExpr<> *>(pde->solution());
            //  ff->set_x_der("((x)^v)*v/(x)-1");
            // }
            break;
        case 2:
            fn_pde = "pde/poisson2d_lshape.xml" ;
            // {
            // gsFunctionExpr<> * ff = 
            //     dynamic_cast<gsFunctionExpr<> *>(pde->solution());
            // ff->set_x_der("(2/3)*sin((2/3)*atan2(y, x)+(1/3)*pi)*x/(x^2+y^2)^(2/3)-(2/3)*(x^2+y^2)^(1/3)*cos((2/3)*atan2(y, x)+(1/3)*pi)*y/(x^2*(1+y^2/x^2))");
            // ff->set_x_der("(2/3)*sin((2/3)*atan2(y, x)+(1/3)*pi)*x/(x^2+y^2)^(2/3)-(2/3)*(x^2+y^2)^(1/3)*cos((2/3)*atan2(y, x)+(1/3)*pi)*y/(x^2*(1+y^2/x^2))");
            // }
            break;
        case 3:
            fn_pde = "pde/poisson3d_sin.xml" ;
            break;
    }
    pde = gsReadFile<>(fn_pde) ;
    gsFunctionExpr<>::uPtr solution = gsReadFile<>(fn_pde);
    
    // Scale the domain
    if ( scaling != 1.0 )
        for (size_t j = 0; j < mp.nPatches(); ++j )
            mp.patch(j).coefs().array() *= scaling;

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    // Create Dirichlet boundary conditions for all boundaries
    for (gsMultiPatch<>::const_biterator 
             bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        BCs.addCondition( *bit, condition_type::dirichlet, solution.get() );
    }

    // Run tests
    gsMatrix<> testsL2(numHref+1,7);
    testsL2.setZero();

    gsMatrix<> testsH1(numHref+1,7);
    testsH1.setZero();

    int i = 0;

    gsSparseSolver<>::CGDiagonal solver;

    do
    {
        // Set up and elevate discretization bases
        gsMultiBasis<> bases;
        for (size_t j = 0; j < mp.nPatches(); ++j )
        {
            gsBasis<> * gb =  makeGradedBasis<real_t>(mp.parDim(),0,1,(1<<i)-1,basisDegree,basisDegree-basisCont,grading, twopointfive);
                //makeQuasiGradedBasis<real_t>(mp.parDim(),0,1,(1<<i)-1,basisDegree,basisDegree-basisCont,grading, twopointfive, 0.5);
            bases.addBasis(gb);
        }
        gsInfo << "Num bases: " << bases.nBases() << "\n";
        
    // Setup the Galerkin assembers
    gsPoissonAssembler<real_t> galerkin(mp,bases,BCs,*pde->rhs(),
                                        dirichlet::elimination,iFace::glue);
    gsPoissonAssembler<real_t> galerkin_wBC(mp,bases,BCs,*pde->rhs(),
                                            dirichlet::nitsche,iFace::glue);
    gsPoissonAssembler<real_t> galerkin_dg(mp,bases,BCs,*pde->rhs(),
                                           dirichlet::nitsche,iFace::dg);
    
        gsInfo<<"Discretization Space for patch 0: \n"<< bases[0] << "\n";
        gsInfo<<"Initial DoFs             : "<< galerkin_dg.numDofs() << "\n";  
        gsInfo<<"Penalty constant         : "<< penalty(bases[0]) << "\n";  
        gsInfo<<"Gauss nodes per direction: " << "\n";  
        gsInfo << "---------------------------------------\n";
        gsInfo<<"System size (elim. BCs)  : "<< galerkin.numDofs() << "\n";  
        gsInfo<<"System size (Nitsche BCs): "<< galerkin_wBC.numDofs() << "\n";  
        gsInfo<<"System size (dg)         : "<< galerkin_dg.numDofs() << "\n";
        gsInfo << "---------------------------------------\n";
        
        gsInfo<< "Computing conforming C^0 solution..\n";
        gsStopwatch time;
        galerkin.assemble();

        solver.compute( galerkin.matrix() );
        gsMatrix<> solVector = solver.solve( galerkin.rhs() );
        const double totalTime = time.stop();
        gsInfo << "residual error: " << solver.error() << "\n";
        gsInfo << "    iterations: " << solver.iterations() << "\n";
        gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
        gsField<> sol = galerkin.constructSolution(solVector);

        gsInfo << "Computing solution with weakly imposed BCs..\n";
        time.restart();
        galerkin_wBC.assemble();
        solver.compute( galerkin_wBC.matrix() );
        solVector = solver.solve( galerkin_wBC.rhs() );
        const double totalTime_wBC = time.stop();
        gsInfo << "residual error: " << solver.error() << "\n";
        gsInfo << "    iterations: " << solver.iterations() << "\n";
        gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
        gsField<> sol_wBC = galerkin_wBC.constructSolution(solVector);
        
        gsField<> sol_dg;
        double totalTime_dg(0);
        gsMultiPatch<> mpsol;
        if ( mp.nPatches() > 1 )
        {
            gsInfo << "Computing solution with patch-wise Disc. Galerkin method..\n";
            time.restart();
            galerkin_dg.assemble();
            solver.compute( galerkin_dg.matrix() );
            solVector = solver.solve( galerkin_dg.rhs() );
            totalTime_dg = time.stop();
            gsInfo << "residual error: " << solver.error() << "\n";
            gsInfo << "    iterations: " << solver.iterations() << "\n";
            gsInfo << "     tolerance: " << solver.tolerance() << "\n"; 
            sol_dg = galerkin_dg.constructSolution(solVector);
            galerkin_dg.constructSolution(solVector, mpsol);
        }

        // Collect data
        testsL2(i,0)  =
        testsH1(i,0)= galerkin_dg.numDofs();
        if ( solution )
        {
            testsL2(i,1)= sol.distanceL2( *solution ) ;
            testsL2(i,3)= sol_wBC.distanceL2( *solution ) ;

            testsH1(i,1) = sol.distanceH1( *solution ) ;
            testsH1(i,3) = sol_wBC.distanceH1( *solution ) ;
        }
                                          
        if ( mp.nPatches() > 1 )
        {
            testsL2(i,5)= sol_dg.distanceL2( *solution ) ;
            testsH1(i,5)= sol_dg.distanceH1( *solution ) ;
        }

        if (i > 0)
        {
            testsL2(i,2)= testsL2(i-1,1) / testsL2(i,1);
            testsL2(i,4)= testsL2(i-1,3) / testsL2(i,3);
            testsL2(i,6)= testsL2(i-1,5) / testsL2(i,5);

            testsH1(i,2)= testsH1(i-1,1) / testsH1(i,1);
            testsH1(i,4)= testsH1(i-1,3) / testsH1(i,3);
            testsH1(i,6)= testsH1(i-1,5) / testsH1(i,5);
        }
        
        gsInfo << "Total time (elim. BCs): " << totalTime        << " s" << "\n";    
        gsInfo << "Total time (weak  BCs): " << totalTime_wBC << " s" << "\n";
        if ( mp.nPatches() > 1 )
            gsInfo << "Total time (full DG)  : " << totalTime_dg     << " s" << "\n";    
        gsInfo << "---------------------------------------\n";    
        gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        CONFORMING      |        WEAK BCs        |         FULL DG            \n";
        gsInfo << "    Dofs   |  L2 error  | err. ratio|  L2 error  | err. ratio|  L2 error  | err. ratio    \n" << testsL2.row(i)  << "\n";  
        gsInfo << "           |  H1 error  | err. ratio|  H1 error  | err. ratio|  DG error  | err. ratio    \n" << testsH1.row(i)  << "\n";  

        if (plot && i== numHref)
        {
            // Write approximate and exact solution to paraview files
            gsInfo<<"Plotting in Paraview...\n";
            //gsWriteParaview<>( *sol_dg, "gradedpoisson_problem", 1000, true);
            gsWriteParaview<>(sol, "gradedpoisson_problem", 1000, true);
            
            // Run paraview
            result = system("paraview gradedpoisson_problem.pvd &");
        }
                	
    } 
    while ( i++ < numHref );
    
    for(i = 1; i<= numHref; ++i)
    {   // Compute convergence rates
        testsL2(i,2)= math::log(testsL2(i,2))/std::log(2.0);
        testsL2(i,4)= math::log(testsL2(i,4))/std::log(2.0);
        testsL2(i,6)= math::log(testsL2(i,6))/std::log(2.0);

        testsH1(i,2)= math::log(testsH1(i,2))/std::log(2.0);
        testsH1(i,4)= math::log(testsH1(i,4))/std::log(2.0);
        testsH1(i,6)= math::log(testsH1(i,6))/std::log(2.0);
    }
    
    gsInfo << "Grading Parameter is " << grading << "\n";
    gsInfo << "Summary:\n\n";
    gsInfo << " (deg= "<<basisDegree <<")  |        CONFORMING      |        WEAK BCs        |         FULL DG            \n";
    gsInfo << "    Dofs   |  L2 error  | conv. rate|  L2 error  | conv. rate|  L2 error  | conv. rate    \n" << testsL2  << "\n";


    gsInfo << " (deg= "<<basisDegree <<")  |        CONFORMING      |        WEAK BCs        |         FULL DG            \n";
    gsInfo << "    Dofs   |  H1 error  | conv. rate|  H1 error  | conv. rate|  DG error  | conv. rate    \n" << testsH1  << "\n";
    
    return result;
    
}
