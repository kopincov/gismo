//
    /*
     * Space-time test example for the heat equation
     *
     * S. Moore
     */
#include <iostream>
#include <fstream> 

#include <gismo.h>
#include <gismo_dev.h>


#include <gsAssembler/gsSpaceTimeAssembler.h>
#include <gsAssembler/gsSpaceTimeNorm.h>
#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsNormH1Boundary.h>
#include <gsAssembler/gsNormL2Boundary.h>

using namespace gismo;


int main(int argc, char *argv[])
{   
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 1;
    index_t basisDegree = 0;
    bool plot       = false;
    index_t dim         = 2;

    // Multipatch object
    gsMultiPatch<> mp;

    // Pde
    memory::unique_ptr< gsPoissonPde<> > pde;
    gsFunctionExpr<>::Ptr solution;
    
    int result = 0;
    std::string fn("planar/square.xml");
    std::string fn_pde("");
    
    gsCmdLine cmd("Testing Space-Time problem.");
    cmd.addInt("r","hRefine","Number of refinements",  numHref);
    cmd.addInt("p","degree", "degree in all direction ", basisDegree); 
    cmd.addInt("d","dom","dimension of problem",dim);
    cmd.addString("g","geometry","File containing Geometry",fn);
    cmd.addString("q","pde","File containing a poisson PDE (.xml)", fn_pde);
    cmd.addInt("e","elevate", "degree elevation",numElevate);
    cmd.addSwitch("plot", "Plot in Paraview", plot);
    
    
    std::ofstream fout("test.txt"); 
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
        
    // Read input
    gsReadFile<>(fn, mp);
    // Read PDE (and known exact solution)  
    if ( fn_pde.empty() )
        switch ( mp.geoDim() )
        {
        case 2:
            fn_pde = "pde/time2d_1.xml"; 
            break;
        case 3:
            fn_pde = "pde/time3d_1.xml";
            break;
        case 4:
            fn_pde = "pde/time4d_1.xml";
            break;
        default:
            gsInfo <<"Got "<< mp;
        }
    pde = gsReadFile<>(fn_pde) ;
    solution  = gsReadFile<>(fn_pde, 100);
    
    // Boundary conditions
    gsBoundaryConditions<> BCs;
    for (gsMultiPatch<>::const_biterator 
          bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
         {
             BCs.addCondition( *bit, condition_type::dirichlet, *solution );
         }        
     
    // Set up and elevate discretization bases
    gsMultiBasis<> bases(mp);

    
    int i = 0;
    if (basisDegree)
            bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate,0);
    
    gsFunctionExpr<real_t>  theta("0.1", dim);

    // Run tests
    gsMatrix<> tests(numHref+1,6);
    tests.setZero();
    
    gsSparseSolver<> ::LU solver;
      
      
    //int i = 0;
    do
    {
        gsInfo<< "SetUp Problem " << "\n";
        // Setup the assemblers
        gsSpaceTimeAssembler<real_t> galerkin(mp,bases,BCs,*pde->rhs(),theta);
        galerkin.options().setInt("DirichletStrategy", dirichlet::homogeneous);
        galerkin.options().setInt("InterfaceStrategy", iFace::glue);
        
        //gsInfo<< "Starting Assembly " << "\n";
        gsStopwatch time;
        galerkin.assemble();
        const double assTime = time.stop();
        gsInfo << "Assembly time : " << assTime << " s" << "\n";
        time.restart();
        
        gsMatrix<> solVector;

        if ( galerkin.numDofs() )
        {
            solver.compute( galerkin.matrix() );
            solVector = solver.solve( galerkin.rhs() );
        }
        gsField<> sol = galerkin.constructSolution(solVector);
            
        tests(i,0)  =  bases.totalSize();
        if ( solution )
        {
            gsNormL2<real_t> L2(sol,*solution);
            tests(i,2) = L2.compute() ;
            gsSpaceTimeNorm<real_t> stnorm(sol, *solution) ;
            tests(i,4)= stnorm.compute() ;
        }

        if (i > 0)
        {
            tests(i,3)= tests(i-1,2) / tests(i,2);
            tests(i,5)= tests(i-1,4) / tests(i,4);
        }
        
        if (plot && i== numHref)
        {
            // Write approximate and exact solution to paraview files
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>(sol, "spacetime_sol", 1000);
            const gsField<> exact( galerkin.patches(), *solution, false );
            gsWriteParaview<>( exact, "spacetime_exact", 1000);
            
            // Run paraview
            result = system("paraview spacetime_sol.pvd &");
        }
        
        bases.uniformRefine();
    } 
    while ( i++ < numHref );
    
    
    for(i = 1; i<= numHref; ++i)
    {   // Compute convergence rates
        tests(i,3)= math::log(tests(i,3))/std::log(2.0);
        tests(i,5)= math::log(tests(i,5))/std::log(2.0);
    }
          
    gsInfo << "Summary:\n\n";
    gsInfo << " (deg= "<<bases[0].minDegree() <<")           TIME-DEPENDENT PROBLEM    \n";
    gsInfo << "    Dofs   |  Iter  |  L2 error | Rate  |Discrete error | Rate \n" << tests << "\n";

    return result;
}
