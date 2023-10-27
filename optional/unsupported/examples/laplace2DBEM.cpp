// example solving the laplace equation with the boundary element method

#include <iostream>
#include <math.h>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsBem/gsBemLaplace.h>



using namespace gismo;

// *******************************************************
// *******************************************************
// *******************************************************


int main(int argc, char *argv[])
{
    gsPlanarDomain<>::uPtr Pdomain;          // defaults to BSplineCube
    index_t numRefine = 2;                       // defaults to 1
    index_t samplingPoints(1000);                // defaults to 1000
    memory::unique_ptr< gsPoissonPde<> > ppde;
    bool plot = false;  

    std::string fn_par("");
    std::string fn("");
    std::string fn_pde("");
    
    gsCmdLine cmd("Solving Laplace problem");
    cmd.addInt("r","uniformRefine", 
               "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addString("g","geometry","File containing Geometry (.axl, .txt)", fn);
    cmd.addString("q","fullgeometry","File containing the parametrized geometry (.axl, .txt)", fn_par);
    cmd.addInt("s","samplingPoints", 
               "Number of sampling points to use for plotting", samplingPoints);
    cmd.addString("p","pde","File containing PDE (.xml)", fn_pde);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    if ( fn_pde.empty() )
    {
      fn_pde = "pde/laplace2d_poly2.xml";
     // fn_pde = "pde/laplace2d_poly.xml";
    }
    ppde = gsReadFile<>(fn_pde);
    gsFunctionExpr<>::Ptr solution = gsReadFile<>(fn_pde);
    
    if ( !ppde )
    {
        gsWarn<< "Did not find any PDE in "<< fn<<", quitting.\n";
        return 1;
    }
    
    if ( fn.empty() )
    {
        fn     = "planar/lake_pd.xml";
        fn_par = "planar/lake.xml";
    }
    else if ( fn_par.empty() )
    {
        //fn_par = "planar/lake.xml";
    }

    Pdomain =  gsReadFile<>( fn ) ;
    if ( !Pdomain )
    {
        gsWarn<< "Did not find any planar domain in "<< fn<<", quitting.\n";
        return 2;
    }
    
    gsInfo<<"Input : "<< *Pdomain <<"\n";
    
    gsStopwatch time;
  
    //Exact solutions on every component
    std::vector<gsFunction <> * > fis;
    gsInfo<<"Exact solution "<< *solution <<".\n" << "\n";
    // allocate boundary conditions
    for( int i = 0; i<Pdomain->numLoops(); i++)
        fis.push_back(solution->clone().release());
    
//    gsFileData<> newdata;
//     newdata << Rdomain ;
//     newdata.dump("amoeba_hole_write");
    
    //******** vector of basis for every component
    std::vector<gsBasis<>* > bb;
    for(int i=0; i<Pdomain->numLoops();i++)
    {
        bb.push_back( Pdomain->loop(i).singleCurve()->basis().clone().release() );
        //dynamic_cast<gsBSplineBasis<>*>(bb.back())->degreeReduce(1);
        //dynamic_cast<gsBSplineBasis<>*>(bb.back())->knots().increaseMultiplicity();
        for (int k = 0; k < numRefine; ++k)
            bb.back()->uniformRefine();
    }
        
    //gsBemLaplace<> Lsolver( cv, fi , b );// old constructor for just one loop
    
    //******** solving Laplace problem
    gsBemLaplace<> Lsolver(Pdomain.get(), fis, bb ); // takes ownership of "fis" pointer vector
//    gsFunction<> * sol_bem = Lsolver.solve() ;
     gsBemSolution<> * sol_bem = Lsolver.solve() ;

    ///*** testing 
    gsMatrix<> u(2,2) ; // 0.3, 0.8

//    u<< 0.9, 1.5,
//        1.2,  2.0;

    u<< 1.0, 1.5, 
        .5 , 1.0;
    

//    u<<  2,3,
//         -3, 2.7;
    gsInfo<<"u\n"<< u<<"\n";
    gsInfo<<"evaluation of the bem solution \n"<< sol_bem->eval(u) <<"\n";
    gsInfo<<"evaluation of the derivatives of the bem solution \n"<<
        sol_bem->deriv(u) <<"\n";
    
    gsInfo << "Total time: " << time << "\n";

    gsInfo<<"Last Basis piece:\n"<< *bb.back()<<"\n";
    
    
    if(plot)
    {
        // Plotting in paraview
        gsInfo<<"Plotting in Paraview...\n";
        // Set a square domain instead of the boundary
        gsMultiPatch<>::uPtr  plotDomain = gsReadFile<>(fn_par);
        gsField<> exact( *plotDomain, *ppde->solution(), false);
        gsField<> bemf ( *plotDomain, *sol_bem, false); // takes ownership of sol_bem
        
        gsWriteParaview<>( exact, "laplace2d_exact", samplingPoints) ; 
        gsWriteParaview<>( bemf, "laplace2d_bem"   , samplingPoints) ;

        freeAll(bb);
        //run paraview
        char cmdParaview[100];
        strcpy(cmdParaview,"paraview laplace2d_bem.pvd&");
        return system(cmdParaview);
    }
    else
   {
       delete sol_bem;
       freeAll(bb);
       return 0;
   }

}
