// example solving the laplace equation with the boundary element method

#include <iostream>

//#include <math.h>//

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
    index_t numRefine = 2;                       // defaults to 2
    index_t samplingPoints(1000);                // defaults to 1000
    memory::unique_ptr< gsPoissonPde<> > ppde;
    bool plot = false;  
    
    int exitCommand = 0;
    
    std::string fn("planar/lake_pd.xml");
    std::string fn_pde("pde/laplace2d_poly1.xml");

    gsCmdLine cmd("Solving Laplace problem");
    cmd.addInt("r","uniformRefine", 
               "Number of Uniform h-refinement steps to perform before solving", 
               numRefine);
    cmd.addString("g","geometry","File containing Geometry (.axl, .txt)", fn);
    cmd.addInt("s","samplingPoints", 
               "Number of sampling points to use for plotting", samplingPoints);
    cmd.addString("p","pde","File containing PDE (.xml)", fn_pde);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);     
        
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    ppde = gsReadFile<>(fn_pde);
    gsFunctionExpr<>::Ptr solution = gsReadFile<>(fn_pde);
    
    if ( !ppde )
    {
        gsWarn<< "Did not find any PDE in "<< fn_pde<<", quitting.\n";
        return 1;
    }
    
    Pdomain = gsReadFile<>( fn ) ;
    if ( !Pdomain )
    {
        gsWarn<< "Did not find any planar domain in "<< fn<<", quitting.\n";
        return 2;
    }
    
    gsInfo<<"Input : "<< *Pdomain <<"\n";
    
    gsStopwatch time;
    
    //Exact solutions on every component
    std::vector<gsFunction <> * > fis(Pdomain->numLoops());
    gsFunction<> * fi = solution.get();
    gsInfo<<"Exact solution "<< *fi <<".\n" << "\n";
    
    //******** vector of basis for every component
    std::vector<gsBasis<>* > bb;
    for(int i=0; i<Pdomain->numLoops();i++)
    {
        bb.push_back( Pdomain->loop(i).singleCurve()->basis().clone().release() );
        //dynamic_cast<gsBSplineBasis<>*>(bb.back())->knots().increaseMultiplicity();
        //dynamic_cast<gsBSplineBasis<>*>(bb.back())->degreeReduce(1);
    }

    gsInfo<<"Last Basis piece:\n"<< *bb.back()<<"\n";
    
    int i = 0;
    do
    {
        for( int r = 0; r<Pdomain->numLoops(); r++)
            fis[r] = fi->clone().release() ;

        //******** solving Laplace problem
        gsBemLaplace<> Lsolver(Pdomain.get(), fis, bb ); // takes ownership of "fis" pointer vector
        time.restart();
        gsBemSolution<> * sol_bem = Lsolver.solve() ;
                        
        gsInfo << "Total time: " << time << "\n";

        // Get a parametrization of the domain
        gsGeometry<>::uPtr sqr = gsReadFile<>("planar/lake.xml");
        gsMultiPatch<> mp( *sqr );

        gsField<> bemf( mp, *sol_bem, false);// taks ownership of sol_bem

        gsInfo << "L2 error: " << sol_bem->distanceL2( *solution ) << "\n";
        
        if ( plot && i== numRefine )
        {
            // Plotting in paraview
            gsInfo<<"Plotting in Paraview...\n";
            gsField<> exact( mp, *solution, false);
            
            //gsMatrix<> col_points, gr_points;
            //sol_bem->flux()[0]->basis().anchors_into(gr_points);
            //gsGeometry<> * outerLoop = sol_bem->getDomain()->loop(0).singleCurve();
            //outerLoop->basis().anchors_into(gr_points);
            //gsInfo << col_points <<"\n";
            //outerLoop->eval_into(gr_points,col_points);
            //gsInfo << col_points <<"\n";
            //gsMatrix<> X = col_points.row(0);
            //gsMatrix<> Y = col_points.row(1);
            //gsWriteParaviewPoints( X, Y, "laplace2d_col_points" );
            //delete outerLoop;

            gsWriteParaview<>( exact, "laplace2d_exact", samplingPoints) ; 
            gsWriteParaview<>( bemf , "laplace2d_bem"  , samplingPoints) ;
            
            //run paraview on exit
            exitCommand = system("paraview laplace2d_bem.pvd&");
        }

        // Refine for next iteration
        for(int k=0; k<Pdomain->numLoops();k++)
        {
            bb[k]->uniformRefine();
        }

    }
    while ( i++ < numRefine );
    
    freeAll(bb);

    return exitCommand;
}
