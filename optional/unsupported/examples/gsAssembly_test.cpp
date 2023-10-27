#include <iostream>

#include <gismo.h>

#include <gsSolver/gsSolverUtils.h>


using namespace gismo;


int main(int argc, char *argv[])
{   
    // Input options
    index_t numElevate  = 0;
    int numIter     = 3;
    index_t numHref     = 0;
    index_t basisDegree = 0;
    //std::string fn("planar/two_squares.xml");

    // Multipatch object
    gsMultiPatch<> mp;

    gsCmdLine cmd("Testing a multipatch problem.");
    cmd.addInt("r","hRefine", 
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("p","degree", 
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    //cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)", fn);
    cmd.addInt("e","degreeElevation", 
               "Number of degree elevation steps to perform on the Geometry's basis before solving", 
               numElevate);
        
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    //gsReadFile<>(fn, mp);
    mp = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);
    
    // Set up and elevate discretization bases
    gsMultiBasis<> bases(mp);

    if (basisDegree)
        for (size_t i = 0; i < mp.nPatches(); ++i)
            bases[i].setDegree(basisDegree);
    else if (numElevate)
        for (size_t i = 0; i < mp.nPatches(); ++i)
            bases[i].degreeElevate(numElevate);

    // Source function
    // K0 = 1; k1 = 2; k2 = 3; k3 = 4
    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",
                              "((pi*3)^2 + (pi*4)^2)*sin(pi*x*3)*sin(pi*y*4)",2);
    // Exact solution
    gsFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)+pi/10",
                              "sin(pi*x*3)*sin(pi*y*4)-pi/10",2);
    // Boundary conditions
    gsBoundaryConditions<> bcInfo;
    //Dirichlet BCs
    bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g);
    bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, &g);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &g);
    bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &g);

    // Neumann BCs
    gsFunctionExpr<> gEast ("1*pi*cos(pi*1)*sin(pi*2*y)", "3*pi*cos(pi*3)*sin(pi*4*y)",2);
    gsFunctionExpr<> gSouth("-pi*2*sin(pi*x*1)","-pi*4*sin(pi*x*3)",2);
    bcInfo.addCondition(3, boundary::east,  condition_type::neumann, &gEast);
    bcInfo.addCondition(2, boundary::east,  condition_type::neumann, &gEast);
    bcInfo.addCondition(0, boundary::south, condition_type::neumann, &gSouth);
    bcInfo.addCondition(2, boundary::south, condition_type::neumann, &gSouth);

    // List of element sizes
    std::vector <real_t> h_list;
    // List of L2-errors
    std::vector <real_t> error_list;

    bases.uniformRefine( (1<<numHref) - 1 );
    for (int i = 0; i < numIter; ++i)
    {
        //gsPoissonAssembler<real_t> pa(mp,bases,bcInfo,f,elimination,glue);
        //gsPoissonAssembler<real_t> pa(mp,bases,bcInfo,f,nitsche,glue   );
        gsPoissonAssembler<real_t> pa(mp,bases,bcInfo,f,dirichlet::nitsche,iFace::dg);
        //gsPoissonAssembler<real_t> pa(mp,bases,bcInfo,f,elimination,dg);
        pa.assemble();

        gsSparseSolver<>::CGDiagonal solver( pa.matrix() );
        gsMatrix<> solVector = solver.solve( pa.rhs() );
    
        gsInfo << "residual error: " << solver.error() << "\n";
        gsInfo << "    iterations: " << solver.iterations() << "\n";
   
        // Find the element size
        const real_t h = math::pow( (real_t) bases.size(0), -1.0 / bases.dim() );
        h_list.push_back(h);

        // Construct the solution
        gsMultiPatch<> mpsol;
        pa.constructSolution(solVector, mpsol);
        gsField<> sol(pa.patches(), mpsol);
    
        // L2 error
        error_list.push_back( sol.distanceL2(g) );

        // Get element errors: pa.estimateResidue()

        // Mark elements plus refine basis
        // at marked elements (bases....)

        gsInfo << "Iteration number "<< i <<"   L2 error: " << error_list.back() 
             << "   mesh size: " << h <<"\n";

        bases.uniformRefine();// todo: refine inplace in assembler
    }


    real_t convrateavg =  gsSolverUtils<>::convergenceRateLS(error_list,h_list);
    gsInfo << "The convergence rate is (fitted line)        : "<< convrateavg<< "\n";
    //TEST(convrateavg > convratetol);


// /*    
    // Generic matrix assembler
    gsGenericAssembler<real_t> assembler(mp, bases);
    
    assembler.assembleMass();

    index_t sz = assembler.matrix().rows();

    if ( sz < 30 )
    {
        //gsInfo<<"Mass Matrix (lower triangular): \n"<< assembler.matrix().toDense() << "\n";
        gsSparseMatrix<> fm = assembler.fullMatrix();
        gsInfo<<"Mass Matrix (full) : \n"<< fm.toDense() << "\n";

    }
    else
        gsInfo<<"Mass Matrix computed.\n";

    assembler.assembleStiffness();

    if ( sz < 30 )
    {
        //gsInfo<<"Stiffness Matrix (lower triangular): \n"<< assembler.matrix().toDense() << "\n";
        gsSparseMatrix<> fm = assembler.fullMatrix();
        gsInfo<<"Stiffness Matrix (full) : \n"<< fm.toDense() << "\n";
    }
    else
        gsInfo<<"Stiffess Matrix computed.\n";
    
//*/


    return 0;
}
