// [Include namespace]

#include <gismo.h>
#include <gsAssembler/gsStokesAssemblerNew.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsBramblePasciakCG.h>
#include <gsMultiGrid/gsMassSmoother.h>

using namespace gismo;
// [Include namespace]

int main(int argc, char *argv[])
{
    // [Parse command line]
    bool plot = false;
    real_t scalP=1;
    real_t alpha =1;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal("P","scalP", "scaling for the preconditioner", scalP);
    cmd.addReal("A","scalAlpha", "scaling for the method", alpha);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    // [Parse command line]

    // [Function data]
    // Define source function
    gsFunctionExpr<> f("0","0",2) ;
    //gsFunctionExpr<> f("0","0","0",3) ;

    // Boundary condition
    //2D_
    gsFunctionExpr<> U0("if(y==0, 1, 0)",2) ;
    gsFunctionExpr<> U1( "0",2) ;
    // For homogeneous term, we can use this (last argument is the dimension of the domain)
    //gsConstantFunction<> f(0.0, 0.0, 2);

    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    // [Function data]

    // [Geometry data]
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;
    patches = gsNurbsCreator<>::BSplineSquareGrid(1, 1, 0.5);
    gsInfo << "The domain is a "<< patches <<"\n";
    // [Geometry data]


    // [Boundary conditions]
    gsBoundaryConditions<> bcInfo;
    // Every patch with a boundary need to be specified. In this
    // there are in total 8 sides (two for each patch)

    // Dirichlet Boundary conditions
    // First argument is the patch number
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &U0,0 );
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &U0,0 );
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &U0,0 );
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &U0,0 );
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &U1,1 );
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &U1,1 );
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &U1,1 );
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &U1,1 );
    // [Boundary conditions]

    /*
  //Alternatively: You can automatically create Dirichlet boundary
  //conditions using one function (the exact solution) for all
  //boundaries like this:

for (gsMultiPatch<>::const_biterator
         bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
{
    bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
}
*/

    // [Refinement]
    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );
    std::vector<gsMultiBasis<> > bases(2,refine_bases);
    bases[0].degreeIncrease();

    // Number for h-refinement of the computational (trial/test) basis.
    int numRefine  = 3;

    // h-refine each basis (4, one for each patch)
    for(int k=0; k<2;k++)
    for (int i = 0; i < numRefine; ++i)
        bases[k].uniformRefine();



    // [Refinement]

    ////////////// Setup solver and solve //////////////
    // Initialize Solver
    // Setup method for handling Dirichlet boundaries, options:
    //
    // * elimination: Eliminate the Dirichlet DoFs from the linear system.
    //
    // * nitsche: Keep the Dirichlet DoFs and enforce the boundary
    //
    // condition weakly by a penalty term.
    // Setup method for handling patch interfaces, options:
    //
    // * glue:Glue patches together by merging DoFs across an interface into one.
    //   This only works for conforming interfaces
    //
    // * dg: Use discontinuous Galerkin-like coupling between adjacent patches.
    //       (This option might not be available yet)
    // [Assemble]
    gsStokesPde<real_t> stokesPDE(patches,bcInfo,&f, NULL,0.001);
    gsStokesAssemblerNew<real_t> StokesAssembler(stokesPDE,bases,dirichlet::elimination);

    // Generate system matrix and load vector
    gsInfo<< "Assembling...\n";
    StokesAssembler.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with "
           << StokesAssembler.numDofs() << " dofs.\n";
    // [Assemble]

    // [Solve]
    // Initialize the conjugate gradient solver
    gsInfo << "Solving...\n";
    gsBlockOp<>::Ptr mat = gsBlockOp<>::make(2,2);
    gsBlockOp<>::Ptr K= gsBlockOp<>::make(2,2);
    gsBlockOp<>::Ptr B= gsBlockOp<>::make(1,2);
    gsBlockOp<>::Ptr BT= gsBlockOp<>::make(2,1);
    gsBlockOp<>::Ptr MG=gsBlockOp<>::make(2,2);
    gsVector<index_t> colBlocks(2);
    colBlocks << StokesAssembler.system().colMapper(0).freeSize()+ StokesAssembler.system().colMapper(1).freeSize(),StokesAssembler.system().colMapper(2).freeSize();
    gsSparseMatrix<real_t>::BlockView view = StokesAssembler.system().blockView();//const_cast<gsSparseMatrix<real_t>* >(&StokesAssembler.matrix())->blockView(colBlocks,colBlocks);
    gsSparseMatrix<real_t>::BlockView view2=   const_cast<gsSparseMatrix<real_t>* >(&StokesAssembler.matrix())->blockView(colBlocks,colBlocks);

  //  K->addOperator(0,0,gsMatrixOp<gsSparseMatrix<real_t> >::make(view(0,0)));
  //  K->addOperator(1,1,gsMatrixOp<gsSparseMatrix<real_t> >::make(view(1,1)));

  //  BT->addOperator(0,0,gsMatrixOp<gsSparseMatrix<real_t> >::make(view(0,2)));
  //  BT->addOperator(1,0,gsMatrixOp<gsSparseMatrix<real_t> >::make(view(1,2)));

  //  B->addOperator(0,0,gsMatrixOp<gsSparseMatrix<real_t> >::make(view(2,0)));
  //  B->addOperator(0,1,gsMatrixOp<gsSparseMatrix<real_t> >::make(view(2,1)));
/*
    mat->addOperator(0,0,K);
    mat->addOperator(1,0,B);
    mat->addOperator(0,1,BT);
*/
   mat->addOperator(0,0,makeMatrixOp(view2(0,0)));
   mat->addOperator(1,0,makeMatrixOp(view2(1,0)));
   mat->addOperator(0,1,makeMatrixOp(view2(0,1)));

    gsBasis<real_t>::uPtr coarseBasis = (patches)[0].basis().clone();
    coarseBasis->degreeIncrease(1);

    std::vector<gsMultiBasis<real_t> > MG_bases;
    std::vector<gsSparseMatrix<real_t, RowMajor> > transferMatrices;

    gsGridHierarchy<>::buildByRefinement(gsMultiBasis<>(*coarseBasis), bcInfo, StokesAssembler.options(), numRefine+1)
        .moveMultiBasesTo(MG_bases)
        .moveTransferMatricesTo(transferMatrices)
        .clear();

    gsMultiGridOp<>::Ptr mg_ptr = gsMultiGridOp<>::make(view(0,0), transferMatrices);

    mg_ptr->setNumPreSmooth( 1);
    mg_ptr->setNumPostSmooth( 1 );
    mg_ptr->setNumCycles( 1 );
    real_t damping = 1;
    for(int i = mg_ptr->numLevels() == 1 ? 0 : 1; i < mg_ptr->numLevels(); ++i)
    {
        switch( 2 ) {
        case 0: mg_ptr->setSmoother(i, makeRichardsonOp(mg_ptr->matrix(i),damping)); break;
        case 1: mg_ptr->setSmoother(i, makeJacobiOp(mg_ptr->matrix(i),damping)); break;
        case 2: mg_ptr->setSmoother(i, makeGaussSeidelOp(mg_ptr->matrix(i))); break;
        case 3: mg_ptr->setSmoother(i, gsPreconditionerFromOp<>::make(mg_ptr->underlyingOp(i), makeMassSmootherOperator(MG_bases[i][0], damping, bcInfo) )); break;
        case 4: mg_ptr->setSmoother(i, gsPreconditionerFromOp<>::make(mg_ptr->underlyingOp(i), makeBoundaryCorrectedMassSmootherOperator(MG_bases[i][0], damping, bcInfo) )); break;
        case 5: mg_ptr->setSmoother(i, gsPreconditionerFromOp<>::make(mg_ptr->underlyingOp(i), makeSubspaceCorrectedMassSmootherOperator(MG_bases[i][0], damping, bcInfo) )); break;
        case 6: mg_ptr->setSmoother(i, gsPreconditionerFromOp<>::make(mg_ptr->underlyingOp(i), makeSubspaceCorrectedMassSmootherOperator(MG_bases[i][0], damping, bcInfo, true) )); break;
        case 7: mg_ptr->setSmoother(i,
                                    gsCompositePrecOp<>::make(
                                        gsPreconditionerFromOp<>::make( mg_ptr->underlyingOp(i), makeSubspaceCorrectedMassSmootherOperator(MG_bases[i][0], damping, bcInfo) ),
                                        makeGaussSeidelOp( mg_ptr->matrix(i) )
                                        )
                                    ); break;
        }
    }
    MG->addOperator(0,0,mg_ptr);
    MG->addOperator(1,1,mg_ptr);
  //  MG->addOperator(0,0,gsIdentityOp<real_t>::make(view(0,0).rows()));
  //  MG->addOperator(1,1,gsIdentityOp<real_t>::make(view(0,0).rows()));

    gsMatrix<real_t> K_mat(view2(0,0));
    gsInfo<<"K: "<<K_mat.eigenvalues().transpose()<<"\n";

    gsMatrix<real_t> MG_mat(view2(0,0));
    gsMatrix<real_t> col;
    for(int i=0; i<MG->rows();++i)
    {
        MG->apply(gsMatrix<real_t>::Identity(MG->rows(),MG->rows()).col(i),col);
        MG_mat.col(i)=col;
    }
  //  gsInfo<<"MG: \n"<<MG_mat<<"\n\n";
  //  gsInfo<<"Mat: \n"<<StokesAssembler.matrix().toDense()<<"\n\n";
    gsInfo<<"MG: "<<MG_mat.eigenvalues().transpose()<<"\n";
    gsInfo<<"K - a * MG eigenvalues: "<<(K_mat-math::min((2*alpha-1),1/(2*alpha-1))*scalP*MG_mat.inverse()).eigenvalues().transpose()<<"\n";

    gsBramblePasciakCG<gsBPCG_Types::BP,real_t> bpcg(mat, MG);
    bpcg.setCombinationScaling(alpha);
    bpcg.setPreconditionerScaling(scalP);
    gsMatrix<real_t> sol;
    sol.setZero(StokesAssembler.numDofs(),1);
    gsMatrix<real_t> history;
    bpcg.solveDetailed(StokesAssembler.rhs(),sol, history);

    gsInfo << "Solved the system with BPCG solver in .\n";
    gsInfo<< bpcg.detail()<<"\n";

    gsInfo<<"\n\n"<<history.transpose()<<"\n";
    // [Solve]
    //

    gsBlockOp<>::Ptr MG_minres=gsBlockOp<>::make(3,3);
    MG_minres->addOperator(0,0,mg_ptr);
    MG_minres->addOperator(1,1,mg_ptr);
    MG_minres->addOperator(2,2,gsIdentityOp<real_t>::make(view(2,2).rows()));
    gsMinimalResidual<> minres(mat,MG_minres);
    sol.setZero();
    minres.solve(StokesAssembler.rhs(),sol);

    gsInfo << "Solved the system with minres solver in .\n";
    gsInfo<< minres.detail()<<"\n";

    return 0;
}
