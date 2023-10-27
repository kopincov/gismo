/** @file KKTPervertedPoisson.cpp

    @brief A optimal controll with PDE constrain example.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/


# include <gismo.h>
# include <gsAssembler/gsKKTPervertedPoisson.h>
# include <gsAssembler/gsBiharmonicAssembler.h>
# include <gsMultiGrid/gsGridHierarchy.h>
# include <gsSolver/gsSolverUtils.h>
# include <gsSolver/gsBlockOp.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine = 3;
    index_t numDegree = 0;
    real_t alpha = 1.0;
    if (argc >= 2)
        numRefine = atoi(argv[1]);
    if (argc >= 3)
        numDegree = atoi(argv[2]);
    if (argc >= 4)
        alpha = atof(argv[3]);

    gsInfo << "Alpha is: " << alpha << "\n";

    /*
    std::map<std::string, std::string> replace;
    replace["alpha"] = "1.0";
    replace["d1"] ="-3.5";
    replace["d2"] ="3.5";*/

    //std::string replace = "u:=-3.5;v:=3.5;"; //put before function string
    //gsFunctionExpr<> solSin("0.5*x*(1.0*x - 1.0)*(1.0*x - 0.5)",2);
    //gsFunctionExpr<> fSin("3.0*(1.0*x - 0.5)",2);


    gsFunctionExpr<> zeroFunc("0.0",2);
    gsFunctionExpr<> solSin ("sin(2*pi*x) * sin(4*pi*y)",2);
    gsFunctionExpr<> fSin("-20*pi**2*sin(2*pi*x)*sin(4*pi*y)",2);

    std::vector<gsFunction<> *> desired;
    desired.push_back(new gsFunctionExpr<>("-2*pi*sin(4*pi*y)*cos(2*pi*x)",2)); //West
    desired.push_back(new gsFunctionExpr<>(" 2*pi*sin(4*pi*y)*cos(2*pi*x)",2)) ; //East
    desired.push_back(new gsFunctionExpr<>("-4*pi*sin(2*pi*x)*cos(4*pi*y)",2)); //South
    desired.push_back(new gsFunctionExpr<>(" 4*pi*sin(2*pi*x)*cos(4*pi*y)",2)); //North


    gsMultiPatch<> * geo;
    geo = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(1)) );
    //geo = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    std::vector<gsMultiBasis<> > bases_vec;
    bases_vec.push_back(gsMultiBasis<>(*geo));
    bases_vec.push_back(gsMultiBasis<>(*geo));
    bases_vec.push_back(gsMultiBasis<>(*geo));

    bases_vec[1].degreeElevate();
    //bases_vec[1].degreeElevate();
    //bases_vec[1].degreeElevate();


    bases_vec[0].degreeReduce();
    bases_vec[2].degreeReduce();

    //const index_t geoDim = geo->geoDim();
    for (int k = 0; k < 3; ++k)
    {
        for (int i = 0; i < numDegree; ++i)
        {
            bases_vec[k].degreeElevate();
        }
        for (int i = 0; i < numRefine; ++i)
        {
            bases_vec[k].uniformRefine();
        }
    }

    gsBoundaryConditions<> bcInfo;
    gsBoundaryConditions<> bcInfo2;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &solSin);//Annulus: small arch lenght
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &solSin);//Annulus: Large arch lenght
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &solSin);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &solSin);
    //bcInfo2.addCondition( boundary::west,  condition_type::neumann, &laplace44);//Annulus: small arch lenght
    //bcInfo2.addCondition( boundary::east,  condition_type::neumann, &laplace44);//Annulus: Large arch lenght
    //bcInfo2.addCondition( boundary::north, condition_type::neumann, &laplace44);
    //bcInfo2.addCondition( boundary::south, condition_type::neumann, &laplace44);

    /////////////////// Setup solver ///////////////////
    //Initilize Solver
    dirichlet::strategy dirStrategy = dirichlet::elimination;
    iFace::strategy intStrategy = iFace::glue;

    //gsBiharmonicAssembler<real_t> BiharmonicAssembler(*geo,basis,bcInfo,bcInfo2,f44,a_const,
    //                                                   dirStrategy, intStrategy);

    gsKKTPervertedPoisson<real_t> KKTPPAssembler(*geo,bases_vec,bcInfo,bcInfo2,desired,alpha,dirStrategy);

    gsBiharmonicAssembler<real_t> BiharmonicAssembler(*geo,bases_vec[1],bcInfo,bcInfo2,zeroFunc,
                                                       dirStrategy, intStrategy);

    gsInfo<<"Assembling KKT..." << "\n";
    KKTPPAssembler.assemble();
    gsInfo<<"Assembling Biharmonic (for multigrid)..." << "\n";
    BiharmonicAssembler.assemble();

    index_t totSz = KKTPPAssembler.numDofs();
    index_t fSz = KKTPPAssembler.dofMapper(0).freeSize();
    index_t uSz = KKTPPAssembler.dofMapper(1).freeSize();
    index_t wSz = KKTPPAssembler.dofMapper(2).freeSize();

    gsInfo << "Number of DOFs: " << totSz << "  Sizes: " << fSz
         << "  " << uSz << "  " << wSz << "\n" << "\n";

    //gsMatrix<> sysDense = KKTPPAssembler.matrix();
    //gsInfo << sysDense << "\n";
    // Initialize the congugate gradient solver

    real_t alphaInv = 1.0/alpha;
    gsMatrix<> mass = KKTPPAssembler.matrix().block(fSz+uSz,0,fSz,fSz);
    gsMatrix<> Amat = KKTPPAssembler.matrix().block(fSz+uSz,fSz,wSz,uSz);
    gsMatrix<> obs = KKTPPAssembler.matrix().block(fSz,fSz,uSz,uSz);

    //gsMatrix<> denseMat = KKTPPAssembler.matrix();
    //gsMatrix<> mass = denseMat.block(fSz+uSz,0,fSz,fSz);
    //gsMatrix<> Amat = denseMat.block(fSz+uSz,fSz,wSz,uSz);
    //gsMatrix<> obs = denseMat.block(fSz,fSz,uSz,uSz);

    gsInfo << "Inverting mass matrix" << "\n";
    gsMatrix<> massInv = mass.inverse();
    gsInfo << "Inverting complete!" << "\n";

    gsSparseMatrix<> massInvAlphaInv((alphaInv*massInv).sparseView());
    //gsSparseMatrix<> massInvAlphaInv(massInv.sparseView());

    gsSparseMatrix<> massInvAlpha((alpha*massInv).sparseView());

    //gsSparseMatrix<> mgMatrix((alpha*(Amat.transpose()*massInv*Amat) + obs).sparseView());
    //gsMatrix <> BiharMat = BiharmonicAssembler.matrix();
    //gsSparseMatrix<> mgMatrix((alpha*BiharMat + obs).sparseView());
    gsMatrix<> mgInverse( ((alpha*(Amat.transpose()*massInv*Amat) + obs).inverse()) );

    //gsSparseMatrix<> mgInverse( ((alpha*(Amat.transpose()*massInv*Amat) + obs).inverse()).spraseView() );


    //gsMatrix<> preconMat;
    //preconMat.setZero(totSz,totSz);
    //preconMat.block(0,0,fSz,fSz) = alphaInv*massInv;
    //preconMat.block(fSz,fSz,uSz,uSz) = (alpha*(Amat.transpose()*massInv*Amat) + obs).inverse();
    //preconMat.block(fSz+uSz,fSz+uSz,wSz,wSz) = alpha*massInv;

    //Setting up Multigrid preconditioner
    gsInfo << "Setting up multigrid" << "\n";
    std::vector< gsMultiBasis<real_t> > basisLvls;
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;
    gsGridHierarchy<>::buildByCoarsening(gsMultiBasis<>(bases_vec[1].basis(0)), bcInfo, gsAssembler<>::defaultOptions(), 100, 60)
        .moveMultiBasesTo(basisLvls)
        .moveTransferMatricesTo(transferMatrices)
        .clear();

    //Create DoF mapper for each level
    const int numTransfer = basisLvls.size() - 1;
    std::vector< gsDofMapper > dofMappers(numTransfer+1);
    const bool conforming = ( intStrategy == iFace::glue );

    for (int k = 0; k < numTransfer; ++k)
    {
        if ( dirStrategy == dirichlet::elimination )
            basisLvls[k].getMapper(conforming, bcInfo, dofMappers[k] );
        else
            basisLvls[k].getMapper(conforming, dofMappers[k] );
    }

    // Memory clean
    basisLvls.clear();

    //Add from the finest discretization
    gsDofMapper tmpDofM = KKTPPAssembler.dofMapper(1);
    dofMappers[numTransfer] = tmpDofM;

    /*
    // incorporate boundary conditions
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatricesWithBC;
    restrictToFreeDofs( transferMatrices, dofMappers, transferMatricesWithBC );
    gsMultiGridOp<> mg(mgMatrix,transferMatricesWithBC);
    for (int i = 0; i < mg.numLevels(); ++i)
        mg.setSmoother(i, makeGaussSeidelOp(mg.matrix(i)));

    mg.setNumPreSmooth(2);
    mg.setNumPostSmooth(2);*/
    //End inite mg

    gsBlockOp<>::Ptr blockPrec = gsBlockOp<>::make(3,3);
    blockPrec->addOperator(0,0,makeMatrixOp(massInvAlphaInv));
    blockPrec->addOperator(1,1,makeMatrixOp(mgInverse));
    //blockPrec->addOperator(1,1,&mg);
    blockPrec->addOperator(2,2,makeMatrixOp(massInvAlpha));

    gsMinimalResidual<> solver(KKTPPAssembler.matrix(), blockPrec);
    solver.setMaxIterations(2000);
    solver.setTolerance(1e-12);

    gsMatrix<> solVector;
    solVector.setRandom(KKTPPAssembler.numDofs(),1);
    
    gsInfo << "Starting solver, number of transfers "<< numTransfer << "\n";

    gsInfo << "\n---------------------" << "\n";
    solver.solve(KKTPPAssembler.rhs(),solVector);
    gsInfo << "Number of iterations " << solver.iterations() << "\n";
    gsInfo << "Residual Error " << solver.error() << "\n";
    gsInfo << "---------------------\n" << "\n";

    //gsInfo << "The condition number is: " << gsSolverUtils<real_t>::conditionNumber(preconMat*denseMat, false) << "\n";

    //std::exit(0);

    /*
    gsSparseSolver<>::QR solver;
    solver.analyzePattern(KKTPPAssembler.matrix() );
    solver.factorize(KKTPPAssembler.matrix());
    gsDebug << "   Eigen  lastErrorMessage: " << solver.lastErrorMessage () << "\n";
    gsMatrix<> solVector= solver.solve(KKTPPAssembler.rhs());
    */


    gsField<>::uPtr solF = memory::make_unique(KKTPPAssembler.constructSolution(solVector,0));
    gsField<>::uPtr solU = memory::make_unique(KKTPPAssembler.constructSolution(solVector,1));
    gsField<>::uPtr solW = memory::make_unique(KKTPPAssembler.constructSolution(solVector,2));

    real_t l2error_f = solF->distanceL2(zeroFunc);
    real_t l2error_u = solU->distanceL2(solSin);

    //real_t norm_f = fSin.distanceL2(zeroFunc);
    //real_t norm_u = solU->distanceL2(solSin);

    gsInfo << "norm u: " << l2error_u << "  norm f: " << l2error_f << "\n";
    //gsInfo << "The minimalization is " << l2error_u*l2error_u + alpha*l2error_f*l2error_f << "\n";
    //gsInfo << "The f norm " << norm_f << "  "<<norm_f*norm_f*alpha << "\n";


    // Plot solution in paraview
    bool plot = true;
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n" <<  "\n";
        gsWriteParaview<>(*solF, "kktControl", 1000);
        gsWriteParaview<>(*solU, "kktState", 1000);
        gsWriteParaview<>(*solW, "kktMultiplier", 1000);
        const gsField<> exactF( *geo, fSin, false );
        const gsField<> exactS( *geo, solSin, false );
        gsWriteParaview<>( exactF, "kktControl_exact", 1000);
        gsWriteParaview<>( exactS, "kktState_exact", 1000);

    }

//    gsInfo << "Basis:\n" << basis[0] << "\n";
    gsInfo << "Test is done: Cleaning up..." << "\n"; //freeAll(m_bconditions);

    delete geo;
    //delete newGeo;
    freeAll(desired);

    gsInfo << "Test is done: Exiting" << "\n";

    return  0;
}

