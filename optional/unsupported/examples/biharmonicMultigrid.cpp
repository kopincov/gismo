/** @file biharmonicMultigrid.cpp

    @brief A Biharmonic example with multigrid.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

# include <gismo.h>
# include <gismo_dev.h>
# include <gsAssembler/gsBiharmonicAssembler.h>
# include <gsAssembler/gsPoissonAssembler.h>
# include <gsMultiGrid/gsGridHierarchy.h>
# include <gsSolver/gsSolverUtils.h>
# include <gsMultiGrid/gsMassSmoother.h>
# include <gsMultiGrid/gsBlockSmoother.h>


using namespace gismo;


int main(int argc, char *argv[])
{

    index_t numRefine = 3;
    index_t numDegree = 0;
    real_t damp = 1;
    gsFunctionExpr<> zeroFunc("0.0",2);
    gsCmdLine cmd("Solves a Biharmonic problem with an isogeometric discretization using a multigrid solver.");
    cmd.addInt("r", "numRefine", "Number of refinement levels", numRefine);
    cmd.addInt("p", "numDegree", "Polynomial degree", numDegree);
    cmd.addReal("", "damp", "Damp", damp);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFunctionExpr<> sol44("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)-y*y*y-x*x*x",2);
    gsFunctionExpr<> f44("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace44 ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))-6*x-6*y",2);
    gsFunctionExpr<> f_po44("16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);

    gsMultiPatch<> geo;
    //real_t a = 4;
    //geo = gsMultiPatch<>( *gsNurbsCreator<>::BSplineRectangle(0,0,1,a) );
    //geo = gsMultiPatch<>( *gsNurbsCreator<>::BSplineSquareDeg(2) );
    geo = gsMultiPatch<>( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    gsMultiBasis<> basis = gsMultiBasis<>(geo);
    gsMultiPatch<> geoPara;
    geoPara = gsMultiPatch<>( *gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(2)) );

    //gsMultiPatch<> geoPara;
    //geoPara = gsMultiPatch<>( *gsNurbsCreator<>::BSplineSquareDeg(2) );
    //gsMultiBasis<> basis = gsMultiBasis<>(geoPara);

    //p-refine to get equal polynomial degree s,t directions (for Annulus)
    basis.degreeElevate(1,0);

    //const index_t geoDim = geo.geoDim();
    for (int i = 0; i < numDegree; ++i)
    {
        basis.degreeElevate();
    }
    for (int i = 0; i < numRefine; ++i)
    {
        basis.uniformRefine();
    }

    gsBoundaryConditions<> bcInfo;
    gsBoundaryConditions<> bcInfo2;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &sol44);//Annulus: small arch lenght
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &sol44);//Annulus: Large arch lenght
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &sol44);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &sol44);
    bcInfo2.addCondition( boundary::west,  condition_type::neumann, &laplace44);//Annulus: small arch lenght
    bcInfo2.addCondition( boundary::east,  condition_type::neumann, &laplace44);//Annulus: Large arch lenght
    bcInfo2.addCondition( boundary::north, condition_type::neumann, &laplace44);
    bcInfo2.addCondition( boundary::south, condition_type::neumann, &laplace44);

    /////////////////// Setup solver ///////////////////
    //Initilize Solver
    dirichlet::strategy dirStrategy = dirichlet::elimination;
    iFace::strategy intStrategy = iFace::glue;

    gsBiharmonicAssembler<real_t> assemblerBi(geo,basis,bcInfo,bcInfo2,f44,dirStrategy,intStrategy);
    gsBiharmonicAssembler<real_t> assemblerBiPara(geoPara,basis,bcInfo,bcInfo2,f44,dirStrategy,intStrategy);

    bool doPoisson = false;
    gsPoissonAssembler<real_t>    assemblerPo(geo, basis, bcInfo, f44, dirStrategy, intStrategy);
    gsPoissonAssembler<real_t>    assemblerPoPara(geoPara, basis, bcInfo, f44, dirStrategy, intStrategy);

    bool doMass = false;
    gsFunctionExpr<real_t>  coeff_A("0.0","0","0","0.0", 2);
    gsFunctionExpr<real_t>  coeff_b("0.0","0.0", 2);
    gsFunctionExpr<real_t>  coeff_c("1.0", 2);
    gsCDRAssembler<real_t> assemblerMass(geo, basis, bcInfo,f44, coeff_A, coeff_b, coeff_c, dirStrategy, intStrategy);
    gsCDRAssembler<real_t> assemblerMassPara(geoPara, basis, bcInfo,f44, coeff_A, coeff_b, coeff_c, dirStrategy, intStrategy);


    gsInfo<<"Assembling..." << "\n";
    assemblerBi.assemble();
    assemblerBiPara.assemble();
    if (doPoisson)
    {
        assemblerPo.assemble();
        assemblerPoPara.assemble();
    }
    if (doMass)
    {
        assemblerMass.assemble();
        assemblerMassPara.assemble();
    }
    gsInfo << "basis degree: " << basis.degree(0,0)<< " and "<< basis.degree(0,1) <<"\n";
    gsInfo << "Number of DOFs: " << assemblerBi.numDofs() << " and " << assemblerBiPara.numDofs()<< "\n";
    if (doPoisson)
    {
        gsInfo << "Number of DOFs: " << assemblerPo.numDofs() << " and " << assemblerPoPara.numDofs()<< "\n";
    }
    if (doMass)
    {
        gsInfo << "Number of DOFs: " << assemblerMass.numDofs() << " and " << assemblerMassPara.numDofs()<< "\n";
    }


    gsSparseMatrix<> Bmat = assemblerBi.matrix();
    gsSparseMatrix<> BmatPara = assemblerBiPara.matrix();

    gsMatrix <> BmatDense = Bmat;
    gsMatrix <> BmatParaDense = BmatPara;
    gsMatrix <> invBmatPara = BmatParaDense.inverse();
    //gsInfo << "The condition number is inv(Bp)*B: "
    //     << gsSolverUtils<real_t>::conditionNumber(BmatDense*invBmatPara, false) << "\n";

    if (doPoisson)
    {
        gsSparseMatrix<> Kmat = assemblerPo.matrix();
        gsSparseMatrix<> KmatPara = assemblerPoPara.matrix();

        gsMatrix <> KmatDense = Kmat;
        gsMatrix <> KmatParaDense = KmatPara;
        gsMatrix <> invKmatPara = KmatParaDense.inverse();
        gsInfo << "The condition number is inv(Pp)*P: "
             << gsSolverUtils<real_t>::conditionNumber(KmatDense*invKmatPara, false) << "\n";
    }

    if (doMass)
    {
        gsSparseMatrix<> Mmat = assemblerMass.matrix();
        gsSparseMatrix<> MmatPara = assemblerMassPara.matrix();
        gsMatrix <> MmatDense = Mmat;
        gsMatrix <> MmatParaDense = MmatPara;
        //gsInfo << MmatPara << "\n";
        gsMatrix <> invMmatPara = MmatParaDense.inverse();
        //gsInfo << invMmatPara<< "\n";
        gsInfo << "The condition number is inv(Mp)*M: "
             << gsSolverUtils<real_t>::conditionNumber(MmatDense*invMmatPara, false) << "\n";
    }




    //Setting up Multigrid preconditioner
    std::vector< gsMultiBasis<real_t> > basisLvls;
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;

    gsGridHierarchy<>::buildByCoarsening(give(basis),bcInfo,gsAssembler<>::defaultOptions(),100,60)
        .moveMultiBasesTo(basisLvls)
        .moveTransferMatricesTo(transferMatrices)
        .clear();

    //Create DoF mapper for each level
    const int numTransfer = basisLvls.size() - 1;
    std::vector< gsDofMapper > dofMappers(numTransfer+1);
    const bool conforming = ( intStrategy == iFace::glue );


    for (int k = 0; k <= numTransfer; ++k)
    {

        if ( dirStrategy == dirichlet::elimination )
            basisLvls[k].getMapper(conforming, bcInfo, dofMappers[k] );
        else
            basisLvls[k].getMapper(conforming, dofMappers[k] );

    }

    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make(assemblerBi.matrix(), transferMatrices);

    std::vector<std::vector<gsVector<index_t> > > blockInfo;
    //mg.setSmoother(i, makeGaussSeidelOp(mg.matrix(i)));

    for (int i = 0; i < mg->numLevels(); ++i)
    {
        blockInfo.push_back(gsBlockOrgaizer::uniformBlocks(1,dofMappers[i].freeSize()));
        //blockInfo.push_back(gsBlockOrgaizer::boundaryBlock(basisLvls[i][0],dofMappers[i], 1));
        mg->setSmoother(i, gsGaussSeidelBlockSmoother::make(mg->matrix(i),blockInfo[i]));
    }

    /*for (int i = 0; i < mg.numLevels(); ++i)
    {
        mg.setSmoother(i,makeMassSmootherTensor<2>(*basisLvls[i], damp, bcInfo));
    }*/


    mg->setNumPreSmooth(1);
    mg->setNumPostSmooth(1);

    //End inite mg

    gsLinearOperator<>::Ptr preConMat = gsIdentityOp<>::make(assemblerBi.numDofs());
    //gsMultiGridOp<>::Ptr mg_ptr = memory::make_shared_not_owned(&mg);

    //gsConjugateGradient<> solver(assemblerBi.matrix(), preConMat);
    gsConjugateGradient<> solver(assemblerBi.matrix(), mg);

    solver.setMaxIterations(2000);
    solver.setTolerance(1e-12);

    gsMatrix<> solVector;
    //solVector.setRandom(BiharmonicAssembler.numDofs(),1);
    solVector.setZero(assemblerBi.numDofs(),1);

    //solver.solve(BiharmonicAssembler.rhs(),solVector,preConMat);
    solver.solve(assemblerBi.rhs(),solVector);


    gsInfo <<"Number of interations :" <<solver.iterations() << "\n";

    gsMultiPatch<> mpsol;
    assemblerBi.constructSolution(solVector, mpsol);
    gsField<> solF(assemblerBi.patches(), mpsol);

    real_t l2error = solF.distanceL2(sol44);
    gsInfo << "The L2 error of the solution: " << l2error << "\n";


    // Plot solution in paraview
    bool plot = false;
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsFunctionExpr<> solConst("1",2);
        gsWriteParaview<>(solF, "Biharmonic2d", 10000);
        const gsField<> exact( geo, solConst, false );
        gsWriteParaview<>( exact, "Biharmonic2d_exact", 10000);
    }

    //gsInfo << "Basis:\n" << basisLvls[0][0] << "\n";
    gsInfo << "Test is done: Cleaning up..." << "\n"; //freeAll(m_bconditions);

    return  0;
}

