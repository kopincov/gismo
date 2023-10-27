/** @file gsParameterDomainPreconditionerTest.cpp

    @brief Checks if the problem on the parameter domain is a preconditioner for the problem on the pysical domain for biharmonic problems.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, S. Takacs
*/

# include <gismo.h>
# include <gismo_dev.h>
# include <gsPde/gsBiharmonicPde.h>

# include <gsRecipeAssembler/gsRecipeAssemblerBiharmonicSogn.h>
# include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>
# include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>
# include <gsAssembler/gsBiharmonicAssembler.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <gsSolver/gsKronecker.h>

using namespace gismo;

gsTensorBSpline<2,real_t>::uPtr BSplineParallelogram(int deg)
{
    gsKnotVector<real_t> KV (0,1,0,2) ;
    gsMatrix<real_t> C(4,2) ;

    C <<  0 , 0 , 1 , 0
        , 0.5 , 1 , 1.5 , 1 ;

    gsTensorBSpline<2,real_t> * result = new gsTensorBSpline<2,real_t>(KV,KV,give(C));
    result->degreeElevate(deg-1);
    return memory::make_unique(result);
}
gsTensorBSpline<2,real_t>::uPtr BSplineTrapezoidal(int deg)
{
    gsKnotVector<real_t> KV (0,1,0,2) ;
    gsMatrix<real_t> C(4,2) ;

    C <<  0 , 0 , 1 , 0
        , 0 , 1 , 0.5 , 1 ;

    gsTensorBSpline<2,real_t> * result = new gsTensorBSpline<2,real_t>(KV,KV,give(C));
    result->degreeElevate(deg-1);
    return memory::make_unique(result);
}

int main(int argc, char *argv[])
{
    index_t numRefine = 1;
    index_t numDegree = 0;
    index_t geoType = 1;
    index_t assemblerType = 1;

    gsCmdLine cmd("Checks if the problem on the parameter domain is a preconditioner for the problem on the pysical domain.");
    cmd.addInt("r", "refine", "Number of refinements", numRefine);
    cmd.addInt("p", "degree", "Polynomial degree", numDegree);
    cmd.addInt("g", "geotype", "1=BSplineFatQuarterAnnulus, 2=BSplineTrapezoidal, 3=BSplineRectangle, 4=BSplineParallelogram", geoType);
    cmd.addInt("a", "assemblertype", "1=gsBiharmonicAssembler, 2=gsRecipeAssemblerBiharmonicSogn, 3=gsPoissonAssembler", assemblerType);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFunctionExpr<> sol44("0",2);
    gsFunctionExpr<> f44("0",2);
    gsFunctionExpr<> laplace44 ("0",2);

    gsGeometry<>::uPtr geoParam = memory::convert_ptr<gsGeometry<> >(gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(2)));
    gsGeometry<>::uPtr geoPhys;
    if ( geoType == 1 )
    {
        geoPhys = gsNurbsCreator<>::BSplineFatQuarterAnnulus();
        geoPhys->degreeElevate(1,0);
    }
    else if ( geoType == 2 )
    {
        geoPhys = BSplineTrapezoidal(2);
    }
    else if ( geoType == 3 )
    {
        geoPhys = gsNurbsCreator<>::BSplineRectangleWithPara(0,0,1,4);
        geoPhys->degreeElevate(1);
    }
    else if ( geoType == 4 )
    {
        geoPhys = BSplineParallelogram(2);
    }
    else
    {
        gsInfo << "Unknown geotype.\n\n";
        return -1;
    }
    
    gsInfo << "Geometry for the parameter domain:\n" << *geoParam << "\n\n";
    gsInfo << "Geometry for the physical domain:\n" << *geoPhys << "\n\n";
    
    gsMultiPatch<> mpParam(*geoParam);
    gsMultiPatch<> mpPhys(*geoPhys);
        
    
    gsMultiBasis<> basisParam = gsMultiBasis<>(mpParam);
    gsMultiBasis<> basisPhys = gsMultiBasis<>(mpPhys);
    basisParam[0].component(0).uniformRefine();
    basisPhys[0].component(0).uniformRefine();

    //const index_t geoDim = geo->geoDim();
    for (int i = 0; i < numRefine; ++i)
    {
        basisParam.uniformRefine();
        basisPhys.uniformRefine();
    }
    for (int i = 2; i < numDegree; ++i)
    {
        basisParam.degreeElevate();
        basisPhys.degreeElevate();
    }

    GISMO_ENSURE( basisParam.minCwiseDegree() == basisParam.maxCwiseDegree(), "Problems with the degrees");
    GISMO_ENSURE( basisPhys.minCwiseDegree()  == basisPhys.maxCwiseDegree(),  "Problems with the degrees");

    gsInfo << "Basis on the parameter domain: degree = " << basisParam.minCwiseDegree() << ", size = " << basisParam.totalSize() << "\n";
    gsInfo << "Basis on the physical  domain: degree = " << basisPhys.minCwiseDegree()  << ", size = " << basisPhys.totalSize()  << "\n";
    gsInfo << "\n\n";

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

    gsSparseMatrix<> matParam;
    gsSparseMatrix<> matPhys;
    
    if ( assemblerType == 1 )
    {
        dirichlet::strategy dirStrategy = dirichlet::elimination;
        iFace::strategy intStrategy = iFace::glue;

        gsBiharmonicAssembler<real_t> BiharmonicAssemblerParam(mpParam,basisParam,bcInfo,bcInfo2,f44,dirStrategy,intStrategy);

        gsBiharmonicAssembler<real_t> BiharmonicAssemblerPhys(mpPhys,basisPhys,bcInfo,bcInfo2,f44,dirStrategy,intStrategy);

        
        gsInfo << "Assembling..." << std::flush;
        BiharmonicAssemblerParam.assemble();
        BiharmonicAssemblerPhys.assemble();
        matParam = BiharmonicAssemblerParam.matrix();
        matPhys = BiharmonicAssemblerPhys.matrix();
        
        gsInfo << "done.\n";
    }
    else if ( assemblerType == 2 )
    {
        gsBiharmonicPde<real_t> pdeBiharmonicParam(mpParam,bcInfo,gsBoundaryConditions<>(),f44);
        gsRecipeAssemblerBiharmonicSogn BiharmonicAssemblerParam(pdeBiharmonicParam);
        BiharmonicAssemblerParam.setDirichletStrategy(dirichlet::elimination);

        std::vector<gsPhysicalSpace*> phySpaceParam;
        phySpaceParam.push_back(new gsPhysicalSpaceScalar(basisParam,mpParam,INVERSE_COMPOSITION));
        BiharmonicAssemblerParam.setSpace(phySpaceParam);


        gsBiharmonicPde<real_t> pdeBiharmonicPhys(mpPhys,bcInfo,gsBoundaryConditions<>(),f44);
        gsRecipeAssemblerBiharmonicSogn BiharmonicAssemblerPhys(pdeBiharmonicPhys);
        BiharmonicAssemblerPhys.setDirichletStrategy(dirichlet::elimination);

        std::vector<gsPhysicalSpace*> phySpacePhys;
        phySpacePhys.push_back(new gsPhysicalSpaceScalar(basisPhys,mpPhys,INVERSE_COMPOSITION));
        BiharmonicAssemblerPhys.setSpace(phySpacePhys);

        gsPoissonPde<real_t> pdePoissonParam(mpParam, bcInfo, f44);
        gsRecipeAssemblerPoisson PoissonAssemblerParam(pdePoissonParam);
        PoissonAssemblerParam.setDirichletStrategy(dirichlet::elimination);
        std::vector<gsPhysicalSpace*> phySpacePoisson;
        phySpacePoisson.push_back(new gsPhysicalSpaceScalar(basisParam,mpParam,INVERSE_COMPOSITION));
        PoissonAssemblerParam.setSpace(phySpacePoisson);
        PoissonAssemblerParam.assemble();



        gsInfo << "Assembling..." << std::flush;
        BiharmonicAssemblerParam.assemble();
        BiharmonicAssemblerPhys.assemble();
        matParam = BiharmonicAssemblerParam.getSystemMatrix();
        matPhys = BiharmonicAssemblerPhys.getSystemMatrix();
        gsInfo << "done.\n";
    
    }
    else if ( assemblerType == 3 )
    {
        dirichlet::strategy dirStrategy = dirichlet::elimination;
        iFace::strategy intStrategy = iFace::glue;

        gsPoissonAssembler<real_t> PoissonAssemblerParam(mpParam,basisParam,bcInfo,f44,dirStrategy,intStrategy);

        gsPoissonAssembler<real_t> PoissonAssemblerPhys(mpPhys,basisPhys,bcInfo,f44,dirStrategy,intStrategy);

        
        gsInfo << "Assembling..." << std::flush;
        PoissonAssemblerParam.assemble();
        PoissonAssemblerPhys.assemble();
        matParam = PoissonAssemblerParam.matrix();
        matPhys = PoissonAssemblerParam.matrix();
        gsInfo << "done.\n";

        gsSparseMatrix<real_t> M1,M2,M3;
        assembleParameterMassForTensorProductSpace(basisParam[0],bcInfo,M1);
        assembleGeneralizedParameterStiffnessForTensorProductSpace(basisParam[0],bcInfo,real_t(0.),real_t(1.),M2);
        gsGenericAssembler<real_t> gen(mpParam, basisParam, gsAssembler<real_t>::defaultOptions(),&bcInfo);
        M3= gen.assembleMass();

        gsKroneckerOp<real_t>::Ptr ptr;
        assembleParameterMassInverseForTensorProductSpace(basisParam[0],bcInfo,ptr);
        gsMatrix<real_t> Minv1;
        ptr->toMatrix(Minv1);
        gsMatrix<real_t> Minv2 = M3.toDense().inverse();

        gsInfo<<"Mat1\n"<<M1.toDense()<<"\n\n\nMat2\n"<<M2.toDense()<<"\n\n\nMat3\n"<<M3.toDense()<<"\n\n\nMinv1\n"<<Minv1<<"\n\n\nMinv2\n"<<Minv2<<"\n\n";

    }
    else
    {
        gsInfo << "Unknown assemblertype.\n\n";
        return -1;
    }
    
    if( numRefine < 4 )
    {
        gsInfo << "\nMatrix on the parameter domain:\n" << gsMatrix<>(matParam) << "\n\n\n";
        gsInfo <<   "Matrix on the physical  domain:\n" << gsMatrix<>(matPhys) << "\n\n\n";
    }

    gsInfo << "Number of DOFs on the parameter domain: " << matParam.rows() << "\n";
    gsInfo << "Number of DOFs on the physical  domain: " << matPhys.rows()  << "\n\n";


    // Initialize the congugate gradient solver
    gsInfo << "Setup the direct solver (sparse Cholesky) for the parameter domain..." << std::flush;

    gsLinearOperator<>::Ptr precond = makeSparseCholeskySolver( matParam );
    
    gsInfo << "done.\n";

    gsInfo << "Setup iterative solver (conjugate gradient) for the physical domain..." << std::flush;
    gsConjugateGradient<> solver(matPhys, precond);
    solver.setMaxIterations(500);
    solver.setTolerance(1e-12);
    solver.setCalcEigenvalues(true);
    gsInfo << "done.\n";

    gsMatrix<> solVector;
    gsMatrix<> rhs;
    
    solVector.setRandom(matParam.rows(),1);
    rhs.setRandom(matParam.rows(),1);
    gsMatrix<> iters, eigs;
    
    gsInfo << "Solve..." << std::flush;
    solver.solveDetailed(rhs,solVector,iters);
    gsInfo << "done.\n\n";


    gsInfo << "Number of interations: " << solver.iterations() << "\n\n";
    solver.getEigenvalues(eigs);
    gsInfo << "Eigenvalues: " << eigs.transpose() << "\n\n";
    gsInfo << "Error: " << iters.transpose() << "\n\n";

    real_t thereshold = iters(0,0) * 1.e-8;
    index_t l = 0;
    for( l=0; iters(l,0) > thereshold && l<iters.size(); ++l ) {}
    if ( iters(l,0) <= thereshold ) 
        gsInfo << "Required steps to reduce initial error by a factor of 1.e-8: " << l << "\n\n";
    else
        gsInfo << "Required steps to reduce initial error by a factor of 1.e-8: More than applied.\n\n";

    gsInfo << "Condition number: " << ( eigs(eigs.rows()-1,0) / eigs(0,0) ) << "\n\n";

    bool plot = false;
    if (plot)
    {
        // Plotting the physical geometry
        gsWriteParaview<>( *geoPhys, "BiharmonicPhysicalGeometry");
    }

    return  0;
}

