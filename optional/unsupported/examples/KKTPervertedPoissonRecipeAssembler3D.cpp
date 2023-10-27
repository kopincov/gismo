/** @file KKTPervertedPoisson.cpp

    @brief A optimal controll with PDE constrain example in 3D.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <gismo.h>
#include <gismo_dev.h>
#include <gsAssembler/gsKKTPervertedPoisson.h>
#include <gsAssembler/gsBiharmonicAssembler.h>
#include <gsMultiGrid/gsGridHierarchy.h>
#include <gsSolver/gsSolverUtils.h>
#include <gsSolver/gsBlockOp.h>
#include <gsMapUtils/gsL2GMapper.h>
#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsRecipeAssembler/gsRecipeAssemblerKKTPervertedPossion.h>
#include <gsRecipeAssembler/gsRecipeAssemblerBiharmonicSogn.h>
#include <gsPde/gsBiharmonicPde.h>
#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsSolver/gsKronecker.h>
#include <gsMultiGrid/gsBlockSmoother.h>




using namespace gismo;


gsMatrix<> solve(const gsSparseMatrix<> &sys, const gsMatrix<> &rhs)
{
    gsSparseSolver<>::QR  solver;
    solver.analyzePattern( sys );
    solver.factorize     ( sys );
    return solver.solve( rhs );
}

template <typename T>
T* Cloner(const T* a)
{
    return memory::convert_ptr<T>(a->clone()).release();
}

template <typename T>
T* Constifier(T* a)
{
    return a;//const_cast<const T*>(a);
}

void DegreeElevate( gsBasis<real_t>* a) {a->degreeElevate();}

void UniformRefine( gsBasis<real_t>* a) {a->uniformRefine();}

void DegreeReduce( gsBasis<real_t>* a) {a->degreeReduce();}

void ReduceContinuity( gsBasis<real_t>* a) {a->reduceContinuity();}


std::vector<gsPhysicalSpace*> constructSpaces (
        const std::vector<gsBasis<real_t> *>       &bases,
        const gsMultiPatch<real_t>                 &geo,
        std::vector<std::vector<gsBasis<>*> >      *outBasis,
        index_t                                    numRefine = 2,
        index_t                                    increasePolydegree = 0,
        bool                                       reduceCont = false
        )
{
    std::vector<gsPhysicalSpace*>        result(3);

    std::vector<gsBasis<real_t>*>       temp,  patchBasisC, patchBasisS, patchBasisM;
    index_t increasePolydegree2 = increasePolydegree +2;
    {   // build control space

        std::transform(bases.begin(),bases.end(),std::back_inserter(patchBasisC), Cloner<gsBasis<real_t> >);

        if (increasePolydegree2 == 0)
        {
            std::for_each(patchBasisC.begin(),patchBasisC.end(), DegreeReduce);
            std::for_each(patchBasisC.begin(),patchBasisC.end(), DegreeReduce);
        }
        else if (increasePolydegree2 == 1)
            std::for_each(patchBasisC.begin(),patchBasisC.end(), DegreeReduce);
        else
        {
            for (index_t k = 0; k < increasePolydegree2-2; ++k)
                std::for_each(patchBasisC.begin(),patchBasisC.end(), DegreeElevate);
        }

        for (index_t k = 0; k < numRefine; ++k)
            std::for_each(patchBasisC.begin(),patchBasisC.end(), UniformRefine);

        //Reduce Continuity
        if (reduceCont)
        {
            std::for_each(patchBasisC.begin(),patchBasisC.end(), ReduceContinuity);
            std::for_each(patchBasisC.begin(),patchBasisC.end(), ReduceContinuity);
        }

        temp.resize(patchBasisC.size());
        cloneAll(patchBasisC.begin(),patchBasisC.end(),temp.begin());
        gsMultiBasis<real_t>                 multipatch (temp, geo);
        gsMapFactoryMultiBasis               mapFactory (multipatch);
        std::vector< gsBasis<real_t>*>       argVec;

        std::transform(patchBasisC.begin(),patchBasisC.end(),std::back_inserter(argVec), Constifier<gsBasis<real_t> >);
        result[controlSpace]=new gsPhysicalSpaceScalar(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
    }

    {   // build state space

        std::transform(bases.begin(),bases.end(),std::back_inserter(patchBasisS), Cloner<gsBasis<real_t> >);


        for (index_t k = 0; k < increasePolydegree; ++k)
            std::for_each(patchBasisS.begin(),patchBasisS.end(), DegreeElevate);

        for (index_t k = 0; k < numRefine; ++k)
            std::for_each(patchBasisS.begin(),patchBasisS.end(), UniformRefine);

        temp.resize(patchBasisS.size());
        cloneAll(patchBasisS.begin(),patchBasisS.end(),temp.begin());
        gsMultiBasis<real_t>                 multipatch (temp,geo);
        gsMapFactoryMultiBasis               mapFactory (multipatch);
        std::vector<gsBasis<real_t>*>        argVec;

        std::transform(patchBasisS.begin(),patchBasisS.end(),std::back_inserter(argVec), Constifier<gsBasis<real_t> >);
        result[stateSpace]=new gsPhysicalSpaceScalar(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
    }
    {   // build multiplier space

        std::transform(bases.begin(),bases.end(),std::back_inserter(patchBasisM), Cloner<gsBasis<real_t> >);

        if (increasePolydegree2 == 0)
        {
            std::for_each(patchBasisM.begin(),patchBasisM.end(), DegreeReduce);
            std::for_each(patchBasisM.begin(),patchBasisM.end(), DegreeReduce);
        }
        else if (increasePolydegree2 == 1)
            std::for_each(patchBasisM.begin(),patchBasisM.end(), DegreeReduce);
        else
        {
            for (index_t k = 0; k < increasePolydegree2-2; ++k)
                std::for_each(patchBasisM.begin(),patchBasisM.end(), DegreeElevate);
        }

        for (index_t k = 0; k < numRefine; ++k)
            std::for_each(patchBasisM.begin(),patchBasisM.end(), UniformRefine);

        //Reduce Continuity
        if (reduceCont)
        {
            std::for_each(patchBasisM.begin(),patchBasisM.end(), ReduceContinuity);
            std::for_each(patchBasisM.begin(),patchBasisM.end(), ReduceContinuity);
        }

        temp.resize(patchBasisM.size());
        cloneAll(patchBasisM.begin(),patchBasisM.end(),temp.begin());
        gsMultiBasis<real_t>                 multipatch (temp,geo);
        gsMapFactoryMultiBasis               mapFactory (multipatch);
        std::vector<gsBasis<real_t>*>        argVec;

        std::transform(patchBasisM.begin(),patchBasisM.end(),std::back_inserter(argVec), Constifier<gsBasis<real_t> >);
        result[multiplierSpace]=new gsPhysicalSpaceScalar(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
    }


    if(outBasis)
    {
        *outBasis= std::vector<std::vector<gsBasis<>*> >(1,patchBasisC);
        outBasis->push_back(patchBasisM);
        outBasis->push_back(patchBasisS);
        //*outBasis= std::vector<std::vector<gsBasis<>*> >(1,patchBasisS);
        //outBasis->push_back(patchBasisM);
        //outBasis->push_back(patchBasisC);
    }
    return result;
}


int main(int argc, char *argv[])
{
    gsStopwatch totalTime;
    gsStopwatch clock;
    gsStopwatch totalTimeAssembler;
    enum 
    {
        controlSpace = 0,
	stateSpace = 2,
	multiplierSpace = 1
    };
    index_t numRefine = 3;
    index_t numDegree = 0;
    real_t alpha = 1.0;
    if (argc >= 2)
        numRefine = atoi(argv[1]);
    if (argc >= 3)
        numDegree = atoi(argv[2]);
    if (argc >= 4)
        alpha = atof(argv[3]);

    gsInfo << "numRefine is: "<< numRefine <<" degree is: "<< numDegree + 3<<" Alpha is: " << alpha << "\n";

    gsInfo << std::setprecision(6);
    //gsInfo<< std::scientific;


    gsFunctionExpr<> zeroFunc("0.0",3);
    gsFunctionExpr<> solSin ("sin(2*pi*x) * sin(4*pi*y) * sin(6*pi*z)",3);
    gsFunctionExpr<> solSinGrad ("2*pi*cos(2*pi*x)*sin(4*pi*y)*sin(6*pi*z)",
                                 "4*pi*sin(2*pi*x)*cos(4*pi*y)*sin(6*pi*z)",
                                 "6*pi*sin(2*pi*x)*sin(4*pi*y)*cos(6*pi*z)",3);

    //gsFunctionExpr<> fSin("-20*pi**2*sin(2*pi*x)*sin(4*pi*y)",3);

    //Container for the desired state at the boundary. (condition type does not matter)
    gsBoundaryConditions<> desired;

    desired.addCondition( boundary::west, condition_type::neumann, &solSinGrad); //West
    desired.addCondition( boundary::east, condition_type::neumann, &solSinGrad); //East
    desired.addCondition( boundary::south,condition_type::neumann, &solSinGrad); //South
    desired.addCondition( boundary::north,condition_type::neumann, &solSinGrad); //North
    desired.addCondition( boundary::front,condition_type::neumann, &solSinGrad); //Front
    desired.addCondition( boundary::back, condition_type::neumann, &solSinGrad); //Back

    //Read 3D geometry
    std::string input("volumes/twistedFlatQuarterAnnulus.xml");
    gsFileData<> fileData(input);
    gsMultiPatch<> geo;
    fileData.getFirst< gsMultiPatch<> >(geo);

    //gsMultiPatch<> geo = gsNurbsCreator<>::BSplineCube(3);

    gsMultiBasis<> geoBases(geo);
    //cubic splines in all directions
    geoBases.degreeElevate(2,0);
    geoBases.degreeElevate(2,2);


    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &zeroFunc);//Annulus: small arch lenght
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &zeroFunc);//Annulus: Large arch lenght
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &zeroFunc);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &zeroFunc);
    bcInfo.addCondition( boundary::front, condition_type::dirichlet, &zeroFunc);
    bcInfo.addCondition( boundary::back,  condition_type::dirichlet, &zeroFunc);


    //PDE info container
    gsPoissonPde<real_t> pde(geo, bcInfo, zeroFunc, solSin.clone().release());

    // Init assembler
    gsRecipeAssemblerKKTPervertedPoisson assembler(pde, desired, alpha);

    //Contruct discretisation spaces
    std::vector<std::vector<gsBasis<>*> > bases;
    std::vector<gsBasis<>*> tmpBasisVec;
    for (size_t i=0; i< pde.domain().nPatches();++i)
        tmpBasisVec.push_back(const_cast<gsBasis<>*>(&geoBases.basis(i)));

    std::vector<gsPhysicalSpace*> phySpace=constructSpaces(tmpBasisVec,pde.domain(),&bases,numRefine,numDegree, false);
    assembler.setSpace(phySpace);
    //gsInfo << phySpace[0]->getMapper()->getNrOfTargets() << " X "<< phySpace[0]->getMapper()->getNrOfSources()<<"\n";
    //gsInfo << phySpace[1]->getMapper()->getNrOfTargets() << " X "<< phySpace[1]->getMapper()->getNrOfSources()<<"\n";
    //gsInfo << phySpace[2]->getMapper()->getNrOfTargets() << " X "<< phySpace[2]->getMapper()->getNrOfSources()<<"\n";
    //std::exit(0);
    clock.restart();
    assembler.assemble();
    gsInfo << "time to assemble: " << clock.stop() << "\n"; 

    gsMatrix<>         eli;

    gsSparseMatrix<> sysSparse = assembler.getSystemMatrix();
    //sysSparse.finalize();
    gsMatrix<> rhs = assembler.getSystemRhs();
    eli  = solve(assembler.getEliminatedMatrix(),assembler.getEliminatedRhs());
    rhs -= assembler.getRhsModMatrix()*eli;

    gsInfo << "Total time for the assembler: " << totalTimeAssembler.stop() <<"sec\n";

    gsInfo << "System size " << assembler.getFreeLimit() << ": ";
    //gsInfo << getElimSize()  << " Number of shifts: " << getShifts().size()<< "\n";

    index_t fSz = assembler.getShifts()[1];
    index_t wSz = assembler.getShifts()[2] - assembler.getShifts()[1];
    index_t uSz = assembler.getShifts()[3] - assembler.getShifts()[2] - assembler.getElimSize();

    index_t totSz = assembler.getFreeLimit();
    //index_t uSz = assembler.getShifts()[1] - assembler.getElimSize();
    //index_t fSz = assembler.getShifts()[2] - assembler.getShifts()[1];
    //index_t wSz = assembler.getShifts()[3] - assembler.getShifts()[2];
    gsInfo <<  fSz << " X " << wSz << " X " << uSz<< "\n";



    /////////////////// Assembling done -> setup solver ///////////////////

    //sysSparse.prune(1,1e-16);
    gsStopwatch totalTimeSolver;
    real_t alphaInv = 1.0/alpha;
    //gsMatrix<> mass = sysSparse.block(fSz,0,wSz,fSz);
    gsSparseMatrix<> massSparse = sysSparse.block(fSz,0,wSz,fSz);
    //gsMatrix<> Amat = sysSparse.block(fSz,fSz+wSz,wSz,uSz);
    gsSparseMatrix<> obs  = sysSparse.block(fSz+wSz,fSz+wSz,uSz,uSz);
    //gsMatrix<> zeroMat = sysSparse.block(fSz+uSz,fSz+uSz,wSz,wSz);
    //gsMatrix<> mass1 = sysSparse.block(0,0,fSz,fSz);

    //gsInfo << mass1<< "\n";
    //return 0;

    //gsMatrix<> denseMat = sysSparse;
    //gsMatrix<> mass = denseMat.block(fSz+uSz,0,fSz,fSz);
    //gsMatrix<> Amat = denseMat.block(fSz+uSz,fSz,wSz,uSz);
    //gsMatrix<> obs = denseMat.block(fSz,fSz,uSz,uSz);

    //Kronecker product mass matrix
    /*
    gsTensorBSplineBasis<3,real_t> * bs = 0;
    bs = dynamic_cast<gsTensorBSplineBasis<3,real_t>*>(bases[0][0]);
    gsSparseMatrix<> Mx, My, Mz;

    assembleParameterMass(bs->component(0), Mx);
    assembleParameterMass(bs->component(1), My);
    assembleParameterMass(bs->component(2), Mz);


    handleDirichletConditions(Mx,desired,boundary::west ,boundary::east);
    handleDirichletConditions(My,desired,boundary::south,boundary::north);
    handleDirichletConditions(Mz,desired,boundary::front,boundary::back);
    gsKroneckerProduct * massKronecker = new gsKroneckerProduct( makeSparseCholeskySolver(Mx), makeSparseCholeskySolver(My), makeSparseCholeskySolver(Mz) );
    */
    //

    //gsInfo << "Inverting mass matrix" << "\n";

    //gsMatrix<> massInv = mass.inverse();
    //gsInfo << "Inverting complete! Time: "<< clock.stop() << "sec"  << "\n";

    //new gsSolverOperator<typename  gsSparseSolver<T>::SimplicialLDLT >(mass)

    //gsSparseMatrix<> massInvAlphaInv((alphaInv*massInv).sparseView());
    //gsMatrix<> massInvAlphaInv(alphaInv*massInv);

    //gsSparseMatrix<> massInvAlphaInv(massInv.sparseView());

    //gsSparseMatrix<> massInvAlpha((alpha*massInv).sparseView());
    //gsMatrix<> massInvAlpha(alpha*massInv);

    /////////////////// Assembling the Biharmonic equation ///////////////////

    gsBiharmonicPde<real_t> pdeBiharmonic(geo,bcInfo,gsBoundaryConditions<>(),zeroFunc);
    gsRecipeAssemblerBiharmonicSogn assemblerBiharmonic(pdeBiharmonic);
    assemblerBiharmonic.setDirichletStrategy(dirichlet::elimination);
    std::vector<gsPhysicalSpace*> phySpaceBiharmonic;
    phySpaceBiharmonic.push_back(phySpace[stateSpace]);
    assemblerBiharmonic.setSpace(phySpaceBiharmonic);

    assemblerBiharmonic.assemble();
    gsSparseMatrix <> BiharMat = assemblerBiharmonic.getSystemMatrix();

    //Condition Number!
    /*
    gsMatrix <> AtInvMA = Amat.transpose()*massInv*Amat;
    gsMatrix <> invAtInvMA = AtInvMA.inverse();
    gsInfo << "The condition number is B*inv(A*inv(M)A): "
         << gsSolverUtils<real_t>::conditionNumber(BiharMat*invAtInvMA, false) << "\n";
    //std::exit(0);*/

    //gsSparseMatrix<> mgMatrix((alpha*(Amat.transpose()*massInv*Amat) + obs).sparseView());
    gsSparseMatrix<> mgMatrix((alpha*BiharMat + obs));
    //clock.restart();
    //gsMatrix<> mgInverse( ((alpha*(Amat.transpose()*massInv*Amat) + obs).inverse()) );
    //gsInfo<< "S_3 inverse time: "<< clock.stop() << "sec"  << "\n";
    /*
    //Setting up Multigrid preconditioner
    gsInfo << "Setting up multigrid" << "\n";
    std::vector< gsMultiBasis<real_t>* > basisLvls;
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;
    gsGridHierarchy<>::buildByCoarsening(gsMultiBasis<>(fineBasis), gsBoundaryConditions<>(), gsAssembler<>::defaultOptions, 100, 60)
        .moveMultiBasesTo(basisLvls)
        .moveTransferMatricesTo(transferMatrices)
        .clear();

    //Create DoF mapper for each level
    const int numTransfer = basisLvls.size() - 1;
    std::vector< gsDofMapper > dofMappers(numTransfer+1);
    const bool conforming =  true ;

    //Get DoF-mappers
    for (int k = 0; k < numTransfer+1; ++k)
    {
        if ( true )
            basisLvls[k]->getMapper(conforming, bcInfo, dofMappers[k] );
        else
            basisLvls[k]->getMapper(conforming, dofMappers[k] );
    }

    // Memory clean
    freeAll(basisLvls);


    // incorporate boundary conditions
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatricesWithBC;
    restrictToFreeDofs( transferMatrices, dofMappers, transferMatricesWithBC );
    gsMultiGridOp<> mg(mgMatrix,transferMatricesWithBC);

    std::vector<std::vector<gsVector<index_t> > > blockInfo;

    for (int i = 0; i < mg.numLevels(); ++i)
    {
        blockInfo.push_back(gsBlockOrgaizer::uniformBlocks(
        int (std::sqrt((double)dofMappers[i].freeSize())), dofMappers[i].freeSize()) ); //std::cbrt(n) c++11
        //mg.setSmoother(i, makeGaussSeidelOp(mg.matrix(i)));
        mg.setSmoother(i, new gsGaussSeidelBlockSmoother(mg.matrix(i),blockInfo[i]));
    }

    mg.setNumPreSmooth(2);
    mg.setNumPostSmooth(2);
    //End inite mg */

    gsBlockOp<>::Ptr blockPrec = gsBlockOp<>::make(3,3);
    // TODO: here is some code missing, something like (KKTPervertedPoissonRecipeAssembler.cpp:391:
    //       blockPrec->addOperator(0,0,makeMatrixOp(massInvAlphaInv));
    //       blockPrec->addOperator(1,1,makeMatrixOp(massInvAlpha));
    //       blockPrec->addOperator(2,2,makeMatrixOp(mgInverse));
    //       blockPrec->rows() must be 3391
    //       maybe be 1109 * 3 + 64
    gsMinimalResidual<> solver(sysSparse, blockPrec);
    solver.setMaxIterations(2000);
    solver.setTolerance(1e-12);
    
    gsMatrix<> solVector;
    solVector.setZero(totSz,1);
    //solVector = solVector/solVector.col(0).norm();

    //blockPrec->addPreconditioner(massKronecker,0,0);
    //blockPrec->addPreconditioner(massKronecker,1,1);
    clock.restart();
    gsLinearOperator<>::Ptr massCholeskyInv = makeSparseCholeskySolver(massSparse);
    //gsLinearOperator * massCholeskyInv = makeSparseCholeskySolver(massSparse);
    gsInfo<< "time to make 1 Cholesky preconditioner " << clock.stop() <<"sec\n";
    //gsScaledOp massAI(*massCholeskyInv, alphaInv);
    //gsScaledLinearOperator massAI(*makeSparseCholeskySolver(massSparse), alphaInv);
    //gsScaledOp massA (*massCholeskyInv, alpha);
    blockPrec->addOperator(controlSpace, controlSpace, gsScaledOp<>::make(massCholeskyInv, alphaInv));
    blockPrec->addOperator(multiplierSpace, multiplierSpace, gsScaledOp<>::make(massCholeskyInv, alpha));

    //blockPrec->addPreconditioner(makeMatrixOperator(massInvAlphaInv),0,0);
    //blockPrec->addPreconditioner(makeMatrixOperator(massInvAlpha),1,1);

    clock.restart();
    blockPrec->addOperator(stateSpace, stateSpace, makeSparseCholeskySolver(mgMatrix));
    gsInfo<< "time to make 3 Cholesky preconditioner " << clock.stop() <<"sec\n";
    //blockPrec->addPreconditioner(&mg,2,2);

    gsInfo << "---------------------" << "\n";
    //gsMatrixOp<gsMatrix<> > precCononical(preconMat);
    clock.restart();
    solver.solve(rhs, solVector);
    gsInfo << "Time to solve with MINRES: " << clock.stop() <<"sec\n";
    gsInfo << "Total time for the solver: " << totalTimeSolver.stop() <<"sec\n";
    gsInfo << "Number of iterations " << solver.iterations() << "\n";
    gsInfo << "Residual Error " << solver.error() << "\n";
    gsInfo << "---------------------" << "\n";



    gsInfo << "Test is done: Cleaning up..." << "\n";
    freeAll(phySpace);
    //delete massCholeskyInv;
    gsInfo << "Exiting, total time: " << totalTime.stop() << " sec\n";
    return  0;

}


