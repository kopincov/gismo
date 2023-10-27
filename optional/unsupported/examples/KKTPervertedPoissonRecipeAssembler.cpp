/** @file KKTPervertedPoisson.cpp

    @brief A optimal controll with PDE constrain example.

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
        result[0]=new gsPhysicalSpaceScalar(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
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
        result[2]=new gsPhysicalSpaceScalar(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
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
        result[1]=new gsPhysicalSpaceScalar(argVec,geo, INVERSE_COMPOSITION, *mapFactory.makeMapper());
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
    index_t numRefine = 3;
    index_t numDegree = 0;
    real_t alpha = 1.0;
    if (argc >= 2)
        numRefine = atoi(argv[1]);
    if (argc >= 3)
        numDegree = atoi(argv[2]);
    if (argc >= 4)
        alpha = atof(argv[3]);

    gsInfo << "numRefine is: "<< numRefine <<" degree is: "<< numDegree + 2<<" Alpha is: " << alpha << "\n";

    gsInfo<< std::setprecision(12);

    gsFunctionExpr<> zeroFunc("0.0",2);
    gsFunctionExpr<> solSin ("sin(2*pi*x) * sin(4*pi*y)",2);
    gsFunctionExpr<> solSinGrad ("2*pi*cos(2*pi*x)*sin(4*pi*y)", "4*pi*sin(2*pi*x)*cos(4*pi*y)",2);
    gsFunctionExpr<> fSin("-20*pi**2*sin(2*pi*x)*sin(4*pi*y)",2);

    //Container for the desired state at the boundary. (condition type does not matter)
    gsBoundaryConditions<> desired;
    desired.addCondition( boundary::west, condition_type::neumann, & solSinGrad ); //West
    desired.addCondition( boundary::east, condition_type::neumann, & solSinGrad ); //East
    desired.addCondition( boundary::south,condition_type::neumann, & solSinGrad ); //South
    desired.addCondition( boundary::north,condition_type::neumann, & solSinGrad ); //North

    gsMultiPatch<> * geo = NULL;

    //geo = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineSquareDeg(2) );
    geo = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );

    gsMultiBasis<> geoBases(*geo);
    geoBases.degreeElevate(1,0);

    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &zeroFunc);//Annulus: small arch lenght
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &zeroFunc);//Annulus: Large arch lenght
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &zeroFunc);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &zeroFunc);



    //PDE info container

    gsPoissonPde<real_t> pde(*geo, bcInfo, zeroFunc, solSin.clone().release());

    // Init assembler
    gsRecipeAssemblerKKTPervertedPoisson assembler(pde, desired, alpha);

    //Contruct discretisation spaces
    std::vector<std::vector<gsBasis<>*> > bases;
    std::vector<gsBasis<>*> tmpBasisVec;
    for (size_t i=0; i< pde.domain().nPatches();++i)
        tmpBasisVec.push_back(const_cast<gsBasis<>*>(&geoBases.basis(i)));

    std::vector<gsPhysicalSpace*> phySpace=constructSpaces(tmpBasisVec,pde.domain(),&bases,numRefine,numDegree, true);
    assembler.setSpace(phySpace);

    assembler.assemble();

    gsMatrix<>         eli;

    gsSparseMatrix<> sysSparse = assembler.getSystemMatrix();
    //sysSparse.finalize();
    gsMatrix<> rhs = assembler.getSystemRhs();
    eli  = solve(assembler.getEliminatedMatrix(),assembler.getEliminatedRhs());
    rhs -= assembler.getRhsModMatrix()*eli;

    gsInfo << "System size " << assembler.getFreeLimit() << "\n";
    //gsInfo << getElimSize()  << " Number of shifts: " << getShifts().size()<< "\n";

    index_t fSz = assembler.getShifts()[1];
    index_t wSz = assembler.getShifts()[2] - assembler.getShifts()[1];
    index_t uSz = assembler.getShifts()[3] - assembler.getShifts()[2] - assembler.getElimSize();

    index_t totSz = assembler.getFreeLimit();
    //index_t uSz = assembler.getShifts()[1] - assembler.getElimSize();
    //index_t fSz = assembler.getShifts()[2] - assembler.getShifts()[1];
    //index_t wSz = assembler.getShifts()[3] - assembler.getShifts()[2];
    gsInfo << "Control    size: " << fSz << "\n";
    gsInfo << "Multiplier size: " << wSz << "\n";
    gsInfo << "State      size: " << uSz<< "\n";



    /////////////////// Assembling done -> setup solver ///////////////////

    //sysSparse.prune(1,1e-16);

    real_t alphaInv = 1.0/alpha;
    gsMatrix<> mass = sysSparse.block(fSz,0,wSz,fSz);
    gsMatrix<> Amat = sysSparse.block(fSz,fSz+wSz,wSz,uSz);
    gsMatrix<> obs  = sysSparse.block(fSz+wSz,fSz+wSz,uSz,uSz);
    //gsMatrix<> zeroMat = sysSparse.block(fSz+uSz,fSz+uSz,wSz,wSz);
    //gsMatrix<> mass1 = sysSparse.block(0,0,fSz,fSz);

    //gsInfo << mass1<< "\n";
    //return 0;

    //gsMatrix<> denseMat = sysSparse;
    //gsMatrix<> mass = denseMat.block(fSz+uSz,0,fSz,fSz);
    //gsMatrix<> Amat = denseMat.block(fSz+uSz,fSz,wSz,uSz);
    //gsMatrix<> obs = denseMat.block(fSz,fSz,uSz,uSz);

    gsInfo << "Inverting mass matrix" << "\n";
    gsMatrix<> massInv = mass.inverse();
    gsInfo << "Inverting complete!" << "\n";

    //gsSparseMatrix<> massInvAlphaInv((alphaInv*massInv).sparseView());
    gsMatrix<> massInvAlphaInv(alphaInv*massInv);

    //gsSparseMatrix<> massInvAlphaInv(massInv.sparseView());

    //gsSparseMatrix<> massInvAlpha((alpha*massInv).sparseView());
    gsMatrix<> massInvAlpha(alpha*massInv);

    /////////////////// Assembling the Biharmonic equation ///////////////////

    gsBiharmonicPde<real_t> pdeBiharmonic(*geo,bcInfo,gsBoundaryConditions<>(),zeroFunc);
    gsRecipeAssemblerBiharmonicSogn assemblerBiharmonic(pdeBiharmonic);
    assemblerBiharmonic.setDirichletStrategy(dirichlet::elimination);
    std::vector<gsPhysicalSpace*> phySpaceBiharmonic;
    phySpaceBiharmonic.push_back(phySpace[2]);
    assemblerBiharmonic.setSpace(phySpaceBiharmonic);

    assemblerBiharmonic.assemble();
    gsMatrix <> BiharMat = assemblerBiharmonic.getSystemMatrix();

    //Condition Number!
    /*
    gsMatrix <> AtInvMA = Amat.transpose()*massInv*Amat;
    gsMatrix <> invAtInvMA = AtInvMA.inverse();
    gsInfo << "The condition number is B*inv(A*inv(M)A): "
         << gsSolverUtils<real_t>::conditionNumber(BiharMat*invAtInvMA, false) << "\n";
    //std::exit(0);*/

    //gsSparseMatrix<> mgMatrix((alpha*(Amat.transpose()*massInv*Amat) + obs).sparseView());
    //gsSparseMatrix<> mgMatrix((alpha*BiharMat + obs).sparseView());
    //gsMatrix<> mgInverse( ((alpha*(Amat.transpose()*massInv*Amat) + obs).inverse()) );
    //gsSparseMatrix<> mgMatrix((alpha*BiharMat + obs).sparseView());
    gsMatrix<> mgInverse((alpha*BiharMat + obs).inverse());


    //Condition Number!

    /*gsMatrix <> schurInv;
    schurInv.setZero(sysSparse.cols(), sysSparse.cols());
    schurInv.block(0,0,fSz,fSz) = massInvAlphaInv;
    schurInv.block(fSz,fSz,wSz,wSz) = massInvAlpha;
    schurInv.block(fSz+wSz,fSz+wSz,uSz,uSz) = mgInverse;
    gsInfo << "Starting kappa calculation\n";
    gsInfo << "The condition number is B*inv(A*inv(M)A): "
         << gsSolverUtils<real_t>::conditionNumber(schurInv*sysSparse, false) << "\n";
    std::exit(0);*/

    /*
    //Setting up Multigrid preconditioner
    gsInfo << "Setting up multigrid" << "\n";
    std::vector< gsMultiBasis<real_t> > basisLvls;
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;
    gsGridHierarchy<>::buildByCoarsening(gsMultiBasis<>(*bases[2][0]), bcInfo, gsAssembler<>::defaultOptions(), 100, 60)
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
            basisLvls[k].getMapper(conforming, bcInfo, dofMappers[k] );
        else
            basisLvls[k].getMapper(conforming, dofMappers[k] );
    }

    // Memory clean
    basisLvls.clear();
   
    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make(mgMatrix,transferMatrices);

    std::vector<std::vector<gsVector<index_t> > > blockInfo;

    //mg->setSmoother(i, makeGaussSeidelOp(mg->matrix(i)));
    for (int i = 0; i < mg->numLevels(); ++i)
    {
        blockInfo.push_back(gsBlockOrgaizer::uniformBlocks(
        int (std::sqrt((double)dofMappers[i].freeSize())), dofMappers[i].freeSize()) );
        mg->setSmoother(i, gsGaussSeidelBlockSmoother::make(mg->matrix(i),blockInfo[i]));
    }

    mg->setNumPreSmooth(2);
    mg->setNumPostSmooth(2);*/
    //End inite mg

    gsMatrix<> solVector;
    solVector.setZero(totSz,1);
    //gsInfo << solVector(5,0)<< ": " <<solVector(3,0) << "\n";
    //gsInfo << "norm random vector: " <<solVector.col(0).norm() << "\n";
    //solVector = solVector/solVector.col(0).norm();
    //gsInfo << "norm random vector: " <<solVector.col(0).norm() << "\n";

    gsBlockOp<>::Ptr blockPrec = gsBlockOp<>::make(3,3);   // blockPrec->rows() := 0
    blockPrec->addOperator(0,0,makeMatrixOp(massInvAlphaInv)); // blockPrec->rows() := 576 (+576)
    //gsMatrixOp<gsMatrix<> > prec11(mgInverse);
    //blockPrec->addOperator(1,1,&prec11);
    blockPrec->addOperator(1,1,makeMatrixOp(massInvAlpha));    // blockPrec->rows() := 1152 (*2)
    //blockPrec->addOperator(2,2,mg);
    blockPrec->addOperator(2,2,makeMatrixOp(mgInverse));       // blockPrec->rows() := 1216 (+64)

    gsMinimalResidual<> solver(sysSparse,blockPrec);
    solver.setMaxIterations(2000);
    solver.setTolerance(1e-12);
    
    gsInfo << "\n---------------------" << "\n";
    //gsMatrixOp<gsMatrix<> > precCononical(preconMat);
    solver.solve(rhs,solVector);
    gsInfo << "Number of iterations " << solver.iterations() << "\n";
    gsInfo << "Residual Error " << solver.error() << "\n";
    gsInfo << "---------------------\n" << "\n";



    //std::exit(0);

    /*
    Eigen::SparseQR< gsSparseMatrix<>, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(sysSparse );
    solver.factorize(sysSparse);
    gsDebug << "   Eigen  lastErrorMessage: " << solver.lastErrorMessage () << "\n";
    gsMatrix<> solVector= solver.solve(KKTPPAssembler.rhs());
    */

    gsInfo << "Test is done: Cleaning up..." << "\n"; //freeAll(m_bconditions);

    delete geo;
    //delete geoTest;
    //freeAll(desired);
    freeAll(phySpace);

    gsInfo << "Test is done: Exiting\n\n" << "\n";

    return  0;

}


