/** @file gsPervertedStokes.cpp

    @brief Tests a perverted Stokes recipe assembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <iostream>
#include <iomanip>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsCore/gsField.h>
#include <gsCore/gsFunctionExpr.h>

#include <gsAssembler/gsAssemblerUtils.h>

#include <gsRecipeAssembler/gsRecipeAssembler.h>
#include <gsRecipeAssembler/gsRecipeAssemblerPervertedStokes.h>
#include <gsRecipeAssembler/gsRecipeAssemblerBiharmonicSogn.h>
#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>

#include <gsNurbs/gsNurbsCreator.h>
#include <gsPde/gsPervertedStokesPde.h>

#include <gsSolver/gsStokesIterativeSolver.h>



#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n");\

bool passed=true;


;
using namespace gismo;

gsMatrix<> solveDirect(const gsSparseMatrix<> &sys, const gsMatrix<> &rhs)
{
    gsSparseSolver<>::QR  solver;
    solver.analyzePattern( sys );
    solver.factorize     ( sys );
    return solver.solve( rhs );
}

class gsRecipeAssemblerMass : public gsRecipeAssembler
{
public:
    gsRecipeAssemblerMass(const gsMultiPatch<real_t> &domain)
        : gsRecipeAssembler(domain)
    {

    }

protected:
    // recipe providing functions
    virtual gsRecipe<real_t>    getPatchRecipe     (index_t patch)
    {
        gsRecipe<real_t>            result;
        gsRecipeIngredient<real_t>  ingr;

        ingr.setOperator(new gsL2ScalarOp<real_t>());
        ingr.setTestSpace(0);
        ingr.setUnknownSpace(0);
        ingr.setRule(getSysWriter());
        result.add(ingr);
        return result;
    }

    virtual gsRecipe<real_t>    getBoundaryRecipe  (patchSide ps)
    {
        gsRecipe<real_t> result;
        return result;
    }
};


int main(int argc, char *argv[])
{

    gsInfo.width(9);
    gsInfo<< std::setprecision(2);
    gsInfo<< std::scientific;

    index_t numRefine = 2;
    index_t numDegree = 1;
    if (argc >= 2)
        numRefine = atoi(argv[1]);
    if (argc >= 3)
        numDegree = atoi(argv[2]);


    gsFunctionExpr<> f("+2*cos(x+y)+2*sin(x-y)","-2*cos(x+y)+2*sin(x-y)",2);

    gsMultiPatch<>        *domain;

    domain = new gsMultiPatch<real_t>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    //      gsReadFile<>("square.xml" )),
    //  new gsMultiPatch<real_t>(*gsNurbsCreator<real_t>::BSplineSquareDeg(1))),
    //      gsReadFile<>("surfaces/multipatch_triangle2.xml")),
    gsMultiBasis<> basis =gsMultiBasis<>(*domain);

    //p-refine to get equal polynomial degree s,t directions (for Annulus)
    basis.degreeElevate(1,0);

    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> zeroFunc("0.0","0.0",2);

    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &zeroFunc);//Annulus: small arch lenght
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &zeroFunc);//Annulus: Large arch lenght
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &zeroFunc);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &zeroFunc);

    //add BC info here!
    //gsPervertedStokesPde<real_t> pde(*domain,gsBoundaryConditions<>(),f.clone());
    gsPervertedStokesPde<real_t> pde(*domain, bcInfo, f.clone().release());

    gsSparseMatrix<real_t> sys;
    gsMatrix<real_t>       rhs;
    gsMatrix<real_t>       eli;
    gsMatrix<real_t>       sol;
    gsVector<index_t>      shifts(5);
    std::vector<gsPhysicalSpace*> phySpace;

    /// Create discretization spaces
    std::vector<std::vector<gsBasis<>*> > bases;
    gsRecipeAssemblerPervertedStokes assembler(pde);
    std::vector<gsBasis<>*> tmpBasisVec;
    for (size_t i=0; i< pde.domain().nPatches();++i)
        tmpBasisVec.push_back(const_cast<gsBasis<>*>(&basis.basis(i)));

    phySpace=constructPervertedTHSpaces(tmpBasisVec,pde.domain(),&bases,numRefine,numDegree);//Put in main

    /// Assemble Perverted stokes system
    assembler.setSpace(phySpace);
    assembler.assemble();

    sys = assembler.getSystemMatrix();
    rhs = assembler.getSystemRhs();
    /*if (pde.boundaryConditions().dirichletSides().size()>0)
    {
        eli  = solveDirect(assembler.getEliminatedMatrix(),assembler.getEliminatedRhs());
        rhs -= assembler.getRhsModMatrix()*eli;
    }
    else
        eli.resize(0,rhs.cols());*/

    gsInfo << "System size: " << assembler.getFreeLimit() << "\n";
    gsInfo << "Eliminated size: " << assembler.getElimSize() << "\n";

    //Print the space sizes
    index_t shiftCounter = 0;
    shifts[0] = assembler.getShifts()[1] - assembler.getElimSize();
    gsInfo << "State velocity size: " << shifts[0] << " (DoFs removed)" << "\n";
    shiftCounter = assembler.getShifts()[1];
    shifts[1] = assembler.getShifts()[2] - shiftCounter;
    gsInfo << "State pressure size: " << shifts[1]<< "\n";
    shiftCounter = assembler.getShifts()[2];
    shifts[2] = assembler.getShifts()[3] - shiftCounter;
    gsInfo << "Multiplier velocity size: " << shifts[2] << "\n";
    shiftCounter = assembler.getShifts()[3];
    shifts[3] = assembler.getShifts()[4] - shiftCounter;
    gsInfo << "Multiplier pressure size: " << shifts[3] << "\n";
    shifts[4] = assembler.getFreeLimit();

    ///------------------Construct preconditioner-------------------///
    //Biharmonic for state velocity
    gsBiharmonicPde<real_t> pdeBiharmonic(*domain,bcInfo,gsBoundaryConditions<>(),zeroFunc);
    gsRecipeAssemblerBiharmonicSogn assemblerBiharmonic(pdeBiharmonic);
    assemblerBiharmonic.setDirichletStrategy(dirichlet::elimination);
    std::vector<gsPhysicalSpace*> phySpaceBiharmonic;
    phySpaceBiharmonic.push_back(phySpace[0]);
    assemblerBiharmonic.setSpace(phySpaceBiharmonic);
    assemblerBiharmonic.assemble();
    gsSparseMatrix<> BiharMat = assemblerBiharmonic.getSystemMatrix();
    gsInfo << "Constructed biharmonic matrix! Dimension: "<< BiharMat.rows()<<" X "  << BiharMat.cols() << "\n";

    //Poisson for state and multiplier pressure
    gsPoissonPde<real_t> pdePoisson(*domain,gsBoundaryConditions<>(),zeroFunc);
    gsRecipeAssemblerPoisson assemblerPoisson(pdePoisson);
    assemblerPoisson.setDirichletStrategy(dirichlet::none);
    std::vector<gsPhysicalSpace*> phySpacePoisson;
    phySpacePoisson.push_back(phySpace[1]);
    assemblerPoisson.setSpace(phySpacePoisson);
    assemblerPoisson.assemble();
    gsSparseMatrix<> poissonMat = assemblerPoisson.getSystemMatrix();
    gsInfo << "Constructed Poisson matrix! Dimension: "<< poissonMat.rows()<<" X "  << poissonMat.cols() << "\n";

    //Mass for multiplier velocity
    gsRecipeAssemblerMass assemblerMass(*domain);
    std::vector<gsPhysicalSpace*> phySpaceMass;
    phySpaceMass.push_back(phySpace[2]);
    assemblerMass.setSpace(phySpaceMass);
    assemblerMass.assemble();
    gsSparseMatrix<> massMat = assemblerMass.getSystemMatrix();
    gsInfo << "Constructed mass matrix! Dimension: "<< massMat.rows()<<" X "  << massMat.cols() << "\n";


    //Clean up
    for (size_t s=0; s<bases.size();++s)
        freeAll(bases[s]);
    freeAll(phySpace);
    delete domain;
    return !passed;
}



