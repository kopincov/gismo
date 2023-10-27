/** @file gsRecipeAssemblerBiharmonicTest.cpp

    @brief This program tests that the gsRecipeAssembler and the gsBiharmonicAssembler
    assemble the same matrixes and right hand side.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, A. Bressan
*/
#include <iostream>
#include <iomanip>

#include <gismo.h>
#include <gismo_dev.h>
#include <gsUtils/gsStopwatch.h>

#include <gsCore/gsField.h>
#include <gsCore/gsFunctionExpr.h>

#include <gsAssembler/gsAssemblerUtils.h>
#include <gsAssembler/gsBiharmonicAssembler.h>
#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsSeminormH1.h>

#include <gsPde/gsBiharmonicPde.h>

#include <gsRecipeAssembler/gsRecipeAssemblerBiharmonicSogn.h>
#include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>

#include <gsNurbs/gsNurbsCreator.h>



;
using namespace gismo;

#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n");\

bool passed=true;

#ifdef NAN
    #undef NAN
#endif
static real_t NAN(std::numeric_limits<real_t>::quiet_NaN());


//Add a first kind of Dirichlet BC and second kind Neumann BC.
void addAllDirichletBoundaries(const gsMultiPatch<>& mp,
                               gsFunction<real_t> *g,
                               gsFunction<real_t> *laplace,
                               gsBoundaryConditions<real_t> & bcInfoFirst,
                               gsBoundaryConditions<real_t> & bcInfoSecond)
{
    std::vector<patchSide> boundaries = mp.boundaries();
    for(unsigned i = 0;i<boundaries.size();++i)
    {
        bcInfoFirst.addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::dirichlet, g);
        bcInfoSecond.addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::neumann, laplace);
    }
}

//gsMatrix<> solve(const gsSparseMatrix<> &sys, const gsMatrix<> &rhs)
//{
//    gsSparseSolver<>::CGDiagonal solver;
//    return solver.compute(sys).solve(rhs);
//}
gsMatrix<> solve(const gsSparseMatrix<> &sys, const gsMatrix<> &rhs)
{
    gsSparseSolver<>::QR  solver;
    solver.analyzePattern( sys );
    solver.factorize     ( sys );
    return solver.solve( rhs );
}

template <typename MatrixT>
real_t diff(const MatrixT &a,const MatrixT &b)
{
    if (a.rows()==b.rows() && a.cols()==b.cols() )
        return (a-b).norm();
    else
        return NAN;
}

class ProblemData
{
protected:
    gsMultiPatch<>         *domain;
    gsFunctionExpr<real_t> source;
    gsFunctionExpr<real_t> laplace;
public:
    gsFunctionExpr<real_t> solution;
    gsMultiBasis<>         basis;
    gsBiharmonicPde<real_t>   pde;

    dirichlet::strategy    Dstrategy;
    iFace::strategy        Istrategy;
public:
    ProblemData()
        :
          //domain(gsNurbsCreator<real_t>::BSplineSquareGrid(1, 1, 1.0)), //Multipatch those not work
          domain(new gsMultiPatch<>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus() )),
          source("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2),
          laplace("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2),
          solution("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2),
          basis(*domain),
          pde(*domain,gsBoundaryConditions<>(),gsBoundaryConditions<>(),source),
          Dstrategy(dirichlet::elimination),
          Istrategy(iFace::glue)

    {

        addAllDirichletBoundaries(*domain, &solution, &laplace, pde.boundaryConditions(),pde.boundaryConditionsSecond());
    }
    ~ProblemData()
    {
        delete domain;
    }
};


// assembling s
typedef void (*assembler) (const ProblemData&, gsSparseMatrix<>&,gsMatrix<>&, gsMatrix<>&, gsMatrix<>&, real_t&, real_t&);

void assembler1 (const ProblemData&, gsSparseMatrix<>&, gsMatrix<>&, gsMatrix<>&, gsMatrix<>&, real_t &, real_t &time);
void assembler2 (const ProblemData&, gsSparseMatrix<>&, gsMatrix<>&, gsMatrix<>&, gsMatrix<>&, real_t &, real_t &time);


assembler      methods[]={assembler1,assembler2};
const int numMet=sizeof(methods)/sizeof(assembler);

std::string    metNames  []={
    std::string("VisitorAssembler"),
    std::string("RecipeAssembler")
};

const int numRefine = 0;     // Lowest number of refinement: numRefine + 1
const int maxIterations = 3; // Highes number of refinement: numRefine + maxIterations

int main(int argc, char *argv[])
{
    gsInfo.width(9);
    gsInfo<< std::setprecision(2);
    gsInfo<< std::scientific;

    ProblemData data;

    data.basis.degreeElevate(1,0);
    data.basis.degreeElevate();

    for (int i = 0; i < numRefine; ++i)
        data.basis.uniformRefine();

    for (int i = 0; i < maxIterations; ++i)
    {
        data.basis.uniformRefine();

        gsMatrix<real_t>        tim;
        tim.setConstant(numMet,1,NAN);
        gsMatrix<real_t>        errL2;
        errL2.setConstant(numMet,1,NAN);
        //gsMatrix<real_t>        errH1;
        //errH1.setConstant(numMet,1,NAN);

        gsMatrix<real_t>        sysDiff, rhsDiff, eliDiff,solDiff;

        gsSparseMatrix<real_t> sys[numMet];
        gsMatrix<real_t>       rhs[numMet];
        gsMatrix<real_t>       eli[numMet];
        gsMatrix<real_t>       sol[numMet];

        for ( int m=0; m<numMet;++m)
        {
            methods[m](data,sys[m],rhs[m],eli[m],sol[m],errL2(m,0),tim(m,0));
        }


        sysDiff.setZero(numMet,numMet);
        rhsDiff.setZero(numMet,numMet);
        eliDiff.setZero(numMet,numMet);
        solDiff.setZero(numMet,numMet);

        for ( int m=0; m<numMet;++m)
        {
            for ( int om=m; om<numMet;++om)
            {
                sysDiff(om,m)= diff(sys[m],sys[om]);
                rhsDiff(om,m)= diff(rhs[m],rhs[om]);
                eliDiff(om,m)= diff(eli[m],eli[om]);
                solDiff(om,m)= diff(sol[m],sol[om]);
            }
        }
//            gsInfo
//                    <<"\nStiff differences::\n"<<sysDiff<<"\n"
//                   <<"\nRhs   differences::\n"<<rhsDiff<<"\n"
//                  <<"\nEli   differences::\n"<<eliDiff<<"\n"
//                 <<"\nSol  differences::\n"<<solDiff<<"\n"<<"\n";


        gsInfo<<"Table entries are in format: (errorL2,time)\n";
        for ( int m=0; m<numMet;++m)
            gsInfo<<metNames[m]<<std::setw(9)<<"   ";
        gsInfo<<"\n";

        for (int m=0;m<numMet;++m)
            gsInfo<<"("<<errL2(m,0)<<", "<<tim(m,0)<<"),    ";
        gsInfo<<"\n";

        real_t minErr=std::numeric_limits<real_t>::infinity();
        for (int m=0; m<numMet;++m)
            minErr= errL2(m,0)==NAN ? minErr : math::min(minErr,errL2(m,0));
        TEST(errL2(1,0)/minErr<=2);
    }
    return !passed;
}

void assembler1 (
        const ProblemData &data,
        gsSparseMatrix<>  &sys,
        gsMatrix<>        &rhs,
        gsMatrix<>        &eli,
        gsMatrix<>        &sol,
        real_t            &errorL2,
        real_t            &time
        )
{
    gsBiharmonicAssembler<real_t> BiharmonicAssembler(
                data.pde.domain(),
                data.basis,
                data.pde.boundaryConditions(),
                data.pde.boundaryConditionsSecond(),
                *data.pde.rhs(),
                data.Dstrategy,
                data.Istrategy
                );
    gsStopwatch clock;
    BiharmonicAssembler.assemble();
    time=clock.stop();

    sys = BiharmonicAssembler.matrix();
    rhs = BiharmonicAssembler.rhs();
    eli = BiharmonicAssembler.fixedDofs();
    sol = solve(sys,rhs);


    gsField<> solField = BiharmonicAssembler.constructSolution(sol);

    if (    !(data.Dstrategy==dirichlet::elimination)
            || ( !(data.pde.boundaryConditions().dirichletSides().size()>0) ) )
        eli.resize(0,rhs.cols());

    gsNormL2<real_t>     normL2(solField, data.solution);
    //gsSeminormH1<real_t> semiH1(solField, data.solution);
    errorL2 = normL2.compute();
    //errorH1 = semiH1.compute();
}


void assembler2 (
        const ProblemData &data,
        gsSparseMatrix<>  &sys,
        gsMatrix<>        &rhs,
        gsMatrix<>        &eli,
        gsMatrix<>        &sol,
        real_t            &errorL2,
        real_t            &time
        )
{
    bool eliminated = data.pde.boundaryConditions().dirichletSides().size()>0
            && data.Dstrategy==dirichlet::elimination;

    gsRecipeAssemblerBiharmonicSogn assembler(data.pde);
    assembler.setDirichletStrategy(data.Dstrategy);
    std::vector<gsPhysicalSpace*> phySpace;
    phySpace.push_back(new gsPhysicalSpaceScalar(data.basis,data.pde.domain(),INVERSE_COMPOSITION));
    assembler.setSpace(phySpace);

    gsStopwatch clock;
    assembler.assemble();
    rhs = assembler.getSystemRhs();

    if (eliminated)
    {
        eli  = solve(assembler.getEliminatedMatrix(),assembler.getEliminatedRhs());
        rhs -= assembler.getRhsModMatrix()*eli;
    }
    else
    {
        eli.resize(0,rhs.cols());
    }
    time=clock.stop();

    sys = assembler.getSystemMatrix();
    sol = solve(sys,rhs);

    // compute the error
    std::vector<gsMatrix<real_t> > coefs;
    for (size_t s=0; s< phySpace.size(); ++s)
        coefs.push_back(assembler.reconstructSolution(s,sol,eli));
    std::vector<gsFunction<real_t>*> svec(1,
    const_cast<gsFunctionExpr<real_t>*>(&data.solution));
    gsRecipeDistance dist(data.pde.domain(),phySpace,coefs,svec);
    gsMatrix<> norms(2,2);
    norms<<0,2,1,2;
    dist.setSeminorms(0,norms);
    dist.assemble();

    gsMatrix<> error=dist.getDistance(0).rowwise().sum();
    errorL2=math::sqrt(error(0,0));
    //errorH1=math::sqrt(error(1,0));
    freeAll(phySpace);
}





