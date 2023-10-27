/** \file
// This program tests that the gsRecipeAssembler and the gsPoissonAssembler
// assemble the same matrixes and right hand side.
//
// Test inclueds:
//
// Inhomogeneous source term (right hand side)
// Inhomogeneous Dirichlet boundary conditions by either Nitsche or Elimination
// Inhomogeneous Neumann boundary conditions
//
// Multiple right hand side
// Multipatch gluing by identification
//
//     x,y \in (0,1)
// Source function:
//
//     f(x,y)_x = ((pi*k0)^2 + (pi*k1)^2)*sin(pi*x*k0)*sin(pi*y*k1)
//     f(x,y)_y = ((pi*k2)^2 + (pi*k3)^2)*sin(pi*x*k2)*sin(pi*y*k3)
//
// Solution:
//
//     u(x,y)_x = sin(pi*x*k0)*sin(pi*y*k1)+pi/10
//     u(x,y)_y = sin(pi*x*k2)*sin(pi*y*k3)-pi/10
//
// This is based on PoissonVerificationTest by Jarle Sogn
//
// The two functions PoissonVolumeContributions and PoissonBoundaryContributions
// together are a Poisson assembler based on gsRecipeAssembler
//
**/

#include <iostream>
#include <iomanip>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsCore/gsField.h>
#include <gsCore/gsFunctionExpr.h>

#include <gsAssembler/gsAssemblerUtils.h>
#include <gsAssembler/gsStokesAssembler.h>
#include <gsAssembler/gsStokesAssembler2.h>

#include <gsRecipeAssembler/gsRecipeAssemblerStokes.h>
#include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>

#include <gsNurbs/gsNurbsCreator.h>


#include <gsSolver/gsStokesIterativeSolver.h>



#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n");\

bool passed=true;


using namespace gismo;


void addAllDirichletBoundaries(const gsMultiPatch<>& mp, gsFunction<real_t> *g,gsBoundaryConditions<real_t> & bcInfo)
{
    std::vector<patchSide> boundaries = mp.boundaries();
    for(unsigned i = 0;i<boundaries.size();++i)
        bcInfo.addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::dirichlet, g);
}


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
        return std::numeric_limits<real_t>::quiet_NaN();
}

class ProblemData
{
public:
    gsFunctionExpr<>      f;
    gsMultiPatch<>        *domain;

    gsFunctionExpr<>      u;
    gsFunctionExpr<>      p;
    gsMultiBasis<>         basis;
    gsStokesPde<real_t>    pde;

    dirichlet::strategy    Dstrategy;
    iFace::strategy        Istrategy;

    ProblemData()
    :
    f("+2*cos(x+y)+2*sin(x-y)","-2*cos(x+y)+2*sin(x-y)",2),
    domain(
        new gsMultiPatch<real_t>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus())),
    //      gsReadFile<>("square.xml" )),
    //  new gsMultiPatch<real_t>(*gsNurbsCreator<real_t>::BSplineSquareDeg(1))),
    //      gsReadFile<>("surfaces/multipatch_triangle2.xml")),

    u("cos(x+y)+sin(x-y)","-1-cos(x+y)+sin(x-y)",2),
    p("0",2),
    
    basis(*domain),
    
    pde(*domain,gsBoundaryConditions<>(),&f),
    Dstrategy(dirichlet::elimination),
    Istrategy(iFace::glue)
    {
        addAllDirichletBoundaries(*domain,&u,pde.boundaryConditions());
    }

    ~ProblemData()
    {
        delete domain;
    }

};


// assembling methods
typedef void (*assembler) (const ProblemData&, gsSparseMatrix<>&,gsMatrix<>&, gsMatrix<>&, gsMatrix<>&, real_t&, real_t&, real_t&);

void assembler1 (const ProblemData&, gsSparseMatrix<>&,gsMatrix<>&, gsMatrix<>&, gsMatrix<>&, real_t&, real_t&, real_t&);
void assembler2 (const ProblemData&, gsSparseMatrix<>&,gsMatrix<>&, gsMatrix<>&, gsMatrix<>&, real_t&, real_t&, real_t&);
void assembler3 (const ProblemData&, gsSparseMatrix<>&,gsMatrix<>&, gsMatrix<>&, gsMatrix<>&, real_t&, real_t&, real_t&);

assembler      methods[]={
    assembler1,
    assembler2,
    assembler3
};
const int numMet=sizeof(methods)/sizeof(assembler);

std::string    metNames  []={
    std::string("Assembler"),
    std::string("RecipeAssembler"),
    std::string("Assembler2")
};

// dirichlet strategy
dirichlet::strategy Dstrategy[]={dirichlet::elimination};
std::string         strName  []={
    std::string("elimination")
};
const int numStr=sizeof(Dstrategy)/sizeof(dirichlet::strategy);

const int numRefine = 0;     // Lowest number of refinement: numRefine + 1
const int maxIterations = 2; // Highes number of refinement: numRefine + maxIterations

int main(int argc, char *argv[])
{
    gsInfo.width(9);
    gsInfo<< std::setprecision(2);
    gsInfo<< std::scientific;

    ProblemData data;

    for (int i = 0; i < numRefine; ++i)
        data.basis.uniformRefine();

    for (int i = 0; i < maxIterations; ++i)
    {
        data.basis.uniformRefine();

        gsMatrix<real_t>        tim;
        tim.setConstant(numMet,numStr,std::numeric_limits<real_t>::quiet_NaN());
        gsMatrix<real_t>        errP;
        gsMatrix<real_t>        errV;

        errP.setConstant(numMet,numStr,std::numeric_limits<real_t>::quiet_NaN());
        errV.setConstant(numMet,numStr,std::numeric_limits<real_t>::quiet_NaN());

        gsMatrix<real_t>        sysDiff, rhsDiff, eliDiff,solDiff;

        for (int s=0; s<numStr;++s)
        {
            data.Dstrategy = Dstrategy[s];
            gsInfo<< ((Dstrategy[s]==dirichlet::elimination) ? "---------Elimination \n" : "------------Nietsche \n") ;

            gsSparseMatrix<real_t> sys[numMet];
            gsMatrix<real_t>       rhs[numMet];
            gsMatrix<real_t>       eli[numMet];
            gsMatrix<real_t>       sol[numMet];

            for ( int m=0; m<numMet;++m)
                methods[m](data,sys[m],rhs[m],eli[m],sol[m],errV(m,s),errP(m,s),tim(m,s));

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
            gsInfo
                    <<"\nStiff differences::\n"<<sysDiff<<"\n"
                   <<"\nRhs   differences::\n"<<rhsDiff<<"\n"
                  <<"\nEli   differences::\n"<<eliDiff<<"\n"
                 <<"\nSol  differences::\n"<<solDiff<<"\n"<<"\n";

        }

        gsInfo<<"Table entries are in format: (errorV,errorP,time)\n";
        gsInfo<<std::setw(18)<<" ";
        for ( int m=0; m<numMet;++m)
            gsInfo<<std::setw(32)<<metNames[m]<<"   ";
        gsInfo<<"\n";
        for ( int s=0; s<numStr;++s)
        {
            gsInfo<<std::setw(15)<<strName[s]<<"     ";
            for (int m=0;m<numMet;++m)
                gsInfo<<"("<<errV(m,s)<<", "<<errP(m,s)<<", "<<tim(m,s)<<"),    ";
            gsInfo<<"\n";
        }

        for (int s=0;s<numStr;++s)
        {
            real_t minErr=std::numeric_limits<real_t>::infinity();
            for (int m=0; m<numMet;++m)
                minErr= math::isnan(errV(m,s)) ? minErr : math::min(minErr,errV(m,s));
            TEST(errV(1,s)/minErr<=2);
        }
    }
    return !passed;
}

void assembler1 (
        const ProblemData  &data,
        gsSparseMatrix<>   &sys,
        gsMatrix<>         &rhs,
        gsMatrix<>         &eli,
        gsMatrix<>         &sol,
        real_t             &errV,
        real_t             &errP,
        real_t             &time
        )
{
    bool eliminated = (data.Dstrategy==dirichlet::elimination && data.pde.boundaryConditions().dirichletSides().size()>0);

    std::vector<std::vector<gsBasis<>*> > bases;

    std::vector<gsBasis<>*> tmpBasisVec;
    for (size_t i=0; i< data.pde.domain().nPatches();++i)
        tmpBasisVec.push_back(const_cast<gsBasis<>*>(&data.basis.basis(i)));

    std::vector<gsPhysicalSpace*> vec=constructTHSpaces(tmpBasisVec,data.pde.domain(),&bases);
    freeAll(vec);

    gsStokesAssembler<real_t> assembler(
                data.pde.domain(),
                data.pde.boundaryConditions(),
                bases,
                *data.pde.force()
                );
    assembler.setDirichletStrategy(data.Dstrategy);
    assembler.setInterfaceStrategy(data.Istrategy);

    gsStopwatch clock;
    assembler.initialize();
    assembler.assemble();
    time=clock.stop();

    sys = assembler.systemMatrix();
    rhs = assembler.systemRhs();
    eli = assembler.fixedValues()[0];

    sol = solve(sys,rhs);
    assembler.setSolutionVector(sol);
    assembler.reconstructSolution();

    if (!eliminated )
        eli.resize(0,rhs.cols());

    const gsField<> &field_V = assembler.solution(0);
    const gsField<> &field_P = assembler.solution(1);

    errV = igaFieldL2Distance( field_V, data.u);
    errP = igaFieldL2Distance( field_P, data.p);
}


void assembler2 (
        const ProblemData  &data,
        gsSparseMatrix<>   &sys,
        gsMatrix<>         &rhs,
        gsMatrix<>         &eli,
        gsMatrix<>         &sol,
        real_t             &errV,
        real_t             &errP,
        real_t             &time
        )
{
    bool eliminated = data.Dstrategy==dirichlet::elimination && data.pde.boundaryConditions().dirichletSides().size()>0;

    std::vector<std::vector<gsBasis<>*> > bases;
    gsRecipeAssemblerStokes              assembler(data.pde);

    std::vector<gsBasis<>*> tmpBasisVec;
    for (size_t i=0; i< data.pde.domain().nPatches();++i)
        tmpBasisVec.push_back(const_cast<gsBasis<>*>(&data.basis.basis(i)));
    std::vector<gsPhysicalSpace*> phySpace=constructTHSpaces(tmpBasisVec,data.pde.domain(),&bases);
    assembler.setSpace(phySpace);
    assembler.setZeroAverage(false);

    gsStopwatch clock;
    assembler.assemble();

    sys = assembler.getSystemMatrix();
    rhs = assembler.getSystemRhs();
    if (eliminated)
    {
        eli  = solve(assembler.getEliminatedMatrix(),assembler.getEliminatedRhs());
        rhs -= assembler.getRhsModMatrix()*eli;
    }
    else
        eli.resize(0,rhs.cols());
    time=clock.stop();

    sol = solve(sys,rhs);
    if (assembler.getZeroAverage())
    {
        sys.conservativeResize(sys.rows()-1,sys.cols()-1);
        rhs.resize(rhs.rows()-1,rhs.cols());
        sol.resize(sol.rows()-1,sol.cols());
    }

    // compute the error
    std::vector<gsMatrix<real_t> > coefs;
    for (size_t s=0; s< phySpace.size(); ++s)
        coefs.push_back(assembler.reconstructSolution(s,sol,eli));
    std::vector<gsFunction<real_t>*> svec(2);
    svec[0] = const_cast<gsFunctionExpr<real_t>*>(&data.u);
    svec[1] = const_cast<gsFunctionExpr<real_t>*>(&data.p);
    gsRecipeDistance dist(data.pde.domain(),phySpace,coefs,svec);
    dist.assemble();

    errV=math::sqrt(dist.getDistance(0).sum());
    errP=math::sqrt(dist.getDistance(1).sum());

    for (size_t s=0; s<bases.size();++s)
        freeAll(bases[s]);
    freeAll(phySpace);
}




void assembler3 (
        const ProblemData  &data,
        gsSparseMatrix<>   &sys,
        gsMatrix<>         &rhs,
        gsMatrix<>         &eli,
        gsMatrix<>         &sol,
        real_t             &errV,
        real_t             &errP,
        real_t             &time
        )
{
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back( gsMultiBasis<>(data.basis) );//Basis for velocity
    discreteBases.push_back( gsMultiBasis<>(data.basis) );//Basis for pressure
    //p-refine the velocity space (Taylor-Hood type)
    discreteBases[0].degreeElevate(1);

    gsStokesAssembler2<real_t> stokes(data.pde.domain(), discreteBases, data.pde.boundaryConditions(), *data.pde.rhs(), data.Dstrategy);
    stokes.setViscosity(1);

    gsStopwatch clock;
    stokes.assemble();
    time = clock.stop();

    sys = stokes.matrix().selfadjointView<Lower>();
    rhs = stokes.rhs();
    eli = stokes.dirValues();
    eli.resize(eli.size(),1);

    sol = solve(sys,rhs);

    gsField<> field_V  = stokes.constructSolution(sol, 0);
    gsField<> field_P  = stokes.constructSolution(sol, 1);

    errV = igaFieldL2Distance( field_V, data.u);
    errP = igaFieldL2Distance( field_P, data.p);
}
