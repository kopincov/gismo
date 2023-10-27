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

#include <gsUtils/gsStopwatch.h>

#include <gsCore/gsField.h>
#include <gsCore/gsFunctionExpr.h>

#include <gsAssembler/gsAssemblerUtils.h>
#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsSeminormH1.h>

#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>
#include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>

#include <gsNurbs/gsNurbsCreator.h>

#include <gsCore/gsTransformedFuncSet.h>


using namespace gismo;

#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n");\

bool passed=true;

#ifdef NAN
#undef NAN
#endif
static real_t NAN(std::numeric_limits<real_t>::quiet_NaN());


void addAllDirichletBoundaries(const gsMultiPatch<>& mp, gsFunction<real_t> *g,gsBoundaryConditions<real_t> & bcInfo)
{
    std::vector<patchSide> boundaries = mp.boundaries();
    for(unsigned i = 0;i<boundaries.size();++i)
        bcInfo.addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::dirichlet, g);
}

gsMatrix<> solve(const gsSparseMatrix<> &sys, const gsMatrix<> &rhs)
{
    gsSparseSolver<>::CGDiagonal solver;
    return solver.compute(sys).solve(rhs);
}

template <typename MatrixT>
real_t diff(const MatrixT &a,const MatrixT &b)
{
    if (
            a.rows()==b.rows() &&
            a.cols()==b.cols() &&
            a.cols()>0 &&
            a.rows()>0 )
        return (a-b).norm();
    else
        return NAN;
}

class ProblemData
{
protected:
    gsMultiPatch<>         *domain;
    gsFunctionExpr<real_t> valu;
    gsFunctionExpr<real_t> f;
    gsFunctionExpr<real_t> gradu;
public:
    gsFunctionWithDerivatives<real_t> u;
    gsMultiBasis<>         basis;
    gsPoissonPde<real_t>   pde;

    dirichlet::strategy    Dstrategy;
    iFace::strategy        Istrategy;
public:
    ProblemData(std::string geoname)
        :
          domain(
              ( (gsMultiPatch<>::uPtr)gsReadFile<real_t>(geoname) ).release()
          ),
          valu(
              "sin(pi*x*1)*sin(pi*y*2)+pi/10",
              /*"sin(pi*x*3)*sin(pi*y*4)-pi/10",*/2),
          f(
              "((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",
              /*"((pi*3)^2 + (pi*4)^2)*sin(pi*x*3)*sin(pi*y*4)",*/2),
          gradu(
              "pi*cos(pi*x*1)*sin(pi*y*2)",
              "2*pi*sin(pi*x*1)*cos(pi*y*2)",
              /*"3*pi*cos(pi*x*3)*sin(pi*y*4)",
              "4*pi*sin(pi*x*3)*cos(pi*y*4)",*/2
              ),
          u(valu,gradu),
          basis(*domain),
          pde(
              *domain,
              gsBoundaryConditions<>(),
              f),
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

void assembleVisitor (const ProblemData&, gsSparseMatrix<>&, gsMatrix<>&, gsMatrix<>&, gsMatrix<>&, real_t &, real_t&, real_t &time);
void assembleRecipe  (const ProblemData&, gsSparseMatrix<>&, gsMatrix<>&, gsMatrix<>&, gsMatrix<>&, real_t &, real_t&, real_t &time);


assembler methods[]={
    assembleRecipe,
    assembleVisitor,
};
const int numMet=sizeof(methods)/sizeof(assembler);

std::string    metNames  []={
    std::string("Recipe"),
    std::string("Visitor2"),
};

// dirichlet strategy
dirichlet::strategy Dstrategy[]={
    dirichlet::elimination,
    dirichlet::nitsche
    };
std::string         strName  []={
    std::string("elimination"),
    std::string("nitsche")
};
const int numStr=sizeof(Dstrategy)/sizeof(dirichlet::strategy);



int main(int argc, char *argv[])
{


    index_t numRefine     = 3; // Lowest number of refinement: numRefine + 1
    index_t maxIterations = 1; // Highes number of refinement: numRefine + maxIterations
    std::string geoFileName="yeti_mp.xml";


    gsInfo.width(9);
    gsInfo<< std::setprecision(2);
    gsInfo<< std::scientific;

    gsCmdLine cmd("Compares two different poisson assemblers");
    cmd.addInt("r","refine","starting refinement level", numRefine);
    cmd.addInt("i","iteration","number of refinement iterations", maxIterations);
    cmd.addString("g","geometry", "file containing the domain", geoFileName);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    ProblemData data(geoFileName);

    data.basis.degreeElevate();
    for (int i = 0; i < numRefine; ++i)
        data.basis.uniformRefine();

    for (int i = 0; i < maxIterations; ++i)
    {
        data.basis.uniformRefine();

        gsMatrix<real_t>        tim;
        tim.setConstant(numMet,numStr,NAN);
        gsMatrix<real_t>        errL2;
        errL2.setConstant(numMet,numStr,NAN);
        gsMatrix<real_t>        errH1;
        errH1.setConstant(numMet,numStr,NAN);

        gsMatrix<real_t>        sysDiff, rhsDiff, eliDiff,solDiff;

        std::vector<gsSparseMatrix<real_t> >  sys(numMet,gsSparseMatrix<>(0,0));
        std::vector<gsMatrix<real_t> >        rhs(numMet,gsMatrix<>(0,0));
        std::vector<gsMatrix<real_t> >        eli(numMet,gsMatrix<>(0,0));
        std::vector<gsMatrix<real_t> >        sol(numMet,gsMatrix<>(0,0));

        for (int s=0; s<numStr;++s)
        {
            data.Dstrategy = Dstrategy[s];

            for ( int m=0; m<numMet;++m)
                methods[m](data,sys[m],rhs[m],eli[m],sol[m],errL2(m,s),errH1(m,s),tim(m,s));

            gsInfo<< ((Dstrategy[s]==dirichlet::elimination) ?
                          "---------Elimination: dofs=" :
                          "------------Nietsche: dofs=")
                  << sys[0].rows() <<std::endl;


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

        gsInfo<<"Table entries are in format: (errorL2,errorH1,time)\n";
        gsInfo<<std::setw(18)<<" ";

        for ( int m=0; m<numMet;++m)
            gsInfo<<std::setw(32)<<metNames[m]<<"   ";
        gsInfo<<"\n";
        for ( int s=0; s<numStr;++s)
        {
            gsInfo<<std::setw(15)<<strName[s]<<"    ";
            for (int m=0;m<numMet;++m)
                gsInfo<<"("<<errL2(m,s)<<", "<<errH1(m,s)<<", "<<tim(m,s)<<"),    ";
            gsInfo<<"\n";
        }

        for (int s=0;s<numStr;++s)
        {
            real_t minErr=std::numeric_limits<real_t>::infinity();
            for (int m=0; m<numMet;++m)
                minErr= errL2(m,s)==NAN ? minErr : math::min(minErr,errL2(m,s));
            TEST(errL2(1,s)/minErr<=2);
        }
    }
    return !passed;
}


void assembleVisitor (
        const ProblemData &data,
        gsSparseMatrix<>  &sys,
        gsMatrix<>        &rhs,
        gsMatrix<>        &eli,
        gsMatrix<>        &sol,
        real_t            &errorL2,
        real_t            &errorH1,
        real_t            &time
        )
{
    gsPoissonAssembler<real_t> PoissonAssembler(
                data.pde.domain(),
                data.basis,
                data.pde.boundaryConditions(),
                *data.pde.rhs(),
                data.Dstrategy,
                data.Istrategy
                );
    gsStopwatch clock;
    PoissonAssembler.assemble();
    time=clock.stop();

    sys = PoissonAssembler.matrix();
    rhs = PoissonAssembler.rhs();
    eli = PoissonAssembler.dirValues();
    sol = solve(sys,rhs);


    gsField<> solField = PoissonAssembler.constructSolution(sol);

    if (
            (!(data.Dstrategy==dirichlet::elimination))
            || ( !(data.pde.boundaryConditions().dirichletSides().size()>0) )
            )
        eli.resize(0,rhs.cols());

    gsNormL2<real_t>     normL2(solField, data.u);

    gsSeminormH1<real_t> semiH1(solField, data.u);
    errorL2 = normL2.compute();
    errorH1 = semiH1.compute();
}


class gsRecipeAssemblerPoisson2 : public gsRecipeAssemblerPoisson
{
public:
    gsRecipeAssemblerPoisson2(const gsPoissonPde<real_t> &pde)
        : gsRecipeAssemblerPoisson(pde)
    {
    }
gsRecipeAssembler::SpaceList getPatchSpaces (index_t patch )
{
    SpaceList result;
    result.push_back(m_space[0]->getPatchSpace(patch));
    static_cast<gsTransformedFuncSet<real_t> *>(result.back().get())->activeShift+=m_shiftsSource[0];
    result.push_back(gsRecipe<real_t>::spacePtr(new gsRestrictTFS(*m_pde.rhs())));
    result.push_back(gsRecipe<real_t>::spacePtr(new gsConstantFunction<real_t>(1.0,m_domain.dim())));
    return result;
}

gsRecipe<real_t>    getPatchRecipe     (index_t patch)
{
    gsRecipe<real_t>            result;
    gsRecipeIngredient<real_t>  ingr;

    ingr.setOperator(new gsGradGradOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(0);
    ingr.setRule(getSysWriter());
    result.add(ingr);

    ingr.setOperator(new gsL2ScalarOp<real_t>());
    ingr.setTestSpace(0);
    ingr.setUnknownSpace(1);
    ingr.setRule(getRhsWriter());
    result.add(ingr);

    if (m_zeroAverage)
    {
        ingr.setOperator(new gsL2ScalarOp<real_t>());
        ingr.setTestSpace(0);
        ingr.setUnknownSpace(2);
        ingr.setRule(new gsL2GMapperRhs<MulW>(*m_map,
                                              MulW(SysW(getFreeLimit(),getFreeLimit(),m_sysM,m_rhsMod),0,getFreeLimit()))
                     );
        result.add(ingr);
    }
    return result;
}
};

void assembleRecipe (
        const ProblemData &data,
        gsSparseMatrix<>  &sys,
        gsMatrix<>        &rhs,
        gsMatrix<>        &eli,
        gsMatrix<>        &sol,
        real_t            &errorL2,
        real_t            &errorH1,
        real_t            &time
        )
{
    bool eliminated = data.pde.boundaryConditions().dirichletSides().size()>0
            && data.Dstrategy==dirichlet::elimination;

    gsRecipeAssemblerPoisson2 assembler(data.pde);
    assembler.setDirichletStrategy(data.Dstrategy);
    std::vector<gsPhysicalSpace*> phySpace;
    phySpace.push_back(
                new gsPhysicalSpaceScalar(data.basis,data.pde.domain(),INVERSE_COMPOSITION)
                );
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
    std::vector<gsFunction<>*> ss(1,const_cast<gsFunctionWithDerivatives<>*>(&data.u));
    gsRecipeDistance dist(data.pde.domain(),phySpace,coefs,ss);
    gsMatrix<> norms(2,2);
    norms<<0,2,1,2;
    dist.setSeminorms(0,norms);
    dist.assemble();

    gsMatrix<> error=dist.getDistance(0).rowwise().sum();
    errorL2=math::sqrt(error(0,0));
    errorH1=math::sqrt(error(1,0));
    freeAll(phySpace);
}


