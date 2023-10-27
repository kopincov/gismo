/** @file gsRecipeAssemblerShowcase.cpp

    @brief Showcase for meny PDEs with the recipe assembler (Currantly under construction).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan J. Sogn
*/


#include <gismo.h>
#include <gismo_dev.h>

#include <gsCore/gsTransformedFuncSet.h>

#include <gsRecipeAssembler/gsRecipe.h>
#include <gsRecipeAssembler/gsPDEOperators.h>
#include <gsRecipeAssembler/gsLocalToGlobal.h>

#include <gsRemappedBasis/gsRemappedBasis.h>

using namespace gismo;

namespace Utils
{
void addAllDirichletBoundaries(const gsMultiPatch<>& mp, gsFunction<real_t> *g,gsBoundaryConditions<real_t> & bcInfo);
void computeError( gsFunctionSet<> &exact, gsFunctionSet<> &approx,
                   const gsMultiPatch<> &domain, const gsElementList<> &elements, const gsVector<index_t> &nGaussP,
                   gsMatrix<> &l2Err, gsMatrix<> &h1Err);


template<typename T, int Major, typename IndexM >
void forceDirichletConditions ( gsSparseMatrix<T,Major>      &M,
                                gsMatrix<T>                  &Rhs,
                                const IndexM                 &dirDofs,
                                const gsMatrix<T>            &dirValues);




void                  computeDirichletValues  (const gsMultiBasis<> &orig,const gsRecipe<>::spacePtr space, const gsMultiPatch<>& domain, const gsRecipe<>::spacePtr dirichletData, std::vector<index_t> &boundaryDofs, gsMatrix<> &boundaryVal);
gsElementList<>       getPatchElements        (const gsMultiBasis<> &basis);
gsElementList<>       getBoundaryElements     (const gsMultiBasis<> &basis);
std::vector<index_t>  getBoundaryDofs         (const gsMultiBasis<> &basis, const gsWeightMapper<real_t>* remap );
void                  getSuggestedIntegration (const gsMultiBasis<> &basis, gsElementList<> &elements, gsVector<index_t> &nGaussP);
}


namespace Poisson {
void solve             (const gsMultiPatch<> &domain, const gsMultiBasis<> &basis);
void solveElimination  (const gsMultiPatch<> &domain, const gsMultiBasis<> &basis);
void solveNitsche      (const gsMultiPatch<> &domain, const gsMultiBasis<> &basis);
}

int main(int argc, char *argv[])
{
    index_t numRefine = 0;
    index_t degree    = 0;

    std::cout.precision(2);
    std::cout<< std::scientific;

    std::string geoFileName("square.xml");

    gsCmdLine cmd("Many examples of assemblers for specific equations");
    cmd.addInt   ("r","refine",   "number of dyadic subdivisions of the geometry basis to obtain the space", numRefine);
    cmd.addInt   ("d","degree",   "minimum degree of the discrete basis",                                    degree);
    cmd.addString("g","geometry", "file describing the domain",                                              geoFileName);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsMultiPatch<>::uPtr domain = gsReadFile<>(geoFileName);
    if (!domain)
    {
        GISMO_ERROR("Error reading "+geoFileName+" or it does not contain a geometry\n");
    }

    gsMultiBasis<> startBasis(*domain);
    if (startBasis.totalSize()==4 && numRefine == 0 && degree == 0)
        gsInfo << "Changed required refinement because otherwise there are no internal dofs. Not it is "<<++numRefine<<std::endl;

    while (startBasis.degree()<degree)
        startBasis.degreeElevate();
    for   (int r=0; r<numRefine;++r)
        startBasis.uniformRefine();



    Poisson::solve            (*domain,startBasis);
    Poisson::solveElimination (*domain,startBasis);
    Poisson::solveNitsche     (*domain,startBasis);

}


namespace Poisson
{
// The Poisson problem is: find u ��� H��:
//
// -��u = f on ��
//   u = g on �����
//
// where f ��� H-��, g ��� H��.
//
// Its weak formulation is:  ��� w ��� H��
//
// ������u.���w = ���f w
//
// where the result of two bilinear operators is compared.
// So there are two operators:
// -gsGradGradOp -> ������u.���w
// -gsL2ScalarOp -> ���f w
// and three spaces involved plus the boundary data
enum
{
    trial     = 0, // the space of u
    test      = 0, // the space of w (we set them to be the same in this case)
    force     = 1, // the space {f}
    dirichlet = 2, // the space {g}
};


void solve (const gsMultiPatch<> &domain,const gsMultiBasis<> &basis)
{
    // data
    gsFunctionExpr<> f("(-2*cos(x^2) + (1 + 4*x^2)*sin(x^2) )*sin(y)",domain.geoDim());
    gsFunctionExpr<> g("sin(x^2)sin(y)",domain.geoDim()); // dirichlet data and exact solution

    // prepare the discretization space (this produces a unique enumerarion of the dofs)
    gsRemappedBasis *paraSpace = gsRemappedBasis::makeMultiPatch(basis);

    // construct the list of spaces mapped to the physical domain
    gsRecipe<>::spaceList spaces;
    spaces.push_back( gsRecipe<>::spacePtr(new gsGradConformingTFS (*paraSpace)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsRestrictTFS (f)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsRestrictTFS (g)));

    // construct integration rule given by a per element quadrature and an element list
    gsVector<index_t>   nGaussP;;
    gsElementList<> elements;
    Utils::getSuggestedIntegration(basis,elements,nGaussP);

    // prepare the destination matrices
    gsSparseMatrix<> A(spaces[test]->size(),spaces[trial]->size());
    A.reserve ( gsVector<index_t>::Constant(spaces[test]->size(),1,(nGaussP*2).prod()) ); //expected number of non zero per rows, depends on spline degree
    gsMatrix<>       y;
    y.setZero (spaces[test]->size(),spaces[force]->size());

    // construct the proper recipe and assemble
    gsRecipe<> poisson;
    poisson.add(gsRecipe<>::opHandle(new gsGradGradOp<>()),getWriter(A) ,trial, test);
    poisson.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()),getWriter(y) ,force, test);
    poisson.assemble(domain,spaces,nGaussP,elements);

    // modify the matrix to force Dirichlet conditions on all boundaries
    std::vector<index_t> fixedDofs=Utils::getBoundaryDofs(basis,&paraSpace->getMapper());
    gsMatrix<>           fixedVals;
    Utils::computeDirichletValues(basis,spaces[trial],domain,spaces[dirichlet],fixedDofs,fixedVals);
    Utils::forceDirichletConditions(A,y,fixedDofs,fixedVals);

    // solve the system
    A.makeCompressed();
    gsMatrix<>       solCoefs = gsSparseSolver<>::LU(A).solve(y);

    // get the solution as function
    gsFunctionSet<> *solParam = paraSpace->makeReduced(solCoefs,false,solCoefs.cols());
    gsFunctionSet<> *solPhys  = new gsGradConformingTFS(*solParam);

    // compare with exact solution
    gsMatrix<> l2Err;
    gsMatrix<> h1Err;
    Utils::computeError(*spaces[dirichlet],*solPhys, domain,elements,nGaussP.array()+2,l2Err,h1Err);
    gsInfo<<"Solved Poisson with matrix modification, the error is (L2,H1): ("<<l2Err<<", "<<h1Err<<")"<<std::endl;

    elements.free();
    delete paraSpace;
    delete solParam;
    delete solPhys;
}


// When the Nitsche method is employed to solve a Poisson problem
// the weak formulation is modified to include boundary integrals
// that penalize the norm of (u-g) on �����.
//
// This method requires that the discretizations of the trial space
// and of the test space are the same. (No Petrov-Galerkin allowed)
//
// Its weak formulation is:  ��� w ��� H��
//
// ������u.���w + s ���_��� u w - ���_��� n.���u w - ���_��� u n.���w = s ���_��� g w +  ���f w - ���_��� g n.���w
//
// where s is a mesh dependent stabilization parameter,
//  ���_��� are boundary integrals
//  n is the outer normal on the boundary
// where the result of two bilinear operators is compared.
// So there are the standard two operators for Poisson on the patches:
// -gsGradGradOp -> ������u.���w
// -gsL2ScalarOp -> ���f w
// plus other integral on the Dirichlet boundary with the following
// operators:
// -gsL2ScalarOp -> s ���_��� u w
// -gsL2ScalarOp -> s ���_��� g w
// -gsBoundaryNormalDerValueOp -> - ���_��� n.���u w - ���_��� u n.���w
// we take care of the coeffients (s and -1) and of the symmetry of
// - ���_��� n.���u w - ���_��� u n.���w in the writing rule.

void solveNitsche      (const gsMultiPatch<> &domain, const gsMultiBasis<> &basis)
{
    // data, these are parametric functions
    gsFunctionExpr<> f("(-2*cos(x^2) + (1 + 4*x^2)*sin(x^2) )*sin(y)",domain.geoDim());
    gsFunctionExpr<> g("sin(x^2)sin(y)",domain.geoDim()); // dirichlet data and exact solution

    // prepare the discretization space (this produces a unique enumerarion of the dofs)
    gsRemappedBasis *paraSpace = gsRemappedBasis::makeMultiPatch(basis);

    // construct the list of spaces mapped to the physical domain
    gsRecipe<>::spaceList spaces;
    spaces.push_back( gsRecipe<>::spacePtr(new gsGradConformingTFS (*paraSpace)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsRestrictTFS (f)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsRestrictTFS (g)));

    // construct integration rule given by a per element quadrature and an element list
    gsVector<index_t>   nGaussP;;
    gsElementList<> elements;
    Utils::getSuggestedIntegration(basis,elements,nGaussP);

    // prepare the destination matrices
    gsSparseMatrix<> A(spaces[test]->size(),spaces[trial]->size());
    gsMatrix<>       y;
    y.setZero (spaces[test]->size(),spaces[force]->size());
    A.reserve ( gsVector<index_t>::Constant(spaces[test]->size(),1,(nGaussP*2).prod()) ); //expected number of non zero per rows, depends on spline degree

    // construct the proper recipe and assemble
    gsRecipe<> poisson;
    poisson.add(gsRecipe<>::opHandle(new gsGradGradOp<>()),getWriter(A) ,trial, test);
    poisson.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()),getWriter(y) ,force, test);
    poisson.assemble(domain,spaces,nGaussP,elements);

    // integrate on the Dirichlet boundary
    real_t nitscheParam = 2*(spaces[trial]->domainDim()+nGaussP.maxCoeff())* nGaussP.maxCoeff()*math::pow(spaces[trial]->size(),1.0/ spaces[trial]->domainDim());
    real_t opCoef=-1;
    gsRecipe<> nitsche;
    nitsche.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()), getWriterWithCoef(A,nitscheParam) ,trial,     test);
    nitsche.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()), getWriterWithCoef(y,nitscheParam) ,dirichlet, test);
    nitsche.add(gsRecipe<>::opHandle(new gsBoundaryNormalDerValueOp<>()), getSymmetricWriterWithCoef(A,opCoef) , trial, test);
    nitsche.add(gsRecipe<>::opHandle(new gsBoundaryNormalDerValueOp<>()), getTransposeWriterWithCoef(y,opCoef) , test, dirichlet);

    gsElementList<> b_elements=Utils::getBoundaryElements(basis);
    nitsche.assemble(domain,spaces,nGaussP,b_elements);

    // solve the system
    A.makeCompressed();
    gsMatrix<>       solCoefs = gsSparseSolver<>::LU(A).solve(y);

    // get the solution as function
    gsFunctionSet<> *solParam = paraSpace->makeReduced(solCoefs,false,solCoefs.cols());
    gsFunctionSet<> *solPhys  = new gsGradConformingTFS(*solParam);

    // compare with exact solution
    gsMatrix<> l2Err;
    gsMatrix<> h1Err;
    Utils::computeError(*spaces[dirichlet],*solPhys, domain,elements,nGaussP,l2Err,h1Err);
    gsInfo<<"Solved Poisson with the Nitsche method,  the error is (L2,H1): ("<<l2Err<<", "<<h1Err<<")"<<std::endl;

    elements.free();
    b_elements.free();

    delete paraSpace;
    delete solParam;
    delete solPhys;
}


// In order to remove the degrees of freedom corresponding to eliminated dofs from
// the system matrix we change the order of the dofs of the space with a
// permutation matrix. Then we can simply split the output into two: over the number of
// eliminated dofs and under the number and write them to different matrices.

void solveElimination (const gsMultiPatch<> &domain,const gsMultiBasis<> &basis)
{
    // data
    gsFunctionExpr<> f("(-2*cos(x^2) + (1 + 4*x^2)*sin(x^2) )*sin(y)",domain.geoDim());
    gsFunctionExpr<> g("sin(x^2)sin(y)",domain.geoDim()); // dirichlet data and exact solution

    // prepare the discretization space (this produces a unique enumerarion of the dofs)
    gsRemappedBasis *paraSpace = gsRemappedBasis::makeMultiPatch(basis);

    // construct the list of spaces mapped to the physical domain
    gsRecipe<>::spaceList spaces;
    spaces.push_back( gsRecipe<>::spacePtr(new gsGradConformingTFS (*paraSpace)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsRestrictTFS (f)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsRestrictTFS (g)));

    // now we collecte the eliminated dofs and we reorder the indices of all dofs
    std::vector<index_t> fixedDofs=Utils::getBoundaryDofs(basis,&paraSpace->getMapper());
    gsMatrix<>           fixedVals;
    Utils::computeDirichletValues(basis,spaces[trial],domain,spaces[dirichlet],fixedDofs,fixedVals);
    const index_t totaSize = paraSpace->size();
    const index_t elimSize = fixedDofs.size();
    const index_t freeSize = totaSize-elimSize;

    gsVector<index_t> newOrder(totaSize);
    std::vector<index_t>::const_iterator it  = fixedDofs.begin();
    std::vector<index_t>::const_iterator end = fixedDofs.end();
    index_t drift=0;
    for (index_t dof=0;dof<totaSize;++dof)
    {
        newOrder(dof-drift)=dof;
        if (it!=end && dof==*it)
        {
            ++it;
            ++drift;
        }
    }
    it=fixedDofs.begin();
    for (index_t dof=freeSize;dof<totaSize && it!= end; ++it,++dof)
        newOrder(dof)=*it;

    paraSpace->reduce(Eigen::PermutationMatrix<-1,-1,index_t>(newOrder)); // the eliminate dofs are the end

    // construct integration rule given by a per element quadrature and an element list
    gsVector<index_t>   nGaussP;;
    gsElementList<> elements;
    Utils::getSuggestedIntegration(basis,elements,nGaussP);

    // the matrices
    gsSparseMatrix<> A(freeSize,freeSize);
    gsSparseMatrix<> R(freeSize,elimSize); // this is the linear operator from the boundary dofs to the
                        // modification of the right hand side
    A.reserve ( gsVector<index_t>::Constant(freeSize,1,(nGaussP*2).prod()) ); //expected number of non zero per rows
    R.reserve ( gsVector<index_t>::Constant(elimSize,1,(nGaussP*2).prod()) ); //expected number of non zero per rows
    gsMatrix<>       y;
    y.setZero (freeSize,spaces[force]->size());

    // construct the proper recipe and assemble
    gsRecipe<> poisson;
    poisson.add(gsRecipe<>::opHandle(new gsGradGradOp<>()),getSplitWriter(A,R,  freeSize) ,trial, test);
    poisson.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()),getSplitWriterRhs(y, freeSize) ,force, test);
    poisson.assemble(domain,spaces,nGaussP,elements);

    A.makeCompressed();
    R.makeCompressed();

    // modify the right hand side according to the Dirichlet conditions
    y-=R*fixedVals; // y can be saved to solve for many different Dirichlet conditions without re-assembling

    // solve the system
    gsMatrix<> solCoefs(totaSize, f.size());
    solCoefs.topRows   (freeSize) = gsSparseSolver<>::LU(A).solve(y);
    solCoefs.bottomRows(elimSize) = fixedVals;

    // get the solution as function
    gsFunctionSet<> *solParam = paraSpace->makeReduced(solCoefs,false,solCoefs.cols());
    gsFunctionSet<> *solPhys  = new gsGradConformingTFS(*solParam);

    // compare with exact solution
    gsMatrix<> l2Err;
    gsMatrix<> h1Err;
    Utils::computeError(*spaces[dirichlet],*solPhys, domain,elements,nGaussP,l2Err,h1Err);
    gsInfo<<"Solved Poisson problem with elimination, the error is (L2,H1): ("<<l2Err<<", "<<h1Err<<")"<<std::endl;

    elements.free();
    delete paraSpace;
    delete solParam;
    delete solPhys;
}

}
namespace KKTPoisson
{
// The optimisation problem: find u ��� H�� and f ��� L�� s. t.:
//
// min(u,f) = {1/2||u-d||�� + a/2 ||f||�� }
//
// subject to
// ��u = f on ��
//  u = g on �����
//
// where d, g ��� L��.
//
// Its KKT formulation is:  ��� v ��� H��_0
//
//  ���u u + ������u.���v = ���d u
// a���f f - ���f v   = 0
//  ������v.���u-���v f   = ���g v
//
// This is a three by three block system.
// There are 2 operators:
// -gsGradGradOp -> ������u.���v
// -gsL2ScalarOp -> ���f v
// and four spaces involved plus the boundary data
enum
{
    state     = 0, // the space of u
    control   = 1, // the space of f
    lagrange  = 2, // the space of w
    desired   = 3, // the space {d}
    dirichlet = 4, // the space {g}
};
/*
void solve (const gsMultiPatch<> &domain,const gsMultiBasis<> &basis)
{
    // data
    //gsFunctionExpr<> f("(-2*cos(x^2) + (1 + 4*x^2)*sin(x^2) )*sin(y)",domain.geoDim());
    gsFunctionExpr<> g("sin(x^2)sin(y)",domain.geoDim()); // dirichlet data and desired solution
    real_t alpha= 0.1; //Tikhonov regularization parameter

    // the matrices
    gsSparseMatrix<> A;
    gsMatrix<>       y;

    // construct the proper recipe
    gsRecipe<> KKTpoisson;

    real_t opCoef=-1;
    KKTpoisson.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()),getWriter(A) ,state, state);
    KKTpoisson.add(gsRecipe<>::opHandle(new gsGradGradOp<>()),getSymmetricWriter(A) ,state, lagrange);
    KKTpoisson.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()),getWriter(y) ,desired, state);

    KKTpoisson.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()),getWriterWithCoef(A,alpha) ,control, control);
    KKTpoisson.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()),getSymmetricWriterWithCoef(A,opCoef) ,control, control);
    KKTpoisson.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()),getWriter(y) ,dirichlet, state);

    // prepare the discretization space (this produces a unique enumerarion of the dofs)
    gsRemappedBasis *paraSpace = gsRemappedBasis::makeMultiPatch(basis);

    // construct the list of spaces mapped to the physical domain
    gsRecipe<>::spaceList spaces;
    spaces.push_back( gsRecipe<>::spacePtr(new gsGradConformingTFS (*paraSpace)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsGradConformingTFS (*paraSpace)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsGradConformingTFS (*paraSpace)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsRestrictTFS (g)));
    spaces.push_back( gsRecipe<>::spacePtr(new gsRestrictTFS (g)));

    // construct integration rule given by a per element quadrature and an element list
    gsVector<index_t>   nGaussP;;
    gsElementList<> elements;
    Utils::getSuggestedIntegration(basis,elements,nGaussP);

    // prepare the destination matrices
    index_t totalSize = spaces[state]->size() + spaces[control]->size() + spaces[lagrange]->size();
    A.resize  (totalSize, totalSize);
    A.reserve ( gsVector<index_t>::Constant(totalSize,(nGaussP*2).prod()) ); //expected number of non zero per rows, depends on spline degree
    y.setZero (totalSize, spaces[g]->size());
    // assemble
    KKTpoisson.assemble(domain,spaces,nGaussP,elements);

    // modify the matrix to force Dirichlet conditions on all boundaries
    std::vector<index_t> fixedDofs=Utils::getBoundaryDofs(basis,&paraSpace->getMapper());
    gsMatrix<>           fixedVals;
    Utils::computeDirichletValues(basis,spaces[trial],domain,spaces[dirichlet],fixedDofs,fixedVals);
    Utils::forceDirichletConditions(A,y,fixedDofs,fixedVals);

    // solve the system
    A.makeCompressed();
    gsMatrix<>       solCoefs = gsSparseSolver<>::LU(A).solve(y);

    // get the solution as function
    gsFunctionSet<> *solParam = paraSpace->makeReduced(solCoefs,false,solCoefs.cols());
    gsFunctionSet<> *solPhys  = new gsGradConformingTFS(*solParam);

    // compare with exact solution
    gsMatrix<> l2Err;
    gsMatrix<> h1Err;
    Utils::computeError(*spaces[dirichlet],*solPhys, domain,elements,nGaussP.array()+2,l2Err,h1Err);
    gsInfo<<"Solved a Poisson problem, the error is (L2,H1): ("<<l2Err<<", "<<h1Err<<")"<<std::endl;

    elements.free();
    delete paraSpace;
    delete solParam;
    delete solPhys;
}*/

}


namespace Utils
{

void computeError(gsFunctionSet<> &exact, gsFunctionSet<> &approx,
                  const gsMultiPatch<> &domain, const gsElementList<> &elements, const gsVector<index_t> &nGaussP,
                  gsMatrix<> &l2Err, gsMatrix<> &h1Err)
{
    gsMatrix<> mass;
    gsMatrix<> stiff;

    mass.setZero (exact.size() , exact.size());
    stiff.setZero(exact.size() , exact.size());

    gsRecipe<>::spaceList errSpace;
    errSpace.push_back( gsRecipe<>::spacePtr(new gsFunctionSetDiff<>(exact ,approx)));
    gsRecipe<> errorsRecipe;
    errorsRecipe.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()),getWriter(mass) ,0, 0);
    errorsRecipe.add(gsRecipe<>::opHandle(new gsGradGradOp<>()),getWriter(stiff) ,0, 0);
    errorsRecipe.assemble(domain,errSpace,nGaussP,elements);

    l2Err = mass.diagonal().transpose().array().sqrt();
    h1Err = (stiff.diagonal().transpose()+mass.diagonal().transpose()).array().sqrt();
}


void  getSuggestedIntegration (const gsMultiBasis<> &basis, gsElementList<> &elements, gsVector<index_t> &nGaussP)
{
    nGaussP.resize(basis.dim());
    elements=getPatchElements(basis);
    for (int i=0;i<basis.dim();++i)
        nGaussP(i)=basis.maxDegree(i)+1;
}

gsElementList<> getPatchElements(const gsMultiBasis<> &basis)
{
    gsElementList<> elements;
    for (size_t b=0;b<basis.nBases();++b)
        elements.push_back(std::make_pair(basis[b].makeDomainIterator().release(),b));
    return elements;
}

gsElementList<> getBoundaryElements(const gsMultiBasis<> &basis)
{
    gsElementList<> elements;
    const gsBoxTopology &top=basis.topology();
    gsBoxTopology::const_biterator bit=top.bBegin(),bend=top.bEnd();

    for (;bit!=bend;++bit)
        elements.push_back(std::make_pair(basis[bit->patch].makeDomainIterator(bit->side()).release(),bit->patch));
    return elements;
}

std::vector<index_t>             getBoundaryDofs(const gsMultiBasis<> &basis, const gsWeightMapper<real_t> *remap )
{
    std::vector<index_t> result;
    std::vector<index_t> tmp;

    const gsBoxTopology &top=basis.topology();
    gsBoxTopology::const_biterator bit=top.bBegin(),bend=top.bEnd();

    for (;bit!=bend;++bit)
    {
        std::swap(tmp,result);
        gsMatrix<index_t> bDofs=basis[bit->patch].boundaryOffset(bit->side(),0);
        if (remap!=NULL)
        {
            std::vector<index_t> source(reinterpret_cast<index_t*>(bDofs.data()),reinterpret_cast<index_t*>(bDofs.data())+bDofs.size());
            std::vector<index_t> target;
            remap->fastSourceToTarget(source,target);
            std::set_union(tmp.begin(),tmp.end(), target.begin(), target.end(), std::back_inserter(result));
        }
        else
            std::set_union(tmp.begin(),tmp.end(), reinterpret_cast<index_t*>(bDofs.data()), reinterpret_cast<index_t*>(bDofs.data())+bDofs.size(), std::back_inserter(result));

        tmp.clear();
    }
    return result;
}

/**
        \brief Modify the SparseMatrix in place to impose Dirichlet type conditions.

        Change the provided matrix in place in order to impose Dirichlet conditions,
        i.e. it replaces the coefficient (i,j) with
             0 if \f$i\neq j\f$ and i or j is a fixed degree of freedom;
             1 if \f$i=j\f$ and i is a fixed degree of freedom.
        It leaves the other coefficients unchanged.

        While replacing the entries it also compute the required modification to the
        right-hand-side matrix and applies them.

        This operation preserves symmetry of the matrix.

        \param[in,out] M                the system matrix
        \param[in,out] Rhs              the right-hand-side matrix
        \param[in]     dirichletDofs    the indexes of the DirichletDofs
        \param[in]     dirichletValues  the coefficients for the fixed degrees of freedom

        \tparam        T                the type of the coefficients
        \tparam        Major            the storage class of the matrixes
    */



template<typename T, int Major,typename II >
void forceDirichletByColumn (
        gsSparseMatrix<T,Major>      &M,
        gsMatrix<T>                  &Rhs,
        const II                     *beg,
        const II                     *end,
        const gsMatrix<T>            &dirichletValues, real_t coef=1)
{
    GISMO_ASSERT(Major==Eigen::ColMajor,"Something is wrong");
    const II *curR=beg;
    const II *curC=beg;

    for (index_t k=0; k<M.outerSize(); ++k)
    {
        if(curC!=end && *curC == k)
        {
            for (typename gsSparseMatrix<T,Major>::InnerIterator it(M,k); it; ++it)
            {
                Rhs.row(it.row())-=dirichletValues.row(curC-beg)*it.value();
                it.valueRef()=0;
            }
            ++curC;
        }
        else
        {
            curR=beg;
            for (typename gsSparseMatrix<T,Major>::InnerIterator it(M,k); it; ++it)
            {
                while (curR!=end && *curR< it.row() )
                    ++curR;
                if (curR==end)
                    break;
                else if (it.row()==*curR)
                    it.valueRef()=0;
            }
        }
    }
}

template<typename T, int Major,typename II >
void forceDirichletByRow (
        gsSparseMatrix<T,Major>      &M,
        gsMatrix<T>                  &Rhs,
        const II                     *beg,
        const II                     *end,
        const gsMatrix<T>            &dirichletValues, real_t coef=1)
{
    GISMO_ASSERT(Major==Eigen::RowMajor,"Something is wrong");

    const II *curR=beg;
    const II *curC=beg;

    for (index_t k=0; k<M.outerSize(); ++k)
    {
        if (curR!=end && *curR == k)
        {
            for (typename gsSparseMatrix<T,Major>::InnerIterator it(M,k); it; ++it)
            {
                it.valueRef()=0;
            }
            ++curR;
        }
        else
        {
            curC=beg;
            for (typename gsSparseMatrix<T,Major>::InnerIterator it(M,k); it; ++it)
            {
                while (curC!=end && *curC<it.col() )
                    ++curC;
                if (curC==end)
                    break;
                if (it.col()==*curC)
                {
                    Rhs.row(it.row())-=dirichletValues.row(curC-beg)*it.value();
                    it.valueRef()=0;
                    ++curC;
                }
            }
        }
    }
}

template<typename T, int Major, typename II >
void forceDirichletConditions (
        gsSparseMatrix<T,Major>      &M,
        gsMatrix<T>                  &Rhs,
        const II                     *beg,
        const II                     *end,
        const gsMatrix<T>            &dirichletValues, real_t coef=1)
{
    M.makeCompressed();

    if ( Major==Eigen::ColMajor )
        forceDirichletByColumn(M,Rhs,beg,end,dirichletValues,coef=1);
    else if ( Major==Eigen::RowMajor )
        forceDirichletByRow(M,Rhs,beg,end,dirichletValues,coef=1);

    // set the rhs on the fixed dofs, and the unit diagonal of the matrix
    // we can not set the diagonal in the previous loop because it only loop
    // over non zero coefficients
    for (const II *cur=beg; cur!=end; ++cur)
    {
        index_t pos = *cur;
        M(pos,pos)=coef;
        Rhs.row(pos).array()=coef*dirichletValues.row(cur-beg);
    }
    return;
}

template<typename T, int Major, typename IndexM >
void forceDirichletConditions (
        gsSparseMatrix<T,Major>      &M,
        gsMatrix<T>                  &Rhs,
        const IndexM                 &dirDofs,
        const gsMatrix<T>            &dirVals)
{
    if (Rhs.cols()!=dirVals.cols())
    {
        GISMO_ERROR("The matrix for the RHS and for the DirichletValues must have the same size");
    }
    forceDirichletConditions(
                M,Rhs,
                dirDofs.data(),dirDofs.data()+dirDofs.size(),
                dirVals
                );
}

void addAllDirichletBoundaries(const gsMultiPatch<>& mp, gsFunction<real_t> *g,gsBoundaryConditions<real_t> & bcInfo)
{
    std::vector<patchSide> boundaries = mp.boundaries();
    for(unsigned i = 0;i<boundaries.size();++i)
        bcInfo.addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::dirichlet, g);
}

void computeDirichletValues(const gsMultiBasis<> &orig,const gsRecipe<>::spacePtr space, const gsMultiPatch<>& domain, const gsRecipe<>::spacePtr dirichletData, std::vector<index_t> &boundaryDofs, gsMatrix<> &boundaryVal)
{
    // prepare list of meshes
    gsElementList<> elements=getBoundaryElements(orig);

    // prepare the space
    gsRecipe<>::spaceList spaces;
    spaces.push_back( space );
    spaces.push_back( dirichletData );

    // prepare the matrices
    gsSparseMatrix<> B(boundaryDofs.size(),boundaryDofs.size());
    gsMatrix<>       b;
    b.setZero(boundaryDofs.size(),dirichletData->size());

    // decide the quadrature
    gsVector<index_t> nGauss;
    nGauss.setConstant(domain.parDim(),orig.maxCwiseDegree()+1);

    // make recipe
    gsRecipe<> boundaryL2Projection;
    boundaryL2Projection.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()), getFilteredWriter<>(B,boundaryDofs), 0, 0 );
    boundaryL2Projection.add(gsRecipe<>::opHandle(new gsL2ScalarOp<>()), getFilteredRowWriter<>(b,boundaryDofs), 1, 0 );

    // assemble
    boundaryL2Projection.assemble(domain,spaces,nGauss,elements);
    elements.free();

    // solve the L2Projection
    B.makeCompressed();
    boundaryVal = gsSparseSolver<>::CGDiagonal(B).solve(b);
}

} // namespace Utils


