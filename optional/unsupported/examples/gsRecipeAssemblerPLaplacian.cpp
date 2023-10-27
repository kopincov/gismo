/** @file gsCompositeH_test.h

    @brief We solve the p-laplacian using a Newton descent method.
  
    The problem is
    
    \f[ div( |grad u|^(p-2) grad u) = f  in Omega \f]
    u = g                            on the boundary
    
    To use the Newton method we write
    
    \f[ F(u)  = div( |grad u|^(p-2) grad u) \f]
    
    DF(u) = jacobian of F
    Then we use Newton iterations.
    
    Author(s): A. Bressan
*/


#include <gismo.h>
#include <gismo_dev.h>
#include <gsRecipeAssembler/gsPoissonCookBook.h> // for the Dirichlet dofs elimination
#include <gsRecipeAssembler/gsDistanceOperators.h>
#include <gsRecipeAssembler/gsLocalToGlobal.h>
#include <gsAssembler/gsAssemblerUtils.h>
#include <gsMapUtils/gsL2GMapper.h>

#include <gsCore/gsTransformedFuncSet.h>

#include <iostream>
#include <iomanip>

#include <vector>
#include <string>


using namespace gismo;

using std::vector;
using std::string;


class gsWeightedGradGrad : public gsBilinearOp<>
{
private:
    const gsMatrix<> &coeffs;
    real_t            exponent;
public:
    gsWeightedGradGrad(const gsMatrix<> &c, real_t e)
        : coeffs(c), exponent(e-2)
    {}

    virtual void pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t>     result
                        ) const
    {
        unsigned au=unknownSpace.actives().rows();
        unsigned gradSize=unknownSpace.data().derivSize();
        gsVector<> grad(gradSize);
        grad.setZero(gradSize);
        for (unsigned i=0; i<au;++i)
        {
            grad+=unknownSpace.deriv().col(i)*coeffs(unknownSpace.actives()(i));
        }
        const real_t gn   = grad.norm(); // grad norm
        const real_t gn_p = math::pow(gn,exponent);

        result += gn_p*testSpace.deriv().transpose()*unknownSpace.deriv();
    }
    virtual unsigned    testSpaceNeeds() const {return NEED_GRAD;}
    virtual unsigned    unknownSpaceNeeds() const {return NEED_GRAD;}
    virtual unsigned    geometryNeeds()    const {return 0;}
};



class PLapDerOp : public gsBilinearOp<>
{
private:
    const gsMatrix<> &coeffs;
    real_t            exponent;
public:
    PLapDerOp(const gsMatrix<> &c, real_t e)
        : coeffs(c), exponent(e-2)
    {}

    virtual void  pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t>    result
                        ) const
    {
        unsigned au=unknownSpace.actives().rows();
        const index_t gradSize=unknownSpace.data().derivSize();
        gsVector<> grad(gradSize);
        grad.setZero(gradSize);
        for (unsigned i=0; i<au;++i)
        {
            grad+=unknownSpace.deriv().col(i)*coeffs(unknownSpace.actives()(i));
        }
        const real_t gn = grad.norm(); // grad norm
        const real_t gn_p = math::pow(gn,exponent);
        const real_t gn_pmo = math::pow(gn,exponent-1);

        gsMatrix<> C=exponent*gn_pmo*grad*grad.transpose();
        C+=gsVector<>::Constant(gradSize,gn_p).asDiagonal();

        result += testSpace.deriv().transpose()*C*unknownSpace.deriv();
    }

    virtual unsigned    testSpaceNeeds() const {return NEED_GRAD;}
    virtual unsigned    unknownSpaceNeeds() const {return NEED_GRAD;}
    virtual unsigned    geometryNeeds()      const {return 0;}
};

class PLapTestOp : public gsBilinearOp<>
{
private:
    const gsMatrix<> &coeffs;
    real_t            exponent;
public:
    PLapTestOp(const gsMatrix<> &c, real_t e)
        : coeffs(c), exponent(e-2)
    {}


    virtual void pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t>     result
                            ) const
    {
        unsigned au=unknownSpace.actives().rows();
        unsigned gradSize=unknownSpace.data().derivSize();
        gsVector<> grad(gradSize);
        grad.setZero(gradSize);
        for (unsigned i=0; i<au;++i)
        {
            grad+=unknownSpace.deriv().col(i)*coeffs(unknownSpace.actives()(i));
        }
        const real_t gn = grad.norm();; // grad norm

        result += geoEval.measure()*math::pow(gn,exponent)*testSpace.deriv().transpose()*grad;
    }

    virtual unsigned    testSpaceNeeds() const {return NEED_GRAD;}
    virtual unsigned    unknownSpaceNeeds() const {return NEED_GRAD;}
    virtual unsigned    geometryNeeds()      const {return NEED_GRAD|NEED_MEASURE;}
    virtual void        outputSize(unsigned &r,unsigned &c)    const {c=1;}
};



class gsNonlinearSolver
{
    typedef gsRecipe<real_t>::spacePtr spacePtr;
protected:
    gsMatrix<>       m_guess;
    gsMatrix<>       m_full_patch_coefs;

    gsSparseMatrix<> m_sys_m;
    gsMatrix<>       m_sys_rhs;
    gsMatrix<>       m_sys_rhs_base;

    gsSparseMatrix<> m_dirichlet_m;
    gsMatrix<>       m_dirichlet_val;

    gsDofMapper   &m_dof_mapper;

    gsVector<index_t>    m_fill_in;

    gsMultiPatch<>  &m_geo;
    gsMultiBasis<>  &m_space;

    gsMatrix<>       m_distance;
    gsFunction<>    &m_exact_solution;

    real_t           m_tol;
    unsigned         m_max_iter;
protected:
    virtual gsRecipe<> getPatchRecipe (index_t patchId) {return gsRecipe<>();}
public:
    gsNonlinearSolver(
            gsMultiPatch<> &geo,
            gsMultiBasis<> &space,
            gsDofMapper  &dom_mapper,
            gsFunction<>   &source,
            gsFunction<>   &solution,
            gsMatrix<>     &dirichlet_val,
            unsigned        max_it,
            real_t          tol
            )
        : m_dirichlet_val(dirichlet_val),
          m_dof_mapper(dom_mapper),
          m_geo(geo), m_space(space),
          m_distance(1,1),
          m_exact_solution(solution),
          m_tol(tol), m_max_iter(max_it)

    {
        const unsigned   freeSize=m_dof_mapper.freeSize();

        m_sys_rhs_base.setZero(freeSize,source.targetDim());
        for (int patch_id=0; patch_id < m_geo.size();++patch_id)
        {
            gsRecipe<>  m_recipe(1);
            m_recipe[0].setOperator(new gsL2TestOp<>(source));
            m_recipe[0].setTestSpace( 0);
            m_recipe[0].setUnknownSpace( 0);
            m_recipe[0].setRule(new gsL2GMappedRhs<>(m_sys_rhs_base,m_dof_mapper,patch_id));


            std::vector<spacePtr> spList;
            spList.push_back(gsRecipe<real_t>::spacePtr(new gsGradConformingTFS(space[patch_id] )));
            gsVector<index_t>  numNodes = gsAssemblerUtils<real_t>::getNumIntNodesFor( m_space[patch_id]);
            gsGaussRule<>  quadrature(numNodes);
            gsBasis<>::domainIter domIt=m_space[patch_id].makeDomainIterator();
            m_recipe.assemble(quadrature, *domIt, m_geo.patch(patch_id), spList);
        }

        short_t             maxD=0;
        for (size_t patchId=0; patchId < geo.nPatches(); ++patchId )
            maxD = std::max(maxD, space[patchId].maxDegree());
        m_fill_in.setConstant(m_dof_mapper.freeSize(), math::ipow( 2*maxD+1, space[0].dim()));

    }

    void assemble()
    {
        const unsigned   freeSize=m_dof_mapper.freeSize();
        const unsigned   boundarySize=m_dof_mapper.boundarySize();

        m_sys_m.resize(freeSize,freeSize);
        m_sys_m.setZero();
        m_sys_m.reserve(m_fill_in);
        m_sys_rhs.setZero(freeSize,1);

        m_dirichlet_m.resize(freeSize,boundarySize);
        m_dirichlet_m.setZero();

        m_distance.setZero(1,1);

        for (index_t patch_id=0; patch_id < m_geo.size();++patch_id)
        {
            m_full_patch_coefs=reconstructFullCoefs(patch_id);
            gsRecipe<> m_problem=getPatchRecipe(patch_id);
            gsRecipe<> l2dist(1);
            l2dist[0].setOperator(new gsDistL2(m_full_patch_coefs, &m_exact_solution));
            l2dist[0].setTestSpace( 0);
            l2dist[0].setUnknownSpace( 0);
            l2dist[0].setRule(new gsL2GPlain< gsMatrix<> > (m_distance));
            m_problem.append(l2dist);

            std::vector<spacePtr> spList;
            spList.push_back(spacePtr(new gsGradConformingTFS(m_space[patch_id])));
            gsVector<index_t>  numNodes = gsAssemblerUtils<real_t>::getNumIntNodesFor( m_space[patch_id]);
            gsGaussRule<>  quadrature(numNodes);
            gsBasis<>::domainIter domIt=m_space[patch_id].makeDomainIterator();
            m_problem.assemble(quadrature, *domIt, m_geo.patch(patch_id), spList);
        }
    }
    gsMatrix<> run(const gsMatrix<> &guess)
    {
        m_guess=guess;
        for (int i=1; true; ++i)
        {
            assemble();
            gsMatrix<> delta_x = solveLinearStep();
            if ( stopCriteria(delta_x, i) )
                break;
            m_guess += delta_x;
        }
        return m_guess;
    }

    gsMatrix<> solveLinearStep ()
    {
        gsSparseSolver<>::BiCGSTABILUT  solver;
        solver.analyzePattern( m_sys_m );
        solver.factorize( m_sys_m );
        return solver.solve ( m_sys_rhs );
    }
    bool stopCriteria(gsMatrix<> &delta_x, unsigned it)
    {
        gsInfo<<std::setprecision(6)<< std::scientific;
        gsInfo
            <<it
            <<"  L2 dist = "<<math::sqrt(m_distance(0,0))
            <<"  delta_x  "<<delta_x.norm()<<"\n";
        return     delta_x.norm() <= m_tol
                || (unsigned)it>m_max_iter
                || !math::isfinite(delta_x.norm());
    }

    gsField<> reconstructSolution()
    {
        gsPiecewiseFunction<> * sols = new gsPiecewiseFunction<>;
        for (size_t np=0; np < m_geo.nPatches(); ++np )
            sols->addPiecePointer( reconstructPatchSolution( np ) );
        return gsField<>( m_geo, gsPiecewiseFunction<>::Ptr(sols), true );
    }

    gsGeometry<>::uPtr reconstructPatchSolution( int p)
    {
        return m_space[p].makeGeometry(reconstructFullCoefs(p));
    }

    gsMatrix<>    reconstructFullCoefs( int p)
    {
        const int sz  = m_space[p].size();
        const int dim = 1;

        gsMatrix<> coeffs(sz, dim);
        for (index_t i = 0; i < sz; ++i)
        {
            if ( m_dof_mapper.is_free(i, p) ) // internal or interface
                coeffs.row(i) = m_guess.row(m_dof_mapper.index(i, p));
            else // eliminated Dof: fill with Dirichlet data
                coeffs.row(i) = m_dirichlet_val.row( m_dof_mapper.bindex(i, p) );
        }
        return coeffs;
    }

};


class PLapNewton
    : public gsNonlinearSolver
{
protected:
    using gsNonlinearSolver::m_full_patch_coefs;
    using gsNonlinearSolver::m_sys_m;
    using gsNonlinearSolver::m_sys_rhs;
    using gsNonlinearSolver::m_dirichlet_m;

    real_t          m_exp;

    virtual gsRecipe<> getPatchRecipe (index_t patchId)
    {
        gsRecipe<> m_recipe(2);
        m_recipe[0].setOperator(new PLapDerOp(m_full_patch_coefs, m_exp));
        m_recipe[0].setTestSpace( 0);
        m_recipe[0].setUnknownSpace( 0);
        m_recipe[0].setRule(new gsL2GMapped<>(m_sys_m,m_dirichlet_m,m_dof_mapper,m_dof_mapper,patchId) );

        m_recipe[1].setOperator(new PLapTestOp(m_full_patch_coefs, m_exp));
        m_recipe[1].setTestSpace( 0);
        m_recipe[1].setUnknownSpace( 0);
        m_recipe[1].setRule(new gsL2GMappedRhs<>(m_sys_rhs,m_dof_mapper,patchId));
        return m_recipe;
    }
public:
    PLapNewton (
            real_t          exp,
            gsMultiPatch<> &geo,
            gsMultiBasis<> &space,
            gsDofMapper  &dom_mapper,
            gsFunction<>   &source,
            gsFunction<>   &solution,
            gsMatrix<>     &dirichlet_val,
            unsigned        max_it,
            real_t          tol
            )
            : gsNonlinearSolver(geo,space,dom_mapper,source,solution,dirichlet_val,max_it,tol),
            m_exp(exp)
    {

    }
    gsMatrix<> run(const gsMatrix<> &guess)
    {
        m_guess=guess;
        for (int i=1; true; ++i)
        {
            assemble();
            m_sys_rhs-=m_sys_rhs_base;
            gsMatrix<> delta_x = solveLinearStep();
            if ( stopCriteria(delta_x, i) )
                break;
            m_guess -= delta_x;
        }
        return m_guess;
    }
};

class PLapSequence
    : public gsNonlinearSolver
{
protected:
    using gsNonlinearSolver::m_full_patch_coefs;
    using gsNonlinearSolver::m_sys_m;
    using gsNonlinearSolver::m_sys_rhs;
    using gsNonlinearSolver::m_dirichlet_m;

    real_t          m_exp;

    virtual gsRecipe<> getPatchRecipe (index_t patchId)
    {
        gsRecipe<> m_recipe(1);
        m_recipe[0].setOperator(new gsWeightedGradGrad(m_full_patch_coefs, m_exp));
        m_recipe[0].setTestSpace( 0);
        m_recipe[0].setUnknownSpace( 0);
        m_recipe[0].setRule(new gsL2GMapped<>
            (m_sys_m,m_dirichlet_m,m_dof_mapper,m_dof_mapper,patchId));
        return m_recipe;
    }
public:
    PLapSequence (
            real_t          exp,
            gsMultiPatch<> &geo,
            gsMultiBasis<> &space,
            gsDofMapper  &dom_mapper,
            gsFunction<>   &source,
            gsFunction<>   &solution,
            gsMatrix<>     &dirichlet_val,
            unsigned        max_it,
            real_t          tol
            )
            : gsNonlinearSolver(geo,space,dom_mapper,source,solution,dirichlet_val,max_it,tol),
            m_exp(exp)
    {
            m_sys_rhs=m_sys_rhs_base;
    }
    gsMatrix<> run(const gsMatrix<> &guess)
    {
        m_guess=guess;
        for (int i=1; true; ++i)
        {
            assemble();
            m_sys_rhs=m_sys_rhs_base-m_sys_m*m_guess-m_dirichlet_m*m_dirichlet_val;
            gsMatrix<> delta_x = solveLinearStep();
            if ( stopCriteria(delta_x, i) )
                break;
            m_guess += delta_x;
        }
        return m_guess;
    }
};

int main (int argc, char **argv)
{
    // configuration variables
    gsMultiPatch<>::uPtr geo;             // the domain geometry
    index_t           h_ref = 3;           // the number of dyadic refinements of the domain map
    index_t           p_ref = 1;           // the number of degree elevations
    real_t            test_tol = 10e-13;
    real_t            norm_tol = 10e-1;
    index_t           maxIter = 75;
    real_t            exponent = 2.5;
    string            inName = "square.xml"; //planar/multipatch_triangle.xml",
    
    gsCmdLine cmd("solving P-laplacian");
    cmd.addPlainString("filename", "File containing the input", inName);
    cmd.addReal("p", "exponent", "the exponenet p in the equation", exponent);
    cmd.addReal("t", "e_tollerance", "the tollerance in the residual", test_tol);
    cmd.addReal("T", "l2_tollerence", "the tollerance in the norm", norm_tol);
    cmd.addInt("r", "refinement", "number of uniform h-refinement steps", h_ref);
    cmd.addInt("e", "elevation", "number of degree elevation steps", p_ref);
    cmd.addInt("i", "iteration", "maximum number of iterations", maxIter);      
        
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    geo = gsReadFile<>(inName);
    if (!geo)
    {
        GISMO_ERROR("No geometry provided in the file "<<inName<<"\n");
    }


    /// Set up the problem data according to the given
    /// geometry and exponent


    std::stringstream f_string;
    f_string<< "(pow(2,"<<exponent<<"/2.)*pow(5,"<<exponent<<")*pow(2-cos(20*x)+cos(20*y),"<<exponent<<"/2.)*((cos(10*x)*("<<exponent<<"-(-1+"<<exponent<<")*cos(20*x)+cos(20*y)))/2.+((-1+"<<exponent<<")*pow(cos(10*y),2)+pow(sin(10*x),2))*sin(10*y)))/pow(pow(cos(10*y),2)+pow(sin(10*x),2),2)";

    gsFunctionExpr<>    u("sin(10*y)+cos(10*x)", 2);
    gsFunctionExpr<>    f(f_string.str(),2);

    gsBoundaryConditions<> bcd;
    gsMultiPatch<>::const_biterator bit;
    for (bit = geo->bBegin(); bit != geo->bEnd(); ++bit)
        bcd.addCondition( *bit, condition_type::dirichlet, &u);
    /// printout data
    gsInfo
            <<"Solving th p-Laplacian equation\n"
           <<"\n\t div ( |grad u|^(p-2) grad u ) = f"
          <<"\n\nfor p = "<<exponent<<" and  f = "<<f<<"\n"
         <<"\non Omega = "<<*geo<<"\n";

    gsInfo
            <<"\nThe exact solution is\n"
           <<u<<"\n";

    gsInfo
            <<"\nDiscretization parameters are:"
           <<"\nelevetate degree  "<<p_ref<<" times"
          <<"\nrefine dyadically "<<h_ref<<" times"
         <<"\nthe tollerence in the euclidean vector norm is < "<<test_tol
        <<"\nthe tollerence in the L2 norm is < "<<test_tol
       <<"\nthe maxximum number of iterations is "<<maxIter<<"\n";


    /// Initialization of the discrete space


    gsMultiBasis<>     mySpace(*geo);
    for (int k = 0; k < p_ref; ++k)
        mySpace.degreeElevate();
    for (int k = 0; k < h_ref; ++k)
        mySpace.uniformRefine();
    gsInfo<<"discreteSpace is"<<mySpace<<"\n";

    gsDofMapper * dofMapper = new gsDofMapper;
    mySpace.getMapper(true, bcd, *dofMapper);

    /// Solve Boundary conditions

    gsSparseMatrix<> dirichlet_sys(dofMapper->boundarySize(),dofMapper->boundarySize());
    gsMatrix<>       dirichlet_rhs;
    gsMatrix<>       dirichlet_val;
    dirichlet_rhs.setZero(dofMapper->boundarySize(),u.targetDim());

    for (unsigned b_id=0; b_id<bcd.dirichletSides().size(); ++b_id)
    {
        unsigned       patchId = bcd.dirichletSides()[b_id].patch();
        boxSide sideId  = bcd.dirichletSides()[b_id].side();
        gsFunction<>  *func    = bcd.dirichletSides()[b_id].function().get();


        gsVector<index_t> numNodes = gsAssemblerUtils<real_t>::getNumIntNodesForSide( mySpace[patchId], sideId.direction() );
        gsBasis<>::domainIter domIt=mySpace[patchId].makeDomainIterator(sideId);
        gsGaussRule<>       quadrature(numNodes);

        std::vector<gsRecipe<real_t>::spacePtr> spList;
        spList.push_back(gsRecipe<real_t>::spacePtr(new gsGradConformingTFS(mySpace[patchId])));
        gsRecipe<> dirichlet_recipe = gsPoissonBoundaryCookBook::makeRecipeForDiricletValues (*func, *dofMapper, dirichlet_sys, dirichlet_rhs, patchId);
        dirichlet_recipe.assemble(quadrature, *domIt, geo->patch(patchId), spList);
    }
    gsSparseSolver<>::CGDiagonal solver;
    dirichlet_val = solver.compute(dirichlet_sys).solve(dirichlet_rhs);

    /// initialize matrixes


    const unsigned size = dofMapper->freeSize();
    gsMatrix<>       initial_guess(size,1);
    // we cannot start from something with 0 gradient as the matrix would be singular
    for (unsigned r=0;r<size;++r)
    {
        initial_guess(r,0)=(real_t)(r)/size;
    }

    PLapNewton my_newton_solver(exponent, *geo, mySpace, *dofMapper, f, u, dirichlet_val, maxIter, test_tol);
    my_newton_solver.run(initial_guess);
    PLapSequence my_seq_solver(exponent, *geo, mySpace, *dofMapper, f, u, dirichlet_val, maxIter, test_tol);
    my_seq_solver.run(initial_guess);



    gsField<> newton_plap_sol = my_newton_solver.reconstructSolution();
    gsField<> lapseq_plap_sol = my_seq_solver.reconstructSolution();
    gsField<> exact_plap_sol(*geo,u,false);

    gsWriteParaview<>( newton_plap_sol, "p-lap-newton", 2500);
    gsWriteParaview<>( lapseq_plap_sol, "p-lap-lapseq", 2500);
    gsWriteParaview<>( exact_plap_sol,  "p-lap-exact", 2500);

    real_t distN = igaFieldL2Distance (newton_plap_sol, u);
    real_t distS = igaFieldL2Distance (lapseq_plap_sol, u);

    gsInfo<<"L2 distance for Newton and Sequence: "<<distN<<",  "<<distS<<"\n";

    delete dofMapper;
    return distN > norm_tol || distS > norm_tol;
}


