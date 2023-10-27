/**
*    Solving singularly perturbed convection-diffusion equations with stabilized IGA
*
*
*
**/


// Stabilized IGA schemes for convection-diffusion in 1D

#include <gismo.h>
#include <gsCore/gsGeometryEvaluator.h>
#include <gsUtils/gsNorms.h>
#include <fstream>

using namespace gismo;
using math::sinh;
using math::cosh;
using math::tanh;



real_t coth(real_t x)
{
    return 1.0 / math::tanh(x);
}

gsGeometry<>::uPtr makeGeo(int p, int numRefine, bool uniformCP)
{
    gsGeometry<>::uPtr geo( gsNurbsCreator<>::BSplineUnitInterval( p ) );
    // we refine the geometry instead of only the discretization basis
    for (int i = 0; i < numRefine; ++i)
        geo->uniformRefine();

    // linearize control points
    if (uniformCP)
    {
        const int n = geo->coefs().rows();
        geo->coefs() = gsVector<>::LinSpaced( n, geo->coefs()(0), geo->coefs()(n-1) );
    }

    return geo;
}

gsGeometry<>::uPtr smoothingInterpolation( const gsBasis<>& g, const gsFunction<>& f )
{
    gsMatrix<> pts = g.anchors();
    gsMatrix<> fpts = f.eval(pts);
    fpts.transposeInPlace();
    return g.makeGeometry( give(fpts) );
}



template<class T>
gsMatrix<T> solve(const gsSparseMatrix<T> &matrix, const gsMatrix<T> &rhs, bool isSymmetric=false )
{

    // Solve linear system
    gsVector<T> result;

    if ( matrix.cols())
    {
        if (isSymmetric)
        {
            gsSparseSolver<>::CGDiagonal solver;
            //gsSparseSolver<>::CGDiagonal solver;
            solver.setMaxIterations(4 * matrix.rows() );
            result = solver.compute( matrix ).solve (  rhs );
            gsInfo << "residual error: " << solver.error() << "\n";
            gsInfo << "    iterations: " << solver.iterations() << "\n";
            gsInfo << "    Tolerance: " << solver.tolerance() << "\n";
        }
        else
        {
            gsSparseSolver<>::BiCGSTABILUT solver;
            //gsSparseSolver<>::BiCGSTABDiagonal solver;
            //gsSparseSolver<>::BiCGSTABDiagonal solver;

            //gsSparseSolver<>::QR solver;
            //solver.analyzePattern(*m_system.matrix());
            //solver.factorize( *m_system.matrix() );
            //res = solver.solve ( * m_system.rhs() );
            result = solver.compute( matrix ).solve ( rhs );
            gsInfo << "residual error: " << solver.error() << "\n";
            gsInfo << "    iterations: " << solver.iterations() << "\n";
            gsInfo << "    Tolerance: " << solver.tolerance() << "\n";
        }
        //gsSparseSolver<>::SimplicialLDLT solver;
        //gsSparseSolver<>::LU solver;
    }

    return result;
}

template<class T>
memory::unique_ptr<gsGeometry<T> >
reconstructSolution(const gsMatrix<T>& data, const gsMatrix<T>& eliDofs, const gsDofMapper &mapper, const gsBasis<T> &basis, const gsGeometry<T> &geo)
{
    // Reconstruct solution coefficients on patch p
    const int sz  = basis.size();
    const int fsz = mapper.freeSize();
    const int dim = 1;

    gsMatrix<T> coeffs(sz, dim);

    for (index_t i = 0; i < sz; ++i)
    {
        if ( mapper.is_free(i, 0) ) // internal or interface
        {
            for (int k = 0; k < dim; ++k)
                coeffs(i,k) = data( k * fsz + mapper.index(i, 0),0);
        }
        else // eliminated Dof: fill with Dirichlet data
        {
            coeffs.row(i) = eliDofs.row( mapper.bindex(i, 0) );
        }
    }
    return basis.makeGeometry( give(coeffs) );
}


template<class T>
void assembleCDR(gsSparseMatrix<T> &matrix, gsMatrix<T> &rhs, const gsBasis<T>& basis, const gsDofMapper& mapper, const gsMatrix<T> & ddof, const gsConvDiffRePde<T> & pde, const gsGeometry<T> &geo)
{
    gsVector<index_t> numNodes;
    numNodes.setConstant(basis.dim(),basis.maxDegree()+2);
    matrix.resize(mapper.freeSize(),mapper.freeSize());
    rhs.setZero(mapper.freeSize(),1);

    gsMatrix<index_t> actives;
    gsMatrix<T> rhsVals;        // values of the right-hand side
    gsMatrix<T> diffVals, convVals, reacVals;   // PDE coefficient values
    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node

    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell
    gsVector<T> localRhs;       // local load vector

    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM, geo));

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element

    gsGaussRule<T> rule(numNodes);

    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsMatrix<T> basisVal;
    gsMatrix<T> basisDeriv;

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        rule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes,quWeights);

        basis.active_into(quNodes.col(0), actives);
        basis.eval_into  ( quNodes,  basisVal);
        basis.deriv_into ( quNodes,  basisDeriv);

        const index_t numActive=actives.rows();

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(quNodes);

        // evaluate right-hand side and PDE coefficients at the geometry points
        if (pde.rhs())              pde.rhs()->eval_into( geoEval->values(), rhsVals );
        if (pde.diffusion())        pde.diffusion()->eval_into( geoEval->values(), diffVals );
        if (pde.convection())       pde.convection()->eval_into( geoEval->values(), convVals );
        if (pde.reaction())         pde.reaction()->eval_into( geoEval->values(), reacVals );

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);
        localRhs.setZero(numActive);

        for (index_t k = 0; k < quNodes.cols(); ++k)      // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval->measure(k); // weight * abs(det J), where J is geometry Jacobian

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval->transformGradients(k, basisDeriv, trf_grads_k);

            if (pde.rhs())
                localRhs += (weight * rhsVals(0,k)) * basisVal.col(k);

            if (pde.diffusion())
                localStiffness.noalias() += (weight * diffVals(0,k)) * (trf_grads_k.transpose() * trf_grads_k);

            if (pde.convection())
                localStiffness.noalias() += weight * (basisVal.col(k) * (convVals.col(k).transpose() * trf_grads_k));

            if (pde.reaction())
                localStiffness.noalias() += (weight * reacVals(0,k)) * (basisVal.col(k) * basisVal.col(k).transpose());
        }  // end loop Gauss nodes

        // add contributions from local stiffness matrix to global stiffness matrix and load vector
        for (index_t i=0; i < numActive; ++i)
        {
            const int ii = mapper.index(actives(i));
            if ( !mapper.is_free_index(ii) )
                continue;

            rhs(ii,0)+=localRhs(i,0);
            for (index_t j=0; j < numActive; ++j)
            {
                const int jj = mapper.index(actives(j));
                if ( !mapper.is_free_index(jj) )
                    continue;
                matrix.coeffRef(ii, jj) += localStiffness(i, j);
            }
        }

    } // loop over all domain elements

    matrix.makeCompressed();
}


int main(int argc, char** argv)
{
    index_t p = 2;
    index_t numRefine = 3;
    bool uniformCP = false;
    real_t eps = 1e-3;
    bool stabilize = false;
    bool plot = false;
    bool log = false;
    
    gsCmdLine cmd("Solving singularly perturbed convection-diffusion equations with stabilized IGA");
    cmd.addInt("p", "degree", "Degree of the B-spline discretization space", p);
    cmd.addInt("r", "refine",
               "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addReal("a", "diffusion", "Diffusion parameter", eps);
    cmd.addSwitch("stabilize", "Use stabilized scheme", stabilize);
    cmd.addSwitch("uniform", "Use uniformly spaced control points (nonlinear parametrization)",
                  uniformCP);
    cmd.addSwitch("plot", "Plot the result", plot);
    cmd.addSwitch("log", "Create log file", log);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Solving convection-diffusion equation -" << eps << " u'' + u' = 1" << "\n";

    gsGeometry<>::uPtr geo = makeGeo(p, numRefine, uniformCP);
    const real_t h = pow(0.5, numRefine);
    const real_t Pe = 1.0 * h / (2 * eps);      // mesh Peclet number
    real_t sigma = Pe * coth(Pe);         // Il'in-Allen-Southwell factor
    switch (p)
    {
        case 1: sigma = Pe * coth(Pe); break;
        case 2: sigma = (Pe * (6*coth(Pe) + sinh(2*Pe))) / (2 * (2 + cosh(2*Pe))); break;
        case 3: sigma = 2./3.*Pe*( cosh(Pe) * (123 + 56*cosh(2*Pe) + cosh(4*Pe))) / (40*sinh(Pe) + 25*sinh(3*Pe) + sinh(5*Pe)); break;
        case 4: sigma = Pe/4.*(15619*cosh(Pe) + 4293*cosh(3*Pe) + 247*cosh(5*Pe) + cosh(7*Pe))/(1225*sinh(Pe) + 119*(9*sinh(3*Pe) + sinh(5*Pe)) + sinh(7*Pe)); break;
    }

    gsInfo << "Peclet number: " << Pe << "\n";
    if (stabilize)  gsInfo << "Fitting factor: " << sigma << "\n";

    gsConstantFunction<> a(stabilize ? sigma*eps : eps, 1);
    gsConstantFunction<> b(1.0, 1);
    gsConstantFunction<> f(1.0, 1);

    // This isn't the exact solution, but it has maximum error
    //   ~ 2E-9  for eps=0.1
    //   ~ 1E-15 for eps=0.05
    // and has the advantage that it can be evaluated without overflows.
    gsFunctionExpr<> exSol("x - exp((x-1)/w) + exp(-1/w) - exp((x-2)/w)",1);
    exSol.set_w(eps);       // use w for diffusion parameters


    // set up problem data
    gsConvDiffRePde<> pde(&a, &b, 0, &f);
    pde.boundaryConditions().addCondition( boundary::left,  condition_type::dirichlet, 0 );
    pde.boundaryConditions().addCondition( boundary::right, condition_type::dirichlet, 0 );

    gsSparseMatrix<> matrix;
    gsMatrix<>       rhs;
    gsMatrix<>       solCoeff;

    gsMultiBasis<>   basis(geo->basis());
    gsDofMapper      dofMapper;
    basis.getMapper(true,pde.boundaryConditions(), 0, dofMapper);

    gsMatrix<> ddofs;
    ddofs.setZero(dofMapper.boundarySize(),1);

    // assembling amd solving
    assembleCDR(matrix,rhs,basis.basis(0), dofMapper, ddofs, pde,*geo);
    solCoeff=solve(matrix,rhs);

    gsGeometry<>::uPtr igaSol = reconstructSolution(solCoeff,ddofs,dofMapper,basis.basis(0),*geo);
    gsField<> sol(*geo, *igaSol);// consumes igaSol

    gsGeometry<>::uPtr interp ( smoothingInterpolation( geo->basis(), exSol ) );

    real_t l2err  = computeL2Distance     (sol, exSol, false, 9*geo->basis().size());
    real_t maxerr = computeMaximumDistance(sol, exSol, false, 9*geo->basis().size());
    real_t maxerr_intp = computeMaximumDistance(sol, *interp, false, 9*geo->basis().size());

    gsInfo << "L2 error:   " << l2err << "\n";
    gsInfo << "Max. error: " << maxerr << "\n";
    gsInfo << "Max. norm:  " << computeMaximumNorm(sol, 9*geo->basis().size()) << "\n";
    gsInfo << "Max. error to interpolant: " << maxerr_intp << "\n";
    gsInfo << "Max. norm of interpolant:  " << computeMaximumNorm(*geo, *interp, false, 9*geo->basis().size()) << "\n";

    if (log)
    {
        std::fstream logf("convdiff.txt", std::fstream::out | std::fstream::app);
        logf << std::setw(4) << p << " "
             << std::setw(15) << std::scientific << h
             << std::setw(15) << std::scientific << l2err
             << std::setw(15) << std::scientific << maxerr
             << "\n";
    }
    if (plot)
    {
        gsInfo << "Plotting..." << "\n";
        gsWriteParaview( sol, "convdiff", 1000 );
    }
}
