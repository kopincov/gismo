
// Finding eigenvalues of the linear Euler-Bernoulli beam problem

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

#include <gismo.h>

#include <gsCore/gsGeometryEvaluator.h>

using namespace std;
using namespace gismo;

void computeEigenvalues(const gsSparseMatrix<>& As, const gsSparseMatrix<>& Ms, gsVector<>& ev)
{
    // convert to dense matrices (sparse eigensolver would be nice!)
    gsMatrix<> A = As;
    gsMatrix<> M = Ms;

    //cout << A << "\n";
    //cout << M << "\n";

    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<>::Base > eig;
    eig.compute( A, M, Eigen::EigenvaluesOnly );

    ev = eig.eigenvalues();
}


// First few eigenvalues, computed by Mathematica.
// For the rest, the formula (2*(i+1)+1) * EIGEN_PI / 2 produces
// a result which is accurate to machine precision.
static const real_t clampedClampedEVs[] = {
    4.730040744862704,
    7.853204624095838,
    10.995607838001671,
    14.137165491257464,
    17.27875965739948,
    20.42035224562606,
    23.561944902040455,
    26.703537555508188,
    29.845130209103253,
    32.98672286269282,
    36.128315516282626,
    39.269908169872416,
    42.411500823462205,
    45.553093477052,
    48.6946861306418
};

real_t clampedClampedEV(int i)
{
    if (static_cast<size_t>(i) < sizeof(clampedClampedEVs) / sizeof(clampedClampedEVs[0]))
        return clampedClampedEVs[i];
    else
        return (2*(i+1)+1) * EIGEN_PI / 2;
}


real_t clampedClampedProblem(const gsSparseMatrix<>& As, const gsSparseMatrix<>& Ms, gsVector<>& ev)
{
    const int N = As.rows();

    gsMatrix<> A = As.block(1, 1, N-2, N-2);
    gsMatrix<> M = Ms.block(1, 1, N-2, N-2);

    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<>::Base > eig;
    eig.compute( A, M, Eigen::EigenvaluesOnly );

    ev = eig.eigenvalues();

    cout << "\n" << "\n" << "   --- CLAMPED-CLAMPED PROBLEM ---";
    cout << "\n EV nr.    exact   approx.   ratio\n" ;
    cout <<   "=============================================\n";
    cout.setf(ios::showpoint);

    real_t relErr1 = 0;
    static const int RelErrRange = 1;
    for (int i = 0; i < ev.size(); ++i)
    {
        real_t ex_ev = clampedClampedEV(i);
        ex_ev *= ex_ev;
        const real_t approx_ev = math::sqrt(ev[i]);

        cout << "  " << setw(5) << (i+1)
             << "  "  << setw(8) << ex_ev
             << "  "  << approx_ev
             << "  "  << approx_ev / ex_ev
             << "\n";
        if (i < RelErrRange)
            relErr1 = math::max(relErr1, math::abs(ex_ev - approx_ev) / ex_ev);
    }

    return relErr1;
}


void assembleMatrices (gsGeometry<> * geo, gsSparseMatrix<> &mass, gsSparseMatrix<> &system, gsMatrix<> &rhsMatrix, gsFunction<>* rhsFunction=NULL)
{
    gsBoundaryConditions<> bvp;
    bvp.addCondition( boundary::left,  condition_type::dirichlet, 0 );
    bvp.addCondition( boundary::right, condition_type::dirichlet, 0 );

    // set up discretization space
    gsBasis<>::uPtr basis = geo->basis().clone();

    cout << "Control points:\n" << geo->coefs().transpose() << "\n";

    cout << "Assembling matrices" << flush;

    gsMultiBasis<> mbasis(*basis);
    gsDofMapper    mapper;
    mbasis.getMapper(true, bvp,mapper,true);
    gsMatrix<> ddof(2,1);
    ddof<<0,0;

    mass.resize(mapper.freeSize(),mapper.freeSize());
    system.resize(mapper.freeSize(),mapper.freeSize());
    rhsMatrix.setZero(mapper.freeSize(),1);

    GISMO_ASSERT(geo->parDim() ==1,"Beam must be 1D");

    gsVector<index_t> numNodes;
    numNodes.setConstant(1, basis->degree(0)+1);
    gsGaussRule<> quad(numNodes);
    typename gsGeometryEvaluator<>::uPtr geoEval(getEvaluator(NEED_VALUE | NEED_JACOBIAN | NEED_2ND_DER, *geo));
    // Make domain element iterator
    gsBasis<>::domainIter domIt = basis->makeDomainIterator();


    // per element locals
    gsMatrix<> rhsVals;        // values of the right-hand side
    gsMatrix<> localStiffness; // (dense) stiffness matrix within one grid cell
    gsMatrix<> localMass;      // (dense) mass matrix within one grid cell
    gsVector<> localRhs;       // local load vector

    gsMatrix<> quNodes;
    gsVector<> quWeights;

    gsMatrix<index_t> basActive;
    vector<gsMatrix<> >basAllValues;

    // per point locals
    gsMatrix<> physDer;

    for (; domIt->good(); domIt->next() )
    {
        // evaluate
        quad.mapTo(domIt->lowerCorner(),domIt->upperCorner(),quNodes,quWeights);
        geoEval->evaluateAt(quNodes);

        basis->active_into(domIt->centerPoint(),basActive);
        const int numActive = basActive.rows();

        basis->evalAllDers_into(quNodes,2,basAllValues);
        const gsMatrix<> & basVal  = basAllValues[0];
        const gsMatrix<> & basDer1 = basAllValues[1];
        const gsMatrix<> & basDer2 = basAllValues[2];

        const gsMatrix<>& geoDer1 = geoEval->jacobians();
        const gsMatrix<>& geoDer2 = geoEval->derivs2();
        // evaluate right-hand side at the geometry points
        if (rhsFunction) rhsFunction->eval_into( geoEval->values(), rhsVals );

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);
        localMass.setZero(numActive, numActive);
        localRhs.setZero(numActive);

        for (index_t k = 0; k < quWeights.size(); ++k)      // loop over quadrature nodes
        {
            const real_t weight = quWeights(k) * geoDer1(0,k);
            // transform parametric 2nd basis derivatives into physical ones
            physDer = (basDer2.col(k) - (geoDer2(0,k)/geoDer1(0,k)) * basDer1.col(k)) / (geoDer1(0,k)*geoDer1(0,k));

            if (rhsFunction) localRhs += (weight * rhsVals(0,k)) * basVal.col(k);
            localStiffness.noalias() += weight * (physDer * physDer.transpose());
            localMass.noalias()      += weight * basVal.col(k)*basVal.col(k).transpose();
        }  // end loop Gauss nodes

        // map to global matrix
        for (index_t i=0; i < numActive; ++i)
        {
            const int ii = mapper.index(basActive(i));
            if ( !mapper.is_free_index(ii) )
                continue;

            rhsMatrix(ii,0)+=localRhs(i,0);
            for (index_t j=0; j < numActive; ++j)
            {
                const int jj = mapper.index(basActive(j));
                if ( !mapper.is_free_index(jj) )
                    continue;
                system.coeffRef(ii, jj) += localStiffness(i, j);
                mass.coeffRef(ii, jj)   += localMass(i, j);
            }
        }

    }
    mass.makeCompressed();
    system.makeCompressed();

//    gsMatrix<> denseMass=mass;
//    gsMatrix<> denseSys=system;
//    cout<<"\n\n"<<denseMass<<"\n\n"<<denseSys<<"\n\n"<<flush;
}




int main(int argc, char** argv)
{
    index_t p = 2;
    index_t numRefine = 3;
    bool uniformCP = false;
    bool scaleMass = false;

    gsCmdLine cmd("Computes eigenvalues of the linear Euler-Bernoulli beam using isogeometric analysis.");
    cmd.addInt("p", "degree",
               "Degree of the B-spline discretization space", p);
    cmd.addInt("r", "uniformRefine",
               "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addSwitch("uniform", "Use uniformly spaced control points (nonlinear parametrization)", uniformCP);
    cmd.addSwitch("scale-mass", "Use scaling of mass matrix", scaleMass);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    cout << "Determining eigenfrequencies of the Euler-Bernoulli beam using B-splines of degree " << p << "..." << "\n";

    gsGeometry<>::uPtr geo( gsNurbsCreator<>::BSplineUnitInterval( p ) );
    // we refine the geometry instead of only the discretization basis
    for (index_t i = 0; i < numRefine; ++i)
        geo->uniformRefine();

    // linearize control points
    if (uniformCP)
    {
        const index_t n = geo->coefs().rows();
        geo->coefs() = gsVector<>::LinSpaced( n, geo->coefs()(0), geo->coefs()(n-1) );
    }

    // set up boundary value problem

    gsSparseMatrix<> massMatrix,systemMatrix;
    gsMatrix<> rhs;

    assembleMatrices(geo.get(), massMatrix,systemMatrix,rhs);

#if 0
    {
        const index_t N = basis->size();
        cout << "Basis size: " << N << "\n";
        gsMatrix<> coeff(N,1);

        gsGeometry<>::uPtr bfun = basis->makeGeometry( coeff );

        for (index_t j = 0; j < N; ++j)
        {
            for (index_t i = 0; i < N; ++i)
                coeff(i,0) = (i==j) ? 1.0 : 0.0;
            bfun->setCoefs( coeff );

            gsMatrix<> x(1,2);
            x(0,0) = 0.0; x(0,1) = 1.0;

            cout << "Basis function " << setw(2) << j << " evaluated at 0 and 1:";
            cout << "  values: " << setw(4) << bfun->eval(x)
                 << "   1st der.: " << setw(4) << bfun->deriv(x)
                 << "   2nd der.: " << setw(4) << bfun->deriv2(x)
                 << "\n";
        }
    }
#endif

    if (scaleMass)
    {
        const index_t n = massMatrix.rows();
        gsVector<> diag(n);

        diag.fill(p+1);
        for (index_t i = 0; i < p; ++i)
            diag[i] = diag[n-i-1] = i+1;
        diag = (diag / (p+1)).array().pow(-1.0);

        massMatrix = diag.asDiagonal() * massMatrix * diag.asDiagonal();
    }

    cout << setprecision(9);
    cout << endl << massMatrix << endl << systemMatrix << endl;


    cout << "\nSolving eigenvalue problem... " << flush;
    gsVector<> ev;
    computeEigenvalues(systemMatrix, massMatrix, ev);
    cout << "done." << "\n";

    cout << "\n EV nr.    exact   approx.   ratio   rel.err.(%)\n" ;
    cout <<   "=================================================\n";
    cout.setf(ios::showpoint);
    real_t relErr = 0;
    real_t relErr1 = 0;
    static const index_t RelErrRange = 1;
    for (index_t i = 0; i < ev.size(); ++i)
    {
        const real_t ex_ev = EIGEN_PI*EIGEN_PI*(i+1)*(i+1);
        const real_t approx_ev = math::sqrt(ev[i]);

        cout << setprecision(6)
             << "  "  << setw(5) << (i+1)
             << "  "  << setw(8) << ex_ev
             << "  "  << approx_ev
             << "  "  << approx_ev / ex_ev
             << "  "  << setprecision(4) << 100.0 * abs(ex_ev - approx_ev) / ex_ev
             << "\n";
        relErr = math::max<real_t>(relErr, math::abs(ex_ev - approx_ev) / ex_ev);
        if (i < RelErrRange)
            relErr1 = math::max(relErr1, math::abs(ex_ev - approx_ev) / ex_ev);
    }
    cout << "\n";
    //cout << "Max. relative error: " << setprecision(1) << fixed << 100*relErr << " %" << "\n";

    //relErr1 = clampedClampedProblem(systemMatrix, massMatrix, ev);
    //cout << "\n" << "Clamped-clamped eigenvalues:\n" << ev.array().sqrt().transpose() << "\n";

    fstream log("out.txt", fstream::out | fstream::app);
    log << p << "\t" << numRefine << "\t" << ev.size() << "\t" << relErr1 << "\n";

}
