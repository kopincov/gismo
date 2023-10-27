
    // Finding eigenvalues of the linear Euler-Bernoulli beam problem

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

#include <gismo.h>


using namespace gismo;

void computeEigenvalues(const gsSparseMatrix<>& As, const gsSparseMatrix<>& Ms, gsVector<>& ev)
{
    // convert to dense matrices (sparse eigensolver would be nice!)
    gsMatrix<> A = As;
    gsMatrix<> M = Ms;

    //gsInfo << A << "\n";
    //gsInfo << M << "\n";

    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<>::Base > eig;
    eig.compute( A, M, Eigen::ComputeEigenvectors );

    ev = eig.eigenvalues();
    gsInfo << eig.eigenvectors().col( A.cols()-2 ) << "\n";
}


real_t laplaceEV(int i)
{
    return (i+1)*(i+1) * EIGEN_PI*EIGEN_PI;
}


int main(int argc, char** argv)
{
    index_t p = 2;
    index_t numRefine = 3;
    bool uniformCP = false;
    
    gsCmdLine cmd("Computes eigenvalues of the Laplace problem using isogeometric analysis.");
    cmd.addInt("p", "degree",
               "Degree of the B-spline discretization space", p);
    cmd.addInt("r", "uniformRefine",
               "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addSwitch("uniform", "Use uniformly spaced control points (nonlinear parametrization)", uniformCP);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    gsInfo << "Determining eigenfrequencies of the Euler-Bernoulli beam using B-splines of degree " << p << "..." << "\n";

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
    gsBoundaryConditions<> BCs;
    BCs.addCondition( boundary::left,  condition_type::dirichlet, 0 );
    BCs.addCondition( boundary::right, condition_type::dirichlet, 0 );

    // set up discretization space
    gsMultiBasis<> basis( geo->basis() );

    gsGenericAssembler<real_t> solver(gsMultiPatch<>(*geo), basis,
                                      gsAssembler<>::defaultOptions(), &BCs);

    gsInfo << "Discretization Space:\n" << basis << "\n";

    //gsInfo << "Control points:\n" << geo->coefs().transpose() << "\n";

    gsInfo << "Assembling stiffness matrix, " << std::flush;
    gsSparseMatrix<> S = solver.assembleStiffness();
    gsInfo << "mass matrix, " << std::flush;
    gsSparseMatrix<> M = solver.assembleMass();
    gsInfo << "done.\n" << "\n";

    //gsInfo << "Stiffness matrix:\n";
    //gsMatrix<> dense = S.toDense();
    //dense = dense.selfadjointView<Lower>();
    //gsInfo << dense << "\n" << "\n";

    //gsInfo << "Mass matrix:\n";
    //dense = M.toDense();
    //dense = dense.selfadjointView<Lower>();
    //gsInfo << dense << "\n" << "\n";

    gsInfo << "\nSolving generalized eigenvalue problem... " << std::flush;
    gsVector<> ev;
    computeEigenvalues( S, M, ev);
    gsInfo << "done." << "\n";

    gsInfo << "\n EV nr.    exact   approx.   ratio\n" ;
    gsInfo <<   "=============================================\n";
    gsInfo.setf(std::ios::showpoint);
    real_t relErr = 0;

    for (index_t i = 0; i < ev.size(); ++i)
    {
        const real_t ex_ev = laplaceEV(i);
        const real_t approx_ev = ev[i];

        gsInfo << "  " << std::setw(5) << (i+1)
            << "  "  << std::setw(8) << ex_ev
            << "  "  << approx_ev
            << "  "  << approx_ev / ex_ev
            << "\n";
        relErr = math::max(relErr, math::abs(ex_ev - approx_ev) / ex_ev);
    }
    gsInfo << "\n";
    gsInfo << "Max. relative error: " << std::setprecision(1) << std::fixed << 100*relErr << " %" << "\n";

    //gsInfo << "\n" << "Laplacian eigenvalues:\n" << ev.array().sqrt().transpose() << "\n";
    std::fstream log("out.txt", std::fstream::out | std::fstream::app);
    log << p << "\t" << numRefine << "\t" << ev.size() << "\t" << relErr << "\n";

    return 0;
}
