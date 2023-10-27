/** @file gsMultiGridExample.cpp

    @brief Provides test examples for multigrid algorithms

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsSolver/gsTimedOp.h>
#include <gsMultiGrid/gsMassSmoother.h>
#include <gsCore/gsLinearCombinationOfFunctionsFunction.h>

//#include <gsIO/gsMatrixIO.h>
#include <gsTensor/gsTensorTools.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsGeneralizedPoissonAssembler.h>

#include <gsIO/gsCmdLineWithEnumSupport.h>
#include <gsCore/gsConnectedBoundaryComponent.h>

#include <gsUtils/gsNorms.h>

using namespace std;
using namespace gismo;

namespace Smoother {
    enum type {
        Richardson,
        Jacobi,
        GaussSeidel,
        MassRichardson,
        MassRichardsonBoundaryCorrection,
        MassRichardsonSubspaceCorrection,
        MassRichardsonSubspaceCorrectionAdditiveMPDDD,
        MassRichardsonSubspaceCorrectionAdditiveMPDDD2,
        MassRichardsonSubspaceCorrectionAdditiveMPDDDRep,
        MassRichardsonSubspaceCorrectionAdditiveMPNDD,
        MassRichardsonSubspaceCorrectionMultiplicativeMPNDD,
        MassRichardsonSubspaceCorrectionGeo,
        MassRichardsonSubspaceCorrectionGS
    };
}

/// A locally used geometry
gsTensorBSpline<2>::uPtr BSplineMySquare(short_t deg)
{
    gsTensorBSpline<2>::uPtr res(gsNurbsCreator<>::BSplineSquareDeg(deg));
    res->insertKnot( 0.5, 0 );
    return res;
}

/// A locally used geometry
gsTensorBSpline<2>::uPtr BSplineQuarterAnnulus()
{
    gsKnotVector<> KV(0,1, 0,3);

    gsMatrix<> C(9,2);
    C  << 1,   0,
          1.5, 0,
          2,   0,
          1,   1,
          1.5, 1.5,
          2,   2,
          0,   1,
          0,   1.5,
          0,   2;

    return memory::make_unique(new gsTensorBSpline<2>(KV,KV, give(C)));
}

/// A locally used geometry
gsGeometry<>::uPtr approximateQuarterAnnulus(int deg)
{
    gsGeometry<>::uPtr quann = gsNurbsCreator<>::NurbsQuarterAnnulus();

    gsKnotVector<> KV1(0,1, 0, deg+1);        // no interior knots in x direction
    gsKnotVector<> KV2(0,1, 1, deg+1);        // 1 interior knot in y direction

    gsTensorBSplineBasis<2> tbsp (new gsBSplineBasis<>(KV1), new gsBSplineBasis<>(KV2));
    return tbsp.interpolateAtAnchors( quann->eval(tbsp.anchors()) );
}

/// Completes the vector of free dofs with the given boundary dofs and stores the result in \a freeDofs.
void completeVector(const gsDofMapper& mapper, const gsMatrix<real_t>& freeDofs, const gsMatrix<real_t>& boundaryDofs, gsMatrix<real_t>& allDofs)
{
    const int nDofs = mapper.size();

    assert (freeDofs.size() == mapper.freeSize());
    assert (boundaryDofs.size() == mapper.boundarySize());
    allDofs.resize(nDofs, 1);

    for (index_t i = 0; i < nDofs; ++i)
    {
        allDofs(i) = mapper.is_free(i)
            ? freeDofs( mapper.index(i) )
            : boundaryDofs( mapper.bindex(i) );
    }
}


/// Completes the vector of free dofs with zero boundary dofs and stores the result in \a freeDofs.
void completeVector(const gsDofMapper& mapper, const gsMatrix<real_t>& freeDofs, gsMatrix<real_t>& allDofs)
{
    const int nDofs = mapper.size();

    assert (freeDofs.size() == mapper.freeSize());

    allDofs.setZero(nDofs, 1);

    for (index_t i = 0; i < nDofs; ++i)
    {
        if (mapper.is_free(i))
            allDofs(i) = freeDofs( mapper.index(i) );
    }
}


/// Extracts the free dofs from the given vector of all dofs and stores the result in \a freeDofs.
void extractFreeDofs(const gsDofMapper& mapper, const gsMatrix<real_t>& allDofs, gsMatrix<real_t>& freeDofs)
{
    const int nFreeDofs = mapper.freeSize();

    assert (allDofs.size() == mapper.size());
    freeDofs.resize(nFreeDofs, 1);

    for (index_t j = 0; j < allDofs.size(); ++j)
    {
        const int gj = mapper.index(j, 0);
        if (mapper.is_free_index(gj))
            freeDofs(gj) = allDofs(j);
    }
}

/// Largest eigenvalue of the self-adjoint generalized EVP   Ax = lambda Bx
real_t lambdaMax(const gsSparseMatrix<>& A, const gsSparseMatrix<>& B)
{
    gsMatrix<> Adense = A;
    gsMatrix<> Bdense = B;
    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<>::Base > eig( Adense, Bdense );

    //gsInfo << Adense << "\n";
    //gsInfo << Bdense << "\n";

    return eig.eigenvalues()[A.rows() - 1];
}

/// TODO: documentation
void spectralTest(gsMultiGridOp<>& mg, const gsMatrix<real_t>& rhs, const gsSparseMatrix<>& M)
{
    if (mg.numLevels() != 2)
    {
        cerr << "Need two-level multigrid!" << "\n";
        exit(1);
    }

    gsMatrix<> Adense = mg.matrix();
    gsMatrix<> Mdense = M;
    gsInfo << setprecision(10) << Mdense << "\n" << "\n";
    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<>::Base > eig( Adense, Mdense );
    gsMatrix<> ev = eig.eigenvectors();
    gsMatrix<> eval = eig.eigenvalues();
    ev.colwise().normalize();

    gsMatrix<> x;

    gsInfo << "\n" << "smoother=[\n";
    for (index_t i = 0; i < ev.cols(); ++i)
    {
        x = ev.col( i );

        mg.smoothingStep(rhs, x);

        gsInfo << " " << eval(i) << "  " << x.norm() << "\n";

        //gsMatrix<> eigCoeffs = ev.transpose() * x;    // expand result in eigenvector basis

        //gsInfo << i << ": " << abs(eigCoeffs(i)) / eigCoeffs.lpNorm<1>() << "\n";
        //gsInfo << eigCoeffs.transpose() << "\n" << "\n";
    }

    gsInfo << "];" << "\n";

    gsInfo << "\ncgcorr=[\n";

    mg.setNumPreSmooth( 0 );
    mg.setNumPostSmooth( 0 );

    for (index_t i = 0; i < ev.cols(); ++i)
    {
        x = ev.col( i );
        mg.step( rhs, x );

        gsInfo << " " << eval(i) << "  " << x.norm() << "\n";

        //gsInfo << (ev.transpose() * x).transpose() << "\n" << "\n";

        //gsMatrix<> eigCoeffs = ev.transpose() * x;    // expand result in eigenvector basis

        //gsInfo << eigCoeffs.transpose() << "\n" << "\n";
        //gsInfo << i << ": " << abs(eigCoeffs(i)) / eigCoeffs.lpNorm<1>() << "\n";
    }
    gsInfo << "];" << "\n";

#if 0
    gsInfo << "\nkaczmarz=[\n";
    for (index_t i = 0; i < ev.cols(); ++i)
    {
        x = ev.col( i );
        for (int j = 0; j < 1; ++j)
            kaczmarzSweepBoundary<>(mg.matrix(), x, rhs, 0, 3);

        gsInfo << " " << eval(i) << "  " << x.norm() << "\n";
    }
    gsInfo << "];" << "\n";
#endif

    // print eigenvectors
    //for (index_t i = 0; i < ev.cols(); ++i)
    //{
    //    x = ev.col( i );
    //    gsInfo << x.transpose() << "\n";
    //}
}

/// Experimental smoother
gsLinearOperator<>::Ptr makeWeightedMassSmoother1D(const gsSparseMatrix<>& M, int p, real_t damping)
{
    const int N = M.rows();

    const real_t h = 1.0 / N;

    gsVector<> diag(N);
    diag.setConstant( 1.0 / ((p+1) * h));
    for (int i = 0; i < p; ++i)
    {
        diag[i] = 1.0 / ((i+2) * h);
        diag[N-1-i] = 1.0 / ((i+2) * h);
    }

    gsSparseMatrix<> sm = (1/damping) * diag.asDiagonal() * M * diag.asDiagonal();
    return makeSparseCholeskySolver(sm);
}

/// Experimental smoother
gsLinearOperator<>::Ptr makeWeightedMassSmoother2D(const gsSparseMatrix<>& M, int p, real_t damping)
{
    const int N2 = M.rows();
    const int N = (int)std::sqrt((double)N2);
    gsInfo << "N2 = " << N2 << "\nN =  " << N << "\n";
    assert(N2 == N*N);

    const real_t h = 1.0 / N;

    gsVector<> diag(N);
    diag.setConstant( 1.0 / ((p+1) * h));
    for (int i = 0; i < p; ++i)
    {
        diag[i] = 1.0 / ((i+2) * h);
        diag[N-1-i] = 1.0 / ((i+2) * h);
    }

    gsVector<> diag2(N2);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
        {
            diag2[i + N*j] = math::sqrt(diag[i] * diag[j]);
        }

    gsSparseMatrix<> sm = (1/damping) * diag2.asDiagonal() * M * diag2.asDiagonal();
    return makeSparseCholeskySolver(sm);
}

/// Experimental smoother
gsLinearOperator<>::Ptr makeWeightedMassSmoother3D(const gsSparseMatrix<>& M, int p, real_t damping)
{
    const int N3 = M.rows();
    // this is independent of the coefficient data type
    const int N = (int)(std::pow((double)N3, (1.0/3.0)) + 0.5);
    gsInfo << "N3 = " << N3 << "\nN =  " << N << "\n";
    assert(N3 == N*N*N);

    const real_t h = 1.0 / N;

    gsVector<> diag(N);
    diag.setConstant( 1.0 / ((p+1) * h));
    for (int i = 0; i < p; ++i)
    {
        diag[i] = 1.0 / ((i+1) * h);
        diag[N-1-i] = 1.0 / ((i+1) * h);
    }

    gsVector<> diag3(N3);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
            {
                diag3[i + N*j + N*N*k] = math::pow(diag[i] * diag[j] * diag[k], 1.0/3.0);
            }

    gsSparseMatrix<> sm = (1/(h*damping)) * diag3.asDiagonal() * M * diag3.asDiagonal();
    return makeSparseCholeskySolver(sm);
}

/// Allows to setup boundary conditions from command line
void bcChoose( char bc, condition_type::type& bc_type, gsFunction<> *& bc_func, gsFunction<> * bc_func_dirichlet, gsFunction<> * bc_func_neumann )
{
    if( bc == 'd' )
    {
        bc_type = condition_type::dirichlet;
        bc_func = bc_func_dirichlet;
    }
    else if( bc == 'n' )
    {
        bc_type = condition_type::neumann;
        bc_func = bc_func_neumann;
    }
    else
    {
        cerr << "Invalid boundary condition. Allowed are: dirichlet (d), neumann (n) and mixed (dd, dn, nd, nn, ddd, ddn, dnd, dnn, ndd, ndn, nnd, nnn).\n";
        exit(-1);
    }
}

void computeMatrixLevels(const std::vector< gsSparseMatrix<real_t,RowMajor> >& transfer, const gsSparseMatrix<real_t>& M, std::vector< gsSparseMatrix<real_t> >& Mlevels)
{
    const index_t sz = transfer.size() + 1;
    Mlevels.resize( sz );
    Mlevels[ sz-1 ] = M;

    for (index_t i = sz-2; i >= 0; --i)
    {
        Mlevels[i] = transfer[i].transpose() * Mlevels[i+1] * transfer[i];
    }
}


index_t maxMRSLevels( index_t p, index_t r )
{
    if(p<3)
        return r+1;
    else if(p<5)
        return r;
    else if(p<9)
        return r-1;
    else
        return r-2;
}

gsBoundaryConditions<real_t> updateBoundaryConditionsAfterSplit(const gsBoundaryConditions<real_t>& bcOld, const gsMultiPatch<real_t>& splittedPatch, index_t dim = -1)
{
    gsBoundaryConditions<> newBC;
    index_t d = (dim < 0) ? splittedPatch.parDim() : dim;
    for(gsBoundaryConditions<>::const_bciterator it = bcOld.beginAll();it!=bcOld.endAll();++it)
    {
        for(gsBoundaryConditions<>::const_iterator itBC = bcOld.begin(it->first);itBC!=bcOld.end(it->first);++itBC)
        {
            const boundary_condition<real_t>& bc= *itBC;
            for(gsMultiPatch<>::const_biterator iitMP = splittedPatch.bBegin();iitMP!=splittedPatch.bEnd();++iitMP)
            {
                const patchSide& boundary = *iitMP;
                if(boundary.patch >= math::exp2(d)*bc.patch() && boundary.patch < math::exp2(d)*(bc.patch()+1) && bc.side() == boundary.side())
                {
                    if(bc.type() == condition_type::unknownType)
                        newBC.add(boundary.patch,bc.side(),bc.ctype(),bc.function(),bc.unknown(),bc.parametric());
                    else
                        newBC.addCondition(boundary.patch,bc.side(),bc.type(),bc.function(),bc.unknown(),bc.parametric());
                }
            }
        }
    }
    return newBC;
}

template<class T>
gsMultiPatch<T> nonUniformSplit( const gsMultiPatch<T> & orig, index_t dir )
{
    const index_t n = 2;
    std::vector<gsGeometry<T>* > result;
    result.reserve(orig.nPatches()*n);

    for(size_t np = 0; np<orig.nPatches();++np)
    {
        std::vector<gsGeometry<T>* > result_temp = orig[np].uniformSplit(dir);
        GISMO_ENSURE( result_temp.size() == 2, "Internal error." );
        result.insert(result.end(),result_temp.begin(),result_temp.end());
    }
    gsMultiPatch<T> mp(result);
    mp.computeTopology();
    return mp;
}

class gsRepeatedParameterDomainSmoother : public gsPreconditionerOp<>
{
public:

    typedef memory::shared_ptr<gsRepeatedParameterDomainSmoother> Ptr;
    typedef memory::unique_ptr<gsRepeatedParameterDomainSmoother> uPtr;
    typedef gsPreconditionerOp<> Base;

    /// Constructor
    gsRepeatedParameterDomainSmoother(const memory::shared_ptr< gsSparseMatrix<> >& A,
                                const gsPreconditionerOp<>::Ptr sm,
                                const gsMultiBasis<>& mb,
                                const gsBoundaryConditions<>& bc,
                                const gsOptionList& opt,
                                index_t nrofsteps,
                                real_t damping = 1.0)
        : m_A(makeMatrixOp(A)), m_sm(sm), m_nrofsteps(nrofsteps), m_damping(damping)
    {

        const index_t nBases = mb.nBases();

        gsSparseMatrix<> Ahat(A->rows(),A->cols());

        gsDofMapper dm;
        mb.getMapper(
            (dirichlet::strategy)opt.askInt("DirichletStrategy",11),
            (iFace    ::strategy)opt.askInt("InterfaceStrategy", 1),
            bc,
            dm,
            0
        );

        for (index_t i=0; i<nBases; ++i)
        {
            const real_t h = mb[i].getMinCellLength();

            gsBoundaryConditions<> localbc;
            gsSparseMatrix<real_t> localstiff;
            assembleGeneralizedParameterStiffnessForTensorProductSpace(mb[i], localbc, (real_t)1, (real_t)(1/(h*h)), localstiff);


            const index_t nDofs = mb[i].size();
            gsSparseMatrix<real_t> transfer( nDofs, A->rows() );
            gsSparseEntries<real_t> tmp;
            tmp.reserve(nDofs);
            for (index_t j=0; j<nDofs; ++j)
            {
                const index_t dof_idx = dm.index(j,i);
                if (dm.is_free_index(dof_idx)) //TODO: check this!
                    tmp.add(j,dof_idx,1);
            }
            transfer.setFromTriplets(tmp.begin(),tmp.end());
            transfer.makeCompressed();
            Ahat += transfer.transpose() * localstiff * transfer;


        }
        m_Ahat = makeMatrixOp(Ahat.moveToPtr());


    }

    /// Make function called by smoother factory
    static uPtr make(const memory::shared_ptr< gsSparseMatrix<> >& A,
                     const gsPreconditionerOp<>::Ptr sm,
                     const gsMultiBasis<>& mb,
                     const gsBoundaryConditions<>& bc,
                     const gsOptionList& opt,
                     index_t nrofsteps,
                     real_t damping = 1.0)
    { return uPtr( new gsRepeatedParameterDomainSmoother( A, sm, mb, bc, opt, nrofsteps, damping ) ); }


    /// Make function called by smoother factory
    static uPtr make(const gsSparseMatrix<>& A,
                     const gsPreconditionerOp<>::Ptr sm,
                     const gsMultiBasis<>& mb,
                     const gsBoundaryConditions<>& bc,
                     const gsOptionList& opt,
                     index_t nrofsteps,
                     real_t damping = 1.0)
    { return uPtr( new gsRepeatedParameterDomainSmoother( memory::make_shared_not_owned(&A), sm, mb, bc, opt, nrofsteps, damping ) ); }

    void step(const gsMatrix<>& f, gsMatrix<>& x) const
    {
        m_A->apply(x,m_res0);
        m_res0 -= f;

        m_tmp.setZero(m_res0.rows(),m_res0.cols());
        m_sm->step(m_res0, m_tmp);
        m_update = - m_tmp;

        for (index_t i=1; i<m_nrofsteps; ++i)
        {
            m_Ahat->apply(m_update,m_tmp);
            m_res = m_res0 + m_tmp;

            m_tmp.setZero(m_res.rows(),m_res.cols());
            m_sm->step(m_res, m_tmp);
            m_update -= m_damping * m_tmp;
        }

        x += m_update;
    }

    index_t rows() const {return m_Ahat->rows();}
    index_t cols() const {return m_Ahat->cols();}

    gsLinearOperator<>::Ptr underlyingOp() const { return m_A; }

private:
    gsLinearOperator<>::Ptr m_A;
    gsLinearOperator<>::Ptr m_Ahat;
    gsPreconditionerOp<>::Ptr m_sm;
    mutable gsMatrix<real_t> m_res0, m_res, m_tmp, m_update;

    index_t m_nrofsteps;
    real_t m_damping;
};

/// @brief Null operator
///
/// \ingroup Solver
template<class T>
class gsNullOp : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsNullOp
    typedef memory::shared_ptr<gsNullOp> Ptr;

    /// Unique pointer for gsNullOp
    typedef memory::unique_ptr<gsNullOp> uPtr;

    /// Constructor taking the dimension of the identity operator
    gsNullOp(index_t dim) : m_dim(dim) {}

    /// Make function returning a smart pointer
    static uPtr make(index_t dim) { return memory::make_unique( new gsNullOp(dim) ); }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        x.setZero(input.rows(),input.cols());
    }

    index_t rows() const {return m_dim;}

    index_t cols() const {return m_dim;}

private:
    const index_t m_dim;
};



int main(int argc, char *argv[])
{
    string geometry("2");
    bool showBasis = false;
    string boundaryCondition("n");
    index_t numRefine = 3;
    index_t mult = 1;
    index_t degree = 2;
    index_t numLevels = -1;
    real_t numSplit = 0;
    Smoother::type smoother = Smoother::GaussSeidel;
    std::string smoother_name;
    index_t numPreSmooth = 1;
    index_t numPostSmooth = 1;
    index_t cycles = 1;
    real_t alpha = 1.;
    bool useCG = false;
    bool compEigs = false;
    bool useFMG = false;
    bool useCascadic = false;
    bool doSpectralTest = false;
    bool convergenceRate = true;
    bool monitorl2 = false, monitorL2 = false;
    bool writeLog = false;
    real_t tol = 1e-8;
    real_t damping = -1.0;
    real_t outerDamping = 1.;
    bool plot = false;
    index_t maxIter = 1000;
    bool tensorAssemble = false;
    std::string tmpDir("");

    index_t par_dom_rep = 1;
    real_t par_dom_rep_damp = 1;

    gsCmdLineWithEnumSupport cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addString("g", "geometry",           "Specification of the geometry (overrides dimension)",             geometry         );
    cmd.addSwitch(     "show-basis",         "Shows the chosen basis and quits.",                               showBasis        );
    cmd.addString("b", "boundary-condition", "Boundary condition",                                              boundaryCondition);
    cmd.addInt   ("r", "uniformRefine",      "Number of uniform h-refinement steps to perform before solving",  numRefine        );
    cmd.addInt   ("m", "multiplicity",       "Multiplicity of knots to insert when refining",                   mult             );
    cmd.addInt   ("p", "degree",             "Degree of the B-spline discretization space",                     degree           );
    cmd.addInt   ("l", "levels",             "Number of levels to use for multigrid iteration",                 numLevels        );
    cmd.addReal  ("",  "split",              "Split every patch uniformly into 2^d patches (default: 0)",       numSplit         );
    cmd.addInt   ("",  "presmooth",          "Number of pre-smoothing steps",                                   numPreSmooth     );
    cmd.addInt   ("",  "postsmooth",         "Number of post-smoothing steps",                                  numPostSmooth    );
    cmd.addEnum  ("s", "smoother",           "Smoothing method",                                                smoother         )
        .add(Smoother::Richardson,                                          "r",      "Richardson smoother"                                                       )
        .add(Smoother::Jacobi,                                              "j",      "Jacobi smoother"                                                           )
        .add(Smoother::GaussSeidel,                                         "gs",     "GaussSeidel smoother"                                                      )
        .add(Smoother::MassRichardson,                                      "mr",     "mass smoother"                                                             )
        .add(Smoother::MassRichardsonBoundaryCorrection,                    "mrb",    "mass smoother with boundary correction"                                    )
        .add(Smoother::MassRichardsonSubspaceCorrection,                    "mrs",    "mass smoother with subspace correction"                                    )
        .add(Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD,       "mrs-ad", "mass smoother with subspace correction based on additive Dirichlet dd"     )
        .add(Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD2,      "mrs-ad2","mass smoother with subspace correction based on additive Dirichlet dd; experimental"     )
        .add(Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDDRep,    "mrs-adr","mass smoother with subspace correction based on additive Dirichlet dd, p times par.-rep.")
        .add(Smoother::MassRichardsonSubspaceCorrectionAdditiveMPNDD,       "mrs-an", "mass smoother with subspace correction based on additive Neumann dd"       )
        .add(Smoother::MassRichardsonSubspaceCorrectionMultiplicativeMPNDD, "mrs-mn", "mass smoother with subspace correction based on multiplicative Neumann dd" )
        .add(Smoother::MassRichardsonSubspaceCorrectionGeo,                 "mrs-g",  "mass smoother with subspace correction for rank 1 geometry approximation"  )
        .add(Smoother::MassRichardsonSubspaceCorrectionGS,                  "mrs-gs", "mass smoother with subspace correction with GaussSeidel"                   )
        .writeDescOfChosenOptTo(smoother_name);
    cmd.addReal  ("",  "damping",            "Damping factor for the smoother (handed over to smoother)",       damping          );
    cmd.addReal  ("",  "outerdamping",       "Damping factor for the smoother (globally)",                      outerDamping     );
    cmd.addInt   ("c", "cycles",             "Number of multi-grid cycles",                                     cycles           );
    cmd.addInt   ("",  "maxiter",            "Maximum number of iterations",                                    maxIter          );
    cmd.addReal  ("a", "alpha",              "alpha in \"- LAPLACE u + alpha u = f\"",                          alpha            );
    cmd.addSwitch(     "log",                "Write results to log file",                                       writeLog         );
    cmd.addSwitch(     "cg",                 "Use CG iteration",                                                useCG            );
    cmd.addSwitch(     "eigs",               "Compute eigenvalues of the preconditioned system. Requires --cg.",compEigs         );
    cmd.addSwitch(     "fmg",                "Use full multigrid cycle",                                        useFMG           );
    cmd.addSwitch(     "cascadic",           "Use cascadic multigrid",                                          useCascadic      );
    cmd.addSwitch(     "tensor",             "Assemble using tensor product (experimental)",                    tensorAssemble   );
    cmd.addReal  ("",  "tol",                "Tolerance for multigrid solver stopping criterion",               tol              );
    cmd.addSwitch(     "spectral-test",      "Perform a numerical spectral test",                               doSpectralTest   );
    cmd.addSwitch(     "noConvergenceRate",  "Do not print convergence rate",                                   convergenceRate  );
    cmd.addSwitch(     "monitor-L2",         "Monitor the L2 errors over the iteration",                        monitorL2        );
    cmd.addSwitch(     "monitor-l2",         "Monitor the discrete l2 errors over the iteration",               monitorl2        );
    cmd.addSwitch(     "plot",               "Plot result in ParaView format",                                  plot             );

    cmd.addString("",  "tmp",  "Write matrices to a temp directory (default: no, otherwise specify path)",      tmpDir           );

    cmd.addInt   ("",  "par_dom_rep",        "Needs -s mrs-adr",                                                par_dom_rep      );
    cmd.addReal  ("",  "par_dom_rep_damp",   "Needs -s mrs-adr",                                                par_dom_rep_damp );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    /******************** Define Geometry ********************/

    gsStopwatch time;
    gsMultiPatch<> mp;

    std::string orig_geometry = geometry;

    if (geometry=="1")
    {
        gsGeometry<>::uPtr geo = gsNurbsCreator<>::BSplineUnitInterval(static_cast<short_t>(1));

        // offset the computational domain by 0.5 in each direction (for historical reasons)
        gsVector<> offs = 0.5 * gsVector<>::Ones( geo->parDim() ); geo->translate( offs );

        mp = gsMultiPatch<>(*geo);
    }
    else if (geometry=="2")
    {
        gsGeometry<>::uPtr geo = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(1));

        // offset the computational domain by 0.5 in each direction (for historical reasons)
        gsVector<> offs = 0.5 * gsVector<>::Ones( geo->parDim() ); geo->translate( offs );

        mp = gsMultiPatch<>(*geo);
    }
    else if (geometry=="3")
    {
        gsGeometry<>::uPtr geo = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(1));

        // offset the computational domain by 0.5 in each direction (for historical reasons)
        gsVector<> offs = 0.5 * gsVector<>::Ones( geo->parDim() ); geo->translate( offs );

        mp = gsMultiPatch<>(*geo);
    }
    else if (geometry=="4")
    {
        gsGeometry<>::uPtr geo = approximateQuarterAnnulus(static_cast<short_t>(2));

        // offset the computational domain by 0.5 in each direction (for historical reasons)
        gsVector<> offs = 0.5 * gsVector<>::Ones( geo->parDim() ); geo->translate( offs );

        mp = gsMultiPatch<>(*geo);
    }
    else if (geometry=="5")
    {
        gsGeometry<>::uPtr geo = BSplineMySquare(static_cast<short_t>(1));

        // offset the computational domain by 0.5 in each direction (for historical reasons)
        gsVector<> offs = 0.5 * gsVector<>::Ones( geo->parDim() ); geo->translate( offs );

        mp = gsMultiPatch<>(*geo);
    }
    else
    {
        if (geometry=="6") geometry = "volumes/twistedFlatQuarterAnnulus.xml";
        if (geometry=="7") geometry = "yeti_mp2.xml";
        if (geometry=="8") geometry = "lshape_3_patches.xml";
        if (geometry=="9") geometry = "volumes/fichera_u7p.xml";

        if ( gsFileManager::fileExists(geometry) )
        {
            geometry = gsFileManager::find(geometry);
            gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
            if (!mpPtr) { cerr << "No multipatch object found in file " << geometry << ".\n"; return 1; }
            mp = *mpPtr;
        }
        else
        {
            cerr << "Invalid geometry. Allowed are:\n"
                << "1: gsNurbsCreator<>::BSplineUnitInterval(1)\n"
                << "2: gsNurbsCreator<>::BSplineSquareDeg(1)\n"
                << "3: gsNurbsCreator<>::BSplineCube(1)\n"
                << "4: approximateQuarterAnnulus(2)\n"
                << "5: BSplineMySquare(1) // unit square with one additional refinement in x-direction\n"
                << "6: volumes/twistedFlatQuarterAnnulus.xml\n"
                << "7: yeti_mp2.xml\n"
                << "8: lshape_3_patches.xml\n"
                << "9: volumes/fichera_u7p.xml\n"
                << "or a valid filename.\n";
            return -1;
        }
    }

    gsInfo << "The geometry consists of " << mp.nPatches() << " patches.\n";

    if ( showBasis )
    {
        gsInfo << mp << endl;
        return 0;
    }

    gsFunction<>::Ptr f0;
    gsFunction<>::Ptr g;

    switch (mp.geoDim())
    {
        case 1:
            f0 = memory::make_shared( new gsFunctionExpr<>("(pi^2 ) * sin(pi*x)",1) );
            g = memory::make_shared( new gsFunctionExpr<>("sin(pi*x)",1) );
            break;
        case 2:
            //f0 = memory::make_shared( new gsFunctionExpr<>("(2*pi^2 ) * sin(pi*x) * sin(pi*y)",2) );
            //g = memory::make_shared( new gsFunctionExpr<>("sin(pi*x) * sin(pi*y)",2) );
            f0 = memory::make_shared( new gsFunctionExpr<>("2*25*pi^2*sin(5*pi*x) * sin(5*pi*y)",2) );
            g = memory::make_shared( new gsFunctionExpr<>("0",2) );
            break;
        case 3:
           // f0 = memory::make_shared( new gsFunctionExpr<>("(3*pi^2 ) * sin(pi*x) * sin(pi*y) * sin(pi*z)",3) );
           // g = memory::make_shared( new gsFunctionExpr<>("sin(pi*x) * sin(pi*y) * sin(pi*z)",3) );
            f0 = memory::make_shared( new gsFunctionExpr<>("3*25*pi^2*sin(5*pi*x) * sin(5*pi*y) * sin(5*pi*z)",3) );
            g = memory::make_shared( new gsFunctionExpr<>("0",3) );
            break;
        default:
            cerr << "Invalid geometry dimension.\n";
            return -1;
    }

    gsFunction<>::Ptr f = gsLinearCombinationOfFunctionsFunction<>::make(1,f0,alpha,g);

    if (numRefine < 1)
    {
        cerr << "Number of refinements must be positive.\n"; return -1;
    }
    if (mult < 1)
    {
        cerr << "Multiplicity must be positive.\n"; return -1;
    }
    if (numLevels < 1)
    {
        numLevels = numRefine + 1;
        if( smoother == Smoother::MassRichardsonSubspaceCorrection
            || smoother == Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD
            || smoother == Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD2
            || smoother == Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDDRep
            || smoother == Smoother::MassRichardsonSubspaceCorrectionAdditiveMPNDD
            || smoother == Smoother::MassRichardsonSubspaceCorrectionMultiplicativeMPNDD
            || smoother == Smoother::MassRichardsonSubspaceCorrectionGeo
            || smoother == Smoother::MassRichardsonSubspaceCorrectionGS
        )
            numLevels = maxMRSLevels(degree,numRefine);
        gsInfo << "The number of levels was chosen to be " << numLevels << ".\n";
    }
    if (numLevels < 1)
    {
        cerr << "Number of levels must be positive.\n"; return -1;
    }
    if (numRefine - numLevels + 1 < 0)
    {
        cerr << "Not enough refinements for the desired number of levels.\n"; return -1;
    }
    if (cycles < 1)
    {
        cerr << "Number of cycles must be positive.\n"; return -1;
    }
    if ( useFMG && useCascadic )
    {
        cerr << "Cannot combine full multigrid and cascadic.\n"; return -1;
    }
    if ( useCG && ( useFMG || useCascadic ) )
    {
        cerr << "Cannot use CG for any full multigrid.\n"; return -1;
    }
    if (compEigs && !useCG)
    {
        cerr << "Cannot compute eigenvalues for the preconditioned system without applying CG.\n"; return -1;
    }
    if (damping < 0)
    {
        // Here, one could add appropriate parameter-choice rules (depending on the degree)
        switch (smoother)
        {
            case Smoother::Richardson:                                     damping = 0.80; break;
            case Smoother::Jacobi:                                         damping = 0.80; break;
            case Smoother::GaussSeidel:                                    break;
            case Smoother::MassRichardson:                                 damping = 0.25 / ( (degree+1.)*(degree+1.) ); break;
            case Smoother::MassRichardsonBoundaryCorrection:               damping = 0.09; break;
            case Smoother::MassRichardsonSubspaceCorrection:               damping = 0.09; break;
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD:       damping = 0.09; break;
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:      damping = 0.09; break;
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDDRep: damping = 0.09; break;
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPNDD:       damping = 0.09; break;
            case Smoother::MassRichardsonSubspaceCorrectionMultiplicativeMPNDD: damping = 0.09; break;
            case Smoother::MassRichardsonSubspaceCorrectionGeo:            damping = 1.00; break;
            case Smoother::MassRichardsonSubspaceCorrectionGS:             damping = 0.09; break;
        }
    }
    if(!convergenceRate && !useCG)
    {
        cerr << "Preventing calculation of convergence rate can only be performed with CG.\n"; return -1;
    }

    // ---------------------------------------------------------------------------------------------------------------------------

    std::string fullFn;
    bool foundSavedMatrix = false;
    if (!tmpDir.empty())
    {
        std::replace(geometry.begin(), geometry.end(), '/', '_');
        fullFn = tmpDir + "/mg_matrix_"
                        + geometry + "_"
                        + util::to_string(numRefine) + "_"
                        + util::to_string(degree) + "_"
                        + util::to_string(mult) + "_"
                        + util::to_string(alpha) + "_"
                        + util::to_string(tensorAssemble) + "_"
                        + boundaryCondition + ".xml";
        if (gsFileManager::fileExists(fullFn))
            foundSavedMatrix = true;
    }

    // ---------------------------------------------------------------------------------------------------------------------------

    if (doSpectralTest)
    {
        // ensure two-grid method with zero exact solution
        numLevels = 2;
        f0 = memory::make_shared( new gsConstantFunction<>(0.0, mp.geoDim()) );
        g = memory::make_shared( new gsConstantFunction<>(0.0, mp.geoDim()) );
    }

    gsInfo << "Source function: " << *f << ".\n" << "\n";
    gsInfo << "Exact solution:  " << *g << ".\n" << "\n";


    // set up boundary conditions

    gsConstantFunction<> zero(0.0, mp.geoDim());
    gsConstantFunction<> one (1.0, mp.geoDim());

    gsBoundaryConditions<> bc;
    condition_type::type bc_type;
    gsFunction<> * bc_func;

    if (mp.nPatches()==1)
    {
        // if only single BC given, use it in all coordinate directions
        if (boundaryCondition.length() == 1)
            boundaryCondition = string(mp.geoDim(), boundaryCondition[0]);

        if( (index_t)boundaryCondition.length() != mp.geoDim() )
            boundaryCondition = "x"; // Let the bcChoose do the work

        bcChoose( boundaryCondition[0], bc_type, bc_func, &*g, &zero );
        bc.addCondition( boundary::west,  bc_type, bc_func );
        bc.addCondition( boundary::east,  bc_type, bc_func );
        if (mp.geoDim() >= 2)
        {
            bcChoose( boundaryCondition[1], bc_type, bc_func, &*g, &zero );
            bc.addCondition( boundary::south, bc_type, bc_func );
            bc.addCondition( boundary::north, bc_type, bc_func );
        }
        if (mp.geoDim() >= 3)
        {
            bcChoose( boundaryCondition[2], bc_type, bc_func, &*g, &zero );
            bc.addCondition( boundary::front, bc_type, bc_func );
            bc.addCondition( boundary::back,  bc_type, bc_func );
        }
    }
    else if (boundaryCondition == "dn" && orig_geometry=="7")
    {
        std::vector<patchSide> outer_bdy = getConnectedBoundaryComponent( mp, *mp.bBegin() ); // this is rather hackisch, we just guess that
                                                                                              // the first is on the outer boundary
        index_t d_nr = 0, n_nr = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
        {
            if ( std::find( outer_bdy.begin(), outer_bdy.end(), *it ) != outer_bdy.end() )
            {
                d_nr++;
                bcChoose( 'd', bc_type, bc_func, &*g, &zero );
            }
            else
            {
                n_nr++;
                bcChoose( 'n', bc_type, bc_func, &*g, &zero );
            }
            bc.addCondition( *it, bc_type, bc_func );
        }
        gsInfo << "Added " << d_nr << " Dirichlet and " << n_nr << " Neumann boundary condtions.\n";

    }
    else if (boundaryCondition == "dn" && orig_geometry=="9")
    {
        index_t d_nr = 0, n_nr = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
        {
            gsMatrix<> midpoint = mp.pointOn(*it);
            if ( midpoint(0,0) > .999 || midpoint(1,0) > .999 || midpoint(2,0) < -.999 )
            {
                d_nr++;
                bcChoose( 'd', bc_type, bc_func, &*g, &zero );
            }
            else
            {
                n_nr++;
                bcChoose( 'n', bc_type, bc_func, &*g, &zero );
            }
            bc.addCondition( *it, bc_type, bc_func );
        }
        gsInfo << "Added " << d_nr << " Dirichlet and " << n_nr << " Neumann boundary condtions.\n";
    }
    else
    {
        if (boundaryCondition.length() != 1)
        {
            gsWarn << "Only one boundary condition is acceptable for multipatch domains.\n";
            return -1;
        }

        index_t i = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
        {
            ++i;
            bcChoose( boundaryCondition[0], bc_type, bc_func, &*g, &zero );
            bc.addCondition( *it, bc_type, bc_func );
        }
        gsInfo << "Added " << i << " boundary conditions.\n";
    }

    for (index_t i=0; i<(index_t)numSplit; ++i)
    {
        gsInfo << "Split multipatch object uniformly... " << flush;
        mp = mp.uniformSplit();

        bc = updateBoundaryConditionsAfterSplit(bc,mp);
        gsInfo << "done." << endl;
    }

    index_t partial_splits = ( numSplit - (index_t)numSplit ) * mp.geoDim() + 0.1;

    for (index_t i=0; i<partial_splits; ++i)
    {
        gsInfo << "Split multipatch object one-sided... " << flush;
        mp = nonUniformSplit(mp,i);

        bc = updateBoundaryConditionsAfterSplit(bc,mp,1);
        gsInfo << "done." << endl;
    }

    bc.print(gsInfo);

    if( numSplit )
        gsInfo << "The geometry consists of " << mp.nPatches() << " patches.\n";


    gsMultiBasis<> mb(mp);

    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    vector< gsMultiBasis<> > bases;
    vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;
    vector< vector< gsSparseMatrix<real_t, RowMajor> > > fullTransferMatrices;

    gsOptionList assemblerOptions = gsAssembler<>::defaultOptions();

    // Define coarse discretization space by refining the basis of the geometry
    for (int i = 0; i < numRefine - numLevels + 1; ++i)
        mb.uniformRefine();       // refine until coarsest level

    if (mult>1)
    {
        gsInfo << "Reduce continuity by " << mult-1 << std::endl;
        mb.reduceContinuity(mult-1);
    }

    const double timeSetupGeo = time.stop(); time.restart();

    gsInfo << "Setup grid hierarchy..." << std::flush;

    // set up the hierarchy of spaces and transfer matrices between them
    const index_t refineKnots = 1;
    gsGridHierarchy<>::buildByRefinement(give(mb), bc, assemblerOptions, numLevels, refineKnots, mult)
        .moveMultiBasesTo(bases)
        .moveTransferMatricesTo(transferMatrices);

    const double timeSetupGH = time.stop(); time.restart();

    gsInfo << "done." << std::endl;

    if ( useFMG || useCascadic )
    {
        // compute full transfer matrices per hand, starting from coarsest level...
        fullTransferMatrices.resize(transferMatrices.size());
        gsMultiBasis<real_t> tmp = bases[0];
        for ( index_t i=1; i<numLevels; ++i )
        {
            fullTransferMatrices[i-1].resize(tmp.nBases());
            for ( size_t j=0; j<tmp.nBases(); ++j )
                tmp[j].uniformRefine_withTransfer(fullTransferMatrices[i-1][j], refineKnots, mult);
        }
    }

    /////////// TEST
    // gsBSplineBasis<> * bsp = dynamic_cast< gsBSplineBasis<>* >(&mb.basis(0));
    // assert(bsp);
    // gsSparseMatrix<> P_tilde, P_compl;
    // tildeSpaceBasis(*bsp, P_tilde, P_compl);
    // cout << P_tilde << endl << P_compl << endl;
    /////////// END TEST

    gsInfo << (useFMG ? "Full multigrid" : (useCascadic ? "Cascadic multigrid" : (useCG ? "CG preconditioned by multigrid" : "Multigrid")));
    gsInfo << " with " << numLevels << " levels using " << smoother_name;
    gsInfo << "(" << numPreSmooth << "," << numPostSmooth << ") smoother" << " and " << cycles << "-cycle.\n";
    for (size_t i=0; i<mb.nBases(); ++i)
        gsInfo << "Coarse discretization space: dim=" << mb[i].dim() << " deg=" << mb[i].degree(0) << " dofs=" << mb[i].size() << "\n";

    gsInfo << "Setup gsGeneralizedPoissonAssembler... " << flush;
    //gsGenericAssembler<real_t> genassm(mp, bases.back(), assemblerOptions, &bc);
    gsGeneralizedPoissonAssembler<real_t> assm(mp, bases.back(), bc, *f, alpha, assemblerOptions);

    const double timeSetupAssembler = time.stop(); time.restart();

    gsInfo << "done." << endl;

    gsSparseMatrix<real_t> Kfine;
    gsMatrix<real_t> rhs;

    if ( ! foundSavedMatrix )
    {
        if ( ! tensorAssemble )
        {
            gsInfo << "Assembling stiffness matrix... " << flush;
            assm.assemble();
            Kfine = assm.matrix();
            rhs = assm.rhs();
            gsInfo << "done, " << Kfine.rows() << " dofs and " << Kfine.nonZeros() << " non-zeros." << endl;
        }
        else
        {
            if ( !(geometry == "1" || geometry == "2" || geometry == "3") )
            {
                gsInfo << "Tensorized assembling is only possible if there is no geometry transformation.\n";
                gsInfo << "This is not the case for geometry = \"" << geometry << "\"\n.";
                return -1;
            }
            if (monitorL2 || plot)
            {
                gsInfo << "These options require the assembler object.\n";
                return -1;
            }
            if ( boundaryCondition != string(mp.geoDim(), 'n') )
            {
                gsInfo << "**********\n* WARNING: The tensor assembler does not treat the inhomogenous Dirichlet boundary conditions properly.\n**********\n";
            }
            gsInfo << "Assemble in a tensorized way... " << flush;
            assembleGeneralizedParameterStiffnessForTensorProductSpace<>(bases.back()[0], bc, (real_t)1., alpha, Kfine);
            assembleParameterMomentsForTensorProduct(bases.back()[0], bc, *f, rhs);
            gsInfo << "done, " << Kfine.rows() << " dofs." << endl;
        }

        if (!fullFn.empty())
        {
            gsInfo << "Write data to file " << fullFn << "..." << std::flush;
            gsFileData<> fd;
            fd << Kfine;
            fd << rhs;
            fd.save(fullFn);
            gsInfo << "done." << std::endl;
        }

    }
    else
    {
        gsInfo << "Found file " << fullFn << ". Read data just from there..." << std::flush;
        gsFileData<> fd(fullFn);
        fd.getFirst(Kfine);
        fd.getFirst(rhs);
        gsInfo << "done, " << Kfine.rows() << " dofs." << endl;
    }

    gsInfo<<"Number of NNZ in K: "<<Kfine.nonZeros()<<"\n";
    //gsInfo << "Assembling mass matrix... " << flush;
    //genassm.assembleMass();
    //gsSparseMatrix<real_t> Mfine = genassm.fullMatrix();
    //gsInfo << "done." << endl;
    //
    //gsInfo << "Setup the problem matrix \"K + " << alpha << " * M\"... " << flush;
    //Kfine += alpha * Mfine;
    //gsInfo << "done." << endl;
    //
    //gsInfo << "Assembling right-hand side... " << flush;
    //genassm.assembleMoments(f);
    //gsInfo << "done." << endl;

    /*************************************************************************/

    const double timeAssembling = time.stop(); time.restart();

    vector< gsMatrix<real_t> > dirichletIntp;
    vector< gsMatrix<real_t> > rhsForAllLevels;

    if ( useFMG || useCascadic )
    {
        gsInfo << "Setup full multi grid hierarchy... " << flush;
        const int lvls = bases.size();
        std::vector< gsDofMapper > dofMappers( lvls );
        rhsForAllLevels.resize( lvls );
        dirichletIntp.resize( lvls );

        rhsForAllLevels[lvls-1] = rhs;

        for (int i = 0; i < lvls; ++i)
        {
            bases[i].getMapper(
                (dirichlet::strategy)assemblerOptions.getInt("DirichletStrategy"),
                (iFace::strategy)assemblerOptions.getInt("InterfaceStrategy"),
                bc, dofMappers[i], 0 );
        }

        for ( int i = lvls-2; i >= 0; --i )
        {
            // We need to recompute the right-hand sides for each level;
            // the easiest way is just to assemble the whole system.
            //gsGenericAssembler<real_t> genassmLevel(mp, bases[i], assemblerOptions, &bc);
            gsGeneralizedPoissonAssembler<real_t> assmLevel(mp, bases[i], bc, *f, alpha, assemblerOptions);
            assmLevel.assemble();
            rhsForAllLevels[i] = assmLevel.rhs();

            // compute prolongation contributions from the Dirichlet boundary conditions
            gsMatrix<real_t> zeroVec = gsMatrix<real_t>::Zero( rhsForAllLevels[i].size(), 1 );
            gsMatrix<real_t> dirValues;
            completeVector( dofMappers[i], zeroVec, assmLevel.fixedDofs(), dirValues );
            dirValues = fullTransferMatrices[i][0] * dirValues;
            extractFreeDofs( dofMappers[i+1], dirValues, dirichletIntp[i+1] );

        }
        dirichletIntp[0].resize(rhsForAllLevels[0].size(),1);
        dirichletIntp[0].setZero();
        gsInfo << "done.\nnorms(dirichletIntp) = [  ";
        for ( int i = 0; i < lvls; ++i )
        {
            gsInfo << dirichletIntp[i].norm() << "  ";
        }
        gsInfo << "]\nnorms(rhs[i]-P rhs[i+1]) = [ ";
        for ( int i = 0; i < lvls-1; ++i )
        {
            gsInfo << ( rhsForAllLevels[i] - transferMatrices[i].transpose() * rhsForAllLevels[i+1] ).norm() << "  ";
        }
        gsInfo << "]\n" << flush;
    }

    // set up the multigrid solver
    gsMultiGridOp<> mg(give(Kfine), transferMatrices);
    //gsMultiGridOp<> mg(give(Kfine), transferMatrices, gsNullOp<real_t>::make(transferMatrices.back().rows()) );

    // Determine coarse solve time:
    gsTimedOp<>::Ptr coarseSolver = gsTimedOp<>::make("CoarseSolver",mg.coarseSolver(),false);
    mg.setCoarseSolver(coarseSolver);

    vector< gsSparseMatrix<real_t> > massMatrices;
    if ( doSpectralTest )
    {
        gsInfo << "Computing mass matrices... " << flush;
        gsGenericAssembler<real_t> genassm(mp, bases.back(), assemblerOptions, &bc);
        genassm.assembleMass();
        gsSparseMatrix<real_t> Mfine = genassm.fullMatrix();
        computeMatrixLevels(transferMatrices, Mfine, massMatrices);
        gsInfo << "done." << "\n";
    }

    mg.setNumPreSmooth( numPreSmooth );
    mg.setNumPostSmooth( numPostSmooth );
    mg.setNumCycles( cycles );

    GISMO_ASSERT( outerDamping == 1 || smoother != Smoother::GaussSeidel, "Gauss-Seidel does not support --outerdamping" );

    const double timeSetupMGObject = time.stop(); time.restart();
    double timeLocalSmoothers = 0;

    gsInfo << "Constructing smoothers... " << flush;
    for (int i = mg.numLevels() == 1 ? 0 : 1; i < mg.numLevels(); ++i)
    {
        switch( smoother ) {
            case Smoother::Richardson:                            mg.setSmoother(i, makeRichardsonOp(mg.matrix(i),damping*outerDamping)); break;
            case Smoother::Jacobi:                                mg.setSmoother(i, makeJacobiOp(mg.matrix(i),damping*outerDamping)); break;
            case Smoother::GaussSeidel:                           mg.setSmoother(i, makeGaussSeidelOp(mg.matrix(i))); break;
            case Smoother::MassRichardson:                        mg.setSmoother(i, gsPreconditionerFromOp<>::make(mg.underlyingOp(i),makeMassSmootherOperator(bases[i][0], damping, bc), outerDamping )); break;
            case Smoother::MassRichardsonBoundaryCorrection:      mg.setSmoother(i, gsPreconditionerFromOp<>::make(mg.underlyingOp(i),makeBoundaryCorrectedMassSmootherOperator(bases[i][0], damping, bc), outerDamping )); break;
            case Smoother::MassRichardsonSubspaceCorrection:      mg.setSmoother(i, gsPreconditionerFromOp<>::make(mg.underlyingOp(i),makeSubspaceCorrectedMassSmootherOperator(bases[i][0], damping, bc), outerDamping )); break;
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD:
            {

                gsStopwatch timer2;
                std::vector<gsLinearOperator<>::Ptr> localSmoothers = makeSubspaceCorrectedMassSmootherOperatorsDirichlet(bases[i],damping);
                timeLocalSmoothers += timer2.stop();

                mg.setSmoother(
                    i,
                    gsAdditiveSmoother::make(
                        mg.underlyingOp(i),
                        setupPiecewisePreconditioner(
                            mg.matrix(i),
                            give(localSmoothers),
                            bases[i],
                            bc,
                            assemblerOptions
                        ),
                        outerDamping
                    )
                );
                break;
            }
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:
            {
                const index_t dim = bases[i].dim();

                std::vector< std::vector< std::pair< typename gsBasis<>::Ptr, gsSparseMatrix<> > > > pieces =
                    constructPieces( bases[i], bc, assemblerOptions );

                std::vector< gsLinearOperator<>::Ptr > localSmoothers;
                std::vector< gsSparseMatrix<> > smootherTransfers;
                const index_t nrPieces = pieces.size();


                real_t h = 1;
                for ( size_t j=0; j<bases[i].nBases(); ++j)
                    h = std::min(h,bases[i][j].getMinCellLength());

                for ( index_t dd = 0; dd<nrPieces; ++dd )
                {
                    gsBoundaryConditions<> localbc;
                    for( index_t ps=0; ps < 2*dim; ++ps )
                        localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

                    const index_t sz = pieces[dd].size();
                    for ( index_t j=0; j<sz; ++j)
                    {
                        gsStopwatch timer2;
                        if (dd == dim)
                        {
                            localSmoothers.push_back( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][j].first, damping, localbc) );
                        }
                        else if (dd > 1)
                        {
                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                            localSmoothers.push_back( gsScaledOp<>::make( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][j].first, damping, localbc, scalingM/scalingK ), 1/scalingK ) );
                        }
                        else if (dd == 1)
                        {
                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg.matrix(i) * pieces[dd][j].second.transpose();
                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );

                            const real_t scalingK = std::pow(2,dim-dd) * (dim-dd) * std::pow( h/(2*degree+1), dim-dd );
                            const real_t scalingM = std::pow(2,dim-dd) * std::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                            gsSparseMatrix<> M, K;
                            assembleParameterMass(*pieces[dd][j].first, M);
                            assembleParameterStiffness(*pieces[dd][j].first, K);
                            M = M.block( 1, 1, M.rows()-2, M.cols()-2 );
                            K = K.block( 1, 1, K.rows()-2, K.cols()-2 );
                            gsSparseMatrix<> mat = scalingK * K + scalingM * M;

                            localSmoothers.push_back( makeSparseCholeskySolver(mat) );
                        }
                        else if (dd == 0)
                        {
                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg.matrix(i) * pieces[dd][j].second.transpose();
                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );
                            const index_t sz = pieces[dd][j].second.rows();
                            const real_t scalingFactor = std::pow(2,dim) * dim * std::pow( h/(1+2*degree), dim-1 ) * degree*degree/(h*(2*degree-1));
                            localSmoothers.push_back( gsScaledOp<>::make( gsIdentityOp<>::make(sz), 1/scalingFactor) );
                        }
                        timeLocalSmoothers += timer2.stop();

                        smootherTransfers.push_back( give( pieces[dd][j].second ) );
                    }
                }
                gsInfo << "[" << localSmoothers.size() << "] ";

                mg.setSmoother(
                    i,
                    gsAdditiveSmoother::make(
                        mg.underlyingOp(i),
                        give(smootherTransfers),
                        give(localSmoothers),
                        outerDamping
                    )
                );
                break;
            }
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDDRep:
            {
                mg.setSmoother(
                    i,
                    gsRepeatedParameterDomainSmoother::make(
                        mg.matrix(i),
                        gsAdditiveSmoother::make(
                            mg.underlyingOp(i),
                            setupPiecewisePreconditioner(
                                mg.matrix(i),
                                makeSubspaceCorrectedMassSmootherOperatorsDirichlet(bases[i],damping),
                                bases[i],
                                bc,
                                assemblerOptions
                            ),
                            outerDamping
                        ),
                        bases[i],
                        bc,
                        assemblerOptions,
                        par_dom_rep,
                        par_dom_rep_damp
                    )
                );
                break;
            }
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPNDD:
            {
                mg.setSmoother(
                    i,
                    gsAdditiveSmoother::make(
                        mg.underlyingOp(i),
                        getPatchwiseTransfers(bases[i], bc, assemblerOptions),
                        makeSubspaceCorrectedMassSmootherOperators(bases[i],damping,bc),
                        outerDamping
                    )
                );
                break;
            }
            case Smoother::MassRichardsonSubspaceCorrectionMultiplicativeMPNDD:
            {
                mg.setSmoother(
                    i,
                    gsMultiplicativeSmoother::make(
                        mg.underlyingOp(i),
                        getPatchwiseTransfers(bases[i], bc, assemblerOptions),
                        makeSubspaceCorrectedMassSmootherOperators(bases[i],damping,bc),
                        outerDamping
                    )
                );
                break;
            }
            case Smoother::MassRichardsonSubspaceCorrectionGeo:   mg.setSmoother(i, gsPreconditionerFromOp<>::make(mg.underlyingOp(i),makeSubspaceCorrectedMassSmootherOperator(mp.patch(0), bases[i][0], damping, bc), outerDamping) ); break;
            case Smoother::MassRichardsonSubspaceCorrectionGS:    mg.setSmoother(i,
                    gsCompositePrecOp<>::make(
                        gsAdditiveSmoother::make(
                            mg.underlyingOp(i),
                            setupPiecewisePreconditioner(
                                mg.matrix(i),
                                makeSubspaceCorrectedMassSmootherOperatorsDirichlet(bases[i],damping),
                                bases[i],
                                bc,
                                assemblerOptions
                            ),
                            outerDamping
                        ),
                        //gsPreconditionerFromOp<>::make( makeMatrixOp( mg.matrix(i) ), makeSubspaceCorrectedMassSmootherOperator(bases[i][0], damping, bc) ),
                        makeGaussSeidelOp( mg.matrix(i) )
                    )
                ); break;
        }
    }
    gsInfo << "done." << "\n";

    if (doSpectralTest)
    {
        spectralTest( mg, rhs, massMatrices[ mg.finestLevel() ] );
        return 0;
    }

    gsMatrix<> x;
    //x.setRandom( mg.nDofs(), 1 );
    x.setZero( mg.nDofs(), 1 );

    const double timeSetupSmoother = time.stop(); time.restart();

    const real_t resNorm0 = (rhs - mg.matrix() * x).norm();
    gsInfo << "Residual norm:     " << resNorm0 << "\n";
    real_t resNorm(0), oldResNorm = resNorm0;

    gsField<> sol;

    real_t l2Err, oldL2Err(0), eucl_error, old_eucl_error(0);
    if (monitorL2)
    {
        sol = assm.constructSolution(x);
        oldL2Err = computeL2Distance(sol, *g, false, 3*mg.nDofs());
        gsInfo << "L2 error: " << oldL2Err << "\n";
    }

    int numIter = 0;
    real_t minReduction = 1e6;

    gsMatrix<> exactDiscreteSol;
    if (monitorl2)
    {
        Eigen::SparseLU< gsSparseMatrix<real_t> > directsolver( mg.matrix() );
        exactDiscreteSol = directsolver.solve( rhs );
        old_eucl_error = (exactDiscreteSol - x).norm();
        gsInfo << "Euclidean error: " << old_eucl_error << "\n";
    }

    gsConjugateGradient<> cg( mg.underlyingOp(), memory::make_shared_not_owned(&mg) );
    cg.setTolerance( tol );

    if (compEigs)
        cg.setCalcEigenvalues(true);

    if (useCG && convergenceRate)
        cg.initIteration( rhs, x );


    if (useFMG) // solve using full multi-grid cycle
    {
        mg.fullMultiGrid(rhsForAllLevels,dirichletIntp,x);
        gsInfo << "Residual norm after full multigrid cycle: " <<
            (rhs - mg.matrix() * x).norm() << std::endl;
    }
    else if (useCascadic)
    {
        mg.cascadicMultiGrid(rhsForAllLevels,dirichletIntp,x);
        gsInfo << "Residual norm after cascadic multigrid cycle: " <<
            (rhs - mg.matrix() * x).norm() << std::endl;
    }
    else if (convergenceRate)       // solve using MG iteration
    {
        do
        {
            if (useCG)
                cg.step(x);
            else if (numLevels == 1)
                mg.smoothingStep(rhs, x);
            else
                mg.step(rhs, x);

            resNorm = (rhs - mg.matrix() * x).norm();
            gsInfo << "Residual norm:     " << left << setw(15) << resNorm << "          reduction:  1 / " << setprecision(3) << (oldResNorm/resNorm) << setprecision(6) << "\n";
            minReduction = math::min(minReduction, oldResNorm/resNorm);
            oldResNorm = resNorm;

            if (monitorL2)
            {
                sol = assm.constructSolution(x);
                l2Err = computeL2Distance(sol, *g, false, 3*mg.nDofs());
                gsInfo << "                                                                   |  L2 error:     "
                     << left << setw(15) << l2Err << "          reduction:  1 / " << setprecision(3) << (oldL2Err/l2Err) << setprecision(6) << "\n";
                oldL2Err = l2Err;
            }

            if (monitorl2)
            {
                eucl_error = (exactDiscreteSol - x).norm();
                gsInfo << "                                                                   |  l2 error:     "
                     << left << setw(15) << eucl_error << "          reduction:  1 / " << setprecision(3) << (old_eucl_error/eucl_error) << setprecision(6) << "\n";
                old_eucl_error = eucl_error;
            }

            ++numIter;
        } while (resNorm / resNorm0 > tol && numIter < maxIter && gsIsfinite(resNorm));
    }
    else
    {
        cg.solve(rhs, x);
        numIter = cg.iterations();

    }

    const double timeSolve = time.stop();

    if (!useFMG && !useCascadic) {
        if (convergenceRate && (resNorm / resNorm0 > tol || ! gsIsfinite(resNorm)))
            gsInfo << "Did not converge.\n";
        else
            gsInfo << "Converged in " << numIter << " iterations.\n";

        if(convergenceRate)
        {
            gsInfo << "Average convergence factor:  1 / " << setprecision(3) << math::pow(resNorm0 / resNorm, 1.0 / numIter) << setprecision(6) << "\n";
            gsInfo << "Worst   convergence factor:  1 / " << setprecision(3) << minReduction << setprecision(6) << "\n";
            gsInfo << "\n";
        }
    }

    const double timeSetup = timeSetupGH+timeSetupMGObject+timeSetupSmoother+timeSetupAssembler;

    gsInfo << "Setup Geometry:  "; formatTime(gsInfo, timeSetupGeo);                 gsInfo << "\n";
    gsInfo << "Assembling time: "; formatTime(gsInfo, timeAssembling);               gsInfo << "\n";
    gsInfo << "MG time setup:   "; formatTime(gsInfo, timeSetup);                    gsInfo << "\n";
    gsInfo << " grid hierarchy: "; formatTime(gsInfo, timeSetupGH);                  gsInfo << "\n";
    gsInfo << " MG object:      "; formatTime(gsInfo, timeSetupMGObject);            gsInfo << "\n";
    gsInfo << " smoother:       "; formatTime(gsInfo, timeSetupSmoother);            gsInfo << "\n";
    gsInfo << " assembler:      "; formatTime(gsInfo, timeSetupAssembler);           gsInfo << "\n";
    if (timeLocalSmoothers>0)
    {
        gsInfo << "   local smooth: "; formatTime(gsInfo, timeLocalSmoothers);           gsInfo << "\n";
    }
    gsInfo << "MG time solving: "; formatTime(gsInfo, timeSolve);
    if (!useFMG && !useCascadic)
    {
        gsInfo << "         (avg. "; formatTime(gsInfo, timeSolve/numIter);          gsInfo << " per iteration)" << "\n";
        gsInfo << " coarse solver:  "; formatTime(gsInfo, coarseSolver->getTime());
    }
    gsInfo << "\n";
    gsInfo << "Total time:      "; formatTime(gsInfo, timeSetupGeo+timeAssembling+timeSetup+timeSolve);  gsInfo << "\n";
    gsInfo << "\n";

    if (monitorL2 || plot)
        sol = assm.constructSolution(x);

    if (monitorL2)
    {
        // Compute L2 error
        const real_t error = computeL2Distance(sol, *g, false, 3*mg.nDofs());
        gsInfo << "L2 error: " << error << "\n";
    }

    if (monitorl2)
    {
        const gsVector<> f_err = exactDiscreteSol - x;
        const real_t eucl_error_f = f_err.norm();
        gsInfo << "l2 error: " << eucl_error_f << "\n";
        gsInfo << "Discrete energy error: " << f_err.dot(mg.matrix()*f_err)  << "\n";
    }

    real_t max = 0., min = 1.e12, essMin = 1.e12;
    if (compEigs)
    {
        //gsInfo << "Condition number: " << cg.getConditionNumber() << std::endl;
        gsMatrix<real_t> eigs;
        cg.getEigenvalues(eigs);
        //gsInfo << "Eigenvalues: " << eigs.transpose() << std::endl;
        for ( index_t i = 0; i < eigs.rows(); ++i )
        {
            if (eigs(i) > max)
                max = eigs(i);
            if (eigs(i) < min)
                min = eigs(i);
            if (eigs(i) < essMin && eigs(i) > 1.e-5)
                essMin = eigs(i);
        }
        gsInfo << "Eigenvalues: Minimum, essentialMinimum, maximum, essentialConditionNumber: " << min << ", " << essMin << ", " << max << ", " << (max/essMin) << std::endl;
    }

    if (writeLog)
    {
        fstream log("out.txt", fstream::out | fstream::app);

        log << "multigrid_example" << (useCG ? "_cg" : "" ) << "\t"
            << geometry << "\t"
            << alpha << "\t"
            << cycles << "\t"
            << degree << "\t"
            << numRefine << "\t"
            << numLevels << "\t"
            << numSplit << "\t"
            << boundaryCondition << "\t"
            << damping << "\t"
            << outerDamping << "\t"
            << (int)smoother << "\t"
            << numPreSmooth << "\t"
            << numPostSmooth << "\t"
            << numIter << "\t"
            << math::pow(resNorm0 / resNorm, 1.0 / numIter) << "\t"
            << timeAssembling << "\t"
            << timeSetup << "\t"
            << timeSolve;

            if (compEigs)
                log << "\t" << min << "\t" << essMin << "\t" << max << "\t" << (max/essMin) ;

            log << "\n";
    }

    if (plot)
    {
        // Plot solution in Paraview
        gsInfo << "Plotting in Paraview: multigrid.pvd.\n";
        gsWriteParaview<>(sol, "multigrid", 1000);
        gsFileManager::open("multigrid.pvd");
    }

    return (resNorm / resNorm0 > tol) ? 1 : 0;
};
