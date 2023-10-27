/** @file gsStokesMultiGridPreconditioner.cpp

    @brief Provides test examples for multigrid preconditioner for Stokes

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsRecipeAssembler/gsRecipeAssemblerStokes.h>
#include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsRankOneAssembler.h>
#include <gsMultiGrid/gsMassSmoother.h>
#include <gsAssembler/gsGeneralizedPoissonAssembler.h>
#include <gsSolver/gsBramblePasciakCG.h>

#include <gsIO/gsCmdLineWithEnumSupport.h>


using namespace gismo;

namespace Smoother {
    enum type {
        Richardson = 0,
        Jacobi = 1,
        GaussSeidel = 2,
        MassRichardsonSubspaceCorrection = 3,
        MassRichardsonSubspaceCorrectionGeo =4
    };
}

namespace PressurePreconditioner {
    enum type {
        pressure_mass = 0,
        pressure_jacobi = 1,
        pressure_sgs = 2,
        pressure_mass_geo = 3
    };
}

namespace StokesDiscrtizationType {
    enum type {
        TaylorHood = 0,
        Nedelec = 1,
        RaviartThomas = 2
    };
}

namespace MultigridSpace {
    enum type {
        param = 0,
        partly = 1,
        phys = 2
    };
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

template <typename MatrixT>
real_t diff(const MatrixT &a,const MatrixT &b)
{
    if (a.rows()==b.rows() && a.cols()==b.cols() )
        return (a-b).norm();
    else
        return std::numeric_limits<real_t>::quiet_NaN();
}


std::vector<gsPhysicalSpace*> constructSpaces(
        const std::vector< gsBasis<>* >& bases,
        const gsMultiPatch<>& geo,
        ValueTransformationType valueTransformationType,
        StokesDiscrtizationType::type discretizationType,
        std::vector< std::vector< gsBasis<>* > >* outBasis
        )
{
    const index_t dim = bases[0]->domainDim();
    std::vector< gsPhysicalSpace* > result(2);
    

    outBasis->clear();
    
    {   // build velocity space
        
        std::vector< gsPhysicalSpaceScalar* > scalarSpaces;
        for (index_t i=0; i<dim; ++i)
        {

            std::vector< gsBasis<>* > patchBasisV;
            cloneAll( bases, patchBasisV );

            for (unsigned j = 0; j < patchBasisV.size(); ++j )
            {
                if( discretizationType == StokesDiscrtizationType::TaylorHood ) {
                    patchBasisV[j]->degreeElevate(); //degreeElevate preserves smoothness, but not the knot multiplicty
                } else if( discretizationType == StokesDiscrtizationType::Nedelec ) {
                    for (index_t d=0; d<dim; ++d)
                    {
                        if (d==i)
                            patchBasisV[j]->degreeIncrease(d); // degreeIncrease preserves knot multiplicty                                
                        else
                            patchBasisV[j]->degreeElevate(d);  // degreeElevate preserves smoothness
                    }
                } else if( discretizationType == StokesDiscrtizationType::RaviartThomas ) {
                    patchBasisV[j]->degreeIncrease(i);         // degreeIncrease preserves knot multiplicty                                
                } else {
                    GISMO_ERROR( "Unknwon discretization type." );
                }
            }

            std::vector< gsBasis<>* > temp;
            cloneAll( patchBasisV, temp );
            
            gsMultiBasis<real_t> multipatch( temp, geo );
            gsMapFactoryMultiBasis mapFactory( multipatch );

            std::vector< gsBasis<>* > argVec = patchBasisV;

            scalarSpaces.push_back( new gsPhysicalSpaceScalar( argVec, geo, valueTransformationType, *memory::make_unique(mapFactory.makeMapper()) ) );
            outBasis->push_back(patchBasisV);
            
        }
        
        result[velocity] = ( new gsPhysicalSpaceVector(scalarSpaces) );
        
        freeAll(scalarSpaces);
        
    }
    
    {   // build pressure space
        
        std::vector< gsBasis<real_t>* > temp;
        cloneAll( bases, temp );
        
        gsMultiBasis<real_t> multipatch( temp, geo );
        gsMapFactoryMultiBasis mapFactory( multipatch );

        std::vector< gsBasis<>* > patchBasesP;
        cloneAll( bases, patchBasesP );
        
        std::vector< gsBasis<real_t>* > argVec = patchBasesP;

        result[pressure] = new gsPhysicalSpaceScalar( argVec, geo, INVERSE_COMPOSITION, *memory::make_unique(mapFactory.makeMapper()) );
        outBasis->push_back(patchBasesP);
    }
    
    return result;
}


template <typename T>
T computeL2DistanceUpToConst(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, const gsFunction<T>& v, bool isParametrized_v, int numEvals)
{
    GISMO_ASSERT( u.targetDim() == v.targetDim(), "Functions need to have same target dimension");

    const int d = geo.parDim();
    assert( d == geo.geoDim() );

    // compute the tensor Gauss rule
    gsMatrix<T> nodes;
    gsVector<T> weights;
    gsMatrix<T> range = geo.basis().support();

    // Number of nodes of the underlying Gauss rule to use
    const index_t nodesPerInterval = 1;
    const int nodesPerElement  = math::ipow(nodesPerInterval, d);
    const int numElements      = (numEvals + nodesPerElement - 1) / nodesPerElement;
    std::vector< std::vector<T> > intervals;
    uniformIntervals<T>(range.col(0), range.col(1), intervals, numElements);

    // perform the quadrature
    gsGaussRule<T> QuRule( gsVector<index_t>::Constant(d,nodesPerInterval) );
    const int numPts = QuRule.numNodes();

    gsTensorDomainIterator<T> domIt(intervals);
    gsMatrix<T> geo_pts, geo_jac, u_val, v_val;
    T u2 = 0.0;
    T u1 = 0.0;
    T u0 = 0.0;

    for (; domIt.good(); domIt.next() )
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), nodes, weights );

        // only compute the geometry points if either function is not parametrized
        geo_pts =  (!isParametrized_u || !isParametrized_v) ?
                    geo.eval(nodes) : gsMatrix<T>();
        geo_jac = geo.jacobian(nodes);
        
        // evaluate u and v
        u_val = isParametrized_u ? u.eval(nodes) : u.eval(geo_pts);
        v_val = isParametrized_v ? v.eval(nodes) : v.eval(geo_pts);
        
        for (index_t k = 0; k < numPts; ++k)
        {
            const T funcDet = math::abs( geo_jac.block(0, k*d, d,d).determinant() ) ;
            const gsVector<T> diff = u_val.col(k) - v_val.col(k);
                  gsVector<T> ones = diff; ones.setOnes();
            
            u2 += weights[k] * funcDet * diff.dot(diff);
            u1 += weights[k] * funcDet * diff.dot(ones);
            u0 += weights[k] * funcDet * ones.dot(ones);
        }
    }

    return math::sqrt(u2-u1*u1/u0);
}

template <typename T>
T computeL2DistanceUpToConst(const gsField<T>& u, const gsFunction<T>& v, bool isParametrized_v, int numEvals)
{
    T dist = T();

    for (index_t i = 0; i < u.nPatches(); ++i)
    {
        T curDist = computeL2DistanceUpToConst( u.patch(i), u.function(i), u.isParametrized(), v, isParametrized_v, numEvals);
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}

template<typename T=real_t>
class gsPressureControlAdvanced : public gsLinearOperator<T>
{
public:
    typedef memory::shared_ptr<gsPressureControlAdvanced> Ptr;
    typedef memory::unique_ptr<gsPressureControlAdvanced> uPtr;

    /// Constructor
    gsPressureControlAdvanced(index_t sz1, index_t sz2, const gsSparseMatrix<T>& mass)
    : m_sz1(sz1), m_sz2(sz2)
    {
        m_ones.setZero( sz2 );
        m_ones.array() += 1;
        m_Mones = mass * m_ones;

        T factor = m_ones.dot(m_Mones);

        m_Mones.array() /= factor;

    }

    static uPtr make(index_t sz1, index_t sz2, const gsSparseMatrix<T>& mass)
    { return uPtr( new gsPressureControlAdvanced( sz1, sz2, mass ) ); }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        x = input;
        const T factor = x.bottomRows(m_sz2).col(0).dot(m_Mones.col(0));
        x.bottomRows(m_sz2) -= factor * m_ones;
    }

    index_t rows() const {
        return m_sz1 + m_sz2;
    }

    index_t cols() const {
        return m_sz1 + m_sz2;
    }

private:
    index_t m_sz1; 
    index_t m_sz2;        
    gsVector<T> m_ones;
    gsVector<T> m_Mones;
};

template<typename T=real_t>
class gsPressureRankOne : public gsLinearOperator<T>
{
public:
    typedef memory::shared_ptr<gsPressureRankOne> Ptr;
    typedef memory::unique_ptr<gsPressureRankOne> uPtr;
    
    /// Constructor
    gsPressureRankOne(T tau, index_t sz1, index_t sz2, const gsSparseMatrix<T>& mass)
    : m_tau(tau), m_sz1(sz1), m_sz2(sz2)
    {
        gsVector<T> tmp;
        tmp.setZero( sz2 );
        tmp.array() += 1;
        
        m_ones = mass * tmp;
        
        T factor = sqrt(m_ones.dot(m_ones));
        
        m_ones.array() /= factor;
        
    }
    
    static uPtr make(T tau, index_t sz1, index_t sz2, const gsSparseMatrix<T>& mass)
    { return uPtr( new gsPressureRankOne( tau, sz1, sz2, mass ) ); }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.cols() == 1, "This only works for a single right-hand side.");
        x.setZero(input.rows(),1);
        const T factor = m_tau * x.bottomRows(m_sz2).col(0).dot(m_ones.col(0));
        x.bottomRows(m_sz2) -= factor * m_ones;
    }

    index_t rows() const {
        return m_sz1 + m_sz2;
        
    }

    index_t cols() const {
        return m_sz1 + m_sz2;
    }

private:
    T m_tau;
    index_t m_sz1; 
    index_t m_sz2;        
    gsVector<T> m_ones;
};

template<class T=real_t>
class gsPreconditionedRichardsonOp : public gsPreconditionerOp<T>
{
public:

    /// Shared pointer for gsPreconditionedRichardsonOp
    typedef memory::shared_ptr<gsPreconditionedRichardsonOp> Ptr;

    /// Unique pointer for gsPreconditionedRichardsonOp
    typedef memory::unique_ptr<gsPreconditionedRichardsonOp> uPtr;

    /// Base class
    typedef gsLinearOperator<T> Base;
    
    /// Base class pointer
    typedef typename Base::Ptr BasePtr;
    
    gsPreconditionedRichardsonOp( const BasePtr& A, const BasePtr& precon )
     : m_A(A), m_precon(precon), m_damping(1)
    {
        GISMO_ASSERT( A->rows() == precon->cols() && A->cols() == precon->rows(), "Dimensions do not agree." );
         
    }
     
    static uPtr make( const BasePtr& A, const BasePtr& precon )
    { return uPtr( new gsPreconditionedRichardsonOp( A, precon) ); }
    
    virtual ~gsPreconditionedRichardsonOp() {}

    virtual void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        GISMO_ASSERT( m_A->rows() == x.rows() && x.rows() == rhs.rows() && m_A->cols() == m_A->rows() && x.cols() == rhs.cols(),
            "The dimensions do not fit." );

        m_A->apply( x, m_res );
        m_res -= rhs;
        m_precon->apply( m_res, m_corr );
        x -= m_damping * m_corr;
    }    
    
    virtual void setDamping( T damping ) { m_damping = damping;     }
    virtual T damping()                  { return m_damping;        }
    virtual index_t rows() const         { return m_precon->rows(); }
    virtual index_t cols() const         { return m_precon->cols(); }
    virtual BasePtr underlyingOp() const { return m_A;              }

protected:
    BasePtr m_A;
    BasePtr m_precon;
    T m_damping;
    mutable gsMatrix<T> m_res;
    mutable gsMatrix<T> m_corr;
    
}; // gsPreconditionedRichardsonOp

template<class T=real_t>
class gsProductOfSteppableOperatorsOp : public gsPreconditionerOp<T>
{
public:

    /// Shared pointer for gsProductOfSteppableOperatorsOp
    typedef memory::shared_ptr<gsProductOfSteppableOperatorsOp> Ptr;

    /// Unique pointer for gsProductOfSteppableOperatorsOp
    typedef memory::unique_ptr<gsProductOfSteppableOperatorsOp> uPtr;

    /// Base class
    typedef gsPreconditionerOp<T> Base;
    
    /// Base class pointer
    typedef typename Base::Ptr BasePtr;
    
    gsProductOfSteppableOperatorsOp( const std::vector<BasePtr>& op )
     : m_op(op)
    {
        //GISMO_ASSERT( A->rows() == precon->cols() && A->cols() == precon->rows(), "Dimensions do not agree." );
    }

    gsProductOfSteppableOperatorsOp( const BasePtr& A, const BasePtr& B )
     : m_op(2)
    {
        //GISMO_ASSERT( A->rows() == precon->cols() && A->cols() == precon->rows(), "Dimensions do not agree." );
        m_op[0] = A;
        m_op[1] = B;
    }

    gsProductOfSteppableOperatorsOp( const BasePtr& A, const BasePtr& B, const BasePtr& C )
     : m_op(3)
    {
        //GISMO_ASSERT( A->rows() == precon->cols() && A->cols() == precon->rows(), "Dimensions do not agree." );
        m_op[0] = A;
        m_op[1] = B;
        m_op[2] = C;
    }
     
    static uPtr make( const std::vector<BasePtr>& op )
    { return uPtr( new gsProductOfSteppableOperatorsOp( op ) ); }

    static uPtr make( const BasePtr& A, const BasePtr& B )
    { return uPtr( new gsProductOfSteppableOperatorsOp( A, B ) ); }

    static uPtr make( const BasePtr& A, const BasePtr& B, const BasePtr& C )
    { return uPtr( new gsProductOfSteppableOperatorsOp( A, B, C ) ); }
    
    virtual ~gsProductOfSteppableOperatorsOp() {}

    virtual void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        //GISMO_ASSERT( m_A->rows() == x.rows() && x.rows() == rhs.rows() && m_A->cols() == m_A->rows() && x.cols() == rhs.cols(),
        //    "The dimensions do not fit." );

        for (size_t i = 0; i<m_op.size(); ++i)
            m_op[i]->step(rhs,x);
    }

    virtual index_t rows() const         { return m_op[0]->rows(); }
    virtual index_t cols() const         { return m_op[0]->cols(); }
    virtual typename gsLinearOperator<T>::Ptr underlyingOp() const { return m_op[0]->underlyingOp(); }

protected:
    std::vector<BasePtr> m_op;
}; // gsProductOfSteppableOperatorsOp

void addAsBlock(const gsSparseMatrix<>& mat, index_t row, index_t col, gsSparseMatrix<>& result)
{
    GISMO_ASSERT( row + mat.rows() <= result.rows(), "Dimensions to not agree" );
    GISMO_ASSERT( col + mat.cols() <= result.cols(), "Dimensions to not agree" );
    //TODO: this is not really good in terms of efficiency
    for (index_t k=0; k < mat.outerSize(); ++k)
        for (gsSparseMatrix<>::InnerIterator it(mat,k); it; ++it)
            result(row+it.row(),col+it.col()) += it.value();
}

void addAsBlock(const gsSparseMatrix<real_t,RowMajor>& mat, index_t row, index_t col, gsSparseMatrix<real_t,RowMajor>& result)
{
    GISMO_ASSERT( row + mat.rows() <= result.rows(), "Dimensions to not agree" );
    GISMO_ASSERT( col + mat.cols() <= result.cols(), "Dimensions to not agree" );
    //TODO: this is not really good in terms of efficiency
    for (index_t k=0; k < mat.outerSize(); ++k)
        for (gsSparseMatrix<real_t,RowMajor>::InnerIterator it(mat,k); it; ++it)
            result(row+it.row(),col+it.col()) += it.value();
}

std::string smootherName( Smoother::type smoother )
{
    switch (smoother)
    {
        case Smoother::Richardson: return "a Richardson smoother";
        case Smoother::Jacobi: return "a Jacobi smoother";
        case Smoother::GaussSeidel: return "a Gauss-Seidel smoother";
        case Smoother::MassRichardsonSubspaceCorrection: return "a mass-Richardson smoother with subspace correction";
        case Smoother::MassRichardsonSubspaceCorrectionGeo: return "a mass-Richardson smoother with subspace correction (with rank-1-geometry approximation)";
    }
    return "an unknwon smoother";
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


int main(int argc, char *argv[])
{

    index_t geoIndex = 1;
    index_t numRefine = 3;
    index_t numLevels = -1;
    index_t degree = 2;
    index_t numPreSmooth = 1;
    index_t numPostSmooth = 1;
    index_t sweeps = 1;
    index_t outerSweeps = 1;
    Smoother::type smoother = Smoother::MassRichardsonSubspaceCorrection;
    real_t damping_v = -1;
    real_t damping_p = 1;
    real_t outerdamping = 1;
    real_t outerdamping_v = 1;
    index_t cycles = 1;
    index_t maxIter = 500;
    bool directSolve = false;
    StokesDiscrtizationType::type discretizationType = StokesDiscrtizationType::TaylorHood;
    std::string discretizationTypeString;
    bool divConforming = false;
    real_t tol = 1.e-6;
    index_t gaussSeidelSweeps = 0;
    PressurePreconditioner::type pprecond = PressurePreconditioner::pressure_mass;
    MultigridSpace::type mg_phys = MultigridSpace::param;
    std::string tmpDir("");
    bool computeError = false;
    bool writeLog = false;

    gsCmdLineWithEnumSupport cmd("Solves the Stokes problem with an isogeometric discretization using a minres solver, preconditioned with multigrid.");
    cmd.addInt("g", "geometry", "Specification of the geometry", geoIndex);
    cmd.addInt("r", "refinements", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("l", "levels", "Number of levels to use for multigrid iteration (defaulted to the value of refinements [r])", numLevels);
    cmd.addInt("p", "degree", "Degree of the B-spline discretization space", degree);
    cmd.addInt("", "presmooth", "Number of pre-smoothing steps", numPreSmooth);
    cmd.addInt("", "postsmooth", "Number of post-smoothing steps", numPostSmooth);
    cmd.addInt("", "sweeps", "Number of multigrid cycles used as preconditioner", sweeps);
    cmd.addInt("", "outersweeps", "Number of outer sweeps for mutligrid on physical domain", outerSweeps);
    cmd.addEnum("s", "smoother", "Smoothing method", smoother)
        .add(Smoother::Richardson,                          "r",     "Richardson"                          )
        .add(Smoother::Jacobi,                              "j",     "Jacobi"                              )
        .add(Smoother::GaussSeidel,                         "gs",    "GaussSeidel"                         )
        .add(Smoother::MassRichardsonSubspaceCorrection,    "mrs",   "MassRichardsonSubspaceCorrection"    )
        .add(Smoother::MassRichardsonSubspaceCorrectionGeo, "mrs-g", "MassRichardsonSubspaceCorrectionGeo" );
    cmd.addReal("", "damping_v", "Damping parameter for the smoother for the velocity v", damping_v);
    cmd.addReal("", "outerdamping_v", "Damping parameter for the smoother for the velocity v", outerdamping_v);
    cmd.addReal("", "damping_p", "Damping parameter for the preconditioner for the pressure p", damping_p);
    cmd.addReal("", "outerdamping", "Outer damping (for outersweeps)", outerdamping);
    cmd.addInt("c", "cycles", "Number of multi-grid cycles", cycles);
    cmd.addInt("", "maxiter", "Maximum number of iterations", maxIter);
    cmd.addSwitch("direct", "Use direct solver", directSolve);
    cmd.addInt("","gs", "Number of GaussSeidel sweeps", gaussSeidelSweeps);
    cmd.addEnum("d", "disc", "Discretization type", discretizationType)
        .add(StokesDiscrtizationType::TaylorHood,       "TH",       "TaylorHood"     )
        .add(StokesDiscrtizationType::Nedelec,          "N",        "Nedelec"        )
        .add(StokesDiscrtizationType::RaviartThomas,    "RT",       "RaviartThomas"  )
        .writeKeyOfChosenOptTo(discretizationTypeString);
    cmd.addSwitch("div", "Use Piola (=div conforming) transformation (if not set: use inverse composition)", divConforming);
    cmd.addReal("", "tol", "Tolerance for multigrid solver stopping criterion", tol);
    cmd.addEnum("","pprecond", "Preconditioner for pressure",pprecond)
        .add(PressurePreconditioner::pressure_mass,     "mass",     "the exact inverse of the mass matrix on the parameter domain"                          )
        .add(PressurePreconditioner::pressure_mass_geo, "mass-geo", "the exact inverse of a rank one approximation of mass matrix on the physical domain"   )
        .add(PressurePreconditioner::pressure_jacobi,   "jacobi",   "the Jacobi preconditioner for the mass matrix on the physical domain"                  )
        .add(PressurePreconditioner::pressure_sgs,      "sgs",      "the Symmetric Gauss Seidel preconditioner for the mass matrix on the physical domain"  );
    cmd.addEnum("","mg_phys", "Where to setup the multigrid preconditioner", mg_phys)
        .add(MultigridSpace::param,      "param",  "Parameter domain"                                      )
        .add(MultigridSpace::partly,     "partly", "On the physical domain, ignoring the Piola transform"  )
        .add(MultigridSpace::phys,       "phys",   "On the physical domain, considering a Piola transform" );
    cmd.addString("", "tmp", "Write matrices to a temp directory (default: no, otherwise specify path)", tmpDir);
    cmd.addSwitch("compute_error", "Compute the final error", computeError);
    cmd.addSwitch("log", "Write to logfile out.txt", writeLog);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (numRefine < 0)
    {
        std::cerr << "Number of refinements must be non-negative.\n"; return -1;
    }

    if (numLevels < 1)
    {
        numLevels = numRefine;
        if( smoother == Smoother::MassRichardsonSubspaceCorrection || smoother == Smoother::MassRichardsonSubspaceCorrectionGeo )
            numLevels = maxMRSLevels(degree,numRefine) - 1;
        gsInfo << "The number of levels was chosen to be " << numLevels << "\n";
    }
 
    if (numRefine - numLevels < 0)
    {
        std::cerr << "Not enough refinements for the desired number of levels.\n"; return -1;
    }

    if (damping_v < 0)
    {
        // Here, one could add appropriate parameter-choice rules (depending on the degree)
        switch (smoother)
        {
            case Smoother::Richardson:                                     damping_v = 0.80; break;
            case Smoother::Jacobi:                                         damping_v = 0.80; break;
            case Smoother::GaussSeidel:                                    break;
            case Smoother::MassRichardsonSubspaceCorrection:               damping_v = 0.04; break;
            case Smoother::MassRichardsonSubspaceCorrectionGeo:            damping_v = 1.00; break;
        }
    }


    // -------------------------------------------------------------------------------------------------------------------------------------

    gsInfo << "Stokes multigrid example.\n";
    gsOptionList opt = cmd.getOptionList();
    gsInfo << opt << "\n\n\n";

    // -------------------------------------------------------------------------------------------------------------------------------------

    gsInfo << "Setup problem specification..." << std::flush;
    gsMultiPatch<> domain;
    switch (geoIndex)
    {
        case 1: domain = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(1)));              break;
        case 2: domain = gsMultiPatch<>(*approximateQuarterAnnulus(2));                    break;
        case 3: domain = gsMultiPatch<>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus());    break;
        default: gsInfo << "Invalid geometry. Allowed are:\n"
            << "1: unit square\n"
            << "2: approximateQuarterAnnulus\n"
            << "3: BSplineFatQuarterAnnulus\n";
        return -1;
    }
    const index_t dim = 2;

    gsInfo << "Chosen domain # " << geoIndex << ":" << std::endl;
    gsInfo << domain << std::endl;

    ValueTransformationType valueTransformationType = divConforming ? DIV_CONFORMING : INVERSE_COMPOSITION;

    gsFunctionExpr<>       f("+2*25*cos(5*x+5*y)+2*25*sin(5*x-5*y)+(1+y)","-2*25*cos(5*x+5*y)+2*25*sin(5*x-5*y)+(1+x)",2);
    gsFunctionExpr<>       u("cos(5*x+5*y)+sin(5*x-5*y)","-1-cos(5*x+5*y)+sin(5*x-5*y)",2);
    gsFunctionExpr<> p;
    if(geoIndex==1)
        p = gsFunctionExpr<>("-(1+x)*(1+y)+9/4",2);
    else if(geoIndex==2)
        p = gsFunctionExpr<>("-(1+x)*(1+y)+8257/2100",2);
    else if(geoIndex==3)
        p = gsFunctionExpr<>("-(1+x)*(1+y)+8257/2100",2);

    gsMultiBasis<>         mb(domain);
    gsStokesPde<real_t>    pde(domain,gsBoundaryConditions<>(),&f);

    dirichlet::strategy    Dstrategy(dirichlet::elimination);  

    std::vector<patchSide> boundaries = domain.boundaries();
    for(unsigned i = 0;i<boundaries.size();++i)
        pde.boundaryConditions().addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::dirichlet, &u);

    bool eliminated = Dstrategy==dirichlet::elimination && pde.boundaryConditions().dirichletSides().size()>0;
    gsInfo << "done." << std::endl;

    // -------------------------------------------------------------------------------------------------------------------------------------

    gsInfo << "The geometry has " << domain.dim() << " dimensions.\n";
    gsInfo << "The geometry has " << domain.nPatches() << " patches.\n";
    GISMO_ENSURE( domain.nPatches() == 1, "Multigrid works only for one patch" );

    // -------------------------------------------------------------------------------------------------------------------------------------

    gsBasis<> * tbasis = &mb.basis(0);

    tbasis->setDegreePreservingMultiplicity(degree);

    gsInfo << "Refine grid " << numRefine << " times..." << std::flush;
    for (int i = 0; i < numRefine; ++i)
        tbasis->uniformRefine();
    gsInfo << "done." << std::endl;


    gsInfo << "Setup Stokes bases..." << std::flush;

    std::vector< std::vector< gsBasis<>* > > bases;
    std::vector<gsPhysicalSpace*> phySpace;
    {
        std::vector<gsBasis<>*> tmpBasisVec;
        tmpBasisVec.push_back(tbasis);
        phySpace = constructSpaces(tmpBasisVec,pde.domain(),valueTransformationType,discretizationType,&bases);
    }
    gsInfo << "done." << std::endl;

    gsInfo << "Got bases:" << bases.size() << std::endl;
    for (unsigned i = 0; i < bases.size(); ++i)
        for (unsigned j = 0; j < bases[j].size(); ++j)
        {
            gsInfo << "Basis(" << i << "," << j << ") has " << bases[i][j]->numElements() << " elements, " 
                << bases[i][j]->size() << " dofs and degree " << bases[i][j]->minDegree() << " to " << bases[i][j]->maxDegree() << ".\n";
        }

    GISMO_ASSERT( bases.size() == dim+1, "Stokes is assumed to have velocity and pressure." );

    // -------------------------------------------------------------------------------------------------------------------------------------

    std::string fullFn;
    bool foundSavedMatrix = false;
    if (!tmpDir.empty())
    {
        fullFn = tmpDir + "/stokes_matrix_"
                        + util::to_string(geoIndex) + "_"
                        + util::to_string(numRefine) + "_"
                        + util::to_string(degree) + "_"
                        + discretizationTypeString + "_"
                        + util::to_string(divConforming) + ".xml";
        if (gsFileManager::fileExists(fullFn) && ! directSolve && ! computeError)
            foundSavedMatrix = true;
    }

    gsRecipeAssemblerStokes assembler(pde);
    gsSparseMatrix<> sys;
    gsMatrix<> rhs;
    gsMatrix<> eli;
    if (!foundSavedMatrix)
    {
        gsInfo << "Setup Stokes assembler..." << std::flush;
        assembler.setSpace(phySpace);
        assembler.setZeroAverage(false);
        gsInfo << "done." << std::endl;

        gsInfo << "Assemble..." << std::flush;
        assembler.assemble();
        sys = assembler.getSystemMatrix();
        rhs = assembler.getSystemRhs();
        gsInfo << "done. System matrix has " << sys.rows() << " rows and " << sys.cols() << " columns." << std::endl;

        eli.resize(0,rhs.cols());
        if (eliminated)
        {
            gsInfo << "Eliminate Dirichlet dofs..." << std::flush;
            gsSparseSolver<>::QR  solver;
            solver.analyzePattern( assembler.getEliminatedMatrix() );
            solver.factorize( assembler.getEliminatedMatrix() );
            eli = solver.solve( assembler.getEliminatedRhs() );
            rhs -= assembler.getRhsModMatrix()*eli;
            gsInfo << "done." << std::endl;
        }

        if (!fullFn.empty())
        {
            gsInfo << "Write data to file " << fullFn << "..." << std::flush;
            gsFileData<> fd;
            fd << sys;
            fd << rhs;
            fd.save(fullFn);
            gsInfo << "done." << std::endl;
        }
    }
    else
    {
        gsInfo << "Found file " << fullFn << ". Read data just from there..." << std::flush;
        gsFileData<> fd(fullFn);
        fd.getFirst(sys);
        fd.getFirst(rhs);
        gsInfo << "done." << std::endl;
    }

    const index_t pdim = bases[dim][0]->size();
    const index_t vdim = sys.rows() - pdim;

    // -------------------------------------------------------------------------------------------------------------------------------------
    gsMatrix<> sol;


    if (directSolve)
    {
        // -------------------------------------------------------------------------------------------------------------------------------------
        gsInfo << "Solve with a direct solver..." << std::flush;
        gsSparseSolver<>::QR  solver;
        solver.analyzePattern( sys );
        solver.factorize     ( sys );
        sol = solver.solve( rhs );
        if (assembler.getZeroAverage())
        {
            sys.conservativeResize(sys.rows()-1,sys.cols()-1);
            rhs.resize(rhs.rows()-1,rhs.cols());
            sol.resize(sol.rows()-1,sol.cols());
        }
        gsMatrix<> tmp;
        gsInfo << "done." << std::endl;
        gsInfo << "Setup mass matrix (for pressure)..." << std::flush;
        gsSparseMatrix<> Mp;
        {
            gsBoundaryConditions<> pBc;
            gsGenericAssembler<> genassm( domain, *bases[dim][0], gsAssembler<>::defaultOptions(), &pBc );
            Mp = genassm.assembleMass();
        }
        gsInfo << "done. System matrix has " << Mp.rows() << " rows and " << Mp.cols() << " columns." << std::endl;
        gsInfo << "Correct mean of pressure..." << std::flush;
        gsPressureControlAdvanced<>::make(vdim,pdim,Mp)->apply(sol,tmp);
        sol.swap(tmp);
        gsInfo << "done. Updated by " << (sol-tmp).norm() << "." << std::endl;
    }
    else
    {

        // -------------------------------------------------------------------------------------------------------------------------------------

        const index_t number_of_velocity_preconders = (mg_phys == MultigridSpace::phys) ? 1 : dim;

        gsBlockOp<>::Ptr blockPrecon = gsBlockOp<>::make(number_of_velocity_preconders+1,number_of_velocity_preconders+1);

        std::vector< std::vector< gsSparseMatrix<real_t, RowMajor> > > transferMatrices(dim);
        std::vector< std::vector< gsMultiBasis<> > > basesLvsls(dim);

        gsInfo << "Coarsen basis...\n" << std::flush;
        for (index_t i=0; i<dim; ++i)
        {
            GISMO_ENSURE( bases[i].size() == 1, "Only the single patch case is covered." );
            gsInfo << "Handle basis(" << i << "," << 0 << ").\n";
            // The velocity has the boundary conditions from pde, the pressure does not have any.
            gsGridHierarchy<>::buildByCoarsening(gsMultiBasis<>(*bases[i][0]), pde.boundaryConditions(), gsAssembler<>::defaultOptions(), numLevels+1 )
                .moveMultiBasesTo( basesLvsls[i] )
                .moveTransferMatricesTo( transferMatrices[i] );
        }
        gsInfo << "done." << std::endl;

        if ( mg_phys == MultigridSpace::phys )
        {
            std::vector< gsSparseMatrix<real_t, RowMajor> > combinedTransferMatrices;
            gsInfo << "Combine transfer matrices...\n" << std::flush;
            const index_t sz = transferMatrices[0].size();
            combinedTransferMatrices.resize(sz);
            for (index_t j=0; j<sz; ++j)
            {
                index_t rows = 0;
                index_t cols = 0;
                for (index_t i=0; i<dim; ++i)
                {
                    rows += transferMatrices[i][j].rows();
                    cols += transferMatrices[i][j].cols();
                }
                combinedTransferMatrices[j].resize( rows, cols );
                index_t row_cnt = 0;
                index_t col_cnt = 0;
                for (index_t i=0; i<dim; ++i)
                {
                    addAsBlock( transferMatrices[i][j], row_cnt, col_cnt, combinedTransferMatrices[j] );
                    row_cnt += transferMatrices[i][j].rows();
                    col_cnt += transferMatrices[i][j].cols();
                }
            }
            gsInfo << "done." << std::endl;

            gsInfo << "Constructing combined multigrid preconditioner for all velocity diretions... " << std::flush;

            // set up the multigrid precondtioner for v
            gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make(sys.block(0,0,vdim,vdim), combinedTransferMatrices);

            mg->setNumPreSmooth( numPreSmooth );
            mg->setNumPostSmooth( numPostSmooth );
            mg->setNumCycles( cycles );

            mg->setNumOfSweeps( sweeps );

            GISMO_ASSERT( numLevels+1 == mg->numLevels(), "Inconsistent.");

            for (index_t i = 1; i < mg->numLevels(); ++i)
            {
                switch( smoother ) {
                    case Smoother::Richardson:                            mg->setSmoother(i, makeRichardsonOp(mg->matrix(i),damping_v*outerdamping_v)); break;
                    case Smoother::Jacobi:                                mg->setSmoother(i, makeJacobiOp(mg->matrix(i),damping_v*outerdamping_v)); break;
                    case Smoother::GaussSeidel:                           mg->setSmoother(i, makeGaussSeidelOp(mg->matrix(i))); break;
                    case Smoother::MassRichardsonSubspaceCorrection:      GISMO_ASSERT( false, "This combination is not possible." ); break;
                    case Smoother::MassRichardsonSubspaceCorrectionGeo:   GISMO_ASSERT( false, "This combination is not possible." ); break;
                }
            }

            blockPrecon->addOperator(0,0,mg);
            gsInfo << "done. System matrix has " << mg->rows() << " rows and " << mg->cols() << " columns." << std::endl;

        }
        else
        {
            index_t pre_dim = 0;
            for (index_t j=0; j<dim; ++j)
            {
                gsInfo << "Constructing multigrid preconditioner for direction " << j << "... " << std::flush;

                gsSparseMatrix<> Kv;
                assembleParameterStiffnessForTensorProductSpace(basesLvsls[j][numLevels][0], pde.boundaryConditions(), Kv);
                if ( mg_phys == MultigridSpace::partly )
                {
                    const index_t curr_dim = Kv.rows();
                    Kv = sys.block(pre_dim,pre_dim,curr_dim,curr_dim);
                    pre_dim += curr_dim;
                }

                gsDebug << Kv.rows() << "x" << Kv.cols() << "\n";

                // set up the multigrid precondtioner for v
                gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make(give(Kv), transferMatrices[j]);

                mg->setNumPreSmooth( numPreSmooth );
                mg->setNumPostSmooth( numPostSmooth );
                mg->setNumCycles( cycles );

                mg->setNumOfSweeps( sweeps );

                GISMO_ASSERT( numLevels+1 == mg->numLevels(), "Inconsistent.");

                for (index_t i = 1; i < mg->numLevels(); ++i)
                {
                    switch( smoother ) {
                        case Smoother::Richardson:                            mg->setSmoother(i, makeRichardsonOp(mg->matrix(i),damping_v*outerdamping_v)); break;
                        case Smoother::Jacobi:                                mg->setSmoother(i, makeJacobiOp(mg->matrix(i),damping_v*outerdamping_v)); break;
                        case Smoother::GaussSeidel:                           mg->setSmoother(i, makeGaussSeidelOp(mg->matrix(i))); break;
                        case Smoother::MassRichardsonSubspaceCorrection:      mg->setSmoother(i, gsPreconditionerFromOp<>::make(mg->underlyingOp(i),makeSubspaceCorrectedMassSmootherOperator(basesLvsls[j][i][0], damping_v, pde.boundaryConditions()), outerdamping_v) ); break;
                        case Smoother::MassRichardsonSubspaceCorrectionGeo:   mg->setSmoother(i, gsPreconditionerFromOp<>::make(mg->underlyingOp(i),makeSubspaceCorrectedMassSmootherOperator(domain.patch(0), basesLvsls[j][i][0], damping_v, pde.boundaryConditions()), outerdamping_v) ); break;
                    }
                }

                blockPrecon->addOperator(j,j,mg);
                gsInfo << "done. System matrix has " << mg->rows() << " rows and " << mg->cols() << " columns." << std::endl;

            }
        }

        gsSparseMatrix<> Mp;
        {
            gsInfo << "Setup mass matrix (for pressure)..." << std::flush;
            gsBoundaryConditions<> pBc;
            gsGenericAssembler<> genassm( domain, gsMultiBasis<>(*bases[dim][0]), gsAssembler<>::defaultOptions(), &pBc );
            Mp = genassm.assembleMass();
            gsInfo << "done. System matrix has " << Mp.rows() << " rows and " << Mp.cols() << " columns." << std::endl;
        }

        if (pprecond==PressurePreconditioner::pressure_mass)
        {
            gsInfo << "Setup inverse of mass matrix operator (for pressure)..." << std::flush;
            gsKroneckerOp<>::Ptr Mp_param_inv;
            assembleParameterMassInverseForTensorProductSpace((*bases[dim][0]),gsBoundaryConditions<>(),Mp_param_inv);
            blockPrecon->addOperator(number_of_velocity_preconders,number_of_velocity_preconders,gsScaledOp<>::make(Mp_param_inv,damping_p));
            gsInfo << "done. System matrix has " << Mp_param_inv->rows() << " rows and " << Mp_param_inv->cols() << " columns." << std::endl;
        }
        else if (pprecond==PressurePreconditioner::pressure_mass_geo)
        {
            gsInfo << "Setup inverse of mass matrix operator (for pressure)..." << std::flush;
            gsKroneckerOp<>::Ptr Mp_param_inv;
            assembleRankOneMassInverse(domain.patch(0),(*bases[dim][0]),gsBoundaryConditions<>(),Mp_param_inv);
            blockPrecon->addOperator(number_of_velocity_preconders,number_of_velocity_preconders,gsScaledOp<>::make(Mp_param_inv,damping_p));
            gsInfo << "done. System matrix has " << Mp_param_inv->rows() << " rows and " << Mp_param_inv->cols() << " columns." << std::endl;
        }
        else if (pprecond==PressurePreconditioner::pressure_jacobi)
        {
            gsInfo << "Setup gsJacobiOp operator (for pressure)..." << std::flush;
            gsJacobiOp< gsSparseMatrix<> >::Ptr Mp_jacobi = gsJacobiOp< gsSparseMatrix<> >::make(Mp);
            Mp_jacobi->setDamping(damping_p);
            gsPressureControlAdvanced<>::Ptr pressCtrl = gsPressureControlAdvanced<>::make(0,pdim,Mp);
            blockPrecon->addOperator(number_of_velocity_preconders,number_of_velocity_preconders,gsProductOp<>::make(Mp_jacobi,pressCtrl));
            gsInfo << "done." << std::endl;
        }
        else if (pprecond==PressurePreconditioner::pressure_sgs)
        {
            gsInfo << "Setup symmetric gsGaussSeidelOp operator (for pressure)..." << std::flush;
            gsGaussSeidelOp< gsSparseMatrix<>, gsGaussSeidel::symmetric >::Ptr Mp_SGS = gsGaussSeidelOp< gsSparseMatrix<>, gsGaussSeidel::symmetric >::make(Mp);
            gsPressureControlAdvanced<>::Ptr pressCtrl = gsPressureControlAdvanced<>::make(0,pdim,Mp);
            blockPrecon->addOperator(number_of_velocity_preconders,number_of_velocity_preconders,gsProductOp<>::make(Mp_SGS,pressCtrl));
            gsInfo << "done." << std::endl;
        }

        gsLinearOperator<>::Ptr precon = blockPrecon;

        if (gaussSeidelSweeps>0)
        {
             gsSparseMatrix<> exactPreconMat(vdim+pdim,vdim+pdim);
             addAsBlock( sys.block(0,0,vdim,vdim), 0, 0, exactPreconMat);
             addAsBlock( Mp, vdim, vdim, exactPreconMat);
             memory::shared_ptr< gsSparseMatrix<> > exactPreconMatPtr = exactPreconMat.moveToPtr();

             gsPreconditionerOp<>::Ptr fgs = gsGaussSeidelOp< gsSparseMatrix<> >::make(exactPreconMatPtr);
             gsPreconditionerOp<>::Ptr bgs = gsGaussSeidelOp< gsSparseMatrix<>, gsGaussSeidel::reverse >::make(exactPreconMatPtr);
             fgs->setNumOfSweeps(gaussSeidelSweeps);
             bgs->setNumOfSweeps(gaussSeidelSweeps);
             //gsPreconditionerOp<>::Ptr sgs = gsGaussSeidelOp< gsSparseMatrix<>, gsGaussSeidel::symmetric >::make(exactPreconMatPtr);sgs->setNumOfSweeps(gaussSeidelSweeps);
             precon = gsProductOfSteppableOperatorsOp<>::make( fgs, gsPreconditionedRichardsonOp<>::make(makeMatrixOp(exactPreconMatPtr),precon), bgs );
        }

        if (outerSweeps>1)
        {
            gsInfo << "Setup the required interfaces for outerSweep..." << std::flush;
            gsBlockOp<>::Ptr exactPrecon = gsBlockOp<>::make(2,2);
            exactPrecon->addOperator(0,0,makeMatrixOp( sys.block(0,0,vdim,vdim) ));
            exactPrecon->addOperator(1,1,makeMatrixOp( Mp ));

            gsPreconditionedRichardsonOp<>::Ptr precon2 = gsPreconditionedRichardsonOp<>::make(exactPrecon,precon);
            precon2->setNumOfSweeps(outerSweeps);
            precon2->setDamping(outerdamping);
            precon = precon2;
            gsInfo << "done.\n";
        }

        // -------------------------------------------------------------------------------------------------------------------------------------

        gsInfo << "Setup Minres solver..." << std::flush;
        gsSumOp<>::Ptr sys_op = gsSumOp<>::make(
            makeMatrixOp(sys),
            gsPressureRankOne<>::make(1,vdim,pdim,Mp)
        );
        gsMinimalResidual<> solver(sys_op,precon);
        gsInfo << "done. System matrix has " << sys_op->rows() << " rows and " << sys_op->cols() << " columns." << std::endl;

        sol.setRandom( sys.rows(), 1 );

        // First project into subspace
        gsMatrix<> tmp;
        gsPressureControlAdvanced<>::make(vdim,pdim,Mp)->apply(sol,tmp); sol.swap(tmp);

        sys_op->apply(sol,tmp); const real_t resNorm0 = (rhs - tmp).norm();
        gsInfo << "Residual norm:     " << resNorm0 << "\n";
        real_t resNorm(0), oldResNorm = resNorm0;

        real_t minReduction = 1e6;
        int numIter = 0;
        solver.setTolerance(0);
        solver.setMaxIterations(maxIter+3);

        solver.initIteration( rhs, sol );

        // -------------------------------------------------------------------------------------------------------------------------------------

        do
        {
            solver.step(sol);

            sys_op->apply(sol,tmp); resNorm = (rhs - tmp).norm();

            gsInfo << "Residual norm:     " << std::left << std::setw(15) << resNorm << "          reduction:  1 / " << std::setprecision(3) << (oldResNorm/resNorm) << std::setprecision(6) << "\n";
            minReduction = math::min(minReduction, oldResNorm/resNorm);
            oldResNorm = resNorm;

            ++numIter;
        } while (resNorm / resNorm0 > tol && numIter < maxIter && gsIsfinite(resNorm));

        //gsPressureControlAdvanced<>::make(vdim,pdim,Mp)->apply(sol,tmp); sol.swap(tmp);

        if (resNorm / resNorm0 > tol || ! gsIsfinite(resNorm))
            gsInfo << "Did not converge.\n";
        else
            gsInfo << "Converged in " << numIter << " iterations.\n";
        gsInfo << "Average convergence factor:  1 / " << std::setprecision(3) << math::pow(resNorm0 / resNorm, 1.0 / numIter) << std::setprecision(6) << "\n";
        gsInfo << "Worst   convergence factor:  1 / " << std::setprecision(3) << minReduction << std::setprecision(6) << "\n";
        gsInfo << "\n\n";

        if (writeLog)
        {
            std::fstream log("out.txt", std::fstream::out | std::fstream::app);

            log << "stokes\t";
            if (resNorm / resNorm0 > tol || ! gsIsfinite(resNorm))
                log << "dnc";
            else
                log << numIter;
            std::vector<gsOptionList::OptionListEntry> entr = opt.getAllEntries();
            for (size_t i=0; i<entr.size(); ++i)
            log << "\t" << entr[i].label << "=" << entr[i].val;
            log << "\n";
        }
    }

    // -------------------------------------------------------------------------------------------------------------------------------------

    if (computeError)
    {
        gsInfo << "Compute error..." << std::flush;
        std::vector<gsMatrix<real_t> > coefs;
        for (size_t s=0; s < phySpace.size(); ++s)
            coefs.push_back(assembler.reconstructSolution(s,sol,eli));
        std::vector<gsFunction<real_t>*> svec(2);
        svec[0] = const_cast<gsFunctionExpr<real_t>*>(&u);
        svec[1] = const_cast<gsFunctionExpr<real_t>*>(&p);
        gsRecipeDistance dist(pde.domain(),phySpace,coefs,svec);
        dist.assemble();
        gsInfo << "done.\n\n" << std::endl;

        gsInfo << "Error in v: " << math::sqrt(dist.getDistance(0).sum()) << "\n";


        if (true)
        {
            gsGeneralizedPoissonAssembler<real_t> assm(domain, *bases[dim][0], gsBoundaryConditions<>(), p, 1, gsAssembler<>::defaultOptions());
            gsInfo << "Error in p: " << computeL2DistanceUpToConst(assm.constructSolution(sol.bottomRows(pdim)), p, false, 3*pdim) << "\n";
        }
        else
            gsInfo << "Error in p: " << math::sqrt(dist.getDistance(1).sum()) << "\n";
    }




    // -------------------------------------------------------------------------------------------------------------------------------------

    gsInfo << "Free variables..." << std::flush;
    for (size_t s=0; s<bases.size();++s)
        freeAll(bases[s]);
    freeAll(phySpace);
    gsInfo << "done." << std::endl;

    // -------------------------------------------------------------------------------------------------------------------------------------
    return 0;
}



