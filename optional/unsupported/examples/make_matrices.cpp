//
// Simple program to assemble IGA matrices and write them to a file
// in Matrix Market format.
//
// Clemens Hofreither, 2014-2016
//

#include <iostream>
#include <iomanip>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsIO/gsMatrixIO.h>

using namespace std;
using namespace gismo;


string makeFilename(int d, int p, int n, const string& bc)
{
    stringstream ss;
    ss << "poisson_";
    ss << ((bc == "D") ? "dir" : "neu");
    ss << "_d" << d << "_p" << p << "_n" << n;
    return ss.str();
}

void argumentError(const char * errmsg)
{
    cerr << errmsg << endl;
    exit(-1);
}

int main(int argc, char *argv[])
{
    std::string geometry("1");
    index_t numIntervals = 50;
    index_t degree = 2;
    string boundaryCond = "D";
    bool refine = false;

    gsCmdLine cmd("Computes mass and stiffness matrices for IGA and saves them to disk.");
    cmd.addString("g", "geometry", "Geometry", geometry);
    cmd.addInt("n", "num-intervals",
            "Number of knot spans per coordinate direction", numIntervals);
    cmd.addInt("p", "degree",
            "Degree of the B-spline discretization space", degree);
    cmd.addString("b", "boundary",
            "Boundary conditions: D (Dirichlet) or N (Neumann)", boundaryCond);
    cmd.addSwitch("r", "refine",
            "If set, uniformly refines the basis and produces a prolongation matrix", refine);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (numIntervals < 1)
        argumentError("Number of dofs must be positive.");

    if (degree < 1)
        argumentError("Degree must be positive.");

    if (boundaryCond != "D" && boundaryCond != "N")
        argumentError("Boundary condition must be either D or N.");

    /******************** Define Geometry ********************/

    gsGeometry<>::uPtr geo;

    //gsGeometry<> * geo = BSplineQuarterAnnulus();
    //gsGeometry<> * geo = approximateQuarterAnnulus( degree );

    if (geometry=="1")
        geo = gsNurbsCreator<>::BSplineUnitInterval(static_cast<short_t>(degree));
    else if (geometry=="2")
        geo = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(degree));
    else if (geometry=="3")
        geo = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(degree));
    else if (geometry=="BSplineQuarterAnnulus")
        geo = gsNurbsCreator<>::BSplineQuarterAnnulus(2);
    else if ( gsFileManager::fileExists(geometry) )
    {
        geometry = gsFileManager::find(geometry);
        cout << "Reading file " << geometry << ".\n";
        gsMultiPatch<> mp;
        gsFileData<> fileData(geometry);
        if (!fileData.has< gsMultiPatch<> >()) argumentError("No multipatch object found.");
        fileData.getFirst< gsMultiPatch<> >(mp);
        if ( mp.nPatches() != 1 )              argumentError("There is not exactly one patch.");
        geo = mp[0].clone();
    }
    else
        argumentError( "Invalid geometry." );

    gsFunctionExpr<> f, g;

    switch (geo->geoDim())
    {
        case 1:
            f = gsFunctionExpr<>("(pi^2 ) * sin(pi*x)", 1);
            g = gsFunctionExpr<>("sin(pi*x)", 1);
            break;
        case 2:
            f = gsFunctionExpr<>("(2*pi^2 ) * sin(pi*x) * sin(pi*y)", 2);
            g = gsFunctionExpr<>("sin(pi*x) * sin(pi*y)", 2);
            break;
        case 3:
            f = gsFunctionExpr<>("(3*pi^2 ) * sin(pi*x) * sin(pi*y) * sin(pi*z)", 3);
            g = gsFunctionExpr<>("sin(pi*x) * sin(pi*y) * sin(pi*z)", 3);
            break;
        default:
            cerr << "Invalid geometry dimension." << endl;
            return -1;
    }

    // set up the boundary value problem
    gsConstantFunction<> zero(0.0, geo->geoDim());

    condition_type::type bc_type;
    gsFunction<> * bc_func;

    if (boundaryCond == "D")
    {
        bc_type = condition_type::dirichlet;
        bc_func = &g;
    }
    else // "N"
    {
        bc_type = condition_type::neumann;
        bc_func = &zero;
    }

#if 0
    //gsConstantFunction<> one (1.0, geo->geoDim());
    //gsBVProblem<> bvp( geo, new gsConvDiffRePde<>(&one, 0, &one, &f) ); // -Delta u + u = f
    gsBVProblem<> bvp( geo, new gsPoissonPde<>(f, geo->geoDim()) ); // -Delta u = f

    bvp.addCondition( boundary::west,  bc_type, bc_func );
    bvp.addCondition( boundary::east,  bc_type, bc_func );
    if (geo->geoDim() >= 2)
    {
        bvp.addCondition( boundary::south, bc_type, bc_func );
        bvp.addCondition( boundary::north, bc_type, bc_func );
    }
    if (geo->geoDim() >= 3)
    {
        bvp.addCondition( boundary::front, bc_type, bc_func );
        bvp.addCondition( boundary::back,  bc_type, bc_func );
    }

#else
    //// new-style boundary conditions
    gsBoundaryConditions<> bc;

    bc.addCondition( boundary::west,  bc_type, bc_func );
    bc.addCondition( boundary::east,  bc_type, bc_func );
    if (geo->geoDim() >= 2)
    {
        bc.addCondition( boundary::south, bc_type, bc_func );
        bc.addCondition( boundary::north, bc_type, bc_func );
    }
    if (geo->geoDim() >= 3)
    {
        bc.addCondition( boundary::front, bc_type, bc_func );
        bc.addCondition( boundary::back,  bc_type, bc_func );
    }
#endif

    gsKnotVector<real_t> KV(0.0, 1.0, numIntervals-1, degree+1);
    cout << "Coarse knots: " << KV << endl;
    gsBasis<real_t> * tbasis;
    switch (geo->geoDim())
    {
        case 1: tbasis = new gsBSplineBasis<real_t>( KV ); break;
        case 2: tbasis = new gsTensorBSplineBasis<2,real_t>( KV, KV ); break;
        case 3: tbasis = new gsTensorBSplineBasis<3,real_t>( KV, KV, KV ); break;
        default: cerr << "Only dimensions 1-3 currently supported."; return -1;
    }

    cout << tbasis->numElements() << endl;

    gsSparseMatrix<real_t, RowMajor> transfer;
    if (refine)
    {
        cout << "Refining... " << flush;
        tbasis->uniformRefine_withTransfer(transfer);
        cout << "done." << endl;
    }

    cout << "Discretization space: dim=" << tbasis->dim() << " deg=" << tbasis->degree(0) << " dofs=" << tbasis->size() << endl;

    cout << "Assembling stiffness matrix... " << flush;

    gsSparseMatrix<real_t> M, A;

#if 0
    gsGalerkinMethod<real_t> galerkin(bvp, *tbasis);

    // assemble linear system
    galerkin.assemble();
    A = galerkin.linearSystem().matrix()->selfadjointView<Lower>();

    // assemble mass matrix
    gsSparseMatrix<real_t> * tmp = galerkin.assembleMass();
    M = tmp->selfadjointView<Lower>();
    delete tmp;
#else
    // new-style assembler (slightly slower)
    gsMultiPatch<> mp(*geo);
    gsMultiBasis<> mb(*tbasis);
    gsPoissonAssembler<> assm(mp, mb, bc, f);
    assm.assemble();

    A = assm.matrix().selfadjointView<Lower>();

    cout << "done.\nAssembling mass matrix... " << flush;

    gsGenericAssembler<real_t> genassm(mp, mb,
    gsGenericAssembler<real_t>::defaultOptions(), &bc);
    genassm.assembleMass();
    M = genassm.fullMatrix();
#endif

    cout << "done.\nSaving... " << flush;

    // export matrices
    const string fname = makeFilename( geo->geoDim(), tbasis->degree(0), tbasis->component(0).numElements(), boundaryCond );
    saveMarket(M, fname + "_mass.mtx");
    saveMarket(A, fname + "_stiff.mtx");
    if (refine) {
        saveMarket(transfer, fname + "_prolong.mtx");
    }
    cout << "done." << endl;

    cout << "Mass matrix was saved to \"" << fname << "_mass.mtx\"\n";
    cout << "Stiffness matrix was saved to \"" << fname << "_stiff.mtx\"\n\n";

    delete tbasis;
    return 0;
};
