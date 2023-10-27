/** @file gsParallelMultiGridExample.cpp

    @brief Provides test examples for parallel multigrid algorithms

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer, S. Takacs
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI

#include <gismo.h>
#include <gismo_dev.h>

#include <gsSolver/gsTimedOp.h>
#include <gsCore/gsLinearCombinationOfFunctionsFunction.h>

//#include <gsIO/gsMatrixIO.h>
#include <gsTensor/gsTensorTools.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsGeneralizedPoissonAssembler.h>

#include <gsIO/gsCmdLineWithEnumSupport.h>

#include <gsMultiGrid/gsParallelGridHierarchy.h>
#include <gsIETI/gsParallelOperator.h>
#include <gsIETI/gsParallelCG.h>
#include <gsMpi/gsMpi.h>
#include <gsMultiGrid/gsParallelMultiPatchPreconditioners.h>
#include <gsSolver/gsParallelPreconditioner.h>
#include <gsCore/gsConnectedBoundaryComponent.h>

#include <gsMultiGrid/gsMassSmoother.h>


using namespace std;
using namespace gismo;

namespace Smoother {
enum type {
    Richardson,
    Jacobi,
    GaussSeidel,
    MassRichardsonSubspaceCorrectionAdditiveMPDDD,
    MassRichardsonSubspaceCorrectionAdditiveMPDDD2,
};
}

/// A locally used geometry
gsTensorBSpline<2>::uPtr BSplineMySquare(int deg)
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

gsBoundaryConditions<real_t> updateBoundaryConditionsAfterSplit(const gsBoundaryConditions<real_t>& bcOld, const gsMultiPatch<real_t>& splittedPatch)
{
    gsBoundaryConditions<> newBC;
    int d = splittedPatch.parDim();
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

int main(int argc, char *argv[])
{
    string geometry("2");
    bool showBasis = false;
    string boundaryCondition("n");
    int numRefine = 3;
    int mult = 1;
    int degree = 2;
    int numLevels = -1;
    int numSplit = 0;
    Smoother::type smoother = Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD2;
    std::string smoother_name;
    int numPreSmooth = 1;
    int numPostSmooth = 1;
    int cycles = 1;
    real_t alpha = 1.;
    bool useCG = false;
    bool compEigs = false;
    bool convergenceRate = true;
    bool monitorl2 = false, monitorL2 = false;
    bool writeLog = false;
    bool detailedLog = false;
    bool noSAOp = false;

    real_t tol = 1e-8;
    real_t damping = -1.0;
    real_t outerDamping = 1.;
    bool plot = false;
    int maxIter = 1000;
    bool tensorAssemble = false;

    gsCmdLineWithEnumSupport cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addString("g", "geometry",           "Specification of the geometry (overrides dimension)",             geometry         );
    cmd.addSwitch(     "show-basis",         "Shows the chosen basis and quits.",                               showBasis        );
    cmd.addString("b", "boundary-condition", "Boundary condition",                                              boundaryCondition);
    cmd.addInt   ("r", "uniformRefine",      "Number of uniform h-refinement steps to perform before solving",  numRefine        );
    cmd.addInt   ("m", "multiplicity",       "Multiplicity of knots to insert when refining",                   mult             );
    cmd.addInt   ("p", "degree",             "Degree of the B-spline discretization space",                     degree           );
    cmd.addInt   ("l", "levels",             "Number of levels to use for multigrid iteration",                 numLevels        );
    cmd.addInt   ("",  "split",              "Split every patch uniformly into 2^d patches (default: 0)",       numSplit         );
    cmd.addInt   ("",  "presmooth",          "Number of pre-smoothing steps",                                   numPreSmooth     );
    cmd.addInt   ("",  "postsmooth",         "Number of post-smoothing steps",                                  numPostSmooth    );
    cmd.addEnum  ("s", "smoother",           "Smoothing method",                                                smoother         )
            .add(Smoother::Richardson,                                          "r",      "Richardson smoother"                                                       )
            .add(Smoother::Jacobi,                                              "j",      "Jacobi smoother"                                                           )
            .add(Smoother::GaussSeidel,                                         "gs",     "GaussSeidel smoother"                                                      )
            .add(Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD,       "mrs-ad", "mass smoother with subspace correction based on additive Dirichlet dd"     )
            .add(Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD2,      "mrs-ad2","mass smoother with subspace correction based on additive Dirichlet dd; experimental"     )
            .writeDescOfChosenOptTo(smoother_name);
    cmd.addReal  ("",  "damping",            "Damping factor for the smoother (handed over to smoother)",       damping          );
    cmd.addReal  ("",  "outerdamping",       "Damping factor for the smoother (globally)",                      outerDamping     );
    cmd.addInt   ("c", "cycles",             "Number of multi-grid cycles",                                     cycles           );
    cmd.addInt   ("",  "maxiter",            "Maximum number of iterations",                                    maxIter          );
    cmd.addReal  ("a", "alpha",              "alpha in \"- LAPLACE u + alpha u = f\"",                          alpha            );
    cmd.addSwitch(     "log",                "Write results to log file",                                       writeLog         );
    cmd.addSwitch(     "detailedLog",        "Writes detailed results to log file",                             detailedLog      );
    cmd.addSwitch(     "cg",                 "Use CG iteration",                                                useCG            );
    cmd.addSwitch(     "eigs",               "Compute eigenvalues of the preconditioned system. Requires --cg.",compEigs         );
    cmd.addSwitch(     "tensor",             "Assemble using tensor product (experimental)",                    tensorAssemble   );
    cmd.addReal  ("",  "tol",                "Tolerance for multigrid solver stopping criterion",               tol              );
    cmd.addSwitch(     "noConvergenceRate",  "Do not print convergence rate",                                   convergenceRate  );
    cmd.addSwitch(     "monitor-L2",         "Monitor the L2 errors over the iteration",                        monitorL2        );
    cmd.addSwitch(     "monitor-l2",         "Monitor the discrete l2 errors over the iteration",               monitorl2        );
    cmd.addSwitch(     "plot",               "Plot result in ParaView format",                                  plot             );
    cmd.addSwitch("","NoSubassembledOperators", "Do not use the subassembledOperator", noSAOp);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //enable log if detailed log is enabled.
    writeLog = detailedLog ? true : writeLog;

    /******************** Init MPI **************************/
    // Initialize the MPI environment
    const gsMpi & mpi = gsMpi::init(argc, argv);

    // Get the world communicator
    gsMpiComm comm = mpi.worldComm();

    //Get size and rank of the processor
    int _size = comm.size();
    int _rank = comm.rank();

    if (0==_rank)
        gsInfo<<"Running on "<<_size<<" processes.\n";
    comm.barrier();


    /******************** Define Geometry ********************/

    gsStopwatch time;
    gsMultiPatch<> mp;

    std::string orig_geometry = geometry;

    if (geometry=="1")
    {
        gsGeometry<>::uPtr geo = gsNurbsCreator<>::BSplineUnitInterval(1);

        // offset the computational domain by 0.5 in each direction (for historical reasons)
        gsVector<> offs = 0.5 * gsVector<>::Ones( geo->parDim() ); geo->translate( offs );

        mp = gsMultiPatch<>(*geo);
    }
    else if (geometry=="2")
    {
        gsGeometry<>::uPtr geo = gsNurbsCreator<>::BSplineSquareDeg(1);

        // offset the computational domain by 0.5 in each direction (for historical reasons)
        gsVector<> offs = 0.5 * gsVector<>::Ones( geo->parDim() ); geo->translate( offs );

        mp = gsMultiPatch<>(*geo);
    }
    else if (geometry=="3")
    {
        gsGeometry<>::uPtr geo = gsNurbsCreator<>::BSplineCube(1);

        // offset the computational domain by 0.5 in each direction (for historical reasons)
        gsVector<> offs = 0.5 * gsVector<>::Ones( geo->parDim() ); geo->translate( offs );

        mp = gsMultiPatch<>(*geo);
    }
    else if (geometry=="4")
    {
        gsGeometry<>::uPtr geo = approximateQuarterAnnulus(2);

        // offset the computational domain by 0.5 in each direction (for historical reasons)
        gsVector<> offs = 0.5 * gsVector<>::Ones( geo->parDim() ); geo->translate( offs );

        mp = gsMultiPatch<>(*geo);
    }
    else if (geometry=="5")
    {
        gsGeometry<>::uPtr geo = BSplineMySquare(1);

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
            gsFileData<> fileData(geometry);
            if (!fileData.has< gsMultiPatch<> >()) { cerr << "No multipatch object found in file " << geometry << ".\n"; return 1; }
            fileData.getFirst< gsMultiPatch<> >(mp);
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

    if(_rank == 0) gsInfo << "The geometry consists of " << mp.nPatches() << " patches.\n";


    if ( showBasis )
    {
        if(_rank == 0) gsInfo << mp << endl;
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
        if( smoother == Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD
                || smoother == Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD2
                )
            numLevels = maxMRSLevels(degree,numRefine);
        if(_rank == 0) gsInfo << "The number of levels was chosen to be " << numLevels << ".\n";
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
        case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD:       damping = 0.09; break;
        case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:      damping = 0.09; break;
        }
    }
    if(!convergenceRate && !useCG)
    {
        cerr << "Preventing calculation of convergence rate can only be performed with CG.\n"; return -1;
    }

    // ---------------------------------------------------------------------------------------------------------------------------

    if(_rank == 0) gsInfo << "Source function: " << *f << ".\n" << "\n";
    if(_rank == 0) gsInfo << "Exact solution:  " << *g << ".\n" << "\n";


    // set up boundary conditions

    gsConstantFunction<> zero(0.0, mp.geoDim());
    //gsConstantFunction<> one (1.0, mp.geoDim());

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
        if(_rank == 0) gsInfo << "Added " << d_nr << " Dirichlet and " << n_nr << " Neumann boundary condtions.\n";
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
        if(_rank == 0) gsInfo << "Added " << d_nr << " Dirichlet and " << n_nr << " Neumann boundary condtions.\n";
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
        if(_rank == 0) gsInfo << "Added " << i << " boundary condtions.\n";
    }

    for (index_t i=0; i<numSplit; ++i)
    {
        if(_rank == 0) gsInfo << "Split multipatch object uniformly... " << flush;
        mp = mp.uniformSplit();

        bc = updateBoundaryConditionsAfterSplit(bc,mp);
        if(_rank == 0) gsInfo << "done." << endl;
    }

    if(_rank == 0) bc.print(gsInfo);


    if( numSplit )
        if(_rank == 0) gsInfo << "The geometry consists of " << mp.nPatches() << " patches.\n";


    gsMultiBasis<> mb(mp);
    for (size_t k = 0; k < mb.nBases(); ++k)
    {
        for(index_t i=0;i<2;++i)
         mb.basis(k).uniformCoarsen();
    }

    //Distribute the patches to the processors
    gsSortedVector<size_t> myPatches;
#if 1
    {
        // Assign based on the # of dofs and try to have similar numbered patches
        // together (assuming them to be close)
        const index_t nPatches   = mp.nPatches();
        const index_t totalDofNr = mb.totalSize();
        const index_t nProc      = comm.size();
        index_t dofs             = 0;
        index_t asignee          = 0;
        index_t ownDofNr         = 0; // for output only
        for (index_t k=0; k<nPatches; ++k)
        {
            if (dofs*nProc > (asignee+1)*totalDofNr)
                asignee++;

            dofs += mb.size(k);
            // assign to processor "asignee"
            if (asignee == _rank)
            {
                myPatches.push_sorted_unique(k);
                ownDofNr += mb.size(k);
            }
        }
        GISMO_ENSURE( asignee < nProc, "Internal error." );

        std::ostringstream out;
        out<<"Patches of proc " << _rank << ":";
        for(size_t i=0; i<myPatches.size();++i )
            out<<" "<<myPatches[i];
        out<<" ("<<(ownDofNr*1000/totalDofNr)/10.<<"% of dofs)\n";

        gsInfo << out.str();
    }
#else
    unsigned ratio = (unsigned)(mp.nPatches()/comm.size());
    for(size_t np =0; np< ratio*comm.size();++np)
        if(np /ratio == (unsigned)_rank)
            myPatches.push_sorted_unique(np);
    for(index_t np=ratio*comm.size(); np<mp.nPatches();++np)
        if(np % comm.size() == (unsigned)_rank)
            myPatches.push_sorted_unique(np);
    gsInfo<<"My Patches: ";
    for(size_t i=0; i<myPatches.size();++i )
        gsInfo<<" "<<myPatches[i];
    gsInfo<<"\n";
#endif

    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    vector< gsMultiBasis<> > bases;
    vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices;
    vector< vector< gsSparseMatrix<real_t, RowMajor> > > fullTransferMatrices;

    gsOptionList assemblerOptions = gsAssembler<>::defaultOptions();
    if(noSAOp) assemblerOptions.addSwitch("NoSubassembledOperators","do not use SAOps",true);
    // Define coarse discretization space by refining the basis of the geometry
    for (int i = 0; i < numRefine - numLevels + 1; ++i)
        mb.uniformRefine();       // refine until coarsest level

    if (mult>1)
    {
        if(_rank == 0) gsInfo << "Reduce continuity by " << mult-1 << std::endl;
        mb.reduceContinuity(mult-1);
    }

    comm.barrier();
    const double timeSetupGeo = time.stop(); time.restart();

    if(_rank == 0) gsInfo << "Setup grid hierarchy..." << std::flush;
    // set up the hierarchy of spaces and transfer matrices between them
    const index_t refineKnots = 1;
    gsParallelGridHierarchy<real_t> GH = gsParallelGridHierarchy<real_t>::buildByRefinement(give(mb), bc, assemblerOptions,myPatches,comm, numLevels,refineKnots, mult);

    comm.barrier();
    const double timeSetupGH = time.stop(); time.restart();

    if(_rank == 0) gsInfo << "done.\nSetup transfer operators" << std::flush;
    std::pair<std::pair<std::vector<gsParallelOperator<real_t>::Ptr>,std::vector<gsParallelOperator<real_t>::Ptr> >,
            std::vector<gsParallelGlobalLocalHandler::Ptr> >
            transferOps = GH.getRestrictionAndProlongationOperators();

    const double timeSetupRestrProlong = time.stop(); time.restart();

    bases = GH.getMultiBases();

    if(_rank == 0) gsInfo << (useCG ? "CG preconditioned by multigrid" : "Multigrid");
    if(_rank == 0) gsInfo << " with " << numLevels << " levels using " << smoother_name;
    if(_rank == 0) gsInfo << "(" << numPreSmooth << "," << numPostSmooth << ") smoother" << " and " << cycles << "-cycle.\n";
    for (size_t i=0; i<mb.nBases(); ++i)
        if(_rank == 0) gsInfo << "Coarse discretization space: dim=" << mb[i].dim() << " deg=" << mb[i].degree(0) << " dofs=" << mb[i].size() << "\n";

    if(_rank == 0) gsInfo << "Setup gsGeneralizedPoissonAssembler and assemble... " << flush;
    //gsGenericAssembler<real_t> genassm(mp, bases.back(), assemblerOptions, &bc);

    gsMatrix<real_t> rhs,rhs_;
    gsGeneralizedPoissonAssembler<real_t> assm(mp, bases.back(), bc, *f, alpha, assemblerOptions); //for construction solution

    const double timeSetupAssembler = time.stop(); time.restart();

    std::vector<gsMatrixOp<gsSparseMatrix<real_t> >::Ptr > _localOps(mp.nPatches()); // for subassembled Ops
    gsMatrixOp<gsSparseMatrix<real_t> >::Ptr _localOp; //no subassembled Op
    if(!noSAOp)
    {
        assm.system().reserve(0,1);
        for(size_t k=0; k<myPatches.size();++k)
        {
            gsBoundaryConditions<real_t> bcLoc;
            bc.getConditionsForPatch(myPatches[k],bcLoc);
            gsGeneralizedPoissonAssembler<real_t> assmL(mp.patch(myPatches[k]), bases.back().basis(myPatches[k]), bcLoc, *f, alpha, assemblerOptions);
            assmL.assemble();
            // gsInfo<<"Local mat: \n"<<i<<": "<<assemblerLoc.system().matrix().toDense()<<"\n";
            _localOps[myPatches[k]] =  makeMatrixOp(assmL.system().matrix().moveToPtr());
            gsMatrix<real_t> rhsL = assmL.system().rhs();

            gsMatrix<unsigned> actives(rhsL.rows(),1);
            size_t count =0;
            for(unsigned i=0; i<(unsigned)bases.back().basis(myPatches[k]).size();++i)
                if(assmL.system().colMapper(0).is_free(i))actives(count++,0)=assm.system().colMapper(0).index(i,myPatches[k]);

            assm.system().pushToRhs(rhsL,actives); //build rhs
        }

        rhs_ = assm.system().rhs();
    }
    else
    {
        gsBoundaryConditions<real_t> reducedBc;
        gsMultiPatch<real_t> reducedPatch;
        gsMultiBasis<real_t> reducedBasis;
        for(size_t i=0; i<myPatches.size();++i)
        {
            gsBoundaryConditions<real_t> bcLoc;
            bc.getConditionsForPatch(myPatches[i],bcLoc);
            gsBoundaryConditions<real_t>::bcContainer container = bcLoc.allConditions();
            for(gsBoundaryConditions<real_t>::const_iterator it = container.begin(); it!=container.end();++it)
                reducedBc.add(i,it->side(),it->ctype(),it->function(),it->unknown(),it->unkComponent(),it->parametric());

            reducedPatch.addPatch(mp.patch(myPatches[i]).clone());
            reducedBasis.addBasis(bases.back().basis(myPatches[i]).clone().release());
        }
        reducedPatch.computeTopology();
        reducedBasis.setTopology(reducedPatch.topology());
        gsGeneralizedPoissonAssembler<real_t> reducedAssembler(reducedPatch,reducedBasis,reducedBc,*f, alpha,assemblerOptions);
        GH.getSubassembledTopology().back()->reorder(assm.system().colMapper(0),reducedAssembler.system().colMapper(0));

        reducedAssembler.assemble();
        _localOp =makeMatrixOp(reducedAssembler.system().matrix().moveToPtr());
        rhs = reducedAssembler.rhs();
    }

    comm.barrier();
    if(_rank == 0) gsInfo << "done." << endl;
    const double timeAssLocalMats = time.stop(); time.restart();

    if(_rank == 0) gsInfo << "Setup system matrix on all levels... " << flush;
    std::vector<gsMatrixOp<gsSparseMatrix<real_t> >::Ptr > c_ops; //for subassembled Op
    gsMatrixOp<gsSparseMatrix<real_t> >::Ptr c_op; //for no subassembled Ops
    std::vector<gsDofMapper> c_locMappers;
    gsDofMapper c_globMapper;
    std::pair<std::vector<typename gsParallelOperator<real_t>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > operators;
    if(!noSAOp)
        operators= GH.generateGalerkinProjection(_localOps,c_ops, c_locMappers,c_globMapper);
    else
        operators= GH.generateGalerkinProjection(_localOp,c_op, c_locMappers,c_globMapper);


    comm.barrier();
    if(_rank == 0) gsInfo << "done." << endl;
    const double timeAssGalerkinProjections= time.stop(); time.restart();



    if(_rank == 0) gsInfo << "Extract rhs... " << flush;
    const gsParallelGlobalLocalHandler & handler =*operators.second.back();

    if(!noSAOp)handler.extractLocalVector(rhs_,rhs);

    comm.barrier();
    if(_rank == 0) gsInfo << "done." << endl;
    const double timeExtractingRHS = time.stop(); time.restart();

    if(_rank == 0) gsInfo << "Setup coarse solver... " << flush;
    gsCoarseSolverAdapter<real_t>::Ptr cs;
    if(!noSAOp)
        cs = give(constructCoarseSolver<real_t>(c_ops,operators.second.front(),myPatches,c_locMappers,c_globMapper));
    else
        cs = give(constructCoarseSolver<real_t>(c_op,operators.second.front()));

    gsTimedOp<>::Ptr coarseSolver = gsTimedOp<>::make("CoarseSolver",cs,false);


    const double timeSetupCoarseSolver = time.stop(); time.restart();
    if(_rank == 0) gsInfo << "done." << endl;

    if(_rank == 0) gsInfo << "Perform casts... " << flush;
    std::vector<gsLinearOperator<real_t>::Ptr> linOps (operators.first.size());
    std::vector<gsLinearOperator<real_t>::Ptr> linOpsR (transferOps.first.first.size());
    std::vector<gsLinearOperator<real_t>::Ptr> linOpsP (transferOps.first.second.size());
    for(size_t l =0; l<operators.first.size();++l)
        linOps[l] = operators.first[l];
    for(size_t l =0; l<transferOps.first.first.size();++l)
    {
        linOpsR[l] = transferOps.first.first[l];
        linOpsP[l] = transferOps.first.second[l];
    }
    if(_rank == 0) gsInfo << "done." << endl;

    if(_rank == 0) gsInfo << "Setup multigrid solver... " << flush;
    gsMultiGridOp<> mg(linOps,linOpsP,linOpsR,coarseSolver );
    //gsMultiGridOp<> mg(give(Kfine), transferMatrices, gsNullOp<real_t>::make(transferMatrices.back().rows()) );

    vector< gsSparseMatrix<real_t> > massMatrices;

    mg.setNumPreSmooth( numPreSmooth );
    mg.setNumPostSmooth( numPostSmooth );
    mg.setNumCycles( cycles );

    comm.barrier();
    if(_rank == 0) gsInfo << "done." << endl;
    const double timeSetupMGClass = time.stop(); time.restart();
    double timeLocalSmoothers = 0;

    GISMO_ASSERT( outerDamping == 1 || smoother != Smoother::GaussSeidel, "Gauss-Seidel does not support --outerdamping" );

    if(_rank == 0) gsInfo << "Constructing smoothers... " << flush; comm.barrier();
    for (int i = mg.numLevels() == 1 ? 0 : 1; i < mg.numLevels(); ++i)
    {
        switch( smoother ) {
            case Smoother::Richardson:                            mg.setSmoother(i, makeRichardsonOp(mg.matrix(i),damping*outerDamping)); break;
            case Smoother::Jacobi:                                mg.setSmoother(i, makeJacobiOp(mg.matrix(i),damping*outerDamping)); break;
            case Smoother::GaussSeidel:                           mg.setSmoother(i, makeGaussSeidelOp(mg.matrix(i))); break;
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD:
            {

                gsStopwatch timer2;
                std::vector<gsLinearOperator<>::Ptr> localSmoothers = makeSubspaceCorrectedMassSmootherOperatorsDirichlet(bases[i],damping);
                timeLocalSmoothers += timer2.stop();

                mg.setSmoother(
                            i,
                            gsParallelAdditivePreconditionerOp<real_t>::make(
                                operators.first[i],
                                setupPiecewisePreconditioner<real_t>(
                                    *operators.first[i],
                                    *operators.second[i],
                                    give(localSmoothers),
                                    mp,
                                    bases[i],
                                    bc,
                                    assemblerOptions,
                                    myPatches,
                                    comm
                                    ),
                                outerDamping
                                )
                            );
                break;
            }
            case Smoother::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:
            {
                const index_t dim = bases[i].dim();
                gsMultiBasis<> proc_loc_mb;
                for (size_t k=0;k<myPatches.size();++k)
                    proc_loc_mb.addBasis(bases[i][myPatches[k]].clone().release());

                const index_t proc_loc_dofs = (*operators.second[i]).localSize();

                gsDofMapper dm; // Global mapper
                bases[i].getMapper(
                   (dirichlet::strategy)assemblerOptions.askInt("DirichletStrategy",11),
                   (iFace    ::strategy)assemblerOptions.askInt("InterfaceStrategy", 1),
                   bc,
                   dm,
                   0
                );

                // Global => Processor local
                gsMatrix<unsigned> glob2loc;
                {
                    gsMatrix<unsigned> inp(proc_loc_dofs,1);
                    inp.col(0) = gsVector<unsigned>::LinSpaced(proc_loc_dofs,0,proc_loc_dofs-1);
                    glob2loc.setZero((*operators.second[i]).globalSize(),1);
                    (*operators.second[i]).addLocalVectorToGlobal( inp, glob2loc );
                }

                std::vector< gsVector<index_t> > proc_loc_maps;
                for (size_t k=0;k<myPatches.size();++k)
                {
                    const index_t sz = proc_loc_mb[k].size();
                    index_t kk = myPatches[k];
                    gsVector<index_t> local(sz);
                    for (index_t j=0; j<sz; ++j)
                    {
                        local[j] = dm.is_free(j,kk) ? glob2loc(dm.index(j, kk),0) : -1u;
                    }
                    proc_loc_maps.push_back(give(local));
                }
                std::vector< std::vector< std::pair< typename gsBasis<>::Ptr, gsSparseMatrix<> > > > pieces =
                    constructPieces( proc_loc_mb, proc_loc_maps, proc_loc_dofs );

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
                            const real_t scalingK = std::pow(2,dim-dd) * (dim-dd) * std::pow( h/(2*degree+1), dim-dd );
                            const real_t scalingM = std::pow(2,dim-dd) * std::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
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
                gsInfo << "[" << _rank << ":" << localSmoothers.size() << "] " << std::flush;

                mg.setSmoother(
                    i,
                        gsParallelAdditivePreconditionerOp<real_t>::make(
                        operators.first[i],
                        give(smootherTransfers),
                        give(localSmoothers),
                        outerDamping
                    )
                );
                break;
            }
        }
    }
    if(_rank == 0) gsInfo << "done." << "\n"<< std::flush; comm.barrier();


    gsMatrix<> x;
    //x.setRandom( mg.nDofs(), 1 );
    x.setZero( mg.nDofs(), 1 );

    comm.barrier();
    const double timeSetupSmoother = time.stop(); time.restart();

    gsMatrix<real_t> tmp, tmp_acc;
    real_t resNorm0 = 0;
    real_t resNorm(0), oldResNorm = resNorm0;
    if(convergenceRate)
    {
        mg.underlyingOp()->apply(x,tmp);
        tmp = rhs - tmp;
        dynamic_cast<gsParallelOperator<real_t>*>(mg.underlyingOp().get())->accumulate(tmp,tmp_acc);
        resNorm0 = (tmp.transpose()*tmp_acc)(0,0);
        resNorm0 = math::sqrt(comm.sum(resNorm0));
        if(_rank == 0) gsInfo << "Residual norm:     " << resNorm0 << "\n";
        oldResNorm = resNorm0;
    }

    gsField<> sol;

    real_t l2Err, oldL2Err(0), eucl_error, old_eucl_error(0);
    if (monitorL2)
    {
        sol = assm.constructSolution(x);
        oldL2Err = computeL2Distance(sol, *g, false, 3*mg.nDofs());
        if(_rank == 0) gsInfo << "L2 error: " << oldL2Err << "\n";
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

    gsParallelPreconditionerOp<real_t>::Ptr mg_prec = gsParallelPreconditionerOp<real_t>::make(memory::make_shared_not_owned(&mg));
    gsParallelCG<real_t> cg( operators.first.back(), mg_prec );
    cg.setTolerance( tol );


    if (compEigs)
        cg.setCalcEigenvalues(true);

    if (useCG && convergenceRate)
        cg.initIteration( rhs, x );

    if (convergenceRate)       // solve using MG iteration
    {
        do
        {
            if (useCG)
                cg.step(x);
            else if (numLevels == 1)
                mg.smoothingStep(rhs, x);
            else
                mg.step(rhs, x);

            mg.underlyingOp()->apply(x,tmp);
            tmp = rhs - tmp;
            dynamic_cast<gsParallelOperator<real_t>*>(mg.underlyingOp().get())->accumulate(tmp,tmp_acc);
            resNorm = (tmp.transpose()*tmp_acc)(0,0);
            resNorm = math::sqrt(comm.sum(resNorm));

            if(_rank == 0) gsInfo << "Residual norm:     " << left << setw(15) << resNorm << "          reduction:  1 / " << setprecision(3) << (oldResNorm/resNorm) << setprecision(6) << "\n";
            minReduction = math::min(minReduction, oldResNorm/resNorm);
            oldResNorm = resNorm;

            if (monitorL2)
            {
                sol = assm.constructSolution(x);
                l2Err = computeL2Distance(sol, *g, false, 3*mg.nDofs());
                if(_rank == 0) gsInfo << "                                                                   |  L2 error:     "
                                      << left << setw(15) << l2Err << "          reduction:  1 / " << setprecision(3) << (oldL2Err/l2Err) << setprecision(6) << "\n";
                oldL2Err = l2Err;
            }

            if (monitorl2)
            {
                eucl_error = (exactDiscreteSol - x).norm();
                if(_rank == 0) gsInfo << "                                                                   |  l2 error:     "
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

    comm.barrier();
    const double timeSolve = time.stop();
    const double timeSetup = timeSetupGH+timeSetupRestrProlong+timeAssGalerkinProjections+timeSetupCoarseSolver+timeSetupMGClass+timeSetupSmoother;
    const double timeAssemble = timeAssLocalMats+timeExtractingRHS;

    if(_rank == 0)
    {
        if (convergenceRate && ( resNorm / resNorm0 > tol || ! gsIsfinite(resNorm)  ) )
            gsInfo << "Did not converge.\n";
        else
            gsInfo << "Converged in " << numIter << " iterations.\n";

        if(convergenceRate)
        {
            gsInfo << "Average convergence factor:  1 / " << setprecision(3) << math::pow(resNorm0 / resNorm, 1.0 / numIter) << setprecision(6) << "\n";
            gsInfo << "Worst   convergence factor:  1 / " << setprecision(3) << minReduction << setprecision(6) << "\n";
            gsInfo << "\n";
        }
        gsInfo << "Setup Geometry:  "; formatTime(gsInfo, timeSetupGeo);                 gsInfo << "\n";
        gsInfo << "Assembling time: "; formatTime(gsInfo, timeAssemble);                 gsInfo << "\n";
        gsInfo << " local matrices: "; formatTime(gsInfo, timeAssLocalMats);             gsInfo << "\n";
        gsInfo << " extract rhs:    "; formatTime(gsInfo, timeExtractingRHS);            gsInfo << "\n";
        gsInfo << "MG time setup:   "; formatTime(gsInfo, timeSetup);                    gsInfo << "\n";
        gsInfo << " grid hierarchy: "; formatTime(gsInfo, timeSetupGH);                  gsInfo << "\n";
        gsInfo << " restr/prol:     "; formatTime(gsInfo, timeSetupRestrProlong);        gsInfo << "\n";
        gsInfo << " Galerkin proj:  "; formatTime(gsInfo, timeAssGalerkinProjections);   gsInfo << "\n";
        gsInfo << " coarse solver:  "; formatTime(gsInfo, timeSetupCoarseSolver);        gsInfo << "\n";
        gsInfo << " MG object:      "; formatTime(gsInfo, timeSetupMGClass);             gsInfo << "\n";
        gsInfo << " smoother:       "; formatTime(gsInfo, timeSetupSmoother);            gsInfo << "\n";
        gsInfo << " assembler:      "; formatTime(gsInfo, timeSetupAssembler);           gsInfo << "\n";

        if (timeLocalSmoothers>0)
        {
            gsInfo << "   local smooth: "; formatTime(gsInfo, timeLocalSmoothers);           gsInfo << "\n";
            gsInfo << "MG time solving: "; formatTime(gsInfo, timeSolve);
        }
        gsInfo << "         (avg. "; formatTime(gsInfo, timeSolve/numIter);          gsInfo << " per iteration)" << "\n";
        gsInfo << " coarse solver:  "; formatTime(gsInfo, coarseSolver->getTime());
        gsInfo << "\n";
        gsInfo << "Total time:      "; formatTime(gsInfo, timeSetupGeo+timeAssemble+timeSetup+timeSolve);  gsInfo << "\n";
        gsInfo << "\n";
    }
    if (monitorL2 || plot)
        sol = assm.constructSolution(x);

    if (monitorL2)
    {
        // Compute L2 error
        const real_t error = computeL2Distance(sol, *g, false, 3*mg.nDofs());
        if(_rank == 0) gsInfo << "L2 error: " << error << "\n";
    }

    if (monitorl2)
    {
        const gsVector<> f_err = exactDiscreteSol - x;
        const real_t eucl_error_f = f_err.norm();
        if(_rank == 0) gsInfo << "l2 error: " << eucl_error_f << "\n";
        if(_rank == 0) gsInfo << "Discrete energy error: " << f_err.dot(mg.matrix()*f_err)  << "\n";
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
        if(_rank == 0) gsInfo << "Eigenvalues: Minimum, essentialMinimum, maximum, essentialConditionNumber: " << min << ", " << essMin << ", " << max << ", " << (max/essMin) << std::endl;
    }

    if (writeLog && _rank == 0)
    {
        fstream log("out.txt", fstream::out | fstream::app);

        log << "parallel_multigrid_example" << (useCG ? "_cg" : "" ) << "\t"
            << geometry << "\t"
            << comm.size() << "\t"
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
            << std::fixed << std::setprecision(6)
            << timeAssemble << "\t"
            << timeSetup << "\t"
            << timeSolve << "\t";
        if(detailedLog)
        {
            log << std::fixed << std::setprecision(6) << std::setfill('0');
            log << std::setw(10) << timeSetupGH<<  "\t"       <<std::setw(10) << timeSetupRestrProlong << "\t"
                << std::setw(10) <<timeAssGalerkinProjections <<"\t"<<std::setw(10) << timeSetupCoarseSolver <<"\t"
                <<std::setw(10) << timeSetupSmoother << "\t"<< std::setw(10) <<timeSetupAssembler<< "\t";
        }
        log << std::setw(10) << assm.numDofs() << "\t";


        if (compEigs)
            log << "\t" << min << "\t" << essMin << "\t" << max << "\t" << (max/essMin) ;

        log << "\n";
    }

    if (plot && _rank == 0)
    {
        // Plot solution in Paraview
        if(_rank == 0) gsInfo << "Plotting in Paraview: multigrid.pvd.\n";
        gsWriteParaview<>(sol, "multigrid", 1000);
        gsFileManager::open("multigrid.pvd");
    }

    return (resNorm / resNorm0 > tol) ? 1 : 0;
}

#else
#include <gismo.h>
using namespace gismo;
int main(int argc, char **argv)
{
    gsInfo<<" No Mpi enabled! \n";
    return 0;
}

#endif
