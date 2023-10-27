/** @file gsCompositeHBasis_AssemblerExamples.h

    @brief File testing the gsMappedBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsSmoothPatches/gsCompositeAssemblerUtils.h>

#include "gsMSplines/gsMappedBasis.h"
#include "gsMSplines/gsMappedSpline.h"

#include "gsSmoothPatches/gsCompositeBSplineBasis.h"
#include "gsSmoothPatches/gsCompositeIncrSmoothnessGeom.h"
#include <gsSmoothPatches/gsCompositeUtils.h>

#include <gsAssembler/gsAdaptiveRefUtils.h>
#include <gsAssembler/gsErrEstPoissonResidual.h>


;
using std::flush;
using std::vector;
using std::pair;
using std::make_pair;
using std::ios;
using namespace gismo;

//========================================== UTILS ========================================//

void printVector(const std::vector<real_t>& vector)
{
    for(unsigned i = 0;i<vector.size();++i)
        gsInfo << vector[i] << " ";
}

//========================================== TESTS ========================================//

bool solvePoissonCollectResults_stepsWithRefine(
        gsCompositeHBasis<2,real_t> * compBasis,
        gsMappedSpline<2,real_t> * domain,
        gsFunctionExpr<real_t>& g,const gsFunctionExpr<real_t>& dg,gsFunctionExpr<real_t>& f,
        int maxIterations,
        int& plotNr,bool plot=false)
{
    std::vector<real_t> l2errors;
    std::vector<real_t> h1errors;

    // loop over different refinements
    for(int iteration = 0;iteration<maxIterations;++iteration)
    {
        gsMultiPatch<real_t> mp = domain->exportToPatches();

        real_t l2error,h1error,time;
        unsigned dof;
        gsFunctionExpr<real_t> vg_copy=g;
        gsFunctionExpr<real_t> dg_copy=dg;
        gsFunctionWithDerivatives<real_t> g_copy(vg_copy,dg_copy);
        gsFunctionExpr<real_t> f_copy=f;
        // Define Boundary conditions
        gsBoundaryConditions<real_t> bcInfo;
        addAllDirichletBoundaries(mp,g_copy,bcInfo);

        gsMatrix<real_t> reconstructedSol;
        solvePoisson(compBasis,bcInfo,mp,f_copy,time,reconstructedSol);
        std::vector<gsBasis<real_t> *> vc_refine_bases;
        for(size_t i = 0;i<compBasis->nPatches();i++)
        {
            vc_refine_bases.push_back(new gsMappedSingleBasis<2,real_t>(compBasis->getMappedSingleBasis(i)));
        }
        gsMappedSpline<2,real_t> geo(*static_cast<gsMappedBasis<2,real_t> *>(compBasis),reconstructedSol);
        gsMultiPatch<real_t> mp2 = geo.exportToPatches();
        gsField<real_t> sol(mp,mp2);

        if(plot)
        {

            gsField<real_t> exact(mp, g_copy, false);
            gsWriteParaview( sol , "PoissonH"+util::to_string(plotNr)+"solutionComposite"+util::to_string(iteration), compBasis->size()*2, true );
            gsWriteParaview( exact , "PoissonH"+util::to_string(plotNr)+"solutionExact"+util::to_string(iteration), compBasis->size()*2, false );

            gsField<real_t> errorPtr = gsFieldCreator<real_t>::absError(sol,g);
            gsWriteParaview(errorPtr,"PoissonH"+util::to_string(plotNr)+"solutionAbsErrComposite"+util::to_string(iteration), compBasis->size()*2, false);

            gsMultiPatch<> patches = domain->exportToPatches();
            for(size_t i =0;i<compBasis->nPatches();++i)
            {
                gsMesh<real_t> mesh;
                makeMesh<>(compBasis->getBase(i),mesh,10);
                patches.patch(i).evaluateMesh(mesh);
                gsWriteParaview(mesh,"PoissonH"+util::to_string(plotNr)+"finalMesh"+util::to_string(i),true);
            }
            gsInfo << "prints done - ";
        }

        l2error = compositeL2Error(sol,g_copy);
        h1error = compositeH1Error(sol,g_copy);

        dof=compBasis->size();

        // Print out the L2 errors and element sizes
        gsInfo << std::scientific;
        gsInfo << "Refinement step: " <<iteration<<"\n";
        gsInfo << "L2 error: "<<l2error;
        gsInfo << "\nH1 error: "<<h1error;
        gsInfo << "\ndofs: "<<dof;
        gsInfo << "\ntimes: "<<time;
        gsInfo << "\n\n";
        gsInfo.unsetf(ios::fixed | ios::scientific);

        h1errors.push_back(h1error);
        l2errors.push_back(l2error);

        //gsInfo << "test passed\n\n" ;
        gsInfo << flush;
        if(iteration<maxIterations-1)
        {
            // Set up and compute the L2-error to the known exact solution...
            gsNormL2<real_t> norm(sol,g);
            // ...and the error estimate, which needs the right-hand-side.
            gsErrEstPoissonResidual<real_t> errEst(sol,f);

            norm.compute(1);
            errEst.compute(1);

            // Get the vector with element-wise local errors...
            const std::vector<real_t> & elError = norm.elementNorms();
            // ...or the vector with element-wise local error estimates.
            const std::vector<real_t> & elErrEst = errEst.elementNorms();

            // Mark elements for refinement, based on the computed local errors and
            // refCriterion and refParameter.
            std::vector<bool> elMarked( elError.size() );

            MarkingStrategy refCriterion = MarkingStrategy::GARU;//1;
            real_t refParameter = 0.5;//0.5;
            // Use the (in this case known) exact error...
            //gsMarkElementsForRef( elError, refCriterion, refParameter, elMarked);
            // ...or the error estimate.
            gsMarkElementsForRef( elErrEst, refCriterion, refParameter, elMarked);

            gsRefineMarkedElements( *compBasis, elMarked);
            std::vector<gsMatrix<real_t> *> coefs;
            for(size_t i = 0; i<compBasis->nPatches(); ++i)
                coefs.push_back(NULL);
            compBasis->repairPatches(coefs);
            compBasis->updateTopol();
        }
        plotNr++;
    }

    std::vector<real_t> convergenceRatesh1;
    std::vector<real_t> convergenceRatesl2;
    getConvergenceRatios(l2errors,convergenceRatesl2);
    getConvergenceRatios(h1errors,convergenceRatesh1);
    gsInfo << "l2 convergence rates: ";
    printVector(convergenceRatesl2);
    gsInfo << "\nh1 convergence rates: ";
    printVector(convergenceRatesh1);
    gsInfo << "\n\n";
    return true;
}

bool test_poissonSolvingExample1(unsigned switch_var,int degSmooth, int & plotNr,bool plot)
{
    // set up problem
    gsInfo<< "\n------------ poissonSolving ------------\n\n";
    std::string path;
    gsFunctionExpr<real_t> g,dg,f;
    int iterations=2;
    switch(switch_var)
    {
    case 1:
        gsInfo << "triangle domain, degree 2\n";
        path="planar/multipatch_tunnel_thb.xml";
        g=gsFunctionExpr<real_t>("tanh( 18-x-2*y )", 2);
        dg=gsFunctionExpr<real_t>("-(1.0/cosh(18 - x - 2*y))^2",
                                   "-2*(1.0/cosh(18 - x - 2*y))^2",2);
        f=gsFunctionExpr<real_t>("10 * tanh( 18-x-2*y ) * ( 1 - tanh(18-x-2*y )*tanh(18-x-2*y) )",2);
        iterations = 2;
        break;
    }
    gsMultiPatch<> * mp = (  (gsMultiPatch<>::uPtr)gsReadFile<>(path)  ).release();
    if(mp==NULL)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    patchCorner pc(4,4);
    std::vector<patchCorner> cornerList;
    cornerList.push_back(pc);
    gsCompositeIncrSmoothnessGeom<2,real_t> domain(*mp,cornerList,degSmooth);
    typename gsCompositeHBasis<2,real_t>::uPtr basis = memory::convert_ptr<gsCompositeHBasis<2,real_t> >(domain.getCompBasis().clone());

    solvePoissonCollectResults_stepsWithRefine(basis.release(),&domain,g,dg,f,iterations,plotNr,plot);

    gsInfo << "test passed\n\n" ;
    delete mp;
    return true;
}

//========================================== MAIN ========================================//

int main(int argc, char *argv[])
{
    bool passed = true;
    bool plot       = false;

    gsCmdLine cmd("Composite basis tests.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    int plotNr = 0;
    test_poissonSolvingExample1(1,-1,plotNr,plot) && passed;
    return 0;
}
