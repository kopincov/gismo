/** @file gsMasonry_test.cpp

    @brief Provides assembler for the Poisson equation of masonry structure.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Xia, A. Mantzaflaris
*/
#include <iostream>
#include <stdio.h>
#include <fstream>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsSelfSuppSurf/gsMasonry_simple.h>
#include <gsSelfSuppSurf/gsMasonry.h>
#include <gsSelfSuppSurf/gsAiryStress.h>
#include <gsSolver/gsSolverUtils.h>

#include <gsPde/gsPointLoads.h>


#include <gsOptimizer/gsOptMasonry.h>


using std::cout;
using std::endl;
using std::vector;
using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine = 0;     // Lowest number of refinement: numRefine + 1
    //int numIter   = 10;
    std::string str_k1, str_k2, str_k3, str_k4;
    std::string fn("surfaces/masonry_bspline1p.xml");

    gsCmdLine cmd("Testing the Masonry optimizer.");
    cmd.addString("g", "geometry", "Input surface", fn);
    cmd.addInt("r", "uniformrefine", "Uniform refinement steps", numRefine);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Read Geometry
    gsMultiPatch<> patches;
    gsReadFile<>(fn, patches);
    patches.computeTopology();

    gsInfo << "Patches: "<< patches <<"\n";

    gsFunctionExpr<> f("1.0",2);

    // Define Boundary conditions coming from the input surface
    gsBoundaryConditions<> bcInfo;// to do: move inside gsOptMasonry
    gsMultiPatch<> bc = patches;
    for (size_t k = 0; k <  bc.nPatches(); ++k)
    {
        bc.patch(k).coefs().col(0).swap(bc.patch(k).coefs().col(2));
        bc.patch(k).embed(1);
    }

    for (gsMultiPatch<>::const_biterator 
             bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
        bcInfo.addCondition( *bit, condition_type::dirichlet, &bc.patch(bit->patch),0,true );  

    gsPointLoads<real_t> pLoads;

    gsOptMasonry<real_t>   masonry(patches, bcInfo, f, pLoads);

    masonry.solve();

    gsMultiPatch<real_t> finalsurf;

    masonry.constructSurface(finalsurf);

    gsWrite(finalsurf               , "selfsupp_bspline"       );
    gsWrite(masonry.airyStress(), "selfsupp_stress_bspline");

    gsWriteParaview(finalsurf, "selfsupp_paraview", 1000);
    gsWriteParaview(patches, "selfsupp_input_paraview", 1000);
    gsWriteParaview(masonry.airyStress(), "selfsupp_stress_paraview", 1000);

    // Constraints:
    //gsVector<> tmp;
    //masonry.evalCon_into( masonry.currentDesign(), tmp);
    //gsInfo<<"Constr:\n" << tmp.transpose() <<"\n";

    gsInfo<<"Target surf residue:\n" << masonry.targetResidue() <<"\n";

    return 0;
}
