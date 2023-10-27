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
#include <gsSelfSuppSurf/gsMasonry_simple.h>
#include <gsSelfSuppSurf/gsMasonry.h>
#include <gsSolver/gsSolverUtils.h>

#include <gsPde/gsPointLoads.h>
#include <time.h>
#include <math.h>


;
using std::vector;
using namespace gismo;


void computeTargetSurfDofs(const gsMultiPatch<real_t> & mp, 
                           const gsDofMapper & mapper,
                           gsMatrix<real_t> & result)
{  
    result.resize( mapper.freeSize(), 1);
    
    for (size_t k = 0; k <  mp.nPatches(); ++k)
        for (unsigned i = 0; i <  mp.patch(k).coefsSize(); ++i)
        {
            const index_t ii = mapper.index(i,k); 
            if ( mapper.is_free_index(ii) )
                result.row( ii ) = mp.patch(k).coef(i);
        }
}


int main(int argc, char *argv[])
{
    index_t numRefine = 1;     // Lowest number of refinement: numRefine + 1
    index_t numIter   = 10;
    std::string str_k1, str_k2, str_k3, str_k4;
    std::string fn("isss/ex_catenary_4p.xml");

    gsCmdLine cmd("Testing the Masonry constructor.");
    cmd.addInt   ("r", "refine", "Number of refinement steps"  , numRefine);
    cmd.addInt   ("i", "iterations", "Max number of Newton iterations", numIter);
    cmd.addString("g", "geometry", "Input 2D geometry"  , fn);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> xmlData(fn);
    // Check if all needed data exist in the file
    if ( !xmlData.has< gsMultiPatch<> >() )
    {
        gsWarn <<"Did not find a multipatch domain in this file.\n";  
        return 0;
    }

    if ( !xmlData.has<gsFunctionExpr<> >() )
    {
        gsWarn <<"Did not find a function (coefficient) in this file.\n"; 
        return 0;
    }

    if ( !xmlData.has<gsBoundaryConditions<> >() )
    {
        gsWarn <<"Did not find boundary conditions in this this file.\n"; 
        return 0;
    }

    gsMultiPatch<> patches;
    gsFunctionExpr<> f("1",2);
    gsFunctionExpr<> k;
    gsBoundaryConditions<> bcInfo;

    xmlData.getFirst(patches);
    xmlData.getFirst(bcInfo);
    xmlData.getFirst(k);
    
    gsInfo << "Patches: "<< patches <<"\n";

    if ( patches.geoDim() == 3 ) // if 3D then project to 2D
        for (size_t i = 0; i<patches.nPatches(); ++i)
            patches.patch(i).coefs().conservativeResize(Eigen::NoChange, 2);

    //Scaling
    //for (size_t i = 0; i<patches.nPatches(); ++i)
    //    patches.patch(i).coefs().array() /= 10.0;
    //gsWrite(patches,"scaledpatches");
    //

    // Point loads
    gsPointLoads<real_t> pLoads;
    gsVector<> point(2); point<< 0.5, 0.5 ;
    gsVector<> load (3); load << 0.0, 0.0, 1.0 ;
    pLoads.addLoad(point, load, 0 ); 
    pLoads.clear();// no point loads
    

    // Refine the initial geometry   
    patches.degreeElevate(1);
    gsInfo << "Degree of basis wrt the first variable: " << patches.basis(0).degree(0)
              << "\n"
              << "Degree of basis wrt the first variable: " << patches.basis(0).degree(1)
              << "\n";

    for (int i = 0; i < numRefine; ++i)
        patches.uniformRefine();
    gsInfo << "The refinement level is " << numRefine << "\n";

    // Initialize Assembler
    gsMasonry<real_t> myMasonry(patches, bcInfo, f, k, pLoads);
    gsNewtonIterator<real_t> newton(myMasonry);
    newton.setMaxIterations(numIter);

    gsStopwatch watch;
    newton.solve();
    const double t = watch.stop();
    gsInfo<< "The time for newton iteration is " << t << " seconds.\n";

    //gsDebugVar( newton.solution().patch(0).basis().numElements() );
    //gsDebug<< "Solution coefs: " <<
    //    newton.solution().patch(0).coefs().transpose() <<"\n";

    gsMatrix<real_t> curCoefs;
    computeTargetSurfDofs(newton.solution(), myMasonry.system().colMapper(0), curCoefs);
    
    // Verify the residue
    myMasonry.assembleSystem( newton.solution() );
    gsDebug<< "Residue Norm: " << ( myMasonry.matrix() * curCoefs -  myMasonry.rhs() ).norm() <<"\n";

    if ( newton.converged() )
    {
        gsInfo <<"Converged after "<<newton.numIterations()
                <<" iterations with tolerance "<<newton.tolerance() <<".\n";
    }
    else
        gsInfo <<"Newton iteration did not converge.\n";

    // Find the l2 error
    gsField<> sol(patches, newton.solution() );
    //
    const gsFunction<> & exactSol = *bcInfo.dirichletBegin()->function();
    gsInfo << "Using exact solution: "<< exactSol <<"\n";
    real_t l2error = sol.distanceL2(exactSol);
    gsInfo << "The L2-error is " << l2error << ".\n";
    
    // Solution after Newton iteration
    for (size_t i = 0; i<patches.nPatches(); ++i)
    {
        const gsGeometry<> & Z = newton.solution().patch(i);
        
        gsGeometry<> & XYZ = patches.patch(i);
        
        GISMO_ASSERT( XYZ.basis().degree(0)==Z.basis().degree(0) &&
                      XYZ.basis().degree(1)==Z.basis().degree(1),
                      "different degree.");
        
        while ( XYZ.coefs().rows() < Z.coefs().rows() ) // refine coefs to match solution
            XYZ.uniformRefine();
        
        //gsInfo<<"Number of coefficients: "<< Z.coefsSize() <<"\n";
        XYZ.coefs().conservativeResize(Eigen::NoChange, 3); // resize coefs to 3D
        XYZ.coefs().col(2) = Z.coefs();// fill z-coordinate with solution
    }

    gsInfo<<"Writing the 3D solution in Paraview file. \n";
    gsWrite(patches, "masonry_bspline");
    gsWriteParaview( patches, "masonry_paraview", 1000);
        
    return 0;
}
