/** @file gsMasonryStress_test.cpp

    @brief Provides assembler for the stress calculation of masonry structure.

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
#include <gsSolver/gsSolverUtils.h>
#include <gsSelfSuppSurf/gsAiryStress.h>
#include <gsPde/gsPointLoads.h>


;
using std::vector;
using namespace gismo;

int main(int argc, char *argv[])
{
    /*
    gsFunctionExpr<> func("-1/2*x*y+1");
    gsKnotVector<> kv(-1, 1, 10, 3, 1, 2);
    gsTensorBSplineBasis<2> bsp(kv,kv);
    gsMatrix<> pointsInt = *bsp.anchors();
//     pointsInt.conservativeResize(Eigen::NoChange, 3);
//     pointsInt.col(3) = func.eval(pointsInt);
//     gsGeometry<> * gt = bsp.interpolate(pointsInt); // , *bsp.anchors()
    gsGeometry<> * gt = bsp.interpolate(func.eval(*bsp.anchors()) ); // , *bsp.anchors()
    gt->embed3d();
    gt->degreeElevate(1);
    for (int i = 0; i < 5; ++i)
        gt->uniformRefine();
        
    gsInfo<< *gt <<"\n";
    gsWrite(*gt, "func");
    gsWriteParaview(*gt, "masonry_binear", 1000);
    delete gt;
    return 0; 
*/
    index_t numRefine = 0;
    index_t TypeDirichBoundary = 0;
    std::string str_f("1");
    std::string fn("");    
    gsCmdLine cmd("Testing the Masonry Stress calculator.");
    cmd.addInt   ("r", "refine", "Number of refinement steps"  , numRefine);
    cmd.addString("g", "geometry", "Input 3D surface geometry"  , fn);
    cmd.addInt   ("t", "boundary","To choose an analytical solution as boundary", TypeDirichBoundary);
    cmd.addString("f", "fourth", "external gravity", str_f );    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Read Geometry 
    gsGeometry<>::uPtr g;
    if ( fn.empty() )
        g = gsReadFile<>("surfaces/masonry_bilinear.xml");//"surfaces/masonary_catenary.xml"
    else
        g = gsReadFile<>(fn); 

    // load first patch:
    gsMultiPatch<> patches(*g);
    gsInfo << "Patches: "<< patches <<"\n";   
    //Refinement
    for (int i = 0; i < numRefine; ++i)
        patches.uniformRefine();

    std::vector<gsMatrix<> > zAllCoors(patches.nPatches());
    for(size_t i = 0;i<patches.nPatches();++i)
    {
        //get coordinates in the sequence of local patch
        zAllCoors[i] = patches.patch(i).coefs().col(2) ; // coefs() cannot called by pointer
        patches.patch(i).coefs().col(2).setZero();
    }

//     std::vector<gsMatrix<> > zAllCoors(patches.nPatches()),xyAllCoors(patches.nPatches());
//     gsMultiPatch<> mp2d;
//     for(size_t i = 0;i<patches.nPatches();++i)
//     {
//         zAllCoors[i] = patches.patch(i).coefs().col(2) ; // coefs() cannot called by pointer
//         xyAllCoors[i] = patches.patch(i).coefs().leftCols(2) ;
//         gsGeometry<>::uPtr g=patches.patch(i).basis().makeGeometry(xyAllCoors[i]);
//         mp2d.addPatch(g);
//     }
//     mp2d.computeTopology();    
//     gsDebugVar(zAllCoor);

    if ( patches.geoDim() == 3 ) // if 3D then project to 2D
    {
        for (size_t i = 0; i<patches.nPatches(); ++i)
        {
            //patches.patch(i).coefs().conservativeResize(Eigen::NoChange, 2);
            patches.patch(i).embed(2);
        }
    }
    
    // Source function
    gsFunctionExpr<> f(str_f,2);
        
    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;
    //Dirichlet BCs
    gsFunctionExpr<> bcCatenary("-cosh(x) + cosh(1.0)",2);
    gsFunctionExpr<> bcPinch("0",2);
    gsFunctionExpr<> bcLinear("-x*y",2);
//     gsFunctionExpr<> bcLinear("-x*y+0.5*x*x+0.5*y*y",2);//this also works
    

    for (gsMultiPatch<>::const_biterator 
             bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
    {       
        bcInfo.addCondition( *bit, condition_type::dirichlet, &bcLinear );          
    }  
//     bcInfo.addCondition(0, boundary::west , condition_type::dirichlet, &bcPinch );
    //bcInfo.addCondition(0, boundary::west , condition_type::dirichlet, &bcPinch );
    // Define discretization space by initial refining the basis of the geometry 
//     patches.degreeElevate(1);

    // Initialize Assembler    
    gsAiryStress<real_t> myMasoStress(patches, bcInfo, f, zAllCoors);

    myMasoStress.assemble();

//     gsSparseSolver<>::CGDiagonal solver( myMasoStress.matrix() );
    gsSparseSolver<>::LU solver( myMasoStress.matrix() );
    gsMatrix<real_t> solVector;
    solVector = solver.solve( myMasoStress.rhs() );
    gsMultiPatch<real_t>  AirStress;
    myMasoStress.constructSolution(solVector, AirStress);
    
    gsField<> result(patches, AirStress ); 
    gsField<> expected(patches, bcLinear, false );

    gsInfo<<"Writing the 3D solution in Paraview file. \n";      
    gsWriteParaview( result  , "airy_paraview", 1000);
    gsWriteParaview( expected, "airy_expected", 1000);
    real_t l2error = result.distanceL2(bcLinear);      
    gsInfo << "The L2-error is " << l2error << ".\n";

    // Solution after Newton iteration
    for (size_t i = 0; i<patches.nPatches(); ++i)
    {
        gsGeometry<> & XYZ = patches.patch(i);
        XYZ.coefs().conservativeResize(Eigen::NoChange, 3); // resize coefs to 3D
        XYZ.coefs().col(2) = zAllCoors[i];// fill z-coordinate with solution
    }
    gsField<> resultOnSurface(patches, AirStress);
    gsField<> expectedOnSurface(patches, bcLinear, true);
    gsWrite(patches, "AiryShapeSurface_bspline");
    gsWriteParaview(resultOnSurface, "AiryShapeSurface_paraview", 1000);
    gsWriteParaview(expectedOnSurface, "AiryExpectedSurface_paraview", 1000);
    
    return 0;
}
