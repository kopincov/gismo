/**  gsShellMixed_test.cpp

    @brief Test file for gsShellMixedAssembler

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/


#include <gismo.h>
#include <gismo_dev.h>
#include <gsPlateMixed/gsShellMixedAssembler.h>
#include <gsPlateMixed/gsNormL2M.h>
#include <gsPde/gsShellMixedPde.h>

#include <gsPlateMixed/gsFunctionPDisplace.h>
#include <gsPlateMixed/gsFunctionSetPDisplace.h>
#include <gsPlateMixed/gsFunctionPDisplaceModulus.h>
#include <gsPlateMixed/gsFunctionSetPDisplaceModulus.h>
#include <gsPlateMixed/gsFunctionM.h>
#include <gsPlateMixed/gsFunctionSetM.h>


using namespace gismo;
using std::cout;
using std::endl;


void rotate4times(gsVector<real_t,3> axis, gsMultiPatch<>& mp)
{
    gsGeometry<>* patch0 = &mp.patch(0);
    gsGeometry<>::uPtr patch1 = patch0->clone();
    patch1->rotate(EIGEN_PI/2, axis);
    gsGeometry<>::uPtr patch2 = patch1->clone();
    patch2->rotate(EIGEN_PI/2, axis);
    gsGeometry<>::uPtr patch3 = patch2->clone();
    patch3->rotate(EIGEN_PI/2, axis);

    mp.addPatch(give(patch1));
    mp.addPatch(give(patch2));
    mp.addPatch(give(patch3));

    mp.computeTopology();
    //cout<<"whole mp" <<mp.detail()<<endl;
}


int main (int argc, char** args)
{
    int result = 0;

    /// Options ///

    index_t testcase = 7;
    // testcase: 1 thin plate
    //           2 Scordelis-Lo roof (1 patch)
    //           3 Scordelis-Lo roof (4 patches)
    //           4 Pinched cylinder (4 patches)
    //           5 Pinched hemisphere (4 patches)

    //           6 Cylindrical shell strip (1 patch)

    //           7 Hyperbolic paraboloid (1 patch)
    //           8 Clamped cylinder (4 patches)


    bool plot = false;
    bool writeFileConvergence = false;
    bool writeFileLocking = false;
    bool writeFileConvergenceNorm = false;

    double mpDeformedScaling = 1;
    bool uniformRefinement = false;
    bool useBSplineApproxGeometry = true; // for plotting deformed goemetry set useBSplineApproxGeometry = true

    bool analytSol = false;
    bool computeFineSolution = false;

    bool outputEvals = false;
    bool outputMatrix = false;
    bool outputBlock = false;

    //int numIncreaseDegree =0;
    index_t numElevateDegree = 0;
    index_t numHRefine = 5;
    index_t numInsertKnots = 8; // knots -3 -elevate
    int numElevateDegreeFine = 3;
    int numHRefineFine = 6;

    int numElevateDegreeGeometry = 3; // 3
    int numHRefineGeometry = 4; // 4

    index_t slenderness = 100; // R/t
    real_t thickness     = 0.05;
    real_t youngsModulus = 1.2E6;
    real_t poissonsRatio = 0.0;

    /// Read command line arguments ///
    gsCmdLine cmd("gsShellMixed_test");
    cmd.addInt("t","testcase","Choosen testcase", testcase);
    cmd.addInt("r","refine","Number of h-refinements", numHRefine);
    cmd.addInt("k","knots","Number of inserted knots", numInsertKnots);
    cmd.addInt("e","elevate","Maximal degree elevation", numElevateDegree);
    cmd.addInt("s","slenderness","Slenderness of the plate", slenderness);

    try { cmd.getValues(argc,args); } catch (int rv) { return rv; }


    /// Set options ///

    gsOptionList options = gsShellMixedAssembler<>::defaultOptions();
    options.setInt("DirichletValues"  , dirichlet::homogeneous);
    options.setInt("DirichletStrategy", dirichlet::elimination );
    options.setInt("InterfaceStrategy", iFace    ::conforming  );
    if(numElevateDegree >= 0)
        options.setReal("quA", 1.0);
    else // p=1
        options.setReal("quA", 2.0);
    options.setInt ("quB", 1  );
    options.setReal("bdA", 2.0);
    options.setInt ("bdB", 1  );

    //options.setInt("phiBcMethodSs", gsShellMixedAssembler<>::nitsche);
    //options.setInt("phiBcMethodF", gsShellMixedAssembler<>::nitsche);
    options.setInt("phiBcMethodSs", gsShellMixedAssembler<>::lagrange);
    options.setInt("phiBcMethodF", gsShellMixedAssembler<>::lagrange);

    options.setSwitch("wHomogenPsiq", false);
    options.setSwitch("wPenaltyPsiq", false);
    options.setSwitch("solveWholeSystem", true);
    options.setSwitch("pBoundary0", false);

    options.setSwitch("Mmixed", true);
    options.setSwitch("Nmixed", false);
    options.setSwitch("lambdaTMean0", false);

    /// Variables ///

    gsMultiPatch<> sol_p, sol_phi, sol_u, solFine_u;
    gsMatrix<> solVector, solVectorFine;

    gsStopwatch watch;
    gsMultiPatch<> mp, mpBSpline;

    gsFunctionExpr<> analytSol_u;
    gsFunctionExpr<> surfForceFun;


    /// Setup geometry ///

    if (testcase == 1)
    {
        gsReadFile<>("surfaces/sphere_u6p.xml", mp);
        //gsReadFile<>("thin_plate2.xml", mp);
    }
    else if (testcase == 2 || testcase == 3)
    {
        thickness = 0.25;
        youngsModulus = 4.32E8;
        poissonsRatio = 0.0;
        gsReadFile<>("scordelis_lo_roof.xml", mp);

        // deg=2 both directions
        mp.patch(0).degreeElevate(1,0);

        mpDeformedScaling = 15;
    }
    else if (testcase == 4)
    {
        thickness = 3;
        youngsModulus = 3E6;
        poissonsRatio = 0.3;
        gsReadFile<>("pinched_cylinder.xml", mp);

        // deg=2 both directions
        mp.patch(0).degreeElevate(1,0);

        mpDeformedScaling = 3e6;

    }
    else if (testcase == 5)
    {
        thickness = 0.04;
        youngsModulus = 6.825E7;
        poissonsRatio = 0.3;
        gsReadFile<>("quarter_hemisphere.xml", mp);

        mpDeformedScaling = 20;
    }
    else if (testcase == 6)
    {

        thickness = 10./slenderness; // t = R / slenderness, slenderness = R / t
        youngsModulus = 1000;
        poissonsRatio = 0;
        gsReadFile<>("cylinder_strip.xml", mp);

        // deg=2 both directions
        mp.patch(0).degreeElevate(1,0);

        mpDeformedScaling = 1;
    }
    else if (testcase == 7)
    {

        thickness = 1./slenderness; // t = L / slenderness, slenderness = L / t
        youngsModulus = 2E11;
        poissonsRatio = 0.3;
        gsReadFile<>("hyperbolic_paraboloid.xml", mp);

        mpDeformedScaling = 1;
    }
    else if (testcase == 8)
    {

        thickness = 1./slenderness; // t = L / slenderness, slenderness = L / t
        youngsModulus = 2E11; // 2E11
        poissonsRatio = 1./3;
        gsReadFile<>("bathe_cylinder.xml", mp);

        // deg=2 both directions
        mp.patch(0).degreeElevate(1,0);

        mpDeformedScaling = 1e-1;
    }

/*
    if (plot)
    {
        gsWriteParaview<>( mp, "geometry", 1000, true, true);

        // Run paraview on exit
        result = system("paraview geometry.pvd &");
    }
*/


    /// Refinement h and p of geometry ///
    /*
    // Elevate degree (= Increase degree if no h-refinement has been done before, i.e. only knots at the boundary)
    mp.patch(0).degreeElevate(numElevateDegree);

    // h-refine
    for(index_t i = 0; i< numHRefine; ++i)
        mp.patch(0).uniformRefine();

    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";
    */


    /// Approximate Nurbs geometry with BSpline geometry ///

    // Approximate NURBS geometry with B-Spline geometry (only for 1 patch)
    gsMultiPatch<> mpNurbs(mp);

    const gsTensorNurbsBasis<2, real_t>* basisNurbs_temp = dynamic_cast< const gsTensorNurbsBasis<2, real_t>* >(&mp.basis(0));
    if(basisNurbs_temp != NULL)
    {
        cout<<"Geometry is Nurbs geometry"<<endl;
        gsTensorBSplineBasis<2> basisBSpline_temp(basisNurbs_temp->source());


        // Store unrefined mp for approximation basis
        gsGeometry<>::uPtr BSpline_temp = basisBSpline_temp.interpolateAtAnchors(mpNurbs.piece(0).eval(basisBSpline_temp.anchors()));
        mpBSpline =  *BSpline_temp;

        // Refine h and p of geometry (for good approximation of NURBS geometry)
        basisBSpline_temp.degreeElevate(numElevateDegreeGeometry);

        for(index_t i = 0; i< numHRefineGeometry; ++i)
            basisBSpline_temp.uniformRefine();

        BSpline_temp = basisBSpline_temp.interpolateAtAnchors(mpNurbs.piece(0).eval(basisBSpline_temp.anchors()));

        if(useBSplineApproxGeometry)
        {
            cout<<"Given Nurbs geometry is approximated with BSpline geometry"<<endl;
            mp =  *BSpline_temp;
        }
    }
    else
    {
        cout<<"Geometry is BSpline geometry"<<endl;
        mpBSpline = mp;
    }


    gsInfo<<"Geometry basis (patch 0): "<< mp.patch(0).basis() << "\n";



    /// Build whole geometry ///

    // Pinched cylinder (4 patches)
    if(testcase == 4 || testcase == 8)
    {
        gsVector<real_t,3> axis;
        axis<<1, 0, 0;

        rotate4times(axis, mp);
        rotate4times(axis, mpBSpline);
    }

    // Pinched hemisphere (4 patches)
    if(testcase == 5)
    {

        gsVector<real_t,3> axis;
        axis<<0, 0, 1;

        rotate4times(axis, mp);
        rotate4times(axis, mpBSpline);

    }



    /// Split geometry into patches ///

    // Scordelis-Lo roof (4 patches)
    if(testcase == 3)
    {
        mp = mp.uniformSplit(); // splits also parameter space !!
        mpBSpline = mpBSpline.uniformSplit();
    }

    gsInfo<<"Geometry " << mp <<"\n";


    /// Boundary conditions and surface forces ///

    enum bcType_plate
    {
        clamped = 0,
        simplysupp   = 1,
        free     = 2
    };

    gsBoundaryConditions<> BCs;
    std::vector<std::vector<gsShellMixedPde<>::corner_type> > corners(mp.nPatches());
    std::vector<std::vector<gsMatrix<> >  > cornerCouplingCoefs(mp.nPatches());
    std::vector<std::vector<int >  > indFreeComp(mp.nPatches());

    gsVector<> tmp(3);
    gsVector<> tmp_lineLoad(3);
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();


    // Scordelis-Lo roof
    if(testcase == 2)
    {

        // u
        // Diaphragm conditions
        BCs.add(0, boundary::west, "Dirichlet", 0, 7 ); // u2
        BCs.add(0, boundary::west, "Dirichlet", 0, 8 ); // u3
        BCs.addCornerValue(boundary::southwest, 0.0, 0, 6); // (corner,value, patch, unknown) // u1

        BCs.add(0, boundary::east, "Dirichlet", 0, 7 ); // u2
        BCs.add(0, boundary::east, "Dirichlet", 0, 8 ); // u3

        // north clamped
        //BCs.add(0, boundary::north, "Dirichlet", 0, 6 ); // u1
        //BCs.add(0, boundary::north, "Dirichlet", 0, 7 ); // u2
        //BCs.add(0, boundary::north, "Dirichlet", 0, 8 ); // u3

        // p (same as u3)
        BCs.add(0, boundary::west, "Dirichlet", 0, 0 );
        BCs.add(0, boundary::east, "Dirichlet", 0, 0 );
        //BCs.add(0, boundary::north, "Dirichlet", 0, 0 ); // north clamped

        // phi
        BCs.add(0, boundary::north, "Free", 0, 1 ); // north clamped
        BCs.add(0, boundary::south, "Free", 0, 1 );

        BCs.add(0, boundary::west, "SimplySupp", 0, 1 );
        BCs.add(0, boundary::east, "SimplySupp", 0, 1 );


        // phi corner values BCs to eliminate kernel
        BCs.addCornerValue(boundary::northeast, 0, 0, 1);
        BCs.addCornerValue(boundary::northeast, 0, 0, 2);

        BCs.addCornerValue(boundary::southwest, 0, 0, 1);
        //BCs.addCornerValue(boundary::southeast, 0, 0, 2);

        //BCs.addCornerValue(boundary::northwest, 0, 0, 1);
        //BCs.addCornerValue(boundary::northwest, 0, 0, 2);



        // lagrange multipliers
        // lambdaN
        // eliminate Clamped
        //BCs.add(0, boundary::north, "Dirichlet", 0, 9 ); // north clamped

        // lambdaT
        // eliminate Clamped
        //BCs.add(0, boundary::north, "Dirichlet", 0, 10 ); // north clamped
        // eliminate SimplySupp
        BCs.add(0, boundary::west, "Dirichlet", 0, 10 );
        BCs.add(0, boundary::east, "Dirichlet", 0, 10 );

        // eliminate corners
        BCs.addCornerValue(boundary::northeast, 0, 0, 10);
        BCs.addCornerValue(boundary::northwest, 0, 0, 10);
        BCs.addCornerValue(boundary::southeast, 0, 0, 10);
        BCs.addCornerValue(boundary::southwest, 0, 0, 10);

        // eliminate corners
        //BCs.addCornerValue(boundary::northeast, 0, 0, 9);
        //BCs.addCornerValue(boundary::northwest, 0, 0, 9);
        //BCs.addCornerValue(boundary::southeast, 0, 0, 9);
        //BCs.addCornerValue(boundary::southwest, 0, 0, 9);


        // set corner coupling for lagrange multipliers
        corners[0].resize(5);
        cornerCouplingCoefs[0].resize(5);

        // N S f, E W ss
        corners[0][1]= gsShellMixedPde<>::sf;
        corners[0][2]= gsShellMixedPde<>::sf;
        corners[0][3]= gsShellMixedPde<>::sf;
        corners[0][4]= gsShellMixedPde<>::sf;
        //corners[0][3]= gsShellMixedPde<>::none;
        //corners[0][4]= gsShellMixedPde<>::none;

        cornerCouplingCoefs[0][1].resize(1,4);
        cornerCouplingCoefs[0][1]<< 0, 1, -1, 0;
        cornerCouplingCoefs[0][2].resize(1,4);
        cornerCouplingCoefs[0][2]<< 1, 0, 0, 1;
        cornerCouplingCoefs[0][3].resize(1,4);
        cornerCouplingCoefs[0][3]<< 1, 0, 0, 1;
        cornerCouplingCoefs[0][4].resize(1,4);
        cornerCouplingCoefs[0][4]<< 0, 1, -1, 0;

        // set index free components for langrange multipliers for mean(lambdaT) = 0
        indFreeComp[0].resize(5);
        indFreeComp[0][boundary::north] = 0; // north clamped
        indFreeComp[0][boundary::south] = 1;
        indFreeComp[0][boundary::west]  = -1;
        indFreeComp[0][boundary::east]  = -1;



        // surface forces
        tmp << 0, 0, -90;
        //tmp.setZero();

        // Point loads
        //gsVector<> point(2); point<< 0.5, 0.5 ;
        //gsVector<> load (3); load << 0.0, 0.0, -90 ;
        //pLoads.addLoad(point, load, 0 );


        /*
        // whole boundary simply supp
        for(int i=1; i!=5 ; ++i)
        {
            BCs.add(0, i, "Dirichlet", 0, 3 ); // u1
            BCs.add(0, i, "Dirichlet", 0, 4 ); // u2
            BCs.add(0, i, "Dirichlet", 0, 5 ); // u3

            BCs.add(0, i, "Dirichlet", 0, 0 ); //p

            BCs.add(0, i, "SimplySupp", 0, 1 ); // phi

            BCs.add(0, i, "Dirichlet", 0, 7 ); // lambdat
       }

        // eliminate corners lambdan
        BCs.addCornerValue(boundary::northeast, 0, 0, 6);
        BCs.addCornerValue(boundary::northwest, 0, 0, 6);
        BCs.addCornerValue(boundary::southeast, 0, 0, 6);
        BCs.addCornerValue(boundary::southwest, 0, 0, 6);
        */


        /*
        // whole boundary clamped
        for(int i=1; i!=5 ; ++i)
        {
            BCs.add(0, i, "Dirichlet", 0, 6 ); // u1
            BCs.add(0, i, "Dirichlet", 0, 7 ); // u2
            BCs.add(0, i, "Dirichlet", 0, 8 ); // u3

            BCs.add(0, i, "Dirichlet", 0, 0 ); //p
       }
       */



    }
    // Scordelis-Lo roof 4 patches
    if(testcase == 3)
    {
        // patch 0
        BCs.add(0, boundary::west, "Dirichlet", 0, 7 );  // u2
        BCs.add(0, boundary::west, "Dirichlet", 0, 8 );  // u3
        BCs.addCornerValue(boundary::southwest, 0.0, 0, 6); // (corner,value, patch, unknown) // u1
        BCs.add(0, boundary::west, "Dirichlet", 0, 0 );  // p
        BCs.add(0, boundary::west, "SimplySupp", 0, 1 ); // phi
        BCs.add(0, boundary::west, "Dirichlet", 0, 10 );  // lambdat

        BCs.add(0, boundary::south, "Free", 0, 1 ); // phi

        // set corner coupling for lagrange multipliers
        corners[0].resize(5);
        cornerCouplingCoefs[0].resize(5);

        corners[0][1]= gsShellMixedPde<>::sf;
        corners[0][2]= gsShellMixedPde<>::none;
        corners[0][3]= gsShellMixedPde<>::none;
        corners[0][4]= gsShellMixedPde<>::none;
        cornerCouplingCoefs[0][1].resize(1,4);
        cornerCouplingCoefs[0][1]<< 0, 1, -1, 0;

        // patch 1
        BCs.add(1, boundary::west, "Dirichlet", 0, 7 );  // u2
        BCs.add(1, boundary::west, "Dirichlet", 0, 8 );  // u3
        BCs.add(1, boundary::west, "Dirichlet", 0, 0 );  // p
        BCs.add(1, boundary::west, "SimplySupp", 0, 1 ); // phi
        BCs.add(1, boundary::west, "Dirichlet", 0, 10 );  // lambdat

        BCs.add(1, boundary::north, "Free", 0, 1 ); // phi

        // set corner coupling for lagrange multipliers
        corners[1].resize(5);
        cornerCouplingCoefs[1].resize(5);

        corners[1][1]= gsShellMixedPde<>::none;
        corners[1][2]= gsShellMixedPde<>::none;
        corners[1][3]= gsShellMixedPde<>::sf;
        corners[1][4]= gsShellMixedPde<>::none;
        cornerCouplingCoefs[1][3].resize(1,4);
        cornerCouplingCoefs[1][3]<< 1, 0, 0, 1;

        // patch 2
        BCs.add(2, boundary::east, "Dirichlet", 0, 7 );  // u2
        BCs.add(2, boundary::east, "Dirichlet", 0, 8 );  // u3
        BCs.add(2, boundary::east, "Dirichlet", 0, 0 );  // p
        BCs.add(2, boundary::east, "SimplySupp", 0, 1 ); // phi
        BCs.add(2, boundary::east, "Dirichlet", 0, 10 );  // lambdat

        BCs.add(2, boundary::south, "Free", 0, 1 ); // phi

        // set corner coupling for lagrange multipliers
        corners[2].resize(5);
        cornerCouplingCoefs[2].resize(5);

        corners[2][1]= gsShellMixedPde<>::none;
        corners[2][2]= gsShellMixedPde<>::sf;
        corners[2][3]= gsShellMixedPde<>::none;
        corners[2][4]= gsShellMixedPde<>::none;
        cornerCouplingCoefs[2][2].resize(1,4);
        cornerCouplingCoefs[2][2]<< 1, 0, 0, 1;

        // patch 3
        BCs.add(3, boundary::east, "Dirichlet", 0, 7 );  // u2
        BCs.add(3, boundary::east, "Dirichlet", 0, 8 );  // u3
        BCs.add(3, boundary::east, "Dirichlet", 0, 0 );  // p
        BCs.add(3, boundary::east, "SimplySupp", 0, 1 ); // phi
        BCs.add(3, boundary::east, "Dirichlet", 0, 10 );  // lambdat

        BCs.add(3, boundary::north, "Free", 0, 1 ); // phi

        // set corner coupling for lagrange multipliers
        corners[3].resize(5);
        cornerCouplingCoefs[3].resize(5);

        corners[3][1]= gsShellMixedPde<>::none;
        corners[3][2]= gsShellMixedPde<>::none;
        corners[3][3]= gsShellMixedPde<>::none;
        corners[3][4]= gsShellMixedPde<>::sf;
        cornerCouplingCoefs[3][4].resize(1,4);
        cornerCouplingCoefs[3][4]<< 0, 1, -1, 0;



        /*
        // whole boundary clamped
        // patch 0
        int boundarySide = boundary::west;
        BCs.add(0, boundarySide, "Dirichlet", 0, 3 ); // u1
        BCs.add(0, boundarySide, "Dirichlet", 0, 4 ); // u2
        BCs.add(0, boundarySide, "Dirichlet", 0, 5 ); // u3
        BCs.add(0, boundarySide, "Dirichlet", 0, 0 ); //p

        boundarySide = boundary::south;
        BCs.add(0, boundarySide, "Dirichlet", 0, 3 ); // u1
        BCs.add(0, boundarySide, "Dirichlet", 0, 4 ); // u2
        BCs.add(0, boundarySide, "Dirichlet", 0, 5 ); // u3
        BCs.add(0, boundarySide, "Dirichlet", 0, 0 ); //p

        // patch 1
        boundarySide = boundary::west;
        BCs.add(1, boundarySide, "Dirichlet", 0, 3 ); // u1
        BCs.add(1, boundarySide, "Dirichlet", 0, 4 ); // u2
        BCs.add(1, boundarySide, "Dirichlet", 0, 5 ); // u3
        BCs.add(1, boundarySide, "Dirichlet", 0, 0 ); //p

        boundarySide = boundary::north;
        BCs.add(1, boundarySide, "Dirichlet", 0, 3 ); // u1
        BCs.add(1, boundarySide, "Dirichlet", 0, 4 ); // u2
        BCs.add(1, boundarySide, "Dirichlet", 0, 5 ); // u3
        BCs.add(1, boundarySide, "Dirichlet", 0, 0 ); //p

        // patch 2
        boundarySide = boundary::south;
        BCs.add(2, boundarySide, "Dirichlet", 0, 3 ); // u1
        BCs.add(2, boundarySide, "Dirichlet", 0, 4 ); // u2
        BCs.add(2, boundarySide, "Dirichlet", 0, 5 ); // u3
        BCs.add(2, boundarySide, "Dirichlet", 0, 0 ); //p

        boundarySide = boundary::east;
        BCs.add(2, boundarySide, "Dirichlet", 0, 3 ); // u1
        BCs.add(2, boundarySide, "Dirichlet", 0, 4 ); // u2
        BCs.add(2, boundarySide, "Dirichlet", 0, 5 ); // u3
        BCs.add(2, boundarySide, "Dirichlet", 0, 0 ); //p

        // patch 3
        boundarySide = boundary::north;
        BCs.add(3, boundarySide, "Dirichlet", 0, 3 ); // u1
        BCs.add(3, boundarySide, "Dirichlet", 0, 4 ); // u2
        BCs.add(3, boundarySide, "Dirichlet", 0, 5 ); // u3
        BCs.add(3, boundarySide, "Dirichlet", 0, 0 ); //p

        boundarySide = boundary::east;
        BCs.add(3, boundarySide, "Dirichlet", 0, 3 ); // u1
        BCs.add(3, boundarySide, "Dirichlet", 0, 4 ); // u2
        BCs.add(3, boundarySide, "Dirichlet", 0, 5 ); // u3
        BCs.add(3, boundarySide, "Dirichlet", 0, 0 ); //p
        */

        // surface forces
        tmp << 0, 0, -90;

    }

    // Pinched cylinder
    if(testcase == 4)
    {
        for(int i =0; i<4; ++i)
        {
            // u
            // Diaphragm conditions
            BCs.add(i, boundary::west, "Dirichlet", 0, 7 ); // u2
            BCs.add(i, boundary::west, "Dirichlet", 0, 8 ); // u3
            //BCs.addCornerValue(boundary::southwest, 0.0, 0, 6); // (corner,value, patch, unknown) // u1

            BCs.add(i, boundary::east, "Dirichlet", 0, 7 ); // u2
            BCs.add(i, boundary::east, "Dirichlet", 0, 8 ); // u3

            // p (same as u3)
            BCs.add(i, boundary::west, "Dirichlet", 0, 0 );
            BCs.add(i, boundary::east, "Dirichlet", 0, 0 );

            // phi
            BCs.add(i, boundary::west, "SimplySupp", 0, 1 );
            BCs.add(i, boundary::east, "SimplySupp", 0, 1 );

            // phi corner values BCs to eliminate kernel
            BCs.addCornerValue(boundary::northeast, 0, 0, 1);
            BCs.addCornerValue(boundary::northeast, 0, 0, 2);
            BCs.addCornerValue(boundary::southwest, 0, 0, 1);

            // lagrange multipliers
            // lambdaT
            // eliminate SimplySupp
            BCs.add(i, boundary::west, "Dirichlet", 0, 10 );
            BCs.add(i, boundary::east, "Dirichlet", 0, 10 );

            // no corner coupling for lagrange multipliers
            corners[i].resize(5);
            cornerCouplingCoefs[i].resize(5);

            corners[i][1]= gsShellMixedPde<>::none;
            corners[i][2]= gsShellMixedPde<>::none;
            corners[i][3]= gsShellMixedPde<>::none;
            corners[i][4]= gsShellMixedPde<>::none;
        }

        /*
        // Symmetry in y-direction
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 4 ); // u2 = uy

        // Symmetry in x-direction
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 3 ); // u1 = ux

        // Symmetry in z-direction
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 5 ); // u3 = uz
        */

        // Surface forces
        tmp.setZero();

        // Point loads
        gsVector<> point(2); point<< 0.5, 1;
        gsVector<> load (3); load << 0.0, 0.0, -1 ;
        pLoads.addLoad(point, load, 0 );

        point<< 0.5, 1 ;
        load << 0.0, 0.0, 1 ;
        pLoads.addLoad(point, load, 2 );

    }

    //Pinched hemisphere
    if(testcase == 5)
    {
        for(int i =0; i<4; ++i)
        {
            // u
            // fix top
            BCs.add(i, boundary::north, "Dirichlet", 0, 6 ); // u1
            BCs.add(i, boundary::north, "Dirichlet", 0, 7 ); // u2
            BCs.add(i, boundary::north, "Dirichlet", 0, 8 ); // u3

            // phi
            BCs.add(i, boundary::south, "Free", 0, 1 );

            // phi corner values BCs to eliminate kernel
            BCs.addCornerValue(boundary::northeast, 0, 0, 1);
            BCs.addCornerValue(boundary::northeast, 0, 0, 2);
            BCs.addCornerValue(boundary::southwest, 0, 0, 1);

            // lagrange multipliers
            // no corner coupling for lagrange multipliers
            corners[i].resize(5);
            cornerCouplingCoefs[i].resize(5);

            corners[i][1]= gsShellMixedPde<>::none;
            corners[i][2]= gsShellMixedPde<>::none;
            corners[i][3]= gsShellMixedPde<>::none;
            corners[i][4]= gsShellMixedPde<>::none;

            // set index free components for langrange multipliers for mean(lambdaT) = 0
            indFreeComp[i].resize(5);
            indFreeComp[i][boundary::south] = 0;
            indFreeComp[i][boundary::north] = -1;
            indFreeComp[i][boundary::west]  = -1;
            indFreeComp[i][boundary::east]  = -1;
        }

        // Surface forces
        tmp.setZero();

        // Point loads
        gsVector<> point(2); point<< 0, 0;
        gsVector<> load (3); load << -2, 0, 0 ;
        pLoads.addLoad(point, load, 0 );

        point<< 0, 0 ;
        load << 0, 2, 0 ;
        pLoads.addLoad(point, load, 1 );

        point<< 0, 0 ;
        load << 2, 0, 0 ;
        pLoads.addLoad(point, load, 2 );

        point<< 0, 0 ;
        load << 0, -2, 0 ;
        pLoads.addLoad(point, load, 3 );

    }

    // Cylindrical shell strip
    if(testcase == 6)
    {

        // u
        // clamped
        BCs.add(0, boundary::north, "Dirichlet", 0, 6 ); // u1
        BCs.add(0, boundary::north, "Dirichlet", 0, 7 ); // u2
        BCs.add(0, boundary::north, "Dirichlet", 0, 8 ); // u3

        // p (same as u3)
        BCs.add(0, boundary::north, "Dirichlet", 0, 0 );

        // phi
        BCs.add(0, boundary::west, "Free", 0, 1 );
        BCs.add(0, boundary::east, "Free", 0, 1 );
        BCs.add(0, boundary::south, "Free", 0, 1 );

        // phi corner values BCs to eliminate kernel
        BCs.addCornerValue(boundary::northeast, 0, 0, 1);
        BCs.addCornerValue(boundary::northeast, 0, 0, 2);
        BCs.addCornerValue(boundary::southwest, 0, 0, 1);


        // lagrange multipliers
        // lambdaN
        // eliminate Clamped
        BCs.add(0, boundary::north, "Dirichlet", 0, 9 );

        // lambdaT
        // eliminate Clamped
        BCs.add(0, boundary::north, "Dirichlet", 0, 10 );


        // set corner coupling for lagrange multipliers
        corners[0].resize(5);
        cornerCouplingCoefs[0].resize(5);

        // N c, E S W f
        corners[0][1]= gsShellMixedPde<>::ff;
        corners[0][2]= gsShellMixedPde<>::ff;
        corners[0][3]= gsShellMixedPde<>::none;
        corners[0][4]= gsShellMixedPde<>::none;
        cornerCouplingCoefs[0][1].resize(2,4);
        cornerCouplingCoefs[0][1]<< 1, 0, 0, 1, 0, 1, -1, 0;
        cornerCouplingCoefs[0][2].resize(2,4);
        cornerCouplingCoefs[0][2]<< 1, 0, 0, 1, 0, 1, -1, 0;

        // set index free components for langrange multipliers for mean(lambdaT) = 0
        indFreeComp[0].resize(5);
        indFreeComp[0][boundary::north] = -1; // north clamped
        indFreeComp[0][boundary::south] = 0;
        indFreeComp[0][boundary::west]  = 0;
        indFreeComp[0][boundary::east]  = 0;

        // Surface forces
        tmp.setZero();

        // Line loads
        tmp_lineLoad << 0, 0.1*pow(thickness,3), 0;
        //tmp_lineLoad << 0, 0, 0;
        gsConstantFunction<> lineLoad(tmp_lineLoad,3);
        BCs.add(0, boundary::south, "LineLoad", lineLoad, 0 );

    }

    // Hyperbolic paraboloid
    if(testcase == 7)
    {

        // u
        // clamped
        BCs.add(0, boundary::west, "Dirichlet", 0, 6 ); // u1
        BCs.add(0, boundary::west, "Dirichlet", 0, 7 ); // u2
        BCs.add(0, boundary::west, "Dirichlet", 0, 8 ); // u3

        // p (same as u3)
        BCs.add(0, boundary::west, "Dirichlet", 0, 0 );

        // phi
        BCs.add(0, boundary::north, "Free", 0, 1 );
        BCs.add(0, boundary::east, "Free", 0, 1 );
        BCs.add(0, boundary::south, "Free", 0, 1 );

        // phi corner values BCs to eliminate kernel
        BCs.addCornerValue(boundary::northeast, 0, 0, 1);
        BCs.addCornerValue(boundary::northeast, 0, 0, 2);
        BCs.addCornerValue(boundary::southwest, 0, 0, 1);


        // lagrange multipliers
        // lambdaN
        // eliminate Clamped
        BCs.add(0, boundary::west, "Dirichlet", 0, 9 );

        // lambdaT
        // eliminate Clamped
        BCs.add(0, boundary::west, "Dirichlet", 0, 10 );


        // set corner coupling for lagrange multipliers
        corners[0].resize(5);
        cornerCouplingCoefs[0].resize(5);

        // N c, E S W f
        corners[0][1]= gsShellMixedPde<>::none;
        corners[0][2]= gsShellMixedPde<>::ff;
        corners[0][3]= gsShellMixedPde<>::none;
        corners[0][4]= gsShellMixedPde<>::ff;
        cornerCouplingCoefs[0][2].resize(2,4);
        cornerCouplingCoefs[0][2]<< 1, 0, 0, 1, 0, 1, -1, 0;
        cornerCouplingCoefs[0][4].resize(2,4);
        cornerCouplingCoefs[0][4]<< 1, 0, 0, 1, 0, 1, -1, 0;

        // set index free components for langrange multipliers for mean(lambdaT) = 0
        indFreeComp[0].resize(5);
        indFreeComp[0][boundary::north] = 0;
        indFreeComp[0][boundary::south] = 0;
        indFreeComp[0][boundary::west]  = -1; // west clamped
        indFreeComp[0][boundary::east]  = 0;

        // Surface forces
        tmp << 0, 0, -8000 * thickness;
        //tmp.setZero();

    }

    // clamped cylinder
    if(testcase == 8)
    {
        for(int i =0; i<4; ++i)
        {
            // u
            // Diaphragm conditions
            BCs.add(i, boundary::west, "Dirichlet", 0, 6 ); // u1
            BCs.add(i, boundary::west, "Dirichlet", 0, 7 ); // u2
            BCs.add(i, boundary::west, "Dirichlet", 0, 8 ); // u3
            //BCs.addCornerValue(boundary::southwest, 0.0, 0, 6); // (corner,value, patch, unknown) // u1

            BCs.add(i, boundary::east, "Dirichlet", 0, 6 ); // u1
            BCs.add(i, boundary::east, "Dirichlet", 0, 7 ); // u2
            BCs.add(i, boundary::east, "Dirichlet", 0, 8 ); // u3

            // p (same as u3)
            BCs.add(i, boundary::west, "Dirichlet", 0, 0 );
            BCs.add(i, boundary::east, "Dirichlet", 0, 0 );

            // phi corner values BCs to eliminate kernel
            BCs.addCornerValue(boundary::northeast, 0, 0, 1);
            BCs.addCornerValue(boundary::northeast, 0, 0, 2);
            BCs.addCornerValue(boundary::southwest, 0, 0, 1);

            // lagrange multipliers
            // lambdaN
            // eliminate Clamped
            BCs.add(i, boundary::west, "Dirichlet", 0, 9 );
            BCs.add(i, boundary::east, "Dirichlet", 0, 9 );

            // lambdaT
            // eliminate Clamped
            BCs.add(i, boundary::west, "Dirichlet", 0, 10 );
            BCs.add(i, boundary::east, "Dirichlet", 0, 10 );


            // no corner coupling for lagrange multipliers
            corners[i].resize(5);
            cornerCouplingCoefs[i].resize(5);

            corners[i][1]= gsShellMixedPde<>::none;
            corners[i][2]= gsShellMixedPde<>::none;
            corners[i][3]= gsShellMixedPde<>::none;
            corners[i][4]= gsShellMixedPde<>::none;
        }


        // Surface forces
        tmp.setZero();
        //surfForceFun = gsFunctionExpr<>("v:="+util::to_string(t)+"; 2*Cos(x^2+y)^2-1 * -x","-y","-z",3);
        surfForceFun = gsFunctionExpr<>("0","v:="+util::to_string(thickness)+";v * 10^6 * (2*y^2-1) * -y","v:="+util::to_string(thickness)+";v * 10^6 * (2*y^2-1) * -z",3);


    }

    /// ShellMixedPde ///

    gsConstantFunction<> surfForceConstant(tmp,3);
    gsFunction<>* surfForce;

    if(testcase != 8)
        surfForce = &surfForceConstant;
    else
        surfForce = &surfForceFun;

    gsShellMixedPde<> pde(mp, BCs, *surfForce, pLoads, youngsModulus, poissonsRatio, thickness, corners, cornerCouplingCoefs, indFreeComp);



    /// Store basis for assembling ///

    // Copy basis from the geometry

    gsMultiBasis<> basis(mpBSpline);
    gsMultiBasis<> basisFine(mpBSpline);
    gsMultiPatch<> mpFine(mp);

    // revert degree elevation for approximation
    //basis.degreeReduce(numElevateDegreeGeometry);

    // Refinement h and p of basis

    // elevate degree (= increase degree if no h-refinement has been done before, i.e. only knots at the boundary)
    if(numElevateDegree >= 0)
    {
        basis.degreeElevate(numElevateDegree);
        basisFine.degreeElevate(numElevateDegreeFine);
        mpFine.degreeElevate(numElevateDegreeFine);
    }
    else
    {
        basis.degreeReduce(-numElevateDegree);
    }

    gsMultiBasis<> basis_degMinus1(basis);
    gsMultiBasis<> basis_degMinus1Comp0(basis);
    gsMultiBasis<> basis_degMinus1Comp1(basis);

    basis_degMinus1.degreeReduce(1);
    basis_degMinus1Comp0.degreeReduce(1);
    basis_degMinus1Comp1.degreeReduce(1);
    basis_degMinus1Comp0.degreeElevate(1,1);
    basis_degMinus1Comp1.degreeElevate(1,0);

    gsMultiBasis<> basisFine_degMinus1(basisFine);
    gsMultiBasis<> basisFine_degMinus1Comp0(basisFine);
    gsMultiBasis<> basisFine_degMinus1Comp1(basisFine);

    basisFine_degMinus1.degreeReduce(1);
    basisFine_degMinus1Comp0.degreeReduce(1);
    basisFine_degMinus1Comp1.degreeReduce(1);
    basisFine_degMinus1Comp0.degreeElevate(1,1);
    basisFine_degMinus1Comp1.degreeElevate(1,0);


    if(uniformRefinement)
    {
        // Scordelis-Lo roof (4 patches) to get same h as for 1 patch
        if(testcase == 3)
        {
            numHRefine -=1;
            numHRefineFine -=1;
        }

        // h-refine
        for (int i = 0; i < numHRefine; ++i)
        {
            basis.uniformRefine();
            basis_degMinus1.uniformRefine();
            basis_degMinus1Comp0.uniformRefine();
            basis_degMinus1Comp1.uniformRefine();
        }


        for (int i = 0; i < numHRefineFine; ++i)
        {
            basisFine.uniformRefine();
            basisFine_degMinus1.uniformRefine();
            basisFine_degMinus1Comp0.uniformRefine();
            basisFine_degMinus1Comp1.uniformRefine();
            mpFine.uniformRefine();
        }
    }
    else
    {

        basis.uniformRefine(numInsertKnots);
        basis_degMinus1.uniformRefine(numInsertKnots);
        basis_degMinus1Comp0.uniformRefine(numInsertKnots);
        basis_degMinus1Comp1.uniformRefine(numInsertKnots);



        // only component 1 refined
        /*
        basis.uniformRefineComponent(1,9);
        basis_degMinus1.uniformRefineComponent(1,9);
        basis_degMinus1Comp0.uniformRefineComponent(1,9);
        basis_degMinus1Comp1.uniformRefineComponent(1,9);
*/

        /*
        basis.uniformRefineComponent(1,20);
        basis_degMinus1.uniformRefineComponent(1,20);
        basis_degMinus1Comp0.uniformRefineComponent(1,20);
        basis_degMinus1Comp1.uniformRefineComponent(1,20);



        basis.uniformRefineComponent(0,10);
        basis_degMinus1.uniformRefineComponent(0,10);
        basis_degMinus1Comp0.uniformRefineComponent(0,10);
        basis_degMinus1Comp1.uniformRefineComponent(0,10);
*/
    }

    /*
    const gsTensorBSplineBasis<2>* basisBSpline_temp_const = dynamic_cast< const gsTensorBSplineBasis<2>* >(&basis.piece(0));
    gsTensorBSplineBasis<2>* basisBSpline_temp = const_cast< gsTensorBSplineBasis<2>* >(basisBSpline_temp_const);

    gsKnotVector<> & knots0 = basisBSpline_temp->knots(0);
    knots0[basis.degree()+1] += 0.05;
    gsKnotVector<> & knots1 = basisBSpline_temp->knots(1);
    knots1[basis.degree()+1] += 0.06;
*/


    for(size_t i = 0; i < basis.nBases(); ++i)
        gsInfo<<"Approximation basis (patch "<<i<<"): " << basis.piece(i) << "\n";
    for(size_t i = 0; i < basis.nBases(); ++i)
        gsInfo<<"Approximation basis_degMinus1 (patch "<<i<<"): " << basis_degMinus1.piece(i) << "\n";
    for(size_t i = 0; i < basis.nBases(); ++i)
        gsInfo<<"Approximation basis_degMinus1Comp0 (patch "<<i<<"): " << basis_degMinus1Comp0.piece(i) << "\n";
    for(size_t i = 0; i < basis.nBases(); ++i)
        gsInfo<<"Approximation basis_degMinus1Comp1 (patch "<<i<<"): " << basis_degMinus1Comp1.piece(i) << "\n";


    // construct bases
    std::vector< gsMultiBasis<> >  bases;
    bases.push_back(basis);
    bases.push_back(basis_degMinus1);
    bases.push_back(basis_degMinus1Comp0);
    bases.push_back(basis_degMinus1Comp1);

    std::vector< gsMultiBasis<> >  basesFine;
    basesFine.push_back(basisFine);
    basesFine.push_back(basisFine_degMinus1);
    basesFine.push_back(basisFine_degMinus1Comp0);
    basesFine.push_back(basisFine_degMinus1Comp1);


    /// Assemble and Solve ///

    gsShellMixedAssembler<real_t> assembler(pde, bases, options);
    gsShellMixedAssembler<real_t> assemblerFine(pde, basesFine, options);

    // Assemble
    watch.restart();
    assembler.assemble();
    cout<<"Assembling time: "<<watch.stop()<<endl;

    if(computeFineSolution == true)
        assemblerFine.assemble();

#if (defined(GISMO_WITH_PARDISO))
    gsSparseSolver<>::PardisoLU solverS;
    //gsSparseSolver<>::PardisoLDLT solverS;
#else
    gsSparseSolver<>::LU solverS;
#endif

    watch.restart();
    solverS.compute(assembler.matrix());
    solVector = solverS.solve(assembler.rhs());
    cout<<"Solving time: "<<watch.stop()<<endl;

    if(computeFineSolution == true)
    {
        solverS.compute(assemblerFine.matrix());
        solVectorFine = solverS.solve(assemblerFine.rhs());
    }


    // Compute eigenvalues

    if(outputEvals)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs(assembler.matrix().toDense());
        cout<<"eigs: \n"<<eigs.eigenvalues()<<"\n";
    }

    if(outputMatrix)
    {
        cout<<assembler.matrix().toDense()<<endl;
    }

    if(outputBlock)
    {
        gsSparseSystem<> spSystem = assembler.system();
        gsVector<index_t> blocks_L(12);
        blocks_L<< 0,0,0,1,1,1,2,2,2,3,3,4;
        gsSparseMatrix<>::BlockView matrixBlockView_L = spSystem.blockView(5, 5, blocks_L, blocks_L);
        gsSparseMatrix<> L = matrixBlockView_L(3,0);

        //gsSparseSystem<> spSystem = assembler.system();
        gsVector<index_t> blocks_LnLt(12);
        blocks_LnLt<< 0,0,0,1,1,1,2,2,2,3,4,5;
        gsSparseMatrix<>::BlockView matrixBlockView_LnLt = spSystem.blockView(6, 6, blocks_LnLt, blocks_LnLt);
        gsSparseMatrix<> Ln = matrixBlockView_LnLt(3,0);
        gsSparseMatrix<> Lt = matrixBlockView_LnLt(4,0);

        /*
        cout<<"Ln \n"<<Ln.toDense()<<endl;
        cout<<"Lt \n"<<Lt.toDense()<<endl;
        */

        gsVector<index_t> blocks_L0L1(12);
        blocks_L0L1<< 0,1,1,2,2,2,3,3,3,4,4,5;
        gsSparseMatrix<>::BlockView matrixBlockView_L0L1 = spSystem.blockView(6, 6, blocks_L0L1, blocks_L0L1);
        gsSparseMatrix<> L0 = matrixBlockView_L0L1(4,0);
        gsSparseMatrix<> L1 = matrixBlockView_L0L1(4,1);

        //gsSparseSystem<> spSystem = assembler.system();
        gsVector<index_t> blocks(12);
        blocks<< 0,1,1,2,2,2,3,3,3,4,5,6;
        gsSparseMatrix<>::BlockView matrixBlockView = spSystem.blockView(7, 7, blocks, blocks);
        gsSparseMatrix<> Ln0 = matrixBlockView(4,0);
        gsSparseMatrix<> Ln1 = matrixBlockView(4,1);
        gsSparseMatrix<> Lt0 = matrixBlockView(5,0);
        gsSparseMatrix<> Lt1 = matrixBlockView(5,1);


        //cout<<"Ln0 \n"<<Ln0.toDense()<<endl;
        //cout<<"Ln1 \n"<<Ln1.toDense()<<endl;
        //cout<<"Lt0 \n"<<Lt0.toDense()<<endl;
        //cout<<"Lt1 \n"<<Lt1.toDense()<<endl;


        Eigen::JacobiSVD<Eigen::MatrixXd> svd_L(L.toDense());
        cout<<"L:"<< " rows "<<L.rows()<<" cols "<<L.cols()<<" rank "<<svd_L.rank()<<"\n";

        if(L0.cols()!=0)
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd_L0(L0.toDense());
            cout<<"L0:"<< " rows "<<L0.rows()<<" cols "<<L0.cols()<<" rank "<<svd_L0.rank()<<"\n";
        }
        else
        {
            cout<<"L0:"<< " rows "<<L0.rows()<<" cols "<<L0.cols()<<" rank "<<0<<"\n";
        }
        Eigen::JacobiSVD<Eigen::MatrixXd> svd_L1(L1.toDense());
        cout<<"L1:"<< " rows "<<L1.rows()<<" cols "<<L1.cols()<<" rank "<<svd_L1.rank()<<"\n";

        Eigen::JacobiSVD<Eigen::MatrixXd> svd_Ln(Ln.toDense());
        cout<<"Ln:"<< " rows "<<Ln.rows()<<" cols "<<Ln.cols()<<" rank "<<svd_Ln.rank()<<"\n";
        if(Lt.rows()!=0)
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd_Lt(Lt.toDense());
            cout<<"Lt:"<< " rows "<<Lt.rows()<<" cols "<<Lt.cols()<<" rank "<<svd_Lt.rank()<<"\n";
        }
        else
        {
            cout<<"Lt:"<< " rows "<<Lt.rows()<<" cols "<<Lt.cols()<<" rank "<<0<<"\n";
        }

        if(L0.cols()!=0)
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd_Ln0(Ln0.toDense());
            cout<<"Ln0:"<< " rows "<<Ln0.rows()<<" cols "<<Ln0.cols()<<" rank "<<svd_Ln0.rank()<<"\n";
        }
        else
        {
            cout<<"Ln0:"<< " rows "<<Ln0.rows()<<" cols "<<Ln0.cols()<<" rank "<<0<<"\n";
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd_Ln1(Ln1.toDense());
        cout<<"Ln1:"<< " rows "<<Ln1.rows()<<" cols "<<Ln1.cols()<<" rank "<<svd_Ln1.rank()<<"\n";

        if(Lt.rows()!=0)
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd_Lt0(Lt0.toDense());
            cout<<"Lt0:"<< " rows "<<Lt0.rows()<<" cols "<<Lt0.cols()<<" rank "<<svd_Lt0.rank()<<"\n";
            Eigen::JacobiSVD<Eigen::MatrixXd> svd_Lt1(Lt1.toDense());
            cout<<"Lt1:"<< " rows "<<Lt1.rows()<<" cols "<<Lt1.cols()<<" rank "<<svd_Lt1.rank()<<"\n";
        }
        else
        {
            cout<<"Lt0:"<< " rows "<<Lt0.rows()<<" cols "<<Lt0.cols()<<" rank "<<0<<"\n";
            cout<<"Lt1:"<< " rows "<<Lt1.rows()<<" cols "<<Lt1.cols()<<" rank "<<0<<"\n";
        }

    }


    /*
    gsSparseSystem<> spSystem = assembler.system();
    gsVector<index_t> blocks(12);
    blocks<< 0,1,1,2,2,2,3,3,3,4,5,6;
    gsSparseMatrix<>::BlockView matrixBlockView = spSystem.blockView(7, 7, blocks, blocks);
    gsMatrix<>::BlockView rhsBlockView = spSystem.blockViewRhs(7, blocks);

    gsSparseMatrix<> A00 = matrixBlockView(0,0);
    gsSparseMatrix<> A10 = matrixBlockView(1,0);
    gsSparseMatrix<> A11 = matrixBlockView(1,1);

    gsSparseMatrix<> B0 = matrixBlockView(3,0);
    gsSparseMatrix<> B1 = matrixBlockView(3,1);

    gsSparseMatrix<> C00 = matrixBlockView(2,2);
    gsSparseMatrix<> C10 = matrixBlockView(3,2);
    gsSparseMatrix<> C11 = matrixBlockView(3,3);

    gsSparseMatrix<> LtMean = matrixBlockView(6,5);
    gsSparseMatrix<> LnMean = matrixBlockView(6,4);


    //cout<<"A00 \n"<<A00<<endl;
    //cout<<"A10 \n"<<A10<<endl;
    //cout<<"A11 \n"<<A11<<endl;

    //cout<<"B0 \n"<<B0<<endl;
    //cout<<"B1 \n"<<B1<<endl;

    cout<<"LtMean \n"<<LtMean<<endl;
    cout<<"LnMean \n"<<LnMean<<endl;

*/


    /// Construct solution ///

    gsVector<index_t> unk_p(1);
    unk_p << 0;
    gsVector<index_t> unk_phi(2);
    unk_phi << 1,2;
    gsVector<index_t> unk_u(3);
    unk_u << 6,7,8;
    //gsVector<index_t> unk_M(3);
    //unk_M << 0,1,2;

    if(options.getSwitch("Mmixed"))
    {
        assembler.constructSolution(solVector, sol_p, unk_p);
        assembler.constructSolution(solVector, sol_phi, unk_phi);
        //assembler.constructSolution(solVector, sol_M, unk_M);
    }

    /*
    else
    {
        gsMultiPatch<> zero_p;
        gsMultiPatch

        for(index_t i = 0; i < mp.nPatches(); ++i)
        {
            gsMatrix<T> coefs_p;
            coefs_p.setZero(basis.size(), 1);
            gsMatrix<T> coefs_phi;
            coefs_phi.setZero(basis.size(), 2);

            gsGeometry<> zero_p(basis, coefs_p);
            gsGeometry<> zero_phi(basis, coefs_phi);
        }


        //gsMultiPatch<> empty;
        sol_p = mp;
        sol_phi = mp;
    }*/

    assembler.constructSolution(solVector, sol_u, unk_u);
    if(computeFineSolution)
        assemblerFine.constructSolution(solVectorFine, solFine_u, unk_u);
    else
        solFine_u = sol_u;

    // Compute physical displacement
    gsFunctionSetPDisplace<real_t>  sol_uPhys(mp, sol_u);
    gsFunctionSetPDisplace<real_t>  solFine_uPhys(mp, solFine_u);
    gsFunctionSetPDisplaceModulus<real_t>  sol_uPhysModulus(mp, sol_u);
    //gsFunctionSetPDisplace<real_t>::Ptr  sol_uPhys = gsFunctionSetPDisplace<real_t>::make(mp, sol_u);
    //gsFunctionSetM<real_t>  sol_M(mp, sol_p, sol_phi);

    // Approximate physical displacement by Bspline multipatch for plotting
    gsMultiPatch<> solApprox_uPhys (mp);
    for(size_t i = 0; i < mp.nPatches(); ++i)
    {
        const gsTensorBSplineBasis<2>& basisBSpline_temp = dynamic_cast< const gsTensorBSplineBasis<2>& >(mp.basis(i).source());
        gsGeometry<>::uPtr BSpline_temp = basisBSpline_temp.interpolateAtAnchors(sol_uPhys.piece(i).eval(basisBSpline_temp.anchors()));
        solApprox_uPhys.patch(i) = *BSpline_temp;
    }


    // Construct solution fields
    // uPhys
    gsField<> solField_u(mp, sol_u);
    const gsField<> solField_uPhysModulus(mp, sol_uPhysModulus, true);
    const gsField<> solField_uPhys(mp, sol_uPhys, true);
    const gsField<> solFieldFine_uPhys(mp, solFine_uPhys, true);
    //const gsField<> solField_M(mp, sol_M, true);

    // uPhys on derformed mp
    gsMultiPatch<> mpDeformed(mp);
    for(size_t i = 0; i < mp.nPatches(); ++i)
        mpDeformed.patch(i).coefs() += mpDeformedScaling*solApprox_uPhys.patch(i).coefs();
    const gsField<> solField_uPhysModulus_mpDeformed(mpDeformed, sol_uPhysModulus, true);
    //const gsField<> solField_M_mpDeformed(mpDeformed, sol_M, true);


    //gsField<> solField_p(mp, sol_p);
    //gsField<> solField_phi(mp, sol_phi);
    //gsField<> solField_M(mp, sol_M);

    /// Compute strain energy ///

    //gsMatrix<> strainEnergy = solVector.transpose() * (assembler.matrix() *  solVector); // wrong way to compute energy for mixed formulation
    //cout<<"strainEnergy "<<strainEnergy<<endl;


    /// Compute discretization errors ///

    real_t errorH1 = 0;
    real_t errorH1_uPhys = 0;

    if(writeFileConvergenceNorm)
    {

        if(analytSol)
        {
            gsNormL2<real_t> L2Norm(solField_uPhys, analytSol_u);
            gsSeminormH1<real_t> H1Seminorm(solField_uPhys, analytSol_u);

            real_t errorL2 = L2Norm.compute();
            real_t errorH1Seminorm = H1Seminorm.compute();

            errorH1 = sqrt(errorL2 * errorL2 + errorH1Seminorm * errorH1Seminorm);
        }
        else
        {
            if(computeFineSolution)
            {
                for (size_t pn=0; pn < mp.nPatches(); ++pn )
                {
                    gsMultiPatch<> mpPatch(mpFine.patch(pn));
                    //gsMultiPatch<> fineSolPatch_u(solFine_uPhys.piece(pn));
                    gsField<> fineFieldPatch_u(mpPatch, solFine_uPhys.piece(pn));

                    gsNormL2<real_t> L2Norm(fineFieldPatch_u, sol_uPhys.piece(pn), true);
                    gsSeminormH1<real_t> H1Seminorm(fineFieldPatch_u, sol_uPhys.piece(pn), true);

                    real_t errorL2 = L2Norm.compute();
                    real_t errorH1Seminorm = H1Seminorm.compute();

                    errorH1_uPhys += errorL2*errorL2 + errorH1Seminorm*errorH1Seminorm;
                }

                for (size_t pn=0; pn < mp.nPatches(); ++pn )
                {
                    gsMultiPatch<> mpPatch(mp.patch(pn));
                    gsMultiPatch<> solFinePatch_u(solFine_u.piece(pn));
                    gsField<> fineFieldPatch_u(mpPatch, solFinePatch_u);

                    gsNormL2<real_t> L2Norm(fineFieldPatch_u, sol_u.piece(pn), true);
                    gsSeminormH1<real_t> H1Seminorm(fineFieldPatch_u, sol_u.piece(pn), true);

                    real_t errorL2 = L2Norm.compute();
                    real_t errorH1Seminorm = H1Seminorm.compute();

                    errorH1 += errorL2*errorL2 + errorH1Seminorm*errorH1Seminorm;
                }


                errorH1_uPhys = sqrt(errorH1_uPhys);
                errorH1 = sqrt(errorH1);


            }
            else
                cout<<"Error: no analytical or fine solution to compare with"<<endl;
        }


        std::ofstream output("dataConvergenceNorm", std::fstream::app|std::fstream::out);
        if(!output.is_open())
            GISMO_ERROR("Opening file failed");

        output<<testcase<<","<<basis.degree()<<","<<basis.piece(0).component(0).size()<<","<<basis.piece(0).getMaxCellLength()<<","<<errorH1_uPhys<<","<<errorH1<<"\n";


    }


    /// Compute point values of solution and write to file ///

    if(writeFileConvergence)
    {

        std::ofstream output("dataConvergence", std::fstream::app|std::fstream::out);
        if(!output.is_open())
            GISMO_ERROR("Opening file failed");

        gsMatrix<> y(2,1);
        gsMatrix<> yN(2,1);
        gsMatrix<> yS(2,1);

        if(testcase == 2)
        {
            yN<< 0.5, 1;
            yS<< 0.5, 0;
            output<<testcase<<","<<basis.degree()<<","<<basis.piece(0).component(0).size()<<","<<solField_uPhysModulus.value(yN,0)(2)<<","<<solField_uPhysModulus.value(yS,0)(2)<<"\n";
        }
        else if(testcase == 3)
        {
            y<< 0.5, 1;
            output<<testcase<<","<<basis.degree()<<","<<basis.piece(0).component(0).size()<<","<<solField_uPhysModulus.value(y,1)(2)<<"\n";
        }
        else if(testcase == 4)
        {
            y<< 0.5, 1;
            output<<testcase<<","<<basis.degree()<<","<<basis.piece(0).component(0).size()<<","<<solField_uPhysModulus.value(y,2)(2)<<"\n";
        }
        else if(testcase == 5)
        {
            y<< 0, 0;
            output<<testcase<<","<<basis.degree()<<","<<basis.piece(0).component(0).size()<<","<<solField_uPhysModulus.value(y,1)(1)<<"\n";
        }
        else if(testcase == 6)
        {
            y<< 0.5, 0;
            output<<testcase<<","<<basis.degree()<<","<<basis.piece(0).component(1).size()<<","<<solField_uPhysModulus.value(y,0)(1)<<"\n";
        }
        else if(testcase == 7)
        {
            y<< 1, 0.5;
            output<<testcase<<","<<basis.degree()<<","<<basis.piece(0).component(1).size()<<","<<solField_uPhysModulus.value(y,0)(2)<<"\n";
        }
        else if(testcase == 8)
        {
            y<< 0.5, 1;
            output<<testcase<<","<<basis.degree()<<","<<basis.piece(0).component(0).size()<<","<<solField_uPhysModulus.value(y,0)(2)<<"\n";
        }
        output.close();
    }

    if(writeFileLocking)
    {
        std::ofstream output("dataLocking", std::fstream::app|std::fstream::out);
        if(!output.is_open())
            GISMO_ERROR("Opening file failed");

        gsMatrix<> y(2,1);

        if(testcase == 6)
        {
            y<< 0.5, 0;
            output<<testcase<<","<<basis.degree()<<","<<3*basis.totalSize()<<","<<10/thickness<<","<<solField_uPhysModulus.value(y,0)(1)<<"\n";
        }
        output.close();
    }


    // Scordelis-Lo roof
    if(testcase == 2)
    {
        cout<<"CP per side"<<endl;
        cout<<basis.piece(0).component(0).size()<<endl;

        gsMatrix<> y(2,1);
        y<< 0.5, 0;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,0).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,0)<<endl;
        cout<<"solField_M"<<endl;
        //cout<<solField_M.value(y,0)<<endl;

        y<< 0.5, 1;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,0).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,0)<<endl;

        /*
        y<< 0.5, 0.5;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,0).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,0)<<endl;
        */

    }
    // Scordelis-Lo roof (4 patches)
    else if (testcase == 3)
    {
        gsMatrix<> y(2,1);
        y<< 0.5, 1;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,1).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,1)<<endl;

        y<<0.5, 0.5;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,3).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,3)<<endl;

    }
    // Pinched cylinder (4 patches)
    else if (testcase == 4)
    {
        gsMatrix<> y(2,1);
        /*
        gsMatrix<> y(2,1);
        y<< 0.5, 1;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,0).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,0)<<endl;

        y<< 0.5, 1;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,2).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,2)<<endl;
        */

        y<< 0,1;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,0).transpose()<<endl;
    }
    // Pinched hemisphere (4 patches)
    else if (testcase == 5)
    {
        cout<<"CP per side"<<endl;
        cout<<basis.piece(0).component(0).size()<<endl;

        gsMatrix<> y(2,1);
        y<< 0, 0;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,0).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,0)<<endl;

        y<< 0, 0;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,1).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,1)<<endl;


        y<< 0, 0;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,2).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,2)<<endl;

        y<< 0, 0;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,3).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,3)<<endl;

        y<< 0, 1;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,3).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,3)<<endl;

    }
    // Cylinder shell strip (1 patches)
    else if (testcase == 6)
    {
        gsMatrix<> y(2,1);
        y<< 0.5, 0;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,0).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,0)<<endl;
    }

    else if (testcase == 7)
    {
        gsMatrix<> y(2,1);
        y<< 1, 0.5;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,0).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,0)<<endl;
    }

    else if (testcase == 8)
    {
        gsMatrix<> y(2,1);
        y<< 0.5, 1;
        cout<<"point_physical"<<endl;
        cout<<solField_uPhysModulus.point(y,0).transpose()<<endl;
        cout<<"solField_u"<<endl;
        cout<<solField_uPhysModulus.value(y,0)<<endl;
    }



    /// Plot solution ///

    if (plot)
    {
        // Write solution to paraview file
        std::cout<<"Plotting in Paraview...\n";
        gsWriteParaview<>( mp, "geometry", 1000, true, true);
        //gsWriteParaview<>( basis.front(), "basis", 1000);


        gsWriteParaview<>( solField_uPhysModulus, "sol_u_mpUndeformed", 5000, false);
        gsWriteParaview<>( solField_uPhysModulus_mpDeformed, "sol_u", 5000, false);
        //gsWriteParaview<>( solField_M_mpDeformed, "sol_M", 10000, false);
        //gsWriteParaview<>( solField_p, "sol_p", 1000);
        //gsWriteParaview<>( solField_phi, "sol_phi", 1000);

        // Run paraview on exit
        result = system("paraview sol_u_mpUndeformed.pvd &");
        //result = system("paraview geometry.pvd &");

    }

    return result;
}




