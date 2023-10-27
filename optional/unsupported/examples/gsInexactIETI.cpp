/**  gsInexactIETI.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):     C. Hofer
    Created on:  2016-07-01

    DESCRIBE THE PURPOSE OF THE EXAMPLE HERE

    WHAT IS HAPPENING

    EXPECTED RESULT
*/




#include <iostream>
#include <gismo.h>
#include <gismo_dev.h>
#include <gsIETI/gsIETIUtils.h>
#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETIdGAssembler.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIScaledDirichlet.h>


#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>

#include <gsSolver/gsConjugateGradient.h>

using namespace gismo;

void makeYETIFootprintNice(gsMultiPatch<real_t>::uPtr& patches)
{
    for(size_t np = 0; np< patches->nPatches(); ++np)
    {
        /*
        if(np==0)
        {
            std::vector<gsGeometry<>* >  ptch_ =  (patches)->patch(np).uniformSplit(0);
            ptch.insert(ptch.end(),ptch_.begin(),ptch_.end());
        }
        else
        {
            ptch.insert(ptch.end(), patches->patch(np).clone());
        }
        */
        if(np==0 || np == 10 || np ==20 )
        {
            gsTensorBSpline<2,real_t>* TB = dynamic_cast<gsTensorBSpline<2,real_t>* >(&patches->patch(np));
            TB->insertKnot(0.25,1);
            TB->insertKnot(0.75,1);
        }
        else if( np == 5 || np == 6  || np == 7)
        {
            gsTensorBSpline<2,real_t>* TB = dynamic_cast<gsTensorBSpline<2,real_t>* >(&patches->patch(np));
            TB->insertKnot(0.25,0);
            TB->insertKnot(0.75,0);
        }
        if(np==19)
        {
            gsTensorBSpline<2,real_t>* TB = dynamic_cast<gsTensorBSpline<2,real_t>* >(&patches->patch(np));
            TB->insertKnot(0.875,1);
        }
        if(np == 18  || np == 16)
        {
            gsTensorBSpline<2,real_t>* TB = dynamic_cast<gsTensorBSpline<2,real_t>* >(&patches->patch(np));
            TB->insertKnot(0.125,1);
            TB->insertKnot(0.375,1);
            TB->insertKnot(0.625,1);
            TB->insertKnot(0.875,1);
        }
    }

}

/// A locally used geometry
gsMultiPatch<real_t> approximateQuarterAnnulus(int deg)
{
    gsGeometry<>::uPtr quann = gsNurbsCreator<>::NurbsQuarterAnnulus();

    gsKnotVector<> KV1(0,1, 0, deg+1);        // no interior knots in x direction
    gsKnotVector<> KV2(0,1, 0, deg+1);        // no interior knot in y direction

    gsTensorBSplineBasis<2> tbsp (give(KV1), give(KV2));
    gsMatrix<real_t> eval = quann->eval(tbsp.anchors());
    gsGeometry<>::uPtr approxGeom = tbsp.interpolateAtAnchors( eval );
    gsMultiPatch<real_t> mp(*approxGeom);

    //gsMultiPatch<real_t> res = mp.uniformSplit();
    return mp;

}

/// A locally used geometry
gsMultiPatch<real_t> approximateTwistedQuarterAnnulus(int deg)
{
    gsTensorBSpline<3>::Ptr mp2;
    gsFileData<> fileData("volumes/twistedFlatQuarterAnnulus.xml");
    if (fileData.has< gsTensorBSpline<3> >())
        mp2 = fileData.getFirst< gsTensorBSpline<3> >();
    else
    {gsInfo<<"Error!!!\n";}


    gsKnotVector<> KV1(0,1, 0, deg+1);        // no interior knots in x direction
    gsKnotVector<> KV2(0,1, 0, deg+1);        // no interior knot in y direction
    gsKnotVector<> KV3(0,1, 0, deg+1);

    gsTensorBSplineBasis<3> tbsp (give(KV1), give(KV2), give(KV3));
    gsGeometry<>::uPtr approxGeom = tbsp.interpolateAtAnchors( mp2->eval(tbsp.anchors()) );
    gsMultiPatch<real_t> mp(*approxGeom);
    return mp;
}


//Important Note: The boundary condition has to be called with
// bcInfo.addCondition(***, *** , *** , g);
// in order that the parallelization works correctly! So no &g !!!
// Otherwise, you have a data race!



real_t getDefaultDamping(std:: string smooth )
{
    if (smooth == "r" || smooth == "Richardson")
        return 0.80;
    else if (smooth == "j" || smooth == "Jacobi")
        return 0.80;
    else if (smooth == "gs" || smooth == "GaussSeidel")
        return 1;
    else if (smooth == "mrb" || smooth == "mass-richardson-boundary-correction")
        return 0.09;
    else if (smooth == "mrs" || smooth == "SubspaceCorrectedMassSmoother")
        return 0.09;
    else if (smooth == "mrst" || smooth == "mass-richardson-subspace-correction-trunc")
        return 0.09;
    else if (smooth == "mrs-gs" || smooth == "mass-richardson-subspace-correction-gauss-seidel")
        return 0.09;
    else
        return 1;
}



int main (int argc, char** args)
{
    index_t testcase =2;
    bool test3D = false;
    bool plot = false;
    bool nitsche = false;
    bool timings = true;
    bool saddlePoint = false;
    bool CalcEigenvals = true;
    //Smoother::type smoother;
    index_t jump = 0;
    std::string out = "";

    //primalDofMethod::strategy pstrat = primalDofMethod::A;
    bool NoMinEnergy = false;
    bool needrescal = false;
    bool nonlinear = false;
    bool calcRhs = false;
    std::string m = "A";
    std::string solvI,solvR,solvKC;
    std::string scaling = "coeff";
    solvI = solvR = solvKC = "D";
    real_t damping = -1;
    real_t coarseDamping = 1;
    index_t levelsKC =-1;
    index_t cyclesKC = 1;
    index_t presmoothKC = 1;
    index_t postsmoothKC = 1;
    real_t dampingKC = -1;
    real_t outerDampingKC = 1;
    std::string smoothKC = "GaussSeidel";
    bool useP1ProjectionKC = false;
    std::string smoothKC_TL = "SubspaceCorrectedMassSmoother";
    index_t presmoothKC_TL = 1;
    index_t postsmoothKC_TL = 1;
    real_t coarseDampingKC = -1;
    real_t dampingKC_TL = -1;
    real_t outerDampingKC_TL = -1;

    index_t levelsKii =-1;
    index_t cyclesKii = 1;
    index_t presmoothKii = 1;
    index_t postsmoothKii = 1;
    real_t dampingKii = -1;
    real_t outerDampingKii = 1;
    std::string smoothKii = "GaussSeidel";
    bool useP1ProjectionKii = false;
    std::string smoothKii_TL = "SubspaceCorrectedMassSmoother";
    index_t presmoothKii_TL = 1;
    index_t postsmoothKii_TL = 1;
    real_t coarseDampingKii = -1;
    real_t dampingKii_TL = -1;
    real_t outerDampingKii_TL = -1;

    real_t reg = 3.e-2;

    index_t maxIterSolveKC = 200;
    real_t tolSolveKC = 1.e-10;
    index_t maxIterSolveRHSKii = 200;
    real_t tolSolveRHSKii = 1.e-10;
    index_t maxIterSolvePrecKii = 2;
    real_t tolSolvePrecKii = 1.e-6;
    real_t tolBasisKC = 1.e-12;
    index_t maxIterBasisKC = 200;

    bool basedOnRankOneCorrectionKii, basedOnRankOneCorrectionKC, basedOnRankOneCorrectionKii_TL, basedOnRankOneCorrectionKC_TL;
    basedOnRankOneCorrectionKii = basedOnRankOneCorrectionKC = basedOnRankOneCorrectionKii_TL = basedOnRankOneCorrectionKC_TL=false;

    real_t BPCG_SZ_scal1 = 1.24;
    real_t BPCG_SZ_scal2 = 0.99/BPCG_SZ_scal1;

    real_t tolerance = 1.e-6;
    index_t maxIterations = 200;

    real_t BPCG_combScal = 1;
    real_t BPCG_Scal = 1;
    real_t BPCG_Scal2 = 1;



    bool dG=false;
    // Number for h-refinement of the computational (trail/test) basis.
    index_t numRefine  = 2;
    // Number for p-refinement of the computational (trail/test) basis.
    index_t numElevate = 0;
    index_t n=4;

    gsCmdLine cmd("Hi, I will test the inexact IETI versions");
    cmd.addSwitch("d", "dim", "test3D", test3D);
    cmd.addSwitch("",  "nitsche","elimination or nitsche",nitsche);
    cmd.addSwitch(     "plot", "Plot result in ParaView format", plot);


    cmd.addSwitch("",  "dG","use dG formulation",dG);
    cmd.addInt("r","refine","Number of refinements",numRefine);
    cmd.addInt("c","case","Choosen testcase",testcase);
    cmd.addInt("e","elevate","maximal elevation",numElevate);
    cmd.addInt("n","number","number of patches in one direction",n);
    cmd.addString("o", "out", "output Filename", out);
    cmd.addInt("j","jumps","jumping coeff from 10^-j to 10^j",jump);

    cmd.addSwitch("",  "IETI.ExtraTiming", "enable extraTimings", timings);
    cmd.addSwitch("",  "IETI.NoMinimumEnergy","choose energy minimizing prim subspaces",NoMinEnergy);
    cmd.addSwitch("",  "IETI.EnforceRescaling","enforces rescaling of local matrices", needrescal);
    cmd.addSwitch("",  "IETI.NonlinearMode","enables the usage for nonlinear PDEs", nonlinear);
    cmd.addSwitch("",  "IETI.CalcRhsNorm","calculates the norm of the original rhs", calcRhs);


    cmd.addString("I", "IETI.SolverKii", "chosen solver for Kii",solvI);
    cmd.addString("C", "IETI.SolverKC", "chosen solver for KC",solvKC);
    cmd.addString("R", "IETI.SolverKrr", "chosen solver for Krr",solvR);
    cmd.addString("",  "IETI.Scaling", "chosen scaling",scaling);
    cmd.addString("s", "IETI.Strategy", "Choosen strategy",m);
    cmd.addSwitch("",  "IETI.SaddlePoint","use saddle point for inexact", saddlePoint);

    cmd.addInt   ("",  "IETI.MG_KC.Levels","Number of levels to use for multigrid iteration", levelsKC);
    cmd.addInt   ("",  "IETI.MG_KC.Cycles",             "Number of multi-grid cycles", cyclesKC);
    cmd.addInt   ("",  "IETI.MG_KC.Presmooth",          "Number of pre-smoothing steps", presmoothKC);
    cmd.addInt   ("",  "IETI.MG_KC.Postsmooth",         "Number of post-smoothing steps", postsmoothKC);
    cmd.addReal  ("",  "IETI.MG_KC.Damping",            "Damping factor for the smoother (handed over to smoother)", dampingKC);
    cmd.addReal  ("",  "IETI.MG_KC.OuterDamping",       "Damping factor for the smoother (globally)", outerDampingKC);
    cmd.addString("",  "IETI.MG_KC.Smoother",            "Smoothing method", smoothKC);
    cmd.addSwitch("",  "IETI.MG_KC.UseP1Projection","enable the P1 projection of coarser grids, then uses TwoLevel", useP1ProjectionKC);
    cmd.addSwitch("",  "IETI.MG_KC.BasedOnRankOneCorrection", "Setup the mass smoother (if chosen) based on a rank one approximation of the matrices", basedOnRankOneCorrectionKC );
    cmd.addInt   ("",  "IETI.MG_KC.TL.Presmooth",          "Number of pre-smoothing steps", presmoothKC_TL);
    cmd.addInt   ("",  "IETI.MG_KC.TL.Postsmooth",         "Number of post-smoothing steps", postsmoothKC_TL);
    cmd.addReal  ("",  "IETI.MG_KC.TL.CoarseDamping",      "Damping factor for the two level smoother (handed over to smoother)", coarseDampingKC);
    cmd.addReal  ("",  "IETI.MG_KC.TL.Damping",        "Damping factor for the smoother (handed over to smoother)", dampingKC_TL);
    cmd.addReal  ("",  "IETI.MG_KC.TL.OuterDamping",       "Damping factor for the smoother (globally)", outerDampingKC_TL);
    cmd.addString("",  "IETI.MG_KC.TL.Smoother",       "Smoothing method", smoothKC_TL);
    cmd.addSwitch("",  "IETI.MG_KC.TL.BasedOnRankOneCorrection", "Setup the mass smoother (if chosen) based on a rank one approximation of the matrices", basedOnRankOneCorrectionKC_TL );

    cmd.addInt   ("",  "IETI.MG_Kii.Levels","Number of levels to use for multigrid iteration", levelsKii);
    cmd.addInt   ("",  "IETI.MG_Kii.Cycles",             "Number of multi-grid cycles", cyclesKii);
    cmd.addInt   ("",  "IETI.MG_Kii.Presmooth",          "Number of pre-smoothing steps", presmoothKii);
    cmd.addInt   ("",  "IETI.MG_Kii.Postsmooth",         "Number of post-smoothing steps", postsmoothKii);
    cmd.addReal  ("",  "IETI.MG_Kii.Damping",            "Damping factor for the smoother (handed over to smoother)", dampingKii);
    cmd.addReal  ("",  "IETI.MG_Kii.OuterDamping",       "Damping factor for the smoother (globally)", outerDampingKii);
    cmd.addString("",  "IETI.MG_Kii.Smoother",             "Smoothing method", smoothKii);
    cmd.addSwitch("",  "IETI.MG_Kii.UseP1Projection","enable the P1 projection of coarser grids, then uses TwoLevel", useP1ProjectionKii);
    cmd.addSwitch("",  "IETI.MG_Kii.BasedOnRankOneCorrection", "Setup the mass smoother (if chosen) based on a rank one approximation of the matrices", basedOnRankOneCorrectionKii );
    cmd.addInt   ("",  "IETI.MG_Kii.TL.Presmooth",          "Number of pre-smoothing steps", presmoothKii_TL);
    cmd.addInt   ("",  "IETI.MG_Kii.TL.Postsmooth",         "Number of post-smoothing steps", postsmoothKii_TL);
    cmd.addReal  ("",  "IETI.MG_Kii.TL.CoarseDamping",      "Damping factor for the two level smoother (handed over to smoother)", coarseDampingKii);
    cmd.addReal  ("",  "IETI.MG_Kii.TL.Damping",         "Damping factor for the smoother (handed over to smoother)", dampingKii_TL);
    cmd.addReal  ("",  "IETI.MG_Kii.TL.OuterDamping",       "Damping factor for the smoother (globally)", outerDampingKii_TL);
    cmd.addString("",  "IETI.MG_Kii.TL.Smoother",       "Smoothing method", smoothKii_TL);
    cmd.addSwitch("",  "IETI.MG_Kii.TL.BasedOnRankOneCorrection", "Setup the mass smoother (if chosen) based on a rank one approximation of the matrices", basedOnRankOneCorrectionKii_TL );

    cmd.addReal ("", "IETI.Regularization", " regularization parameter for local neumann problems (inexact problems) ", reg);

    cmd.addReal ("", "IETI.MG_KC.Solve.Tolerance",           "Stopping criterion for cg", tolSolveKC);
    cmd.addInt  ("", "IETI.MG_KC.Solve.MaxIterations",       "Stopping criterion for cg", maxIterSolveKC);
    cmd.addReal ("", "IETI.MG_KC.Basis.Tolerance",           "Stopping criterion for cg", tolBasisKC);
    cmd.addInt  ("", "IETI.MG_KC.Basis.MaxIterations",       "Stopping criterion for cg", maxIterBasisKC);
    cmd.addReal ("", "IETI.MG_KC.Basis.PreconditionScaling1","K-Scaling for SchöberlZuLehner", BPCG_SZ_scal1);
    cmd.addReal ("", "IETI.MG_KC.Basis.PreconditionScaling2","S-Scaling for SchöberlZuLehner", BPCG_SZ_scal2);
    cmd.addReal ("", "IETI.MG_Kii.RHS.Tolerance",            "Stopping criterion for cg", tolSolveRHSKii);
    cmd.addInt  ("", "IETI.MG_Kii.RHS.MaxIterations",        "Stopping criterion for cg", maxIterSolveRHSKii);
    cmd.addReal ("", "IETI.MG_Kii.Prec.Tolerance",           "Stopping criterion for cg", tolSolvePrecKii);
    cmd.addInt  ("", "IETI.MG_Kii.Prec.MaxIterations",       "Stopping criterion for cg", maxIterSolvePrecKii);


    cmd.addReal  ("", "Tolerance",           "Stopping criterion for cg", tolerance);
    cmd.addInt   ("",  "MaxIterations",      "Stopping criterion for cg", maxIterations);
    cmd.addSwitch("", "CalcEigenvalues",    "calculation of Eigenvalues in cg", CalcEigenvals);
    cmd.addReal  ("", "CombinationScaling",           "Stopping criterion for cg", BPCG_combScal);
    cmd.addReal  ("", "PreconditionScaling1",           "Stopping criterion for cg", BPCG_Scal);
    cmd.addReal  ("", "PreconditionScaling2",           "Stopping criterion for cg", BPCG_Scal2);


    //  cmd.addString("R", "solverKrr", "chosen solver for Krr",solvR);
    try { cmd.getValues(argc,args); } catch (int rv) { return rv; }

    gsOptionList option = cmd.getOptionList();

    dirichlet::strategy dstrat = nitsche ? dirichlet::nitsche : dirichlet::elimination;
    iFace::strategy fstrat =  dG ? iFace::dg : iFace::glue;

    //Handle trivial cases
    if (option.getInt("IETI.MG_Kii.Levels") <0)  option.setInt( "IETI.MG_Kii.Levels", numRefine -1 );
    if (option.getInt("IETI.MG_KC.Levels") <0)   option.setInt( "IETI.MG_KC.Levels", numRefine -1 );

    if(option.getReal("IETI.MG_KC.Damping")<0) option.setReal("IETI.MG_KC.Damping",getDefaultDamping(option.getString("IETI.MG_KC.Smoother")));
    if(option.getReal("IETI.MG_KC.TL.Damping")<0) option.setReal("IETI.MG_KC.TL.Damping",getDefaultDamping(option.getString("IETI.MG_KC.TL.Smoother")));
    if(option.getReal("IETI.MG_Kii.Damping")<0) option.setReal("IETI.MG_Kii.Damping",getDefaultDamping(option.getString("IETI.MG_Kii.Smoother")));
    if(option.getReal("IETI.MG_Kii.TL.Damping")<0) option.setReal("IETI.MG_Kii.TL.Damping",getDefaultDamping(option.getString("IETI.MG_Kii.TL.Smoother")));

    gsInfo << "Run gsInexactIETI with options:\n" << option << std::endl;
    
    int dim;
    (test3D) ? dim =3:dim =2;

    if(testcase == 9) test3D=false;

    gsInfo<<"Dimension of Domain: "<<dim<<"\n";
    gsInfo<<"Using Nitsche: "<<nitsche<<"\n";

    /*
    if (smooth == "r" || smooth == "richardson")
        smoother = Smoother::Richardson;
    else if (smooth == "j" || smooth == "jacobi")
        smoother = Smoother::Jacobi;
    else if (smooth == "gs" || smooth == "gauss-seidel")
        smoother = Smoother::GaussSeidel;
    else if (smooth == "mr" || smooth == "mass-richardson")
        smoother = Smoother::MassRichardson;
    else if (smooth == "mrb" || smooth == "mass-richardson-boundary-correction")
        smoother = Smoother::MassRichardsonBoundaryCorrection;
    else if (smooth == "mrs" || smooth == "mass-richardson-subspace-correction")
        smoother = Smoother::MassRichardsonSubspaceCorrection;
    else if (smooth == "mrst" || smooth == "mass-richardson-subspace-correction-trunc")
        smoother = Smoother::MassRichardsonSubspaceCorrectionTrunc;
    else if (smooth == "mrs-gs" || smooth == "mass-richardson-subspace-correction-gauss-seidel")
        smoother = Smoother::MassRichardsonSubspaceCorrectionGS;
    else if (smooth == "g" || smooth == "generic")
        smoother = Smoother::Generic;
    else
    {
        GISMO_ERROR( "Unknown smoother \"" << smooth << "\".\n"
                     << "Allowed are: richardson (r), jacobi (j), and gauss-seidel (gs).\n"
                     << "Moreover, the following experimental smoothers exist: mass-richardson (mr), mass-richardson-boundary-correction (mrb), mass-richardson-subspace-correction (mrs), and generic (g).\n");
    }
*/
    //////// Right-hand side and analytical solution ////////
    // Define source function7
    gsFunctionExpr<> f,g;
    if(testcase == 10)
    {
        f=gsFunctionExpr<>("0",dim);
        g=gsFunctionExpr<>("0",dim);
    }
    else
    {
        // f=gsFunctionExpr<>("((0.2*pi*1)^2 + (pi*0.4)^2)*sin(0.2*pi*(x+0.4))*sin(pi*(y+0.3)*0.4)",dim);
        // g=gsFunctionExpr<>("sin(0.2*pi*(x+0.4))*sin(pi*(y+0.3)*0.4)+x+y",dim);
        f=gsFunctionExpr<>("((pi*8)^2 + (pi*9)^2)*sin(8*x*pi)*sin(pi*y*9)",dim);
        g=gsFunctionExpr<>("sin(8*x*pi)*sin(pi*y*9)",dim);

    }

    gsPiecewiseFunction<real_t>  alpha;
    gsFunctionExpr<> a1("1.e-"+util::to_string(jump),dim);
    gsFunctionExpr<> a2("1.e+"+util::to_string(jump),dim);

    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n" << "\n";

    /////////////////// Setup geometry ///////////////////
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> mp, patches;


    switch(testcase){
    case 1:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,1, 0.5);
        patches.degreeElevate(1);
        break;
    case 2:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,2, 0.5);
        patches.degreeElevate(1);
        break;
    case 3:
        patches = gsNurbsCreator<>::BSplineSquareGrid(3,1, 1);
        patches.degreeElevate(1);
        break;
    case 4:
    case 5:
        patches = gsNurbsCreator<>::BSplineSquareGrid(1,2, 1);
        patches.degreeElevate(1);
        break;
    case 6:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,2, 1);
        patches.degreeElevate(1);
        break;
    case 7:
    case 10:
    {
        gsFileData<> fileData("yeti_mp2.xml");
        //gsFileData<> fileData("gaps/chairStandTopFace.xml");
        if (fileData.has< gsMultiPatch<> >())
            fileData.getFirst(patches);
        else
            return 1;

      //  makeYETIFootprintNice(patches);

        /*
        mp = patches->uniformSplit();
       // mp= mp.uniformSplit();
        mp.uniformRefine();
        for(index_t np = 0; np< mp.nPatches(); ++np)
            (mp).patch(np).coefs()*=4;
        // delete patches;
        *patches = mp; //just change the pointer for compatibility
        //*/
        break;
    }
    case 8:
    {
        mp= approximateQuarterAnnulus(2);
        {
            
            std::vector<gsGeometry<>*> ptch;
            for(size_t np = 0; np< mp.nPatches(); ++np)
            {
                std::vector<gsGeometry<>* >  ptch_ =  (mp).patch(np).uniformSplit(1);

                ptch.insert(ptch.end(),ptch_.begin(),ptch_.end());
            }
            mp = gsMultiPatch<real_t>(ptch);
            
            // delete patches;
            for(int i=0; i<log2(n);++i)
                mp= mp.uniformSplit();
            patches = mp; //just change the pointer for compatibility
            patches.computeTopology();
        }

        break;
    }
    case 9:
    {
        mp = approximateTwistedQuarterAnnulus(2);

        std::vector<gsGeometry<>*> ptch;
        for(size_t np = 0; np< mp.nPatches(); ++np)
        {
            std::vector<gsGeometry<>* >  ptch_ =  (mp).patch(np).uniformSplit(1);

            ptch.insert(ptch.end(),ptch_.begin(),ptch_.end());
        }

        mp= gsMultiPatch<real_t>(ptch);

        for(int i=0; i<log2(n);++i)
            mp= mp.uniformSplit();
        patches = mp;
        patches.computeTopology();
    }
        break;
    default:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,1, 1);
        patches.degreeElevate(1);
        break;
    }

#ifdef _OPENMP
    Eigen::initParallel();
    // To get the number of threads from the environment
    int nThreads;
    const char * nProcs = getenv("OMP_NUM_THREADS");
    if(nProcs != NULL)
        sscanf( nProcs, "%d", &nThreads );
    else
        nThreads = 1; //one OpenMP thread
    nThreads = math::min(nThreads,(int)patches.nPatches());

    std::ostringstream oss;
    oss << nThreads;

    setenv("OMP_NUM_THREADS", oss.str().c_str(),1);
#endif


    /////////////////// Setup boundary conditions ///////////////////
    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    gsFunctionExpr<> hEast, hSouth,hNorth, hWest;

    hWest = gsFunctionExpr<>("-pi*cos(pi*0.4)*sin(2*pi*(y+0.3))-1",dim);
    hSouth = gsFunctionExpr<>("-pi*2*cos(2*pi*0.3)*sin(pi*(x+0.4))-1",dim);

    switch(testcase){
    case 1:
        hEast = gsFunctionExpr<>("1*pi*cos(pi*1.4)*sin(2*pi*(y+0.3))+1",dim);

        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(0.5+0.3))*sin(pi*(x+0.4))+1",dim);

        bcInfo.addCondition(1, boundary::east, condition_type::neumann, hEast);

        bcInfo.addCondition(0, boundary::north, condition_type::neumann, hNorth);
        bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);

        bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);
        bcInfo.addCondition(1, boundary::south, condition_type::neumann, hSouth);

        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);

        break;
    case 2:
        hEast = gsFunctionExpr<>("1*pi*cos(pi*1.4)*sin(2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(1+0.3))*sin(pi*(x+0.4))+1",dim);

        //     bcInfo.addCondition(2, boundary::east, condition_type::neumann, hEast);
        //     bcInfo.addCondition(3, boundary::east, condition_type::neumann, hEast);

        bcInfo.addCondition(2, boundary::east,  condition_type::dirichlet, g);
        bcInfo.addCondition(3, boundary::east,  condition_type::dirichlet, g);

        //     bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);
        //    bcInfo.addCondition(3, boundary::north, condition_type::neumann, hNorth);

        //     bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);
        //    bcInfo.addCondition(2, boundary::south, condition_type::neumann, hSouth);
        bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, g);
        bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, g);

        bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, g);
        bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, g);


        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
        bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);
        break;
    case 3:
        hEast = gsFunctionExpr<>("1*pi*cos(pi*(3+0.4))*sin(2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(1+0.3))*sin(pi*(x+0.4))+1",dim);


        bcInfo.addCondition(2, boundary::east, condition_type::neumann, hEast);

        bcInfo.addCondition(0, boundary::north, condition_type::neumann, hNorth);
        bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);
        bcInfo.addCondition(2, boundary::north, condition_type::neumann, hNorth);

        bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);
        bcInfo.addCondition(1, boundary::south, condition_type::neumann, hSouth);
        bcInfo.addCondition(2, boundary::south, condition_type::neumann, hSouth);

        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
        break;
    case 4:
        hEast = gsFunctionExpr<>("1*pi*cos(pi*(1+0.4))*sin(2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(2+0.3))*sin(pi*(x+0.4))+1",dim);

        bcInfo.addCondition(0, boundary::east, condition_type::neumann, hEast);
        bcInfo.addCondition(1, boundary::east, condition_type::neumann, hEast);

        bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);

        bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);

        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
        bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);
        break;
    case 5:
        hEast = gsFunctionExpr<>("1*pi*cos(pi*(1+0.4))*sin(2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(2+0.3))*sin(pi*(x+0.4))+1",dim);


        bcInfo.addCondition(0, boundary::east, condition_type::neumann, hEast);
        bcInfo.addCondition(1, boundary::east, condition_type::neumann, hEast);

        bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);

        bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);

        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
        bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);
        //bcInfo.addCondition(1, boundary::west,  condition_type::neumann, &hWest);
        break;
    case 6:

        hEast = gsFunctionExpr<>("1*pi*cos(pi*(4+0.4))*sin(2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(2+0.3))*sin(pi*(x+0.4))+1",dim);

        bcInfo.addCondition(6, boundary::east, condition_type::neumann, hEast);
        bcInfo.addCondition(7, boundary::east, condition_type::neumann, hEast);

        bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);
        bcInfo.addCondition(3, boundary::north, condition_type::neumann, hNorth);
        bcInfo.addCondition(5, boundary::north, condition_type::neumann, hNorth);
        bcInfo.addCondition(7, boundary::north, condition_type::neumann, hNorth);

        bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);
        bcInfo.addCondition(2, boundary::south, condition_type::neumann, hSouth);
        bcInfo.addCondition(4, boundary::south, condition_type::neumann, hSouth);
        bcInfo.addCondition(6, boundary::south, condition_type::neumann, hSouth);

        bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);
        // bcInfo.addCondition(0, boundary::west,  condition_type::neumann, hWest);
        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
        break;

    case 7:
    case 10:
        //bcInfo.addCondition(0,boundary::east, condition_type::dirichlet, &g);

        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
            if((*it).patch == 4)
                bcInfo.addCondition(*it, condition_type::dirichlet, g);
        /*
        bcInfo.addCondition(0,boundary::west, condition_type::dirichlet, g);
        bcInfo.addCondition(5,boundary::north, condition_type::dirichlet, g);
        bcInfo.addCondition(6,boundary::south, condition_type::dirichlet, g);
        bcInfo.addCondition(10,boundary::east, condition_type::dirichlet, g);
        */
        break;
    case 8:
    case 9:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
            bcInfo.addCondition(*it, condition_type::dirichlet, g);
        break;
    default:
        hEast = gsFunctionExpr<>("1*pi*cos(pi*(2+0.4))*sin( primalDofMethod::Y2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(1+0.3))*sin(pi*(x+0.4))+1",dim);

        bcInfo.addCondition(1, boundary::east, condition_type::neumann, hEast);

        bcInfo.addCondition(0, boundary::north, condition_type::neumann, hNorth);
        bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);

        bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);
        bcInfo.addCondition(1, boundary::south, condition_type::neumann, hSouth);

        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
        break;
    }


    if(test3D)
    {
        gsMultiPatch<> patches3D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&patches.patch(i));

            patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb));
        }
        if(testcase==8)
        {
            mp= patches3D.uniformSplit();
            std::vector<gsGeometry<>*> ptch1, ptch2;
            for(size_t np = 0; np< mp.nPatches(); ++np)
            {
                mp.patch(np).degreeElevate(1,2);
                std::vector<gsGeometry<>* >  ptch_ =  (mp).patch(np).uniformSplit(1);
                ptch1.insert(ptch1.end(),ptch_.begin(),ptch_.end());
            }
            for(size_t np = 0; np< mp.nPatches(); ++np)
            {
                std::vector<gsGeometry<>* >  ptch_ =  ptch1[np]->uniformSplit(2);
                ptch2.insert(ptch2.end(),ptch_.begin(),ptch_.end());
            }
            freeAll(ptch1);
            patches = gsMultiPatch<>(ptch2); //just change for compatibility
            patches.computeTopology();
        }
        else
        {
            patches3D.computeTopology();
            patches = give(patches3D);
        }

    }

    int result = 0;
    if(plot)
    {
        gsWriteParaview<>( patches, "IETI_patch", 1000, true);
            result = system("paraview IETI_patch.pvd");
    }
    ////////////////////// Refinement h and p //////////////////////
    // Refinement

    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        refine_bases.degreeIncrease(numElevate);
        gsInfo<<"Degree set to: "<<refine_bases.maxDegree(0)<<"\n";
    }

    // h-refine each basis
    for (int i = 0; i < numRefine; ++i)
        refine_bases.uniformRefine();

    ///////////////////////ASSEMBLE////////////////////////////////////

    gsInfo<<bcInfo<<"\n";

    switch(testcase)
    {
    case 2:
        for(size_t np = 0; np<patches.nPatches();np++)
            if(np == 1 || np == 2)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
        break;
    case 8:
        for(size_t np = 0; np<patches.nPatches();np++)
            if ((np / 4 + np / 8) % 2 == 0)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
        break;

    case 9:
        for(size_t np = 0; np<patches.nPatches();np++)
            if ((np / 8 + np / 16 + np / 32) % 2 == 0)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
        break;
    default:
        for(size_t np = 0; np<patches.nPatches();np++)
            if(np%2 == 0)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
        break;
    }

    // gsPoissonPde<real_t>ppde(*patches,bcInfo,f);
    // gsPoissonAssembler<real_t>assembler (ppde,refine_bases,dstrat);

    gsPoissonHeterogeneousPde<real_t>ppde(patches,bcInfo,f,alpha);
    gsPoissonHeterogeneousAssembler<real_t>assembler (ppde,refine_bases,dstrat,fstrat);

    option.update(assembler.options());

    gsStopwatch time;

    //IETI Setup:
    /*
    gsIETIOptions opt;
    opt.extraTimings=timings;
    opt.strat.strat = pstrat;
    opt.isMinimumEnergy= MinEnergy;
    opt.scal =  IETIPrecondScaling::coefficient;
    opt.solvingSaddlePoint = saddlePoint;
   // opt.KCSolver = chooseSolver(solvKC);
   // opt.KiiSolver = chooseSolver(solvI);

    damping = damping == -1? getDefaultDamping(smoother,refine_bases.maxDegree(0)) : damping;

    if(opt.KCSolver == IETILocalSolver::Multigrid)
    {
        opt.KC_optionsMG.numLevels = numRefine;
        opt.KC_optionsMG.damping = damping;
        opt.KC_optionsMG.smoother = smoother;
        opt.KC_optionsMG.coarseDamping = coarseDamping;
        opt.KC_optionsMG.numPreSmooth = 1;
        opt.KC_optionsMG.numPostSmooth = 1;
        if(!opt.solvingSaddlePoint)
        {
            opt.KC_optionsMG.max_iter = 99;
            opt.KC_optionsMG.tol = 1.e-10;
        }
        else
        {
            opt.KC_optionsMG.max_iter = dim == 2 ? 2 : 3;
            opt.KC_optionsMG.tol = 1.e-2;
        }
        //fill up here

    }
    else if(opt.KCSolver == IETILocalSolver::FastDiagonalization)
    {
        if(!opt.solvingSaddlePoint)
        {
            opt.KC_optionsBS.max_iter = 10000;
            opt.KC_optionsBS.tol = 1.e-10;
        }
        else
        {
            opt.KC_optionsBS.max_iter = 2;
            opt.KC_optionsBS.tol = 1.e-2;
        }
    }
    //--------------------------
    if(opt.KiiSolver == IETILocalSolver::Multigrid)
    {
        opt.Kii_optionsMG.numLevels = numRefine;
        opt.Kii_optionsMG.damping = damping;
        opt.Kii_optionsMG.coarseDamping = coarseDamping;
        opt.Kii_optionsMG.smoother = smoother;
        opt.Kii_optionsMG.max_iter = 2;
        opt.Kii_optionsMG.tol = 1.e-2;
        opt.Kii_optionsMG.numPreSmooth = 1;
        opt.Kii_optionsMG.numPostSmooth = 1;
        //fill up here
    }
    else if(opt.KiiSolver == IETILocalSolver::FastDiagonalization)
    {

        opt.Kii_optionsBS.max_iter = 3;
        opt.Kii_optionsBS.tol = 1.e-2;

    }
    */
    gsIETIAssembler<real_t>* iass;
    if(dG)
        iass = new gsIETIdGAssembler<real_t>(assembler);
    else
        iass = new gsIETIAssembler<real_t>(assembler);
    iass->setOptions(option.getGroup("IETI"));
    iass->init();

    time.restart();
    iass->assemble();
    real_t t1 = time.stop();

    gsInfo<<"\n";
    gsInfo << "Time for assembling: "<< t1<<"\n"<<std::flush;

    gsMatrix<> isolVector, isolution;
    unsigned it=0;
    real_t t2;
    if(!saddlePoint)
    {
        isolVector.setZero(iass->systemSize(),iass->numberRhs());
        if(testcase ==10)
            isolVector = gsMatrix<real_t>::Constant(isolVector.rows(), isolVector.cols(), 1);

        gsIETISolver<real_t>::Ptr isolv = gsIETISolver<real_t>::make(*iass);
        isolv->init();
        gsScaledDirichletPrecond<real_t>::Ptr iprecDir =
                gsScaledDirichletPrecond<real_t>::make(*iass);
        gsConjugateGradient<> iPCG(isolv,iprecDir);
        iPCG.setOptions(option);

        time.restart();
        iPCG.solve(isolv->getRhs(),isolVector);
        t2 = time.stop();
        gsInfo<< "Time for Solving: "<<t2<<"\n"<<std::flush;

        gsMatrix<> ieigs;
        iPCG.getEigenvalues(ieigs);
        gsInfo<<"Eigenvalues: "<<ieigs.minCoeff() <<" - "<<ieigs.maxCoeff()<<"\n";
        //   gsInfo<<"eigenvalues: "<<"\n"<<eigs<<"\n";
        //   gsInfo<<"condition number: "<<"\n"<<PCG.getConditionNumber()<<"\n";

        gsInfo<<"number of CG iterations: " <<iPCG.iterations()<<"\n";
        gsInfo<<"residual error: " <<iPCG.error()<<"\n";
        gsInfo<<"\n";
        iass->printTiming();
        it = iPCG.iterations();
        isolv->calculateSolution(isolVector,isolution);
    }
    else
    {
        isolVector.setZero(iass->systemSize(),iass->numberRhs());
        if(testcase ==10)
            isolVector = gsMatrix<real_t>::Constant(isolVector.rows(), isolVector.cols(), 1);
        gsInexactIETISolver<real_t>::Ptr isolv = gsInexactIETISolver<real_t>::make(*iass);
        isolv->init();
        gsInexactIETIPrecond<real_t>::Ptr iprecDir = gsInexactIETIPrecond<real_t>::make(*iass);
        gsLinearOperator<real_t>::Ptr opK = iprecDir->getPrecForK();
        gsLinearOperator<real_t>::Ptr opS = iprecDir->getPrecForS();
        //gsMinimalResidual<> MinRes(isolv,iprecDir);
        gsBramblePasciakCG<gsBPCG_Types::BP_Schur> MinRes(isolv->makeBlockOp(),opK,opS);
        MinRes.setOptions(option);
        //      gsMatrix<real_t> F,Prec;
        //     isolv.calculateMatrixForm(F);
        //     iprecDir.calculateMatrixForm(Prec);
        
        //   gsInfo<<"F:\n"<<F<<"\n";
        //   gsInfo<<"Prec:\n"<<Prec<<"\n";
        //   gsInfo<<"PF:\n"<<Prec*F<<"\n";
        //   gsInfo<<"rhs:\n "<<isolv.getRhs()<<"\n";
        time.restart();
	gsMatrix<real_t> errHist;
        MinRes.solveDetailed(isolv->getRhs(),isolVector,errHist);
	gsInfo<<errHist<<"\n\n";
        //  MinRes.solve(isolv.getRhs(),isolVector);
        t2 = time.stop();
        gsInfo<< "Time for Solving: "<<t2<<"\n"<<std::flush;

        gsInfo<<"number of MinRes iterations: " <<MinRes.iterations()<<"\n";
        gsInfo<<"residual error: " <<MinRes.error()<<"\n";
        gsInfo<<"\n";
        iass->printTiming();
        it = MinRes.iterations();
        //     gsSolverOp< Eigen::FullPivLU< Eigen::MatrixXd > >::Ptr solv = makeFullPivLUSolver(F);
        //    solv->apply(isolv.getRhs(),isolVector);

        isolv->calculateSolution(isolVector,isolution);
    }

    gsInfo<<"Number of Dofs: "<<assembler.numDofs()<<"\n";
    gsField<> isol = assembler.constructSolution(isolution);
    real_t il2error = isol.distanceL2(g);

    gsInfo<<"L2 error ||u_iIETI - g||_L2: "<<il2error<<"\n";

    //Comparism with exact IETI --------------------------------------------------------------
    /*
    opt.KCSolver = chooseSolver("D");
    opt.KiiSolver = chooseSolver("D");
    opt.solvingSaddlePoint = false;

    gsIETIAssembler<real_t> ass(assembler, opt);
    ass.init();

    time.restart();
    ass.assemble();
    time.stop();

    gsInfo<<"\n";
    gsInfo << "Time for assembling: "<< time<<"\n"<<std::flush;

    gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(ass);
    solv->init();


            gsScaledDirichletPrecond<real_t>::Ptr precDir =
                gsScaledDirichletPrecond<real_t>::make(ass);
    gsConjugateGradient<> PCG(solv, precDir);
    
    PCG.setMaxIterations(50);
    PCG.setTolerance(1.e-6);
    PCG.setCalcEigenvalues(true);
    gsMatrix<> solVector(ass.systemSize(),ass.numberRhs());
    if(testcase ==10)
        solVector = gsMatrix<real_t>::Constant(solVector.rows(), solVector.cols(), 1);
    else
        solVector.setZero(solVector.rows(), solVector.cols());

    time.restart();
    PCG.solve(solv->getRhs(),solVector);
    time.stop();
    gsInfo<< "Time for Solving: "<<time<<"\n";

    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo<<"Eigenvalues: "<<eigs.minCoeff() <<" - "<<eigs.maxCoeff()<<"\n";
    //   gsInfo<<"eigenvalues: "<<"\n"<<eigs<<"\n";
    //   gsInfo<<"condition number: "<<"\n"<<PCG.getConditionNumber()<<"\n";

    gsInfo<<"number of CG iterations: " <<PCG.iterations()<<"\n";
    gsInfo<<"residual error: " <<PCG.error()<<"\n";
    gsInfo<<"\n";



    gsMatrix<> solution;
    solv->calculateSolution(solVector,solution);
    gsField<> SolEx = assembler.constructSolution(solution);
  //  real_t l2error = sol->distanceL2(g);
 //   gsInfo<<"L2 error ||u_iIETI - g||_L2: "<<l2error<<"\n";

    // gsInfo<<"lambda: "<<"\n"<<solVector<<"\n";
    //  gsInfo<<"solIETI: \n"<< solution.transpose()<<"\n\n";
    //   gsInfo<<"solDirect: \n"<< solEx.transpose()<<"\n\n";
    gsMatrix<> diff = (isolution-solution);
    
  //  gsField<> SolEx = (assembler.constructSolution(sol));
    assembler.homogenizeFixedDofs();
    gsField<> Diff = (assembler.constructSolution(diff));
     real_t l2error2 = SolEx.distanceL2(g);
      real_t l2error3 = SolEx.distanceL2(isol);
     gsInfo<<"L2 error ||u_direct - g||: "<<l2error2<<"\n";
     gsInfo<<"L2 error ||u_direct - u_IETI||: "<<l2error3<<"\n";
*/


    // Plot solution in paraview
    if(out!="")
    {
        std::string filename(out + ".dat");
        std::fstream output(filename.c_str(), std::fstream::app|std::fstream::out);
        if(!output.is_open())
            GISMO_ERROR("Opening file failed");
        output<<"#dofs: "<<assembler.numDofs() <<"    ,   "<<iass->getInfo().dofTotalP<<"   ,    "<<iass->getInfo().numberPatches<<"   ,    "<< it<<"   ,    "<<iass->getTiming().totalAssemble.sum()<<"   ,   "<<t1-iass->getTiming().totalAssemble.sum()<<"   ,   "<<t2<<"   ,   "<<t1+t2<<
                " ||  "<<iass->getTiming().averageItNumKCBasis.sum()/iass->getTiming().numKCItBasis.sum()<<" |  "<< iass->getTiming().averageItNumKCSpace.sum()/iass->getTiming().numKCItSpace.sum()<<
                "  |  "<<iass->getTiming().averageItNumKii.sum()/iass->getTiming().numKiiIt.sum()<<"\n";

        output.close();
    }

    if (plot)
    {
        gsPiecewiseFunction<real_t> iter;
        for(size_t np = 0; np<patches.nPatches();np++)
            iter.addPiecePointer(new gsConstantFunction<real_t>(iass->getTiming().averageItNumKCBasis[np]/iass->getTiming().numKCItBasis[np],dim));
        std::string str = "IETI_iter_"+util::to_string(refine_bases.maxDegree(0))+"_"+util::to_string(numRefine)+"_"+util::to_string(damping)+"_"+util::to_string(coarseDamping);
        const gsField<> ITER( patches, iter, false );
        gsWriteParaview<>( ITER, str, 1000);
        /*
        std::vector<gsFunction<>*> fun(patches->nPatches());
        for(size_t np = 0; np <patches->nPatches();np++)
            fun[np] = &alpha[np];
        const gsField<> alph( *patches, fun, false );
        gsWriteParaview<>( alph, "IETI_alpha", 1000);
        */

        //    gsWriteParaview<>( *Diff, "IETI_diff", 1000);
        //    gsWriteParaview<>( *SolEx, "IETI_SolEx", 1000, true);

        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>(isol, "IETI", 1000);
        const gsField<> exact( patches, g, false );
        gsWriteParaview<>( exact, "IETI_ex", 1000);
        gsField<> alph( patches, alpha, false );
        gsWriteParaview<>( alph, "IETI_alpha", 1000);

        // Run paraview
        result = system(("paraview "+str+".pvd  &").c_str());

        //  result = system("paraview IETI_ex.pvd  &");
    }

    delete iass;
    return result;
}
