/**  gsIETI_experimental.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):     C. Hofer
    Created on:  2014-12-03

    DESCRIBE THE PURPOSE OF THE EXAMPLE HERE

    WHAT IS HAPPENING

    EXPECTED RESULT
*/

#include <iostream>
#include <gismo.h>
#include <gismo_dev.h>
#include <gsIETI/gsIETIUtils.h>
#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIScaledDirichlet.h>

#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>

#include <gsSolver/gsConjugateGradient.h>


using namespace gismo;


//Important Note: The boundary condition has to be called with
// bcInfo.addCondition(***, *** , *** , g);
// in order that the parallelization works correctly! So no &g !!!
// Otherwise, you have a data race!


int main (int argc, char** args)
{
    index_t testcase =1;
    bool test3D = false;
    bool plot = false;
    bool nitsche = false;
    bool timings = false;

    dirichlet::strategy dstrat = dirichlet::nitsche;
    bool NoMinEnergy = false;
    std::string m = "A";
    std::string scaling = "coeff";

    // Number for h-refinement of the computational (trail/test) basis.
    index_t numRefine  = 1;
    // Number for p-refinement of the computational (trail/test) basis.
    index_t numElevate = 0;

    gsCmdLine cmd("Hi, I will test IETIdG");
    cmd.addSwitch("d", "dim", "test3D", test3D);
    cmd.addSwitch("n", "nitsche", "Use Nitsche if chosen; otherwise elimination", nitsche);
    cmd.addSwitch("p", "plot", "Plot solution", plot);
    cmd.addSwitch("",  "IETI.ExtraTiming", "enable extraTimings", timings);
    cmd.addSwitch("",  "IETI.NoMinimumEnergy","choose energy minimizing prim subspaces",NoMinEnergy);
    cmd.addInt("r","refine","Number of refinements",numRefine);
    cmd.addInt("c","case","Choosen testcase",testcase);
    cmd.addInt("e","elevate","maximal elevation",numElevate);
    cmd.addString("s", "IETI.Strategy", "Choosen strategy",m);
    cmd.addString("",  "IETI.Scaling", "chosen scaling",scaling);
    try { cmd.getValues(argc,args); } catch (int rv) { return rv; }

    gsOptionList option = cmd.getOptionList();
    //test3D = true;
    //nitsche = true;

    if(nitsche)
        dstrat = dirichlet::nitsche;
    else
        dstrat = dirichlet::elimination;

    int dim;
    (test3D ||testcase == 8) ? dim =3:dim =2;
    gsInfo<<"Dimension of Domain: "<<dim<<"\n";
    gsInfo<<"Using Nitsche: "<<nitsche<<"\n";

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
        f=gsFunctionExpr<>("((pi*1)^2 + (pi*2)^2)*sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)",dim);
        g=gsFunctionExpr<>("sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)+x+y",dim);
    }
    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n" << "\n";

    /////////////////// Setup geometry ///////////////////
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;


    switch(testcase){
    case 1:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,1, 0.5);
        break;
    case 2:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,2, 0.5);
        break;
    case 3:
        patches = gsNurbsCreator<>::BSplineSquareGrid(3,1, 1);
        break;
    case 4:
    case 5:
        patches = gsNurbsCreator<>::BSplineSquareGrid(1,2, 1);
        break;
    case 6:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,2, 1);
        break;
    case 7:
    case 10:
    {
        gsFileData<> fileData("yeti_mp2.xml");
        //gsFileData<> fileData("gaps/chairStandTopFace.xml");
        if (fileData.has< gsMultiPatch<> >())
        {
            fileData.getFirst(patches);
        }
        else
            return 1;
        // patches->computeTopology();
        break;
    }
    case 8:
    {
        gsFileData<> fileData("yeti3d_u21p.xml");
        if (fileData.has< gsMultiPatch<> >())
        {
            fileData.getFirst(patches);
        }
        else
            return 1;
        break;
    }
    default:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,1, 1);
        break;
    }


#ifdef _OPENMP
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

        bcInfo.addCondition(2, boundary::east, condition_type::neumann, hEast);
        bcInfo.addCondition(3, boundary::east, condition_type::neumann, hEast);

        bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);
        bcInfo.addCondition(3, boundary::north, condition_type::neumann, hNorth);

        bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);
        bcInfo.addCondition(2, boundary::south, condition_type::neumann, hSouth);

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
        bcInfo.addCondition(1, boundary::west,  condition_type::neumann, hWest);
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

        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
        //bcInfo.addCondition(1, boundary::west,  condition_type::neumann, hWest);
        bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);
        break;

    case 7:
    case 8:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
           // if((*it).patch == 1)
                bcInfo.addCondition(*it, condition_type::dirichlet, g);
        /*
        bcInfo.addCondition(0,boundary::west, condition_type::dirichlet, g);
        bcInfo.addCondition(5,boundary::north, condition_type::dirichlet, g);
        bcInfo.addCondition(6,boundary::south, condition_type::dirichlet, g);
        bcInfo.addCondition(10,boundary::east, condition_type::dirichlet, g);
        */
        break;
    default:
        hEast = gsFunctionExpr<>("1*pi*cos(pi*(2+0.4))*sin(2*pi*(y+0.3))+1",dim);
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
        gsMultiPatch<real_t>  patches3D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&patches.patch(i));

            patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb));
        }
        patches3D.computeTopology();

        patches = give(patches3D);

    }



    ////////////////////// Refinement h and p //////////////////////
    // Refinement

    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int max_tmp = refine_bases.maxDegree(0);
        for (int j = 1; j < patches.parDim(); ++j )
            if ( max_tmp < refine_bases.maxDegree(j) )
                max_tmp = refine_bases.maxDegree(j);

        // Elevate all degrees uniformly
        max_tmp += numElevate;

        gsInfo<<"Degree set to: "<<max_tmp<<"\n";
        refine_bases.setDegree(max_tmp);
    }

    // h-refine each basis (4, one for each patch)
    for (int i = 0; i < numRefine; ++i)
        refine_bases.uniformRefine();


    ///////////////////////ASSEMBLE////////////////////////////////////

    gsInfo<<bcInfo<<"\n";


    gsPiecewiseFunction<real_t>  alpha;
    gsFunctionExpr<> a1("1.e-0",dim);
    gsFunctionExpr<> a2("1.e+0",dim);
    for(size_t np = 0; np<patches.nPatches();np++)
        if(np%2 == 0)
            alpha.addPiece(a1);
        else
            alpha.addPiece(a2);


    gsPoissonHeterogeneousPde<real_t> pde(patches, bcInfo,f,alpha);
    gsPoissonHeterogeneousAssembler<real_t> assembler(pde,refine_bases,dstrat);
    assembler.system().rhs().conservativeResize(Eigen::NoChange, 1); // Now necessary
    gsPoissonPde<real_t>ppde(patches,bcInfo,f);
   // gsPoissonAssembler<real_t>assembler (ppde,refine_bases,dstrat);
    gsPoissonAssembler<real_t> assembler2(patches,refine_bases,bcInfo,f, dstrat);
    assembler2.system().rhs().conservativeResize(Eigen::NoChange, 1); // Now necessary
    gsStopwatch time;

    // For comparism, the exact solution

    time.restart();
    assembler.assemble();
    gsSparseMatrix<> K = assembler.fullMatrix();
    K.makeCompressed();

#if defined(GISMO_WITH_PARDISO)
    gsSparseSolver<>::PardisoLU solver;
#elif defined(GISMO_WITH_SUPERLU)
    gsSparseSolver<>::SuperLU solver;
#else
    gsSparseSolver<>::LU solver;
#endif

    solver.compute(K);
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eiK(K.toDense());
    //gsInfo<<"eigenvalues K: "<<eiK.eigenvalues().transpose()<<"\n";
    gsMatrix<> solEx = solver.solve( assembler.rhs() );
    time.stop();

    gsInfo<<"Time direct solve: "<<time<<"\n";

    assembler2.assemble();
    gsSparseMatrix<> K2 = assembler2.fullMatrix();
    K2.makeCompressed();

#if defined(GISMO_WITH_PARDISO)
    gsSparseSolver<>::PardisoLU solver2;
#elif defined(GISMO_WITH_SUPERLU)
    gsSparseSolver<>::SuperLU solver2;
#else
    gsSparseSolver<>::LU solver2;
#endif

    solver2.compute(K);
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eiK(K.toDense());
    //gsInfo<<"eigenvalues K: "<<eiK.eigenvalues().transpose()<<"\n";
    gsMatrix<> solEx2 = solver2.solve( assembler2.rhs() );

    gsInfo<<"K Diff: "<< (K-K2).norm()<<"\n";
    gsInfo<<"rhs Diff: "<< (assembler.rhs()-assembler2.rhs()).norm()<<"\n";
    gsInfo<<"solEx Diff: "<< (solEx-solEx2).norm()<<"\n";

    time.restart();

    gsIETIAssembler<real_t> ass(assembler);
    ass.setOptions(option.getGroup("IETI"));
    ass.init();
    std::vector<real_t>  scal;
    for(size_t np = 0; np<patches.nPatches();np++)
        if(np%2 == 0)
            scal.push_back(1.e-2);
        else
            scal.push_back(1.e+2);

    ass.assemble();
    time.stop();


    gsIETIInfo inf = ass.getInfo();
    inf.print();

    gsInfo<<"\n";
    gsInfo << "Time for assembling: "<< time<<"\n"<<std::flush;

    gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(ass);
    solv->init();
    gsScaledDirichletPrecond<real_t>::Ptr precDir =
        gsScaledDirichletPrecond<real_t>::make(ass);
    gsConjugateGradient<> PCG(solv,precDir);
    PCG.setMaxIterations(100);
    PCG.setCalcEigenvalues(true);

    gsMatrix<> solVector(ass.systemSize(),ass.numberRhs());
    if(testcase == 10)
        solVector = gsMatrix<real_t>::Constant(solVector.rows(), solVector.cols(), 1);
    else
        solVector.setZero(solVector.rows(), solVector.cols());


    gsMatrix<> rhs = solv->getRhs();


    //---------------------------
    /*
    //matrix forms for debugging
    gsMatrix<> F,Md;
    solv.calculateMatrixForm(F);
    precDir.calculateMatrixForm(Md);

    gsMatrix<> B,BT;
    //solv.calculateMatrixForm(B,BT);
    //gsInfo<<B<<"\n"<<"\n";
    //gsInfo<<BT<<"\n";

    //gsIETIJumpOperator<real_t> jump(ass,IETIPrecondScaling::stiffness);
    //jump.calculateMatrixForm(B,BT);

    //gsInfo<<B<<"\n"<<"\n";
    //gsInfo<<BT<<"\n";

    gsMatrix<> MdF = Md*F;
    Eigen::EigenSolver<Eigen::Matrix<real_t,Dynamic,Dynamic> > ei(MdF);
    gsInfo<<"eigenvalues MdF: "<<ei.eigenvalues().transpose()<<"\n";

    //Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real_t,Dynamic,Dynamic> > eiM(Md);
    //gsInfo<<"eigenvalues Md: "<<eiM.eigenvalues().transpose()<<"\n";


    //Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real_t,Dynamic,Dynamic> > eiF(F);
    //gsInfo<<"eigenvalues F: "<<eiF.eigenvalues().transpose()<<"\n";


    //gsInfo << "F: "<<"\n"<<F<<"\n"<<"\n";
   // gsInfo << "Md: "<<"\n"<<Md<<"\n"<<"\n";
    //gsInfo << "MdF: "<<"\n"<<MdF<<"\n"<<"\n";

*/
    //---------------------------


    time.restart();
    PCG.solve(rhs,solVector);
    time.stop();

    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);

    gsInfo<<"eigenvalues: "<<"\n"<<eigs<<"\n";

    //  gsInfo<<"condition number: "<<"\n"<<PCG.getConditionNumber()<<"\n";

    gsInfo<< "Time for Solving: "<<time<<"\n";

    gsInfo<<"number of CG iterations: " <<PCG.iterations()<<"\n";
    gsInfo<<"residual error: " <<PCG.error()<<"\n";
    gsInfo<<"\n";

    //gsInfo<<"lambda: "<<"\n"<<solVector<<"\n";


    gsMatrix<> solution;
    solv->calculateSolution(solVector,solution);
    gsField<> sol = assembler.constructSolution(solution);
    real_t l2error = sol.distanceL2(g);
    gsInfo<<"L2 error ||u_IETI - g||: "<<l2error<<"\n";

    gsMatrix<> diff = (solution-solEx);
    gsField<> SolEx = assembler.constructSolution(solEx);
    assembler.homogenizeFixedDofs();
    gsField<> Diff = assembler.constructSolution(diff);
    real_t l2error2 = SolEx.distanceL2(g);
    real_t l2error3 = SolEx.distanceL2(sol);
    gsInfo<<"L2 error ||u_direct - g||: "<<l2error2<<"\n";
    gsInfo<<"L2 error ||u_direct - u_IETI||: "<<l2error3<<"\n";

    // Plot solution in paraview
    int result = 0;
    if (plot)
    {
        /*
        std::vector<gsFunction<>* > fun(patches.nPatches());
        for(index_t np = 0; np <patches.nPatches();np++)
            fun[np] = const_cast<gsFunction<>*>(&alpha[np]);
        const gsField<> alph( patches, fun, false );
        gsWriteParaview<>( alph, "IETI_alpha", 1000);
*/

        gsWriteParaview<>( Diff, "IETI_diff", 1000);
        gsWriteParaview<>( SolEx, "IETI_SolEx", 1000);

        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>(sol, "IETI", 1000, true);
        const gsField<> exact( patches, g, false );
        gsWriteParaview<>( exact, "IETI_ex", 1000);

        // Run paraview
        result = system("paraview IETI.pvd  &");

        result = system("paraview IETI_diff.pvd  &");
    }

    //return  1 for failures and 0 for success
    return result;
}


