/**  gsParallelIETI.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):      C. Hofer
    Created on:  2014-12-03

    This file tests the MPI version of the PCG method with IETI

    WHAT IS HAPPENING

    EXPECTED RESULT
*/

#include <gismo.h>
#include <gismo_dev.h>

#ifdef GISMO_WITH_MPI
using namespace gismo;

#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>

#include <fstream>

#include <gsIETI/gsParallelCG.h>
#include <gsIETI/gsIETIAssemblerMPI.h>
#include <gsIETI/gsIETIdGAssemblerMPI.h>
#include <gsIETI/gsIETISolverMPI.h>
#include <gsIETI/gsIETIScaledDirichletMPI.h>

#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETIdGAssembler.h>
#include <gsIETI/gsIETIScaledDirichlet.h>
#include <gsIETI/gsIETISolver.h>

/// A locally used geometry
gsMultiPatch<real_t> approximateQuarterAnnulus(int deg)
{
    gsGeometry<>::uPtr quann = gsNurbsCreator<>::NurbsQuarterAnnulus();
    gsKnotVector<> KV1(0,1, 0, deg+1);        // no interior knots in x direction
    gsKnotVector<> KV2(0,1, 0, deg+1);        // no interior knot in y direction

    gsTensorBSplineBasis<2> tbsp (new gsBSplineBasis<>(KV1), new gsBSplineBasis<>(KV2));
    gsGeometry<>::uPtr approxGeom = tbsp.interpolateAtAnchors( quann->eval(tbsp.anchors()) );
    gsMultiPatch<real_t> mp(*approxGeom);
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

    gsTensorBSplineBasis<3> tbsp (new gsBSplineBasis<>(KV1), new gsBSplineBasis<>(KV2), new gsBSplineBasis<>(KV3));
    gsGeometry<>::uPtr approxGeom = tbsp.interpolateAtAnchors( mp2->eval(tbsp.anchors()) );
    gsMultiPatch<real_t> mp(*approxGeom);
    return mp;
}


int main(int argc, char **argv)
{


    int testcase =1;
    int max_ref = 2;
    int max_elev = 0;
    int max_inc = 0;
    int jump = 0;
    real_t dom_size =1;
    // bool plot = false;
    bool timings = false;
    bool test3D=false;
    bool comp=false;
    bool dg= false;
    std::string out_str="";
    bool out = false;
    bool plot = false;
    bool mixedInterfaces = false;
    int result=1;
    int n = 8;
    int nSppHolder=1;
    std::string m = "B";
    std::string scaling = "coeff";


    iFace::strategy face = iFace::glue;

    //#if HAVE_MPI
    gsMpi& mpi = gsMpi::init(argc, argv);
    gsMpiComm comm(mpi.worldComm());

    gsCmdLine cmd("Hi, I will test the convergence behaviour of IETI");
    cmd.addSwitch("d", "dim", test3D);
    cmd.addSwitch("comp", "comparism with serial IETI",comp);
    cmd.addSwitch("dg", "use dg",dg);
    cmd.addString("o", "out","file output",out_str);
    cmd.addInt("r","refine","Number of refinements",max_ref);
    cmd.addInt("c","case","Choosen testcase",testcase);
    cmd.addInt("e","elevate","maximal elevation",max_elev);
    cmd.addInt("i","increase","maximal elevation (keeping multiplicity)",max_inc);
    cmd.addInt("j","jumps","jumping coeff from 10^-j to 10^j",jump);
    cmd.addInt("n", "number","number of domain in one direction",n);


    cmd.addSwitch("plot","plot the domain and error", plot);
    cmd.addSwitch("",  "IETI.ExtraTiming", "enable extraTimings", timings);
    cmd.addSwitch("",  "IETI.EnableMixedInterfaceCoupling", "enables the method to handle mixed interface coupling strategies", mixedInterfaces);
    cmd.addString("s", "IETI.Strategy", "Choosen strategy",m);
    cmd.addString("",  "IETI.Scaling", "chosen scaling",scaling);
    cmd.addInt   ("p", "IETI.NumberSppHolder","number of Spp holder",nSppHolder);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsOptionList option = cmd.getOptionList();

    if(nSppHolder==0)
        nSppHolder=1;
    else if(nSppHolder<0)
        nSppHolder=comm.size();

    if(!plot)
        result--; //prevent stupid unused warning

    face = dg ? iFace::dg : iFace::glue;

    int dim;
    test3D ? dim =3:dim =2;
    if(testcase == 8) dim =3;
    if(testcase == 8) test3D=false;
    std::stringstream ss;
    ss << jump;
    std::string jump_str = ss.str();

    if(out_str == "")
        out = false;
    else
        out= true;

    real_t a,b,c;
    // a=70;
    // b=100;
    a=4;
    b=2;
    c=3;

    std::string str("v:="+internal::to_string(a)+";w:="+internal::to_string(b)+";u:="+internal::to_string(c)+";");
    gsFunctionExpr<> f,g;
    if(dim==2)
    {
        f = gsFunctionExpr<>(str+"((pi*v)^2 + (pi*w)^2)*sin(v*pi*(x+0.4))*sin(w*pi*(y+0.3))",dim);
        g = gsFunctionExpr<>(str+"sin(v*pi*(x+0.4))*sin(pi*(y+0.3)*w)+x+y",dim);
    }
    else
    {
        f = gsFunctionExpr<>(str+"((pi*v)^2 + (pi*w)^2+ (pi*u)^2)*sin(v*pi*(x+0.4))*sin(w*pi*(y+0.3))*sin(u*pi*(z+0.8))",dim);
        g = gsFunctionExpr<>(str+"sin(v*pi*(x+0.4))*sin(pi*(y+0.3)*w)*sin(u*pi*(z+0.8))+x+y",dim);
    }

    gsMultiPatch<> patches;
    if(comm.rank()==0)
    {
        gsInfo<<"Using jumping coeff: 10^"+jump_str+" - 10^-"+jump_str<<"\n";
        gsInfo<<"Using method "+m<<"\n";
        if(dg || mixedInterfaces)
            gsInfo<<"Going with dG fomulation\n";
        else
            gsInfo<<"Going with cG fomulation\n";
        gsInfo<<"Using "<<nSppHolder<<" processors for holding the Spp\n";

        gsInfo<< "Right hand side: "<<f<<"\n";
        gsInfo<< "Solution: "<<g<<"\n";
    }

    switch(testcase){
    case 1:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,2, dom_size);
        break;
    case 2:
    default:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,2, dom_size);
        break;
    case 3:
        if(dim==2)
            patches = gsNurbsCreator<>::BSplineSquareGrid(n,n, 1./n);
        else
            patches = gsNurbsCreator<>::BSplineCubeGrid(n,n,2*n,1./n);
        break;
    case 4:
    {
        gsFileData<> fileData("yeti_mp2.xml");
        if (!fileData.getFirst(patches))
            return 1;
        //patches->computeTopology();
        break;
    }
    case 5:
    {
        gsFileData<> fileData("gaps/chairStandTopFace.xml");
        if (!fileData.getFirst(patches))
            return 1;
        patches.computeTopology();
        break;
    }
    case 6:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,1, dom_size);
        break;
    case 7:
    {
        if(false && n==32)
        {
            gsFileData<> fileData("quarter_annulus_1024p.xml");
            if (!fileData.getFirst(patches))
                return 1;
        }
        else
        {
            gsMultiPatch<real_t> mp= approximateQuarterAnnulus(2);

            for(int i=0; i<log2(n);++i)
                mp= mp.uniformSplit();
            patches = give(mp);
            patches.computeTopology(1.e-6);
        }

        break;
    }
    case 8:
    {
        if(n==8)
        {
            gsFileData<> fileData("twisted_quarter_annulus_1024p.xml");
            if (!fileData.getFirst(patches))
                return 1;
        }
        else
        {
            gsMultiPatch<> mp = approximateTwistedQuarterAnnulus(2);

            std::vector<gsGeometry<>*> ptch;
            for(size_t np = 0; np< mp.nPatches(); ++np)
            {
                std::vector<gsGeometry<>* >  ptch_ =  (mp).patch(np).uniformSplit(1);

                ptch.insert(ptch.end(),ptch_.begin(),ptch_.end());
            }

            mp= gsMultiPatch<real_t>(ptch);

            for(int i=0; i<log2(n);++i)
                mp= mp.uniformSplit();
            patches = give(mp);
            patches.computeTopology();
        }
        break;
    }



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


    gsPiecewiseFunction<real_t>  alpha;
    gsFunctionExpr<> a1("1.e-"+jump_str,dim);
    gsFunctionExpr<> a2("1.e+"+jump_str,dim);
    if(testcase==7)
    {
        /*
        for(index_t np = 0; np<patches.nPatches();np++)
            if((np/(patches.nPatches()/4)) == 0 || (np/(patches.nPatches()/4)) == 3)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
        */

        for(size_t np = 0; np<patches.nPatches();np++)
            if((np/(patches.nPatches()/4)) == 0 || (np/(patches.nPatches()/4)) == 3)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
    }
    else if(testcase ==8 )
    {
        for(size_t np = 0; np<patches.nPatches();np++)
            if((((unsigned)np/8)+((unsigned)np/16)+((unsigned)np/32))%2  == 0)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
    }
    else
    {
        for(size_t np = 0; np<patches.nPatches();np++)
            if(np%2  == 0)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
    }

    //Handle interface strategies
    std::vector<iFace::strategy> iFaceStrategies;
    if(mixedInterfaces)
    {
        for(gsMultiPatch<real_t>::iiterator it = patches.iBegin(); it!=patches.iEnd();++it)
        {
            boundaryInterface bI = *it;
            real_t val1= (alpha.piece(bI.first().patch).eval(gsMatrix<real_t>::Constant(dim,1,0)))(0,0);
            real_t val2 = (alpha.piece(bI.second().patch).eval(gsMatrix<real_t>::Constant(dim,1,0)))(0,0);
            if(val1!=val2)
                iFaceStrategies.push_back(iFace::dg);
            else
                iFaceStrategies.push_back(iFace::conforming);
        }
    }

    /////////////////// Setup boundary conditions ///////////////////
    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    gsFunctionExpr<> hEast, hSouth,hNorth, hWest;

    hWest = gsFunctionExpr<>(str+"-v*pi*cos(v*pi*0.4)*sin(w*pi*(y+0.3))-1",dim);
    hSouth = gsFunctionExpr<>(str+"-pi*w*cos(w*pi*0.3)*sin(v*pi*(x+0.4))-1",dim);

    switch(testcase)
    {
    case 1:

        hEast = gsFunctionExpr<>(str+"v*pi*cos(v*pi*(4+0.4))*sin(w*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>(str+"pi*w*cos(w*pi*(2+0.3))*sin(v*pi*(x+0.4))+1",dim);

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
        bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);
        break;
    case 2:

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
        bcInfo.addCondition(1, boundary::west,  condition_type::neumann, hWest);
        break;
    case 3:

        hEast = gsFunctionExpr<>(str+"v*pi*cos(v*pi*("+internal::to_string(1)+"+0.4))*sin(w*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>(str+"pi*w*cos(w*pi*("+internal::to_string(1)+"+0.3))*sin(v*pi*(x+0.4))+1",dim);

        for(int i=0;i<n;i++)
        {
            bcInfo.addCondition(i+n*(n-1), boundary::east, condition_type::neumann, hEast);
            bcInfo.addCondition(i*n+n-1, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(i*n, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(i, boundary::west,  condition_type::dirichlet, g);
        }
        break;
    case 4:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
            // if((*it).patch==1)
            bcInfo.addCondition(*it, condition_type::dirichlet, g);
        break;
    case 5:
        bcInfo.addCondition(0,boundary::east, condition_type::dirichlet, g);
        break;
    case 6:

        hEast = gsFunctionExpr<>("1*pi*cos(pi*(3+0.4))*sin(2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(3+0.3))*sin(pi*(x+0.4))+1",dim);

        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
            bcInfo.addCondition(*it, condition_type::dirichlet, g);

        break;
    default:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
            bcInfo.addCondition(*it, condition_type::dirichlet, g);
        break;
    }


    if(test3D && (testcase!=3 || testcase !=8))
    {
        // gsInfo<<"Making Domain 3D"<<"\n";
        gsMultiPatch<real_t> patches3D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&patches.patch(i));

            patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb));
        }
        patches3D.computeTopology(1.e-4,true);

        patches = give(patches3D);

    }
    else
    {
        //  gsInfo<<"2D Domain"<<"\n";
    }
    if(comm.rank()==0 && plot)
    {
        gsWriteParaview<>( patches, "IETI_patch", 1000, true);
        result = system("paraview IETI_patch.pvd &");
    }
    ///////////////////// Refinement h and p //////////////////////
    // Refinement

    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );
    // Number for h-refinement of the computational (trail/test) basis.



    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( max_elev > -1 || max_inc> -1 )
    {
        // Find maximum degree with respect to all the variables
        int max_tmp = refine_bases.maxDegree(0);
        for (int j = 1; j < patches.parDim(); ++j )
            if ( max_tmp < refine_bases.maxDegree(j) )
                max_tmp = refine_bases.maxDegree(j);

        if(comm.rank()==0)
            gsInfo<<"Degree set to: "<<max_tmp+max_inc<<"\n";
        refine_bases.setDegree(max_tmp);

        refine_bases.degreeIncrease(max_inc);

    }

    // h-refine each basis (4, one for each patch)
    for (int i = 0; i < max_ref; ++i)
        refine_bases.uniformRefine();

    // refine_bases.reduceContinuity(max_elev);

    ///////////////////////ASSEMBLE////////////////////////////////////

    gsPoissonHeterogeneousPde<real_t>ppde(patches,bcInfo,f,alpha);
    gsPoissonHeterogeneousAssembler<real_t>* assembler;

    if(mixedInterfaces)
        assembler = new gsPoissonHeterogeneousAssembler<real_t>(ppde,refine_bases,dirichlet::elimination,iFaceStrategies);
    else
        assembler = new gsPoissonHeterogeneousAssembler<real_t>(ppde,refine_bases,dirichlet::elimination,face);

    //This class should parse all the available patchwise classes for the parallelization
    gsIETIAssemblerMPI<real_t>* ass_ptr;

    if(dg || mixedInterfaces)
        ass_ptr = new gsIETIdGAssemblerMPI<real_t>(*assembler);
    else
        ass_ptr = new gsIETIAssemblerMPI<real_t>(*assembler);

    ass_ptr->setOptions(option.getGroup("IETI"));
    ass_ptr->init();

    //Assemble

    comm.barrier();

    double t0 = MPI_Wtime();
    ass_ptr->assemble();
    //  gsMPIHelper::getCollectiveCommunication().barrier();
    //time2.stop();
    double t1 = MPI_Wtime();

    ////////////////////SOLVE//////////////////////////////////////////

    // gsInfo<<"Init IETISolver "<<mpiHelper.rank()<<"\n"<<std::flush;
    gsIETISolverMPI<real_t>::Ptr solv = gsIETISolverMPI<real_t>::make(*ass_ptr);
    solv->init();
    //   gsInfo<<"Init Precond "<<mpiHelper.rank()<<"\n"
    gsScaledDirichletPrecondMPI<real_t>::Ptr precDir = gsScaledDirichletPrecondMPI<real_t>::make(*ass_ptr);

    //   gsInfo<<"Init PCG "<<mpiHelper.rank()<<"\n";
    //   gsParallelCG<real_t> PCG( memory::make_shared_not_owned(&solv), precDir);
    gsParallelCG<real_t> PCG(solv, precDir);
    PCG.setMaxIterations(200);
    PCG.setTolerance(1.e-8); //1.e-11
    gsMatrix<> solVector = gsMatrix<>::Zero(ass_ptr->getInfoMPI().lagrangeMultReduce,ass_ptr->numberRhs());


    if(comm.rank() == 0)
        PCG.setCalcEigenvalues(true);

    //PCG2.solve(solv2.getRhs(),solVector,precDir);

    //  solVector.setZero();
    //   gsInfo<<"start solving "<<mpiHelper.rank()<<"\n nDofs"<<refine_bases[ass_ptr->m_patchIdx[0]].size()<<"\n";
    comm.barrier();
    double t2 = MPI_Wtime();
    PCG.solve(solv->getRhs(),solVector);
    double t3 = MPI_Wtime();
    //  gsInfo<<"Done solving "<<mpiHelper.rank()<<"\n";

    gsMatrix<> solution;
    solv->calculateSolution(solVector,solution);
    ass_ptr->combineToCommonSolution(solution);
    // gsInfo<<"Done combining "<<mpiHelper.rank()<<"\n";
    gsField<> sol ;
    if(comm.rank()==0 && out)
    {
        sol = assembler->constructSolution(solution);
        real_t l2error, h1error;
        l2error=h1error=0;

        //l2error= sol.distanceL2(g);
        gsInfo<<"L2 error ||u_IETI_MPI - g||: "<<l2error<<"\n";
        // h1error = sol.distanceH1(g);
        gsInfo<<"H1 error ||u_IETI_MPI - g||: "<<h1error<<"\n";

        gsInfo<<" nLag: "<<ass_ptr->systemSize()<<"\n nDofs: "<<ass_ptr->getInfo().dofTotal<<"\n nPrimalDofs: "<<ass_ptr->getInfo().dofTotalP<<"\n";
        gsInfo<<"Assembling time: "<< t1-t0<<"\nSolving time: "<<t3-t2<<"\n"<<"Solver Iterations: "<<PCG.iterations()<<"\n\n"<<std::flush;
        std::string filename(out_str + ".dat");
        if(dg)
            filename = out_str +"_dg"+ ".dat";
        std::fstream output(filename, std::fstream::app|std::fstream::out);
        if(!output.is_open())
            GISMO_ERROR("Opening file failed");
        output<<comm.size()<<" , "<<t1-t0<<" , "<<t3-t2<<", "<<ass_ptr->getInfo().dofTotal<<", "<<l2error<<", "<<h1error<<" , " <<  PCG.iterations()<<"\n";
        output.close();
    }


    ////////////////////BUILD SOLUTION//////////////////////////////////////////
    gsMatrix<> eigs;

    if(comm.rank() == 0)
    {

        //    gsInfo<<(solVector.transpose())<<"\n\n";
        //    gsInfo<<(solution.transpose())<<"\n\n";

        //    gsField<>::uPtr sol = memory::make_unique(assembler.constructSolution(solution));
        //    real_t l2error = sol->distanceL2(g);
        //     gsInfo<<"L2 error ||u_IETI_MPI - g||: "<<l2error<<"\n";


        if(comp)
        {
            gsIETIAssembler<real_t>* ass2_ptr;
            if(dg || mixedInterfaces )
                ass2_ptr = new gsIETIdGAssembler<real_t>(*assembler);
            else
                ass2_ptr = new gsIETIAssembler<real_t>(*assembler);

            gsIETIAssembler<real_t>& ass2 = *ass2_ptr;
            ass2.setOptions(option.getGroup("IETI"));
            ass2.init();
            ass2.assemble();
            gsIETISolver<real_t> solv2(ass2);
            solv2.init();
            gsScaledDirichletPrecond<real_t>::Ptr precDir2 =
                    gsScaledDirichletPrecond<real_t>::make(ass2);
            gsConjugateGradient<> PCG2(memory::make_shared_not_owned(&solv2), precDir2);
            PCG2.setMaxIterations(100);
            PCG2.setTolerance(1.e-8);

            gsMatrix<> solVector2 = gsMatrix<>::Zero(ass2.getInfo().lagrangeMult,ass2.numberRhs());
            PCG2.solve(solv2.getRhs(),solVector2);
            gsMatrix<> solution2;
            solv2.calculateSolution(solVector2,solution2);
            gsInfo<<(solVector2.transpose())<<"\n\n";
            gsInfo<<(solution2.transpose())<<"\n\n";
            gsField<> sol2 = assembler->constructSolution(solution2);

            gsInfo<<"Solver Iterations: "<<PCG2.iterations()<<"\n\n"<<std::flush;
            gsInfo<<"finished comp\n";

            real_t norm=0;
            for(size_t pn =0; pn<patches.nPatches();pn++)
            {
                gsMultiPatch<> mpPatch(patches.patch(pn));
                gsMultiPatch<> solMP(sol.igaFunction(pn));
                gsField<> solMP_field(mpPatch, solMP);

                gsNormL2<real_t> L2Norm(solMP_field, sol2.igaFunction(pn), true);
                norm+=L2Norm.compute();
            }
            gsInfo<<"L2 error ||u_IETI_MPI - u_IETI||: "<<norm<<"\n";


            gsMatrix<> diff = solution2 - solution;
            assembler->homogenizeFixedDofs(-1);
            gsField<> diff_ptr = assembler->constructSolution(diff);
            if(plot)
            {
                gsWriteParaview<>(diff_ptr, "IETI_MPI_diff", 1000);
                gsFileManager::open("IETI_MPI_diff.pvd");
            }

            delete ass2_ptr;
            //   gsWriteParaview<>(*sol2, "IETI", 1000);
        }

        // Write approximate and exact solution to paraview files
        if(plot)
        {
            gsInfo<<"Plotting in Paraview...\n";
            //      gsWriteParaview<>(sol, "IETI_MPI", 1000);
            const gsField<> exact( patches, g, false );
            //     gsWriteParaview<>( exact, "IETI_ex", 1000);

            gsPiecewiseFunction<real_t>  procs;
            for(size_t np = 0; np<patches.nPatches();++np)
                procs.addPiecePointer(new gsFunctionExpr<>(util::to_string(ass_ptr->m_patch2proc[np]),dim));

            gsField<> proc( patches, procs, false );
            //   gsWriteParaview<>( proc, "IETI_procs", 1000);

            gsField<> alph( patches, alpha, false );
            gsWriteParaview<>( alph, "IETI_alpha", 1000);


            // Run paraview
            //  result = system("paraview IETI_MPI.pvd & ");

            //   result = system("paraview IETI_ex.pvd & ");
        }
        gsInfo<<"Done Everything\n";
    }
    delete ass_ptr;
    delete assembler;
    return 0;
}

#else
int main(int argc, char **argv)
{
    gsInfo<<" No Mpi enabled! \n";
    return 0;
}

#endif
