/**  gsIETI_testConvergence.cpp
  
    This file is a test of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s):     C. Hofer
    Created on:  2015-01-20
    
    DESCRIBE THE PURPOSE OF THE EXAMPLE HERE
    
    WHAT IS HAPPENING
    
    EXPECTED RESULT
*/

#include <iomanip>
#include <iostream>
#include <gismo.h>
#include <gismo_dev.h>
#include <gsIETI/gsIETIUtils.h>
#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIScaledDirichlet.h>

#include <gsSolver/gsConjugateGradient.h>


#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>

using namespace gismo;

template<typename T> void printField(T t, const int& width, int prec = 3)
{
    gsInfo << std::left << std::setw(width) << std::setfill(' ') <<std::setprecision(prec)<< t;
}
void printElement(unsigned globDof, unsigned lagMult, unsigned locDofs, real_t Hh,const gsMatrix<unsigned>& iter,const gsMatrix<real_t>& cond)
{
    printField(globDof,8);
    printField(lagMult,6);
    printField(locDofs,6);
    printField( cast<real_t,int>(Hh),6);
    for(index_t i=0; i<iter.cols();i++)
        printField(iter(0,i),12);

    for(index_t i=0; i<cond.cols();i++)
        printField(cond(0,i),12);

    gsInfo<<"\n";

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

//Important Note: The boundary condition has to be called with
// bcInfo.addCondition(***, *** , *** , g);
// in order that the parallelization works correctly! So no &g !!!
// Otherwise, you have a data race!

enum testType
{
    hScaling = 0,
    pScaling,
    weakScaling
};

int main (int argc, char** args)
{
    index_t testcase =4;
    index_t max_ref = 2;
    index_t max_elev = 0;
    index_t max_inc = 0;
    index_t jump = 0;
    index_t max_patchRef = 1;
    real_t dom_size =1;
    index_t nPrec =4;
    index_t testType=0;
    bool plot = false;
    bool timings = false;
    bool test3D=false;

    std::string m = "B";
    // std::string scaling = "coeff";


    gsCmdLine cmd("Hi, I will test the convergence behaviour of IETI");
    cmd.addSwitch("d", "dim", "test3D", test3D);
    cmd.addInt("","testType","h-ref test (0), p-ref test (1) or h/H test(2)",testType);
    cmd.addSwitch("",  "IETI.ExtraTiming", "enable extraTimings", timings);
    cmd.addInt("r","refine","Number of refinements",max_ref);
    cmd.addInt("n","patchRefine", "number of uniform patch splittings", max_patchRef);
    cmd.addInt("c","case","Choosen testcase",testcase);
    cmd.addInt("e","elevate","maximal elevation",max_elev);
    cmd.addInt("i","increase","maximal elevation (keeping multiplicity)",max_inc);
    cmd.addInt("j","jumps","jumping coeff from 10^-j to 10^j",jump);
    cmd.addInt("", "nPrec","number of tests scalings",nPrec);
    cmd.addString("s", "IETI.Strategy", "Choosen strategy",m);
    //  cmd.addString("",  "IETI.Scaling", "chosen scaling",scaling);
    try { cmd.getValues(argc,args); } catch (int rv) { return rv; }

    gsOptionList option = cmd.getOptionList();

    int dim;
    test3D ? dim =3:dim =2;
    std::stringstream ss;
    ss << jump;
    std::string jump_str = ss.str();
    gsInfo<<"Using jumping coeff: 10^"+jump_str+" - 10^-"+jump_str<<"\n";
    gsInfo<<"Using method "+m+"\n";
    gsStopwatch time;
    //////// Right-hand side and analytical solution ////////
    // Define source function

    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)",dim);
    gsFunctionExpr<> g("sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)+x+y",dim);

    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n" << "\n";

    /////////////////// Setup geometry ///////////////////
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;

    switch(testcase){
    case 1:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,2, dom_size);
        break;
    case 2:
    default:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,2, dom_size);
        break;
    case 3:
        patches = gsNurbsCreator<>::BSplineSquareGrid(8,8, dom_size);
        break;
    case 4:
    {
        gsFileData<> fileData("yeti_mp2.xml");
        if (fileData.has< gsMultiPatch<> >())
        {
            fileData.getFirst(patches);
        }
        else
            return 1;
        break;
    }
    case 5:
    {
        gsFileData<> fileData("gaps/chairStandTopFace.xml");
        if (fileData.has< gsMultiPatch<> >())
            fileData.getFirst(patches);
        else
            return 1;
        patches.computeTopology();
        break;
    }
    case 6:
        patches = gsNurbsCreator<>::BSplineSquareGrid(3,3, dom_size);
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

    gsPiecewiseFunction<real_t>  alpha;
    gsFunctionExpr<> a1("1.e-"+jump_str,dim);
    gsFunctionExpr<> a2("1.e+"+jump_str,dim);
    for(size_t np = 0; np<patches.nPatches();np++)
        if(np%2 == 0)
            alpha.addPiece(a1);
        else
            alpha.addPiece(a2);



    /////////////////// Setup boundary conditions ///////////////////
    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    gsFunctionExpr<> hEast, hSouth,hNorth, hWest;

    hWest = gsFunctionExpr<>("-pi*cos(pi*0.4)*sin(2*pi*(y+0.3))-1",dim);
    hSouth = gsFunctionExpr<>("-pi*2*cos(2*pi*0.3)*sin(pi*(x+0.4))-1",dim);

    switch(testcase)
    {
    case 1:

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
        bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);
        break;
    case 2:
    default:

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

        hEast = gsFunctionExpr<>("1*pi*cos(pi*(8+0.4))*sin(2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(8+0.3))*sin(pi*(x+0.4))+1",dim);

        for(int i=0;i<8;i++)
        {
            bcInfo.addCondition(i+56, boundary::east, condition_type::neumann, hEast);
            bcInfo.addCondition(i*8+7, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(i*8, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(i, boundary::west,  condition_type::dirichlet, g);
        }
        break;
    case 4:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
            if((*it).patch==1)
                bcInfo.addCondition(*it, condition_type::dirichlet, g);
        break;
    case 5:
        bcInfo.addCondition(0,boundary::east, condition_type::dirichlet, g);
        break;
    case 6:

        hEast = gsFunctionExpr<>("1*pi*cos(pi*(3+0.4))*sin(2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(3+0.3))*sin(pi*(x+0.4))+1",dim);

        for(int i=0;i<3;i++)
        {
            bcInfo.addCondition(i+2*3, boundary::east, condition_type::neumann, hEast);
            bcInfo.addCondition(i*3+2, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(i*3, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(i, boundary::west,  condition_type::dirichlet, g);
        }
        break;

    }

    if(test3D)
    {
        gsInfo<<"Making Domain 3D"<<"\n";
        gsMultiPatch<real_t> patches3D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&patches.patch(i));

            patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb));
        }
        patches3D.computeTopology();

        patches = give(patches3D);

    }
    else
    {
        gsInfo<<"2D Domain"<<"\n";
    }

    int result = 0;
    if(plot)
    {
        gsWriteParaview<>( patches, "IETI_patch", 1000, true);
        gsWriteParaview<>( patches, "IETI_patch2_", 1000);
        result = system("paraview IETI_patch.pvd &");
    }

    gsMatrix<> eigs,solVector,rhs,solution;
    gsField<> sol;

    gsInfo<<testType<<"\n";
    if(testType == hScaling)
    {
        ////////////////////// Refinement h and p //////////////////////
        // Refinement


        // Copy basis from the geometry
        gsMultiBasis<> refine_bases( patches );
	

        // Elevate and p-refine the basis to order k + numElevate
        // where k is the highest degree in the bases
        if ( max_elev > -1 || max_inc> -1 )
        {
            // Find maximum degree with respect to all the variables
            int max_tmp = refine_bases.maxDegree(0);
            for (int j = 1; j < patches.parDim(); ++j )
                if ( max_tmp < refine_bases.maxDegree(j) )
                    max_tmp = refine_bases.maxDegree(j);

            // Elevate all degrees uniformly
            max_tmp += max_elev;

            gsInfo<<"Degree set to: "<<max_tmp+max_inc<<"\n";
	    for ( size_t i = 0; i < refine_bases.nBases(); ++ i )
	       refine_bases[i].setDegreePreservingMultiplicity(max_tmp);

            refine_bases.degreeIncrease(max_inc);
        }
        
        if(test3D)
            refine_bases.uniformRefineComponent(2);


        gsMatrix<unsigned> iter(max_ref,nPrec);
        gsMatrix<real_t> l2Error(max_ref,nPrec+1);
        gsMatrix<real_t> Hh(max_ref,1);
        gsMatrix<unsigned> globDof(max_ref,1);
        gsMatrix<unsigned> lagMult(max_ref,1);
        gsMatrix<unsigned>  locDof(max_ref,1);

        gsMatrix<real_t> cond(max_ref,nPrec);
        cond.setZero();
        iter.setZero();
        l2Error.setZero();


        for(int r = 0;r < max_ref;r++)
        {
            gsInfo<<"refing "<<r<<"\n";
            //gsPoissonPde<real_t>ppde(patches,bcInfo,f);
            //gsPoissonAssembler<real_t>assembler (ppde,refine_bases,dirichlet::elimination);

            gsPoissonHeterogeneousPde<real_t>ppde(patches,bcInfo,f,alpha);
            gsPoissonHeterogeneousAssembler<real_t>assembler (ppde,refine_bases,dirichlet::elimination);

            /*
            //for the comparism of the convergence rate
            assembler.assemble();
            gsSparseMatrix<> K = assembler.fullMatrix();
            gsSparseSolver<>::LU solver( K );
            gsMatrix<> solEx = solver.solve( assembler.rhs() );
            gsField<> SolEx = assembler.constructSolution(solEx);
            l2Error(r,nPrec) = SolEx.distanceL2(g);
            */

            time.restart();
            gsIETIAssembler<real_t> ass(assembler);
            ass.setOptions(option.getGroup("IETI"));

            ass.init();
            time.stop();
            gsInfo<<"Time for preparing the bookkeeping: "<<time<<"\n";
	    
	    gsIETIInfo info = ass.getInfo();
	    if(!test3D)
                Hh(r,0) = math::sqrt( (real_t)(refine_bases.size(0)) )-1; //Hh domsize cancels
            else
                Hh(r,0) = math::pow( (real_t)(refine_bases.size(0)), 1.0/3)-1;

            locDof(r,0) = refine_bases.size(0);
            globDof(r,0)=info.dofTotal;
            lagMult(r,0)=info.lagrangeMult;
	    
	    
	     if(nPrec == 1) {refine_bases.uniformRefine();continue;}
	    
            time.restart();
            ass.assemble();
            time.stop();
            gsInfo<<"Time for assembling everything: "<<time<<"\n";

            gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(ass);
            solv->init();


            for(int pre = 1;pre<nPrec;pre++)
            {
                gsInfo<<"starting preconditioner "<<pre<<"\n";

                solVector.setZero(ass.systemSize(),ass.numberRhs());

                rhs = solv->getRhs();


                // gsConjugateGradient<> PCG(solv,precDir);
                //  PCG.setMaxIterations(200);
                //  PCG.setCalcEigenvalues(true);
                IETIPrecondScaling::strategy preEnum = static_cast<IETIPrecondScaling::strategy>(pre+1);
                ass.getOptions().scal = preEnum;
                gsScaledDirichletPrecond<real_t>::Ptr precDir = gsScaledDirichletPrecond<real_t>::make(ass);
                gsConjugateGradient<> PCG(solv,precDir);
                PCG.setCalcEigenvalues(true);
                PCG.setMaxIterations(200);
                PCG.setTolerance(1.e-8);

                /*
                if(pre == 0 && false)
                {
                    continue;
                    PCG.solve(rhs,solVector);

                    iter(r,pre)=PCG.iterations();
                    PCG.getEigenvalues(eigs);
                    cond(r,pre)=eigs.maxCoeff()/eigs.minCoeff();
                    continue;
                }
                //*/
                //  IETIPrecondScaling::strategy preEnum = static_cast<IETIPrecondScaling::strategy>(pre+1);
                //  ass.getOptions().scal = preEnum;
                time.restart();
                PCG.solve(rhs,solVector);
                time.stop();
                gsInfo<<"Time for solving PCG total: "<<time<<"\n";
                gsInfo<<"Time per PCG-Iteration: "<<time<<"/"<<PCG.iterations()<<"\n";

                iter(r,pre)=PCG.iterations();

                PCG.getEigenvalues(eigs);
                cond(r,pre)=eigs.maxCoeff()/eigs.minCoeff();

                time.restart();
                //solv->calculateSolution(solVector,solution);
                time.stop();
                gsInfo<<"Time for reconstructing solution: "<<time<<"\n";

                //sol = assembler.constructSolution(solution);
                //l2Error(r,pre)= sol.distanceL2(g);
            }

            refine_bases.uniformRefine();


        }

        gsMatrix<real_t> l2Factor(max_ref,l2Error.cols());
        l2Factor.setZero();
        for(int i=1;i<max_ref;i++)
            for(int j=0;j<l2Error.cols();j++)
                l2Factor(i,j) = math::log(l2Error(i-1,j)/l2Error(i,j))/math::log(2.0);

        /*
        gsInfo << "\n";
        printField("gDof",8);printField("lagM",6);printField("lDof",6); printField("H/h",6); printField("iF",10); printField("iPF-non",10); printField("iPF-mult",10);printField("iPF-coeff",10);printField("iPF-stiff",10);printField("cF",12); printField("cPF-non",12); printField("cPF-mult",12);printField("cPF-coeff",12);printField("cPF-stiff",12);
        gsInfo<<"\n";
        for(int i=0;i<max_ref;i++)
            printElement(globDof(i,0),lagMult(i,0),locDof(i,0),Hh(i,0),iter.row(i),cond.row(i));
*/
        gsInfo << "\n";
        printField("gDof",8);printField("lagM",6);printField("lDof",6); printField("H/h",6); printField("iF",12);printField("iPF-coeff",12);printField("iPF-stiff",12);printField("iPF-stiffM",12);printField("cF",12); printField("cPF-coeff",12);printField("cPF-stiff",12);printField("cPF-stiffM",12);
        gsInfo<<"\n";
        for(int i=0;i<max_ref;i++)
            printElement(globDof(i,0),lagMult(i,0),locDof(i,0),Hh(i,0),iter.row(i),cond.row(i));


        gsInfo<<"\n";
        gsInfo << l2Factor <<"\n"<<"\n";
        gsInfo<<"\n";
        gsInfo <<l2Error<<"\n"<<"\n";


    }
    else if(testType == pScaling)
    {
        int max_degRef= math::max(max_elev,max_inc);
        // Copy basis from the geometry
        gsMultiBasis<> refine_bases( patches );
        //std::vector<gsBasis<> *> refine_bases = patches.basesCopy();

        // Number for h-refinement of the computational (trail/test) basis.
        int numRefine  = max_ref;
        // Number for p-refinement of the computational (trail/test) basis.
        // int numElevate = 0;

        // h-refine each basis (4, one for each patch)


        // Elevate and p-refine the basis to order k + numElevate
        // where k is the highest degree in the bases
        if ( max_elev > -1 || max_inc> -1 )
        {
            // Find maximum degree with respect to all the variables
            int max_tmp = refine_bases.maxDegree(0);
            for (int j = 1; j < patches.parDim(); ++j )
                if ( max_tmp < refine_bases.maxDegree(j) )
                    max_tmp = refine_bases.maxDegree(j);

            // Elevate all degrees uniformly
            max_tmp += max_elev;

            gsInfo<<"Degree set to: "<<max_tmp+max_inc<<"\n";
            for ( size_t i = 0; i < refine_bases.nBases(); ++ i )
	       refine_bases[i].setDegreePreservingMultiplicity(max_tmp);

            refine_bases.degreeIncrease(max_inc);
        }
        if(test3D)
            refine_bases.uniformRefineComponent(2);

        for (int i = 0; i < numRefine; ++i)
            refine_bases.uniformRefine();

        gsMatrix<unsigned> iter(max_degRef,nPrec);
        gsMatrix<real_t> l2Error(max_degRef,1);
        gsMatrix<unsigned> globDof(max_degRef,1);
        gsMatrix<unsigned> lagMult(max_degRef,1);
        gsMatrix<unsigned>  locDof(max_degRef,1);
        gsMatrix<real_t> Hh(max_degRef,1);
        gsMatrix<real_t> deg(max_degRef,1);
        gsMatrix<real_t> cond(max_degRef,nPrec);
        cond.setZero();
        iter.setZero();

        for(int r = 0;r < max_degRef;r++)
        {
            gsInfo<<"refing degree "<<r<<"\n";

            gsPoissonHeterogeneousPde<real_t>ppde(patches,bcInfo,f,alpha);
            gsPoissonHeterogeneousAssembler<real_t>assembler (ppde,refine_bases,dirichlet::elimination);

            time.restart();
            gsIETIAssembler<real_t> ass(assembler);
            ass.setOptions(option.getGroup("IETI"));

            ass.init();
            time.stop();
            //   gsInfo<<"Time for preparing the bookkeeping: "<<time<<"\n";
            time.restart();
            ass.assemble();
            time.stop();
            //    gsInfo<<"Time for assembling everything: "<<time<<"\n";

            gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(ass);
            solv->init();
            gsIETIInfo info = ass.getInfo();

            locDof(r,0) = refine_bases.size(0);
            globDof(r,0)=info.dofTotal;
            lagMult(r,0)=info.lagrangeMult;
            deg(r,0) = refine_bases.degree();

            if(!test3D)
                Hh(r,0) = math::sqrt( (real_t)(refine_bases.size(0)) )-1; //Hh domsize cancels
            else
                Hh(r,0) = math::pow( (real_t)(refine_bases.size(0)), 1.0/3)-1;

            if(max_elev>max_inc)
                refine_bases.degreeElevate(1);
            else
                refine_bases.degreeIncrease(1);

            for(int pre = 1;pre<nPrec;pre++)
            {
                solVector.setZero(ass.systemSize(),ass.numberRhs());
                rhs = solv->getRhs();


                gsConjugateGradient<> PCG(solv);
                PCG.setCalcEigenvalues(true);
                PCG.setMaxIterations(200);
                //*
                if(pre == 0 && false)
                {
                    PCG.solve(rhs,solVector);

                    iter(r,pre+1)=PCG.iterations();

                    PCG.getEigenvalues(eigs);
                    cond(r,pre+1)=eigs.maxCoeff()/eigs.minCoeff();
                    continue;
                }
                //*/
                IETIPrecondScaling::strategy preEnum = static_cast<IETIPrecondScaling::strategy>(pre+1);
                ass.getOptions().scal = preEnum;
                gsScaledDirichletPrecond<real_t>::Ptr precDir = gsScaledDirichletPrecond<real_t>::make(ass);

                time.restart();
                PCG.setPreconditioner(precDir);
                PCG.solve(rhs,solVector);
                time.stop();
                gsInfo<<"Time for solving PCG total: "<<time<<"\n";
                iter(r,pre)=PCG.iterations();
                PCG.getEigenvalues(eigs);
                cond(r,pre)=eigs.maxCoeff()/eigs.minCoeff();

                if(pre==nPrec-1 && false)
                {

                    solv->calculateSolution(solVector,solution);
                    sol = assembler.constructSolution(solution);
                    l2Error(r,0)= sol.distanceL2(g);
                }

            }



        }


        gsInfo << "\n";
        printField("gDof",8);printField("lagM",6);printField("lDof",6); printField("deg",6); printField("iF",12);printField("iPF-coeff",12);printField("iPF-stiff",12);printField("iPF-stiffM",12);printField("cF",12); printField("cPF-coeff",12);printField("cPF-stiff",12);printField("cPF-stiffM",12);
        gsInfo<<"\n";
        for(int i=0;i<max_degRef;i++)
            printElement(globDof(i,0),lagMult(i,0),locDof(i,0),deg(i,0),iter.row(i),cond.row(i));

        gsInfo<<"\n";
        gsInfo <<l2Error<<"\n"<<"\n";

    }
    else if(testType == weakScaling)
    {
        gsMatrix<unsigned> iter(max_patchRef,nPrec);
        gsMatrix<real_t> l2Error(max_patchRef,1);
        gsMatrix<unsigned> globDof(max_patchRef,1);
        gsMatrix<unsigned> lagMult(max_patchRef,1);
        gsMatrix<unsigned>  nPatches(max_patchRef,1);
        gsMatrix<real_t> Hh(max_patchRef,1);
        gsMatrix<real_t> dofsP(max_patchRef,1);
        gsMatrix<real_t> cond(max_patchRef,nPrec);
        cond.setZero();
        iter.setZero();

        // Number for h-refinement of the computational (trail/test) basis.
        int numRefine  = max_ref;

        int d = patches.parDim();
        for(int r = 0;r < max_patchRef;r++)
        {
	    if(r!=0)
	    {
                patches = patches.uniformSplit();   
	        bcInfo = updateBoundaryConditionsAfterSplit(bcInfo,patches);
	    }
	    
            gsPiecewiseFunction<real_t>  alphaN;
            for(size_t np = 0; np<patches.nPatches();np++)
                if((np/math::ipow(2*d,(unsigned)r))%2 == 0)
                    alphaN.addPiece(a1);
                else
                    alphaN.addPiece(a2);


            // Copy basis from the geometry
            gsMultiBasis<> refine_bases( patches );

            // Elevate and p-refine the basis to order k + numElevate
            // where k is the highest degree in the bases
	    if ( max_elev > -1 || max_inc> -1 )
	    {
		// Find maximum degree with respect to all the variables
		int max_tmp = refine_bases.maxDegree(0);
		for (int j = 1; j < patches.parDim(); ++j )
		    if ( max_tmp < refine_bases.maxDegree(j) )
			max_tmp = refine_bases.maxDegree(j);

		// Elevate all degrees uniformly
		max_tmp += max_elev;

		gsInfo<<"Degree set to: "<<max_tmp+max_inc<<"\n";
                for ( size_t i = 0; i < refine_bases.nBases(); ++ i )
	           refine_bases[i].setDegreePreservingMultiplicity(max_tmp);

		refine_bases.degreeIncrease(max_inc);
	    }

            if(test3D && r==0)
                refine_bases.uniformRefineComponent(2);
            if(r != 0)
                refine_bases.uniformRefine();

            for (int i = 0; i < numRefine; ++i)
                refine_bases.uniformRefine();

            gsInfo<<"refing patch "<<r<<"\n";

            //BC not handled properly TODO!

	    /*
            for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
                for(int i=0; i<math::ipow(2*d,(unsigned)r);++i)
                    if((*it).patch==math::ipow(2*d,(unsigned)r)+i)
                        bcInfoN.addCondition(*it, condition_type::dirichlet, g);
	    */

            gsPoissonHeterogeneousPde<real_t>ppde(patches,bcInfo,f,alphaN);
            gsPoissonHeterogeneousAssembler<real_t>assembler (ppde,refine_bases,dirichlet::elimination);

            time.restart();
            gsIETIAssembler<real_t> ass(assembler);
            ass.setOptions(option.getGroup("IETI"));

            ass.init();
            time.stop();
            gsInfo<<"Time for preparing the bookkeeping: "<<time<<"\n";
	    gsIETIInfo info = ass.getInfo();

            nPatches(r,0) = patches.nPatches();
	    gsInfo<<"locDof: "<<refine_bases.size(0)<<"\n";
            globDof(r,0)=info.dofTotal;
            lagMult(r,0)=info.lagrangeMult;
            dofsP(r,0) = info.dofTotalP;

            if(!test3D)
                Hh(r,0) = math::sqrt( (real_t)(refine_bases.size(0)) )-1; //Hh domsize cancels
            else
                Hh(r,0) = math::pow( (real_t)(refine_bases.size(0)), 1.0/3)-1;
	    
	    if(nPrec == 1) {continue;}
	    
            time.restart();
            ass.assemble();
            time.stop();
            gsInfo<<"Time for assembling everything: "<<time<<"\n";

            gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(ass);
            solv->init();


            for(int pre = 1;pre<nPrec;pre++)
            {
                solVector.setZero(ass.systemSize(),ass.numberRhs());
                rhs = solv->getRhs();


                gsConjugateGradient<> PCG(solv);
                PCG.setCalcEigenvalues(true);
                PCG.setMaxIterations(200);
                //*
                if(pre == 0 && false)
                {
                    PCG.solve(rhs,solVector);

                    iter(r,pre+1)=PCG.iterations();

                    PCG.getEigenvalues(eigs);
                    cond(r,pre+1)=eigs.maxCoeff()/eigs.minCoeff();
                    continue;
                }
                //*/
                IETIPrecondScaling::strategy preEnum = static_cast<IETIPrecondScaling::strategy>(pre+1);
                ass.getOptions().scal = preEnum;
                gsScaledDirichletPrecond<real_t>::Ptr precDir = gsScaledDirichletPrecond<real_t>::make(ass);

                time.restart();
                PCG.setPreconditioner(precDir);
                PCG.solve(rhs,solVector);
                time.stop();
                gsInfo<<"Time for solving PCG total: "<<time<<"\n";
                iter(r,pre)=PCG.iterations();
                PCG.getEigenvalues(eigs);
                cond(r,pre)=eigs.maxCoeff()/eigs.minCoeff();
                //gsInfo<<eigs<<"\n\n";
                if(pre==nPrec-1 && false)
                {
                    solv->calculateSolution(solVector,solution);
                    sol = assembler.constructSolution(solution);
                    l2Error(r,0)= sol.distanceL2(g);
                    gsWriteParaview<>(sol, "IETI", 1000, true);
                    int unused = system("paraview IETI.pvd");
                    GISMO_UNUSED(unused);
                }

            }


        }


        gsInfo << "\n";
        printField("gDof",8);printField("lagM",6);printField("lDof",6); printField("nPi",6); printField("iF",12);printField("iPF-coeff",12);printField("iPF-stiff",12);printField("iPF-stiffM",12);printField("cF",12); printField("cPF-coeff",12);printField("cPF-stiff",12);printField("cPF-stiffM",12);
        gsInfo<<"\n";
        for(int i=0;i<max_patchRef;i++)
            printElement(globDof(i,0),lagMult(i,0),nPatches(i,0),dofsP(i,0),iter.row(i),cond.row(i));

        gsInfo<<"\n";
        gsInfo <<l2Error<<"\n"<<"\n";

    }


    return result;
}

