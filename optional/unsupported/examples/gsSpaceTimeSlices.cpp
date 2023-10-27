/**  gsSpaceTimeSlices.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):     C. Hofer
    Created on:  2017-27-07

    DESCRIBE THE PURPOSE OF THE EXAMPLE HERE

    WHAT IS HAPPENING

    EXPECTED RESULT
*/

#include <iostream>
#include <gismo.h>

#include <gsAssembler/gsSpaceTimesliceAssembler.h>
#include <gsPde/gsSpaceTimePoissonPde.h>
#include <gsAssembler/gsSpaceTimeNorm.h>
#include <gsIETI/gsIETIAssembler.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIScaledDirichlet.h>


#include <gsSolver/gsConjugateGradient.h>

using namespace gismo;

void printNonZeros(gsSparseMatrix<real_t>& mat, std::string filename)
{
    std::fstream out(filename.c_str(),  std::ofstream::out | std::ofstream::trunc);
    for (int k=0; k<mat.outerSize(); ++k)
      for (gsSparseMatrix<real_t>::InnerIterator it(mat,k); it; ++it)
      {
          if(math::abs(it.value())>1.e-10)
               out<<util::to_string(it.col())+" "+util::to_string(-it.row())+" "+util::to_string(it.value())+"\n";
      }
    out.close();
}

int doExperiment(const gsSpaceTimePoissonPde<real_t>& pde,gsMultiBasis<real_t> refine_bases, int numRefine,const gsOptionList& opt, gsFunction<real_t>& exact)
{
    int result = 0;
    // h-refine each basis (4, one for each patch)
    for (int i = 0; i < numRefine; ++i)
        //refine_bases.uniformRefine();
        refine_bases.uniformRefineComponent(refine_bases.dim()-1);

    gsSpaceTimesliceAssembler<real_t> assembler(pde,refine_bases,opt.getReal("theta"), dirichlet::elimination);

 //   typedef gsSpaceTimesliceAssembler<real_t>::Permutation Permutation;
 //   Permutation perm = assembler.getTensorPermutation();
    // For comparism, the exact solution
    assembler.assemble(false,false);
    assembler.printTimings();

    gsSparseMatrix<> K = assembler.matrix();
    gsMatrix<> rhs = assembler.rhs();
    //K=K.twistedBy(perm);

  //  gsMatrix<real_t> mat1,mat2;
 //   assembler.computeTransferForTwoSlices(0,mat1,mat2);

    //gsInfo<<"K: "<<K.toDense()<<"\n";
  //  gsInfo<<"f_: "<<assembler.rhs()<<"\n";
  //  gsInfo<<"f: "<<perm*assembler.rhs()<<"\n";

    assembler.assemble(true,true);
    typedef gsSpaceTimesliceAssembler<real_t>::TPData TPData;
    TPData data = assembler.getTensorProductData();

        /*
    gsVector<index_t> sizes(data.Kt.size());
    gsVector<index_t> size1(1);
    size1[0] =1 ;
    for(size_t i=0; i<data.Kt.size();++i)
        sizes[i] = data.Kt[i]->rows()*data.Mx[i]->rows();


    gsSparseMatrix<real_t>::BlockView viewK = K.blockView(sizes,sizes);
    gsMatrix<real_t>::BlockView viewF = rhs.blockView(sizes,size1);
    for(size_t i=0; i<data.Kt.size();++i)
    {
        gsSparseMatrix<real_t> A,B, temp;
        temp=data.Kt[i]->kron(*data.Mx[i]);
        A=data.Mt[i]->kron(*data.Kx[i]);
        if(i>0)
            B=data.Nt[i]->kron(*data.MWx[i]);
        A+=temp;

        gsInfo<<"A on "<<i<<":\n"<<A.toDense()<<"\n\n";
        if(i>0)
       gsInfo<<"-B on "<<i<<":\n"<<-B.toDense()<<"\n\n";
       gsInfo<< "F on "<<i<<":\n"<< data.rhs[i]<<"\n\n"<<std::flush;

        gsInfo<<"A-K: "<<(A-viewK(i,i)).norm()<<" , f-F: "<<(viewF(i)-data.rhs[i]).norm()<<"\n\n";
        if(i==0)
        {
            gsInfo<<"f classic: "<<viewF(i).transpose()<<"\n\n"<<"f TP: "<<data.rhs[i].transpose()<<"\n";
            gsInfo<<(viewF(i)-data.rhs[i]).transpose()<<"\n\n";
        }
        if(i>0)
            gsInfo<<"B-K: "<<(-B-viewK(i,i-1)).norm()<<"\n\n";
    }
*/
    gsMatrix<real_t> Mt =data.Mt[0]->toDense(); gsMatrix<real_t> Kt =data.Kt[0]->toDense();

    gsMatrix<real_t>::EigenSolver esolvMt,esolvKt;
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd > esolvKtMt;
    esolvMt.compute(Mt); esolvKt.compute(Kt); esolvKtMt.compute(Kt,Mt);
    index_t i;
    gsInfo<<"p= "<<refine_bases.degree()<<"\t tau= "<<assembler.getTau(0)<<"\t theta= "<<assembler.getTheta()<<"\n";
    gsInfo<<"Min Eigenvalue of Mt: "<<esolvMt.eigenvalues().real().minCoeff(&i)<<"\n"; //"\t ("<<esolvMt.eigenvalues()(i)<<")\n";
    gsInfo<<"Min Eigenvalue of Kt: "<<esolvKt.eigenvalues().real().minCoeff(&i)<<"\n"; //"\t ("<<esolvKt.eigenvalues()(i)<<"\n";
    gsInfo<<"Min Eigenvalue of Mt^-1Kt: "<<esolvKtMt.eigenvalues().real().minCoeff(&i)<<"\n"; //"\t ("<<esolvKtMt.eigenvalues()(i)<<"\n\n";
    std::fstream out("EigenvaluesMtKt.dat",  std::ofstream::out | std::ofstream::app);
    out<<"("<<numRefine<<","<<assembler.getTimeBases().degree()<<","<<assembler.getTheta()<<")\t\t MtKt: "<<esolvKtMt.eigenvalues().real().minCoeff(&i)<<
      "\t\t Mt:"<<esolvMt.eigenvalues().real().minCoeff(&i)<<"\t\t Kt:"<<esolvKt.eigenvalues().real().minCoeff(&i)<<"\n";

/*
    gsMatrix<std::complex<real_t> > X =Mt*esolvKt.eigenvectors();
    gsInfo<<X<<"\n\n"<<esolvKtMt.eigenvectors()<<"\n\n";
    for(index_t i=0; i<X.cols();++i)
    {

        X.col(i)*=esolvKt.eigenvectors()(0,i)/X(0,i);
        gsInfo<<"difference: "<<(X.col(i)-esolvKt.eigenvectors().col(i)).norm()<<"\n";
        gsInfo<<esolvKt.eigenvectors().col(i).transpose()<<"\n";
        gsInfo<<X.col(i).transpose()<<"\n\n";

      //  gsMatrix<real_t> x = esolvKtMt.eigenvectors().col(i).real();
      //  gsMatrix<real_t> y = esolvKtMt.eigenvectors().col(i).imag();
        gsMatrix<real_t> x = gsMatrix<real_t>::Random(Mt.rows(),1);
        gsMatrix<real_t> y = gsMatrix<real_t>::Random(Mt.rows(),1);
        real_t norm = math::sqrt(x.squaredNorm()+y.squaredNorm());
      //  x/=norm; y/=norm;


        real_t a = (x.transpose()*Kt*x+ y.transpose()*Kt*y)(0,0);
        real_t b = (x.transpose()*Mt*x+ y.transpose()*Mt*y)(0,0);
        real_t c = (x.transpose()*Mt*y- y.transpose()*Mt*x)(0,0);
        real_t d = (x.transpose()*Kt*y- y.transpose()*Kt*x)(0,0);

        real_t a = (x.transpose()*Kt*x)(0,0);
        real_t b = (x.transpose()*Mt*x)(0,0);
        real_t c = (x.transpose()*Mt*y)(0,0);
        real_t d = (x.transpose()*Kt*y)(0,0);


        gsInfo<<"a: "<<a<<"\t b:"<<b<<"\t c:"<<c<<"\t d:"<<d<<"\n";
        gsInfo<<"c*d:"<<c*d<<"\t sqrt(-cd):"<<math::sqrt(-c*d)<<"\n";
        gsInfo<<"realPart: "<<(a*b)/(b*b+c*c)<<" + "<<c*d/(b*b+c*c)<<"= "<<a*b/(b*b+c*c)+c*d/(b*b+c*c)<<" compared to "<<esolvKtMt.eigenvalues()(i)<<"\n\n";
      //  d = ((x.transpose()*Kt*y).array().abs()+ (y.transpose()*Kt*x).array().abs())(0,0);
      //  gsInfo<<"\t realPart: "<<a*b/(b*b+c*c)<<" + "<<-assembler.getTau(0)*assembler.getTheta()*d*d/(b*b+c*c)<<"= "<<a*b/(b*b+c*c)-assembler.getTau(0)*assembler.getTheta()*d*d/(b*b+c*c)<<" compared to "<<esolvKtMt.eigenvalues()(i)<<"\n";
    }
    */
   // gsInfo<<"Mt:\n "<<Mt<<"\n\n"<<"Kt:\n "<<Kt<<"\n\n";


    gsInfo<<"Mt: "<<esolvMt.eigenvalues()<<"\n\n";
    gsInfo<<"Kt: "<<esolvKt.eigenvalues()<<"\n\n";
    gsInfo<<"Mt^{-1}Kt: "<<esolvKtMt.eigenvalues()<<"\n\n";
    return 0;
    gsIETIAssembler<real_t> IETI_ass(assembler);
    IETI_ass.setOptions(opt.getGroup("IETI"));
    IETI_ass.init();
    IETI_ass.assemble();

    gsIETISolver<real_t>::Ptr IETI_solv = gsIETISolver<real_t>::make(IETI_ass);
    IETI_solv->init();

    gsScaledDirichletPrecond<real_t>::Ptr prec = gsScaledDirichletPrecond<real_t>::make(IETI_ass);

    gsMatrix<real_t> x(IETI_solv->rows(),1);
    gsGMRes<> gmres(IETI_solv,prec);
    gmres.setMaxIterations(100);
    gmres.setTolerance(1.e-8);
    gmres.solve(IETI_solv->getRhs(),x);
    gsInfo<<"Converged after "<<gmres.iterations()<< " iterations with residuum "<<gmres.error()<<".\n"<<std::flush;
    gsMatrix<real_t> vecIETI;
    IETI_solv->calculateSolution(x,vecIETI);

#if defined(GISMO_WITH_PARDISO)
    gsSparseSolver<>::PardisoLU solver;
    //gsSparseSolver<>::PardisoLU::ParameterType& param =  solver.pardisoParameterArray();

   // param[10] = 0;
  //  param[12] = 0;
#elif defined(GISMO_WITH_SUPERLU)
    gsSparseSolver<>::SuperLU solver;
#else
    gsSparseSolver<>::LU solver;
#endif

    solver.compute(K);
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eiK(K.toDense());
    //gsInfo<<"eigenvalues K: "<<eiK.eigenvalues().transpose()<<"\n";
    gsMatrix<> solV = solver.solve( rhs );
   // gsInfo<<"solution:\n"<<solV<<"\n";

    gsField<real_t> sol = assembler.constructSolution(solV);
    gsField<real_t> solIETI = assembler.constructSolution(vecIETI);
    gsInfo<<"Dofs: "<<assembler.numDofs()<<"\n";
    //   gsFunctionExpr<> gder("pi*cos(pi*(x))*sin(pi*(y))","pi*sin(pi*x)*cos(pi*y)",2);
    gsSpaceTimeSliceNorm<real_t> norm(sol,exact,assembler.getdGInterfaces(), assembler.getInitialBoundary(),assembler.getTopBoundary());
    norm.setTheta(assembler.getTheta());
    gsSpaceTimeSliceNorm<real_t> normIETI(solIETI,exact,assembler.getdGInterfaces(), assembler.getInitialBoundary(),assembler.getTopBoundary());
    normIETI.setTheta(assembler.getTheta());
    real_t l2error = sol.distanceL2(exact);
    gsInfo<<"L2 error ||uh - g||: "<<l2error<<"\n";
    gsInfo<<"dG error ||uh - g||: "<<norm.compute()<<"\n";
    gsInfo<<"L2 error ||uhIETI - g||: "<<solIETI.distanceL2(exact)<<"\n";
    gsInfo<<"dG error ||uhIETI - g||: "<<normIETI.compute()<<"\n";

    // Plot solution in paraview

    if (opt.getSwitch("p"))
    {
        gsWriteParaview<>( sol, "ST_Sol", 1000,true);
    //    gsWriteParaview<>( solIETI, "ST_IETI", 1000,true);

        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        const gsField<real_t> exactF( pde.domain(), exact, false );
        const gsField<real_t> rhsF( pde.domain(), *pde.rhs(), false );
      //  gsWriteParaview<real_t>( exactF, "ST_Ex", 1000);
      //  gsWriteParaview<real_t>( rhs, "ST_rhs", 1000);

        // Run paraview
        result = system("paraview ST_Sol.pvd  &");


        gsField<real_t> field = gsFieldCreator<real_t>::absError(solIETI,exact);
        gsWriteParaview<>( field, "ST_diff"+util::to_string(0), 1000);

        // Run paraview
        result += system(("paraview ST_diff"+util::to_string(0)+".pvd  &").c_str());

    }
    return result;
}

int doConvergenceTest(const gsSpaceTimePoissonPde<real_t>& pde, gsMultiBasis<real_t> refine_bases, int numRefine, const gsOptionList& opt, gsFunction<real_t>& exact)
{
    int result = 0;
    gsStopwatch time;
    gsMatrix<real_t> error(numRefine,1);
    gsMatrix<real_t> factor(numRefine,1);
    error.setZero();
    factor.setZero();

    gsInfo<<"chosen theta:"<<opt.getReal("theta")<<"\n";
    for (int i = 0; i < numRefine; ++i)
    {
        refine_bases.uniformRefine();

        gsSpaceTimesliceAssembler<real_t> assembler(pde,refine_bases,opt.getReal("theta"),dirichlet::elimination);
        time.restart();
        assembler.assemble();
        assembler.printTimings();
        gsInfo<<"Time for assembling matrix: "<<time.stop()<<" \n";

        const gsSparseMatrix<>& K = assembler.matrix();

        #if defined(GISMO_WITH_PARDISO)
        gsInfo<<"Using Pardiso!\n";
            gsSparseSolver<>::PardisoLU solver;
            //gsSparseSolver<>::PardisoLU::ParameterType& param =  solver.pardisoParameterArray();
           // gsInfo<<"parameters: "<<param<<"\n";

         //   param[10] = 0;
         //   param[12] = 0;
        #elif defined(GISMO_WITH_SUPERLU)
        gsInfo<<"Using SuperLU!\n";
            gsSparseSolver<>::SuperLU solver;
        #else
        gsInfo<<"Using SparseLU!\n";
            gsSparseSolver<>::LU solver;
        #endif


        time.restart();
        solver.compute(K);
        gsMatrix<> solV = solver.solve( assembler.rhs() );
        gsInfo<<"Time for  solving system: "<<time.stop()<<" \n";

        gsField<> sol = assembler.constructSolution(solV);
        //   gsFunctionExpr<> gder("pi*cos(pi*(x))*sin(pi*(y))","pi*sin(pi*x)*cos(pi*y)",2);
        gsSpaceTimeSliceNorm<real_t> norm(sol,exact,assembler.getdGInterfaces(), assembler.getInitialBoundary(),assembler.getTopBoundary());
        norm.setTheta(assembler.getTheta());
        // error(i,0) = sol.distanceH1(exact);
        time.restart();
        error(i,0) = norm.compute();
        gsInfo<<"Time to calculate norm: "<<time.stop()<<" \n";

        gsInfo<<"Dofs: "<<assembler.numDofs()<<"\n";
        gsInfo<<"Finished:"<<i<<"\n";

        if (opt.getSwitch("p"))
        {
            if(i==0)
            {
                gsWriteParaview<>( sol, "ST_Sol", 1000,true);

                // Write approximate and exact solution to paraview files
                gsInfo<<"Plotting in Paraview...\n";
                const gsField<> exactF( pde.domain(), exact, false );
              //  const gsField<> rhs( pde.domain(), *pde.rhs(), false );
                gsWriteParaview<>( exactF, "ST_Ex", 1000);
              //  gsWriteParaview<>( rhs, "ST_rhs", 1000);

                // Run paraview
                result = system("paraview ST_Sol.pvd  &");
            }

            gsField<real_t> field = gsFieldCreator<real_t>::absError(sol,exact);
            gsWriteParaview<>( field, "ST_diff"+util::to_string(i), 1000);

            // Run paraview
            result += system(("paraview ST_diff"+util::to_string(i)+".pvd  &").c_str());

        }

    }
    for(int i=1;i<numRefine;i++)
    {
        factor(i,0) = math::log(error(i-1,0)/error(i,0))/math::log(2.0);
    }
    gsInfo<<"Error:\n"<<error<<"\n";
    gsInfo<<"Convergence factor:\n"<<factor<<"\n";
    return result;
}



//Important Note: The boundary condition has to be called with
// bcInfo.addCondition(***, *** , *** , g);
// in order that the parallelization works correctly! So no &g !!!
// Otherwise, you have a data race!


int main (int argc, char** args)
{
    index_t testcase =1;
    bool test3D = false;
    bool plot = false;
    index_t test=0;
    real_t theta = 1;
    real_t Tlen = 1;

    std::string m = "A";

    // Number for h-refinement of the computational (trail/test) basis.
    index_t numRefine  = 1;
    // Number for p-refinement of the computational (trail/test) basis.
    index_t numElevate = 0;

    gsCmdLine cmd("Hi, I will test IETIdG");
    cmd.addSwitch("p","plot","Plot solution",plot);
    cmd.addInt("r","refine","Number of refinements",numRefine);
    cmd.addInt("c","case","Choosen testcase",testcase);
    cmd.addInt("e","elevate","maximal elevation",numElevate);
    cmd.addInt("","test","(0) experiment, (1) convergence test", test);
    cmd.addString("s", "IETI.Strategy", "Choosen strategy",m);
    cmd.addReal("","theta","specify the parameter theta",theta);
    cmd.addReal("","Tlen","length of one timeslice", Tlen);

    try { cmd.getValues(argc,args); } catch (int rv) { return rv; }

    gsOptionList option = cmd.getOptionList();
    gsOptionList optIETI = gsIETIAssembler<real_t>::defaultOptions();
    optIETI = optIETI.wrapIntoGroup("IETI");
     option.update(optIETI,gsOptionList::addIfUnknown);
     option.setSwitch("IETI.NonSymmetric",true);
    //test3D = true;
    //nitsche = true;

    index_t result = 0;

    int dim;
    if(testcase==6) test3D=true;
    (test3D ||testcase == 8 || testcase == 7) ? dim =3:dim =2;
    if(testcase == 10)
        dim = 4;
    if(dim == 4)
        plot  = false;
    gsInfo<<"Dimension of Domain: "<<dim<<"\n";

    //////// Right-hand side and analytical solution ////////
    // Define source function
    gsFunctionExpr<> f,f_x,f_t,g;

    if(dim==2)
    {

        f= gsFunctionExpr<>("pi*sin(pi*x)*(0.5*cos(0.5*pi*y+pi/2)+pi*sin(0.5*pi*y+pi/2)) ",dim);
        f_x= gsFunctionExpr<>("pi*sin(pi*x)",dim-1);
        f_t=gsFunctionExpr<>("0.5*cos(0.5*pi*x+pi/2)+pi*sin(0.5*pi*x+pi/2)",1);
        g=gsFunctionExpr<>("sin(pi*x)*sin(0.5*pi*y+pi/2)",dim);


    //          f= gsFunctionExpr<>("0",dim);
   //           f_x= gsFunctionExpr<>("0",dim-1);
   //         f_t=gsFunctionExpr<>("0",1);
   //         g=gsFunctionExpr<>("1",dim);



    }
    else if(dim==3)
    {
        f= gsFunctionExpr<>("pi*sin(pi*x)*sin(pi*y)*(cos(pi*z)+2*pi*sin(pi*z)) ",dim);
        f_x= gsFunctionExpr<>("pi*sin(pi*x)*sin(pi*y)",dim-1);
        f_t=gsFunctionExpr<>("cos(pi*x)+2*pi*sin(pi*x)",1);
        g=gsFunctionExpr<>("sin(pi*x)*sin(pi*y)*sin(pi*z)",dim);

      //  f= gsFunctionExpr<>("if((x-0.7)^2+(y-4.5)^2<0.04,5,0)",dim);
     //   f_x= gsFunctionExpr<>("if((x-0.4)^2+(y-4.6)^2<0.04,5,0)",dim-1);
    //      f= gsFunctionExpr<>("0",dim);
    //      f_x= gsFunctionExpr<>("0",dim-1);
    //    f_t=gsFunctionExpr<>("0",1);
    //    g=gsFunctionExpr<>("10",dim);
    }
    else if(dim==4)
    {
        f= gsFunctionExpr<>("pi*sin(pi*x)*sin(pi*y)*sin(pi*z)*(cos(pi*w)+3*pi*sin(pi*w)) ",dim);
        f_x= gsFunctionExpr<>("pi*sin(pi*x)*sin(pi*y)*sin(pi*z)",dim-1);
        f_t=gsFunctionExpr<>("cos(pi*x)+3*pi*sin(pi*x)",1);
        g=gsFunctionExpr<>("sin(pi*x)*sin(pi*y)*sin(pi*z)*sin(pi*w)",dim);
    }

    //  gsFunctionExpr<> dn_g("pi*sin(pi*x)",dim-1);//Neumann at x=1
    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n" << "\n";

    /////////////////// Setup geometry ///////////////////
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;


    switch(testcase){
    case 0:
        patches = gsMultiPatch<>(*gsNurbsCreator<>::BSplineRectangle(0,0,1,1));
        patches.patch(0).degreeReduce();
        break;
    case 1:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,1, 1);
        break;
    case 2:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,2, 0.5);
        break;
    case 3:
        patches = gsNurbsCreator<>::BSplineSquareGrid(3,1, 1);
        break;
    case 4:
        patches = gsNurbsCreator<>::BSplineSquareGrid(1,4, 1);
        break;
    case 5:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,4, 0.5);
        break;
    case 6:
    {
        patches = gsMultiPatch<>(*gsNurbsCreator<>::BSplineRectangle(0,0,1,1));
        break;
    }
    case 7:
    {
        patches = gsNurbsCreator<>::BSplineCubeGrid(4,4,1,0.2,0.5,0.5,0.5);
        break;
    }
    case 8:
    {
        gsFileData<> fileData("yeti3d_u21p.xml");
        if (fileData.has< gsMultiPatch<> >())
        {
            fileData.getFirst(patches);
            size_t nPatches= patches.nPatches();
            int max_slices = 1;
            for(int slice = 1; slice<max_slices;++slice)
            {
                for(size_t np = (slice-1)*nPatches; np< slice*nPatches;++np)
                {
                    gsGeometry<real_t>::uPtr newPatch = patches.patch(np).clone();
                    gsMatrix<real_t>& coefs = newPatch->coefs();
                    real_t diff = coefs.col(2).maxCoeff() - coefs.col(2).minCoeff();
                    for(int j=0; j<coefs.rows();++j)
                        coefs.col(2)(j)+=diff;
                    patches.addPatch(give(newPatch));
                }
            }
            for(size_t np =0; np<patches.nPatches();++np)
                patches.patch(np).degreeElevate(1,2);
            patches.computeTopology();
        }

        else
            return 1;
        break;
    }
    case 10:
    {
        patches = gsNurbsCreator<>::BSplineCubeGrid(2,2,2,1);
        gsMultiPatch<real_t> patches4D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<3,real_t>* tb = static_cast<gsTensorBSpline<3,real_t>* >(&patches.patch(i));
            patches4D.addPatch(gsNurbsCreator<>::lift4D(*tb));
        }

            int nPatches= patches4D.nPatches();
            int max_slices = 2;
            for(int slice = 1; slice<max_slices;++slice)
            {
                for(int np = (slice-1)*nPatches; np< slice*nPatches;++np)
                {
                    gsGeometry<real_t>::uPtr newPatch = patches4D.patch(np).clone();
                    gsMatrix<real_t>& coefs = newPatch->coefs();
                    real_t diff = coefs.col(3).maxCoeff() - coefs.col(3).minCoeff();
                    for(int j=0; j<coefs.rows();++j)
                        coefs.col(3)(j)+=diff;
                    patches4D.addPatch(give(newPatch));
                }
            }
            patches4D.computeTopology();
            patches = give(patches4D);

        break;
    }
    default:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,1, 1);
        break;
    }


#ifdef _OPENMP
    Eigen::initParallel();
    // To get the number of threads from the environment
#endif


    /////////////////// Setup boundary conditions ///////////////////
    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    gsFunctionExpr<> hEast, hSouth,hNorth, hWest;

    hWest = gsFunctionExpr<>("-pi*cos(pi*0.4)*sin(2*pi*(y+0.3))-1",dim);

    switch(testcase){
    case 0:

        bcInfo.addCondition(0, boundary::east,  condition_type::dirichlet, g);
        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);

        bcInfo.add(0, boundary::south, "Initial", g);

        break;
    case 1:
        bcInfo.addCondition(1, boundary::east,  condition_type::dirichlet, g);
        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);

        bcInfo.add(0, boundary::south, "Initial", g);
        bcInfo.add(1, boundary::south, "Initial", g);

        break;
    case 2:
        bcInfo.addCondition(2, boundary::east,  condition_type::dirichlet, g);
        bcInfo.addCondition(3, boundary::east,  condition_type::dirichlet, g);

        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
        bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);

        bcInfo.add(0, boundary::south, "Initial", g);
        bcInfo.add(2, boundary::south, "Initial", g);

        break;
    case 3:
        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, g);
        bcInfo.addCondition(2, boundary::east,  condition_type::dirichlet, g);

        bcInfo.add(0, boundary::south, "Initial", g);
        bcInfo.add(1, boundary::south, "Initial", g);
        bcInfo.add(2, boundary::south, "Initial", g);

        break;
    case 4:

        for(size_t i=0; i<patches.nPatches();++i)
        {
            bcInfo.addCondition(i, boundary::east,  condition_type::dirichlet, g);
            bcInfo.addCondition(i, boundary::west,  condition_type::dirichlet, g);
        }

        bcInfo.add(0, boundary::south, "Initial", g);

        break;
    case 5:

        bcInfo.addCondition(12, boundary::east,  condition_type::dirichlet, g);
        bcInfo.addCondition(13, boundary::east,  condition_type::dirichlet, g);
        bcInfo.addCondition(14, boundary::east,  condition_type::dirichlet, g);
        bcInfo.addCondition(15, boundary::east,  condition_type::dirichlet, g);
        //*/

        /*
        bcInfo.addCondition(12, boundary::east,  condition_type::neumann, dn_g);
        bcInfo.addCondition(13, boundary::east,  condition_type::neumann, dn_g);
        bcInfo.addCondition(14, boundary::east,  condition_type::neumann, dn_g);
        bcInfo.addCondition(15, boundary::east,  condition_type::neumann, dn_g);
        //*/

        bcInfo.add(0, boundary::south, "Initial", g);
        bcInfo.add(4, boundary::south, "Initial", g);
        bcInfo.add(8, boundary::south, "Initial", g);
        bcInfo.add(12, boundary::south, "Initial", g);

        bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
        bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);
        bcInfo.addCondition(2, boundary::west,  condition_type::dirichlet, g);
        bcInfo.addCondition(3, boundary::west,  condition_type::dirichlet, g);
        break;

    case 6:

        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
            bcInfo.addCondition(*it, condition_type::dirichlet, g);
        break;
    case 7:
    case 8:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
        {
            if(it->side() == boundary::front)
            {
                bcInfo.add(it->patch, boundary::front, "Initial", g);
               // bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, f);
            }
            else if(it->side() != boundary::back)
                bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, g);
        }
        break;
    case 10:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
        {
            if(it->side() == boundary::stime)
            {
                bcInfo.add(it->patch, it->side(), "Initial", g);
             //   bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, g);
            }
            else if(it->side() != boundary::etime)
                bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, g);
        }
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
        gsMultiPatch<real_t> patches3D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&patches.patch(i));
            patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb));
            bcInfo.add(i, boundary::front, "Initial", g);
            if(testcase == 6)
                patches3D.patch(i).degreeElevate(1,2);

        }
        patches3D.computeTopology();

        patches = give(patches3D);

    }
    
    if (plot)
    {
        gsWriteParaview<>( patches, "IETI_patch", 1000, true);
        result = system("paraview IETI_patch.pvd &");
    }

    real_t len = patches.patch(0).coefs().col(dim-1).maxCoeff() - patches.patch(0).coefs().col(dim-1).minCoeff();
    for(size_t np = 0; np<patches.nPatches();++np)
        patches.patch(np).scale(Tlen/len,dim-1);

    ////////////////////// Refinement h and p //////////////////////
    // Refinement

    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );


    // Find maximum degree with respect to all the variables
   // refine_bases.setDegree(refine_bases.maxCwiseDegree());
    refine_bases.degreeIncrease(numElevate);
    gsInfo<<"Degree set to: "<<refine_bases.maxCwiseDegree() <<" -- "<<refine_bases.minCwiseDegree()<<"\n";

    if(true && (testcase == 4 || testcase == 5 || testcase == 7))
    {
         for(size_t i = 0;i<(patches).nPatches();++i)
         {
             if(i%2 == 0)
                 refine_bases[i].uniformRefine();
         }
    }

    ///////////////////////ASSEMBLE////////////////////////////////////

    gsInfo<<bcInfo<<"\n";


    gsPiecewiseFunction<real_t>  alpha;
    gsFunctionExpr<> a1("1.e+0",dim-1);
    gsFunctionExpr<> a2("1.e+0",dim-1);
    for(size_t np = 0; np<patches.nPatches();np++)
        if(np%2 == 0)
            alpha.addPiece(a1);
        else
            alpha.addPiece(a2);


    gsSpaceTimePoissonPde<real_t> pde(patches,bcInfo,f,f_x,f_t,alpha,&g);

    switch(test)
    {
    case 0:
        result += doExperiment(pde, refine_bases, numRefine,option,g);
        break;
    case 1:
        result += doConvergenceTest(pde, refine_bases, numRefine,option,g);
        break;
    }





    //return  1 for failures and 0 for success
    return result;
}


