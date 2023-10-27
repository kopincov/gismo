// This file is not a part of the library!

// NOT IMPLEMENTED YET

// Implementation of SSP RK method for time dependent convdiff problems, AFC stabilization used

// Author: A.Jaeschke, parts copied from ConvDiff.cpp
#define _USE_MATH_DEFINES

#include <iostream>

#include <gismo.h>
#include <gsAssembler/gsCDRAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <gsAssembler/gsNormL2.h>
#include <iomanip>


using namespace gismo;

int main(int argc, char *argv[])
{
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 0;
    index_t basisDegree = 0;
    bool plot       = false;
    // Flag for SUPG-Stabilization
    //bool Flag_Stabilization = 0;
    //Flag_Stabilization = 1;

    gsCmdLine cmd("Testing a multipatch problem.");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("p","degree",
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    cmd.addInt("e","degreeElevation",
               "Number of degree elevation steps to perform on the Geometry's basis before solving",
               numElevate);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //gsFunctionExpr<real_t>  init("sin(x)*sin(y)", 2);
    //gsFunctionExpr<real_t>  init("if((sqrt((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5))/0.15<=1),0.25*(1+cos(pi*min(sqrt((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5))/0.15,1))),0) + if((sqrt((x-0.5)*(x-0.5)+(y-0.25)*(y-0.25))/0.15<=1),1-(sqrt((x-0.5)*(x-0.5)+(y-0.25)*(y-0.25))/0.15),0)+if((sqrt((x-0.5)*(x-0.5)+(y-0.75)*(y-0.75))/0.15<=1) and((abs(x-0.5)>=0.025) or (y>=0.85)),1,0)", 2);
    gsFunctionExpr<real_t>  init("if((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)<=0.35*0.35,1.0,0.0)", 2);
    gsFunctionExpr<real_t>  coeff_cM("1.0", 2);
    gsFunctionExpr<real_t>  coeff_AM("0.0","0","0","0.0", 2);
    gsFunctionExpr<real_t>  coeff_bM("0","0", 2);
    gsMultiPatch<real_t> mp( *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0));

    gsMultiBasis<> bases(mp);

    if (basisDegree)
        bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate);

    gsSparseSolver<>::LU solver;
    bases.uniformRefine();
    bases.uniformRefine();
    gsMatrix<> testsL2(numHref+1,7);
    testsL2.setZero();

    for (int i=2; i<=numHref ; i++)
    {


        gsBoundaryConditions<> BCs;
        gsCDRAssembler<real_t> galerkinM(mp,bases,BCs,init,coeff_AM,coeff_bM,coeff_cM,
                                         dirichlet::none, iFace::glue, stabilizerCDR::none);
        galerkinM.assemble();
        gsSparseMatrix<real_t> M (galerkinM.matrix());
        for (int e=0; e<M.outerSize(); ++e)
            for (gsSparseMatrix<real_t>::InnerIterator it(M,e); it; ++it)
            {
                int k = it.row();
                int l = it.col();
                if((M.at(k,l) != 0)&& k!=l)
                {
                    M.coeffRef(k,k) = M.at(k,k) + M.at(k,l);
                    M.coeffRef(k,l) = 0.0;
                }
            }


        //Constrained L2 data projection
        gsMatrix<real_t> u_l (galerkinM.rhs());
        gsMatrix<real_t> u_h (galerkinM.rhs());
        solver.compute(galerkinM.matrix());
        u_h = solver.solve(galerkinM.rhs());
        solver.compute(M);
        u_l = solver.solve(galerkinM.rhs());
        gsMatrix<real_t> u_init (u_l);
        //----------------------------------------------------------------------------------------

        {
            gsMatrix<real_t> Pp (galerkinM.rhs());
            gsMatrix<real_t> Pm (galerkinM.rhs());
            gsMatrix<real_t> Qp (galerkinM.rhs());
            gsMatrix<real_t> Qm (galerkinM.rhs());
            gsMatrix<real_t> ff (galerkinM.rhs());
            for (int k = 0; k < M.outerSize(); k++)
            {
                Pp.coeffRef(k) = 0.0;
                Pm.coeffRef(k) = 0.0;
                Qp.coeffRef(k) = 0.0;
                Qm.coeffRef(k) = 0.0;
                ff.coeffRef(k) = 0.0;
            }
            gsSparseMatrix<real_t> f (galerkinM.matrix());
            for (int e=0; e<M.outerSize(); ++e)
                for (gsSparseMatrix<real_t>::InnerIterator it(f,e); it; ++it)
                {
                    int k = it.row();
                    int l = it.col();
                    if ((u_h(k)-u_h(l))*(u_l(k)-u_l(l)) > 0)
                        f.coeffRef(k,l) = galerkinM.matrix().at(k,l)*(u_h.at(k)-u_h.at(l));
                    else
                        f.coeffRef(k,l) = 0.0;
                    if (l!= k)
                    {
                        Pp.coeffRef(k) += math::max<real_t>(0.0,f.at(k,l));
                        Pm.coeffRef(k) += math::min<real_t>(0.0,f.at(k,l));
                        Qp.coeffRef(k)  = math::max<real_t>(Qp.coeffRef(k),u_l(l)-u_l(k));
                        Qm.coeffRef(k)  = math::min<real_t>(Qm.coeffRef(k),u_l(l)-u_l(k));
                    }
                }
            gsMatrix<real_t> Rp (galerkinM.rhs());
            gsMatrix<real_t> Rm (galerkinM.rhs());
            for (int k = 0; k < M.outerSize(); k++)
            {
                if (Pp.at(k)==0)
                    Rp.coeffRef(k) = 1.0;
                else
                    Rp.coeffRef(k) = math::min<real_t>(1.0,M.at(k,k)*Qp.at(k)/(Pp.at(k)));
                if (Pm.at(k)==0)
                    Rm.coeffRef(k) = 1.0;
                else
                    Rm.coeffRef(k) = math::min<real_t>(1.0,M.at(k,k)*Qm.at(k)/(Pm.at(k)));
            }
            for (int e=0; e<f.outerSize(); ++e)
                for (gsSparseMatrix<real_t>::InnerIterator it(f,e); it; ++it)
                {
                    int k = it.row();
                    int l = it.col();
                    if (f.at(k,l) > 0.0)
                    {
                        ff.coeffRef(k) += math::min(Rp.at(k),Rm.at(l))*f.at(k,l);
                    }
                    else
                    {
                        ff.coeffRef(k) += math::min(Rm.at(k),Rp.at(l))*f.at(k,l);
                    }
                }
            for (int k=0; k<M.outerSize();k++)
                ff.coeffRef(k) = ff.at(k)/M.at(k,k);
            u_init = u_l + ff;

            testsL2(i,0) =bases.totalSize();
            gsField<> solH = galerkinM.constructSolution(u_h);
            gsField<> solL = galerkinM.constructSolution(u_l);
            gsField<> solAFC = galerkinM.constructSolution(u_init);

            if (plot && i== numHref)
            {
                gsInfo<<"Plotting in Paraview...\n";
                gsWriteParaview<>( solH, "aaaHIGH", 100000);
                gsWriteParaview<>( solL, "aaaLOW", 100000);
                gsWriteParaview<>( solAFC, "aaaAFCinit", 100000);
                //result = system("paraview poisson_problem.pvd &");

            }
            gsNormL2<real_t>* L2Norm;
            L2Norm=new gsNormL2<real_t>(solH,init);
            //L2Norm->compute();
            testsL2(i,1)= L2Norm->compute() ;
            L2Norm=new gsNormL2<real_t>(solL,init);
            //L2Norm->compute();
            testsL2(i,3)= L2Norm->compute() ;
            L2Norm=new gsNormL2<real_t>(solAFC,init);
            //L2Norm->compute();
            testsL2(i,5)= L2Norm->compute() ;
            delete L2Norm;
            if (i > 2)
            {
                testsL2(i,2)= testsL2(i-1,1) / testsL2(i,1);
                testsL2(i,4)= testsL2(i-1,3) / testsL2(i,3);
                testsL2(i,6)= testsL2(i-1,5) / testsL2(i,5);
            }
        }
        bases.uniformRefine();
    }
    for(int i = 1; i<= numHref; ++i)
    {   // Compute convergence rates
        testsL2(i,2)= math::log(testsL2(i,2))/math::log(2.0);
        testsL2(i,4)= math::log(testsL2(i,4))/math::log(2.0);
        testsL2(i,6)= math::log(testsL2(i,6))/math::log(2.0);
    }

    gsInfo << "Summary:\n\n";
    gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        High          |        Low        |   Constrained     \n";
    gsInfo << "    Dofs   |  L2 error  | conv. rate   \n" << testsL2  << "\n";


    return 0;
}
