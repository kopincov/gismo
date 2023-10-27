// This file is not a part of the library!

// NOT IMPLEMENTED YET

// Implementation of SSP RK method for time dependent convdiff problems, AFC stabilization used

// Author: A.Jaeschke, parts copied from ConvDiff.cpp


#include <iostream>

#include <gismo.h>
#include <gsAssembler/gsCDRAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <iomanip>
#include <gismo.h>

#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsSeminormH1.h>

using namespace gismo;

gsMatrix<real_t> AFC_Correction(gsSparseMatrix<real_t>* L, gsMatrix <real_t> sol, int bases1d, int bases2d,
                                gsCDRAssembler<real_t> *galerkinK, gsCDRAssembler<real_t> *galerkinM, gsSparseMatrix<real_t>* Ml,
                                gsSparseSolver<>::LU * solver, real_t tau,
                                gsMultiPatch<real_t>* mp, gsMultiBasis<>* bases, gsFunctionExpr<real_t>  g);

int main(int argc, char *argv[])
{
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 0;
    index_t basisDegree = 0;
    bool plot       = false;
    index_t ntsteps = 300;
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
    cmd.addInt("n","nTimesteps", "Number of time steps", ntsteps);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //Time dependent problem definition:
    //Partial
    //gsFunctionExpr<real_t>  init("if((sqrt((x-0.5)*(x-0.5)+(y-0.75)*(y-0.75))/0.15<=1) and((abs(x-0.5)>=0.025) or (y>=0.85)),1,0)", 2);
    //gsFunctionExpr<real_t>  init("if((sqrt((x-0.5)*(x-0.5)+(y-0.25)*(y-0.25))/0.15<=1),1-(sqrt((x-0.5)*(x-0.5)+(y-0.25)*(y-0.25))/0.15),0)", 2);
    //gsFunctionExpr<real_t>  init("if((sqrt((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5))/0.15<=1),0.25*(1+cos(pi*min(sqrt((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5))/0.15,1))),0)", 2);
    //Total - solid body rotation
    gsFunctionExpr<real_t>  init("if((sqrt((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5))/0.15<=1),0.25*(1+cos(pi*min(sqrt((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5))/0.15,1))),0) + if((sqrt((x-0.5)*(x-0.5)+(y-0.25)*(y-0.25))/0.15<=1),1-(sqrt((x-0.5)*(x-0.5)+(y-0.25)*(y-0.25))/0.15),0)+if((sqrt((x-0.5)*(x-0.5)+(y-0.75)*(y-0.75))/0.15<=1) and((abs(x-0.5)>=0.025) or (y>=0.85)),1,0)", 2);
    //Simplified
    //gsFunctionExpr<real_t>  init("if((sqrt((x-0.6)*(x-0.6)+(y-0.6)*(y-0.6))/0.35<=1),1,0)", 2);
    //sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))/0.15


    gsFunctionExpr<real_t>  f("0.0", 2);
    gsFunctionExpr<real_t>  g("0.0", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
    gsFunctionExpr<real_t>  coeff_cM("1.0", 2);
    gsFunctionExpr<real_t>  coeff_AS("0.0","0","0","0.0", 2);
    gsFunctionExpr<real_t>  coeff_bS("0","0", 2);
    gsFunctionExpr<real_t>  coeff_AK("0","0","0","0", 2);

    gsFunctionExpr<real_t>  coeff_bK("0.5-y","x-0.5", 2);
    gsMultiPatch<real_t> mp( *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0) );


    gsMultiBasis<> bases(mp);

    if (basisDegree)
        bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate);

    gsSparseSolver<>::LU solver;
    for (int i=0; i<numHref ; i++)
        bases.uniformRefine();


    // Basic version, all parameters are NOT time dependent
    gsBoundaryConditions<> BCs;
    gsCDRAssembler<real_t> galerkinS(mp,bases,BCs,f,coeff_AS,coeff_bS,coeff_c,
                                     dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinS.assemble();
    gsCDRAssembler<real_t> galerkinK(mp,bases,BCs,f,coeff_AK,coeff_bK,coeff_c,
                                     dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinK.assemble();
    gsCDRAssembler<real_t> galerkinM(mp,bases,BCs,f,coeff_AK,coeff_bS,coeff_cM,
                                     dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinM.assemble();
    gsCDRAssembler<real_t> galerkinIN(mp,bases,BCs,init,coeff_AK,coeff_bS,coeff_cM,
                                      dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinIN.assemble();

    gsSparseMatrix<real_t> M (galerkinM.matrix());

    ///*
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
    //*/
    gsMatrix<real_t> rhs (galerkinK.rhs());
    gsSparseMatrix<real_t> K (galerkinK.matrix());
    ///*
    real_t delta;
    for (int e=0; e<K.outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
        {
            int k = it.row();
            int l = it.col();
            real_t temp1 = K.at(k,l);
            real_t temp2 = K.at(l,k);
            if((temp1 > 0 || temp2 > 0 )&& k!=l)
            {
                delta = (math::max(temp1, temp2));
                K.coeffRef(k,l) = temp1-delta;
                K.coeffRef(l,k) = temp2-delta;
                K.coeffRef(k,k) = K.at(k,k)+delta;
                K.coeffRef(l,l) = K.at(l,l)+delta;
            }
        }
    //*/
    K = K + galerkinS.matrix();
    K = -K;
    //Evaluation of the initial condition
    //Standard L2
    /*gsMatrix<real_t> anch;
      gsMatrix<real_t> fpts;
      gsGeometry<real_t> & geom = mp.patch(0);
      gsBasis<real_t> & basis = bases[0];
      anch = basis.anchors();
      fpts = init.eval(  geom.eval(anch) );
      gsGeometry<real_t> * geo = basis.interpolateAtAnchors(fpts);
      const gsMatrix<real_t> Vals =  geo->coefs();
      gsMatrix<real_t> u_init (Vals);
      delete geo;
    //*/
    //Constrained L2 data projection
    gsMatrix<real_t> u_l (galerkinK.rhs());
    gsMatrix<real_t> u_h (galerkinK.rhs());
    solver.compute(galerkinM.matrix());
    u_h = solver.solve(galerkinIN.rhs());
    solver.compute(M);
    u_l = solver.solve(galerkinIN.rhs());
    gsMatrix<real_t> u_init (u_l);
    //----------------------------------------------------------------------------------------

    {
        gsMatrix<real_t> Pp (galerkinK.rhs());
        gsMatrix<real_t> Pm (galerkinK.rhs());
        gsMatrix<real_t> Qp (galerkinK.rhs());
        gsMatrix<real_t> Qm (galerkinK.rhs());
        gsMatrix<real_t> ff (galerkinK.rhs());
        for (int k = 0; k < M.outerSize(); k++)
        {
            Pp.coeffRef(k) = 0.0;
            Pm.coeffRef(k) = 0.0;
            Qp.coeffRef(k) = 0.0;
            Qm.coeffRef(k) = 0.0;
            ff.coeffRef(k) = 0.0;
        }
        gsSparseMatrix<real_t> fMatrix (galerkinM.matrix());
        for (int e=0; e<K.outerSize(); ++e)
            for (gsSparseMatrix<real_t>::InnerIterator it(fMatrix,e); it; ++it)
            {
                int k = it.row();
                int l = it.col();
                if ((u_h(k)-u_h(l))*(u_l(k)-u_l(l)) > 0)
                    fMatrix.coeffRef(k,l) = galerkinM.matrix().at(k,l)*(u_h.at(k)-u_h.at(l));
                else
                    fMatrix.coeffRef(k,l) = 0.0;
                if (l!= k)
                {
                    Pp.coeffRef(k) += math::max<real_t>(0.0,fMatrix.at(k,l));
                    Pm.coeffRef(k) += math::min<real_t>(0.0,fMatrix.at(k,l));
                    Qp.coeffRef(k)  = math::max<real_t>(Qp.coeffRef(k),u_l(l)-u_l(k));
                    Qm.coeffRef(k)  = math::min<real_t>(Qm.coeffRef(k),u_l(l)-u_l(k));
                }
            }
        gsMatrix<real_t> Rp (galerkinK.rhs());
        gsMatrix<real_t> Rm (galerkinK.rhs());
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
        for (int e=0; e<fMatrix.outerSize(); ++e)
            for (gsSparseMatrix<real_t>::InnerIterator it(fMatrix,e); it; ++it)
            {
                int k = it.row();
                int l = it.col();
                if (fMatrix.at(k,l) > 0.0)
                {
                    ff.coeffRef(k) += math::min(Rp.at(k),Rm.at(l))*fMatrix.at(k,l);
                }
                else
                {
                    ff.coeffRef(k) += math::min(Rm.at(k),Rp.at(l))*fMatrix.at(k,l);
                }
            }
        for (int k=0; k<M.outerSize();k++)
            ff.coeffRef(k) = ff.at(k)/M.at(k,k);
        //u_init = u_l + ff;
        u_init = u_h;
    }
    //----------------------------------------------------------------------------------------
    //TEST
    if (plot)
    {
        gsField<> sol4 = galerkinS.constructSolution(u_init);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>(sol4, "aaaINIT", 10000);//1000000);
    }
    //Time stepping loop
    gsParaviewCollection collection("aaaTIMEL");
    //FOR TEST ONLY
    /*
      BCs.addCondition(0,1, condition_type::dirichlet, &g);
      BCs.addCondition(0,2, condition_type::neumann, &g);
      BCs.addCondition(0,3, condition_type::neumann, &g);
      BCs.addCondition(0,4, condition_type::neumann, &g);

      gsCDRAssembler<real_t> galerkintry(mp,bases,BCs,f,coeff_AS,coeff_bK,coeff_c,
      dirichlet::elimination, iFace::glue, 0);
      galerkintry.assemble();
      K = galerkintry.matrix();
      rhs = galerkintry.rhs();
    //*/

    // Basic version, all parameters are NOT time dependent


    double tau = 2.0*EIGEN_PI / ntsteps;
    int orderFlag = 2;
    //gsSparseSolver<>::LU solver;

    gsMatrix<real_t> sol (u_init);
    //gsMatrix<real_t> M_inv = M.toDense().inverse();
    if (orderFlag == 2)
    {
        gsMatrix<real_t> rhs_it1 (galerkinK.rhs());
        gsMatrix<real_t> rhs_it2 (galerkinK.rhs());
        gsSparseMatrix<real_t> M_it (galerkinM.matrix());
        gsMatrix<real_t> u1 (u_init);
        int bases1d = bases[0].component(0).size();
        int bases2d = bases[0].component(1).size();
        for (int i=0; i< ntsteps; i++)
        {
            //u1 = sol - tau*M_inv*K*sol + tau*M_inv*rhs;
            //sol = 0.5*sol + 0.5*u1 - 0.5*M_inv*tau*K*u1 + 0.5*M_inv*tau*rhs;
            M_it = M;



            rhs_it1 = M*sol + tau*K*sol + tau*rhs;
            //BCs
            for (int k=1; k<=4;k++)
            {
                //this works on the parametric domain !!!!!!!!
                gsMatrix<real_t> anch;
                gsMatrix<real_t> fpts;
                gsGeometry<real_t>::uPtr bd = mp.patch(0).boundary(k);
                gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(k);
                anch = bd_basis->anchors();
                fpts = g.eval(  bd->eval(anch) );
                //gsDebugVar(fpts);
                gsGeometry<real_t>::uPtr geo = bd_basis->interpolateAtAnchors(fpts);
                //gsDebugVar(*geo);

                const gsMatrix<real_t> & dVals =  geo->coefs();
                for (int o = 0; o<bases[0].component(boxSide(k).direction()).size()/2;o++)
                {
                    if (k==3)
                    {
                        ///*
                        int l = bases1d-o-1;
                        rhs_it1.coeffRef(l) = dVals.at(l);
                        for (int e=0; e<M_it.outerSize(); ++e)
                            for (gsSparseMatrix<real_t>::InnerIterator it(M_it,e); it; ++it)
                            {
                                if (it.row() == l)
                                    M_it.coeffRef(l,it.col()) = 0.0;
                            }
                        M_it.coeffRef(l,l) = 1.0;
                        //*/
                    }
                    else if (k==4)
                    {
                        ///*
                        int l = bases1d-o-1;
                        rhs_it1.coeffRef(K.outerSize()-l-1) = dVals.at(bases1d-l-1);
                        for (int e=0; e<M_it.outerSize(); ++e)
                            for (gsSparseMatrix<real_t>::InnerIterator it(M_it,e); it; ++it)
                            {
                                if (it.row() == K.outerSize()-l-1)
                                    M_it.coeffRef(K.outerSize()-l-1,it.col()) = 0.0;
                            }
                        M_it.coeffRef(K.outerSize()-l-1,K.outerSize()-l-1) = 1.0;
                        //*/
                    }
                    else if (k==1)
                    {
                        int l = o;
                        rhs_it1.coeffRef(bases1d*l) = dVals.at(l);
                        for (int e=0; e<M_it.outerSize(); ++e)
                            for (gsSparseMatrix<real_t>::InnerIterator it(M_it,e); it; ++it)
                            {
                                if (it.row() == bases1d*l)
                                    M_it.coeffRef(bases1d*l,it.col()) = 0.0;
                            }
                        M_it.coeffRef(bases1d*l,bases1d*l) = 1.0;

                    }
                    else
                    {
                        ///*
                        int l =bases2d-o-1;
                        rhs_it1.coeffRef(bases1d*(l+1)-1) = dVals.at(l);
                        for (int e=0; e<M_it.outerSize(); ++e)
                            for (gsSparseMatrix<real_t>::InnerIterator it(M_it,e); it; ++it)
                            {
                                if (it.row() == bases1d*(l+1)-1)
                                    M_it.coeffRef(bases1d*(l+1)-1,it.col()) = 0.0;
                            }
                        M_it.coeffRef(bases1d*(l+1)-1,bases1d*(l+1)-1) = 1.0;
                        //*/
                    }
                }
            }

            M_it.makeCompressed();
            solver.compute(M_it);
            u1 = solver.solve(rhs_it1);
            u1 = AFC_Correction(&K, u1, bases1d, bases2d, &galerkinK, &galerkinM, &M_it,
                                &solver, tau, &mp, &bases, g);
            rhs_it2 = 0.5*M*sol + 0.5*M*u1 + 0.5*tau*K*u1 + 0.5*tau*rhs;

            //BCs
            for (int k=1; k<=4;k++)
            {
                //this works on the parametric domain !!!!!!!!
                gsMatrix<real_t> anch;
                gsMatrix<real_t> fpts;
                gsGeometry<real_t>::uPtr bd = mp.patch(0).boundary(k);
                gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(k);
                anch = bd_basis->anchors();
                fpts = g.eval(  bd->eval(anch) );
                //gsDebugVar(fpts);
                gsGeometry<real_t>::uPtr geo = bd_basis->interpolateAtAnchors(fpts);
                //gsDebugVar(*geo);

                const gsMatrix<real_t> & dVals =  geo->coefs();
                for (int o = 0; o<bases[0].component(boxSide(k).direction()).size()/2;o++)
                {
                    if (k==3)
                    {
                        ///*
                        int l = bases1d-o-1;
                        rhs_it2.coeffRef(l) = dVals.at(l);
                        //*/
                    }
                    else if (k==4)
                    {
                        ///*
                        int l = bases1d-o-1;
                        rhs_it2.coeffRef(K.outerSize()-l-1) = dVals.at(bases1d-l-1);
                        //*/
                    }
                    else if (k==1)
                    {
                        int l = o;
                        rhs_it2.coeffRef(bases1d*l) = dVals.at(l);
                    }
                    else
                    {
                        ///*
                        int l = bases2d-o-1;
                        rhs_it2.coeffRef(bases1d*(l+1)-1) = dVals.at(l);
                        //*/
                    }
                }
            }
            sol = solver.solve(rhs_it2);

            ///*
            // defect correction scheme
            sol = AFC_Correction(&K, sol, bases1d, bases2d, &galerkinK, &galerkinM, &M_it,
                                 &solver, tau/2.0, &mp, &bases, g); // change tau
            // plotting in ParaView
            //sol = u1;
            gsInfo << i << "/" << ntsteps << "\n";
            if (plot)
            {
                std::string fileName = "aaaTIME" + util::to_string(i);
                gsField<> sol5 = galerkinS.constructSolution(sol);
                gsWriteParaview<>(sol5, fileName, 1000);
                fileName = fileName + "0";
                collection.addTimestep(fileName,i,".vts");
            }
        }
    }
    else if (orderFlag == 3)
    {
        gsMatrix<real_t> u1 (u_init);
        gsMatrix<real_t> u2 (u_init);
        for (int i=0; i< ntsteps; i++)
        {
            //u1 = sol - tau*M_inv*K*sol + tau*M_inv*rhs;
            //u2 = 0.75*sol + 0.25*u1 - 0.25*tau*M_inv*K*u1 + 0.25*tau*M_inv*rhs;
            //sol = sol/3.0 + 2.0/3.0*u2 - 2.0/3.0*tau*M_inv*K*u2 + 2.0/3.0*tau*M_inv*rhs;
        }
    }
    else
    {
        gsInfo << "NOT IMPLEMENTED ORDER"<< "\n";
        return -1;
    }


    //Postprocessing
    if (plot)
    {
        collection.save();
        gsField<> sol5 = galerkinK.constructSolution(sol);
        gsInfo << "p= " << bases.degree() << "     h = "
               <<bases[0].component(0).size()*bases[0].component(1).size() << "\n";
        gsNormL2<real_t>* L2Norm;
        L2Norm=new gsNormL2<real_t>(sol5,init);
        gsInfo << "Errors L2:" << "\n"
               << "AFC  : " << L2Norm->compute()<< "\n";
        delete L2Norm;
        gsSeminormH1<real_t>* h1Seminorm;
        h1Seminorm = new gsSeminormH1 <real_t> (sol5,init);
        gsInfo << "Errors H1:" << "\n"
               << "AFC  : " << h1Seminorm->compute() << "\n";

        gsNormL<1,real_t>* L1Norm;
        L1Norm=new gsNormL<1,real_t>(sol5,init);
        real_t errorL1 = L1Norm->compute();
        gsInfo << "Errors L1:" << "\n"
               << "AFC  : " << errorL1*errorL1 << "\n"   ;
        delete L1Norm;


        real_t error = 0.0;
        for (int c = 0;c<bases[0].component(0).size()*bases[0].component(1).size();c++)
        {
            gsMatrix<real_t> point (2,1);
            point << c%bases[0].component(0).size() * 1.0/bases[0].component(0).size(),
                c/bases[0].component(0).size() * 1.0/bases[0].component(1).size();
            gsMatrix<real_t> res;
            g.eval_into(point,res);
            error = error + M.at(c,c)*(res.at(0)-sol.at(c))*(res.at(0)-sol.at(c));
        }
        error = math::sqrt(error);
        gsInfo << "Errors L2 Kuzmin approximation:" << "\n"
               << "AFC  : " << error << "\n";

        error = 0.0;
        for (int c = 0;c<bases[0].component(0).size()*bases[0].component(1).size();c++)
        {
            gsMatrix<real_t> point (2,1);
            point << c%bases[0].component(0).size() * 1.0/bases[0].component(0).size(),
                c/bases[0].component(0).size() * 1.0/bases[0].component(1).size();
            gsMatrix<real_t> res;
            g.eval_into(point,res);
            error = error + M.at(c,c)*math::abs(res.at(0)-sol.at(c));
        }
        gsInfo << "Errors L1 Kuzmin approximation:" << "\n"
               << "AFC  : " << error << "\n";

            
        delete h1Seminorm;
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( sol5, "aaaSOL", 10000);//1000000);
    }

    //gsInfo << M;

    return 0;
}


gsMatrix<real_t> AFC_Correction(gsSparseMatrix<real_t>* L, gsMatrix <real_t> sol, int bases1d, int bases2d,
                                gsCDRAssembler<real_t> *galerkinK, gsCDRAssembler<real_t> *galerkinM, gsSparseMatrix<real_t>* Ml,
                                gsSparseSolver<>::LU * solver, real_t tau,
                                gsMultiPatch<real_t>* mp, gsMultiBasis<>* bases, gsFunctionExpr<real_t>  g)
{

    gsMatrix<real_t> Pp (sol);
    gsMatrix<real_t> Pm (sol);
    gsMatrix<real_t> Qp (sol);
    gsMatrix<real_t> Qm (sol);
    gsMatrix<real_t> ff (sol);
    for (int k = 0; k < L->outerSize(); k++)
    {
        Pp.coeffRef(k) = 0.0;
        Pm.coeffRef(k) = 0.0;
        Qp.coeffRef(k) = 0.0;
        Qm.coeffRef(k) = 0.0;
        ff.coeffRef(k) = 0.0;
    }
    // low order solution - lumped mass approximation
    gsMatrix<real_t> u_dot_l (sol);
    gsMatrix<real_t> rhsloc (sol);
    rhsloc = *L*sol;
    //BCs

    for (int o = 0; (o<bases1d/2 || o<bases2d/2);o++)
    {
        if (o<bases1d/2)
        {
            int l = bases1d-o-1;
            rhsloc.coeffRef(l) = 0.0;

            l = bases1d-o-1;
            rhsloc.coeffRef(L->outerSize()-l-1) = 0.0;
        }
        if (o<bases2d/2)
        {
            int l = o;
            rhsloc.coeffRef(bases1d*l) = 0.0;

            l = bases2d-o-1;
            rhsloc.coeffRef(bases1d*(l+1)-1) = 0.0;
        }

    }
    gsSparseMatrix<real_t> f (galerkinK->matrix());
    f = f + galerkinM->matrix();
    u_dot_l = solver->solve(rhsloc);
    //for (int k = 0; k < K.outerSize(); k++)
    //{
    //  for (int l = 0; l < K.innerSize(); l++)
    //  {
    for (int e=0; e<L->outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(f,e); it; ++it)
        {
            int k = it.row();
            int l = it.col();
            //                              f.coeffRef(k,l) = M_it.at(k,l)*(u_dot_l.at(k)-u_dot_l.at(l))+
            //  (math::max(real_t(0.0),math::max(-galerkinK.matrix().at(k,l),
            //                              -galerkinK.matrix().at(l,k)))) * (sol.at(k)-sol.at(l));
            f.coeffRef(k,l) = galerkinM->matrix().at(k,l)*(u_dot_l.at(k)-u_dot_l.at(l))+
                (math::max<real_t>(0.0,math::max(-galerkinK->matrix().at(k,l),
                                                 -galerkinK->matrix().at(l,k)))) * (sol.at(k)-sol.at(l));
            if (l!= k)
            {
                Pp.coeffRef(k) += math::max<real_t>(0.0,f.at(k,l));
                Pm.coeffRef(k) += math::min<real_t>(0.0,f.at(k,l));

                Qp.coeffRef(k) = math::max(Qp.coeffRef(k),sol.at(l)-sol.at(k));
                Qm.coeffRef(k) = math::min(Qm.coeffRef(k),sol.at(l)-sol.at(k));

                //mod
                //Qp.coeffRef(l) = math::max(Qp.coeffRef(l),sol(k)-sol(l));
                //Qm.coeffRef(l) = math::min(Qm.coeffRef(l),sol(k)-sol(l));
            }
        }
    gsMatrix<real_t> Rp (sol);
    gsMatrix<real_t> Rm (sol);
    for (int k = 0; k < L->outerSize(); k++)
    {
        if (tau*Pp.at(k)==0)
            Rp.coeffRef(k) = 1.0;
        else
            Rp.coeffRef(k) = math::min<real_t>(1.0,Ml->at(k,k)*Qp.at(k)/(tau*Pp.at(k)));
        if (tau*Pm.at(k)==0)
            Rm.coeffRef(k) = 1.0;
        else
            Rm.coeffRef(k) = math::min<real_t>(1.0,Ml->at(k,k)*Qm.at(k)/(tau*Pm.at(k)));
    }
    //for (int k = 0; k < K.outerSize(); k++)
    //{
    //  for (int l = 0; l < K.innerSize(); l++)
    //  {
    for (int e=0; e<f.outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(f,e); it; ++it)
        {
            int k = it.row();
            int l = it.col();
            if (f.at(k,l) > 0.0)
            {
                ff.coeffRef(k) += math::min<real_t>(Rp.at(k),Rm.at(l))*f.at(k,l);
                //ff.coeffRef(l) -= f.at(k,l);//math::min(Rp.at(k),Rm.at(l))*f.at(k,l);
            }
            else if (f.at(k,l) < 0.0)
            {
                ff.coeffRef(k) += math::min<real_t>(Rm.at(k),Rp.at(l))*f.at(k,l);
                //ff.coeffRef(l) -= f.at(k,l);//math::min(Rm.at(k),Rp.at(l))*f.at(k,l);
            }
        }
    gsMatrix<real_t> rhs_loc (sol);
    rhs_loc = *Ml*sol+tau*ff;
    //BCs
    for (int k=1; k<=4;k++)
    {
        //this works on the parametric domain !!!!!!!!
        gsMatrix<real_t> anch;
        gsMatrix<real_t> fpts;
        gsGeometry<real_t>::uPtr bd = mp->patch(0).boundary(k);
        gsBasis<>::uPtr bd_basis = (*bases)[0].boundaryBasis(k);
        anch = bd_basis->anchors();
        fpts = g.eval(  bd->eval(anch) );
        //gsDebugVar(fpts);
        gsGeometry<real_t>::uPtr geo = bd_basis->interpolateAtAnchors(fpts);
        //gsDebugVar(*geo);

        const gsMatrix<real_t> & dVals =  geo->coefs();
        for (int o = 0; o<(*bases)[0].component(boxSide(k).direction()).size()/2;o++)
        {
            if (k==3)
            {
                int l = bases1d-o-1;
                rhs_loc.coeffRef(l) = dVals.at(l);
            }
            else if (k==4)
            {
                int l = bases1d-o-1;
                rhs_loc.coeffRef(L->outerSize()-l-1) = dVals.at(bases1d-l-1);
            }
            else if (k==1)
            {
                int l = o;
                rhs_loc.coeffRef(bases1d*l) = dVals.at(l);
            }
            else
            {
                int l = bases2d-o-1;
                rhs_loc.coeffRef(bases1d*(l+1)-1) = dVals.at(l);
            }
        }
    }
    sol = solver->solve(rhs_loc);
    return sol;
}
