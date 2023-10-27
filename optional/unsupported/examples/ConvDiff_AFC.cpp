// This file is not a part of the library!

// First attempt to implement AFC stabilization
// for convection diffusion problem with GISMO

// requires cleaning up and further work, not working correctly yet
// in the compile else part the test of the artificial diffusion
// stablization and of the different assembly strategies

// Author: A.Jaeschke, parts copied from ConvDiff.cpp


#include <gismo.h>

#include <iostream>
#include <iomanip>

using namespace gismo;

int main(int argc, char *argv[])
{
	// Input options
	index_t numElevate  = 0;
	index_t numHref     = 2;
	index_t basisDegree = 0;
	bool plot       = false;
	// Flag for SUPG-Stabilization
	//bool Flag_Stabilization = 0;
	//Flag_Stabilization = 1;

    int result = 0;
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
	cmd.addSwitch("plot", "Plot result in ParaView format", plot);

	try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

	// --- Example 3: Unit square, advection-diffusion, internal layer skew to mesh --------------------------------------------

    gsFunctionExpr<real_t>  funExprF("0", 2);
    gsFunctionExpr<real_t>  g("if( y<=0.2-0.2*x,1,0)", 2);

    gsFunctionExpr<real_t>  coeff_A("0.0001","0","0","0.0001", 2);
    gsFunctionExpr<real_t>  coeff_b("sqrt(2.0)","sqrt(2.0)", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
    gsFunctionExpr<real_t>  coeff_AS("0.0001","0","0","0.0001", 2);
    gsFunctionExpr<real_t>  coeff_bS("0","0", 2);
    gsFunctionExpr<real_t>  coeff_AK("0","0","0","0", 2);
    gsFunctionExpr<real_t>  coeff_bK("sqrt(2.0)","sqrt(2.0)", 2);
    gsFunctionExpr<real_t>  coeff_cM("1.0", 2);

	//test of more complex domain
	/*
    gsMultiPatch<real_t> mp;
    //gsReadFile<>("planar/lake.xml", mp);
    //gsReadFile<>("planar/deformedSquare_copy.xml", mp);
    //gsReadFile<>("planar/Square_m_neq_n.xml", mp);
    gsReadFile<>("planar/chanel_test.xml", mp);
    //gsReadFile<>("planar/blade_test.xml", mp);
    //gsReadFile<>("planar/deformedSquare.xml", mp);
	//*/

	gsMultiPatch<real_t> mp( *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0) );
	gsBoundaryConditions<> BCs;
	for (gsMultiPatch<>::const_biterator
             bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
	{
		BCs.addCondition( *bit, condition_type::dirichlet, &g);
	}
	//
	// Set up and elevate discretization bases
	gsMultiBasis<> bases(mp);

	if (basisDegree)
		bases.setDegree(basisDegree);
	else if (numElevate)
		bases.degreeElevate(numElevate);
	//               bases.degreeReduce(1);

    // Construct the problem data
    gsBoundaryConditions<> BCs2;
    gsConvDiffRePde<> pde5(mp, BCs, &coeff_A,&coeff_b,&coeff_c, &funExprF);
    gsConvDiffRePde<> pdeS(mp, BCs, &coeff_AS,&coeff_bS,&coeff_c, &g);
    gsConvDiffRePde<> pdeK(mp, BCs, &coeff_AK,&coeff_bK,&coeff_c, &g);
    gsConvDiffRePde<> pdeM(mp, BCs, &coeff_AK,&coeff_bS,&coeff_cM, &g);

	// Run tests
	gsMatrix<> testsL2(numHref+1,7);
	testsL2.setZero();

	gsMatrix<> testsH1(numHref+1,7);
	testsH1.setZero();
	gsSparseSolver<>::LU solver4;
	gsSparseSolver<>::LU solver5;

	gsField<> sol4h ;
	gsField<> sol2h ;
	gsField<> solh ;
	gsField<> sol4hSUPG ;
	gsField<> sol2hSUPG ;
	gsField<> solhSUPG ;
	for (int i=0; i<numHref-2; i++)
		bases.uniformRefine();
	int i = numHref-2;
	do
	{
		// Setup the assemblers
		gsCDRAssembler<real_t> galerkin5(pde5, bases, dirichlet::elimination, iFace::glue, stabilizerCDR::SUPG);
		//------------------------------------
		//gsCDRAssembler<real_t> galerkinTRY(..)
		gsCDRAssembler<real_t> galerkinS(pdeS, bases, dirichlet::none, iFace::glue, stabilizerCDR::none);
		gsCDRAssembler<real_t> galerkinK(pdeK, bases, dirichlet::none, iFace::glue, stabilizerCDR::none);
		gsCDRAssembler<real_t> galerkinM(pdeM, bases, dirichlet::none, iFace::glue, stabilizerCDR::none);

        

		//------------------------------------
		galerkin5.assemble();
		//galerkinTRY.assemble();
		galerkinS.assemble();
		galerkinK.assemble();
		galerkinM.assemble();

		gsMatrix<> solVector4;
		gsMatrix<> solVector5;
        if ( galerkin5.numDofs() )
		{
			// version with the AFC
			gsMatrix<real_t> rhs (galerkinK.rhs());
			gsSparseMatrix<real_t> K (galerkinK.matrix());
            real_t delta;
			// outer is col
			///*
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
            //if (i== numHref)
            //{
            K = K + galerkinS.matrix();
            /* //BCs v1
            // Works only for m=n
            int bases1d = (int)math::sqrt((double)bases.totalSize());
            for (int k=0;k<K.outerSize();k++)
            {
            if (((k<bases1d || k>=K.outerSize()-bases1d)||(k % bases1d ==0)||(k % bases1d ==bases1d-1)))
            {
            for (int l=0;l<K.innerSize();l++)
            {
            K.coeffRef(k,l) = 0.0;
            }
            K.coeffRef(k,k) = 1.0;
            // this will work only for uniform mesh
            double ya = 1.0/((double)bases1d-1.0) * (double)(k / bases1d);
            double xa = 1.0/((double)bases1d-1.0) * (double)(k % bases1d);
            if ( ya<=0.2-0.2*xa)
            {
            rhs.coeffRef(k) = 1.0;
            }
            else
            {
            rhs.coeffRef(k) = 0.0;
            }
            }
            }*/
            // BCs v2
            int bases1d = bases[0].component(0).size();
            int bases2d = bases[0].component(1).size();
            /*
              for (int k=1; k<=4;k++)
              {
              //this works on the parametric domain !!!!!!!!
              gsMatrix<real_t> anch;
              gsMatrix<real_t> fpts;

              gsGeometry<real_t> bd = mp.patch(0).boundary(k);
              gsBasis<> * bd_basis = bases[0].boundaryBasis(k);
              anch = bd_basis->anchors();
              fpts = g.eval(  bd->eval(anch) );
              //gsDebugVar(fpts);
              gsGeometry<real_t> * geo = bd_basis->interpolateAtAnchors(fpts);
              //gsDebugVar(*geo);
              // free memory
              const gsMatrix<real_t> & dVals =  geo->coefs();

              for (int l = 0; l<bases[0].boundaryBasis(k)->size();l++)
              {
              if (k==3)
              {
              rhs.coeffRef(l) = dVals.at(l);
              for (int e=0; e<K.outerSize(); ++e)
              for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
              {
              if (it.row() == l)
              K.coeffRef(l,it.col()) = 0.0;
              }
              K.coeffRef(l,l) = 1.0;
              }
              else if (k==4)
              {
              rhs.coeffRef(K.outerSize()-l-1) = dVals.at(bases1d-l-1);
              for (int e=0; e<K.outerSize(); ++e)
              for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
              {
              if (it.row() == K.outerSize()-l-1)
              K.coeffRef(K.outerSize()-l-1,it.col()) = 0.0;
              }
              K.coeffRef(K.outerSize()-l-1,K.outerSize()-l-1) = 1.0;
              }
              else if (k==1)
              {
              rhs.coeffRef(bases1d*l) = dVals.at(l);
              for (int e=0; e<K.outerSize(); ++e)
              for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
              {
              if (it.row() == bases1d*l)
              K.coeffRef(bases1d*l,it.col()) = 0.0;
              }
              K.coeffRef(bases1d*l,bases1d*l) = 1.0;
              }
              else
              {
              rhs.coeffRef(bases1d*(l+1)-1) = dVals.at(l);
              for (int e=0; e<K.outerSize(); ++e)
              for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
              {
              if (it.row() == bases1d*(l+1)-1)
              K.coeffRef(bases1d*(l+1)-1,it.col()) = 0.0;
              }
              K.coeffRef(bases1d*(l+1)-1,bases1d*(l+1)-1) = 1.0;
              }
              }
              //free memory
              delete geo;
              delete bd_basis;
              }//*/
            ///*
            //BC Constrained L2
            //Constrained L2 data projection
            gsMatrix<real_t> u_l (galerkinK.rhs());
            gsMatrix<real_t> u_h (galerkinK.rhs());
            solver4.compute(galerkinM.matrix());
            u_h = solver4.solve(galerkinM.rhs());
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
            solver4.compute(M);
            u_l = solver4.solve(galerkinM.rhs());
            gsMatrix<real_t> u_boundary (u_l);
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
                u_boundary = u_l + ff;
            }
            for (int k=1; k<=4;k++)
            {
                //this works on the parametric domain !!!!!!!!

                for (int l = 0; l<bases[0].component(boxSide(k).direction()).size();l++)
                {
                    if (k==3)
                    {
                        rhs.coeffRef(l) = 1.0;//u_boundary.at(l);
                        for (int e=0; e<K.outerSize(); ++e)
                            for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
                            {
                                if (it.row() == l)
                                    K.coeffRef(l,it.col()) = 0.0;
                            }
                        K.coeffRef(l,l) = 1.0;
                    }
                    else if (k==4)
                    {
                        rhs.coeffRef(K.outerSize()-l-1) = u_boundary.at(K.outerSize()-l-1);
                        for (int e=0; e<K.outerSize(); ++e)
                            for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
                            {
                                if (it.row() == K.outerSize()-l-1)
                                    K.coeffRef(K.outerSize()-l-1,it.col()) = 0.0;
                            }
                        K.coeffRef(K.outerSize()-l-1,K.outerSize()-l-1) = 1.0;
                    }
                    else if (k==1)
                    {
                        rhs.coeffRef(bases1d*l) = u_boundary.at(bases1d*l);
                        for (int e=0; e<K.outerSize(); ++e)
                            for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
                            {
                                if (it.row() == bases1d*l)
                                    K.coeffRef(bases1d*l,it.col()) = 0.0;
                            }
                        K.coeffRef(bases1d*l,bases1d*l) = 1.0;
                    }
                    else
                    {
                        rhs.coeffRef(bases1d*(l+1)-1) = u_boundary.at(bases1d*(l+1)-1);
                        for (int e=0; e<K.outerSize(); ++e)
                            for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
                            {
                                if (it.row() == bases1d*(l+1)-1)
                                    K.coeffRef(bases1d*(l+1)-1,it.col()) = 0.0;
                            }
                        K.coeffRef(bases1d*(l+1)-1,bases1d*(l+1)-1) = 1.0;
                    }
                }
            }
            //	gsField<> * sol4 = galerkinS.constructSolution(u_boundary);
            //     gsInfo<<"Plotting in Paraview...\n";
            //    gsWriteParaview<>( *sol4, "aaaINIT",100000);
            //delete sol4;//*/
            //----------------------------------------------------------------------------------------
            K.makeCompressed();
            solver4.compute(K);
            solVector4 = solver4.solve(rhs);
            // defect correction scheme
            //*
            int niter = 100;
            double tol = 0.0001;
            gsMatrix<real_t> soltemp (galerkinK.rhs());
            for (int o = 0; o < niter; o++)
            {
                gsInfo<<"Iteration" << o << "\n";

                gsSparseMatrix<real_t> f (galerkinK.matrix());
                gsMatrix<real_t> Pp (galerkinK.rhs());
                gsMatrix<real_t> Pm (galerkinK.rhs());
                gsMatrix<real_t> Qp (galerkinK.rhs());
                gsMatrix<real_t> Qm (galerkinK.rhs());
                gsMatrix<real_t> ff (galerkinK.rhs());
                for (int k = 0; k < K.outerSize(); k++)
                {
                    Pp.coeffRef(k) = 0.0;
                    Pm.coeffRef(k) = 0.0;
                    Qp.coeffRef(k) = 0.0;
                    Qm.coeffRef(k) = 0.0;
                    ff.coeffRef(k) = 0.0;
                }

                for (int e=0; e<K.outerSize(); ++e)
                    for (gsSparseMatrix<real_t>::InnerIterator it(galerkinK.matrix(),e); it; ++it)
                    {
                        int k = it.row();
                        int l = it.col();
                        if (l > k || (l<k && galerkinK.matrix().at(l,k)==0)) //make sure that it still works if the matrix is not symmetric (sparsity)
                        {
                            real_t temp = -(math::max((real_t)(0.0),math::max(galerkinK.matrix().at(k,l),
                                                                              galerkinK.matrix().at(l,k)))) * (solVector4.at(l)-solVector4.at(k));
                            f.coeffRef(k,l) = temp;
                            Pp.coeffRef(k) += math::max((real_t)(0.0),temp);
                            Pm.coeffRef(k) += math::min((real_t)(0.0),temp);
                            Qp.coeffRef(k) += math::max((real_t)(0.0),-temp);
                            Qm.coeffRef(k) += math::min((real_t)(0.0),-temp);
                            Qp.coeffRef(l) += math::max((real_t)(0.0),temp);
                            Qm.coeffRef(l) += math::min((real_t)(0.0),temp);
                        }
                        else
                        {
                            f.coeffRef(k,l) = 0.0;
                        }
                    }

                gsMatrix<real_t> Rp (galerkinK.rhs());
                gsMatrix<real_t> Rm (galerkinK.rhs());
                for (int k = 0; k < K.outerSize(); k++)
                {
                    if (Pp.at(k)==0)
                        Rp.coeffRef(k) = 1.0;
                    else
                        Rp.coeffRef(k) = math::min((real_t)(1.0),Qp.at(k)/Pp.at(k));
                    if (Pm.at(k)==0)
                        Rm.coeffRef(k) = 1.0;
                    else
                        Rm.coeffRef(k) = math::min((real_t)(1.0),Qm.at(k)/Pm.at(k));
                }

                for (int e=0; e<f.outerSize(); ++e)
                    for (gsSparseMatrix<real_t>::InnerIterator it(f,e); it; ++it)
                    {
                        int k = it.row();
                        int l = it.col();
                        if (l > k)
                        {
                            real_t temp = f.at(k,l);
                            if (temp > 0.0)
                            {
                                ff.coeffRef(k) += Rp.at(k)*temp;
                                ff.coeffRef(l) -= Rp.at(k)*temp;
                            }
                            else
                            {
                                ff.coeffRef(k) += Rm.at(k)*temp;
                                ff.coeffRef(l) -= Rm.at(k)*temp;
                            }
                        }
                    }

                rhs = galerkinK.rhs()+ff-K*solVector4;

                //BCs
                //int bases1d = bases[0].boundaryBasis(1)->size();
                for (int k=0;k<K.outerSize();k++)
                {
                    if (((k<bases1d || k>=K.outerSize()-bases1d)||
                         (k % bases1d ==0)||(k % bases1d == bases1d-1)))
                    {
                        rhs.coeffRef(k) = 0.0;
                    }
                }
                soltemp = solver4.solve(rhs);
                solVector4 = solVector4 + soltemp;

                gsInfo << soltemp.norm() <<"  "<< solVector4.norm() << "\n";
                if (soltemp.norm() / solVector4.norm() < tol)
                    break;
            }
            //*/
            gsInfo << "Degrees of freedom: " << bases1d*bases2d << "\n";
            gsInfo << "Degree: " << basisDegree << "\n";
            solver5.compute(galerkin5.matrix());
            solVector5 = solver5.solve( galerkin5.rhs() );
            //}


		}
		//postprocessing ----------------------------------------
		if (i == numHref -2)
		{
			sol4h = galerkinK.constructSolution(solVector4);
			sol4hSUPG = galerkin5.constructSolution(solVector5);
		}
		if (i == numHref -1)
		{
			sol2h = galerkinK.constructSolution(solVector4);
			sol2hSUPG = galerkin5.constructSolution(solVector5);
		}
		if (i == numHref)
		{
			solh = galerkinK.constructSolution(solVector4);
			solhSUPG = galerkin5.constructSolution(solVector5);
			gsInfo << "p= " << bases.degree() << "      (finest) h = "
                   <<bases[0].component(0).size()*bases[0].component(1).size() << "\n";
            //this is not very reliable
			real_t rate = math::log(sol2h.distanceL2(sol4h)/solh.distanceL2(sol2h))/log(2.0);
			real_t rateSUPG = math::log(sol2hSUPG.distanceL2(sol4hSUPG)/
                                        solhSUPG.distanceL2(sol2hSUPG))/log(2.0);
			gsInfo << "Rates of convergence" << "\n"
                   << "AFC  : " << rate << "\n"
                   << "SUPG : " << rateSUPG << "\n";

		}
		if (plot && i== numHref)
		{

			gsField<> sol4 = galerkinK.constructSolution(solVector4);
			gsField<> sol5 = galerkin5.constructSolution(solVector5);
			// Write approximate and exact solution to paraview files
			gsInfo<<"Plotting in Paraview...\n";
			gsWriteParaview<>( sol5, "aaa_SUPG_CD", 100000);
			gsInfo<<"Plotting in Paraview...\n";
			gsWriteParaview<>( sol4, "aaa_AFC__CD", 100000);
			// Run paraview
			//result = system("paraview poisson_problem.pvd &");
		}
		bases.uniformRefine();
	}
	while ( i++ < numHref );
	return result;
}
