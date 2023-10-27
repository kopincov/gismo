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


using namespace gismo;

int main(int argc, char *argv[])
{
	// Input options
	index_t ntsteps = 250;
	real_t tau = 0.001;
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
	cmd.addInt("n","nTimesteps", "Number of time steps", ntsteps);
	cmd.addReal("t","Timestep","Time step size",tau);
	cmd.addSwitch("plot", "Plot result in ParaView format", plot);

	try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

	//Time dependent problem definition:
    //gsFunctionExpr<real_t>  init("if( (x-0.5)*(x-0.5)+ (y-0.5)*(y-0.5) <= 0.0625,1,0)", 2);
    gsFunctionExpr<real_t>  init("if(x<1.0,10*sin(pi*x)*sin(pi*y),0)", 2);
    //gsFunctionExpr<real_t>  init("if(x<1.0,1,0)", 2);
    gsFunctionExpr<real_t>  f("0.0", 2);
    gsFunctionExpr<real_t>  g("0.0", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
    gsFunctionExpr<real_t>  coeff_cM("1.0", 2);
    gsFunctionExpr<real_t>  coeff_AS("0","0","0","0", 2);
    gsFunctionExpr<real_t>  coeff_bS("0","0", 2);
    gsFunctionExpr<real_t>  coeff_AK("0","0","0","0", 2);
    gsFunctionExpr<real_t>  coeff_bK("10.0","0.0", 2);

	//test of more complex domain
	/*
      gsMultiPatch<real_t> mp;
      gsReadFile<>("planar/lake.xml", mp);
      //gsReadFile<>("planar/deformedSquare_copy.xml", mp);
      //gsReadFile<>("planar/Square_m_neq_n.xml", mp);
      //gsReadFile<>("planar/chanel_test.xml", mp);
      //mp_p = gsReadFile<>("planar/blade_test.xml", mp);
      //mp_p = gsReadFile<>("planar/deformedSquare.xml", mp);
	//*/


	gsMultiPatch<real_t> mp( *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,3.0,1.0) );

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

	gsSparseMatrix<real_t> M (galerkinM.matrix());
	gsMatrix<real_t> rhs (galerkinK.rhs());
	gsSparseMatrix<real_t> K (galerkinK.matrix());
	K = K + galerkinS.matrix();

	//Evaluation of the initial condition
	gsMatrix<real_t> anch;
	gsMatrix<real_t> fpts;
	gsGeometry<real_t> & geom = mp.patch(0);
	gsBasis<real_t> & basis = bases[0];
	anch = basis.anchors();
	fpts = init.eval(  geom.eval(anch) );
	gsGeometry<real_t>::uPtr geo = basis.interpolateAtAnchors(fpts);
	const gsMatrix<real_t> Vals =  geo->coefs();
	gsMatrix<real_t> u_init (Vals);

	//TEST
	if (plot)
	{
		gsField<> sol4 = galerkinS.constructSolution(u_init);
		gsInfo<<"Plotting in Paraview...\n";
		gsWriteParaview<>( sol4, "aaaINIT", 1000, true);
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
		//int bases2d = bases[0].boundaryBasis(1)->size();
		for (int i=0; i< ntsteps; i++)
		{
			//u1 = sol - tau*M_inv*K*sol + tau*M_inv*rhs;
			//sol = 0.5*sol + 0.5*u1 - 0.5*M_inv*tau*K*u1 + 0.5*M_inv*tau*rhs;
			M_it = M;
			rhs_it1 = M*sol - tau*K*sol + tau*rhs;
			//BCs
			for (int k=1; k<=4;k++)
			{
				//this works on the parametric domain !!!!!!!!
				gsMatrix<real_t> anch2;
				gsMatrix<real_t> fpts2;
				gsGeometry<real_t>::uPtr bd = mp.patch(0).boundary(k);
				gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(k);
				anch2 = bd_basis->anchors();
				fpts2 = g.eval(  bd->eval(anch2) );
				//gsDebugVar(fpts);
				gsGeometry<real_t>::uPtr geo2 = bd_basis->interpolateAtAnchors(fpts2);
				//gsDebugVar(*geo);

				const gsMatrix<real_t> & dVals =  geo2->coefs();
				for (int l = 0; l<bases[0].component(boxSide(k).direction()).size();l++)
				{
					if (k==3)
					{
						/*
                          rhs_it1.coeffRef(l) = dVals.at(l);

                          for (int e=0; e<K.outerSize(); ++e)
                          for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
                          {
                          if (it.row() == l)
                          K.coeffRef(l,it.col()) = 0.0;
                          }
                          M_it.coeffRef(l,l) = 1.0;
						//*/
					}
					else if (k==4)
					{
						/*
                          rhs_it1.coeffRef(K.outerSize()-l-1) = dVals.at(bases1d-l-1);
                          for (int e=0; e<K.outerSize(); ++e)
                          for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
                          {
                          if (it.row() == K.outerSize()-l-1)
                          K.coeffRef(K.outerSize()-l-1,it.col()) = 0.0;
                          }
                          M_it.coeffRef(K.outerSize()-l-1,K.outerSize()-l-1) = 1.0;
						//*/
					}
					else if (k==1)
					{

						rhs_it1.coeffRef(bases1d*l) = dVals.at(l);
						for (int e=0; e<K.outerSize(); ++e)
							for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
							{
								if (it.row() == bases1d*l)
									K.coeffRef(bases1d*l,it.col()) = 0.0;
							}
                        M_it.coeffRef(bases1d*l,bases1d*l) = 1.0;

					}
					else
					{
						/*
                          rhs_it1.coeffRef(bases1d*(l+1)-1) = dVals.at(l);
                          for (int e=0; e<K.outerSize(); ++e)
                          for (gsSparseMatrix<real_t>::InnerIterator it(K,e); it; ++it)
                          {
                          if (it.row() == bases1d*(l+1)-1)
                          K.coeffRef(bases1d*(l+1)-1,it.col()) = 0.0;
                          }
                          M_it.coeffRef(bases1d*(l+1)-1,bases1d*(l+1)-1) = 1.0;
						//*/
					}
				}
			}

			M_it.makeCompressed();
			solver.compute(M_it);
			u1 = solver.solve(rhs_it1);
			rhs_it2 = 0.5*M*sol + 0.5*M*u1 - 0.5*tau*K*u1 + 0.5*tau*rhs;
			//BCs
			for (int k=1; k<=4;k++)
			{
				//this works on the parametric domain !!!!!!!!
				gsMatrix<real_t> anch3;
				gsMatrix<real_t> fpts3;
				gsGeometry<real_t>::uPtr bd = mp.patch(0).boundary(k);
				gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(k);
				anch3 = bd_basis->anchors();
				fpts3 = g.eval(  bd->eval(anch3) );
				//gsDebugVar(fpts);
				gsGeometry<real_t>::uPtr geo3 = bd_basis->interpolateAtAnchors(fpts3);
				//gsDebugVar(*geo);

				const gsMatrix<real_t> & dVals =  geo3->coefs();
				for (int l = 0; l<bases[0].component(boxSide(k).direction()).size();l++)
				{
					if (k==3)
					{
						/*
                          rhs_it1.coeffRef(l) = dVals.at(l);
						//*/
					}
					else if (k==4)
					{
						/*
                          rhs_it1.coeffRef(K.outerSize()-l-1) = dVals.at(bases1d-l-1);
						//*/
					}
					else if (k==1)
					{
						rhs_it2.coeffRef(bases1d*l) = dVals.at(l);
					}
					else
					{
						/*
                          rhs_it1.coeffRef(bases1d*(l+1)-1) = dVals.at(l);
						//*/
					}
				}
			}
			sol = solver.solve(rhs_it2);
			// plotting in ParaView
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
		gsField<> sol5 = galerkinS.constructSolution(sol);
		gsInfo<<"Plotting in Paraview...\n";
		gsWriteParaview<>( sol5, "aaaSOL", 1000);//, true);
		collection.save();
	}

	return 0;
}
