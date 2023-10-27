// This file is not a part of the library!

// First attempt to implement AFC stabilization
// for convection diffusion problem with GISMO

// Convergence test for stable problem

// requires cleaning up and further work, not working correctly yet

// Author: A.Jaeschke, parts copied from ConvDiff.cpp

#include <iostream>

#include <gismo.h>
#include <gsAssembler/gsCDRAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsSeminormH1.h>
#include <gsUtils/gsPointGrid.h>
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

	// --- Example 5:Unit square, advection-diffusion, test of convergence ------------------------------------
      gsFunctionExpr<real_t>  f("-1*(60*x^2+24*x+6)+0.2*(20*x^3+12*x^2+6*x+4)", 2);
      gsFunctionExpr<real_t>  g("5*x^4+4*x^3+3*x^2+2*x+y", 2);
      gsFunctionExpr<real_t>  coeff_A("1.0","0","0","1.0", 2);
      gsFunctionExpr<real_t>  coeff_b("0.2","0.4", 2);
      gsFunctionExpr<real_t>  coeff_c("0", 2);
      gsFunctionExpr<real_t>  coeff_AS("1.0","0","0","1.0", 2);
      gsFunctionExpr<real_t>  coeff_bS("0","0", 2);
      gsFunctionExpr<real_t>  coeff_AK("0","0","0","0", 2);
      gsFunctionExpr<real_t>  coeff_bK("0.2","0.4", 2);
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

	gsSparseSolver<>::LU solver4;
	gsSparseSolver<>::LU solver5;

	//Added one refinement at the beginning
	bases.uniformRefine();
	bases.uniformRefine();
    int i = 2;

    gsMatrix<> testsL2,  testsH1;
    testsL2.setZero(i+numHref+1,5);
    testsH1.setZero(i+numHref+1,5);

    // Run tests
    do
    {
        // Setup the assemblers
		gsCDRAssembler<real_t> galerkin5(mp,bases,BCs,f,coeff_A,coeff_b,coeff_c,
                                        dirichlet::elimination, iFace::glue, stabilizerCDR::none);
		//------------------------------------
		gsBoundaryConditions<> BCs2;
	    //gsCDRAssembler<real_t> galerkinTRY(mp,bases,BCs2,f,coeff_A,coeff_b,coeff_c,
		//	dirichlet::none, iFace::glue, stabilizerCDR::none);
		gsCDRAssembler<real_t> galerkinS(mp,bases,BCs2,f,coeff_AS,coeff_bS,coeff_c,
			dirichlet::none, iFace::glue, stabilizerCDR::none);
		gsCDRAssembler<real_t> galerkinK(mp,bases,BCs2,f,coeff_AK,coeff_bK,coeff_c,
			dirichlet::none, iFace::glue, stabilizerCDR::none);

		//------------------------------------
		galerkin5.assemble();
		//galerkinTRY.assemble();
		galerkinS.assemble();
		galerkinK.assemble();

		gsMatrix<> solVector4;
		gsMatrix<> solVector5;

        if ( galerkin5.numDofs() )
        {
			// version with the AFC
			//int bases1d = (int)math::sqrt((double)bases.totalSize());
			gsMatrix<real_t> rhs (galerkinK.rhs());
			gsSparseMatrix<real_t> K (galerkinK.matrix());
			real_t delta;
			///*
			// Perform discrete upwinding (K:=K+D)
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
			if (i > 1)
			{
				K = K + galerkinS.matrix();
				int bases1d = bases[0].component(0).size();
				//int bases2d = bases[0].boundaryBasis(1)->size();
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
                        // free memory
                        const gsMatrix<real_t> & dVals =  geo->coefs();
						for (int l = 0; l<bases[0].component(boxSide(k).direction()).size();l++)
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
					}

					K.makeCompressed();
					solver4.compute(K);
					solVector4 = solver4.solve(rhs);
					// defect correction scheme
					///*
					int niter = 100;
					double tol = 0.00000001;
					gsMatrix<real_t> soltemp (galerkinK.rhs());
					for (int o = 0; o < niter; o++)
					{
						gsInfo<<"Iteration" << o << "\n";

						gsSparseMatrix<real_t> fMatrix (galerkinK.matrix());
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
								real_t temp = -(math::max<real_t>(0.0,math::max(galerkinK.matrix().at(k,l),
																galerkinK.matrix().at(l,k)))) * (solVector4.at(l)-solVector4.at(k));
								fMatrix.coeffRef(k,l) = temp;
								Pp.coeffRef(k) += math::max<real_t>(0.0,temp);
								Pm.coeffRef(k) += math::min<real_t>(0.0,temp);
								Qp.coeffRef(k) += math::max<real_t>(0.0,-temp);
								Qm.coeffRef(k) += math::min<real_t>(0.0,-temp);
								Qp.coeffRef(l) += math::max<real_t>(0.0,temp);
								Qm.coeffRef(l) += math::min<real_t>(0.0,temp);
								}
								else
								{
									fMatrix.coeffRef(k,l) = 0.0;
								}
							}
						gsMatrix<real_t> Rp (galerkinK.rhs());
						gsMatrix<real_t> Rm (galerkinK.rhs());
						for (int k = 0; k < K.outerSize(); k++)
						{
							if (Pp.at(k)==0)
								Rp.coeffRef(k) = 1.0;
							else
								Rp.coeffRef(k) = math::min<real_t>(1.0,Qp.at(k)/Pp.at(k));
							if (Pm.at(k)==0)
								Rm.coeffRef(k) = 1.0;
							else
								Rm.coeffRef(k) = math::min<real_t>(1.0,Qm.at(k)/Pm.at(k));
						}
						for (int e=0; e<fMatrix.outerSize(); ++e)
							for (gsSparseMatrix<real_t>::InnerIterator it(fMatrix,e); it; ++it)
							{
								int k = it.row();
								int l = it.col();
								if (l > k)
								{
								real_t temp = fMatrix.at(k,l);
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
						// BCs
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
					solver5.compute(galerkin5.matrix());
				solVector5 = solver5.solve( galerkin5.rhs() );

				}


			}
        //postprocessing ----------------------------------------

		gsField<> sol4 = galerkinK.constructSolution(solVector4);
		gsField<> sol5 = galerkin5.constructSolution(solVector5);

		if (plot && i== numHref)
        {
            // Write approximate and exact solution to paraview files
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>( sol5, "aaaSUPG__CD", 1000,true);
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>( sol4, "aaaAFC__CD", 1000,true);

            //result = system("paraview poisson_problem.pvd &");

        }

        /*gsDebugVar(i);
        gsDebugVar(testsL2.size());
        gsDebugVar(testsH1.size());
        */
		// Collect data
        testsL2(i,0)  =
            testsH1(i,0)= bases.totalSize();

        gsNormL2<real_t>* L2Norm;
        L2Norm=new gsNormL2<real_t>(sol4,g);
        //L2Norm->compute();
        testsL2(i,1)= L2Norm->compute() ;
        L2Norm=new gsNormL2<real_t>(sol5,g);
        //L2Norm->compute();
        testsL2(i,3)= L2Norm->compute() ;
        delete L2Norm;
        gsSeminormH1<real_t>* h1Seminorm;
        h1Seminorm = new gsSeminormH1 <real_t> (sol4,g);
        //h1Seminorm->compute();
        testsH1(i,1) = h1Seminorm->compute() ;
        h1Seminorm = new gsSeminormH1 <real_t> (sol5,g);
        h1Seminorm->compute();
		testsH1(i,3) = h1Seminorm->value();
		delete h1Seminorm;

        if (i > 0)
        {
            testsL2(i,2)= testsL2(i-1,1) / testsL2(i,1);
			testsL2(i,4)= testsL2(i-1,3) / testsL2(i,3);
            testsH1(i,2)= testsH1(i-1,1) / testsH1(i,1);
			testsH1(i,4)= testsH1(i-1,3) / testsH1(i,3);
        }

        bases.uniformRefine();
    }
    while ( i++ < numHref );


	    for(i = 1; i<= numHref; ++i)
    {   // Compute convergence rates
        testsL2(i,2)= math::log(testsL2(i,2))/math::log(2.0);
		testsL2(i,4)= math::log(testsL2(i,4))/math::log(2.0);
        testsH1(i,2)= math::log(testsH1(i,2))/math::log(2.0);
		testsH1(i,4)= math::log(testsH1(i,4))/math::log(2.0);
    }

    gsInfo << "Summary:\n\n";
    gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        AFC          |     NO STABILIZATION    \n";
    gsInfo << "    Dofs   |  L2 error  | conv. rate   \n" << testsL2  << "\n";


    gsInfo << " (deg= "<<bases[0].minDegree() <<")  |        AFC          |     NO STABILIZATION    \\n";
    gsInfo << "    Dofs   |  H1 error  | conv. rate   \n" << testsH1  << "\n";


    return result;
}


