// This file is not a part of the library!

// NOT IMPLEMENTED YET

// Implementation of SSP RK method for the Euler equations, AFC stabilization used

// Author: A.Jaeschke, parts copied from ConvDiff.cpp

#define _USE_MATH_DEFINES

#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>
#include <gsAssembler/gsCDRAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <iomanip>
#include <gismo.h>
#include <time.h>


using namespace gismo;
gsSparseMatrix<real_t> AddABlock(int sizei, int sizej, int ibegin,
                                 int jbegin, gsSparseMatrix<real_t> A);
gsMatrix<real_t> A_jx(gsMatrix<real_t> u);
gsMatrix<real_t> A_jy(gsMatrix<real_t> u);
gsMatrix<real_t> F_x(gsMatrix<real_t> u);
gsMatrix<real_t> F_y(gsMatrix<real_t> u);
gsMatrix<real_t> UtoW (gsMatrix<real_t> n, gsMatrix<real_t> u);
gsMatrix<real_t> WtoU (gsMatrix<real_t> n, gsMatrix<real_t> w);
gsMatrix<real_t> nA (gsMatrix<real_t> n, gsMatrix<real_t> u, gsMatrix<real_t> u_inf);
gsMatrix<real_t> F_n_(gsMatrix<real_t> n, gsMatrix<real_t> u, gsMatrix<real_t> u_inf);
gsMatrix<real_t> F_n_wall (gsMatrix<real_t> n, gsMatrix<real_t> u);
gsMatrix<real_t> F_n_subin (gsMatrix<real_t> n, gsMatrix<real_t> u,
                            real_t rho, real_t p, real_t v_tangential);
gsMatrix<real_t> F_n_subout (gsMatrix<real_t> n, gsMatrix<real_t> u, real_t p_out);
gsMatrix<real_t> F_n_supin (gsMatrix<real_t> n, real_t rho, real_t v_x,
                            real_t v_y, real_t E);
gsMatrix<real_t> F_n_supout (gsMatrix<real_t> n, gsMatrix<real_t> u);
gsMatrix<real_t> F_n_type (gsMatrix<real_t> n, gsMatrix<real_t> u, std::string bc);
gsMatrix<real_t> F_boundary (gsMultiPatch<real_t> mp, gsMultiBasis<> bases, gsField<> & U1,
                             gsField<> & U2, gsField<> & U3, gsField<> & U4, std::string bc1,
                             std::string bc2, std::string bc3, std::string bc4);
gsMatrix<real_t> ConstL2(gsMatrix <real_t> INrhs, gsSparseMatrix <real_t> M, gsSparseMatrix <real_t> Ml );
gsMatrix<real_t> outer_normal(gsMatrix<real_t> points, gsMultiPatch<real_t> mp, boxSide a);
gsMatrix<real_t> mapping(gsMatrix<real_t> points, gsMultiPatch<real_t> mp);

int flag = 0;

int main(int argc, char *argv[])
{
    /// begin: Input options --------------------------------------------
    index_t numHref     = 0;
    index_t basisDegree = 0;
    bool plot       = false;
    index_t nsamples    = 10000;
    // ----------------------------
    index_t ntsteps = 300;
    real_t tau = 0.001;
    index_t SSPRKorder = 1;
    index_t plot_sparsity = 1;
    std::string path = "planar/Ubend.xml";
    std::string i1 = "if(x<0,1.0,0.125)";
    std::string i2 = "0.0";
    std::string i3 = "0.0";
    std::string i4 = "if(x<0,1.0,0.1)";
    std::string bc1 = "wall";
    std::string bc2 = "wall";
    std::string bc3 = "wall";
    std::string bc4 = "wall";
    index_t n_threads = 1;
  
    /// Possible types of BCs -------------------------------------------
    // "wall"
    // "subin rho p v_tangential"
    // "subout p" 
    // "supin rho v_x v_y p"
    // "supout"
    // "periodic"
    /// -----------------------------------------------------------------
	
    gsCmdLine cmd("Testing compressible Euler problem.");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("p","degree",
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    cmd.addInt("S","samples",
               "Number of samples to use for plotting",
               nsamples);
    cmd.addInt("n","nTimesteps", "Number of time steps", ntsteps);
    cmd.addReal("t","Timestep","Time step size",tau);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addInt("o","SSPRKorder", "Order of used SSP-RK method, 1-Euler, 2-SSPRK order 2, 3-SSPRK order 3", SSPRKorder);
    cmd.addInt("s","plotSparsity", "Plotting sparsity, number of time-steps between consecutive plots", plot_sparsity);
    cmd.addString("g","geometry","Path to the geometry",path);
    cmd.addString("1","initial_density","Inital condition for the density",i1);
    cmd.addString("2","initial_velocity_x","Inital condition for the x component of the velocity",i2);
    cmd.addString("3","initial_velocity_y","Inital condition for the y component of the velocity",i3);
    cmd.addString("4","initial_pressure","Inital condition for the pressure",i4);
    cmd.addString("a","bc_west","East boundary condition - see specification",bc1);
    cmd.addString("b","bc_east","West boundary condition - see specification",bc2);
    cmd.addString("c","bc_south","South boundary condition - see specification",bc3);
    cmd.addString("d","bc_north","North boundary condition - see specification",bc4);
    cmd.addInt("j","nThreads", "Number of threads used", n_threads);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    /// end: Input options ----------------------------------------------
#ifdef _OPENMP
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(n_threads);

    double t;
#endif
    real_t gamma = 1.4; // Be careful it is redefined several times!!!!c
  
    /// begin: Configuration of the matrix coefficients -----------------
    /// Basic version, all parameters are NOT time dependent
    gsFunctionExpr<real_t>  f("0.0", 2);
    gsFunctionExpr<real_t>  g("0.0", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
    gsFunctionExpr<real_t>  coeff_cM("1.0", 2);
    gsFunctionExpr<real_t>  coeff_AS("0.0","0","0","0.0", 2);
    gsFunctionExpr<real_t>  coeff_bS("0","0", 2);
    gsFunctionExpr<real_t>  coeff_AK("0","0","0","0", 2);
    gsFunctionExpr<real_t>  coeff_bCx("1.0","0.0", 2);
    gsFunctionExpr<real_t>  coeff_bCy("0.0","1.0", 2);
    /// end: Configuration of the matrix coefficients -------------------
    
    /// begin: Geometry configuration -----------------------------------
    gsMultiPatch<real_t> mp;
    gsReadFile<>(path, mp);
    gsMultiBasis<> bases(mp);
    /// end: Geometry configuration -------------------------------------
	
    /// begin: Refinement -----------------------------------------------
    if (basisDegree)
        bases.setDegree(basisDegree);

    gsSparseSolver<>::LU solver;
    for (index_t i=0; i<numHref ; i++)
        bases.uniformRefine();
    /// end: Refinement -------------------------------------------------
#ifdef _OPENMP
    t = omp_get_wtime();
#endif
    /// begin: Basic matrices construction ------------------------------
    gsBoundaryConditions<> BCs;
    gsCDRAssembler<real_t> galerkinCx(mp,bases,BCs,f,coeff_AK,coeff_bCx,coeff_c,
                                      dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinCx.assemble();
    gsCDRAssembler<real_t> galerkinCy(mp,bases,BCs,f,coeff_AK,coeff_bCy,coeff_c,
                                      dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinCy.assemble();
    gsCDRAssembler<real_t> galerkinM(mp,bases,BCs,f,coeff_AK,coeff_bS,coeff_cM,
                                     dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinM.assemble();
    /// end: Basic matrices construction --------------------------------

	
    /// begin: Basic global matrices construction + space delcaration----
    index_t basisSize = bases.totalSize();
    gsSparseMatrix<real_t> M (4*basisSize, 4*basisSize);
    M = AddABlock(4*basisSize, 4*basisSize, 0,0,galerkinM.matrix()) +
        AddABlock(4*basisSize, 4*basisSize, basisSize,basisSize ,galerkinM.matrix()) +
        AddABlock(4*basisSize, 4*basisSize, 2*basisSize, 2* basisSize ,galerkinM.matrix()) +
        AddABlock(4*basisSize, 4*basisSize, 3*basisSize, 3* basisSize ,galerkinM.matrix());
    gsSparseMatrix<real_t> K (4*basisSize, 4*basisSize);
    gsSparseMatrix<real_t> L (4*basisSize, 4*basisSize);
    gsSparseMatrix<real_t> Ml_ (galerkinM.matrix());
#ifdef _OPENMP
    t = omp_get_wtime() - t;
	std::cout << "Basic matrices construction took " << t << "s" << std::endl;
    /// end: Basic global matrices construction + space delcaration------
	
	/// begin: Mass lumping ---------------------------------------------
	t = omp_get_wtime(); 
	double t0 = omp_get_wtime();
#endif
	#pragma omp parallel
	{
	#pragma omp for 
    for (int e=0; e<Ml_.outerSize(); ++e) //cols
        for (gsSparseMatrix<real_t>::InnerIterator it(Ml_,e); it; ++it)
        {
            int k = it.row();
            //int l = it.col();
            if((it.value() != 0)&& k!=e)
            {
                Ml_.coeffRef(k,k) +=it.value();
                it.valueRef() = 0.0;
            }
        }
	}
#ifdef _OPENMP
	    t0 = omp_get_wtime() - t0;
	std::cout << "SUB Mass lumping took " << t0 << "s" << std::endl;
#endif
    gsSparseMatrix<real_t> Ml (4*basisSize, 4*basisSize);
    Ml = AddABlock(4*basisSize, 4*basisSize, 0,0,Ml_) +
        AddABlock(4*basisSize, 4*basisSize, basisSize,basisSize ,Ml_) +
        AddABlock(4*basisSize, 4*basisSize, 2*basisSize, 2* basisSize ,Ml_) +
        AddABlock(4*basisSize, 4*basisSize, 3*basisSize, 3* basisSize ,Ml_);
        
    solver.compute(Ml);
#ifdef _OPENMP
    t = omp_get_wtime() - t;
	std::cout << "Mass lumping took " << t << "s" << std::endl;
    /// end: Mass lumping -----------------------------------------------

	   
    /// begin: Initial condition with constrained L2 appied--------------
    t = omp_get_wtime();
#endif
    i1 = "(" + i1 + ")";
    i2 = "(" + i2 + ")";
    i3 = "(" + i3 + ")";
    i4 = "(" + i4 + ")";
    gsFunctionExpr<real_t>  init1(i1, 2);
    gsCDRAssembler<real_t> galerkinIN1(mp,bases,BCs,init1,coeff_AK,coeff_bS,coeff_cM,
                                       dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinIN1.assemble();
    gsMatrix<real_t> u_init1 = ConstL2(galerkinIN1.rhs(), galerkinM.matrix(), Ml_);
  
    i2 = i2 + "*" + i1;
    gsFunctionExpr<real_t>  init2(i2, 2);
    gsCDRAssembler<real_t> galerkinIN2(mp,bases,BCs,init2,coeff_AK,coeff_bS,coeff_cM,
                                       dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinIN2.assemble();
    gsMatrix<real_t> u_init2 = ConstL2(galerkinIN2.rhs(), galerkinM.matrix(), Ml_);
  
    i3 = i3 + "*" + i1;
    gsFunctionExpr<real_t>  init3(i3, 2);
    gsCDRAssembler<real_t> galerkinIN3(mp,bases,BCs,init3,coeff_AK,coeff_bS,coeff_cM,
                                       dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinIN3.assemble();
    gsMatrix<real_t> u_init3 = ConstL2(galerkinIN3.rhs(), galerkinM.matrix(), Ml_);
    
    i4 = i4 + "/(" + std::to_string(gamma) + "-1.0)+0.5*(" + i2 + "*" + i2 + "+" + i3 + "*" + i3 + ")/" + i1;
    gsFunctionExpr<real_t>  init4(i4, 2);
    gsCDRAssembler<real_t> galerkinIN4(mp,bases,BCs,init4,coeff_AK,coeff_bS,coeff_cM,
                                       dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinIN4.assemble();
    gsMatrix<real_t> u_init4 = ConstL2(galerkinIN4.rhs(), galerkinM.matrix(), Ml_);
    /// end: Initial condition with constrained L2 appied----------------
        
    /// begin: Creating a common vector of initial conditions -----------
    gsMatrix <real_t> u_init (4*basisSize, 1);
    for (index_t k=0; k<4*basisSize; k++)
    {
        if (k<basisSize)
            u_init.coeffRef(k) = u_init1.at(k);
        else if (k<2*basisSize)
            u_init.coeffRef(k) = u_init2.at(k-basisSize);
        else if (k<3*basisSize)
            u_init.coeffRef(k) = u_init3.at(k-2*basisSize);
        else
            u_init.coeffRef(k) = u_init4.at(k-3*basisSize);
    }
#ifdef _OPENMP
    t = omp_get_wtime() - t;
	std::cout << "Initial condition evaluation took " << t << "s" << std::endl;
#endif
    /// end: Creating a common vector of initial conditions -------------

    /// begin: Plot of the initial conditions ---------------------------
    if (plot)      
    {  
        for (int j=0;j<int(basisSize);j++)
        {
            u_init2.coeffRef(j) = u_init2.at(j)/u_init1.at(j);
            u_init3.coeffRef(j) = u_init3.at(j)/u_init1.at(j);
            u_init4.coeffRef(j) = (gamma-1.0)*(u_init4.at(j)-0.5*u_init1.at(j)
                                               *(u_init2.at(j)*u_init2.at(j)+u_init3.at(j)*u_init3.at(j)));
        }
        gsField<> sol1 = galerkinIN1.constructSolution(u_init1);
        gsWriteParaview<>( sol1, "Initial_density", nsamples);
        gsField<> sol2 = galerkinIN2.constructSolution(u_init2);
        gsWriteParaview<>( sol2, "Initial_velocity_x", nsamples);
        gsField<> sol3 = galerkinIN3.constructSolution(u_init3);
        gsWriteParaview<>( sol3, "Initial_velocity_y", nsamples);
        gsField<> sol4 = galerkinIN4.constructSolution(u_init4);
        gsWriteParaview<>( sol4, "Initial_pressure", nsamples);
    }
    /// end: Plot of the initial conditions -----------------------------
  
    /// begin: Plotting collection definitions --------------------------  
    gsParaviewCollection collection("Time_density");
    gsParaviewCollection collection2("Time_velocity_x");
    gsParaviewCollection collection3("Time_velocity_y");
    gsParaviewCollection collection4("Time_pressure");
    /// end: Plotting collection definitions ----------------------------
    
    /// begin: Declaration of necessary vectors for time stepping -------
    gsMatrix<real_t> sol (u_init);
    // We assume one order of SSP-RK method - EULER!!!, further orders to be implemented
    gsMatrix<real_t> rhs_it1 (galerkinCx.rhs());
    gsMatrix<real_t> rhs_it2 (galerkinCx.rhs());
    gsMatrix<real_t> u1 (u_init);
    /// end: Declaration of necessary vectors for time stepping ---------
  
    /// begin: Time stepping loop  --------------------------------------
    for (int i=0; i< ntsteps; i++)
    {
	
	
        /// begin: Construction of the K matrix ---------------------
#ifdef _OPENMP
        t = omp_get_wtime();
#endif
        std::vector<Eigen::Triplet<real_t> > tripletList;
        tripletList.reserve(galerkinCx.matrix().nonZeros());
		//#pragma omp parallel
		//{
		//#pragma omp for
        for (int e=0; e<galerkinCx.matrix().outerSize(); ++e)
            for (gsSparseMatrix<real_t>::InnerIterator it(galerkinCx.matrix(),e); it; ++it)
            {
                int k = it.row();
                int l = it.col();
                gsMatrix<real_t> solat (4,1);
                solat << sol.at(l),sol.at(l+basisSize),sol.at(l+2*basisSize),sol.at(l+3*basisSize);
                gsMatrix<real_t> k_kl = galerkinCx.matrix().at(l,k)*A_jx(solat) +
                    galerkinCy.matrix().at(l,k)*A_jy(solat);
                for (int o=0;o<4;o++)
                    for (int oo=0;oo<4;oo++)
                        tripletList.push_back(Eigen::Triplet<real_t> (k+o*basisSize,l+oo*basisSize,k_kl.coeff(o,oo)));
            }
		//}
        K.setFromTriplets(tripletList.begin(), tripletList.end());
#ifdef _OPENMP
    t = omp_get_wtime() - t;
	std::cout << "K matrix construction took " << t << "s" << std::endl;
        /// end: Construction of the K matrix -----------------------

	              
        /// begin: Construction of the L matrix ---------------------
        t = omp_get_wtime();
#endif
        L = K;
        #pragma omp parallel
		{
		#pragma omp for
        for (int e=0; e<galerkinCy.matrix().outerSize(); ++e)
        {
            for (gsSparseMatrix<real_t>::InnerIterator it(galerkinCy.matrix(),e); it; ++it)
            {
                int k = it.row();
                int l = it.col();
                if (l>k)
                {
                    real_t e_klx = (galerkinCx.matrix().coeff(l,k) - galerkinCx.matrix().coeff(k,l))/2.0;
                    real_t e_kly = (galerkinCy.matrix().coeff(l,k) - galerkinCy.matrix().coeff(k,l))/2.0;
                    real_t v_klx = (math::sqrt(sol.at(k))*sol.at(k+basisSize)/sol.at(k) +
                                    math::sqrt(sol.at(l))*sol.at(l+basisSize)/sol.at(l))/
                        (math::sqrt(sol.at(k))+math::sqrt(sol.at(l)));
                    real_t v_kly = (math::sqrt(sol.at(k))*sol.at(k+2*basisSize)/sol.at(k) +
                                    math::sqrt(sol.at(l))*sol.at(l+2*basisSize)/sol.at(l))/
                        (math::sqrt(sol.at(k))+math::sqrt(sol.at(l)));
                    real_t v_kl2 = v_klx*v_klx+v_kly*v_kly;
                    real_t Hk = (gamma*sol.at(k+3*basisSize) + (1.0-gamma)*(sol.at(k+basisSize)*
                                                                            sol.at(k+basisSize)+sol.at(k+2*basisSize)*sol.at(k+2*basisSize))/(2.0*sol.at(k)))/sol.at(k);
                    real_t Hl = (gamma*sol.at(l+3*basisSize) + (1.0-gamma)*(sol.at(l+basisSize)*
                                                                            sol.at(l+basisSize)+sol.at(l+2*basisSize)*sol.at(l+2*basisSize))/(2.0*sol.at(l)))/sol.at(l);
                    real_t Hkl = (math::sqrt(sol.at(k))*Hk + math::sqrt(sol.at(l))*Hl)/
                        (math::sqrt(sol.at(k))+math::sqrt(sol.at(l)));
                    real_t ckl = math::sqrt((gamma-1.0)*(Hkl-v_kl2/2.0));
                    real_t v_kl = (e_klx*v_klx+e_kly*v_kly)/math::sqrt(e_klx*e_klx+e_kly*e_kly);
                    real_t dij = math::sqrt(e_klx*e_klx+e_kly*e_kly)*(math::abs(v_kl)+ckl);
                    for (int o = 0; o < 4; o++)
                    {
                        L.coeffRef(k+o*basisSize,l+o*basisSize) += dij;
                        L.coeffRef(l+o*basisSize,k+o*basisSize) += dij;
                        L.coeffRef(k+o*basisSize,k+o*basisSize) -= dij;
                        L.coeffRef(l+o*basisSize,l+o*basisSize) -= dij;
                    }
                }
            }
		}
		}
#ifdef _OPENMP
	t = omp_get_wtime() - t;
	std::cout << "L matrix construction took " << t << "s" << std::endl;
        /// end: Construction of the L matrix -----------------------
	        
        /// begin: Construction of the boundary fluxes --------------
        t = omp_get_wtime();
#endif
        gsMatrix <real_t> sol1 (basisSize,1);
        gsMatrix <real_t> sol2 (basisSize,1);
        gsMatrix <real_t> sol3 (basisSize,1);
        gsMatrix <real_t> sol4 (basisSize,1);
        for (int j=0;j<int(basisSize);j++)
        {
            sol1.coeffRef(j) = sol.at(j);
            sol2.coeffRef(j) = sol.at(j+basisSize);
            sol3.coeffRef(j) = sol.at(j+2*basisSize);
            sol4.coeffRef(j) = sol.at(j+3*basisSize);
        }
        gsField<> U1 = galerkinM.constructSolution(sol1);
        gsField<> U2 = galerkinM.constructSolution(sol2);
        gsField<> U3 = galerkinM.constructSolution(sol3);
        gsField<> U4 = galerkinM.constructSolution(sol4);

        gsMatrix<real_t> S = F_boundary (mp, bases, U1, U2, U3, U4, bc1, bc2, bc3, bc4);
#ifdef _OPENMP
        t = omp_get_wtime() - t;
		std::cout << "Boundary treatment took" << t << "s" << std::endl;
        /// end: Construction of the boundary fluxes ----------------
	        
        /// begin: Performing a step solution -----------------------
        t = omp_get_wtime();
#endif
        gsMatrix<real_t> sol0 (sol);
        rhs_it1 = Ml*sol + tau*L*sol - tau*S;
        sol = solver.solve(rhs_it1);
#ifdef _OPENMP
        t = omp_get_wtime() - t;
		std::cout << "Solving the system took" << t << "s" << std::endl;
        /// end: Performing a step solution -------------------------
            

		t = omp_get_wtime();
#endif
        gsMatrix<real_t> u_dot_l (rhs_it1);
        gsMatrix<real_t> rhsloc (rhs_it1);
        rhsloc = L*sol-S;
        u_dot_l = solver.solve(rhsloc);
           
        /// begin: Performing flux correction -----------------------
        gsSparseMatrix<real_t> fMatrix (K);
        fMatrix = fMatrix + M;
        gsMatrix<real_t> Pp (rhs_it1);
        gsMatrix<real_t> Pm (rhs_it1);
        gsMatrix<real_t> Qp (rhs_it1);
        gsMatrix<real_t> Qm (rhs_it1);
        gsMatrix<real_t> pPp (rhs_it1);
        gsMatrix<real_t> pPm (rhs_it1);
        gsMatrix<real_t> pQp (rhs_it1);
        gsMatrix<real_t> pQm (rhs_it1);
        gsMatrix<real_t> ff (rhs_it1);
        #pragma omp parallel
		{
		#pragma omp for
        for (int k = 0; k < K.outerSize(); k++)
        {
            Pp.coeffRef(k) = 0.0;
            Pm.coeffRef(k) = 0.0;
            Qp.coeffRef(k) = 0.0;
            Qm.coeffRef(k) = 0.0;
            pPp.coeffRef(k) = 0.0;
            pPm.coeffRef(k) = 0.0;
            pQp.coeffRef(k) = 0.0;
            pQm.coeffRef(k) = 0.0;
            ff.coeffRef(k) = 0.0;
        }
		}
		#pragma omp parallel
		{
		#pragma omp for
        for (int e=0; e<K.outerSize(); ++e)
        {
            gsSparseMatrix<real_t>::InnerIterator itK(K,e);
            gsSparseMatrix<real_t>::InnerIterator itL(L,e);
            for (gsSparseMatrix<real_t>::InnerIterator it(fMatrix,e); it; ++it)
            {
                index_t k = it.row();
                index_t l = it.col();
                if (k/basisSize == l/basisSize)
                    it.valueRef() = galerkinM.matrix().at(k%basisSize,l%basisSize)*(u_dot_l.at(k)-u_dot_l.at(l))+ (itL.value()-itK.value()) * (sol.at(k)-sol.at(l)); //change
                else
                    it.valueRef() = 0.0;
                if (l != k && k < basisSize && l < basisSize)
                {
                    Pp.coeffRef(k) += math::max<real_t>(0.0,it.value());
                    Pm.coeffRef(k) += math::min<real_t>(0.0,it.value());

                    Qp.coeffRef(k) = math::max(Qp.coeffRef(k),sol.at(l)-sol.at(k));
                    Qm.coeffRef(k) = math::min(Qm.coeffRef(k),sol.at(l)-sol.at(k));
                }
                ++itK;
                ++itL;
            }
        }
		}
		#pragma omp parallel
		{
		#pragma omp for
        for (int e=0; e<K.outerSize(); ++e)
        {
            for (gsSparseMatrix<real_t>::InnerIterator it(fMatrix,e); it; ++it)
            {
                index_t k = it.row();
                index_t l = it.col();
                if (k != l && k < basisSize && l < basisSize)
                {
                    real_t f_p_kl = (gamma-1.0)*(fMatrix.at(k+3*basisSize,l+3*basisSize)+
                                                 (sol.at(k+basisSize)*sol.at(k+basisSize)+sol.at(k+2*basisSize)*sol.at(k+2*basisSize))/(2.0*sol.at(k)*sol.at(k))*fMatrix.at(k,l)
                                                 - (sol.at(k+basisSize)/sol.at(k)*fMatrix.at(k+basisSize,l+basisSize) + sol.at(k+2*basisSize)/sol.at(k)*fMatrix.at(k+2*basisSize,l+2*basisSize)));
                    real_t p_k = (gamma-1.0)*(sol.at(k+3*basisSize)-
                                              (sol.at(k+basisSize)*sol.at(k+basisSize)+sol.at(k+2*basisSize)*sol.at(k+2*basisSize))/(2.0*sol.at(k)));
                    real_t p_l = (gamma-1.0)*(sol.at(l+3*basisSize)-
                                              (sol.at(l+basisSize)*sol.at(l+basisSize)+sol.at(l+2*basisSize)*sol.at(l+2*basisSize))/(2.0*sol.at(l)));

                    pPp.coeffRef(k) += math::max<real_t>(0.0,f_p_kl);
                    pPm.coeffRef(k) += math::min<real_t>(0.0,f_p_kl);

                    pQp.coeffRef(k) = math::max(pQp.coeffRef(k),p_l-p_k);
                    pQm.coeffRef(k) = math::min(pQm.coeffRef(k),p_l-p_k);
                }
            }
        }
		}
        gsMatrix<real_t> Rp (sol);
        gsMatrix<real_t> Rm (sol);
        gsMatrix<real_t> pRp (sol);
        gsMatrix<real_t> pRm (sol);
        #pragma omp parallel
		{
		#pragma omp for
        for (int k = 0; k < K.outerSize(); k++)
        {
            if (tau*Pp.at(k)==0)
                Rp.coeffRef(k) = 1.0;
            else
                Rp.coeffRef(k) = math::min<real_t>(1.0,Ml_.at(k,k)*Qp.at(k)/(tau*Pp.at(k)));
            if (tau*Pm.at(k)==0)
                Rm.coeffRef(k) = 1.0;
            else
                Rm.coeffRef(k) = math::min<real_t>(1.0,Ml_.at(k,k)*Qm.at(k)/(tau*Pm.at(k)));

            if (tau*pPp.at(k)==0)
                pRp.coeffRef(k) = 1.0;
            else
                pRp.coeffRef(k) = math::min<real_t>(1.0,Ml_.at(k,k)*pQp.at(k)/(tau*pPp.at(k)));
            if (tau*pPm.at(k)==0)
                pRm.coeffRef(k) = 1.0;
            else
                pRm.coeffRef(k) = math::min<real_t>(1.0,Ml_.at(k,k)*pQm.at(k)/(tau*pPm.at(k)));
        }
		}
		#pragma omp parallel
		{
		#pragma omp for
        for (int e=0; e<fMatrix.outerSize(); ++e)
            for (gsSparseMatrix<real_t>::InnerIterator it(fMatrix,e); it; ++it)
            {
                int k = it.row();
                int l = it.col();
                real_t f_p_kl = (gamma-1.0)*(fMatrix.at(k%basisSize+3*basisSize,l%basisSize+3*basisSize)+
                                             (sol.at(k%basisSize+basisSize)*sol.at(k%basisSize+basisSize)+sol.at(k%basisSize+2*basisSize)*sol.at(k%basisSize+2*basisSize))/(2.0*sol.at(k%basisSize)*sol.at(k%basisSize))*fMatrix.at(k%basisSize,l%basisSize)
                                             - (sol.at(k%basisSize+basisSize)/sol.at(k%basisSize)*fMatrix.at(k%basisSize+basisSize,l%basisSize+basisSize) + sol.at(k%basisSize+2*basisSize)/sol.at(k%basisSize)*fMatrix.at(k%basisSize+2*basisSize,l%basisSize+2*basisSize)));
                real_t f_p_lk = (gamma-1.0)*(fMatrix.at(l%basisSize+3*basisSize,k%basisSize+3*basisSize)+
                                             (sol.at(l%basisSize+basisSize)*sol.at(l%basisSize+basisSize)+sol.at(l%basisSize+2*basisSize)*sol.at(l%basisSize+2*basisSize))/(2.0*sol.at(l%basisSize)*sol.at(l%basisSize))*fMatrix.at(l%basisSize,k%basisSize)
                                             - (sol.at(l%basisSize+basisSize)/sol.at(l%basisSize)*fMatrix.at(l%basisSize+basisSize,k%basisSize+basisSize) + sol.at(l%basisSize+2*basisSize)/sol.at(l%basisSize)*fMatrix.at(l%basisSize+2*basisSize,k%basisSize+2*basisSize)));
                real_t R_kl;
                real_t R_lk;
                real_t pR_kl;
                real_t pR_lk;
                if (fMatrix.at(k,l) >= 0.0)
                {
                    R_kl = Rp.at(k%basisSize);
                    R_lk = Rm.at(l%basisSize);
                }
                else
                {
                    R_kl = Rm.at(k%basisSize);
                    R_lk = Rp.at(l%basisSize);
                }
                if (f_p_kl >= 0.0 )
                    pR_kl = pRp.at(k%basisSize);
                else
                    pR_kl = pRm.at(k%basisSize);
                if (f_p_lk >= 0.0 )
                    pR_lk = pRp.at(l%basisSize);
                else
                    pR_lk = pRm.at(l%basisSize);
                ff.coeffRef(k) += math::min<real_t>(math::min<real_t>(R_kl,R_lk),
                                                    math::min<real_t>(pR_kl,pR_lk))*it.value();
            }
		}
        /// end: Performing flux correction -------------------------
#ifdef _OPENMP
	t = omp_get_wtime() - t;
	std::cout << "Flux correction took" << t << "s" << std::endl;
	t = omp_get_wtime();
#endif
        /// begin: Solving the corrected system ---------------------
        gsMatrix<real_t> rhs_loc (rhs_it1);
        rhs_loc = Ml*sol+tau*ff;
        sol = solver.solve(rhs_loc);
#ifdef  _OPENMP
    t = omp_get_wtime() - t;
	std::cout << "Solving the corrected system took" << t << "s" << std::endl;
#endif
        /// end: Solving the corrected system -----------------------
        /// begin: Step postprocessing ------------------------------
        gsInfo << i+1 << "/" << ntsteps << "\n";
        if (plot && i%plot_sparsity == 0)
        {
            for (int j=0;j<int(basisSize);j++)
            {
                sol1.coeffRef(j) = sol.at(j);
                sol2.coeffRef(j) = sol.at(j+basisSize)/sol.at(j);
                sol3.coeffRef(j) = sol.at(j+2*basisSize)/sol.at(j);
                sol4.coeffRef(j) = (gamma-1.0)*(sol.at(j+3*basisSize)-0.5*sol.at(j)*
                                                (sol.at(j+basisSize)/sol.at(j)*sol.at(j+basisSize)/sol.at(j)+
                                                 sol.at(j+2*basisSize)/sol.at(j)*sol.at(j+2*basisSize)/sol.at(j))); 
            }
            gsField<> U11 = galerkinM.constructSolution(sol1);
            gsField<> U21 = galerkinM.constructSolution(sol2);
            gsField<> U31 = galerkinM.constructSolution(sol3);
            gsField<> U41 = galerkinM.constructSolution(sol4);
            std::string fileName = "T_density" + util::to_string(std::ceil(i/10));
            gsWriteParaview<>( U11, fileName, nsamples);
            fileName = fileName + "0";
            collection.addTimestep(fileName,std::ceil(i/10),".vts");
            fileName = "T_velocity_x" + util::to_string(std::ceil(i/10));
            gsWriteParaview<>( U21, fileName, nsamples);
            fileName = fileName + "0";
            collection2.addTimestep(fileName,std::ceil(i/10),".vts");
            fileName = "T_velocity_y" + util::to_string(std::ceil(i/10));
            gsWriteParaview<>( U31, fileName, nsamples);
            fileName = fileName + "0";
            collection3.addTimestep(fileName,std::ceil(i/10),".vts");
            fileName = "T_pressure" + util::to_string(std::ceil(i/10));
            gsWriteParaview<>( U41, fileName, nsamples);
            fileName = fileName + "0";
            collection4.addTimestep(fileName,std::ceil(i/10),".vts");
        }
        /// end: Step postprocessing --------------------------------
    }
    /// end: Time stepping loop -----------------------------------------

    /// begin: Post-processing ------------------------------------------
    if (plot)
    {
        /// begin: Saving the collection of time steps ------------------
        collection.save();
        collection2.save();
        collection3.save();
        collection4.save();
        /// end: Saving the collection of time steps --------------------
    
        /// begin: Plotting the final solution --------------------------
        gsMatrix <real_t> sol1 (basisSize,1);
        gsMatrix <real_t> sol2 (basisSize,1);
        gsMatrix <real_t> sol3 (basisSize,1);
        gsMatrix <real_t> sol4 (basisSize,1);
        for (int j=0;j<int(basisSize);j++)
        {
            sol1.coeffRef(j) = sol.at(j);
            sol2.coeffRef(j) = sol.at(j+basisSize)/sol.at(j);
            sol3.coeffRef(j) = sol.at(j+2*basisSize)/sol.at(j);
            sol4.coeffRef(j) = (gamma-1.0)*(sol.at(j+3*basisSize)-0.5*sol.at(j)*
                                            (sol.at(j+basisSize)/sol.at(j)*sol.at(j+basisSize)/sol.at(j)+
                                             sol.at(j+2*basisSize)/sol.at(j)*sol.at(j+2*basisSize)/sol.at(j))); 
        }

        gsField<> U1 = galerkinM.constructSolution(sol1);
        gsField<> U2 = galerkinM.constructSolution(sol2);
        gsField<> U3 = galerkinM.constructSolution(sol3);
        gsField<> U4 = galerkinM.constructSolution(sol4);
        gsWriteParaview<>( U1, "Final_density", nsamples);
        gsWriteParaview<>( U2, "Final_velocity_x", nsamples);
        gsWriteParaview<>( U3, "Final_velocity_y", nsamples);
        gsWriteParaview<>( U4, "Final_pressure", nsamples);
        /// end: Plotting the final solution --------------------------
    }
    /// end: Prost-processing -------------------------------------------
	gsInfo << "DOF 1D" << sqrt(basisSize) << "\n";
	gsInfo << "basisSize" << basisSize << "\n";
	gsInfo << "bases.totalSize()" << bases.totalSize() << "\n";
	gsInfo << "bases_x" << bases[0].component(0).size() << "\n";
	gsInfo << "bases_y" << bases[0].component(1).size() << "\n";
    return 0;
}

///Function adding a block to the matrix
gsSparseMatrix<real_t> AddABlock(int sizei, int sizej, int ibegin,
                                 int jbegin, gsSparseMatrix<real_t> A){
    //Correct only for ColMajor Sparse Matrix!
    int i,j;
    std::vector<Eigen::Triplet<real_t> > tripletList;
    tripletList.reserve(A.nonZeros());
    for (int e=0; e<A.outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(A,e); it; ++it)
        {
            i = it.row() + ibegin;
            j = it.col() + jbegin;

            tripletList.push_back(Eigen::Triplet<real_t> (i,j,it.value()));
        }

    gsSparseMatrix<real_t> res(sizei,sizej);
    res.setFromTriplets(tripletList.begin(), tripletList.end());
    return res;
}

/// Function computing the x component of the Jacobian tensor
gsMatrix<real_t> A_jx(gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,4);
    real_t gamma = 1.4;
    real_t gamma1 = (gamma-1.0)/2.0;
    real_t gamma2 = (gamma-3.0)/2.0;
    real_t v1 = u.at(1)/u.at(0);
    real_t v2 = u.at(2)/u.at(0);
    real_t H = (gamma*u.at(3) + (1.0-gamma)*(u.at(1)*u.at(1)+u.at(2)*u.at(2))/
                (2.0*u.at(0)))/u.at(0);
    res << 0, 1.0, 0, 0,
        gamma2*v1*v1+gamma1*v2*v2, (3.0-gamma)*v1, (1.0-gamma)*v2, gamma-1.0,
        -v1*v2, v2, v1, 0,
        gamma1*(v1*v1*v1+v2*v2*v1)-H*v1, H-(gamma-1.0)*v1*v1, (1.0-gamma)*v1*v2, gamma*v1;
    return res;
}

/// Function computing the y component of the Jacobian tensor
gsMatrix<real_t> A_jy(gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,4);
    real_t gamma = 1.4;
    real_t gamma1 = (gamma-1.0)/2.0;
    real_t gamma2 = (gamma-3.0)/2.0;
    real_t v1 = u.at(1)/u.at(0);
    real_t v2 = u.at(2)/u.at(0);
    real_t H = (gamma*u.at(3) + (1.0-gamma)*(u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)))/
        u.at(0);
    res << 0, 0, 1, 0,
        -v1*v2, v2, v1, 0,
        gamma2*v2*v2+gamma1*v1*v1, (1.0-gamma)*v1, (3.0-gamma)*v2, gamma-1.0,
        gamma1*(v1*v1*v2+v2*v2*v2)-H*v2, (1.0-gamma)*v1*v2, H-(gamma-1.0)*v2*v2, gamma*v2;
    return res;
}

/// Function computing the x component of the inviscid flux
gsMatrix<real_t> F_x(gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    real_t p = (gamma-1.0)*(u.at(3) - (u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)));
    res << u.at(1),
        u.at(1)*u.at(1)/u.at(0) + p,
        u.at(1)*u.at(2)/u.at(0),
        u.at(3)*u.at(1)/u.at(0) + p*u.at(1)/u.at(0);
    return res;
}

/// Function computing the y component of the inviscid flux
gsMatrix<real_t> F_y(gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    real_t p = (gamma-1.0)*(u.at(3) - (u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)));
    res << u.at(2),
        u.at(1)*u.at(2)/u.at(0),
        u.at(2)*u.at(2)/u.at(0)+p,
        u.at(3)*u.at(2)/u.at(0) + p*u.at(2)/u.at(0);
    return res;
}

/// Function transforming state variables to Riemann invariants
gsMatrix<real_t> UtoW (gsMatrix<real_t> n, gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    real_t c_v = 7.0/2.0;
    real_t gamma = 1.4;
    real_t p = (gamma-1.0)*(u.at(3) - (u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)));
    real_t vn = n.at(0)*u.at(1)/u.at(0) + n.at(1)*u.at(2)/u.at(0);
    real_t vxi = n.at(1)*u.at(1)/u.at(0) - n.at(0)*u.at(2)/u.at(0);
    real_t c = math::sqrt(gamma*p/u.at(0));
    res << vn - (2.0*c)/(gamma-1.0),
        c_v * math::log(p / math::pow(u.at(0), gamma)),
        vxi,
        vn + (2.0*c)/(gamma-1.0);
    return res;
}

/// Function transforming Riemann invariants to state variables
gsMatrix<real_t> WtoU (gsMatrix<real_t> n, gsMatrix<real_t> w){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    real_t c_v = 7.0/2.0;
    real_t vn = (w.at(3)+w.at(0))/2.0;
    real_t vxi = w.at(2);
    real_t c = (gamma-1.0)/4.0*(w.at(3)-w.at(0));
    real_t rho = math::pow(c*c/gamma*exp(-w.at(1)/c_v),1.0/(gamma-1.0));
    real_t p = rho*c*c/gamma;
    res << rho,
        rho*(vn*n.at(0)+vxi*n.at(1)),
        rho*(vn*n.at(1)-vxi*n.at(0)),
        p/(gamma-1.0)+rho/2.0*(vn*vn+vxi*vxi);
    return res;
}

/// Function computing the result of the scalar product of normal vector and the Jacobi matrix
gsMatrix<real_t> nA (gsMatrix<real_t> n, gsMatrix<real_t> u, gsMatrix<real_t> u_inf){
    gsMatrix<real_t> res (4,4);
    real_t gamma = 1.4;
    real_t v_ix = u.at(1)/u.at(0);
    real_t v_iy = u.at(2)/u.at(0);
    real_t v_jx = u_inf.at(1)/u_inf.at(0);
    real_t v_jy = u_inf.at(2)/u_inf.at(0);
    real_t v_ijx = (math::sqrt(u.at(0))*v_ix + math::sqrt(u_inf.at(0))*v_jx)/
        (math::sqrt(u.at(0)) + math::sqrt(u_inf.at(0)));
    real_t v_ijy = (math::sqrt(u.at(0))*v_iy + math::sqrt(u_inf.at(0))*v_jy)/
        (math::sqrt(u.at(0)) + math::sqrt(u_inf.at(0)));
    real_t v_ijn = n.at(0)*v_ijx + n.at(1)*v_ijy;
    real_t H_i = (gamma*u.at(3) + (1.0-gamma)*(u.at(1)*u.at(1)+u.at(2)*u.at(2))/
                  (2.0*u.at(0)))/u.at(0);
    real_t H_j = (gamma*u_inf.at(3) + (1.0-gamma)*
                  (u_inf.at(1)*u_inf.at(1)+u_inf.at(2)*u_inf.at(2))/(2*u_inf.at(0)))/u_inf.at(0);
    real_t H_ij = (math::sqrt(u.at(0))*H_i + math::sqrt(u_inf.at(0))*H_j)/
        (math::sqrt(u.at(0)) + math::sqrt(u_inf.at(0)));
    real_t c_ij = math::sqrt((gamma-1.0)*(H_ij-(v_ijx*v_ijx+v_ijy*v_ijy)/2.0));
    real_t q = math::sqrt(v_ijx*v_ijx+v_ijy*v_ijy)/2.0;
    real_t b2 = (gamma-1.0)/(c_ij*c_ij);
    real_t b1 = b2*q;
    gsMatrix<real_t> lambda (4,4);
    lambda << math::abs(v_ijn - c_ij), 0, 0, 0,
        0, math::abs(v_ijn), 0, 0,
        0, 0, math::abs(v_ijn + c_ij), 0,
        0, 0, 0, math::abs(v_ijn);
    gsMatrix<real_t> R(4,4);
    R << 1.0, 1.0, 1.0, 0,
        v_ijx-c_ij*n.at(0), v_ijx, v_ijx+c_ij*n.at(0), n.at(1),
        v_ijy-c_ij*n.at(1), v_ijy, v_ijy+c_ij*n.at(1), -n.at(0),
        H_ij-c_ij*v_ijn, q, H_ij+c_ij*v_ijn, v_ijx*n.at(1)-v_ijy*n.at(0);
    gsMatrix<real_t> L(4,4);
    L << 0.5*(b1+v_ijn/c_ij), 0.5*(-b2*v_ijx-n.at(0)/c_ij), 0.5*(-b2*v_ijy - n.at(1)/c_ij), 0.5*b2,
        1.0-b1, b2*v_ijx, b2*v_ijy, -b2,
        0.5*(b1 - v_ijn/c_ij), 0.5*(-b2*v_ijx+n.at(0)/c_ij), 0.5*(-b2*v_ijy+ n.at(1)/c_ij), 0.5*b2,
        n.at(0)*v_ijn-n.at(1)*v_ijx, n.at(1), -n.at(0), 0;
    res = R*lambda*L;
    return res;
}

/// Neumann boundary flux evaluation
gsMatrix<real_t> F_n_ (gsMatrix<real_t> n, gsMatrix<real_t> u, gsMatrix<real_t> u_inf){
    gsMatrix<real_t> res (4,1);
    res = n.at(0)*(F_x(u)+F_x(u_inf))/2.0+n.at(1)*(F_y(u)+F_y(u_inf))/2.0
        - 0.5*nA(n,u,u_inf)*(u_inf-u);
    return res;
}

/// Fluxes vector for wall boundary condition
gsMatrix<real_t> F_n_wall (gsMatrix<real_t> n, gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    gsMatrix<real_t> u_inf (4,1);
    real_t vx = u.at(1)/u.at(0);
    real_t vy = u.at(2)/u.at(0);
    u_inf << u.at(0),
        u.at(0)*(vx - 2*n.at(0)*(n.at(0)*vx+n.at(1)*vy)),
        u.at(0)*(vy - 2*n.at(1)*(n.at(0)*vx+n.at(1)*vy)),
        u.at(3);
    res = F_n_(n,u,u_inf);
    return res;
}

/// Fluxes vector for subsonic outlet
gsMatrix<real_t> F_n_subin (gsMatrix<real_t> n, gsMatrix<real_t> u,
                            real_t rho, real_t p, real_t v_tangential){
    real_t gamma = 1.4;
    real_t c_v = 7.0/2.0;
    gsMatrix<real_t> res (4,1);
    gsMatrix<real_t> w = UtoW(n,u);
    w.coeffRef(2) = v_tangential;
    w.coeffRef(1) = c_v*math::log(p/math::pow(rho,gamma));
    w.coeffRef(0) = w.at(3)-4.0/(gamma-1.0)*math::sqrt(gamma*p/rho);
    gsMatrix<real_t> u_inf = WtoU(n,w);
    res = F_n_(n,u,u_inf);
    return res;
}

/// Fluxes vector for subsonic outlet
gsMatrix<real_t> F_n_subout (gsMatrix<real_t> n, gsMatrix<real_t> u, real_t p_out){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    real_t p = (gamma-1.0)*(u.at(3) - (u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)));
    gsMatrix<real_t> w = UtoW(n,u);
    w.coeffRef(0) = w.at(3) - 4.0/(gamma-1)*math::sqrt(gamma*p_out/u.at(0)*math::pow(p/p_out,1/gamma)) ;
    gsMatrix<real_t> u_inf = WtoU(n,w);
    res = F_n_(n,u,u_inf);
    return res;
}

/// Fluxes vector for supersonic inlet
gsMatrix<real_t> F_n_supin (gsMatrix<real_t> n, real_t rho, real_t v_x,
                            real_t v_y, real_t p_in){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    gsMatrix<real_t> u_inf (4,1);
    real_t rhoE = p_in/(gamma-1.0)+rho*(v_x*v_x+v_y*v_y)/2.0;
    u_inf << rho,
        rho*v_x,
        rho*v_y,
        rhoE;
    res = F_n_(n,u_inf,u_inf);
    return res;
}

/// Fluxes vector for supersonic outlet
gsMatrix<real_t> F_n_supout (gsMatrix<real_t> n, gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    res = F_n_(n,u,u);
    return res;
}
// "wall"
// "subin rho p v_tangential"
// "subout p" 
// "supin rho v_x v_y p"
// "supout"
// "periodic"


/// Wrapper for boundary condition types
gsMatrix<real_t> F_n_type (gsMatrix<real_t> n, gsMatrix<real_t> u, std::string bc){
    gsMatrix<real_t> res (4,1);
    if (bc.substr(0,4) == "wall")
        res = F_n_wall(n,u);
    else if (bc.substr(0,5) == "subin")
    {
        bc = bc.substr(bc.find(" "));
        std::string::size_type pos; 
        real_t rho = std::stod(bc,&pos);
        bc = bc.substr(pos+1);
        real_t p = std::stod(bc,&pos);
        bc = bc.substr(pos+1);
        real_t v_tan = std::stod(bc,&pos);
        res = F_n_subin(n,u, rho, p, v_tan);
    }
    else if (bc.substr(0,6) == "subout")
    {
        bc = bc.substr(bc.find(" "));
        std::string::size_type pos;
        real_t p = std::stod(bc,&pos);
        res = F_n_subout(n,u, p);
    }
    else if (bc.substr(0,5) == "supin")
    {
        bc = bc.substr(bc.find(" "));
        std::string::size_type pos;
        real_t rho = std::stod(bc,&pos);
        bc = bc.substr(pos+1);
        real_t v_x = std::stod(bc,&pos);
        bc = bc.substr(pos+1);
        real_t v_y = std::stod(bc,&pos);
        bc = bc.substr(pos+1);
        real_t p = std::stod(bc,&pos);
        res = F_n_supin(n,rho, v_x, v_y, p);
    }
    else if (bc.substr(0,6) == "supout")
        res = F_n_supout(n,u);
    else if (bc.substr(0,4) == "periodic")
        std::cout << "TODO" << std::endl;
    else
        std::cout << "ERROR: Not implemented type of boundary condition - see specification." << std::endl;
    return res;
}

/// Function building the boundary conditions
gsMatrix<real_t> F_boundary (gsMultiPatch<real_t> mp, gsMultiBasis<> bases, gsField<> & U1,
                             gsField<> & U2, gsField<> & U3, gsField<> & U4, std::string bc1,
                             std::string bc2, std::string bc3, std::string bc4)
{
    gsMatrix <real_t> res (4*bases.totalSize(),1);
    for (int k = 0; k < 4*int(bases.totalSize()); k++)
        res(k) = 0.0;
    int p = bases.degree();
    gsMatrix <real_t> x (14,1);
    x << -0.9862838086968121,
        -0.9284348836635741,
        -0.8272013150697653,
        -0.6872929048116853,
        -0.5152486363581548,
        -0.3191123689278905,
        -0.1080549487073448,
        0.108054948707344,
        0.319112368927890,
        0.515248636358154,
        0.687292904811686,
        0.827201315069765,
        0.928434883663574,
        0.986283808696812 ;
    gsMatrix <real_t> w (14,1);
    w << 0.035119460331751,
        0.080158087159761,
        0.121518570687903,
        0.157203167158193,
        0.185538397477938,
        0.205198463721295,
        0.215263853463158,
        0.215263853463158,
        0.205198463721296,
        0.185538397477938,
        0.157203167158194,
        0.121518570687903,
        0.080158087159761,
        0.035119460331751;
    //west boundary - subsonic inlet
    {
        //real_t rho = 1.0;
        //real_t press = 10.0;
        //real_t v_tangential = 0.0;
        gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(1);
        #pragma omp parallel
		{
		#pragma omp for
        for (int i=0; i< bd_basis->size(); i++)              //Disabled: Not repeating corners
        {
            gsMatrix <real_t> interval = bd_basis->function(i).support();
            real_t subint = (interval.at(1) - interval.at(0))/(p+1.0);
            gsMatrix <real_t> pres (4,1);
            pres << 0.0, 0.0, 0.0, 0.0;
            for (int k=0; k<p+1;k++)
            {
                gsMatrix <real_t> pres0 (4,1);
                pres0 << 0.0, 0.0, 0.0, 0.0;
                gsMatrix <real_t> quad (2,14);
                for (int l=0; l<14;l++)
                {
                    quad.coeffRef(1,l) = ((interval.at(0)+subint*k)*(1.0-x.at(l))+
                                          (interval.at(0)+subint*(k+1.0))*(1.0+x.at(l)))/2.0; //Maybe an error here
                    quad.coeffRef(0,l) = 0.0;
                }
                gsMatrix <real_t> u (4,14);
                gsMatrix <real_t> u_val;
                u_val = U1.value(quad);
                u.row(0) = u_val.row(0);
                u_val = U2.value(quad);
                u.row(1) = u_val.row(0);
                u_val = U3.value(quad);
                u.row(2) = u_val.row(0);
                u_val = U4.value(quad);
                u.row(3) = u_val.row(0);
                gsMatrix <real_t> normal (2,1);
                normal = outer_normal(quad, mp, 1);
                gsMatrix <real_t> funct = bd_basis->function(i).eval(quad.row(1));
                gsMatrix <real_t> F_n (4,14);
                quad = mapping(quad, mp);
                real_t len_of = 0.0;
                for (int l=0; l<14; l++)
                {
                    F_n.col(l) = F_n_type(normal.col(l),u.col(l),bc1);
                }
                for (int l=0; l<14; l++)
                {
                    pres0 = pres0 + w.at(l)*funct.at(l)*F_n.col(l)*0.5;
                }
                for (int l=0; l<13; l++)
                {
                    len_of = len_of + math::sqrt((quad.coeffRef(0,l+1)-quad.coeffRef(0,l))*(quad.coeffRef(0,l+1)-quad.coeffRef(0,l))
                                                 + (quad.coeffRef(1,l+1)-quad.coeffRef(1,l))*(quad.coeffRef(1,l+1)-quad.coeffRef(1,l)));
                }
                pres = pres + pres0 * len_of;
            }
            for (int k=0; k<4; k++)
            {
                res.coeffRef(i*bases[0].component(0).size()+k*bases.totalSize()) += pres.at(k);
            }
        }
	}
    }
    //east boundary 
    {
        //real_t press = 1.0;
        gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(2);
        #pragma omp parallel
		{
		#pragma omp for
        for (int i=0; i< bd_basis->size(); i++)          //Disabled: Not repeating corners
        {
            gsMatrix <real_t> interval = bd_basis->function(i).support();
            real_t subint = (interval.at(1) - interval.at(0))/(p+1.0);
            gsMatrix <real_t> pres (4,1);
            pres << 0.0, 0.0, 0.0, 0.0;
            for (int k=0; k<p+1;k++)
            {
                gsMatrix <real_t> pres0 (4,1);
                pres0 << 0.0, 0.0, 0.0, 0.0;
                gsMatrix <real_t> quad (2,14);
                for (int l=0; l<14;l++)
                {
                    quad.coeffRef(1,l) = ((interval.at(0)+subint*k)*(1.0-x.at(l))+
                                          (interval.at(0)+subint*(k+1))*(1.0+x.at(l)))/2.0; //Maybe an error here
                    quad.coeffRef(0,l) = 1.0;
                }
                gsMatrix <real_t> quad_mapped = mp.patch(0).eval(quad);
                gsMatrix <real_t> u (4,14);
                gsMatrix <real_t> u_val;
                u_val = U1.value(quad);
                u.row(0) = u_val.row(0);
                u_val = U2.value(quad);
                u.row(1) = u_val.row(0);
                u_val = U3.value(quad);
                u.row(2) = u_val.row(0);
                u_val = U4.value(quad);
                u.row(3) = u_val.row(0);
                gsMatrix <real_t> normal (2,1);
                normal = outer_normal(quad, mp, 2);
                gsMatrix <real_t> funct = bd_basis->function(i).eval(quad.row(1)); // this was not working
                gsMatrix <real_t> F_n (4,14);
                quad = mapping(quad, mp);
                real_t len_of = 0.0;
                for (int l=0; l<14; l++)
                {
                    F_n.col(l) = F_n_type(normal.col(l),u.col(l),bc2);
                }
                for (int l=0; l<14; l++)
                {
                    pres0 = pres0 + w.at(l)*funct.at(l)*F_n.col(l)*0.5;
                }
                for (int l=0; l<13; l++)
                {
                    len_of = len_of + math::sqrt((quad.coeffRef(0,l+1)-quad.coeffRef(0,l))*(quad.coeffRef(0,l+1)-quad.coeffRef(0,l))
                                                 + (quad.coeffRef(1,l+1)-quad.coeffRef(1,l))*(quad.coeffRef(1,l+1)-quad.coeffRef(1,l)));
                }
                pres = pres + pres0 * len_of;

            }
            for (int k=0; k<4; k++)
            {
                res.coeffRef((i+1)*bases[0].component(0).size()-1+k*bases.totalSize()) += pres.at(k);
            }
        }
	}
    }
    //south boundary 
    {
        gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(3);
        #pragma omp parallel
		{
		#pragma omp for
        for (int i=0; i< bd_basis->size(); i++)
        {
            gsMatrix <real_t> interval = bd_basis->function(i).support();
            real_t subint = (interval.at(1) - interval.at(0))/(p+1.0);
            gsMatrix <real_t> pres (4,1);
            pres << 0.0, 0.0, 0.0, 0.0;
            for (int k=0; k<p+1;k++)
            {
                gsMatrix <real_t> pres0 (4,1);
                pres0 << 0.0, 0.0, 0.0, 0.0;
                gsMatrix <real_t> quad (2,14);
                for (int l=0; l<14;l++)
                {
                    quad.coeffRef(0,l) = ((interval.at(0)+subint*k)*(1.0-x.at(l))+
                                          (interval.at(0)+subint*(k+1))*(1.0+x.at(l)))/2.0; //Maybe an error here
                    quad.coeffRef(1,l) = 0.0;
                }
                gsMatrix <real_t> u (4,14);
                gsMatrix <real_t> u_val;
                u_val = U1.value(quad);      
                u.row(0) = u_val.row(0);
                u_val = U2.value(quad);
                u.row(1) = u_val.row(0);
                u_val = U3.value(quad);
                u.row(2) = u_val.row(0);
                u_val = U4.value(quad);
                u.row(3) = u_val.row(0);
                gsMatrix <real_t> normal (2,1);
                normal = outer_normal(quad, mp, 3);
                gsMatrix <real_t> funct = bd_basis->function(i).eval(quad.row(0));
                gsMatrix <real_t> F_n (4,14);
                quad = mapping(quad, mp);
                real_t len_of = 0.0;
                for (int l=0; l<14; l++)
                {
                    F_n.col(l) = F_n_type(normal.col(l),u.col(l),bc3);
                }
                for (int l=0; l<14; l++)
                {
                    pres0 = pres0 + w.at(l)*funct.at(l)*F_n.col(l)*0.5;
                }
                for (int l=0; l<13; l++)
                {
                    len_of = len_of + math::sqrt((quad.coeffRef(0,l+1)-quad.coeffRef(0,l))*(quad.coeffRef(0,l+1)-quad.coeffRef(0,l))
                                                 + (quad.coeffRef(1,l+1)-quad.coeffRef(1,l))*(quad.coeffRef(1,l+1)-quad.coeffRef(1,l)));
                }
                pres = pres + pres0 * len_of;
            }
            for (int k=0; k<4; k++)
            {
                res.coeffRef(i+k*bases.totalSize()) += pres.at(k);
            }
        }
	}
    }

    //north boundary 
    {
        gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(4);
        #pragma omp parallel
		{
		#pragma omp for
        for (int i=0; i< bd_basis->size(); i++)
        {
            gsMatrix <real_t> interval = bd_basis->function(i).support();
            real_t subint = (interval.at(1) - interval.at(0))/(p+1.0);
            gsMatrix <real_t> pres (4,1);
            pres << 0.0, 0.0, 0.0, 0.0;
            for (int k=0; k<p+1;k++)
            {
                gsMatrix <real_t> pres0 (4,1);
                pres0 << 0.0, 0.0, 0.0, 0.0;
                gsMatrix <real_t> quad (2,14);
                for (int l=0; l<14;l++)
                {
                    quad.coeffRef(0,l) = ((interval.at(0)+subint*k)*(1.0-x.at(l))+
                                          (interval.at(0)+subint*(k+1))*(1.0+x.at(l)))/2.0; //Maybe an error here
                    quad.coeffRef(1,l) = 1.0;
                }
                gsMatrix <real_t> u (4,14);
                gsMatrix <real_t> u_val;
                u_val = U1.value(quad);
                u.row(0) = u_val.row(0);
                u_val = U2.value(quad);
                u.row(1) = u_val.row(0);
                u_val = U3.value(quad);
                u.row(2) = u_val.row(0);
                u_val = U4.value(quad);
                u.row(3) = u_val.row(0);

                gsMatrix <real_t> normal (2,1);
                normal = outer_normal(quad, mp, 4);
                gsMatrix <real_t> funct = bd_basis->function(i).eval(quad.row(0));
                gsMatrix <real_t> F_n (4,14);
                quad = mapping(quad, mp);
                real_t len_of = 0.0;
                for (int l=0; l<14; l++)
                {
                    F_n.col(l) = F_n_type(normal.col(l),u.col(l),bc4);
                }
                for (int l=0; l<14; l++)
                {
                    pres0 = pres0 + w.at(l)*funct.at(l)*F_n.col(l)*0.5;
                }
                for (int l=0; l<13; l++)
                {
                    len_of = len_of + math::sqrt((quad.coeffRef(0,l+1)-quad.coeffRef(0,l))*(quad.coeffRef(0,l+1)-quad.coeffRef(0,l))
                                                 + (quad.coeffRef(1,l+1)-quad.coeffRef(1,l))*(quad.coeffRef(1,l+1)-quad.coeffRef(1,l)));
                }
                pres = pres + pres0 * len_of;
            }
            for (int k=0; k<4; k++)
            {
                res.coeffRef(bases.totalSize()-bd_basis->size()+i+k*bases.totalSize()) += pres.at(k);
            }
        }
	}
    }
    return res;
}

// Function computing initial condition with constrained L2 projection
gsMatrix<real_t> ConstL2(gsMatrix<real_t>        INrhs,
                         gsSparseMatrix<real_t>  Mc,
                         gsSparseMatrix <real_t> Ml)
{
    gsSparseSolver<>::LU solver;
    gsMatrix<real_t> u_l (INrhs);
    gsMatrix<real_t> u_h (INrhs);
    solver.compute(Mc);
    u_h = solver.solve(INrhs);
    solver.compute(Ml);
    u_l = solver.solve(INrhs);
    gsMatrix<real_t> u_init (u_l);
    gsMatrix<real_t> Pp (INrhs);
    gsMatrix<real_t> Pm (INrhs);
    gsMatrix<real_t> Qp (INrhs);
    gsMatrix<real_t> Qm (INrhs);
    gsMatrix<real_t> ff (INrhs);
    for (int k = 0; k < Mc.outerSize(); k++)
    {
        Pp.coeffRef(k) = 0.0;
        Pm.coeffRef(k) = 0.0;
        Qp.coeffRef(k) = 0.0;
        Qm.coeffRef(k) = 0.0;
        ff.coeffRef(k) = 0.0;
    }
    gsSparseMatrix<real_t> f (Mc);
    #pragma omp parallel
	{
	#pragma omp for
    for (int e=0; e<f.outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(f,e); it; ++it)
        {
            int k = it.row();
            int l = it.col();
            if ((u_h(k)-u_h(l))*(u_l(k)-u_l(l)) > 0)
                it.valueRef() = Mc.at(k,l)*(u_h.at(k)-u_h.at(l));
            else
                it.valueRef() = 0.0;
            if (l!= k)
            {
                Pp.coeffRef(k) += math::max<real_t>(0.0,it.value());
                Pm.coeffRef(k) += math::min<real_t>(0.0,it.value());
                Qp.coeffRef(k)  = math::max<real_t>(Qp.coeffRef(k),u_l(l)-u_l(k));
                Qm.coeffRef(k)  = math::min<real_t>(Qm.coeffRef(k),u_l(l)-u_l(k));
            }
        }
    }
    gsMatrix<real_t> Rp (INrhs);
    gsMatrix<real_t> Rm (INrhs);
    for (int k = 0; k < f.outerSize(); k++)
    {
        if (Pp.at(k)==0)
            Rp.coeffRef(k) = 1.0;
        else
            Rp.coeffRef(k) = math::min<real_t>(1.0,Ml.at(k,k)*Qp.at(k)/(Pp.at(k)));
        if (Pm.at(k)==0)
            Rm.coeffRef(k) = 1.0;
        else
            Rm.coeffRef(k) = math::min<real_t>(1.0,Ml.at(k,k)*Qm.at(k)/(Pm.at(k)));
    }
    #pragma omp parallel
	{
	#pragma omp for    
    for (int e=0; e<f.outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(f,e); it; ++it)
        {
            int k = it.row();
            int l = it.col();
            if (f.at(k,l) > 0.0)
            {
                ff.coeffRef(k) += math::min(Rp.at(k),Rm.at(l))*it.value();
            }
            else
            {
                ff.coeffRef(k) += math::min(Rm.at(k),Rp.at(l))*it.value();
            }
        }
	}
    for (int k=0; k<f.outerSize();k++)
        ff.coeffRef(k) = ff.at(k)/Ml.at(k,k);
    u_init = u_l + ff;
    return u_init;
}

/// Function returning the outer normal vector
gsMatrix<real_t> outer_normal(gsMatrix<real_t> points, gsMultiPatch<real_t> mp, boxSide a){
    gsMapData<real_t> data;
    data.addFlags(NEED_OUTER_NORMAL);
    data.side = a;
    data.points = points;
    mp.patch(0).computeMap(data);
    data.outNormals.colwise().normalize();
    return data.outNormals;
}

/// Function returning mapping of parametric domain points to physical domain
gsMatrix<real_t> mapping(gsMatrix<real_t> points, gsMultiPatch<real_t> mp){
    gsMapData<real_t> data;
    data.addFlags(NEED_VALUE);
    data.points = points;
    mp.patch(0).computeMap(data);
    return data.values.at(0);
}
