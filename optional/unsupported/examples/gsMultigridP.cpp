#include <gismo.h>
#include <string>
#include <gsAssembler/gsNormL2.h>
#include <gsUtils/gsNorms.h>

using namespace std;
using namespace gismo;

/** @brief The p-multigrid base class provides the basic
 *  methods (smoothing, prolongation, restriction) for
 *  implementing p-multigrid methods
 */

template<class T>
struct pMultigridBase
{

public:

  /// @brief Apply p-multigrid solver to given right-hand side on level l
  virtual void solve(const gsMatrix<T> & rhs, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis , gsMatrix<T>& x, const int& numLevels, const int& numSmoothing, const int& typeCycle, const int& 		typeSolver, const int& typeBCHandling, gsBoundaryConditions<T> bcInfo, gsGeometry<>::Ptr geo, const int& lumping)
  {
  	if( numLevels == 1)
  	{
  		solvecoarse(rhs, x, numLevels);
  		return;
  	}
  	else
  	{
        	gsMatrix<T> fineRes;
        	gsMatrix<T> coarseRes;
        	gsMatrix<T> fineCorr;
        	gsMatrix<T> coarseCorr;

    		presmoothing(rhs, x, numLevels, numSmoothing, fineRes);
        	restriction(fineRes, coarseRes, numLevels, m_basis, lumping, typeBCHandling, bcInfo, geo);
                //coarseRes.setZero(coarseRes.rows(),1);
             	coarseCorr.setZero(coarseRes.rows(),1);
             	for( int j = 0 ; j < (typeCycle == 2 ? 2 : 1) ; j++)
             	{
             		solve(coarseRes, m_basis, coarseCorr, numLevels-1, numSmoothing, typeCycle, typeSolver, typeBCHandling, bcInfo, geo, lumping);
             	}
		prolongation(coarseCorr, fineCorr, numLevels, m_basis, lumping, typeBCHandling, bcInfo, geo);
    		postsmoothing(rhs,x, numLevels, numSmoothing, fineCorr, typeSolver);
    	}
  }

  /// @brief Apply fixed number of smoothing steps (pure virtual method)
  virtual void presmoothing(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels, const int& numSmoothing, gsMatrix<T> & fineRes ) = 0;

  /// @brief Apply fixed number of smoothing steps (pure virtual method)
  virtual void postsmoothing(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels, const int& numSmoothing, gsMatrix<T> & fineCorr, const int& typeSolver) = 0;

  /// @brief Apply coarse solver (pure virtual method)
  virtual void solvecoarse(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels) = 0;

  /// @brief Apply Galerkin Condition (pure virtual method)
  virtual void galerkincondition(gsSparseMatrix<T> & Acoarse,
  	                       const gsSparseMatrix<T> & Afine,
  	                       const int               & numLevels,
  	 vector<memory::shared_ptr<gsMultiBasis<T> > >   m_basis,
  	                       //const int               & lumping,
  	                       const int               & typeBCHandling,
  	                       gsBoundaryConditions<T>   bcInfo,
  	                       gsGeometry<>::Ptr         geo) = 0;

  /// @brief Prolongate coarse space function to fine space
  virtual void prolongation(const gsMatrix<T>& Xcoarse, gsMatrix<T>& Xfine, const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& lumping, const int& typeBCHandling,    	gsBoundaryConditions<T> bcInfo, gsGeometry<>::Ptr geo)
  {
  	// Define the low and high order basis
    	gsMultiBasis<> basisL = *m_basis[numLevels-2];
    	gsMultiBasis<> basisH = *m_basis[numLevels-1];

     	// Determine matrix P (high_order * low_order)
    	gsMultiPatch<> mp(*geo);
    	typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    	gsExprAssembler<real_t> ex(1,1);
    	geometryMap G = ex.getMap(mp);
    	typedef gsExprAssembler<real_t>::variable variable;
    	typedef gsExprAssembler<real_t>::space    space;
    	space v_n = ex.getSpace(basisH ,1, 0);
    	space u_n = ex.getTestSpace(v_n , basisL);
    	if(typeBCHandling == 1)
    	{
    		v_n.addBc(bcInfo.get("Dirichlet"));
    		u_n.addBc(bcInfo.get("Dirichlet"));
    	}
    	ex.setIntegrationElements(basisH);
    	ex.initSystem();
    	ex.assemble(u_n*meas(G) * v_n.tr()); //*meas(G));
    	gsSparseMatrix<> P = ex.matrix().transpose();

    	// Determine matrix M (high_order * high_order)
    	gsExprAssembler<real_t> ex2(1,1);
    	geometryMap G2 = ex2.getMap(mp);
    	space w_n = ex2.getSpace(basisH ,1, 0);
    	if(typeBCHandling == 1)
    	{
    		w_n.addBc(bcInfo.get("Dirichlet"));
    	}
    	ex2.setIntegrationElements(basisH);
    	ex2.initSystem();
    	ex2.assemble(w_n*meas(G2) * w_n.tr()); //*meas(G2));
    	gsSparseMatrix<> M = ex2.matrix();

    	// Prolongate Xcoarse to Xfine
    	gsVector<> temp = P*Xcoarse;

	if (lumping == 1)
   	{
        	// Lump the mass matrix (should be a gsSparseMatrix!)
        	gsVector<> one = gsVector<>::Ones(M.rows());
        	gsMatrix<> M_L = (M*one).asDiagonal();
        	//gsConjugateGradient<> CGSolver(M_L);
    		//CGSolver.solve(temp,Xfine);

		// Alternative approach (no solve necessary)
                gsMatrix<> M_L_inv_d = (M*one).array().inverse();
		gsMatrix<> M_L_inv = (M_L_inv_d).asDiagonal();
                Xfine = M_L_inv*temp;


    	}
    	else
    	{
    		gsConjugateGradient<> CGSolver(M);
    		CGSolver.solve(temp,Xfine);
    	}
  }

  /// @brief Restrict fine space function to coarse space
  virtual void restriction(const gsMatrix<T>& Xfine, gsMatrix<T>& Xcoarse, const int& numLevels, vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis, const int& lumping, const int& typeBCHandling, gsBoundaryConditions<T> bcInfo, gsGeometry<>::Ptr geo)
  {
    	// Define the low and high order basis
    	gsMultiBasis<> basisL = *m_basis[numLevels-2];
    	gsMultiBasis<> basisH = *m_basis[numLevels-1];

    	// Determine matrix P (high_order * low_order)
    	gsExprAssembler<real_t> ex(1,1);
    	gsMultiPatch<> mp(*geo);
    	typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    	geometryMap G = ex.getMap(mp);

    	typedef gsExprAssembler<real_t>::variable variable;
    	typedef gsExprAssembler<real_t>::space    space;
    	space v_n = ex.getSpace(basisH ,1, 0);
    	space u_n = ex.getTestSpace(v_n , basisL);
    	if( typeBCHandling == 1)
    	{
    		u_n.addBc(bcInfo.get("Dirichlet"));
    		v_n.addBc(bcInfo.get("Dirichlet"));
    	}
    	ex.setIntegrationElements(basisH);
    	ex.initSystem();
    	ex.assemble(u_n * meas(G)* v_n.tr()); // * meas(G));
    	gsSparseMatrix<> P = ex.matrix();

    	// Determine matrix M (low_order * low_order)
    	gsExprAssembler<real_t> ex2(1,1);
    	geometryMap G2 = ex2.getMap(mp);
    	space w_n = ex2.getSpace(basisL ,1, 0);
    	if(typeBCHandling == 1)
    	{
    		w_n.addBc(bcInfo.get("Dirichlet"));
    	}
    	ex2.setIntegrationElements(basisL);
    	ex2.initSystem();
    	ex2.assemble(w_n * meas(G2) * w_n.tr()); // * meas(G2));
    	gsSparseMatrix<> M = ex2.matrix();


    	// Restrict Xfine to Xcoarse
    	gsMatrix<> temp = P*Xfine;

        if (lumping == 1)
    	{
        	// Lump the mass matrix (should be a gsSparseMatrix!)
        	gsVector<> one = gsVector<>::Ones(M.rows());
        	gsMatrix<> M_L = (M*one).asDiagonal();
    		//gsConjugateGradient<> CGSolver(M_L);
    		//CGSolver.solve(temp,Xcoarse);
		// Alternative approach (no solve necessary)
                gsMatrix<> M_L_inv_d = (M*one).array().inverse();
                gsMatrix<> M_L_inv = (M_L_inv_d).asDiagonal();
                Xcoarse = M_L_inv*temp;

    	}
    	else
    	{
    		gsConjugateGradient<> CGSolver(M);
        	CGSolver.setTolerance(1e-10);
    		CGSolver.solve(temp,Xcoarse);
    	}
  }
};

/** @brief The p-multigrid class implements a generic p-multigrid solver
 *  that can be customized by passing assembler and coarse
 *  solver as template arguments.
 *
 *  @note: This implementation assumes that all required prolongation/
 *  restriction operators are generated internally. Therefore, a
 *  problem-specific assembler has to be passed as template argument.
 */
template<class T, class CoarseSolver, class Assembler = void>
struct pMultigrid : public pMultigridBase<T>
{
private:

  /// Base class type
  typedef pMultigridBase<T> Base;

  /// Shared pointer to multi-patch geometry
  memory::shared_ptr<gsMultiPatch<T> > m_mp_ptr;

  /// Shared pointer to boundary conditions
  memory::shared_ptr<gsBoundaryConditions<T> > m_bcInfo_ptr;

  /// Vector of multi-basis objects
  vector<memory::shared_ptr<gsMultiBasis<T> > > m_basis;

  /// Vector of operator objects
  vector< gsSparseMatrix<T> > m_operator;

  /// Vector of assembler objects
  vector<Assembler> m_assembler;

public:

  // Constructor
  pMultigrid(const gsMultiPatch<T> & mp,
	     const gsMultiBasis<T> & basis,
	     const gsBoundaryConditions<T> & bcInfo)
  {
    m_mp_ptr = memory::make_shared_not_owned(&mp);
    m_bcInfo_ptr = memory::make_shared_not_owned(&bcInfo);
    m_basis.push_back(memory::make_shared_not_owned(&basis));
  }

public:

  ///  @brief Apply p-multigrid solver to given right-hand side on level l
  void solve(const gsFunctionExpr<T> & rhs,
  	         const gsFunctionExpr<T> & sol_exact,
  	               gsMatrix<T>       & x,
  	         const int               & numSmoothing,
  	               gsMatrix<T>         f,
  	         const int               & typeSolver,
  	               int               & iterTot,
  	         const int               & typeCycle,
  	         const int               & numLevels,
  	         //const int               & numRefine,
  	         //const int               & typeMultigrid,
  	         const int               & typeBCHandling,
  	               gsGeometry<>::Ptr   geo,
  	         const int               & lumping,
  	         const gsMatrix<>        & hp)
  {
      	// Generate sequence of bases on all levels
      	for (int i = 1; i < numLevels; i++)
      	{
          	m_basis.push_back(give(m_basis.back()->clone()));

		switch((int) hp(i-1,0) )
		{
			case 0 : m_basis.back()->degreeIncrease(); break;

			case 1 : m_basis.back()->uniformRefine();  break;

			case 2:  m_basis.back()->uniformRefine();
			m_basis.back()->degreeIncrease(); break;

		}
      	}

      	// Define coefficients CDR equation (= Poisson)
      	gsFunctionExpr<> coeff_diff, coeff_conv, coeff_reac;
      	(sol_exact.domainDim() == 1 ? coeff_diff = gsFunctionExpr<>("1",1) : coeff_diff = gsFunctionExpr<>("1","0","0","1",2));
      	(sol_exact.domainDim() == 1 ? coeff_conv = gsFunctionExpr<>("0",1) : coeff_conv = gsFunctionExpr<>("0","0",2));
      	(sol_exact.domainDim() == 1 ? coeff_reac = gsFunctionExpr<>("0",1) : coeff_reac = gsFunctionExpr<>("0",2));

      	// Generate sequence of assembler objects and assemble
      	for (typename vector<memory::shared_ptr<gsMultiBasis<T> > >::iterator it = m_basis.begin();
           it != m_basis.end(); ++it)
      	{
      		m_assembler.push_back(Assembler(*m_mp_ptr,*(*it).get(),*m_bcInfo_ptr,rhs, coeff_diff , coeff_conv , coeff_reac ,(typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche)));
	}

        // Obtain fine grid operator
        m_operator.resize(numLevels);

    	// Print hierarchy information and assemble
    	for (int i = 0; i < numLevels; i++)
    	{
        	m_assembler[i].assemble();

		// REDISCRETIZATION APPROACH
                //m_operator[i] = m_assembler[i].matrix();

		if(typeSolver == 1)
        	{
			gsInfo << "\nLevel: " << i+1 << endl;
			gsInfo << "Order of discretization: " << m_basis[i]->degree() << endl;
	    		gsInfo << "Degrees of freedom: " << m_basis[i]->totalSize() << endl;
		}
    	}

        m_operator[numLevels-1] = m_assembler[numLevels-1].matrix();

        // ALTERNATIVE: APPLY GALERKIN CONDITION!!
      	for (int i = numLevels-1; i > 0; i--)
      	{
	        galerkincondition(m_operator[i-1], m_operator[i], numLevels-(numLevels-1-i), m_basis, /*lumping,*/ typeBCHandling, *m_bcInfo_ptr, geo);


	}

    	gsMatrix<> b;
    	typeSolver == 1 ? b = m_assembler.back().rhs() : b = f;
    	real_t r0 = (m_operator[numLevels-1]*x - b).norm();
    	real_t r = r0;
        if(typeSolver == 1)
    	{
        	gsInfo << "\nInitial residual norm: " << r0 << endl;
    	}
        real_t tol = 1e-8;
    	int iter = 1;

    	while(typeSolver == 1 ? r/r0 > tol : iter < 2)
    	{
        	// Call solver from base class
        	Base::solve(b, m_basis,  x, numLevels, numSmoothing, typeCycle, typeSolver, typeBCHandling, *m_bcInfo_ptr, geo, lumping);

        	r = (m_operator[numLevels-1]*x - b).norm();
        	iter++;
        	iterTot++;
    	}
    	if( typeSolver == 1)
    	{
	    	// Determine residual and L2 errpr
	    	gsField<> solMG = m_assembler.back().constructSolution(x);
	    	gsNormL2<real_t> L2Norm(solMG,sol_exact);
	    	real_t errorL2 = L2Norm.compute();
	    	gsInfo << "Residual after cycle " << iter-1 << " : "  << r << endl;
	    	gsInfo << "L2 error: " << errorL2 << endl;
/*

                // Write V-cycles to file
		stringstream myString;
	        myString << "RD" << numLevels << numRefine << numSmoothing << hp << ".txt";
	        ofstream myfile;
	        myfile.open ("Results//ILUK//" + myString.str());
	        myfile << iter-1;
	        myfile.close();

*/
	    	// Plot solution in Paraview
	    	gsInfo << "Plotting in Paraview...\n";
	    	gsWriteParaview<>(solMG, "multigridMG", 100*x.rows());
	    	gsField<> Exact( *m_mp_ptr, sol_exact, false );
	    	gsWriteParaview<>( Exact, "multigrid_exact", 100*x.rows());
    	}
  }

private:

  /// @brief Apply coarse solver
  virtual void solvecoarse(const gsMatrix<T>& rhs, gsMatrix<T>& x, const int& numLevels)
  {
        gsInfo << " Coarse solver is applied at level " << numLevels << "\n" << endl;
    	CoarseSolver solver(m_operator[0]);
    	solver.setTolerance(1e-10);
    	solver.solve(rhs, x);
  }

  /// @brief Apply fixed number of presmoothing steps
  virtual void presmoothing(const gsMatrix<T>& rhs,  gsMatrix<T>& x, const int& numLevels, const int& numSmoothing, gsMatrix<T> & fineRes)
  {
    	gsInfo << "Presmoothing is applied at level " << numLevels  <<  "\n" << endl;
    	gsInfo << "Norm residual before smoothing: " << (rhs - m_operator[numLevels-1]*x).norm() << endl;

    	for(int i = 0 ; i < numSmoothing ; i++)
    	{
	/*
		// Implementation of a overlapping, multiplicative Schwarz smoother
		real_t Blocksize = 3;
		real_t NumBlocks = m_operator[numLevels-1].rows() - (Blocksize-1);

		for(int i = 0 ; i < NumBlocks ; i++)
		{
			// Obtain block
			gsMatrix<> B = m_operator[numLevels-1].block(i,i,Blocksize, Blocksize);

		        // Find the inverse of B (should be adjusted)
			gsMatrix<> B_inv = B.inverse();

		        // Define projection matrix
			gsMatrix<> V = gsMatrix<>::Zero(Blocksize, m_operator[numLevels-1].rows());
		        for(int j = 0; j < Blocksize ; j++)
			{
				V(j,i+j) = 1;
			}

			// Apply smoothing
			x = x + V.transpose()*B_inv*V*(rhs-m_operator[numLevels-1]*x);

		}

		// Apply ILU(k) as smoother
                Eigen::SuperILU<gsSparseMatrix<real_t>> super_ilu;
                super_ilu.analyzePattern(m_operator[numLevels-1]);
		super_ilu.factorize(m_operator[numLevels-1]);
		gsMatrix<> e;
                gsMatrix<> d = rhs-m_operator[numLevels-1]*x;
		e = super_ilu.solve(d);
                x = x + e;
*/
		// Try to apply ILUT as solver/preconditioner
		Eigen::IncompleteLUT<real_t, index_t> ilu;
		ilu.setFillfactor(1);
		ilu.compute(m_operator[numLevels-1]);
		//gsMatrix<> P = ilu.m_P;
                gsMatrix<> e;
                gsMatrix<> d = rhs-m_operator[numLevels-1]*x;
		ilu._solve_impl(d, e);
                x = x + e;


		//internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs);
    	}

    	gsInfo << "Norm residual after smoothing: " << (rhs - m_operator[numLevels-1]*x).norm() << endl;
    	fineRes = m_operator[numLevels-1]*x - rhs;
  }

  /// @brief Apply fixed number of postsmoothing steps
  virtual void postsmoothing(const gsMatrix<T>& rhs,  gsMatrix<T>& x, const int& numLevels, const int& numSmoothing, gsMatrix<T> & fineCorr, const int& typeSolver)
  {
  	    GISMO_UNUSED(typeSolver); // candidate to remove
      	gsInfo << "Postsmoothing is applied at level " << numLevels  <<  "\n" << endl;
    	real_t omega = 1.0;
    	x = x - omega*fineCorr;

    	gsInfo << "Norm residual before postsmoothing: " << (rhs - m_operator[numLevels-1]*x).norm() << endl;

    	for(int i = 0 ; i < numSmoothing ; i++)
    	{
/*
		// Implementation of a overlapping, multiplicative Schwarz smoother
		real_t Blocksize = 3;
		real_t NumBlocks = m_operator[numLevels-1].rows() - (Blocksize-1);

		for(int i = 0 ; i < NumBlocks ; i++)
		{
			// Obtain block
			gsMatrix<> B = m_operator[numLevels-1].block(i,i,Blocksize, Blocksize);

		        // Find the inverse of B (should be adjusted)
			gsMatrix<> B_inv = B.inverse();

		        // Define projection matrix
			gsMatrix<> V = gsMatrix<>::Zero(Blocksize, m_operator[numLevels-1].rows());
		        for(int j = 0; j < Blocksize ; j++)
			{
				V(j,i+j) = 1;
			}

			// Apply smoothing
			x = x + V.transpose()*B_inv*V*(rhs-m_operator[numLevels-1]*x);

		}

		// Apply ILU(k) as smoother
                Eigen::SuperILU<gsSparseMatrix<real_t>> super_ilu;
		super_ilu.analyzePattern(m_operator[numLevels-1]);
		super_ilu.factorize(m_operator[numLevels-1]);
		gsMatrix<> e;
                gsMatrix<> d = rhs-m_operator[numLevels-1]*x;
		e = super_ilu.solve(d);
                x = x + e;
*/
		// Try to apply ILUT as smoother
		Eigen::IncompleteLUT<real_t, index_t> ilu;
		ilu.setFillfactor(1);
		ilu.compute(m_operator[numLevels-1]);
		gsMatrix<> e;
                gsMatrix<> d = rhs-m_operator[numLevels-1]*x;
		ilu._solve_impl(d, e);
                x = x + e;

		//( typeSolver == 3 ? internal::reverseGaussSeidelSweep(m_operator[numLevels-1],x,rhs) : internal::gaussSeidelSweep(m_operator[numLevels-1],x,rhs));
    	}
    	gsInfo << "Norm residual after postsmoothing: " << (rhs - m_operator[numLevels-1]*x).norm() << endl;
  }


  /// @brief Apply Galerkin Condition
  virtual void galerkincondition(gsSparseMatrix<T> & Acoarse,
  	                       const gsSparseMatrix<T> & Afine,
  	                       const int               & Level,
  	 vector<memory::shared_ptr<gsMultiBasis<T> > >   m_basis,
  	                       //const int               & lumping,
  	                       const int               & typeBCHandling,
  	                       gsBoundaryConditions<T>   bcInfo,
  	                       gsGeometry<>::Ptr         geo)
  {
	// Define the low and high order basis
        gsInfo << "Level: " << Level << endl;
    	gsMultiBasis<> basisL = *m_basis[Level-2];
    	gsMultiBasis<> basisH = *m_basis[Level-1];

     	// Determine matrix P (high_order * low_order)
    	gsMultiPatch<> mp(*geo);
    	typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    	gsExprAssembler<real_t> ex(1,1);
    	geometryMap G = ex.getMap(mp);
    	typedef gsExprAssembler<real_t>::variable variable;
    	typedef gsExprAssembler<real_t>::space    space;
    	space v_n = ex.getSpace(basisH ,1, 0);
    	space u_n = ex.getTestSpace(v_n , basisL);
    	if(typeBCHandling == 1)
    	{
    		v_n.addBc(bcInfo.get("Dirichlet"));
    		u_n.addBc(bcInfo.get("Dirichlet"));
    	}
    	ex.setIntegrationElements(basisH);
    	ex.initSystem();
    	ex.assemble(u_n*meas(G) * v_n.tr()); //*meas(G));
    	gsSparseMatrix<> P = ex.matrix().transpose();

    	// Determine matrix M (high_order * high_order)
    	gsExprAssembler<real_t> ex2(1,1);
    	geometryMap G2 = ex2.getMap(mp);
    	space w_n = ex2.getSpace(basisH ,1, 0);
    	if(typeBCHandling == 1)
    	{
    		w_n.addBc(bcInfo.get("Dirichlet"));
    	}
    	ex2.setIntegrationElements(basisH);
    	ex2.initSystem();
    	ex2.assemble(w_n*meas(G2) * w_n.tr()); //*meas(G2));
    	gsSparseMatrix<> M_H = ex2.matrix();

	// Determine matrix M (low_order * low_order)
    	gsExprAssembler<real_t> ex3(1,1);
    	geometryMap G3 = ex3.getMap(mp);
    	space x_n = ex3.getSpace(basisL,1, 0);
    	if(typeBCHandling == 1)
    	{
    		x_n.addBc(bcInfo.get("Dirichlet"));
    	}
    	ex3.setIntegrationElements(basisL);
    	ex3.initSystem();
    	ex3.assemble(x_n * meas(G3) * x_n.tr()); // * meas(G2));
    	gsSparseMatrix<> M_L = ex3.matrix();

        // Obtain inverse lumped mass matrices
        gsVector<> one_L = gsVector<>::Ones(M_L.rows());
        gsMatrix<> M_L_inv_d = (M_L*one_L).array().inverse();
	gsMatrix<> M_L_inv = (M_L_inv_d).asDiagonal();

	gsVector<> one_H = gsVector<>::Ones(M_H.rows());
        gsMatrix<> M_H_inv_d = (M_H*one_H).array().inverse();
	gsMatrix<> M_H_inv = (M_H_inv_d).asDiagonal();

	// Determine coarse grid operator
	gsMatrix<> temp =  M_L_inv*P.transpose()*Afine*M_H_inv*P;
        Acoarse = temp.sparseView();


  }

};

/** @brief The p-multigrid class implements a generic p-multigrid solver
 *  that can be customized by passing assembler and coarse
 *  solver as template arguments.
 *
 *  @note: This implementation assumes that all required prolongation/
 *  restriction operators are generated externally and provided as
 *  constant references through the constructor. Therefore, no assembler
 *  is passed as template parameter.
 */
template<class T, class CoarseSolver>
struct pMultigrid<T, CoarseSolver, void> : public pMultigridBase<T>
{
  // Default constructor
  pMultigrid()
  {
    	gsInfo << "The specific case";
  }
};

int main(int argc, char* argv[])
{
	index_t numDegree = 2;
  	index_t numRefine = 4;
  	index_t numSmoothing = 1;
  	index_t numLevels = 2;
  	index_t numBenchmark = 3;
  	index_t typeSolver = 1;
  	index_t typeCycle = 1;
  	index_t typeMultigrid = 1;
  	index_t typeBCHandling = 2;
  	index_t lumping = 1;
        std::string coarsening = "h";

  	// Command line argument parser
  	gsCmdLine cmd("This file solves the Poisson equation with a p-multigrid / h-multigrid method");

  	// Add command line arguments
  	cmd.addInt("p", "Degree", "Number of order elevation steps", numDegree);
  	cmd.addInt("r", "Refinement", "Number of global refinements", numRefine);
  	cmd.addInt("v", "Smoothing", "Number of pre/post smoothing steps", numSmoothing);
  	cmd.addInt("l", "Levels", "Number of levels in multigrid method", numLevels);
  	cmd.addInt("b", "Benchmark", "Number of the benchmark",numBenchmark);
  	cmd.addInt("s", "Solver", "Type of solver: (1) p-mg as stand-alone solver (2) BiCGStab prec. with p-mg (3) CG prec. with p-mg", typeSolver);
  	cmd.addInt("t", "Multigridtype", "p-multigrid (1) or h-multigrid (2)", typeMultigrid);
  	cmd.addInt("m", "Multigridcycle", "Type of cycle, eather V-cycle (1) or W-cycle (2)", typeCycle);
  	cmd.addInt("d", "BoundaryConditionHandling", "Handles Dirichlet BC's by elimination (1) or Nitsche's method (2)", typeBCHandling);
  	cmd.addInt("L", "Lumpedprojection", "Restriction and Prolongation performed with the lumped (1) or consistent (2) mass matrix", lumping);
        cmd.addString("z", "Coarsening", "Expression that defines coarsening strategy", coarsening);

        // Read parameters from command line
  	try { cmd.getValues(argc,argv);  } catch (int rv) { return rv; }

  	// Initialize solution, rhs and geometry
  	std::string solution_exact,rhs_exact;
  	gsGeometry<>::Ptr geo;

  	switch(numBenchmark)
  	{
	  	case 1 : gsInfo << "Poisson equation the unit square (1)" << endl;
		         solution_exact = "sin(pi*x)*sin(pi*y)";
		         rhs_exact = "pi*pi*sin(pi*x)*sin(pi*y)+ pi*pi*sin(pi*x)*sin(pi*y)";
		         geo = gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0); break;

		case 2 : gsInfo << "Poisson equation on the unit square (2)" << endl;
		         solution_exact = "sin(pi*(x+1/2))*sin(pi*(y+1/2))";
		         rhs_exact = "pi*pi*sin(pi*(x+1/2))*sin(pi*(y+1/2))+ pi*pi*sin(pi*(x+1/2))*sin(pi*(y+1/2))";
			 geo = gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0); break;

		case 3:  gsInfo << "Poisson equation on the quarter annulus (1)" << endl;
		         solution_exact = "-(x*x+y*y-1)*(x*x+y*y-4)*x*y*y";
		         rhs_exact = "2*x*(22*x*x*y*y+21*y*y*y*y-45*y*y+x*x*x*x-5*x*x+4)";
		         geo = gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0); break;

		case 4:  gsInfo << "Poisson equation on the quarter annulus (2)" << endl;
		         solution_exact = "(x^2+y^2-3*sqrt(x^2+y^2)+2)*sin(2*atan(y/x))";
		         rhs_exact = "(8-9*sqrt(x^2 + y^2))*sin(2*atan(y/x))/(x^2+y^2)";
			 geo = gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0); break;
		case 5:  gsInfo << "Poisson equation on the line segment [0,1] " << endl;
		         solution_exact = "sin(pi*x)";
		         rhs_exact = "pi*pi*sin(pi*x)";
		         geo = gsNurbsCreator<>::BSplineUnitInterval(numDegree); break;
 	}

  	// Assign geometry, solution and rhs
  	gsMultiPatch<> mp(*geo);
  	gsFunctionExpr<> sol_exact(solution_exact,mp.geoDim());
  	gsFunctionExpr<> f(rhs_exact,mp.geoDim());

  	// Construct two bases (low and high order/refinement level)
  	gsMultiBasis<> basisL(mp);
  	gsMultiBasis<> basisH(mp);
        gsMatrix<> hp = gsMatrix<>::Zero(numLevels-1,1);

        // Read string from command line
	gsInfo << coarsening << endl;
	real_t numRefH = 0;
	real_t numRefP = 0;
	real_t numRefZ = 0;

	for( int i = 0; i < numLevels-1 ; ++i)
        {
		if( coarsening[i] == 'h')
                {
			gsInfo << "There is a h on position " << i << endl;
			hp(i,0) = 1;
                        numRefH = numRefH + 1;
		}
		else if( coarsening[i] == 'p')
		{
			gsInfo << "There is a p on position " << i << endl;
			hp(i,0) = 0;
                        numRefP = numRefP + 1;
		}
		else
		{
			gsInfo << "There is a z on position " << i << endl;
			hp(i,0) = 2;
                        numRefZ = numRefZ + 1;
		}
	}

        for (int i = 0; i < numRefine - numRefH - numRefZ; ++i)
  	{
        	basisL.uniformRefine();
        }
        for (int i = 0; i < numRefine ; ++i)
  	{
        	basisH.uniformRefine();
        }
        basisL.degreeIncrease(numDegree-numRefP-numRefZ-1);
        basisH.degreeIncrease(numDegree-1);

  	// Define boundary conditions
  	gsBoundaryConditions<> bcInfo;

  	if(numBenchmark != 5)
  	{
  	  	bcInfo.addCondition(boundary::west,  condition_type::dirichlet, &sol_exact, 0,false);
	  	bcInfo.addCondition(boundary::east,  condition_type::dirichlet, &sol_exact, 0,false);
	  	bcInfo.addCondition(boundary::south, condition_type::dirichlet, &sol_exact, 0,false);
	  	bcInfo.addCondition(boundary::north, condition_type::dirichlet, &sol_exact, 0,false);
  	}
  	else
  	{
  	  	bcInfo.addCondition(boundary::west,  condition_type::dirichlet, &sol_exact, 0,false);
	  	bcInfo.addCondition(boundary::east,  condition_type::dirichlet, &sol_exact, 0,false);
  	}

  	int iterTot = 1;
  	pMultigrid<real_t, gsConjugateGradient<real_t>,gsCDRAssembler<real_t> > My_MG(mp, basisL, bcInfo);

  	// CDR assembler (needed for Krylov solver)
  	gsFunctionExpr<> coeff_diff, coeff_conv, coeff_reac;
  	(sol_exact.domainDim() == 1 ? coeff_diff = gsFunctionExpr<>("1",1) : coeff_diff = gsFunctionExpr<>("1","0","0","1",2));
  	(sol_exact.domainDim() == 1 ? coeff_conv = gsFunctionExpr<>("0",1) : coeff_conv = gsFunctionExpr<>("0","0",2));
  	(sol_exact.domainDim() == 1 ? coeff_reac = gsFunctionExpr<>("0",1) : coeff_reac = gsFunctionExpr<>("0",2));

  	// Construct assembler and assemble
  	gsCDRAssembler<real_t> pa(mp, basisH, bcInfo, f, coeff_diff, coeff_conv, coeff_reac, (typeBCHandling == 1 ? dirichlet::elimination : dirichlet::nitsche));
  	pa.assemble();

  	// Apply p-Multigrid as stand-alone solver
  	gsMatrix<real_t> x = gsMatrix<>::Random(pa.matrix().rows(),1);

/*
	  // Mass matrix needed
	  gsExprAssembler<real_t> ex2(1,1);
	  typedef gsExprAssembler<real_t>::variable variable;
	  typedef gsExprAssembler<real_t>::space    space;
	  space w_n = ex2.getSpace(basisH ,1, 0);    // trial function
	  //w_n.addBc(bcInfo.get("Dirichlet"));
	  ex2.setIntegrationElements(basisH);
	  ex2.initSystem();
	  ex2.assemble(w_n * w_n.tr());
	  gsSparseMatrix<> M = ex2.matrix();


	  // Define matrix with eigenvectors
	  gsMatrix<> A = pa.matrix();
	  Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<>::Base > eig(A, M);
	  gsMatrix<> ev = eig.eigenvectors();
	  gsMatrix<> eval = eig.eigenvalues();
	  ev.colwise().normalize();

	  // Define matrix with eigenvectors
	  ofstream myfile2;
	  myfile2.open ("SSOR45Initial.txt");
	  myfile2 << ev;
	  myfile2.close();

	  gsMatrix <real_t> ResultMatrix = gsMatrix<>::Zero(pa.matrix().rows(),pa.matrix().rows());

	  for(int i = 0 ; i < pa.matrix().rows(); i++)
	  {
		// Make p-Multigrid structure
		pMultigrid<real_t, gsConjugateGradient<real_t>, gsCDRAssembler<real_t> > My_MG1(mp, basisL, bcInfo);
		gsMatrix<> Test = ev.col(i);
	 	// Apply p-Multigrid (with not needed last parameter)
	  	My_MG1.solve(f, sol_exact, Test, numSmoothing, x, typeSolver, iterTot, typeCycle, numLevels, typeBCHandling, geo, lumping, hp);
		ResultMatrix.col(i) = Test;
	  }

	  ofstream myfile1;
	  myfile1.open ("SSOR45Result.txt");
	  myfile1 << ResultMatrix;
	  myfile1.close();
	  exit(1);
*/
/*
	gsMatrix<> IterMatrix;
        for( int i = 0; i < pa.matrix().rows() ; i++)
	{
		gsMatrix<> e = gsMatrix<>::Zero(pa.matrix().rows(),1);
                e(i,0) = 1;
		My_MG.solve(f, sol_exact, e, numSmoothing, x, typeSolver, iterTot, typeCycle, numLevels, typeBCHandling, geo, lumping, hp);
		IterMatrix.col(i) = e;
	}

	ofstream myfile1;
	myfile1.open ("IterationMatrix25.txt");
	myfile1 << IterMatrix;
	myfile1.close();
	exit(1);
*/

 	if(typeSolver == 1)
  	{
        	My_MG.solve(f, sol_exact, x, numSmoothing, x, typeSolver, iterTot, typeCycle, numLevels, typeBCHandling, geo, lumping, hp);
        	return 0;
  	}

  	// Apply BiCGStab or CG
  	gsVector <> r = pa.rhs() - pa.matrix() * x;
  	gsVector<> r0 = r;
  	real_t maxIter = pa.matrix().rows();
  	real_t tol = 1e-8;
  	int i = 1;

  	// Set up, determine initial L2 error
  	gsField<> solKrylov;
  	real_t oldL2Err;
  	gsField<> sol = pa.constructSolution(x);
  	oldL2Err = computeL2Distance(sol, sol_exact, false, 3*pa.matrix().rows());
  	real_t oldResNorm = r0.norm();

  	// Perform BiCGStab
  	if(typeSolver == 2)
  	{
                gsInfo << "BiCGStab is applied as solver\n" << endl;
	  	// Define vectors needed in BiCGStab
	  	gsVector<> t = gsVector<>::Zero(pa.matrix().rows());
	  	gsVector<> s = gsVector<>::Zero(pa.matrix().rows());
	  	gsVector<> p = gsVector<>::Zero(pa.matrix().rows());
	  	gsVector<> v = gsVector<>::Zero(pa.matrix().rows());
	  	gsMatrix<> y = gsMatrix<>::Zero(pa.matrix().rows(),1);
	  	gsMatrix<> z = gsMatrix<>::Zero(pa.matrix().rows(),1);

	  	real_t alp = 1;
	  	real_t rho = 1;
	  	real_t w = 1;

	  	// Set residual norm and #restarts
	  	real_t r0_sqnorm = r0.dot(r0);
	  	int restarts = 0;

	  	// Perform BiCGStab
	  	while(r.norm()/r0.norm() > tol && i < maxIter)
	  	{
			// Construct P-Multigrid objects
			pMultigrid<real_t, gsConjugateGradient<real_t>,gsCDRAssembler<real_t> > My_MG1(mp, basisL, bcInfo);
			pMultigrid<real_t, gsConjugateGradient<real_t>,gsCDRAssembler<real_t> > My_MG2(mp, basisL, bcInfo);

		  	real_t rho_old = rho;
		  	rho = r0.dot(r);
		  	if (abs(rho) < 1e-32*r0_sqnorm)
		  	{
		      		gsInfo << "Residual too orthogonal, restart with new r0 \n";
		      		r = pa.rhs() - pa.matrix()*x;
		      		r0 = r;
		      		rho = r0_sqnorm = r.dot(r);
		      		if (restarts++ == 0)
		      		{
		    			i=0;
		      		}
		  	}
		  	real_t beta = (rho/rho_old)*(alp/w);
		  	p = r + beta*(p - w*v);

			// Apply preconditioning by solving Ay = p
			y.setZero();
			My_MG1.solve(f, sol_exact, y, numSmoothing, p, typeSolver, iterTot, typeCycle, numLevels, typeBCHandling, geo, lumping, hp);
			v = pa.matrix()*y;
			alp = rho/(r0.dot(v));
			s = r - alp*v;

			// Apply preconditioning by solving Az = s
			z.setZero();
			My_MG2.solve(f, sol_exact, z, numSmoothing, s, typeSolver, iterTot, typeCycle, numLevels, typeBCHandling, geo, lumping, hp);
			t = pa.matrix()*z;
			if (t.dot(t) > 0)
				w = t.dot(s)/t.dot(t);
			else
				w = 0;

			x= x + alp*y + w*z;
			r = s - w*t;

		        // Print information about residual and L2 error
		  	gsInfo << "BiCGStab iteration: " << i << "      |   Residual norm: "   << left << setw(15) << r.norm() << "           reduction:  1 / " << setprecision(3) << (oldResNorm/r.norm()) <<   			setprecision	(6) << "\n";
		  	oldResNorm = r.norm();
		  	solKrylov = pa.constructSolution(x);
		  	gsNormL2<real_t> L2Norm(solKrylov,sol_exact);
			real_t l2Err = L2Norm.compute();
		        gsInfo << "                           |   L2 error:      "
			  << left << setw(15) << l2Err << "           reduction:  1 / "  << setprecision(3) << (oldL2Err/l2Err) << setprecision(6) << "\n";
		  	oldL2Err = l2Err;

			++i;
	  }
  }
  else if(typeSolver == 3)
  {
                  gsInfo << "CG is applied as solver\n" << endl;
		  // Apply preconditioner
		  pMultigrid<real_t, gsConjugateGradient<real_t>,gsCDRAssembler<real_t> > My_MG2(mp, basisL, bcInfo);
		  gsMatrix<> z1 = gsMatrix<>::Zero(pa.matrix().rows(),1);
		  My_MG2.solve(f, sol_exact, z1, numSmoothing, r0, typeSolver, iterTot, typeCycle, numLevels, typeBCHandling, geo, lumping, hp);
		  gsVector<> z = z1;
		  gsVector<> p = z;
		  real_t alpha, beta;

		  while(r.norm()/r0.norm() >  tol && i < maxIter)
		  {
		    	// Determine alpha
		    	alpha = r.transpose()*z;
		    	alpha = alpha/(p.transpose()*pa.matrix()*p);

		    	// Update solution and residual
		    	x = x + alpha*p;
		    	gsVector<> r_new = r - alpha*pa.matrix()*p;

		    	// Obtain new values
		    	pMultigrid<real_t, gsConjugateGradient<real_t>,gsCDRAssembler<real_t> > My_MG3(mp, basisL, bcInfo);
		    	gsMatrix<> z2 = gsMatrix<>::Zero(pa.matrix().rows(),1);
		    	My_MG3.solve(f, sol_exact, z2, numSmoothing,r_new, typeSolver, iterTot, typeCycle, numLevels, typeBCHandling, geo,lumping, hp);
		    	gsVector<> z3 = z2;

		    	// Determine beta
		    	beta = z3.transpose()*r_new;
		    	beta = beta/(z.transpose()*r) ;
		    	p = z3 + beta*p;
		    	z = z3;
		    	r = r_new;

		    	// Print information about residual and L2 error
		    	gsInfo << "CG iteration: " << i << "       |  Residual norm: "   << left << setw(15) << r.norm() << "            reduction:  1 / " << setprecision(3) << (oldResNorm/r.norm()) <<   		    			    setprecision	(6) << "\n";
		    	oldResNorm = r.norm();
		    	solKrylov = pa.constructSolution(x);
		    	gsNormL2<real_t> L2Norm(solKrylov,sol_exact);
		    	real_t l2Err = L2Norm.compute();
		    	gsInfo << "                      |  L2 error:      "
		    	<< left << setw(15) << l2Err << "            reduction:  1 / "  << setprecision(3) << (oldL2Err/l2Err) << setprecision(6) << "\n";
		    	oldL2Err = l2Err;

		    	++i;
		   }
   	}

   	// Plot solution in Paraview
   	gsInfo << "Plotting in Paraview...\n";
   	gsWriteParaview<>(solKrylov,"Krylov", 3*x.rows());
   	gsField<> Exact( mp, sol_exact, false );
   	gsWriteParaview<>(Exact, "Krylov_exact", 3*x.rows());

   	return 0;

}
