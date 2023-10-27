#include <gismo.h>
#include <gsTrilinos/gsTrilinos.h>
#include <gsSolver/gsSolverUtils.h>
#include <fstream>
#include <algorithm>  
#include <iostream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <cmath>
#include <gsSolver/gsConjugateGradient.h>
#include <gsCore/gsLinearCombinationOfFunctionsFunction.h>
#include <gsTensor/gsTensorTools.h>
#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsGeneralizedPoissonAssembler.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsGaussSurfaceAssembler.h>

using namespace std;
using namespace gismo;

/** @brief The MPM base class provides the basic methods 
 *  (timestep, assembleMatrix, assembleVector) for implementing 
 *   different versions of the MPM. 
 */

template <class T>
struct MPMBase
{

public:

	// @brief Apply all steps in a single timestep (Assuming m_p = 1...)
	virtual void timestep(gsMatrix<T> & u_i, gsMatrix<T> & v_i, T dt, memory::shared_ptr<gsMultiBasis<T> > & m_basis, memory::shared_ptr<gsGeometry<T> > & m_geo)   
	{
		// Declare stress vector
		//gsMatrix<T> stress;

		// Determine Mass matrix		
		gsSparseMatrix<T> M;
		assembleMatrix(M, m_basis, m_geo);
                
		// Determine Force vectors
		gsSparseMatrix<T> F_grav, F_int;
		assembleVector(F_grav, F_int, /*stress,*/ m_basis, m_geo);
		gsVector<T> Id = gsVector<T>::Ones(m_basis->totalSize(),1);				
		gsVector<T> F_total = F_grav*Id;
				
		// Solve the equation of motion
		gsInfo << "Solving the equation of motion" << endl;
		gsSparseSolver<>::CGIdentity solverCGI;
    		solverCGI.compute(M);
    		gsMatrix<T> a_i = solverCGI.solve(F_total);
    		gsInfo << "Equation of motion solved " << endl;
	
		// Determine velocity via L2 projection (not used yet)
		gsInfo << "Update velocity via momentum (not being used yet!)" << endl;
		
		// Update properties based on acceleration
		gsInfo << "Update particle properties / time integration" << endl;

		// Obtain DOF which are nonzero on boundary
		gsMatrix<index_t> allBound = m_basis->basis(0).allBoundary();

		// Apply homogeneous Dirichlet boundary condition (with loop... :( ) 
		for (int i = 0 ; i < allBound.rows() ; i++)
		{
			a_i(allBound(i,0),0) = 0;
		}

		// Update properties at DOF (Modified Euler)		
		v_i = v_i + dt*a_i;
		u_i = u_i + dt*v_i;
		
		// Update stresses 
	
		//dstrain = v_i*dt;
		//stress = stress + E*dstrain;
	}		
	

	// @brief Assembly procedure for Mass matrix
	virtual void assembleMatrix(gsSparseMatrix<T> & Mass, memory::shared_ptr<gsMultiBasis<T> > & m_basis, memory::shared_ptr<gsGeometry<T> > & m_geo)
	{
		gsInfo << "Matrix is assembled..." << endl;

		// Define geometry evaluator
		typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_VALUE|NEED_MEASURE, *m_geo));
		
		// Define global system
		gsDofMapper dof(m_basis->basis(0));
		dof.finalize();
		gsSparseSystem<T> system (dof);

		gsGaussRule<T> quRule(gsGaussSurfaceAssembler<T>::getNumIntNodesFor(m_basis->basis(0)));
		gsMatrix<T> quNodes;
		gsVector<T> quWeights;
		gsFuncData<T> bdata(NEED_VALUE|NEED_ACTIVE|SAME_ELEMENT);
		
		typename gsBasis<T>::domainIter domIter = m_basis->basis(0).makeDomainIterator();

		for (; domIter->good(); domIter->next() )
		{
			// Compute quadrature rule 
			quRule.mapTo(domIter->lowerCorner(), domIter->upperCorner(), quNodes, quWeights);
			
			// Determine active functions
			m_basis->basis(0).compute(quNodes, bdata);

			// Determine local matrices
			gsMatrix<T> localMat;
			localMat.setZero(bdata.actives.rows(), bdata.actives.rows());

			// Values at quadrature points
			geoEval->evaluateAt(quNodes);
			gsMatrix<T> values = bdata.values[0];

			// Loop over quadrature points for contributions
			for (index_t k = 0 ; k < values.cols(); k++)
			{
				// Obtain integration weight with geometry measure
				const T weight = quWeights[k]*geoEval->measure(k);
			
				// Add contribution to local Matrix
				localMat.noalias() += weight*(values.col(k)*values.col(k).transpose());
			} 
			
			// Push contribution to system
			gsMatrix<T> eliminatedDofs;
			system.pushToMatrix(localMat, bdata.actives, eliminatedDofs);
		}
		gsInfo << "Assembly done! " << endl;
		Mass = system.matrix();
	}

	// @brief Assembly procedure Force vectors
	virtual void assembleVector(gsSparseMatrix<T> & F_grav,
		                        gsSparseMatrix<T> & F_int,
		                        // const gsMatrix<T> & stress,
		     memory::shared_ptr<gsMultiBasis<T> > & m_basis,
		     memory::shared_ptr<gsGeometry<T> >   & m_geo)
	{
		gsInfo << "Vector is assembled..." << endl;
		
		// Define constants
		const T g = 9.81;
			
		// Define geometry evaluator
		typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_VALUE|NEED_MEASURE, *m_geo));
		
		// Define global system
		gsDofMapper dof(m_basis->basis(0));
		dof.finalize();
		gsSparseSystem<T> system_grav (dof);
		dof.finalize(); 
		gsSparseSystem<T> system_int (dof);

		gsGaussRule<T> quRule(gsGaussSurfaceAssembler<T>::getNumIntNodesFor(m_basis->basis(0)));
		gsMatrix<T> quNodes;
		gsVector<T> quWeights;
		gsFuncData<T> bdata(NEED_VALUE|NEED_ACTIVE|SAME_ELEMENT);

		typename gsBasis<T>::domainIter domIter = m_basis->basis(0).makeDomainIterator();
		for (; domIter->good(); domIter->next() )
		{			
			// Apply general quadrature rule
			quRule.mapTo(domIter->lowerCorner(), domIter->upperCorner(), quNodes, quWeights);
	
			// Determine active functions
			m_basis->basis(0).compute(quNodes, bdata);

			// Determine local matrices
			gsMatrix<T> localF_grav;
			gsMatrix<T> localF_int; 
			localF_grav.setZero(bdata.actives.rows(), bdata.actives.rows());
			localF_int.setZero(bdata.actives.rows(), bdata.actives.rows());
			
			// Values at quadrature points.
			geoEval->evaluateAt(quNodes);
			gsMatrix<T> & values = bdata.values[0];
			//gsMatrix<T> & derivs = bdata.values[1];
            
			// Loop over quadrature points for contributions
			for (index_t k = 0 ; k < values.cols(); k++)
			{
				const T weight = quWeights[k]*geoEval->measure(k);
				localF_grav.noalias() += weight*values.col(k)*values.col(k).transpose()*g;
				localF_int.noalias()  += weight*values.col(k)*values.col(k).transpose()*g; // How to do this in 2D?				
			} 
			
			// Push contribution rhs to system
			gsMatrix<T> eliminatedDofs;
			system_grav.pushToMatrix(localF_grav, bdata.actives);
			system_int.pushToMatrix(localF_int, bdata.actives);
		}
		gsInfo << "Assembly done! " << endl;
		F_grav = system_grav.matrix();
		F_int = system_int.matrix();
	}
};

/** @brief The MPM class provides the basic methods 
 *  solve for solving with the MPM 
 */


template <typename T>
struct MPM : public MPMBase<T>
{

private:  
	// Base class type
	typedef MPMBase<T> Base;

	// Shared pointer to multi-patch geometry
	memory::shared_ptr<gsMultiPatch<T> > m_mp_ptr;

	// Vector of multi-basis objects
	memory::shared_ptr<gsMultiBasis<T> > m_basis;

	// Shared pointer to geometry
	memory::shared_ptr<gsGeometry<T> > m_geo;	
	
public:
	// Constructor
  	MPM(const gsMultiPatch<T> & mp, const gsMultiBasis<T> & basis, const gsGeometry<T>& geo)
	{ 	
        	m_mp_ptr = memory::make_shared_not_owned(&mp);
		m_basis = memory::make_shared_not_owned(&basis);
		m_geo = memory::make_shared_not_owned(&geo);
		
	}

public:

	//  @brief Apply p-multigrid solver to given right-hand side on level l
        void solve() 
        {
		gsMatrix<real_t> u_i = gsMatrix<real_t>::Zero(m_basis->totalSize(),1);
		gsMatrix<real_t> v_i = gsMatrix<real_t>::Zero(m_basis->totalSize(),1);

		// Determine timestep
		real_t dt = 0.01;

		for(int i = 1 ; i < 100 ; i++)
		{
			Base::timestep(u_i, v_i, dt, m_basis, m_geo);
			gsInfo << "Velocity norm in solve: " << v_i.norm() << endl;
		}

		// Plot solution in Paraview 
		gsGenericAssembler<T> assm(*m_mp_ptr,*m_basis);
		gsField<T> solMPM = assm.constructSolution(u_i);
  		gsInfo << "Plotting in Paraview...\n";
  		gsWriteParaview<>(solMPM,"MPM", 3*u_i.rows());
	}
};
 

int main()
{
	int numDegree = 2;	
	int numRefine = 2;
	
	gsInfo << "This file solves a PDE with MPM \n" << endl;
	
	// Exact solution and rhs 
	gsFunctionExpr<> u("sin(pi*(x+1/2))*sin(pi*(y+1/2))",2);
 	gsFunctionExpr<> f("pi*pi*sin(pi*(x+1/2))*sin(pi*(y+1/2))+ pi*pi*sin(pi*(x+1/2))*sin(pi*(y+1/2))",2);
  
	//gsGeometry<>::Ptr geo = gsNurbsCreator<>::BSplineSquare(2.0, 0.0, 0.0);
	gsGeometry<>::Ptr geo = gsNurbsCreator<>::BSplineFatQuarterAnnulus(1.0, 2.0); 	
	gsMultiPatch<> mp(*geo);
	gsMultiBasis<> basis(mp);

  	for (int i = 0; i < numRefine; ++i)
  	{
      		basis.uniformRefine();
  	}
	basis.degreeIncrease(numDegree-1);

	// Print some general info
	gsInfo << "Degree of the basis functions: " << basis.degree() << endl;
	gsInfo << "Number of basis functions: " << basis.totalSize() << endl;

  	// Make MPM structure
  	MPM<real_t> My_MPM(mp, basis, *geo);
	My_MPM.solve();
	
	return 0;
}







