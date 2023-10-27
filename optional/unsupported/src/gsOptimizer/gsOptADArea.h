



#pragma once

#include <gsIpopt/gsOptProblem.h>
#include <gsAssembler/gsCDRAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <gsModeling/gsSpringPatch.h>
#include <gsModeling/gsCoonsPatch.h>
#include <gsModeling/gsCrossApPatch.h>

namespace gismo
{


/**
 * @brief
 * Simple shape-optimization example, to demonstrate the AD in reverse mode with use of CoDiPack,
 * Optimization done with IPOpt.
 *
 *  Goal: maximize area preserving the original circumference.
 *
 * 	!Works for single patch geometries!
 * 	!Constraint implemented as penalty function, penalty parameter to be determined by user!
 *
 *	Author: A. Jaeschke
 */


class gsOptADArea : public gsOptProblem<real_t>
{
private:
//Efficient system solving
static void solveSystem_b(codi::RealReverse::TapeType* tape, codi::DataStore* data, codi::AdjointInterface<double>* interface) {

  gsSparseMatrix<codi::RealReverse>* matrix;
  gsMatrix<codi::RealReverse>* rhs;
  gsMatrix<codi::RealReverse>* sol;

  data->getData(matrix);
  data->getData(rhs);
  data->getData(sol);

  gsSparseSolver<codi::RealReverse>::LU solver;


  gsMatrix<codi::RealReverse> solAdj(*sol);
  gsMatrix<codi::RealReverse> rhsAdj(*sol);
  for(index_t i = 0; i < sol->size(); ++i) {
      solAdj[i] = (*sol)[i].getGradient();
  }

  solver.compute(matrix->transpose());
  rhsAdj = solver.solve(solAdj);

  for(index_t i = 0; i < sol->size(); ++i) {
      auto index = (*rhs)[i].getGradientData();
      tape->gradient(index) += rhsAdj[i].getValue();
  }
  for (int e=0; e<matrix->outerSize(); ++e) {
      for (gsSparseMatrix<codi::RealReverse>::InnerIterator it(*matrix,e); it; ++it) {
          int k = it.row();
          int l = it.col();
          codi::RealReverse& temp1 = matrix->at(k,l);
          tape->gradient(temp1.getGradientData()) += -rhsAdj[l].getValue() * (*sol)[k].getValue();
      }
  }
}


static void solveSystem_delete(codi::RealReverse::TapeType* tape, codi::DataStore* data) {}

static void solveSystem(gsSparseSolver<codi::RealReverse>::LU&  solver, const gsSparseMatrix<codi::RealReverse>& matrix,
	const gsMatrix<codi::RealReverse>& rhs, gsMatrix<codi::RealReverse>& sol) {

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setPassive();
  codi::DataStore* dataHandler = new codi::DataStore();
  dataHandler->addData(&matrix);
  dataHandler->addData(&rhs);
  solver.compute( matrix );
  sol = solver.solve( rhs );

  dataHandler->addData(&sol);

  tape.pushExternalFunction(&solveSystem_b, dataHandler, &solveSystem_delete);
  tape.setActive();
  for(index_t i = 0; i < sol.size(); ++i) {
      tape.registerInput(sol.at(i));
  }
}

public:

    gsOptADArea(bool flag_plot, std::string path, real_t penalty_parameter,
		int numElevate_,int numHref_,int basisDegree_, int info_dom_, int info_diff_)
    {
		//Defining the size of the problem and the bounds for IpOpt
        m_numDesignVars  = 0;
        gsMultiPatch<real_t> mp_paa;
		gsReadFile<real_t>(path, mp_paa);
		int size_temp = 0;
		for(index_t ip=0; ip < mp_paa.nPatches();ip++)
		{
			gsBSpline<real_t> * g = dynamic_cast<gsBSpline<real_t> *>(&mp_paa.patch(ip));
			m_numDesignVars  = m_numDesignVars + g->coefs().size();
			m_numConstraints = 0;
			m_numConJacNonZero = 0;
		}
		m_desLowerBounds.resize(m_numDesignVars);
		m_desUpperBounds.resize(m_numDesignVars);
		m_curDesign.resize(m_numDesignVars,1);
		//Restricting domain to boundary only
		for(index_t ip=0; ip < mp_paa.nPatches();ip++)
		{
			gsBSpline<real_t> * g = dynamic_cast<gsBSpline<real_t> *>(&mp_paa.patch(ip));
			gsMatrix<codi::RealReverse> coef(g->coefs().rows(), g->coefs().cols());
			for (int i=0; i<g->coefs().rows(); i++)
			{
				for (int j=0; j<g->coefs().cols(); j++)
					coef(i,j).setValue(g->coefs()(i,j));
			}
			const gsKnotVector<real_t> & k0 = g->basis().knots(0);
			std::vector<codi::RealReverse> v0 (k0.size());

			for (unsigned int i=0; i<k0.size(); i++)
			{
				v0[i].setValue(k0[i]);
			}
			gsKnotVector<codi::RealReverse> nk0(give(v0));
			//std::cout << nk0 << std::endl;
			//std::cout << coef << std::endl;
			gsBSpline<codi::RealReverse> vv(nk0,coef);
			//std::cout << vv << std::endl;

			mp.addPatch(vv);
		//Defining the bounds for IpOpt
        for (int i=0; i<vv.coefs().rows();i++)
        {
			for (int j=0; j<vv.coefs().cols();j++)
			{
				m_desLowerBounds[size_temp + j+i*vv.coefs().cols()] = -100.0;
				m_desUpperBounds[size_temp + j+i*vv.coefs().cols()] =  100.0;
				m_curDesign(size_temp + j+i*vv.coefs().cols(),0) = vv.coefs()(i,j).getValue();
			}
		}
		size_temp  = size_temp + vv.coefs().size();
	}

        // ----------------------- Essentials -----------
        p_resultGrad = gsVector<real_t> (m_numDesignVars);
        info_diff = info_diff_;
        info_dom = info_dom_;
        plot = flag_plot;
        if (plot)
			collection = new gsParaviewCollection ("OptSteps");

		iter = 0;
		numElevate  = numElevate_;
		numHref     = numHref_;
		basisDegree = basisDegree_;
		f = gsFunctionExpr<codi::RealReverse>("0.0", 2);
		gf = gsFunctionExpr<codi::RealReverse>("1.0", 2);
		coeff_A = gsFunctionExpr<codi::RealReverse>("1.0","0","0","1.0", 2);
		coeff_b = gsFunctionExpr<codi::RealReverse>("0.2","0.4", 2);
		coeff_c = gsFunctionExpr<codi::RealReverse>("0", 2);
		//Determination of original circumference (+plot)
		gsSpringPatch<codi::RealReverse> spring(mp);
		if (info_dom)
			gsInfo<<"Created a " << spring.compute() <<"\n";
		else
			spring.compute();
        gsMultiPatch<codi::RealReverse> mp_full = spring.result();
			for (gsMultiPatch<codi::RealReverse>::const_biterator bit = mp_full.bBegin(); bit != mp_full.bEnd(); ++bit)
			{
				BCs.addCondition( *bit, condition_type::dirichlet, &gf);
			}

			gsMultiBasis<codi::RealReverse> bases(mp_full);
			if (basisDegree)
				bases.setDegree(basisDegree);
			else if (numElevate)
				bases.degreeElevate(numElevate);
			for (int i=0; i<numHref ; i++)
				bases.uniformRefine();
			gsSparseSolver<codi::RealReverse>::LU solver;
			gsCDRAssembler<codi::RealReverse> galerkin(mp_full,bases,BCs,f,coeff_A,coeff_b,coeff_c,
				dirichlet::elimination, iFace::glue, stabilizerCDR::none);
			galerkin.assemble();
			gsMatrix<codi::RealReverse> solVector;
			//Efficient solving of the system
			solveSystem(solver, galerkin.matrix(), galerkin.rhs(), solVector);
			//Ineficient way
			//solver.compute( galerkin.matrix() );
			//solVector = solver.solve( galerkin.rhs() );
			gsField<codi::RealReverse> sol = galerkin.constructSolution(solVector);
			gsConstantFunction<codi::RealReverse> zero (0.0,2);
			// old version not supported anymore
			//gsNormL2Boundary<codi::RealReverse> nor = gsNormL2Boundary<codi::RealReverse>(sol,zero, false);
			//circ0 = nor.compute()*nor.compute();
			// new version
			 // Initiate the expression evaluator
			gsExprEvaluator<codi::RealReverse> ev;
			// Set the parameter mesh as the integration mesh
            ev.setIntegrationElements(bases);
			// Define integrant variables
			typedef gsExprEvaluator<codi::RealReverse>::geometryMap geometryMap;
			typedef gsExprEvaluator<codi::RealReverse>::variable variable;
			geometryMap G = ev.getMap(mp_full);
			variable u = ev.getVariable(zero, G);
			variable s = ev.getVariable(sol.function());
			circ0 = ev.integralBdr( (u - s).sqNorm() * nv(G).norm() ) ;

			
			std::cout << "Original Circumference " << circ0 << std::endl;

			penalty.setValue(penalty_parameter);
			if (plot)
			{
				gsInfo<<"Plotting in Paraview...Initial \n";
				gsWriteParaview<>( mp_full, "Initial", 1000);
				
			}
    }

    ~gsOptADArea()
    {
		if(plot)
		{
			collection->save();
			delete collection;
		}
	}


public:

    real_t evalObj( const gsAsConstVector<real_t> & u ) const
    {
        // return the value of the objective function
        compute(u, false);
        return p_result;
    }

    void gradObj_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
    {
		compute(u, true);
		result = p_resultGrad;
    }
    void evalCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
    {
		GISMO_NO_IMPLEMENTATION
    }

    void jacobCon_into( const gsAsConstVector<real_t> & u, gsAsVector<real_t> & result) const
    {
		GISMO_NO_IMPLEMENTATION
    }

private:

    using gsOptProblem<real_t>::m_numDesignVars;
    using gsOptProblem<real_t>::m_numConstraints;
    using gsOptProblem<real_t>::m_numConJacNonZero;

    using gsOptProblem<real_t>::m_desLowerBounds;
    using gsOptProblem<real_t>::m_desUpperBounds;

    using gsOptProblem<real_t>::m_conLowerBounds;
    using gsOptProblem<real_t>::m_conUpperBounds;

    using gsOptProblem<real_t>::m_conJacRows;
    using gsOptProblem<real_t>::m_conJacCols;

    using gsOptProblem<real_t>::m_curDesign;

    mutable gsVector<real_t> p_resultGrad;
	gsParaviewCollection * collection;

	codi::RealReverse penalty;

	codi::RealReverse circ0;

	bool plot;
	bool info_diff;
	bool info_dom;

    mutable real_t p_result;

    mutable int iter;

    mutable gsMultiPatch<codi::RealReverse> mp;

    int numElevate  = 0;
	int numHref     = 1;
	int basisDegree = 0;

	gsFunctionExpr<codi::RealReverse>  f;
	mutable gsFunctionExpr<codi::RealReverse>  gf;
	gsFunctionExpr<codi::RealReverse>  coeff_A;
	gsFunctionExpr<codi::RealReverse>  coeff_b;
	gsFunctionExpr<codi::RealReverse>  coeff_c;

	gsBoundaryConditions<codi::RealReverse> BCs;

    void compute( const gsAsConstVector<real_t> & u, bool flag_grad) const
    {
			//Creating tape inputs from active control points
			codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
			gsVector<codi::RealReverse> ctrlpoints (m_numDesignVars);
			for (int i=0; i<m_numDesignVars;i++)
			{
				ctrlpoints(i).setValue(u(i));
			}
			if (flag_grad)
			{
				tape.setActive();
				for (int i=0; i<m_numDesignVars;i++)
				{
					tape.registerInput(ctrlpoints(i));
				}
			}
			//Reconstruction of the geometry from control points
			int size_temp = 0;
			for(index_t ip=0; ip < mp.nPatches();ip++)
			{
				gsGeometry<codi::RealReverse>* pGeom = &mp.patch(ip);
				for (int i=0; i<pGeom->coefs().rows();i++)
				{
					for (int j=0; j<pGeom->coefs().cols();j++)
						pGeom->coefs()(i,j) = ctrlpoints(size_temp + i*pGeom->coefs().cols()+j);
				}
				size_temp = size_temp + pGeom->coefs().size();
			}
			gsSpringPatch<codi::RealReverse> spring(mp);
			if (info_dom)
				gsInfo<<"Created a " << spring.compute() <<"\n";
			else
				spring.compute();
			gsMultiPatch<codi::RealReverse> mp_full = spring.result();

			gsMultiBasis<codi::RealReverse> bases(mp_full);
			if (basisDegree)
				bases.setDegree(basisDegree);
			else if (numElevate)
				bases.degreeElevate(numElevate);
			for (int i=0; i<numHref ; i++)
				bases.uniformRefine();

			// Solver
			gsSparseSolver<codi::RealReverse>::LU solver;

			gsCDRAssembler<codi::RealReverse> galerkin(mp_full,bases,BCs,f,coeff_A,coeff_b,coeff_c,
				dirichlet::elimination, iFace::glue, stabilizerCDR::none);

			galerkin.assemble();
			gsMatrix<codi::RealReverse> solVector;
			//Efficient solving of the system
			solveSystem(solver, galerkin.matrix(), galerkin.rhs(), solVector);
			//Ineficient way
			//solver.compute( galerkin.matrix() );
			//solVector = solver.solve( galerkin.rhs() );

			gsField<codi::RealReverse> sol = galerkin.constructSolution(solVector);

			// Post-processing
			gsConstantFunction<codi::RealReverse> zero (0.0,2);
			codi::RealReverse result = sol.distanceL2(zero,true);

			result = result * result;
			//gsNormL2Boundary<codi::RealReverse> nor = gsNormL2Boundary<codi::RealReverse>(sol,zero,false);
			//codi::RealReverse circ =  nor.compute()*nor.compute();
			gsExprEvaluator<codi::RealReverse> ev;
			// Set the parameter mesh as the integration mesh
            ev.setIntegrationElements(bases);
			// Define integrant variables
			typedef gsExprEvaluator<codi::RealReverse>::geometryMap geometryMap;
			typedef gsExprEvaluator<codi::RealReverse>::variable variable;
			geometryMap G = ev.getMap(mp);
			variable uu = ev.getVariable(zero, G);
			variable s = ev.getVariable(sol.function());
			codi::RealReverse circ = ev.integralBdr( (uu - s).sqNorm() * nv(G).norm() );
			codi::RealReverse result2 = -result+ penalty*(circ-circ0)*(circ-circ0)*(circ-circ0)*(circ-circ0);
			std::cout << circ << std::endl;
			if (flag_grad)
			{
				tape.registerOutput(result2);
				if (info_diff)
					tape.printStatistics();
				tape.setPassive();
				result2.setGradient(1.0);
				tape.evaluate();

				for (int i=0; i<m_numDesignVars;i++)
					p_resultGrad(i) = ctrlpoints(i).getGradient();

				iter++;
			}

			p_result = result2.getValue();

			if(plot && flag_grad)
			{
				std::string fileName = "STEP" + util::to_string(iter);
				gsWriteParaview<>( sol, fileName, 1000);
				for (index_t ptc = 0; ptc < mp_full.nPatches(); ptc++)
				{
					collection->addTimestep(fileName,ptc,iter,".vts");
				}
			}
			if(plot && flag_grad)
			{
			gsInfo<<"Plotting in Paraview...\n";
			gsWriteParaview<>( mp_full, "FINAL", 1000);
			}
    }
};

} // end namespace gismo
