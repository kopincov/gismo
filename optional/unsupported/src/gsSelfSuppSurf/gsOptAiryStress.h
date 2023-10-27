
#pragma once

#include <gsIpopt/gsOptProblem.h>


namespace gismo
{



/** 
 * @brief Optimize given vault structure by transforming it to
 self-supporting. This is done by finding a convex stess potential and
 a near-by surface that admits it.

 This is done by minimizing the residue

 \f$ ||  u - u_0   ||^2   + r  || \phi ||^2 \f$

 under the constraints

 - \f$ -\Delta u = area(u) \f$ is a convex function

 - \f$ \phi \f$ is a convex function

 - r is a regularization parameter

 The first constraint can be imposed weekly by an isogeometric
 Galerkin discretization.

 The second constraint is equivalent to 

 \f$ \det H_\phi  \geq 0 \f$

 \f$ trace(  H_\phi  ) \geq 0 \f$
  
 and a sufficient condition is 
  
 \f$c_i \geq  0\f$
   
 \f$t_i \geq  0\f$

 where

 \f$\det H_\phi  =  \sum_i   c_i * \hat N_i\f$


 \f$trace(  H_\phi  ) \sum_j   t_j * \hat N_j\f$
   
*/
template <typename T>
class gsOptAiryStress : public gsOptProblem<T>
{
public:

    gsOptAiryStress(const gsMultiPatch<T> & patches, 
                    const gsBoundaryConditions<> & bcs,
                    const gsFunction<T> & f,
                    const gsPointLoads<T> & pLoads, 
                    T reg = 1e-2)
    : m_reg(reg)
    {
        // 
        gsMultiPatch<T> projPatches = patches;
        for (index_t i = 0; i<projPatches.nPatches(); ++i)
            projPatches.patch(i).embed(2);

        // Global airy potential
        gsMatrix<T> bbox;        
        projPatches.boundingBox(bbox);
        /*
          m_airBasis = dynamic_cast<gsTensorBSplineBasis<2,T>&>(projPatches.patch(0).basis());
          m_airBasis.knots(0).transform(bbox(0,0),bbox(0,1));
          m_airBasis.knots(1).transform(bbox(1,0),bbox(1,1));        

          // Define bases for determinant and trace
          m_trBasis  = m_airBasis; // trace has cont. reduced by 2
          m_trBasis.reduceContinuity(2);
          m_detBasis = m_trBasis;  // det has cont. reduced by 2 and double polynomial degree
          m_detBasis.component(0).degreeElevate( m_detBasis.degree(0) );
          m_detBasis.component(1).degreeElevate( m_detBasis.degree(1) );

          m_detPoints = m_detBasis.anchors(); // Greville points for interpolation
          m_trPoints  = m_trBasis .anchors();
        
          // Factorizations of col matrices
          gsSparseMatrix<T> colMat;
          m_detBasis .collocationMatrix(m_detPoints, colMat);
          m_detFact  .analyzePattern(colMat);
          m_detFact  .factorize(colMat);
          m_trBasis  .collocationMatrix(m_trPoints, colMat);
          m_trFact   .analyzePattern(colMat);
          m_trFact   .factorize(colMat);
        */

        gsVector<T> v(4);
        v << 1, 0, 0, 1;
        m_curAiry = gsConstantFunction<T>(v, 2);
        m_assembler = new gsMasonry<T>(projPatches, bcs, f, m_curAiry, pLoads);

        computeTargetSurfDofs(patches, 2, m_tarSurf);

        setupProblem();
    }

    ~gsOptAiryStress()
    {
        delete m_assembler;
    }

protected:

    void setupProblem()
    {
        // Our design variables are the coefficients of the Airy stress and the surface
        m_surfDofs = m_assembler->dofMapper().freeSize();
        m_airyDofs = 4;
        m_numDesignVars = m_surfDofs + m_airyDofs;
        
        computeInitialDesign();
        
        // Set number of constraints to as many as the pde dofs plus
        // determinant rep. plus trace rep.
        m_numConstraints = m_surfDofs + 2;
        
        m_desLowerBounds.setConstant(m_numDesignVars, -1.0e3);
        m_desUpperBounds.setConstant(m_numDesignVars,  1.0e3);

        // NOTE: if lower and upper bounds are set to the same number
        // then we have equality constraints.
        m_conLowerBounds.resize(m_numConstraints);
        m_conUpperBounds.resize(m_numConstraints);

        // Constraints set 1: PDE
        m_conLowerBounds.topRows(m_surfDofs).setConstant( -1e-9 );
        m_conUpperBounds.topRows(m_surfDofs).setConstant(  1e-9 );

        // Constraints set 2: Airy convexity
        m_conLowerBounds.bottomRows(m_numConstraints-m_surfDofs).setConstant(0.0   );
        m_conUpperBounds.bottomRows(m_numConstraints-m_surfDofs).setConstant( 1.0e3);

        // Dense Jacobian for now
        this->computeJacStructure();

        //gsDebugVar( m_airBasis );
        //gsDebugVar( m_detBasis );
        //gsDebugVar( m_trBasis  );
        //gsDebugVar(m_numDesignVars );
        //gsDebugVar(m_numConstraints);
    }


public:

    // Evaluate the residue
    T evalObj( const gsAsConstVector<T> & u ) const
    {
        return (u.topRows(m_surfDofs) - m_tarSurf).squaredNorm() //;
        + m_reg * u.bottomRows(m_airyDofs).squaredNorm() ;
    }

    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const  
    {
        //result.resize( m_numDesignVars );
        result.topRows   (m_surfDofs) =  2 * ( u.topRows(m_surfDofs) - m_tarSurf );
        result.bottomRows(m_airyDofs) = (2 * m_reg) * u.bottomRows(m_airyDofs);
        //result.bottomRows(m_airyDofs).setZero();
    }

    // Constraints
    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        //result.resize( m_numConstraints );

        m_curAiry.setValue(u.bottomRows(m_airyDofs));

        const typename gsAsConstVector<T>::ConstRowsBlockXpr surf = u.topRows(m_surfDofs);

        gsMultiPatch<T>  curr;
        m_assembler->constructSolution(surf, curr);
        m_assembler->assembleSystem(curr);
        
        // 1 -- Pde constraints
        result.head(m_surfDofs) = m_assembler->matrix() * surf -  m_assembler->rhs();
        
        // 2 -- Airy convexity constraints
        result(m_surfDofs  ) = m_curAiry.value().reshapeCol(0,2,2).determinant();
        result(m_surfDofs+1) = m_curAiry.value().reshapeCol(0,2,2).trace();
    
    }

    // not used yet
    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        m_curAiry.setValue(u.bottomRows(m_airyDofs)); 
            
        //1. PDE
            
        //2. Convexity
            
        // Determinant of Hessian matrix
            
        // Trace of Hessian matrix
            
        //gsDebugVar( btr.transpose() );
            
        gsDebug<< "CALLED\n";
    }
        
    T targetResidue() const
    {
        return (m_curDesign.topRows(m_surfDofs) - m_tarSurf).squaredNorm();
    }
        
    void computeTargetSurfDofs(const gsMultiPatch<T> & mp,
                               index_t c, // component to pick
                               gsMatrix<T> & result)
    {
        const gsDofMapper & mapper = m_assembler->dofMapper();
        result.resize( mapper.freeSize(), 1);
            
        for (index_t k = 0; k <  mp.nPatches(); ++k)
            for (unsigned i = 0; i <  mp.patch(k).coefsSize(); ++i)
            {
                const index_t ii = mapper.index(i,k); 
                if ( mapper.is_free_index(ii) )
                    result(ii,0) = mp.patch(k).coef(i,c);
            }
    }
        
    void computeInitialDesign()
    {
        m_curDesign.resize(m_numDesignVars,1);
        
        // Initial Airy stress potential
        /*
          gsMatrix<T> tmp = m_airBasis.anchors();
          gsFunctionExpr<T> initialPhi(" 0.5 * ( x^2 + y^2)", 2);
          m_curAiry = m_airBasis.interpolateAtAnchors(initialPhi.eval(tmp));
        */
        
        m_curDesign.bottomRows(m_airyDofs) = m_curAiry.value();
        
        //m_assembler->computeDirichletDofs();
        gsInfo << ( m_assembler->dirValues().transpose() );

        gsNewtonIterator<T> newton(*m_assembler);
        newton.setMaxIterations(60);
        newton.solve();

        gsMatrix<T> tmp;
        computeTargetSurfDofs(newton.solution(), 0, tmp);

        /*
        // Verify the residue
        //m_assembler->assembleSystem( newton.solution() );
        gsDebug<< "Residue Norm: " << 
        ( m_assembler->matrix() * tmp -  m_assembler->rhs() ).norm() <<"\n";
        gsDebug<< "Target dist: " << 
        ( m_assembler->matrix() * m_tarSurf -  m_assembler->rhs() ).norm() <<"\n";
        */

        m_curDesign.topRows(m_surfDofs).swap(tmp);
    }

    // Reconstruct solution from computed solution vector
    void constructSurface(gsMultiPatch<T> & result)
    {
        m_assembler->constructSolution(m_curDesign.topRows(m_surfDofs), result);

        for (index_t k = 0; k <  result.nPatches(); ++k)
        {
            gsGeometry<T> & p = result.patch(k);
            p.embed3d();
            p.coefs().col(0).swap( p.coefs().col(2) );
            p.coefs().leftCols(2) = m_assembler->patches().patch(k).coefs();
        }
    }


    const gsFunction<T> & airyStress() const {return m_curAiry;}

protected:

    T m_reg;

    gsMatrix<T> m_tarSurf;
    
    mutable gsConstantFunction<T> m_curAiry;

    index_t m_surfDofs;

    index_t m_airyDofs;

private:

    /// Basis for the Hessian determinant
    gsTensorBSplineBasis<2,T> m_detBasis;

    /// Basis for the trace of the Hessian
    gsTensorBSplineBasis<2,T> m_trBasis;

    /// Basis for the airy stress potential
    gsTensorBSplineBasis<2,T> m_airBasis;

    //gsAiryStress<T> m_assembler;
    gsMasonry<T> * m_assembler;

    /// Collocation points for determinant
    gsMatrix<T> m_detPoints;

    /// Collocation points for trace
    gsMatrix<T> m_trPoints;

    /// Linear solver
    gsSparseSolver<>::LU  m_detFact;

    /// Linear solver
    gsSparseSolver<>::LU  m_trFact;

    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;
};

} // end namespace gismo
