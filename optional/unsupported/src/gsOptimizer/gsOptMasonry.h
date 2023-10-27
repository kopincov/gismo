
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
class gsOptMasonry : public gsOptProblem<T>
{
public:

    gsOptMasonry(const gsMultiPatch<T> & patches, 
                 const gsBoundaryConditions<> & bcs,
                 const gsFunction<T> & f,
                 const gsPointLoads<T> & pLoads, 
                 T reg = 1e-2)
    : m_reg(reg)
    {
        // 
        gsMultiPatch<T> projPatches = patches;
        for (size_t i = 0; i< projPatches.nPatches(); ++i)
            projPatches.patch(i).embed(2);

         // Global airy potential
        gsMatrix<T> bbox;        
        projPatches.boundingBox(bbox);
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

        m_assembler = new gsMasonry<T>(projPatches, bcs, f, m_phi, pLoads);

        computeTargetSurfDofs(patches, 2, m_tarSurf);

        setupProblem();
    }

    ~gsOptMasonry()
    {
        delete m_curAiry; 
        delete m_assembler;
    }

protected:

    void setupProblem()
    {
        // Our design variables are the coefficients of the Airy stress and the surface
        m_surfDofs = m_assembler->system().colMapper(0).freeSize();
        m_airyDofs = m_airBasis.size();
        m_numDesignVars = m_surfDofs + m_airyDofs;
        
        computeInitialDesign();
        
        // Set number of constraints to as many as the pde dofs plus
        // determinant rep. plus trace rep.
        m_numConstraints = m_surfDofs + m_trBasis.size() + m_detBasis.size();
        
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
        return (u.topRows(m_surfDofs) - m_tarSurf).squaredNorm()
                + m_reg * u.bottomRows(m_airyDofs).squaredNorm() ;
    }

    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const  
    {
        //result.resize( m_numDesignVars );
        result.topRows   (m_surfDofs) =  2 * ( u.topRows(m_surfDofs) - m_tarSurf );
        result.bottomRows(m_airyDofs) = (2 * m_reg) * u.bottomRows(m_airyDofs);
    }

    // Constraints
    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        //result.resize( m_numConstraints );

        m_curAiry->coefs() = u.bottomRows(m_airyDofs);

        const typename gsAsConstVector<T>::ConstRowsBlockXpr surf = u.topRows(m_surfDofs);

        gsMultiPatch<T>  curr;
        m_assembler->constructSolution(surf, curr);
        m_assembler->assembleSystem(curr);
        
        // 1 -- Pde constraints
        result.head(m_surfDofs) = m_assembler->matrix() * surf -  m_assembler->rhs();
        
        // 2 -- Airy convexity constraints
        gsVector<T> bdet( m_detBasis.size() );
        gsVector<T> btr ( m_trBasis .size() );
        gsMatrix<T> tmpHes;
        
        for ( index_t i = 0; i != m_detBasis.size(); ++i)
        {
            // evaluation of determinant
            m_curAiry->deriv2_into(m_detPoints.col(i), tmpHes);
            bdet[i] = tmpHes(0,0) * tmpHes(1,0) - tmpHes(2,0)*tmpHes(2,0);
        }

        for ( index_t i = 0; i != m_trBasis.size(); ++i)
        {
            // evaluation of trace
            m_curAiry->deriv2_into(m_trPoints.col(i), tmpHes);
            btr [i] = tmpHes(0,0) + tmpHes(1,0);
        }

        // Solve for the constraints
        result.middleRows(m_surfDofs, bdet.size() ) = m_detFact.solve(bdet);
        result.tail( btr.size() ) = m_trFact.solve(btr);

        // gsDebugVar(bdet.transpose() );
        // gsDebugVar(btr .transpose() );
    }
    
    // not used yet
    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        gsMatrix<T> tmpHes, der2N;
        gsMatrix<index_t> activeB;

        m_curAiry->coefs() = u.bottomRows(m_airyDofs);
        
        //1. PDE

        //2. Convexity
        const index_t detSz = m_detBasis.size();
        const index_t  trSz = m_trBasis .size();

        // Initialize right-hand side to zero
        gsMatrix<T> bdet(detSz, m_detPoints.cols() ),
                    btr (trSz , m_trPoints .cols() );

        // Determinant of Hessian matrix
        for ( index_t i = 0; i != detSz; ++i)
        {
            // Get Hessian
            m_curAiry->deriv2_into(m_detPoints.col(i), tmpHes);

            // Active basis function at point i
            m_airBasis.active_into(m_detPoints.col(i), activeB);
            const int numActive = activeB.rows();

            // Basis second derivative matrix
            m_airBasis.deriv2_into(m_detPoints.col(i), der2N);

            for ( index_t k = 0; k != numActive; ++k) // For all active
            {
                // get index of the design control point index
                const index_t kk = activeB(k,0); //m_mapper.index(activeB(k,0),0);
                
                //if ( m_mapper.is_free_index(kk) ) // interior node?
                {
                    bdet(i,kk) =      der2N(3*k+1, 0) * tmpHes(0,0) 
                               +      der2N(3*k  , 0) * tmpHes(1,0) 
                               -  2 * der2N(3*k+2, 0) * tmpHes(2,0);
                }
            }
        }
        
        // Trace of Hessian matrix
        for ( index_t i = 0; i != trSz; ++i)
        {
            m_airBasis.active_into(m_trPoints.col(i), activeB);
            const int numActive = activeB.rows();
            // Basis Gradients matrix
            m_airBasis.deriv2_into(m_trPoints.col(i), der2N);
            for ( index_t k = 0; k != numActive; ++k) // For all active
            {
                // get index of the design control point index
                const index_t kk = activeB(k,0); //m_mapper.index(activeB(k,0),0);
                
                //if ( m_mapper.is_free_index(kk) ) // interior node?
                {
                    btr(i,kk) = der2N(3*k, 0) + der2N(3*k+1, 0);
                }
            }
        }

        //gsDebugVar( btr.transpose() );

        // Solve for the Constraint Jacobian (sensitivities)
        gsMatrix<T> sensDet = m_detFact.solve(bdet);
        gsMatrix<T> sensTr  = m_trFact .solve(btr );
        const index_t sz1 = sensDet.size();
        const index_t sz2 = sensTr .size();
        //result.resize( sz1 + sz2, 1);

        gsDebugVar( m_numConstraints );
        gsDebugVar( m_numDesignVars  );
        gsDebugVar( m_numConJacNonZero );
        gsDebugVar( result.rows() );
        gsDebugVar( result.cols() );
        gsDebugVar( sz1 );
        gsDebugVar( sz2 );

        result.head(sz1) = sensDet.asVector();// !!!!!
        result.tail(sz2) = sensTr .asVector();
        
        //gsDebugVar( result.transpose() );

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
        const gsDofMapper & mapper = m_assembler->system().colMapper(0);
        result.resize( mapper.freeSize(), 1);
        
        for (size_t k = 0; k <  mp.nPatches(); ++k)
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
        gsMatrix<T> tmp = m_airBasis.anchors();

        // TO DO: move to gsNurbsCreator
        gsFunctionExpr<T> initialPhi(" 0.5 * ( x^2 + y^2)", 2);
        m_curAiry = m_airBasis.interpolateAtAnchors(initialPhi.eval(tmp)).release();
        m_phi.setFunction(*m_curAiry);
        
        m_curDesign.bottomRows(m_airyDofs) = m_curAiry->coefs();

        gsNewtonIterator<T> newton(*m_assembler);
        newton.setMaxIterations(60);
        newton.solve();
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

        for (size_t k = 0; k <  result.nPatches(); ++k)
        {
            gsGeometry<T> & p = result.patch(k);
            p.embed3d();
            p.coefs().col(0).swap( p.coefs().col(2) );
            p.coefs().leftCols(2) = m_assembler->patches().patch(k).coefs();
        }
    }


    gsGeometry<T> & airyStress() {return *m_curAiry;}

protected:

    T m_reg;

    gsMatrix<T> m_tarSurf;
    
    mutable gsAiryStressCoeff<T> m_phi;
    
    gsGeometry<T> * m_curAiry;

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
