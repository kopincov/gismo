/** @file gsOptParam.h

    @brief Base class for a problem of optimization of geometry parameterization

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsIpopt/gsOptProblem.h>

#define JAC_STRUCTURE // undefine for full structure

namespace gismo
{

/** 
 * \brief Base class for a problem of optimization of geometry parameterization
 *
 */

template <typename T>
class gsOptParam : public gsOptProblem<T>
{
public:

    gsOptParam() : m_geoBasis(NULL), m_jacBasis(NULL) { }

 gsOptParam(const gsGeometry<T> & geo) : m_geoBasis( geo.basis().clone().release() )
    {
        //to do: test cast to bspline
        m_dim = m_geoBasis->dim();

        // Dimension of the control points
        m_geoDim = geo.geoDim();

        // Get the right basis for the Jacobian
        // Note: degree of the JacDet is d*p-1 in each direction for dimension d, degree=p
        m_jacBasis =  m_geoBasis->clone().release();
        m_jacBasis->reduceContinuity(1);

        for ( int k = 0; k<m_dim; k++)
            m_jacBasis->component(k).degreeElevate( (m_dim-1) * m_geoBasis->degree(k) - 1);
        
        // Store corner coordinates
        m_corner.resize(1<<m_dim);
        for (index_t i = 0; i!= m_corner.size(); ++i)
            m_corner[i] = m_jacBasis->functionAtCorner(i);

        gsDebugVar( (*m_geoBasis).detail() );
        gsDebugVar( (*m_jacBasis).detail() );

        // Set number of constraints to as many as the Jacobian
        // coefficients
        m_numConstraints = m_jacBasis->size();
        //gsDebugVar(m_numConstraints);

        // Get design control points, fixing the boundary control
        // points
        m_mapper = gsDofMapper(*m_geoBasis);
        // boundary control point indices
        gsMatrix<index_t> bnd = m_geoBasis->allBoundary();
        m_mapper.markBoundary(0,bnd);
        m_mapper.finalize();

        GISMO_ENSURE( m_mapper.freeSize() > 0, "No design variables!");

        // Number of design variables (inner control point coordinates)
        m_numDesignVars = m_geoDim * m_mapper.freeSize();
        //gsDebugVar(m_numDesignVars );

        // Compute and factorize the collocation matrix. This is used
        // to evaluate the Geometry's Jacobian coefficients (and their
        // derivatives w.r.t. the design variables
        // TO DO: in the case of tensor-product, exploit the structure
        // by solving d small linear systems instead
        m_jacBasis->anchors_into(m_colPoints);
        gsSparseMatrix<T> colMat;
        m_jacBasis->collocationMatrix(m_colPoints, colMat);
        m_jacSolver.analyzePattern(colMat);
        m_jacSolver.factorize(colMat);
        //gsDebug<<"Factorization OK\n";

        // Compute the Jacobian structure
        //gsDebugVar(m_numConstraints * m_numDesignVars);
        #ifdef JAC_STRUCTURE
        this->computeJacStructure();// non-virtual call
        #else
        gsOptProblem<T>::computeJacStructure();// dense version
#       endif

        // Set boundary coefficients equal to the input once and
        // forever
        const gsMatrix<T> & coefs = geo.coefs();
        m_tmpCoefs.resize( m_geoBasis->size(), m_geoDim);
        for (index_t i = 0; i < bnd.rows(); ++i) // copy to m_tmpCoefs
            m_tmpCoefs.row( bnd(i,0) ) = coefs.row( bnd(i,0) );

        // Design bounds
        m_desLowerBounds.setConstant(m_numDesignVars, -1.0e19);
        m_desUpperBounds.setConstant(m_numDesignVars,  1.0e19);
        
        // Constraint bounds
        m_conLowerBounds.setConstant(m_numConstraints,  1.0e-3);
        m_conUpperBounds.setConstant(m_numConstraints,  1.0e19);

        // Current design -- used by gsOptProblem as starting point
        m_curDesign.resize(m_numDesignVars, 1);
        // view current design as a control point matrix
        gsAsMatrix<T> curDesign( m_curDesign.data(), m_mapper.freeSize(), m_geoDim);
        for ( int i = 0; i< m_tmpCoefs.rows(); i++ ) // TO DO: more efficient
        {
            // get index of the design control point
            const index_t ii = m_mapper.index(i,0);
            
            if ( m_mapper.is_free_index(ii) ) // interior node?
                curDesign.row(ii) = coefs.row(i);
        }
      
        gsDebugVar(m_numConJacNonZero);
        gsDebugVar( m_mapper.freeSize() );

    }

    ~gsOptParam()
    {
        delete m_geoBasis;
        delete m_jacBasis ;
    }
    
public:

    // Look at gsOptProblem for a description of this function
    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        // get the coefficients
        jacobianDetCoefs(u, result);

        // Disactivate the corners, since they have constant coefficient
        for (index_t i = 0; i!= m_corner.size(); ++i)
            result[m_corner[i]] = T(1.0);
    }
    
    // Look at gsOptProblem for a description of this function
    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        jacobianDetSensitivities(u, result);
    }


    const gsMatrix<T> & currentCoefs() const 
    {
        gsAsConstVector<T> curDes (m_curDesign.data(), m_numDesignVars);
        updateTempDesign(curDes);// todo: add constructor in gsMatrix for implicit cast
        return m_tmpCoefs; 
    }

    void currentJacobianCoefs(gsVector<T> & result) const 
    {
        result.resize(m_numConstraints);
        gsAsConstVector<T> curDes(m_curDesign.data(), m_numDesignVars);
        gsAsVector<T>      res   (result.data()     , m_numConstraints);
        
        jacobianDetCoefs(curDes, res);
    }

    const gsBasis<T> & jacBasis() const 
    {
        return *m_jacBasis; 
    }

protected:

    inline void jacobianDetCoefs(const gsAsConstVector<T> & u, 
                                 gsAsVector<T> & result) const
    {
        updateTempDesign(u);
        // Evaluate right hand side
        gsVector<T> b(m_colPoints.cols());
        gsMatrix<T> tmpJac;
        for ( index_t i = 0; i != m_colPoints.cols(); ++i)
        {
            m_geoBasis->jacobianFunc_into(m_colPoints.col(i), m_tmpCoefs, tmpJac);
            b[i] = tmpJac.determinant();
        }

        // Solve for the constraints
        // note: result is expected to be allocated at least
        // m_numConstraints numbers
        result.topRows(m_colPoints.cols()) = m_jacSolver.solve(b);
    }

    inline void jacobianDetSensitivities(const gsAsConstVector<T> & u, 
                                         gsAsVector<T> & result) const
    {
        updateTempDesign(u);

        #ifdef JAC_STRUCTURE

        gsMatrix<T> JacAdj(m_dim,m_dim); //m_geoDim == m_dim
        gsMatrix<T> gradB;

        const index_t sz = m_mapper.freeSize();
 
        // Space for right-and side and solution vector
        gsMatrix<T> currDV(m_colPoints.cols(), m_geoDim),
                         b(m_colPoints.cols(), m_geoDim);
        
        for ( index_t k = 0; k != m_geoBasis->size(); ++k) // for all control points
        {
            // get index of the design control point
            const index_t kk = m_mapper.index(k, 0);
            
            if ( m_mapper.is_free_index(kk) ) // interior ?
            {
                for ( index_t i = 0; i != m_colPoints.cols(); ++i) // for all constraints on detJ
                {
                    // Skip corner constraints
                    if ( (m_corner.array() == i).any() )
                    {
                        b.row(i).setZero();
                        continue; // for
                    }
   
                    if ( m_geoBasis->isActive(k, m_colPoints.col(i)) )
                    {
                        // Get Jacobian
                        m_geoBasis->jacobianFunc_into(m_colPoints.col(i), m_tmpCoefs, JacAdj);

                        // Adjugate (cofactor transpose)
                        JacAdj.adjugateInPlace();
                        
                        // Gradient of B_k
                        m_geoBasis->derivSingle_into(k, m_colPoints.col(i), gradB);
                        
                        for ( index_t s = 0; s != m_geoDim; ++s) // for all components
                            b(i, s) = ( JacAdj.col(s) * gradB.transpose() ).trace();
                    }
                }

                // Solve for the kk-th sensitivities
                currDV = m_jacSolver.solve(b);
                
                // Fill in sparse result
                for (index_t s = 0; s != m_geoDim; ++s) // for all components
                {
                    const index_t rr = s*sz + kk;

                    for (index_t m = m_marker[rr]; m != m_marker[rr+1]; ++m)// for all entries
                        result[m] = currDV(m_conJacRows[m],s);
                }
            }
        }

        #else

        gsMatrix<T> JacAdj(m_dim,m_dim); //m_geoDim == m_dim
        gsMatrix<T> gradB;
        gsMatrix<index_t> activeB;

        const index_t sz = m_mapper.freeSize();

        // Initialize right-hand side to zero
        gsMatrix<T> b = gsMatrix<T>::Zero(m_colPoints.cols(), m_geoDim*m_mapper.freeSize());

        // Evaluate right hand side b
        for ( index_t i = 0; i != m_colPoints.cols(); ++i)
        {
            // m_geoBasis->compute(active, deriv, jacob, m_tmpCoefs);
            // gsBasis<T>::linearComb(active, evals, m_tmpCoefs, result);
            // gsBasis<T>::jacobianFromGradients(active, grads, m_tmpCoefs, result);

            // Get Jacobian determinant and inverse
            m_geoBasis->jacobianFunc_into(m_colPoints.col(i), m_tmpCoefs, JacAdj);
            //Adjugate (cofactor transpose)
            JacAdj.adjugateInPlace();

            // Active basis function at point i
            //const int numActive = 
            m_geoBasis->active_into(m_colPoints.col(i), activeB);
            const int numActive = activeB.rows();

            // Basis Gradients matrix
            m_geoBasis->deriv_into(m_colPoints.col(i), gradB);

            for ( index_t k = 0; k != numActive; ++k) // For all active
            {
                // get index of the design control point index
                const index_t kk = m_mapper.index(activeB(k,0),0);
                
                if ( m_mapper.is_free_index(kk) ) // interior node?
                {
                    for ( index_t s = 0; s != m_geoDim; ++s) // for all components
                    {
                        b(i, s*sz + kk) =
                            // JacDet*( JacInv.col(s) * gradB.template block<d,1>(k*d,0).transpose() ).trace();
                            ( JacAdj.col(s) * gradB.block(k*m_dim,0,m_dim,1).transpose() ).trace();
                    }
                }
            }
        }

        // Solve for the Constraint Jacobian (sensitivities)
        gsAsMatrix<T> JC(result.data(), m_numConstraints, b.cols() );
        // Note: result may contain extra data
        // ( ie. m_numConstraints>b.cols() ), so we fill the
        // principal sub-matrix
        JC.block(0,0, b.rows(),b.cols() ) = m_jacSolver.solve(b);

        // Disactivate the corners by setting their sensitivities to 0
        for (index_t i = 0; i!= m_corner.size(); ++i)
            //JC(m_corner[i], m_numDesignVars-1) = 0.0;
            JC.row(m_corner[i]).setZero();

#endif//JAC_STRUCTURE
    }


    void updateTempDesign(const gsAsConstVector<T> & u) const // ..
    {
        // Reorder design variable values
        gsAsConstMatrix<T> u_mat( u.data(), m_mapper.freeSize(), m_geoDim);

        // Update interior coefficients
        for ( index_t i = 0; i< m_tmpCoefs.rows(); i++ ) // TO DO: more efficient
        {
            const index_t ii = m_mapper.index(i,0);
            
            if ( m_mapper.is_free_index(ii) ) // interior node?
                m_tmpCoefs.row(i) = u_mat.row(ii);
        }
    }

protected:


    void computeJacStructure()
    {
        //to try: typedef gsSortedVector<std::pair<index_t,index_t> > pairList;
        typedef std::set<std::pair<index_t,index_t> > pairList;
        gsMatrix<index_t> actRow, actCol;
        pairList pairs; // pair: (DesVar,Constr)
        
        typename gsBasis<T>::domainIter domIt = m_geoBasis->makeDomainIterator();
        int sz = m_mapper.freeSize();

        for (; domIt->good(); domIt->next() )
        {
            m_jacBasis->active_into(domIt->center, actRow);
            m_geoBasis->active_into(domIt->center, actCol);

            for (index_t j = 0; j!= actCol.rows(); ++j)
            {
                const index_t jj = m_mapper.index(actCol(j),0);

                if ( m_mapper.is_free_index(jj) )
                    for (index_t i = 0; i!= actRow.rows(); ++i)
                        for ( index_t s = 0; s != m_geoDim; ++s) // for all components
                            pairs.insert( std::make_pair(s*sz + jj, actRow(i) ) );
            }
        }

        // sz now bacomes the size of the constraint jacobian related to detJ
        sz = pairs.size();
        m_conJacRows.reserve(sz);
        m_conJacCols.reserve(sz);
        
        // the marker encodes the columns: we don't need really need
        // m_conJacCols anymore, due to sorting wrt design vars!
        m_marker.setZero(m_numDesignVars+1);

        for (pairList::const_iterator it = pairs.begin(); it != pairs.end(); ++it)
        {
            if ( (m_corner.array() == it->second).any() )
                continue;

            m_conJacRows.push_back(it->second);//constraint
            m_conJacCols.push_back(it->first );//design var

            ++m_marker[it->first+1];
        }

        // Complete the marking
        for (index_t k = 1; k<= m_numDesignVars; ++k)
            m_marker[k] += m_marker[k-1];

        m_numConJacNonZero = m_conJacRows.size();
        
        // gsDebugVar( gsAsVector<index_t>(m_conJacRows).transpose());
        // gsDebugVar( gsAsVector<index_t>(m_conJacCols).transpose());
        //gsDebugVar( m_marker.transpose());

/*
        //More sophisticated impl. ?
        gsTensorBSplineBasis<2,T> * geoBasis =  
            dynamic_cast<gsTensorBSplineBasis<2,T>*>(m_geoBasis);
        gsTensorBSplineBasis<2,T> * jacBasis =  
            dynamic_cast<gsTensorBSplineBasis<2,T>*>(m_jacBasis);

        gsMatrix<index_t,-1,2> sup(m_dim,2);
        gsMatrix<index_t> act;

        const index_t sz = m_jacBasis->size();
        gsVector<index_t> perm = 
            gsVector<index_t>::LinSpaced(sz,0,sz-1);
        index_t c = m_colPoints.cols() - (1<<m_dim) - 1;
        for (int i = 0; i!= (1<<m_dim); ++i)
            perm.row(m_jacBasis->functionAtCorner(i))
                .swap( perm.row(c++) );
        c = sz - (1<<m_dim);
    
        for ( index_t i = 0; i != c; ++i)
        {
            // Get influence area of i-th constraint
            jacBasis->elementSupport_into(perm[i], sup);
            
            // Get all active design variables in the area
            geoBasis->elementActive_into(sup, act);

            for ( index_t j = 0; j < act.rows(); ++j)
            {
                const index_t jj = m_mapper.index(act(j,0), 0);

                if ( m_mapper.is_free_index(jj) )
                {
                    m_conJacRows.push_back( perm[i] );
                    m_conJacCols.push_back( jj      );
                }
            }
            // dummy variable
            m_conJacRows.push_back( perm[i] );
            m_conJacCols.push_back( m_numDesignVars-1 );
        }

        m_numConJacNonZero = m_conJacRows.size();

        gsDebugVar( gsAsVector<index_t>(m_conJacRows).transpose());
        gsDebugVar( gsAsVector<index_t>(m_conJacCols).transpose());
        //*/

    }


protected:

    int m_dim;

    /// Dimension of the ambient space
    int m_geoDim;

    /// Collocation points for Jacobian determinant
    gsMatrix<T> m_colPoints; // to do: maybe the basis values on those ?

    /// Mapper from control points to design variables
    gsDofMapper m_mapper;

    /// Basis of the input parameterization
    //gsTensorBSplineBasis<d,T> m_geoBasis;
    gsBasis<T> * m_geoBasis;

    /// Basis for the Jacobian determinant
    //gsTensorBSplineBasis<d,T> m_jacBasis;
    gsBasis<T> * m_jacBasis;

    /// Keeps the indices of the corner coefficients
    gsVector<T> m_corner;

    gsVector<index_t> m_marker;

    /// Sparse factorization of the collocation matrix
    typename gsSparseSolver<T>::LU  m_jacSolver;
    //typename gsSparseSolver<T>::QR  m_jacSolver;// slow
    //typename gsSparseSolver<T>::SuperLU  m_jacSolver;
    //typename gsSparseSolver<T>::SimplicialLDLT  m_jacSolver;
    
    /// Temporary used to translate back designs as control points
    mutable gsMatrix<T> m_tmpCoefs;

protected: // Inherited members from gsOptProblem

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
