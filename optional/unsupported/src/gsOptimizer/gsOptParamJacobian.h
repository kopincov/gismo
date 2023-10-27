/** @file gsOptParamJacobian.h

    @brief Optimization of geometry parameterization: singulatity
    obiviation by Jacobian positivity constraints

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsOptimizer/gsOptParam.h>

#include <gsNurbs/gsTensorBSplineBasis.h>


namespace gismo
{

/** 
    \brief Optimization of geometry parameterization: singulatity
    obiviation by Jacobian positivity constraints
 */

template <typename T>
class gsOptParamJacobian : public gsOptParam<T>
{
public:
    typedef gsOptParam<T> Base;
public:

    gsOptParamJacobian(const gsGeometry<T> & geo)
    {
        m_geoBasis = geo.basis().clone().release();// not calling Base constructor!

        //to do: test cast to bspline
        m_dim = m_geoBasis->dim();

        // Shortcut to the input coefficients
        const gsMatrix<T> & coefs = geo.coefs();

        // Dimension of the control points
        m_geoDim = geo.geoDim();

        //const int sz = m_geoBasis->size() * m_geoDim ;

        // Get the right  basis for the Jacobian
        m_jacBasis =  m_geoBasis->clone().release();
        m_jacBasis->reduceContinuity(1);
        // degree of the JacDet is d*p-1 in each direction for dimension d, degree=p
        for ( int k = 0; k<m_dim; k++)
            m_jacBasis->component(k).degreeElevate( (m_dim-1) * m_geoBasis->degree(k) - 1);
            //m_jacBasis->component(k).degreeElevate(1);

        //m_jacBasis->uniformRefine();

        // Store corner indices
        m_corner.resize(1<<m_dim);
        for (index_t i = 0; i!= m_corner.size(); ++i)
            m_corner[i] = m_jacBasis->functionAtCorner(i);

        gsDebugVar(*m_geoBasis);
        gsDebugVar(*m_jacBasis);

        // Set number of constraints to as many as the Jacobian
        // coefficients
        m_numConstraints = m_jacBasis->size();
        //gsDebugVar(m_numConstraints);

        // Get design control points, fixing the boundary control
        // points
        m_mapper = gsDofMapper(*m_geoBasis);
        // boundary control point indices
        gsMatrix<unsigned> bnd = m_geoBasis->allBoundary();
        m_mapper.markBoundary(0,bnd);
        m_mapper.finalize();

        GISMO_ENSURE( m_mapper.freeSize() > 0, "No design variables!");

        // ****  Number of design variables (inner control point coordinates
        // plus a dummy)
        m_numDesignVars = m_geoDim * m_mapper.freeSize() + 1;
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

        gsDebugVar(m_numConstraints * m_numDesignVars);
        
        #ifdef JAC_STRUCTURE
        this->computeJacStructure();// non-virtual call
        #else
        gsOptProblem<T>::computeJacStructure();// dense version
#       endif

        gsDebugVar(m_numConJacNonZero );
        gsDebugVar(m_mapper.freeSize());

        // Set boundary coefficients equal to the input once and
        // forever
        m_tmpCoefs.resize( m_geoBasis->size(), m_geoDim);
        for (index_t i = 0; i < bnd.rows(); ++i)
            m_tmpCoefs.row( bnd(i,0) ) = coefs.row( bnd(i,0) );

        // Design bounds
        m_desLowerBounds.setConstant(m_numDesignVars, -1.0e19);
        m_desUpperBounds.setConstant(m_numDesignVars,  1.0e19);
        
        // Constraint bounds
        //m_conLowerBounds.setConstant(m_numConstraints,  1.0e-3);
        m_conLowerBounds.setConstant(m_numConstraints,  0.0);
        m_conUpperBounds.setConstant(m_numConstraints,  1.0e19);

        // **** Current design -- used by gsOptProblem as starting point
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

        // Set initial dummy value to JacMinCoeff
        gsVector<T> tmp;
        m_curDesign(m_numDesignVars-1,0) = T(0.0);
        this->currentJacobianCoefs(tmp);
        m_curDesign(m_numDesignVars-1,0) = tmp.minCoeff(); // Initial dummy value
        // = tmp.minCoeff() - 0.01;
    }

public:

    T evalObj( const gsAsConstVector<T> & u ) const
    {
        //updateTempDesign(u);
        return -u[m_numDesignVars-1]; // dummy is maximized (as -dummy is minimized)
    }
    
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const  
    {
        //updateTempDesign(u);
        result.topRows(m_numDesignVars-1).setZero();
        result[m_numDesignVars-1] = -1.0; 
    }
    
    // Look at gsOptProblem for a description of this function
    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        // get the coefficients
        this->jacobianDetCoefs(u, result);

        // substruct dummy variable 
        result.array() -= u(m_numDesignVars-1,0);

        // disactivate the corners, since they have constant coefficient
        for (index_t i = 0; i!= m_corner.size(); ++i)
        {
            //gsDebug<<"Corner: "<< result[m_corner[i]] + u(m_numDesignVars-1,0) <<"\n";
            result[m_corner[i]] = T(1.0);//big values spoil the result due to scaling of constraints in IPOPT
        }
    }

    void evalCon1_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        // Constraint aggregation based on Kreisselmeier���Steinhauser function

        gsVector<T>   JacCoef(m_numConstraints);
        gsAsVector<T> JacData(JacCoef.data(), m_numConstraints);
        this->jacobianDetCoefs(u, JacData);

        // substruct dummy variable 
        JacCoef.array() -= u(m_numDesignVars-1,0);

        // Disactivate the corners by setting them to a big value
        for (index_t i = 0; i!= m_corner.size(); ++i)
            JacCoef[m_corner[i]] = T(1.0);

        const T dummy = u(m_numDesignVars-1, 0);
        T val = 0.0;
        for ( index_t i = 0; i!= JacCoef.size(); ++i)
        {
            const T diff = JacCoef[i] - dummy;
            if ( diff < T(0.0) )
                val += diff * diff * diff ;
        }
        gsDebugVar(val);
        result.setConstant(val);
    }
    

    // Look at gsOptProblem for a description of this function
    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        this->jacobianDetSensitivities(u, result);

        #ifdef JAC_STRUCTURE        
        // dummy partial derivs
        result.bottomRows(m_numConstraints-m_corner.size()).setConstant(-1.0); 
        #else
        // dummy partial derivs
        gsAsVector<T> Jlast(result.bottomRows(m_numConstraints).data(), m_numConstraints);
        Jlast.setConstant(-1.0);

        // Set corners back to zero as they ought to be
       for (index_t i = 0; i!= m_corner.size(); ++i)
           Jlast[m_corner[i]] = T(0.0);
        #endif


/* // Kreisselmeier���Steinhauser constraint aggregation

   
   
 */

    }

    T dummy () const { return m_curDesign(m_numDesignVars-1,0); };

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
        m_conJacRows.reserve(sz+m_numConstraints); // ****
        m_conJacCols.reserve(sz+m_numConstraints);
        
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

        // *** dummy gradient
        for (index_t k = 0; k!= m_numConstraints; ++k)
        {
            if ( (m_corner.array() == k).any() )
                continue;

            m_conJacRows.push_back(k);
            m_conJacCols.push_back(m_numDesignVars-1);
        }
        m_marker[m_numDesignVars] = m_numConstraints - m_corner.size(); // *** dummy

        // Complete the marking
        for (index_t k = 1; k<= m_numDesignVars; ++k)
            m_marker[k] += m_marker[k-1];

        m_numConJacNonZero = m_conJacRows.size();
        
        //gsDebugVar( gsAsVector<index_t>(m_conJacRows).transpose());
        //gsDebugVar( gsAsVector<index_t>(m_conJacCols).transpose());
        //gsDebugVar( m_marker.transpose());
    }

private:

    using Base::m_dim;
    
    /// Dimension of the ambient space
    using Base::m_geoDim;

    /// Collocation points for Jacobian determinant
    using Base::m_colPoints; // to do: maybe the basis values on those ?

    /// Mapper from control points to design variables
    using Base::m_mapper;

    /// Basis of the input parameterization
    //gsTensorBSplineBasis<d,T> m_geoBasis;
    using Base::m_geoBasis;

    /// Basis for the Jacobian determinant
    //gsTensorBSplineBasis<d,T> m_jacBasis;
    using Base::m_jacBasis;

    /// Keeps the indices of the corner coefficients
    using Base::m_corner;
    using Base::m_marker;

    /// Sparse QR factorization of the collocation matrix
    using Base::m_jacSolver;
    
    /// Temporary used to translate back designs as control points
    using Base::m_tmpCoefs;

private:// Inherited members from gsOptProblem

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
