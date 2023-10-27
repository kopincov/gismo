/** @file gsOptParamLiao.h

    @brief Optimization of geometry parameterization using Liao
    functional under non-singularity constraints

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsOptimizer/gsOptParam.h>

namespace gismo
{

/** 
 * \brief
 * Optimization of Parameterization
 *
 */

template <typename T>
class gsOptParamLiao : public gsOptParam<T>
{
public:
    typedef gsOptParam<T> Base;
public:

    gsOptParamLiao() : gsOptParam<T>() 
    { 

    }

    gsOptParamLiao(const gsGeometry<T> & geo) : gsOptParam<T>(geo)
    { 
        // Setup a quadrature rule
        gsVector<index_t> numQuadNodes( m_dim );
        for (int i = 0; i < m_geoBasis->dim(); ++i)
            numQuadNodes[i] = 2 * m_jacBasis->degree(i) + 1;
        quRule = gsGaussRule<T>(numQuadNodes);
    }

    ~gsOptParamLiao() { }
    
public:

    T evalObj( const gsAsConstVector<T> & u ) const
    {
        // Liao functional: trace( (J*J^T)^2 ) = norm(J*J^T)^2
        this->updateTempDesign(u);

        gsMatrix<T> Jac(m_dim,m_dim);
        T val = 0.0;

        // Initialize domain element iteration
        typename gsBasis<T>::domainIter domIt = m_geoBasis->makeDomainIterator();
        // Start iteration over integration elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the integration element
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

            for (index_t i = 0; i!= quWeights.size(); ++i)
            {
                m_geoBasis->jacobianFunc_into(quNodes.col(i), m_tmpCoefs, Jac);
                
                val += quWeights[i] * (Jac*Jac.transpose()).squaredNorm();
                //val += quWeights[i] * (Jac*Jac.transpose()*Jac*Jac.transpose()).trace();
            }
        }
        
        return val;
    }
    
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const  
    {
        // Liao deriv: 2*trace( (J*J^T)*(d_cJ*J^T +(d_cJ*J^T)^T ) )
        // where , d_{c_ij,k}Jac = e_{k} * grad(B_{ij})^T
        // or 2*(Jac*Jac.transpose()*(dJac*Jac.transpose()+Jac*dJac.transpose())).trace();
        this->updateTempDesign(u);

        const int sz = m_mapper.freeSize();
        gsMatrix<T> Jac, dMetric(m_dim,m_dim);
        gsMatrix<T> gradB;
        gsMatrix<index_t> activeB;
        result.setZero();

        // Initialize domain element iteration
        typename gsBasis<T>::domainIter domIt = m_geoBasis->makeDomainIterator();
        // Start iteration over integration elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the integration element
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

            // Active basis functions at current element
            m_geoBasis->active_into(quNodes.col(0), activeB);
            const int numActive = activeB.rows();

            // Basis Gradients matrix
            m_geoBasis->deriv_into(quNodes, gradB);

            m_geoBasis->jacobianFunc_into(quNodes, m_tmpCoefs, Jac);

            for ( index_t k = 0; k != numActive; ++k) // For all active
            {
                // get index of the design control point index
                const index_t kk = m_mapper.index(activeB(k,0),0);
                
                if ( m_mapper.is_free_index(kk) ) // interior node?
                {
                    for ( index_t s = 0; s != m_geoDim; ++s) // for all components
                    {
                        T & val = result[s*sz + kk];
                        for (index_t i = 0; i!= quWeights.size(); ++i)// for all quadrature points
                        {
                            const typename gsMatrix<T>::ColsBlockXpr & Jaci = Jac.middleCols(i*m_dim,m_dim);
                            
                            dMetric = ( gsVector<T>::Unit(m_dim,s) * gradB.block(k*m_dim,i,m_dim,1).transpose() )
                                * Jaci.transpose();
                            
                            val += 2.0 * quWeights[i] * 
                              ((Jaci * Jaci.transpose()).array() * (dMetric+dMetric.transpose()).array() ).sum();
                            //( Jaci * Jaci.transpose() * (dMetric+dMetric.transpose()) ).trace();
                        }
                    }
                }
            }
        }// for domIt
    }

protected:

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

    /// Sparse QR factorization of the collocation matrix
    using Base::m_jacSolver;
    
    /// Temporary used to translate back designs as control points
    using Base::m_tmpCoefs;

    gsGaussRule<T> quRule;
    mutable gsMatrix<T> quNodes;
    mutable gsVector<T> quWeights;

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
