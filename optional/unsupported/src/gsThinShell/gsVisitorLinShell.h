/** @file gsVisitorLinShell.h

    @brief Element visitor for linear elasticity on thin shell.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Goyal, A. Mantzaflaris
*/

#pragma once

namespace gismo
{

/** 
    @brief The visitor computes element system matrix and right hand side for linear elasticity
            on thin shells.
    
    \tparam T coefficient type
    
    \ingroup ThinShell
*/
template <class T>
class gsVisitorLinShell
{
public:

    /// Constructor with thickness, material parameters and the surface force as inputs.
    gsVisitorLinShell(T thickness, T lambda, T mu, const gsFunction<T> & srf_force) : 
    m_thickness(thickness),
    m_lambda(lambda),
    m_mu(mu),
    m_surfaceForce_ptr(&srf_force)
    { }

    /// Function to initialize the assembly procedure.
    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T>    & rule, 
                    unsigned         & evFlags )
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_JACOBIAN | NEED_MEASURE |  NEED_2ND_DER;
    }

    /// Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T> const      & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(quNodes.col(0), actives);
        numActive = actives.rows();
        
        // Evaluate basis functions on element
        basis.evalAllDers_into( quNodes, 2, basisData);
        
        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval.evaluateAt(quNodes);
        
        // Evaluate right-hand side at the geometry points
        m_surfaceForce_ptr->eval_into( geoEval.values(), forceVals );
        
        // Initialize local matrix/rhs
        localMat.setZero(3*numActive, 3*numActive);
        localRhs.setZero(3*numActive, 1          );

        // Initialize auxiliary matrices
        E_m_der.resize(3,3*numActive);
        E_f_der.resize(3,3*numActive);
    }
    
    /// Assembles the local matrix and right hand side using only the undeformed geometry.
    inline void assemble(gsDomainIterator<T>    & /*element*/,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {           
            // Compute needed quantities...
            computeMaterialMatrix(geoEval, k);
            computeStrainDers(geoEval, bGrads, k);

            // Multiply weight by the geometry measure and the thickness
            const T weight = quWeights[k] * geoEval.measure(k) * m_thickness;

            localMat.noalias() += weight *  ( E_m_der.transpose() * m_C * E_m_der + 
                                              E_f_der.transpose() * m_C * E_f_der * 
                                              (m_thickness*m_thickness/3.0) );
            
            for (index_t j = 0; j!= 3; ++j)
                localRhs.middleRows(j*numActive,numActive).noalias() += 
                    weight * forceVals(j,k) * bVals.col(k) / m_thickness ;
        }
        //gsDebug<< "local Mat: \n"<< localMat << "\n";
    }
    
    /// Putting local information in the global matrices and rhs.
    inline void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                              const gsMatrix<T>     & eliminatedDofs,
                              const index_t           /*patchIndex*/,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
        for (index_t ci = 0; ci!= 3; ++ci)
            for (index_t ai=0; ai < numActive; ++ai)
            {
                const index_t gi = ci * numActive +  ai; // row index
                const index_t ii = mappers[ci].index( actives(ai) );

                if ( mappers[ci].is_free_index(ii) )
                {
                    rhsMatrix.row(ii) += localRhs.row(gi);
                    
                    for (index_t cj = 0; cj!= 3; ++cj)
                        for (index_t aj=0; aj < numActive; ++aj)
                        {
                            const index_t gj = cj * numActive +  aj; // column index
                            const index_t jj = mappers[cj].index( actives(aj) ); 
                            
                            if ( mappers[cj].is_free_index(jj) )
                            {
                                sysMatrix.coeffRef(ii, jj) += localMat(gi, gj);
                            }
                            else // Fixed DoF ?
                            {
                                rhsMatrix.row(ii).noalias() -= localMat(gi, gj) * 
                                    eliminatedDofs.row( mappers[cj].global_to_bindex(jj) );
                            }
                        }
                }
            }
    }
    

    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

    /// Computes the material matrix \em m_C at quadrature point \em k
    void computeMaterialMatrix(const gsGeometryEvaluator<T> & geoEval,
                               const index_t k)
    {
        // ---------------  Material matrix
        gsMatrix<T,3,3> F0;
        geoEval.normal(k,normal);
        normal.normalize();
        F0.leftCols(2) = geoEval.jacobian(k);
        F0.col(2)      = normal;
        
        //F0 = F0.inverse(); F0 = F0 * F0.transpose();
        F0 = F0.inverse() * F0.inverse().transpose();

        const T C_constant = 4*m_lambda*m_mu/(m_lambda+2*m_mu);

        m_C(0,0) = C_constant*F0(0,0)*F0(0,0) + 2*m_mu*(2*F0(0,0)*F0(0,0));
        m_C(1,1) = C_constant*F0(1,1)*F0(1,1) + 2*m_mu*(2*F0(1,1)*F0(1,1));
        m_C(2,2) = C_constant*F0(0,1)*F0(0,1) + 2*m_mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
        m_C(1,0) = 
        m_C(0,1) = C_constant*F0(0,0)*F0(1,1) + 2*m_mu*(2*F0(0,1)*F0(0,1));
        m_C(2,0) = 
        m_C(0,2) = C_constant*F0(0,0)*F0(0,1) + 2*m_mu*(2*F0(0,0)*F0(0,1));
        m_C(2,1) = m_C(1,2) = C_constant*F0(0,1)*F0(1,1) + 2*m_mu*(2*F0(0,1)*F0(1,1)); 
        //gsDebug<< "C: \n"<< m_C << "\n";
    }

    /// Computes the membrane and flexural strain first derivatives at quadrature point \em k
    void computeStrainDers(const gsGeometryEvaluator<T> & geoEval,
                           const gsMatrix<T> & bGrads,
                           const index_t k)
    {
        gsMatrix<T> & bGrads2 = basisData[2];

        gsAsConstMatrix<T,3,3> GsecDer( geoEval.deriv2(k).data(),3,3 );

        gsVector<T,3> m_v, n_der;

        geoEval.normal(k,normal); //geoEval or defShell
        normal.normalize();

        const typename gsMatrix<T>::constColumns & Jac = geoEval.jacobian(k);
        for (index_t j = 0; j!= 3; ++j)
        {
            const index_t s = j*numActive;
            for (index_t i = 0; i!= numActive; ++i)
            {
                // ---------------  Membrane strain derivative
                E_m_der(0,i+s) = bGrads(2*i  ,k) * Jac(j,0) ;
                E_m_der(1,i+s) = bGrads(2*i+1,k) * Jac(j,1) ;
                E_m_der(2,i+s) = bGrads(2*i  ,k) * Jac(j,1) + 
                                 bGrads(2*i+1,k) * Jac(j,0) ;

                m_v.noalias() = vecFun(j,bGrads(2*i,k)).cross( 
                                geoEval.jacobian(k).template block<3,1>(0,1) )
                                - vecFun(j,bGrads(2*i+1,k)).cross( 
                                  geoEval.jacobian(k).template block<3,1>(0,0) );               

                // ---------------  First variation of the normal
                n_der.noalias() = (m_v - ( normal.dot(m_v) ) * normal) / geoEval.measure(k);

                // ---------------  Bending strain derivative
                E_f_der.col(i+s) = bGrads2.template block<3,1>(3*i,k) * normal[j] 
                                   + GsecDer * n_der  ;
                E_f_der(2, i+s) *= 2.0 ;
            }
        }
        //gsDebug<< "E_m_der: \n"<< E_m_der << "\n";
        //gsDebug<< "E_f_der: \n"<< E_f_der << "\n";
    }

    /// Extends a scalar to a 3D vector by putting \em val at index \em pos of a zero vector
    static inline gsVector<T,3> vecFun(index_t pos, T val) 
    { 
        gsVector<T,3> result = gsVector<T,3>::Zero();
        result[pos] = val;
        return result;
    }

protected:

    /// Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<index_t> actives;
    index_t numActive;

    /// Normal to the shell centerline
    gsVector<T> normal;

    /// Derivatives of the membrane strain
    gsMatrix<T,3> E_m_der;

    /// Derivatives of the bending (or flexural) strain
    gsMatrix<T,3> E_f_der;

    /// Material matrix
    gsMatrix<T,3,3> m_C;
    

protected:

    /// Half the shell thickness
    T m_thickness;

    /// Lame's constants
    T m_lambda, m_mu;

protected:

    /// Pointer to the surface forces
    const gsFunction<T> * m_surfaceForce_ptr;

    /// Local values of the surface forces
    gsMatrix<T> forceVals;
    
protected:
  
    /// Local matrix
    gsMatrix<T> localMat;
    
    /// Local right hand side
    gsMatrix<T> localRhs;
};


} // namespace gismo

