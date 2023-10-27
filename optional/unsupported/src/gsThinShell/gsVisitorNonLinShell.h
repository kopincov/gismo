/** @file gsVisitorNonLinShell.h

    @brief Element visitor for nonlinear elasticity on thin shell.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Goyal, A. Mantzaflaris
*/

#pragma once

#include <gsThinShell/gsVisitorLinShell.h>

namespace gismo
{

/** 
    @brief The visitor computes element system matrix and right hand side for nonlinear 
            elasticity on thin shells.
    
    \tparam T coefficient type
    
    \ingroup ThinShell
*/
template <class T>
class gsVisitorNonLinShell : public gsVisitorLinShell<T>
{
public:
    typedef gsVisitorLinShell<T> Base;
public:

    /// Constructor with thickness, material parameters and the surface force as inputs.
    gsVisitorNonLinShell(T thickness, T lambda, T mu, 
                         const gsFunction<T> & srf_force, 
                         const gsGeometry<T> & deformed) 
    : Base(thickness, lambda, mu, srf_force),
      defShell(getEvaluator(NEED_MEASURE|NEED_2ND_DER|NEED_JACOBIAN, deformed))
    { 

    }
    
    /// Sets the gsGeometryEvaluator \em defShell using \em deformed
    void setDeformed(const gsGeometry<T> & deformed)
    {
        defShell = memory::make_unique(getEvaluator(NEED_MEASURE|NEED_2ND_DER|NEED_JACOBIAN, deformed));
    }

    /// Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis,
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T> const      & quNodes)
    {
        Base::evaluate(basis,geoEval,quNodes);

        // Evaluate deformed shell 
        defShell->evaluateAt(quNodes);

        //E_m_der2_comp.resize(3,   numActive*numActive);
        E_m_der2_comp.resize(3, 9*numActive*numActive);
        E_f_der2     .resize(3, 9*numActive*numActive);
    }

    /// Assembles the local matrix and right hand side using both deformed and undeformed geometries.
    inline void assemble(gsDomainIterator<T>    & /*element*/,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];

        //gsDebug<< quWeights.size() <<"\n";

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {           
            // Compute needed things...
            Base::computeMaterialMatrix( geoEval , k);
            computeStrainsAndSecDer    ( geoEval , bGrads, k);
            Base::computeStrainDers    (*defShell, bGrads, k);

            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k) * m_thickness;

            localMat.noalias() += weight *  ( E_m_der.transpose() * m_C * E_m_der + 
                                              E_f_der.transpose() * m_C * E_f_der * 
                                              (m_thickness*m_thickness/3.0) );

            localMat.reshape(1,9*numActive*numActive).noalias() -= 
                (weight*m_thickness*m_thickness/3.0) * (E_f.transpose() * m_C * E_f_der2) ;

            localMat.reshape(1,9*numActive*numActive).noalias() += 
                weight * ( E_m.transpose() * m_C * E_m_der2_comp );

            localRhs.transpose().noalias() -= weight * (
                E_m.transpose() * m_C * E_m_der -
               (E_f.transpose() * m_C * E_f_der )* (m_thickness*m_thickness/3.0)
                );
            
            for (index_t j = 0; j!= 3; ++j)
            {
                localRhs.middleRows(j*numActive,numActive).noalias() += 
                     weight * forceVals(j,k) * bVals.col(k) / m_thickness;

                //localMat.block(j*numActive,j*numActive,numActive,numActive) += weight *
                //    gsMatrix<T>( E_m.transpose() * m_C * E_m_der2_comp ).reshape(numActive,numActive);
            }
        }

        //gsDebug<< "local Mat: \n"<< localMat << "\n";
    }


    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

    /// Computes the membrane and flexural strains and their second derivatives at quadrature point \em k
    void computeStrainsAndSecDer(const gsGeometryEvaluator<T> & geoEval,
                                 const gsMatrix<T> & bGrads,
                                 const index_t k)
    {
        gsMatrix<T> & bGrads2 = basisData[2];

        // 2nd derivatives of undeformed shell
        gsAsConstMatrix<T,3,3> GsecDer  (geoEval.deriv2(k).data(),3,3   );
        // 2nd derivatives of deformed shell
        gsAsConstMatrix<T,3,3> defSecDer(defShell->deriv2(k).data(),3,3 );

        geoEval.normal(k, normal);// already computed 
        normal.normalize();

        defShell->normal(k, defNormal);
        defNormal.normalize();

        const typename gsMatrix<T>::constColumns & oriJac = geoEval.jacobian(k);
        const typename gsMatrix<T>::constColumns & defJac = defShell->jacobian(k);

        // ---------------  E_m
        E_m[0] = (defJac.col(0).squaredNorm() - oriJac.col(0).squaredNorm())/2.0;
        E_m[1] = (defJac.col(1).squaredNorm() - oriJac.col(1).squaredNorm())/2.0;
        E_m[2] = (defJac.col(0).transpose() * defJac.col(1) ).value() - 
                 (oriJac.col(0).transpose() * oriJac.col(1) ).value();

        // ---------------  E_f
        E_f = GsecDer * normal - defSecDer * defNormal;
        E_f[2] *= 2.0;

        // --------------- First variation of the normal 
        gsMatrix<T> m_v(3,3*numActive), n_der(3,3*numActive);
        index_t c = 0;
        for (index_t j = 0; j!= 3; ++j)
            for (index_t i = 0; i!= numActive; ++i)
            {
                // m_v
                m_v.col(c).noalias() = ( vecFun(j,bGrads(2*i,k)).cross( 
                                         defJac.template block<3,1>(0,1) ) - 
                                         vecFun(j,bGrads(2*i+1,k)).cross( 
                                         defJac.template block<3,1>(0,0) )
                                        ) / defShell->measure(k);

                // n_der
                n_der.col(c).noalias() = ( m_v.col(c) - normal.dot(m_v.col(c)) * normal );

                c++;
            }

        //gsDebug<<"m_v\n:"<< m_v <<"\n"; 
        //gsDebug<<"n_der:\n"<< n_der <<"\n";

/* check
        E_m_der2_comp.reshape(3*numActive,3*numActive)
        E_m_der2_comp.row(0) = gsMatrix<T>( bGrads.row(0).transpose() * bGrads.row(0) )
                               .reshape(1,numActive*numActive);

        E_m_der2_comp.row(1) = gsMatrix<T>( bGrads.row(1).transpose() * bGrads.row(1) )
                               .reshape(1,numActive*numActive);
        
        const gsMatrix<T> temp = bGrads.row(0).transpose() *  bGrads.row(1);
        E_m_der2_comp.row(2) = gsMatrix<T>( temp + temp.transpose() )
                               .reshape(1,numActive*numActive);
*/        
        gsVector<T,3> m_v_der, n_der2;
        c = 0;
        for (index_t ci = 0; ci!= 3; ++ci)
        {
            for (index_t ai=0; ai < numActive; ++ai)
            {
                const index_t gi = ci*numActive + ai; // v 

                for (index_t cj = 0; cj!= 3; ++cj)
                {
                    for (index_t aj=0; aj < numActive; ++aj)
                    {
                        const index_t gj = cj*numActive + aj; // w

                        // ---------------  E_m_der2
                        E_m_der2_comp(0,c) =  vecFun(ci, bGrads(2*ai  ,k) ).dot( 
                                              vecFun(cj, bGrads(2*aj  ,k) ) );
                        E_m_der2_comp(1,c) =  vecFun(ci, bGrads(2*ai+1,k) ).dot( 
                                              vecFun(cj, bGrads(2*aj+1,k) ) );
                        E_m_der2_comp(2,c) =  vecFun(ci, bGrads(2*ai ,k) ).dot( 
                                              vecFun(cj, bGrads(2*aj+1,k) ) )
                                              +
                                              vecFun(ci, bGrads(2*ai+1,k) ).dot( 
                                              vecFun(cj, bGrads(2*aj,k) ) )
                                              ;

                        // m_v_der
                        m_v_der = ( vecFun(ci, bGrads(2*ai  ,k) ).cross( 
                                    vecFun(cj, bGrads(2*aj+1,k))  ) + 
                                    vecFun(cj, bGrads(2*aj  ,k) ).cross( 
                                    vecFun(ci, bGrads(2*ai+1,k))  ) 
                                  ) / defShell->measure(k)
                                  - normal.dot( m_v.col(gj) ) *  m_v.col(gi);

                        // n_der2
                        n_der2.noalias() = m_v_der - ( m_v.col(gi).dot(n_der.col(gj)) + 
                                             normal.dot(m_v_der) )   * normal
                                         -   normal.dot(m_v.col(gi)) * n_der.col(gj);

                        // ---------------  E_f_der2                        
                        E_f_der2.col(c) = bGrads2.template block<3,1>(3*ai,k) * n_der(ci,gj) + 
                                          bGrads2.template block<3,1>(3*aj,k) * n_der(cj,gi)
                                          + defSecDer * n_der2;

                        E_f_der2(2,c) *= 2.0 ;

                        c++;
                    }
                }
            }
        }

/*
gsDebug<< "0E_m_der2: \n"<< gsAsMatrix<T>(E_m_der2_comp.row(0).data(),3*numActive,3*numActive) << "\n";
gsDebug<< "1E_m_der2: \n"<< gsAsMatrix<T>(E_m_der2_comp.row(1).data(),3*numActive,3*numActive) << "\n";
gsDebug<< "2E_m_der2: \n"<< gsAsMatrix<T>(E_m_der2_comp.row(2).data(),3*numActive,3*numActive) << "\n";
*/       
        //gsDebug<< "E_f_der2: \n"<< E_f_der2 << "\n";

    } // end

        
protected:

    /// Contains the geometry evaluations for the deformed configuration
    typename gsGeometryEvaluator<T>::uPtr defShell;

    /// Deformed normal at the shell centerline
    gsVector<T> defNormal;

protected:

    using Base::vecFun;

    // Basis values
    using Base::basisData;
    using Base::actives;
    using Base::numActive;
    using Base::normal;

    /// Membrane strain
    gsVector<T,3> E_m;

    /// Bending (or flexural) strain
    gsVector<T,3> E_f;

    // Derivative of the membrane strain
    using Base::E_m_der;

    // Derivative of the bending (or flexural) strain
    using Base::E_f_der;

    // Material matrix
    using Base::m_C;

    /// Second derivative of membrane strain
    gsMatrix<T,3> E_m_der2_comp;

    /// Second derivative of bending (or flexural) strain
    gsMatrix<T,3> E_f_der2;
    

protected:

    // thickness
    using Base::m_thickness;

    // Lambda and mu
    using Base::m_lambda;
    using Base::m_mu;

protected:

    // Right hand side
    using Base::m_surfaceForce_ptr;

    
protected:
    // Local values of the surface forces
    using Base::forceVals;

protected:
    // Local matrices
    using Base::localMat;
    using Base::localRhs;
};


} // namespace gismo

