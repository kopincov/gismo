/** @file gsVisitorDgSpaceTime.h

    @brief A DG interface visitor for the Parabolic problem .

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Moore
*/

#pragma once

#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{
/** @brief
    Implementation of a interface condition for the 
    discontinuous Galerkin Space-Time Assembler.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can also be imposed weakly (i.e Nitsche ) 
*/

template <class T>
class gsVisitorDgSpaceTime
{
public:
/** \brief Visitor for adding the interface conditions for the 
 * interior penalty method of the parabolic problem.
 * 
 * This visitor adds the following term to the left-hand side (bilinear form).
 * \f[ -\{\nabla u \cdot \mathbf{n}\}[ \partial_t v\cdot \mathbf{n} ] -
 * \{\nabla v \cdot \mathbf{n}\}[\ u \cdot \mathbf{n}] 
 * + \mu [ \partial_t u ][ \partial_t v]
 *  \f].
 * Where \f[v\f] is the test function and \f[ u \f] is trial function.
 */

    gsVisitorDgSpaceTime(T _penalty, boxSide s) :
    penalty(_penalty), side1(s), d(0)
    { }

    void initialize(const gsBasis<T> & basis1, 
                    const gsBasis<T> & basis2, 
                    gsQuadRule<T> & rule,
                    unsigned & evFlags
        )
    {
        d = basis1.dim();
        const int dir = side1.direction();
        gsVector<int> numQuadNodes ( d );
        for (int i = 0; i < basis1.dim(); ++i)
            numQuadNodes[i] = 2*basis1.degree(i) + 1;
        numQuadNodes[dir] = 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T> const       & B1, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval1,
                         gsBasis<T> const       & B2, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval2,
                         gsMatrix<T>            & quNodes1,
                         gsMatrix<T>            & quNodes2)
    {
        // Compute the active basis functions
        B1.active_into(quNodes1.col(0), actives1);
        B2.active_into(quNodes2.col(0), actives2);
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        // Evaluate basis functions and their first derivatives
        B1.evalAllDers_into( quNodes1, 1, basisData1);
        B2.evalAllDers_into( quNodes2, 1, basisData2);
        
        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval1.evaluateAt(quNodes1);
        geoEval2.evaluateAt(quNodes2);
        
        // Initialize local matrices
        B11.setZero(numActive1, numActive1); B12.setZero(numActive1, numActive2);
        E11.setZero(numActive1, numActive1); E12.setZero(numActive1, numActive2);
        B22.setZero(numActive2, numActive2); B21.setZero(numActive2, numActive1);
        E22.setZero(numActive2, numActive2); E21.setZero(numActive2, numActive1);
    }

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
                         gsDomainIterator<T>    & element2,
                         gsGeometryEvaluator<T> & geoEval1,
                         gsGeometryEvaluator<T> & geoEval2,
                         gsVector<T>            & quWeights)
    {
        //const unsigned d = element1.dim();
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();
        
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            // Compute the outer normal vector from patch1
            // unormal : dim x 1
            geoEval1.outerNormal(k, side1, unormal);

            // Compute the outer normal for space and time
            unormalSpace = unormal.head(d-1);
            unormalTime  = unormal.tail(1);
            
            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T Sweight = quWeights[k] * unormalSpace.norm();
            //const T Tweight = quWeights[k] * unormalTime.norm();
        
            // Compute the outer unit normal vector from patch1 in place
            unormalSpace.normalize();
            // unormalTime.normalize();
                    
            // Take blocks of values and derivatives of basis functions
            // const typename gsMatrix<T>::Column val1 = basisData1.col(k);//bug
            const typename gsMatrix<T>::Block val1 = basisData1[0].block(0,k,numActive1,1);
            gsMatrix<T> & grads1 = basisData1[1];// all grads
            //const typename gsMatrix<T>::Column val2 = basisData2.col(k);//bug
            const typename gsMatrix<T>::Block val2 = basisData2[0].block(0,k,numActive2,1);
            gsMatrix<T> & grads2 = basisData2[1];// all grads
        
            // Transform the basis gradients
            geoEval1.transformGradients(k, grads1, phGrad1);
            geoEval2.transformGradients(k, grads2, phGrad2);
            
            // Gradient in Space : (Dim - 1) X NumActive matrix
            physGradSpace1   = phGrad1.topRows(d-1);
            physGradSpace2   = phGrad2.topRows(d-1);
            
            // Gradient in Time  : 1 X NumActive matrix
            physGradTime1    = phGrad1.bottomRows(1);
            physGradTime2    = phGrad2.bottomRows(1);                    
        
            // Compute element matrices
            const T c1     = Sweight * T(0.5);
            N1.noalias()   = unormalSpace.transpose() * physGradSpace1;
            N2.noalias()   = unormalSpace.transpose() * physGradSpace2;
            
            S11.noalias() +=  physGradSpace1.transpose() * physGradSpace1;
            S12.noalias() +=  physGradSpace1.transpose() * physGradSpace2;
            S21.noalias() +=  physGradSpace2.transpose() * physGradSpace1;
            S22.noalias() +=  physGradSpace2.transpose() * physGradSpace2;
             
            B11.noalias() += c1 * ( physGradTime1.transpose() * N1 );
            B12.noalias() += c1 * ( physGradTime1.transpose() * N2 );
            B22.noalias() -= c1 * ( physGradTime2.transpose() * N2 );
            B21.noalias() -= c1 * ( physGradTime2.transpose() * N1 ); 
            
            //
            const T h1     = element1.getCellSize();
            const T h2     = element2.getCellSize();
            const T c2     = Sweight * penalty*(1/h1 + 1/h2)  ;
            E11.noalias() += c2 * ( physGradTime1.transpose() * physGradTime1 );
            E12.noalias() += c2 * ( physGradTime1.transpose() * physGradTime2 );
            E22.noalias() += c2 * ( physGradTime2.transpose() * physGradTime2 );
            E21.noalias() += c2 * ( physGradTime2.transpose() * physGradTime1 );
            
        }
    }
    
    void localToGlobal(const gsDofMapper     & mapper,
                       const gsMatrix<T>     & eliminatedDofs,
                       const index_t patch1,
                       const index_t patch2,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        // Local DoFs to global DoFs
        mapper.localToGlobal(actives1, patch1, actives1);
        mapper.localToGlobal(actives2, patch2, actives2);
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();
        
//         if (unormalTime >= 0)
//         {
//             // Push element contributions 1-2 to the global matrix and load vector
//             for (index_t j=0; j!=numActive1; ++j)
//             {
//                 const index_t jj1 = actives1(j); // N1_j
// 
//                 for (index_t i=0; i!=numActive1; ++i)
//                 {
//                     const index_t  ii1 = actives1(i); // N1_i
//                     
//                     if ( jj1 <= ii1 )
//                         sysMatrix( ii1, jj1 ) -=  B11(i,j) - B11(j,i) - E11(i,j) - S11(i,j);
//                 }
//                 
//                 for (index_t i=0; i!=numActive2; ++i)
//                 {
//                     const index_t  ii2 = actives2(i); // N2_i
//                     
//                     if ( jj1 <= ii2 ) 
//                         sysMatrix( ii2, jj1)  -=  B21(i,j) - B12(j,i) + E21(i,j) + S21(i,j);
//                 }
//             }
//             // Push element contributions 2-1 to the global matrix and load vector
//             for (index_t j=0; j!=numActive2; ++j)
//             {
//                 const index_t jj2 = actives2(j); // N2_j
//                 
//                 for (index_t i=0; i!=numActive2; ++i)
//                 {
//                     const index_t  ii2 = actives2(i); // N2_i
//                     
//                     if ( jj2 <= ii2 ) 
//                         sysMatrix( ii2, jj2 ) -=  B22(i,j) - B22(j,i) - E22(i,j) ;
//                 }
//                 
//                 for (index_t i=0; i!=numActive1; ++i)
//                 {
//                     const index_t  ii1 = actives1(i); // N1_i
//                     
//                     if ( jj2 <= ii1 ) 
//                         sysMatrix( ii1, jj2)  -=  B12(i,j) - B21(j,i) + E12(i,j) ;
//                 }
//             }                  
//         }
//         else 
//         {
            // Push element contributions 1-2 to the global matrix and load vector
            for (index_t j=0; j!=numActive1; ++j)
            {
                const index_t jj1 = actives1(j); // N1_j

                for (index_t i=0; i!=numActive1; ++i)
                {
                    const index_t  ii1 = actives1(i); // N1_i
                    
                    if ( jj1 <= ii1 )
                        sysMatrix( ii1, jj1 ) -=  B11(i,j) - B11(j,i) - E11(i,j);
                }
                
                for (index_t i=0; i!=numActive2; ++i)
                {
                    const index_t  ii2 = actives2(i); // N2_i
                    
                    if ( jj1 <= ii2 ) 
                        sysMatrix( ii2, jj1)  -=  B21(i,j) - B12(j,i) + E21(i,j);
                }
            }
        // Push element contributions 2-1 to the global matrix and load vector
            for (index_t j=0; j!=numActive2; ++j)
            {
                const index_t jj2 = actives2(j); // N2_j
                
                for (index_t i=0; i!=numActive2; ++i)
                {
                    const index_t  ii2 = actives2(i); // N2_i
                    
                    if ( jj2 <= ii2 ) 
                        sysMatrix( ii2, jj2 ) -=  B22(i,j) - B22(j,i) - E22(i,j) - S22(i,j);
                }
                
                for (index_t i=0; i!=numActive1; ++i)
                {
                    const index_t  ii1 = actives1(i); // N1_i
                    
                    if ( jj2 <= ii1 ) 
                        sysMatrix( ii1, jj2)  -=  B12(i,j) - B21(j,i) + E12(i,j) - S12(i,j);
                }
            }   
        } 
//     }

private:

    // Penalty constant
    T penalty;

    // Side
    boxSide side1;

    // dimension of the problem
    unsigned d;

private:

    // Basis values etc
    std::vector<gsMatrix<T> > basisData1, basisData2;
    gsMatrix<T>        phGrad1   , phGrad2;
    gsMatrix<T>    physGradSpace1, physGradSpace2;
    gsMatrix<T>    physGradTime1, physGradTime2;
    gsMatrix<index_t> actives1  , actives2;



    // Outer normal
    gsVector<T> unormal, unormalSpace, unormalTime;

    // Auxiliary element matrices
    gsMatrix<T> B11, B12, E11, E12, N1,
                B22, B21, E22, E21, N2,
                S11, S12, S21, S22;
};
} // namespace gismo
