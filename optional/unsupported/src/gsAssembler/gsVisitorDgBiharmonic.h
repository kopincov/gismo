/** @file gsVisitorDgBiharmonic.h

    @brief A first-type dG interface visitor for 
    the Biharmonic problem .

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Moore
*/

#pragma once

namespace gismo
{
/** @brief
    Implementation of the interface condition for the 
    discontinuous Galerkin Biharmonic Assembler of 
    the first type.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can also be weakly enforced (i.e Nitsche )
*/
template <class T>
class gsVisitorDgBiharmonic
{
public:
/** \brief Visitor for adding the interface conditions for the 
 * interior penalty method of the biharmonic problem.
 * 
 * This visitor adds the following term to the bilinear form.
 * \f[ -\{\nabla^2 u \cdot \mathbf{n}\}[\nabla v\cdot \mathbf{n} ] -
 * \{\nabla^2 v \cdot \mathbf{n}\}[\nabla u \cdot \mathbf{n}] 
 * + \alpha [\nabla u \cdot \mathbf{n}][\nabla v \cdot \mathbf{n}] \f].
 * Where \f[ v$ \f]  is the test function and \f[u \f] is trial function.
 */

    gsVisitorDgBiharmonic(T _penalty, boxSide s) :
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
            numQuadNodes[i] = basis1.degree(i) + 1;
            numQuadNodes[dir] = 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
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

       
        // Evaluate basis functions and their derivatives
        B1.evalAllDers_into( quNodes1, 2, basisData1);
        B2.evalAllDers_into( quNodes2, 2, basisData2);        
        
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
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsDomainIterator<T>    & element2,
                         gsGeometryEvaluator<T> & geoEval1,
                         gsGeometryEvaluator<T> & geoEval2,
                         gsVector<T>            & quWeights)
    {
        const unsigned d = element.dim();
        
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();
        
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            // Compute the outer normal vector from patch1
            geoEval1.outerNormal(k, side1, unormal);
        
            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T weight = quWeights[k] * unormal.norm();
        
            // Compute the outer unit normal vector from patch1 in place
            unormal.normalize();
        
            // Take blocks of values and derivatives of basis functions
            //const typename gsMatrix<T>::Block val1 = basisData1.block(0,k,numActive1,1);
            const typename gsMatrix<T>::Block basisGrad1 = // all grads
                basisData1.middleRows(numActive1, numActive1 * d );
            const typename gsMatrix<T>::Block basis2ndDerivs1 = // all second derivs
                basisData1.bottomRows( numActive1 * ((d+1) *d )/2 );     
                
            //const typename gsMatrix<T>::Block val2 = basisData2.block(0,k,numActive2,1);
            const typename gsMatrix<T>::Block basisGrad2 = // all grads
                basisData2.middleRows(numActive2, numActive2 * d );
            const typename gsMatrix<T>::Block basis2ndDerivs2 = // all second derivs
                basisData2.bottomRows( numActive2 * ((d+1) *d )/2 );                
        
           // Compute physical gradients at k as a Dim x NumActive matrix
            geoEval1.transformGradients(k, basisGrad1, physBasisGrad1);
            geoEval2.transformGradients(k, basisGrad2, physBasisGrad2);
            
           // Compute physical laplacian at k as a 1 x numActive matrix
            geoEval1.transformLaplaceHgrad (k, basisGrad1, basis2ndDerivs1, physBasisLaplace1);           
            geoEval2.transformLaplaceHgrad (k, basisGrad2, basis2ndDerivs2, physBasisLaplace2);           
            
            // Compute element matrices
            const T c1     = weight * T(0.5);
            //N1.noalias()   = unormal.transpose() * physBasisLaplace1;
            //N2.noalias()   = unormal.transpose() * physBasisLaplace2;
            N1.noalias()   = physBasisGrad1.transpose()*unormal;
            N2.noalias()   = physBasisGrad2.transpose()*unormal;
            B11.noalias() += c1 * ( N1 * physBasisLaplace1 );
            B12.noalias() += c1 * ( N1 * physBasisLaplace2 );
            B22.noalias() -= c1 * ( N2 * physBasisLaplace2 );
            B21.noalias() -= c1 * ( N2 * physBasisLaplace1 );
        
            const T c2     = weight * penalty / element.getCellSize(); //TODO: fix this like in VisitorDG
            E11.noalias() += c2 * ( physBasisGrad1.transpose() * physBasisGrad1 );
            E12.noalias() += c2 * ( physBasisGrad1.transpose() * physBasisGrad2 );
            E22.noalias() += c2 * ( physBasisGrad2.transpose() * physBasisGrad2 );
            E21.noalias() += c2 * ( physBasisGrad2.transpose() * physBasisGrad1 );
        }
    }
    
    void localToGlobal(const gsDofMapper  & mapper,
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
        
        // Push element contributions 1-2 to the global matrix and load vector
        for (index_t j=0; j!=numActive1; ++j)
        {
            const index_t jj1 = actives1(j); // N1_j

            for (index_t i=0; i!=numActive1; ++i)
            {
                const index_t  ii1 = actives1(i); // N1_i
                
                if ( jj1 <= ii1 )
                    sysMatrix( ii1, jj1 ) -=  B11(i,j) + B11(j,i) - E11(i,j);
            }
            
            for (index_t i=0; i!=numActive2; ++i)
            {
                const index_t  ii2 = actives2(i); // N2_i
                
                if ( jj1 <= ii2 ) 
                    sysMatrix( ii2, jj1)  -=  B21(i,j) + B12(j,i) + E21(i,j);
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
                    sysMatrix( ii2, jj2 ) -=  B22(i,j) + B22(j,i) - E22(i,j);
            }
            
            for (index_t i=0; i!=numActive1; ++i)
            {
                const index_t  ii1 = actives1(i); // N1_i
                
                if ( jj2 <= ii1 ) 
                    sysMatrix( ii1, jj2)  -=  B12(i,j) + B21(j,i) + E12(i,j);
            }
        }    
    }

private:

    // Penalty constant
    T penalty;

    // Side
    boxSide side1;

    // dimension of the problem
    unsigned d;

private:

    // Basis values etc
    gsMatrix<T>        basisData1, basisData2;
    gsMatrix<T>        physBasisGrad1    , physBasisGrad2;
    gsMatrix<T>        physBasisLaplace1 , physBasisLaplace2;
    gsMatrix<index_t> actives1  , actives2;

    // Outer normal
    gsVector<T> unormal;

    // Auxiliary element matrices
    gsMatrix<T> B11, B12, E11, E12, N1,
        B22, B21, E22, E21, N2;
// protected:
//     void basisData2(numActive2*(d+2)*(d+1)*d/6 arg1);
};


} // namespace gismo
