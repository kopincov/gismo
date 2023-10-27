/** @file gsVisitorGradGradBoundary.h

    @brief Assembles the Grad Grad for the Boundary.

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
    Implementation of a Grad-Grad BC for elliptic Assembler.

    It sets up an assembler and adds the following term to the bilinear term.
   \f[ (\nabla u \cdot \nabla u)_{\partial \Omega} \f]
*/

template <class T>
class gsVisitorGradGradBoundary
{
public:
/** @brief
 * Constructor of the assembler object 
 * 
   \param[in] s are the sides of the geometry where grad-grad BC are prescribed.
   
   \f[ (\nabla u \cdot \nabla u)_{\partial \Omega} \f]
*/
    gsVisitorGradGradBoundary(const boxSide s) : 
    side(s)
    { }

    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T>    & rule, 
                    unsigned         & evFlags )
    {
        const int dir = side.direction();
        gsVector<int> numQuadNodes ( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
            numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_JACOBIAN|NEED_GRAD_TRANSFORM;;
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         // todo: add element here for efficiency
                         gsMatrix<T>            & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(quNodes.col(0) , actives);
        const index_t numActive = actives.rows();
 
        // Evaluate basis functions on element
        basis.evalAllDers_into(quNodes, 1, basisData);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate the Neumann data

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
    }

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        const unsigned d = element.dim();
        const index_t numActive = actives.rows();
        gsMatrix<T> basisPhGrads ;
        
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);
            
            // Compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisData, basisPhGrads);
            
            localMat.noalias() += weight * ( basisPhGrads.transpose() * basisPhGrads );
        }
    }
    
    void localToGlobal(const gsDofMapper     & mapper,
                       const gsMatrix<T>     & eliminatedDofs,
                       const index_t           patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();

        // Push element contribution to the global matrix
        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = actives(j);
            //rhsMatrix.row(jj) -= localRhs.row(j);
            for (index_t i=0; i!=numActive; ++i)
            {
                const unsigned ii = actives(i);
                if ( jj <= ii ) // assuming symmetric problem (!) probably we must not.
                    sysMatrix( ii, jj ) -= localMat(i,j);
            }
        }
    }

protected:
    
    // Neumann function
    boxSide side;

    // Basis values
    gsMatrix<T>      basisData;
    gsMatrix<T>      basisPhGrads;
    gsMatrix<index_t> actives;
    index_t numActive;    

    // Normal 
    gsVector<T> unormal;

    // Local matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

};


} // namespace gismo
