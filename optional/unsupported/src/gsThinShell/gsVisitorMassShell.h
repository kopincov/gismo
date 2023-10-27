/** @file gsVisitorMassShell.h

    @brief Mass visitor for assembling element mass matrix for a shell

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
    @brief The visitor computes element mass integrals for thin shells
    
    \ingroup ThinShell
*/
template <class T>
class gsVisitorMassShell
{
public:

    /// Constructor with thickness and material density as inputs.
    gsVisitorMassShell(T rho, T thickness) : 
    m_rho(rho),
    m_thickness(thickness)
    { }

    /// Function to initialize the assembly procedure.
    static void initialize(const gsBasis<T> & basis, 
                           gsQuadRule<T> & rule, 
                           unsigned & evFlags )
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE;
    }

    /// Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T>            & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(quNodes.col(0) , actives);
        numActive = actives.rows();
 
        // Evaluate basis functions on element
        basis.eval_into(quNodes, basisData);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Initialize local matrix
        localMat.setZero(numActive,numActive);
    }

    /// Assembles the local mass matrix.
    inline void assemble(gsDomainIterator<T>    & /*element*/,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply quadrature weight by the geometry measure and material data
            const T weight = 2 * m_thickness * m_rho * quWeights[k] * geoEval.measure(k);
        
            localMat += weight * ( basisData.col(k) * basisData.col(k).transpose() );
        }
    }
    
    /// Putting local information in the global mass matrix.
    void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                       const gsMatrix<T>     & /*eliminatedDofs*/,
                       const index_t           /*patchIndex*/,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & /*rhsMatrix*/ )
    {
        for (index_t ci = 0; ci!= 3; ++ci)
            for (index_t ai=0; ai < numActive; ++ai)
            {
                const index_t ii = mappers[ci].index( actives(ai) );

                if ( mappers[ci].is_free_index(ii) )
                {
                    for (index_t cj = 0; cj!= 3; ++cj)
                        for (index_t aj=0; aj < numActive; ++aj)
                        {
                            const index_t jj = mappers[cj].index( actives(aj) ); 
                            
                            if ( mappers[cj].is_free_index(jj) )
                            {
                                sysMatrix.coeffRef(ii, jj) += localMat(ai,aj);
                            }
                        }
                }
            }
    }
    
    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

    /// Basis values
    gsMatrix<T>      basisData;
    gsMatrix<index_t> actives;
    index_t numActive;

    /// Local matrix
    gsMatrix<T> localMat;
    
protected:

    /// Material density
    T m_rho;
    
    /// Half the shell thickness
    T m_thickness;
};


} // namespace gismo

