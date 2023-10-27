/** @file gsVisitorNLElasticityNeumann.h

    @brief Neumann conditions visitor for thin shells.

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
    @brief The visitor adds the Neumann force contribution to element right hand side for thin shells.
    
    \tparam T coefficient type
    
    \ingroup ThinShell
*/
template <class T>
class gsVisitorShellNeumann
{
public:

    /// Constructor with boundary force and boxside object as inputs.
    gsVisitorShellNeumann(const gsFunction<T> & neudata, boxSide s) : 
    neudata_ptr(&neudata), side(s)
    { }

    /// Function to initialize the assembly procedure.
    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T>    & rule, 
                    unsigned         & evFlags )
    {
        const int dir = side.direction();
        gsVector<index_t> numQuadNodes ( basis.dim() );
        for (short_t i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_JACOBIAN;
    }

    /// Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         // todo: add element here for efficiency
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

        // Evaluate the Neumann data
        neudata_ptr->eval_into(geoEval.values(), neuData);

        // Initialize local matrix/rhs
        localRhs.setZero(3*numActive, 1 );
    }

    /// Assembles the local right hand side using only the undeformed geometry.
    inline void assemble(gsDomainIterator<T>    & /*element*/,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector on the side
            //geoEval.outerNormal(k, side, unormal); // (!)
            //geoEval.normal(k,unormal);  
            
            // Multiply quadrature weight by the measure of normal on the boundary
            const T weight = quWeights[k] * geoEval.jacobian(k).col( ! side.direction() ).norm();
            
            for (index_t j = 0; j!= 3; ++j)
                localRhs.middleRows(j*numActive,numActive).noalias() += 
                    weight * neuData(j,k) * basisData.col(k) ;
        }
    }
    
    /// Putting local information in the global rhs.
    void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                       const gsMatrix<T>                 & /*eliminatedDofs*/,
                       const index_t                       /*patchIndex*/,
                       gsSparseMatrix<T>                 & /*sysMatrix*/,
                       gsMatrix<T>                       & rhsMatrix )
    {
        for (index_t ci = 0; ci!= 3; ++ci)
        {
            for (index_t ai=0; ai < numActive; ++ai)
            {
                const index_t gi = ci * numActive +  ai; // row index

                const int ii = mappers[ci].index(ai);

                if ( mappers[ci].is_free_index(ii) )
                {
                    rhsMatrix.row(ii) += localRhs.row(gi);
                }
            }
        }
    }

protected:

    index_t numActive;
    
    /// Pointer to the Neumann forces
    const gsFunction<T> * neudata_ptr;
    boxSide side;

    /// Basis values
    gsMatrix<T>      basisData;
    gsMatrix<index_t> actives;

    /// Neumann force values
    gsMatrix<T> neuData;

    /// Local matrix
    gsMatrix<T> localMat;
    
    /// Local right hand side
    gsMatrix<T> localRhs;

};


} // namespace gismo
