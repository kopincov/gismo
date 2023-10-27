/** @file gsVisitorNeumann2.h

    @brief Neumann conditions visitor for elliptic problems.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>

namespace gismo
{
/** @brief
    Implementation of a Neumann BC for elliptic Assembler.

    It sets up an assembler and adds the following term to the linear term.
    \f[ \nabla u \cdot \mathbf{n} = g_N  \f]
*/

template <class T>
class gsVisitorNeumann2
{
public:

    gsVisitorNeumann2(const gsPde<T> & , const boundary_condition<T> & s)
    : neudata_ptr( s.function().get() ), side(s.side())
    { }

/** @brief
 * Constructor of the assembler object 
 * 
   \param[in] neudata is the Neumann boundary data.
   \param[in] s are the sides of the geometry where neumann BC are prescribed.
   
   \f[ \nabla u \cdot \mathbf{n} = g_N  \f]
*/
    gsVisitorNeumann2(const gsFunction<T> & neudata, boxSide s) : 
    neudata_ptr(&neudata), side(s)
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
        evFlags = NEED_VALUE|NEED_JACOBIAN;
    }

    void initialize(const gsBasis<T> & basis,
                    const index_t ,
                    const gsOptionList & options, 
                    gsQuadRule<T>    & rule,
                    unsigned         & evFlags )
    {
        // Setup Quadrature (harmless slicing occurs)
        rule = gsQuadrature::get(basis, options, side.direction());

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
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
        basis.eval_into(quNodes, basisData);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate the Neumann data
        neudata_ptr->eval_into(geoEval.values(), neuData);

        // Initialize local matrix/rhs
        localRhs.setZero(numActive, neudata_ptr->targetDim() );
    }

    inline void assemble(gsDomainIterator<T>    & ,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, side, unormal);
            
            // Multiply quadrature weight by the measure of normal
            const T weight = quWeights[k] * unormal.norm();
            
            localRhs.noalias() += weight * basisData.col(k) * neuData.col(k).transpose() ;
        }
    }
    
    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> >   & ,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.pushToRhs(localRhs, actives, 0);
    }

    void localToGlobal(const gsDofMapper  & mapper,
                       const gsMatrix<T>     & eliminatedDofs,
                       const index_t patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();

        // Push element contribution to the global load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            // convert local DoF index to global DoF index
            const unsigned jj = actives(j);
            if (mapper.is_free_index(jj))
                rhsMatrix.row(jj) += localRhs.row(j);
        }
    }

protected:

    
    // Neumann function
    const gsFunction<T> * neudata_ptr;
    boxSide side;

    // Basis values
    gsMatrix<T>      basisData;
    gsMatrix<index_t> actives;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> neuData;

    // Local matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

};


} // namespace gismo
