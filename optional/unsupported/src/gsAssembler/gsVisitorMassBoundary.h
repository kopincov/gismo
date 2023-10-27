/** @file gsVisitorMassBoundary.h

    @brief Mass visitor for assmbling element mass matrix on the boundary

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J.Sogn
*/

#pragma once

#include <gsAssembler/gsVisitorMass.h>

namespace gismo
{

template <class T>
class gsVisitorMassBoundary : public gsVisitorMass<T>
{
public:

    gsVisitorMassBoundary(boxSide s) :
    side(s)
    { }

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned & evFlags  )
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

    inline void assemble(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        const index_t numActive = actives.rows();
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const typename gsMatrix<T>::Block bVals  = basisData.block(0,k,numActive,1);

            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, side, unormal);

            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] * unormal.norm();
            //fixme: is this NOT the same as measure(k) ..
            
            // Sum up quadrature point evaluations
            localMat.noalias() += weight * ( bVals * bVals.transpose() );
        }
    }

    // Not same as gsVisitorMass (I don't use symmertri properties:( )
    void localToGlobal(const gsDofMapper  & mapper,
                       const gsMatrix<T>     & eliminatedDofs,
                       const index_t patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        // Local Dofs to global dofs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();
        // Push element contribution to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = actives(j);
            for (index_t i=0; i!=numActive; ++i)
            {
                const unsigned ii = actives(i);
                sysMatrix.coeffRef(ii, jj) += localMat(i, j);
            }
        }
    }

private:
    // Side
    boxSide side;

private:
    // Basis values
    using gsVisitorMass<T>::basisData;
    using gsVisitorMass<T>::actives;

    // Normal vector
    gsVector<T> unormal;

    // Local  matrix
    using gsVisitorMass<T>::localMat;

};


} // namespace gismo
