
/** @file gsVisitorPoissonHeterogeneous.h

    @brief Poisson equation element visitor.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#pragma once

#include<gsAssembler/gsVisitorPoisson.h>
#include <gsPde/gsPoissonHeterogeneousPde.h>

namespace gismo
{

template <class T>
class gsVisitorPoissonHeterogeneous
{
public:

    gsVisitorPoissonHeterogeneous(const gsPde<T> & pde)
    {
        const gsPoissonHeterogeneousPde<T> *ppde = static_cast<const gsPoissonHeterogeneousPde<T> * >(&pde);

        rhs_ptr = ppde->rhs();
        m_alpha = ppde->getAlpha();
    }

    /// Constructor with the right hand side function of the Poisson equation
    gsVisitorPoissonHeterogeneous(const gsFunction<T> & rhs, const gsPiecewiseFunction<T> & alpha) :
        rhs_ptr(&rhs), m_alpha(&alpha)
    { }

    gsVisitorPoissonHeterogeneous() :     rhs_ptr(NULL), m_alpha(NULL) { }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        // Setup Quadrature
        const real_t quA = options.getReal("quA");
        const index_t quB = options.getInt("quB");
        rule = gsGaussRule<T>(basis, quA, quB);// harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>    & basis, // to do: more unknowns
                         const gsGeometry<T> & geo,
                         const gsMatrix<T>   & quNodes)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);
        numActive = actives.rows();

        // Evaluate basis functions on element
        basis.evalAllDers_into(md.points, 1, basisData);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo.computeMap(md);


        // Evaluate right-hand side at the geometry points
        rhs_ptr->eval_into(md.values[0], rhsVals); // to do: parametric rhs ?

        m_alpha->piece(geo.id()).eval_into(md.values[0], alphaVals);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
    }

    inline void assemble(      gsDomainIterator<T> & ,
                         const gsVector<T>         & quWeights)
    {
        const gsMatrix<T> & bVals = basisData[0];
        const gsMatrix<T> & bGrads = basisData[1];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);

            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, physGrad);

            localRhs.noalias() += weight * (bVals.col(k) * rhsVals.col(k).transpose());
            localMat.noalias() += weight * (alphaVals(0, k) * physGrad.transpose() * physGrad);
        }
    }


    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, eliminatedDofs.front(), 0, 0);
    }

    inline void localToGlobal(const gsDofMapper     & mapper,
                              const gsMatrix<T>     & eliminatedDofs,
                              const index_t           patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
        //Assert eliminatedDofs.rows() == mapper.boundarySize()

        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        //const int numActive = actives.rows();

        for (index_t i = 0; i < numActive; ++i)
        {
            const int ii = actives(i);
            if (mapper.is_free_index(ii))
            {
                rhsMatrix.row(ii) += localRhs.row(i);

                for (index_t j = 0; j < numActive; ++j)
                {
                    const int jj = actives(j);
                    if (mapper.is_free_index(jj))
                    {
                        // Matrix is symmetric, store only lower triangular part
                        // if ( jj <= ii )
                        sysMatrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else // if ( mapper.is_boundary_index(jj) ) // Fixed DoF?
                    {
                        rhsMatrix.row(ii).noalias() -= localMat(i, j) *
                            eliminatedDofs.row(mapper.global_to_bindex(jj));
                    }
                }
            }
        }
    }

private:
    // Right hand side
    const gsFunction<T> * rhs_ptr;
    const gsPiecewiseFunction<T>* m_alpha;

    // Basis values
    std::vector<gsMatrix<T> >      basisData;
    gsMatrix<T>        physGrad;
    gsMatrix<index_t> actives;
    index_t numActive;

    // Local values of the right hand side
    gsMatrix<T> rhsVals;
    gsMatrix<T> alphaVals;

    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    gsMapData<T> md;
};


} // namespace gismo



