/** @file gsVisitorNitsche2.h

    @brief Weak (Nitsche-type) BC imposition visitor for elliptic problems.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris , S. Moore, C.Hofer, J.Vogl
*/

#include <gsPde/gsPoissonHeterogeneousPde.h>
#pragma once


namespace gismo
{
/** \brief Visitor for the weak imposition of the dirichlet boundary condition.
 * This Visitor includes the use of inhomogeneous diffusion coefficient
 *
 * Adds this term to the bilinear terms
 * \f[ (\nabla u, v)_{\partial \Omega} + (u, \nabla v )_{\partial \Omega}
 *                                     + (\mu*u, v)_{\partial \Omega} \f]
 *
 * The following term is also added to the linear form
 * \f[ (g_D, \mu*v + \nabla v)_{\partial \Omega} \f],
 * where the dirichlet term is given as \f[ g_D \f].
 */

template <class T>
class gsVisitorNitsche2
{
public:

    gsVisitorNitsche2(const gsPde<T> & pde, const boundary_condition<T> & s)
        : m_domain(&pde.domain()), dirdata_ptr( s.function().get() ), side(s.side()), isHomogeneous(false)
    {
        const gsPoissonHeterogeneousPde<T>* ppde = static_cast<const gsPoissonHeterogeneousPde<T>*>(&pde);
        if(ppde == NULL)
            isHomogeneous = true;
        else
            m_alpha = ppde->getAlpha();
    }

    /** @brief
    Constructor of the assembler object.

    \param[in] basis a multi-basis that contains patch-wise basis
    \param[in] dirData  is a gsBoundaryConditions object that holds boundary conditions of the form:
    \f[ \text{Dirichlet: } u = g_D \text{ on } \Gamma.\f]
    \f$ v \f$ is the test function and \f$ \Gamma\f$ is the boundary side.
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary in the \em bconditions object.
    \param[in] penalty for inputing the penalty choice
*/
    gsVisitorNitsche2(const gsFunction<T> & dirdata, T _penalty, boxSide s, const gsPiecewiseFunction<T> * alpha) :
        dirdata_ptr(&dirdata),penalty(_penalty), side(s), m_alpha(alpha), isHomogeneous(false)
    { }

    gsVisitorNitsche2(const gsFunction<T> & dirdata, T _penalty, boxSide s) :
        dirdata_ptr(&dirdata),penalty(_penalty), side(s), isHomogeneous(true)
    {
    }

    void initialize(const gsBasis<T>    & basis,
                    const index_t         patchIndex,
                    const gsOptionList  & options,
                          gsQuadRule<T> & rule)
    {
        d = basis.dim();
        const T quA = options.getReal("quA") + 1;
        const index_t quB = options.getInt ("quB");
        rule = gsGaussRule<T>(basis, quA, quB, side.direction() );

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        // Compute penalty parameter
        const int deg = basis.maxDegree();
        penalty = (deg + basis.dim()) * (deg + 1) * T(2.5);

        if(m_domain==NULL)
            m_patchDiameter = 1;
        else
             m_patchDiameter = (m_domain->patch(patchIndex).coefAtCorner(boxCorner::getFirst(1)) - m_domain->patch(patchIndex).coefAtCorner(boxCorner::getLast(1))).norm();
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>    & basis,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T>   & quNodes)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(md.points.col(0), actives);
        const index_t numActive = actives.rows();

        // Evaluate basis values and derivatives on element
        basis.evalAllDers_into(md.points, 1, basisData);

        // Compute geometry related values
        geo.computeMap(md);

        // Evaluate the Dirichlet data
        dirdata_ptr->eval_into(md.values[0], dirData);

        if (!isHomogeneous)
            m_alpha->piece(geo.id()).eval_into(md.values[0], alphaVals);
        else
            alphaVals = Eigen::Matrix<T, Dynamic, Dynamic>::Constant(1, md.points.cols(), 1);
        // penalty=std::max(penalty,penalty*((alphaVals.minCoeff()*alphaVals.minCoeff())/alphaVals.maxCoeff() ));

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, dirdata_ptr->targetDim());
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T>   & quWeights)
    {
        gsMatrix<T> & bGrads = basisData[1];
        const index_t numActive = actives.rows();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            const typename gsMatrix<T>::Block bVals =
                basisData[0].block(0, k, numActive, 1);

            // Compute the outer normal vector on the side
            outerNormal(md, k, side, unormal);

            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] * unormal.norm();

            // Compute the unit normal vector
            unormal.normalize();

            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, pGrads);

            // Get penalty parameter
            const T mu = (penalty / (element.getCellSize() * m_patchDiameter)) * alphaVals(0, k);

            // Sum up quadrature point evaluations
            localRhs.noalias() -= weight * ((alphaVals(0, k) * pGrads.transpose() * unormal - mu * bVals)
                * dirData.col(k).transpose());

            localMat.noalias() -= weight * (alphaVals(0, k) * bVals * unormal.transpose() * pGrads
                + (alphaVals(0, k) * bVals * unormal.transpose() * pGrads).transpose()
                - mu * bVals * bVals.transpose());
        }
    }

    inline void localToGlobal(const index_t                     patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                                    gsSparseSystem<T>         & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.pushAllFree(localMat, localRhs, actives, 0);
    }

    void localToGlobal(const gsDofMapper       & mapper,
                       const gsMatrix<T>       & eliminatedDofs,
                       const index_t             patchIndex,
                             gsSparseMatrix<T> & sysMatrix,
                             gsMatrix<T>       & rhsMatrix )
    {
        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();

        // Push element contribution to the global matrix and load vector
        for (index_t j = 0; j != numActive; ++j)
        {
            const unsigned jj = actives(j);
            rhsMatrix.row(jj) += localRhs.row(j);
            for (index_t i = 0; i != numActive; ++i)
            {
                const unsigned ii = actives(i);
                // if ( jj <= ii ) // assuming symmetric problem
                sysMatrix(ii, jj) += localMat(i, j);
            }
        }
    }

private:
    const gsMultiPatch<T>* m_domain;

    // Dirichlet function
    const gsFunction<T> * dirdata_ptr;

    // Penalty constant
    T penalty;

    // Side
    boxSide side;

    const gsPiecewiseFunction<T> * m_alpha;

    bool isHomogeneous;

    int d;

    T m_patchDiameter;

    // Basis values
    std::vector<gsMatrix<T> >     basisData;
    gsMatrix<T>      pGrads;
    gsMatrix<index_t> actives;
    gsMatrix<T>      alphaVals;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> dirData;

    // Local matrix, rhs and md
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    gsMapData<T> md;
};


} // namespace gismo
