/** @file gsVisitorBiharmonicModified.h

    @brief Visitor for a simple modified Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

// INCORRECT, WORK IN PROGRESS!!!

#pragma once

namespace gismo
{

/** \brief Visitor for the convection-diffusion-reaction equation.
 *
 * Visitor for PDEs of the form\n
 * Find \f$ u: \mathbb R^d \rightarrow \mathbb R\f$
 * \f[ \alpha (\Delta - 1)^2 u = f \quad in \quad \Omega,\f]
 * \f[ \frac{\partial u}{\partial \mathbf{n}} = 0 \quad on \quad \partial \Omega, \f]
 * \f[ \alpha\frac{\partial \nabla u}{\partial \mathbf{n}} = 0 \quad on \quad \partial \Omega \f], where\n
 * \f$ \alpha \f$ (reaction coefficient) is a scalar.
 *
 */

template <class T>
class gsVisitorBiharmonicModified
{
public:

    /** \brief Constructor for gsVisitorBiharmonicModified.
     *
     * \param[in] rhs Given right-hand-side function/source term that, for
     * \param[in] coeff_a Given coefficient
     */
    gsVisitorBiharmonicModified(const gsFunction<T> & rhs,
                                const gsFunction<T> & coeff_a) :
        rhs_ptr(&rhs),
        coeff_a_ptr( & coeff_a)
    {
        GISMO_ASSERT( rhs.targetDim() == 1 ,"Not yet tested for multiple right-hand-sides");
    }

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T>    & rule,
                    unsigned         & evFlags )
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// NB!

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T>            & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(quNodes.col(0), actives);
        numActive = actives.rows();

        //deriv2_into()
        //col(point) = B1_xx B2_yy B1_zz B_xy B1_xz B1_xy B2_xx ...

        // Evaluate basis functions on element
        basis.evalAllDers_into( quNodes, 2, basisData);

        //gsInfo << "Evaluation points :\n" << quNodes << std::endl;

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval.evaluateAt(quNodes);// is this generic ??

        // Evaluate the coefficients
        coeff_a_ptr->eval_into( geoEval.values(), coeff_a_vals ); // Dim: 1 X NumPts

        // Evaluate right-hand side at the geometry points
        rhs_ptr->eval_into( geoEval.values(), rhsVals ); // Dim: 1 X NumPts

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand sides
    }


    inline void assemble(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {

        const typename gsMatrix<T>::Block basisVals  = basisData.topRows(numActive);
        const typename gsMatrix<T>::Block basisGrads =
            basisData.middleRows( numActive, numActive * element.dim() );
        const typename gsMatrix<T>::Block basis2ndDerivs =
            basisData.bottomRows( numActive * (element.dim() * (element.dim()+1 ))/2 );

         //(ParDim + (ParDim*(ParDim-1))/2);
        const unsigned d = element.dim();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);

            // Compute physical gradients at k as a Dim x numActive matrix
            geoEval.transformGradients(k, basisGrads, physBasisGrad);
            // Compute physical laplacian at k as a 1 x numActive matrix
            geoEval.transformLaplaceHgrad(k, basisGrads, basis2ndDerivs, physBasisLaplace);

            // tmp_A         : d x d
            gsMatrix<T> tmp_a = coeff_a_vals.col(k);
            tmp_a.resize(d,d);

            // (u,v)
            // ( N x 1 ) * ( 1 x N) = N x N
            localMat.noalias() += weight * (basisVals.col(k) * basisVals.col(k).transpose());

            // (\nabla u, \nabla v)
            // ( N x d ) * ( d x N ) = N x N
            localMat.noalias() += 2 * weight * (physBasisGrad.transpose() * physBasisGrad);

            // (\Delta u, \Delta v)
            localMat.noalias() += weight * (physBasisLaplace.transpose() * physBasisLaplace);

            //(f,v)
            localRhs.noalias() += weight * ( basisVals.col(k) * rhsVals.col(k).transpose() ) ;
        }
    }

    inline void localToGlobal(const gsDofMapper     & mapper,
                              const gsMatrix<T>     & eliminatedDofs,
                              const index_t           patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
        // Local Dofs to global dofs
        mapper.localToGlobal(actives, patchIndex, actives);

        for (index_t i=0; i < numActive; ++i)
        {
            const int ii = actives(i);
            if ( mapper.is_free_index(ii) )
            {
                rhsMatrix.row(ii) += localRhs.row(i);

                for (index_t j=0; j < numActive; ++j)
                {
                    const int jj = actives(j);
                    if ( mapper.is_free_index(jj) )
                    {
                        sysMatrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else
                    {
                        rhsMatrix.row(ii).noalias() -= localMat(i, j) *
                            eliminatedDofs.row( mapper.global_to_bindex(jj) );
                    }
                }
            }
        }
    }


protected:
    // Right hand side
    const gsFunction<T> * rhs_ptr;
    // PDE Coefficient
    const gsFunction<T> * coeff_a_ptr;
    // flag for stabilization method
    unsigned flagStabType;

protected:
    // Basis values
    gsMatrix<T>        basisData;
    gsMatrix<T>        physBasisGrad;
    gsMatrix<T>        physBasisLaplace;
    gsMatrix<index_t> actives;
    index_t numActive;

    gsMatrix<T> coeff_a_vals;

protected:
    // Local values of the right hand side
    gsMatrix<T> rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
};


} // namespace gismo

