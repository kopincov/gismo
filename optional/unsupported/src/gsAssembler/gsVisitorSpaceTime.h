/** @file gsVisitorSpaceTime.h

    @brief Heat equation element visitor.

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

/** \brief Visitor for the Heat equation.
 *
 * Assembles the bilinear terms
 * \f[ (\partial_t u- \Delta u = f  \f]
 * For \f[ u = g \quad on  \quad \Sigma \f],
 *     \f[ u(x,0) = u_0(x) \quad on \quad \Sigma_0 \f],
 */

template <class T, bool paramCoef = false>
class gsVisitorSpaceTime
{
public:
    
    /** \brief Constructor for gsVisitorSpaceTime.
     *
     * \param[in] rhs Given right-hand-side function/source term that, for
     * \param[in] theta Given coefficient
     */
    /// Constructor with the right hand side function of the Poisson equation
    gsVisitorSpaceTime(const gsPde<T> & pde)
    { 
        const gsPoissonHeterogeneousPde<T>* pde_ptr
            = static_cast<const gsPoissonHeterogeneousPde<T>*>(&pde);
        rhs_ptr = pde_ptr->rhs();
        theta_ptr = pde_ptr->diffusion();
    }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        // Setup Quadrature
        rule = gsGaussRule<T>(basis, options);// harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
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
        basis.evalAllDers_into( md.points, 2, basisData );

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo.computeMap(md);

        // Evaluate right-hand side at the geometry points paramCoef
        // specifies whether the right hand side function should be
        // evaluated in parametric(true) or physical (false)
        rhs_ptr->eval_into( (paramCoef ?  md.points :  md.values[0] ), rhsVals );
        // Evaluate the coefficient
        theta_ptr->eval_into( md.values[0], theta_vals );

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand sides
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        const unsigned d = element.dim();

        const gsMatrix<T> & basisVals      = basisData[0];
        const gsMatrix<T> & basisGrads     = basisData[1];
        const gsMatrix<T> & basis2ndDerivs = basisData[2];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);

            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, basisGrads, physGrad);
            transformDeriv2Hgrad(md, k, basisGrads, basis2ndDerivs, physBasisLaplace);


            // Physical gradients in only (Dim)-coordinates at k 
            // is a (Dim-1) x NumActive matrix    
            physGradSpace   = physGrad.topRows(d-1);
            physGradTime    = physGrad.bottomRows(1);
            physMixedDeriv  = physBasisLaplace.rightCols(d-1);
            const T h     = element.getCellSize();
            //const T theta = 1.0; 

            localRhs.noalias() += weight * (theta_vals(0,k) *  h* ( physGradTime.transpose() * rhsVals.col(k).transpose() )
                + basisVals.col(k) * rhsVals.col(k).transpose() ) ;
            localMat.noalias() += weight *  ( theta_vals(0,k) * h* (physGradTime.transpose() * physGradTime
                -  (physMixedDeriv*physGradSpace).transpose() )
                +  physGradSpace.transpose() * physGradSpace
                +  basisVals.col(k) * physGradTime);

        }
    }

    /* gsAssembler2.h */
//    void initialize(const gsBasis<T> & basis,
//                    const index_t patchIndex,
//                    const gsOptionList & options,
//                    gsQuadRule<T>    & rule,
//                    unsigned         & evFlags )
//    {
//        // Setup Quadrature
//        rule = gsGaussRule<T>(basis, options);// harmless slicing occurs here
//
//        // Set Geometry evaluation flags
//        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
//    }
//
//    // Evaluate on element.
//    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
//                         gsGeometryEvaluator<T> & geoEval,
//                         gsMatrix<T> const      & quNodes)
//    {
//        // Compute the active basis functions
//        // Assumes actives are the same for all quadrature points on the elements
//        basis.active_into(quNodes.col(0), actives);
//        numActive = actives.rows();
//
//        // Evaluate basis functions on element
//        basis.evalAllDers_into( quNodes, 2, basisData );
//
//        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
//        geoEval.evaluateAt(quNodes);// is this generic ??
//
//        // Evaluate right-hand side at the geometry points paramCoef
//        // specifies whether the right hand side function should be
//        // evaluated in parametric(true) or physical (false)
//        rhs_ptr->eval_into( (paramCoef ?  quNodes :  geoEval.values() ), rhsVals );
//        // Evaluate the coefficient
//        theta_ptr->eval_into( geoEval.values(), theta_vals );
//
//        // Initialize local matrix/rhs
//        localMat.setZero(numActive, numActive      );
//        localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand sides
//    }
//
//    inline void assemble(gsDomainIterator<T>    & element,
//                         gsGeometryEvaluator<T> & geoEval,
//                         gsVector<T> const      & quWeights)
//    {
//        const unsigned d = element.dim();
//
//        const gsMatrix<T> & basisVals      = basisData[0];
//        const gsMatrix<T> & basisGrads     = basisData[1];
//        const gsMatrix<T> & basis2ndDerivs = basisData[2];
//
//        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
//        {
//            // Multiply weight by the geometry measure
//            const T weight = quWeights[k] * geoEval.measure(k);
//
//            // Compute physical gradients at k as a Dim x NumActive matrix
//            geoEval.transformGradients(k, basisGrads, physGrad);
//            geoEval.transformDeriv2Hgrad(k, basisGrads, basis2ndDerivs, physBasisLaplace);
//
//
//            // Physical gradients in only (Dim)-coordinates at k
//            // is a (Dim-1) x NumActive matrix
//            physGradSpace   = physGrad.topRows(d-1);
//            physGradTime    = physGrad.bottomRows(1);
//            physMixedDeriv  = physBasisLaplace.rightCols(d-1);
//            const T h     = element.getCellSize();
//            //const T theta = 1.0;
//
//            localRhs.noalias() += weight * (theta_vals(0,k) *  h* ( physGradTime.transpose() * rhsVals.col(k).transpose() )
//                + basisVals.col(k) * rhsVals.col(k).transpose() ) ;
//            localMat.noalias() += weight *  ( theta_vals(0,k) * h* (physGradTime.transpose() * physGradTime
//                -  (physMixedDeriv*physGradSpace).transpose() )
//                +  physGradSpace.transpose() * physGradSpace
//                +  basisVals.col(k) * physGradTime);
//
//        }
//    }

    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, eliminatedDofs.front(), 0, 0);
    }

protected:
    // Right hand side
    const gsFunction<T> * rhs_ptr;
    const gsFunction<T> * theta_ptr;

protected:
    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T>        physGrad, physBasisLaplace;
    gsMatrix<index_t> actives;
    index_t numActive;

    gsMatrix<T> physGradTime, physGradSpace, physMixedDeriv;

protected:
    // Local values of the right hand side
    gsMatrix<T> rhsVals;
    gsMatrix<T> theta_vals;
    boxSide side1;

protected:
    // Local matrices
    gsMatrix<T> localMat, localMass;
    gsMatrix<T> localRhs;
    gsMapData<T> md;
};

} // namespace gismo
