/** @file gsVisitorFlux.h

    @brief Visitor for a nonlinear FluxPde

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#pragma once
#include <gsPde/gsFluxPde.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/** \brief Visitor for a generic Flux bilinearform
 *
 * Assembles the bilinear terms
 * \f[ (F'(\nabla u)\nabla z,\nabla v)_\Omega \text{ and } (f,v)_\Omega - (F(\nabla u),\nabla v)\f]
 * For \f[ u = g \quad on \quad \partial \Omega \f], where u is the given current solution.
 *
 */

template <class T, bool paramCoef = false>
class gsVisitorFlux
{
public:

    gsVisitorFlux(const gsPde<T> & pde): flux(static_cast<const gsFluxPde<T>& >(pde).getFlux())
    {
        const gsFluxPde<T>& fpde = static_cast<const gsFluxPde<T>& >(pde);
        d = fpde.dim();
        rhs_ptr = fpde.rhs();
    }

    gsVisitorFlux() :     rhs_ptr(NULL) { }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        // Setup Quadrature
        const real_t quA = options.getReal("quA");
        const index_t quB = options.getInt ("quB");
        rule = gsGaussRule<T>(basis, quA, quB);// harmless slicing occurs here

        patch = patchIndex;

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>    & basis,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T>   & quNodes)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);
        numActive = actives.rows();

        // Evaluate basis functions on element
        basis.evalAllDers_into( md.points, 1, basisData);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo.computeMap(md);

        //Give the flux the MapData (for pyhsical quadrature points and transformation)
        flux[patch]->setMapData(md);

        // Calculate Value of flux
        flux[patch]->eval_into(md.points, fluxVal);

        //calculate Linearization of Flux //geoEval.id()
        flux[patch]->deriv_into(md.points, linFlux);

        // Evaluate right-hand side at the geometry points paramCoef
        // specifies whether the right hand side function should be
        // evaluated in parametric(true) or physical (false)
        rhs_ptr->eval_into( (paramCoef ?  md.points :  md.values[0] ), rhsVals );

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand sides
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);

            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, physGrad);

            //rhs  = -a(u_old,phi)+(f,phi)
            const gsMatrix<T>& a_old = physGrad.transpose()*fluxVal.col(k);
            const gsMatrix<T>& rhs = bVals.col(k) * rhsVals.col(k).transpose();
            localRhs.noalias() += weight * (-a_old + rhs) ;

            //K = F'(u)\grad Phi * \grad \Phi
            //The linearized Flux at a quadrature point must bereshaped into a square dxd matrix
            const gsMatrix<T>& FgradZ=linFlux.reshapeCol(k,d,d).transpose()*physGrad;
            localMat.noalias() += weight * (FgradZ.transpose() * physGrad);
        }

        //  gsInfo<<"rhs:\n "<<localRhs.transpose()<<"\n";
        //  gsInfo<<"mat:\n "<<localMat<<"\n";
    }

    /* start old assembler */
//    void initialize(const gsBasis<T> & basis,
//                    const index_t patchIndex,
//                    const gsOptionList & options,
//                    gsQuadRule<T>    & rule,
//                    unsigned         & evFlags )
//    {
//        // Setup Quadrature
//        const real_t quA = options.getReal("quA");
//        const index_t quB = options.getInt ("quB");
//        rule = gsGaussRule<T>(basis, quA, quB);// harmless slicing occurs here
//
//        patch = patchIndex;
//
//        // Set Geometry evaluation flags
//        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
//    }
//
//    // Evaluate on element.
//    inline void evaluate(gsBasis<T> const       & basis,
//                         gsGeometryEvaluator<T> & geoEval,
//                         gsMatrix<T> const      & quNodes)
//    {
//        // Compute the active basis functions
//        // Assumes actives are the same for all quadrature points on the elements
//        basis.active_into(quNodes.col(0), actives);
//        numActive = actives.rows();
//
//        // Evaluate basis functions on element
//        basis.evalAllDers_into( quNodes, 1, basisData);
//
//        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
//        geoEval.evaluateAt(quNodes);// is this generic ??
//
//        //Give the flux the geometry Evaluator (for pyhsical quadrature points and transformation)
//        flux[patch]->setGeometryEvaluator(geoEval);
//
//        // Calculate Value of flux
//        flux[patch]->eval_into(quNodes, fluxVal);
//
//        //calculate Linearization of Flux //geoEval.id()
//        flux[patch]->deriv_into(quNodes, linFlux);
//
//        // Evaluate right-hand side at the geometry points paramCoef
//        // specifies whether the right hand side function should be
//        // evaluated in parametric(true) or physical (false)
//        rhs_ptr->eval_into( (paramCoef ?  quNodes :  geoEval.values() ), rhsVals );
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
//        gsMatrix<T> & bVals  = basisData[0];
//        gsMatrix<T> & bGrads = basisData[1];
//
//        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
//        {
//            // Multiply weight by the geometry measure
//            const T weight = quWeights[k] * geoEval.measure(k);
//
//            // Compute physical gradients at k as a Dim x NumActive matrix
//            geoEval.transformGradients(k, bGrads, physGrad);
//
//            //rhs  = -a(u_old,phi)+(f,phi)
//            const gsMatrix<T>& a_old = physGrad.transpose()*fluxVal.col(k);
//            const gsMatrix<T>& rhs = bVals.col(k) * rhsVals.col(k).transpose();
//            localRhs.noalias() += weight * (-a_old + rhs) ;
//
//            //K = F'(u)\grad Phi * \grad \Phi
//            //The linearized Flux at a quadrature point must bereshaped into a square dxd matrix
//            const gsMatrix<T>& FgradZ=linFlux.reshapeCol(k,d,d).transpose()*physGrad;
//            localMat.noalias() += weight * (FgradZ.transpose() * physGrad);
//        }
//
//      //  gsInfo<<"rhs:\n "<<localRhs.transpose()<<"\n";
//      //  gsInfo<<"mat:\n "<<localMat<<"\n";
//    }
    /* end old assembler */

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
    int d;

    size_t patch;
    // Right hand side
    const gsFunction<T> * rhs_ptr;
    const std::vector<gsFlux<T>* >& flux;

protected:
    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T>        physGrad;
    gsMatrix<index_t> actives;
    gsMatrix<T>        linFlux;
    gsMatrix<T>        fluxVal;
    index_t numActive;

protected:
    // Local values of the right hand side
    gsMatrix<T> rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    gsMapData<T> md;
};


} // namespace gismo

