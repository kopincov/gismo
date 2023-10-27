/** @file gsVisitorNitscheSpaceTime.h

    @brief Weak (Nitsche-type) BC imposition visitor for heat equation.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Moore
*/

#pragma once

namespace gismo
{
/** \brief Visitor for the weak imposition of the dirichlet boundary condition.
 *
 * Adds this term to the bilinear terms
 * \f[ (\nabla_x u, v)_{\partial \Omega} + (u, \nabla_x v )_{\partial \Omega} 
 *                                     + (\mu*u, v)_{\partial \Omega} \f]
 * 
 * The following term is also added to the linear form
 * \f[ (g_D, \mu*v + \nabla_x v)_{\partial \Omega} \f],
 * where the dirichlet term is given as \f[ g_D \f].
 */
template <class T>
class gsVisitorNitscheSpaceTime
{
public:
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
    gsVisitorNitscheSpaceTime(const gsFunction<T> & dirdata, T _penalty, boxSide s) : 
    dirdata_ptr(&dirdata),penalty(_penalty), side(s)
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
        evFlags = NEED_VALUE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;
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

        // Evaluate basis values and derivatives on element
        basis.evalAllDers_into( quNodes, 1, basisData);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate the Dirichlet data
        dirdata_ptr->eval_into(geoEval.values(), dirData);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, dirdata_ptr->targetDim() );
    }

    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        const unsigned d = element.dim();
        
        const gsMatrix<T> & basisVals  = basisData[0];
        const gsMatrix<T> & basisGrads = basisData[1];
        
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, side, unormal);

            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] *unormal.norm();   

            // Compute the unit normal vector 
            // Dim x Num Active Basis
            unormal.normalize();
            
            // Compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGrads, physGrad);
        
            physGradSpace   = physGrad.topRows(d-1);
            physGradTime    = physGrad.bottomRows(1);
            unormalSpace    = unormal.topRows(d-1);

     
            // Get penalty parameter
            const T mu = penalty / element.getCellSize();

            // Sum up quadrature point evaluations
            localRhs.noalias() += weight * (( physGradSpace.transpose() * unormalSpace - mu *basisVals )
                                            * dirData.col(k).transpose() );

            localMat.noalias() += weight * (basisVals * unormalSpace.transpose() * physGradSpace
                            +  (basisVals * unormalSpace.transpose() * physGradSpace).transpose()
                            -  mu *basisVals *basisVals.transpose() );
            }
    }
    
    void localToGlobal(const gsDofMapper  & mapper,
                       const gsMatrix<T>  & eliminatedDofs,
                       const index_t        patchIndex,
                       gsSparseMatrix<T>  & sysMatrix,
                       gsMatrix<T>        & rhsMatrix )
    {
        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();

        // Push element contribution to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = actives(j);
            rhsMatrix.row(jj) -= localRhs.row(j);
            for (index_t i=0; i!= numActive; ++i)
            {
                const unsigned ii = actives(i);
//                if ( jj <= ii ) // assuming symmetric problem
                    sysMatrix( ii, jj ) -= localMat(i,j);
            }
        }

    }

private:
    // Dirichlet function
    const gsFunction<T> * dirdata_ptr;

    // Penalty constant
    T penalty;

    // Side
    boxSide side;

private:
    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T>      physGrad;
    gsMatrix<index_t> actives;

    // Normal and Neumann values
    gsVector<T> unormal, unormalSpace;
    gsMatrix<T> dirData;

     gsMatrix<T> physGradTime, physGradSpace;
    // Local  matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

};


} // namespace gismo
