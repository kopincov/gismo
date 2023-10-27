
#pragma once

#include<gsPde/gsStokesPde.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/** @brief Visitor for the Stokes problem using Taylor-Hood elements
 * and a common basis for all the components of the velocity.

   ( A   0  0   B1^T )      (rhs_u1)
   ( 0   A  0   B2^T )  and (rhs_u2)
   ( 0   0  A   B3^T )      (rhs_u3)
   ( B1  B2  B3  0   )      (rhs_p )
*/
template <class T>
class gsVisitorStokes
{
public:

    gsVisitorStokes(const gsPde<T> & pde)
    {
        const gsStokesPde<T>* spde = static_cast<const gsStokesPde<T>* >(&pde);
        rhs_ptr =spde->rhs() ;
        m_viscosity =spde->viscocity();

    }


    /// Constructor with the right and side function of the poisson equation
    gsVisitorStokes(const gsFunction<T> & rhs, const T viscosity) :
    rhs_ptr(&rhs),
    m_viscosity(viscosity)
    { }

    /* start code for stable/src/gsAssembler/gsAssemblerBase.h */
    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        d = basis.dim();

        // Setup Quadrature
        rule = gsGaussRule<T>(basis, options);// harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        // Create the d local B_i matrices
        localMatB.resize(d);
    }

    // initialized for each patch
    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T>    & rule)
    {
        d = basis.dim();
        gsVector<index_t> numQuadNodes(d);
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            //numQuadNodes[i] = basis.degree(i) + 1;
            numQuadNodes[i] = basis.maxDegree() + 1; // take quadrature from highest degree

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        // Create the d local B_i matrices
        localMatB.resize(d);
    }

    // Evaluate on element.
    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T>  & geo,
                         const gsMatrix<T>    & quNodes)
    {
        md.points = quNodes;

        // Evaluate basis functions on element  nodes
        basisRefs.front().active_into(md.points.col(0), actives_u);
        basisRefs.back().active_into(md.points.col(0), actives_p);
        const index_t numActU = actives_u.rows();
        const index_t numActP = actives_p.rows();
        basisRefs.front().evalAllDers_into(md.points, 1, basisData_u);
        basisRefs.back().eval_into(md.points, basisVals_p);

        // Evaluate Geometry on element nodes
        geo.computeMap(md);

        // Evaluate right-hand side at the geometry points
        rhs_ptr->eval_into(md.values[0], rhsVals);

        localMatA.setZero(numActU, numActU);//local_A
        localRhs_u.setZero(numActU, rhsVals.rows());
        for (index_t i = 0; i != d; ++i)
            localMatB[i].setZero(numActP, numActU);//local_B_i
    }
    
    inline void assemble(gsDomainIterator<T>    & element,
                         gsVector<T> const      & quWeights)
    {
        //const index_t numActP = actives_p.rows();
        const gsMatrix<T> & basisVals_u = basisData_u[0];
        const gsMatrix<T> & bGrads_u = basisData_u[1];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * md.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads_u, physGrad_u);

            // Right-hand side
            localRhs_u += weight * (basisVals_u.col(k) * rhsVals.col(k).transpose());

            // Local block A
            localMatA.noalias() += weight * m_viscosity * (physGrad_u.transpose() * physGrad_u);

            // Local blocks B_i
            for (index_t i = 0; i != d; ++i)
                localMatB[i].noalias() += weight * (basisVals_p.col(k) * physGrad_u.row(i));
        }
    }
    /* end code for stable/src/gsAssembler/gsAssemblerBase.h */

    /* start code for devel/src/gsAssembler/gsAssemblerBase2.h */
    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule,
                    unsigned         & evFlags )
    {
        d  = basis.dim();

        // Setup Quadrature
        rule = gsGaussRule<T>(basis, options);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        // Create the d local B_i matrices
        localMatB.resize(d);
    }

    // initialized for each patch
    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T>    & rule,
                    unsigned         & evFlags )
    {
        d  = basis.dim();
        gsVector<index_t> numQuadNodes( d );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            //numQuadNodes[i] = basis.degree(i) + 1;
            numQuadNodes[i] = basis.maxDegree() + 1; // take quadrature from highest degree

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_MEASURE|NEED_GRAD_TRANSFORM;

        // Create the d local B_i matrices
        localMatB.resize(d);
    }

    // Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        basisRefs.front().active_into(quNodes.col(0), actives_u);
        basisRefs.back ().active_into(quNodes.col(0), actives_p);
        const index_t numActU = actives_u.rows();
        const index_t numActP = actives_p.rows();
        basisRefs.front().evalAllDers_into(quNodes, 1, basisData_u);
        basisRefs.back() .eval_into       (quNodes,    basisVals_p);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate right-hand side at the geometry points
        rhs_ptr->eval_into( geoEval.values(), rhsVals );

        localMatA .setZero(numActU, numActU);//local_A
        localRhs_u.setZero(numActU, rhsVals.rows());
        for ( index_t i = 0; i!=d; ++i)
            localMatB[i].setZero(numActP, numActU);//local_B_i
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        //const index_t numActP = actives_p.rows();
        const gsMatrix<T> & basisVals_u  = basisData_u[0];
        const gsMatrix<T> & bGrads_u     = basisData_u[1];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_u, physGrad_u);

            // Right-hand side
            localRhs_u +=  weight * ( basisVals_u.col(k) *  rhsVals.col(k).transpose() );

            // Local block A
            localMatA.noalias()  += weight * m_viscosity * (physGrad_u.transpose() * physGrad_u);

            // Local blocks B_i
            for ( index_t i = 0; i!=d; ++i)
                localMatB[i].noalias() += weight * ( basisVals_p.col(k) * physGrad_u.row(i) );
        }
    }
    /* end code for devel/src/gsAssembler/gsAssemblerBase2.h */
    
    inline void localToGlobal(const index_t                     patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                                    gsSparseSystem<T>         & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives_u, patchIndex, actives_u,0);
        system.mapColIndices(actives_p, patchIndex, actives_p,d);

        // Add contributions to the system matrix and right-hand side
        for(index_t i=0; i!=d;++i)
        {
            system.push(localMatA, localRhs_u.col(i), actives_u, eliminatedDofs[i], i, i);
            system.pushToMatrix(localMatB[i], actives_p, actives_u, eliminatedDofs[i], d, i);
            system.pushToMatrix(localMatB[i].transpose(), actives_u, actives_p, eliminatedDofs[d], i, d);

        }
    }

    inline void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                              const gsMatrix<T>     & eliminatedDofs,
                              const index_t           patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
        GISMO_ASSERT(mappers.size() ==2,
                     "Expecting one Dof mapper for pressure and one for (all components of) velocity.");

        const gsDofMapper & Umap = mappers[0];
        const index_t usz = Umap.freeSize();
        const index_t ps = d*usz;

        // Local Dofs to global dofs
        mappers.back().localToGlobal(actives_p, patchIndex, actives_p);
        const index_t numActP = actives_p.rows();
        Umap.localToGlobal(actives_u, patchIndex, actives_u);
        const index_t numActU = actives_u.rows();

        for (index_t i=0; i < numActU; ++i)
        {
            const int ii = actives_u(i);
            if ( Umap.is_free_index(ii) )
            {
                for (index_t s = 0; s!=d; ++s) // todo, take out
                    rhsMatrix(ii+s*usz, 0) += localRhs_u(i,s);
                
                for (index_t j=0; j < numActP; ++j) // Build B-part of the matrix
                {
                    const int jj = actives_p(j);

                    for (index_t s = 0; s!=d; ++s)
                    {
                        sysMatrix.coeffRef(ps+jj   , ii+s*usz ) += localMatB[s](j, i);
                        sysMatrix.coeffRef(ii+s*usz, ps+jj    ) += localMatB[s](j, i);
                    }
                }

                for (index_t j=0; j < numActU; ++j) // Build A-part of the matrix
                {
                    const int jj = actives_u(j);
                    if ( Umap.is_free_index(jj) )
                    {
                        //   if ( jj <= ii ) // store only lower triangular part
                        for (index_t s = 0; s!=d; ++s)
                            sysMatrix.coeffRef(ii+s*usz, jj+s*usz) += localMatA(i, j);
                    }
                    else // Umap.is_boundary_index(jj)
                    {
                        const int bb = Umap.global_to_bindex(jj);
                        for (index_t s = 0; s!=d; ++s)
                            rhsMatrix(ii+s*usz,0) -= // assuming single rhs
                                localMatA(i,j) * eliminatedDofs(bb,s);
                    }
                }
            }
            else // Umap.is_boundary_index(ii)
            {
                const int bb = Umap.global_to_bindex(ii);
                for (index_t k=0; k < numActP; ++k)
                {
                    const int kk = actives_p(k);
                    T tmp = localMatB[0](k, i)*eliminatedDofs(bb,0);
                    for (index_t s = 1; s!=d; ++s)
                        tmp += localMatB[s](k, i) * eliminatedDofs(bb,s);
                    rhsMatrix(ps+kk,0) -= tmp;// assuming single rhs
               }
            }
        }
    }

protected:

    // Velocity vector dimension
    index_t d;

    // Right hand side
    const gsFunction<T> * rhs_ptr;

    T m_viscosity;

protected:
    // Basis values
    std::vector<gsMatrix<T> > basisData_u;
    gsMatrix<T> basisVals_p;
    gsMatrix<T> physGrad_u;
    gsMatrix<index_t> actives_u, actives_p;

protected:
    // Local matrices
    gsMatrix<T> localMatA;
    std::vector<gsMatrix<T> > localMatB;
    gsMatrix<T> localRhs_u; //, localRhs_p;

    // Local values of the right hand side
    gsMatrix<T> rhsVals;
    gsMapData<T> md;
};


} // namespace gismo

