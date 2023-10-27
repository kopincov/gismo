/** @file gsKKTPervertedPoissonBoundary.h

    @brief KKT Perverted Poisson equation element visitor.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once


namespace gismo
{

/** @brief Visitor for the KKT formulated perverted Poisson optimal control problem.
   ( aM   0  M  )  and (0  )
   (  0   Md A^T)      (Mdd)
   (  M   A  0  )      (0  )
   This visitor creats Md and Mdd.
*/
template <class T>
class gsVisitorKKTPervertedPoissonBoundary
{
public:

    /// Constructor with the right hand side function of the poisson equation
    gsVisitorKKTPervertedPoissonBoundary(std::vector<gsFunction<T> * > const & desired, boxSide s) :
    m_desired(desired), side(s)
    {

    }

    // initialized for each patch
    void initialize(const gsBasisRefs<T> & basis,
                    gsQuadRule<T>        & rule,
                    unsigned             & evFlags )
    {
        if ( (int) m_desired.size() != (int) 2*basis[1].dim())
            GISMO_ERROR("Assumes one function for each side");

        const int dir = side.direction();
        gsVector<index_t> numQuadNodes( basis[1].dim() );
        for (int i = 0; i < basis[1].dim(); ++i) // to do: improve
            numQuadNodes[i] = basis[1].maxDegree() + 1; // take quadrature from highest degree
        numQuadNodes[dir] = 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_MEASURE|NEED_GRAD_TRANSFORM;

    }

    // Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        basisRefs[1].active_into(quNodes.col(0), actives_u);
        const index_t numActU = actives_u.rows();
        //basisRefs[1].eval_into(quNodes,  basisData_u);
        basisRefs[1].evalAllDers_into(quNodes, 1, basisData_u);
        
        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate the desired data
        m_desired[side-1]->eval_into(geoEval.values(), rhsData);

        localMatMb .setZero(numActU,numActU);
        localRhs_ub.setZero(numActU, 1);
    }
    
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        //const typename gsMatrix<T>::Block basisVals_u  = basisData_u.topRows(numActU);
        gsMatrix<T> & basisGrads_u = basisData_u[1];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, side, unormal);
            geoEval.transformGradients(k, basisGrads_u, physGrad_u);

            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] *unormal.norm();

            // Compute the unit normal vector
            unormal.normalize();

            //gsMatrix<T> test123 = unormal;//physGrad_u;//*unormal;
            //gsInfo << test123.rows() <<" X " << test123.cols() << std::endl;
            // Local block A
            localMatMb.noalias() += weight * ((physGrad_u.transpose()*unormal)*(unormal.transpose()*physGrad_u)) ;
            //localMatMb.noalias() += weight * (basisVals_u.col(k) * basisVals_u.col(k).transpose());
            //localRhs_ub.noalias()   += weight * ( rhsData.col(k)*(unormal.transpose()*physGrad_u) );
            localRhs_ub.noalias()   += weight * (( physGrad_u.transpose() * unormal)* rhsData.col(k).transpose());

            //localRhs.noalias() += weight *(( physBasisGrad.transpose() * unormal )* neuData.col(k).transpose());

        }
    }
    
    inline void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                              const gsMatrix<T>     & eliminatedDofs,
                              const index_t           patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
        GISMO_ASSERT(mappers.size() ==3,
                     "Expecting three dof mappers, f,u and w.");

        const gsDofMapper & uMap = mappers[1];

        const index_t fsz = mappers[0].freeSize();

        const index_t numActU = actives_u.rows();

        // Local Dofs to global dofs
        uMap.localToGlobal(actives_u, patchIndex, actives_u);


        for (index_t iu=0; iu < numActU; ++iu)
        {
            const int iiu = actives_u(iu);
            if ( uMap.is_free_index(iiu) )
            {
                rhsMatrix.row(iiu+fsz) += localRhs_ub.row(iu);

                for (index_t ju=0; ju < numActU; ++ju) // Build Mb
                {
                    const int jju = actives_u(ju);
                    if ( uMap.is_free_index(jju) )
                    {
                        sysMatrix.coeffRef(iiu+fsz, jju+fsz) += localMatMb(iu, ju);
                    }
                    else
                    {
                        rhsMatrix.coeffRef(iiu+fsz,0) -= localMatMb(iu, ju)*eliminatedDofs( uMap.global_to_bindex(jju),0 );
                    }
                }
            }
        }
    }

protected:

    // Right hand side
    std::vector<gsFunction<T> * > m_desired;

    // Side
    boxSide side;

protected:
    // Basis values
    std::vector<gsMatrix<T> > basisData_u;
    gsMatrix<T> physGrad_u;
    gsMatrix<index_t> actives_u;

protected:
    // Local matrices
    gsMatrix<T> localMatMb;
    gsMatrix<T> localRhs_ub;

    gsMatrix<T> rhsVals;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> rhsData;
};


} // namespace gismo

