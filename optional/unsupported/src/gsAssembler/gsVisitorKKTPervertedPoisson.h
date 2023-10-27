/** @file gsKKTPervertedPoisson.h

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
*/
template <class T>
class gsVisitorKKTPervertedPoisson
{
public:

    /// Constructor with the right hand side function of the poisson equation
    gsVisitorKKTPervertedPoisson(const T alpha) :
    m_alpha(alpha)
    { }

    // initialized for each patch
    void initialize(const gsBasisRefs<T> & basis,
                    gsQuadRule<T>        & rule,
                    unsigned             & evFlags )
    {
        gsVector<index_t> numQuadNodes( basis[1].dim() );
        for (int i = 0; i < basis[1].dim(); ++i) // to do: improve
            numQuadNodes[i] = basis[1].maxDegree() + 1; // take quadrature from highest degree
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_MEASURE|NEED_GRAD_TRANSFORM|NEED_2ND_DER;

    }

    // Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        basisRefs[0].active_into(quNodes.col(0), actives_f);
        basisRefs[1].active_into(quNodes.col(0), actives_u);
        basisRefs[2].active_into(quNodes.col(0), actives_w);
        const index_t numActF = actives_f.rows();
        const index_t numActU = actives_u.rows();
        const index_t numActW = actives_w.rows();
        basisRefs[0].eval_into       (quNodes,    basisVals_f);
        basisRefs[1].evalAllDers_into(quNodes, 2, basisData_u);
        basisRefs[2].eval_into       (quNodes,    basisVals_w);
        
        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate the Dirichlet data
        //rhs_ptr->eval_into(geoEval.values(), rhsData);
        //dis_ptr->eval_into(geoEval.values(), disData);


        GISMO_ASSERT(numActF == numActW, "Assume that some absis is used for f and w");

        localMatA .setZero(numActW,numActU);
        localMatM .setZero(numActW,numActF);
        localMatMu.setZero(numActU,numActU);
        localMataM.setZero(numActF,numActF);
        localRhs_f.setZero(numActF, 1);
        localRhs_d.setZero(numActU, 1);
        localRhs_w.setZero(numActW, 1);
    }
    
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        //const index_t numActF = actives_u.rows();
        //const index_t numActW = actives_w.rows();

        gsMatrix<T> & basisVals_u  = basisData_u[0];
        gsMatrix<T> & basisGrads_u = basisData_u[1];
        gsMatrix<T> & basis2ndDerivs_u = basisData_u[2];
        //basisVals_f basisVals_w

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            //geoEval.transformGradients(k, basisGrads_u, physGrad_u);
            // Compute physical laplacian at k as a 1 x numActive matrix
            geoEval.transformLaplaceHgrad (k, basisGrads_u, basis2ndDerivs_u, physBasisLaplace_u);


            // Local block A
            //gsInfo << " dim basisVals_u " << basisVals_u.rows() << " X " << basisVals_u.cols() << std::endl;
            //gsInfo << " dim physBasisLaplace_u " << physBasisLaplace_u.rows() << " X " << physBasisLaplace_u.cols() << std::endl;
            //gsInfo << " dim basisVals_w " << basisVals_w.rows() << " X " << basisVals_w.cols() << std::endl;
            //gsMatrix<T> tmptest = physBasisLaplace_u.transpose() * basisVals_w.col(k).transpose();
            //gsInfo << tmptest << std::endl;
            //localMatA.noalias()  += weight * ((basisVals_u.col(k) - physBasisLaplace_u.transpose()) * basisVals_w.col(k).transpose()).transpose();
            localMatA.noalias()  += -weight * (physBasisLaplace_u.transpose() * basisVals_w.col(k).transpose()).transpose();

            //gsInfo << "test12" << std::endl;

            localMataM.noalias() += m_alpha * weight * (basisVals_f.col(k) * basisVals_f.col(k).transpose());
            localMatMu.noalias() += weight * (basisVals_u.col(k) * basisVals_u.col(k).transpose());
            //gsInfo << "test13" << std::endl;
            localMatM.noalias()  += weight *(basisVals_f.col(k) * basisVals_w.col(k).transpose());
            //gsInfo << "test15" << std::endl;
            //localRhs_d.noalias() += weight * (basisVals_u.col(k) * disData.col(k).transpose() );
            //localRhs_f.noalias() += m_alpha * weight * (basisVals_f.col(k) * rhsData.col(k).transpose() );


        }
        //gsInfo << "alpha: " << m_alpha << std::endl;
    }
    
    inline void localToGlobal(const gsStdVectorRef<gsDofMapper> & mappers,
                              const gsMatrix<T>     & eliminatedDofs,
                              const index_t patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
        GISMO_ASSERT(mappers.size() ==3,
                     "Expecting three dof mappers, f,u and w.");

        const gsDofMapper & fMap = mappers[0];
        const gsDofMapper & uMap = mappers[1];
        const gsDofMapper & wMap = mappers[2];

        const index_t fsz = fMap.freeSize();
        const index_t usz = uMap.freeSize();
        //const index_t wsz = wMap.freeSize();


        // Local Dofs to global dofs
        fMap.localToGlobal(actives_f, patchIndex, actives_f);
        uMap.localToGlobal(actives_u, patchIndex, actives_u);
        wMap.localToGlobal(actives_w, patchIndex, actives_w);

        const index_t numActF = actives_f.rows();
        const index_t numActU = actives_u.rows();
        const index_t numActW = actives_w.rows();

        for (index_t i_f=0; i_f < numActF; ++i_f)
        {
            const int iif = actives_f(i_f);
            if (! fMap.is_free_index(iif) )
            {
                GISMO_ERROR("Something is wrong!");
            }

            //rhsMatrix.row(iif) += localRhs_f.row(i_f);

            for (index_t jf=0; jf < numActF; ++jf) // Build \alpha M
            {
                const int jjf = actives_f(jf);
                sysMatrix.coeffRef(iif, jjf) += localMataM(i_f,jf);
            }
            //gsInfo << " test L2G 03 " << std::endl;

            for (index_t jw=0; jw < numActW; ++jw) // Build M
            {
                const int jjw = actives_w(jw);
                sysMatrix.coeffRef(iif, jjw+fsz+usz) += localMatM(i_f, jw);

            }
            //gsInfo << " test L2G 04 " << std::endl;

        }
        //gsInfo << " test L2G 100 " << std::endl;

        for (index_t iu=0; iu < numActU; ++iu)
        {
            const int iiu = actives_u(iu);

            if ( uMap.is_free_index(iiu) )
            {

                for (index_t jw=0; jw < numActW; ++jw) // Build A^T
                {
                    const int jjw = actives_w(jw);
                    sysMatrix.coeffRef(iiu+fsz, jjw+fsz+usz) += localMatA(jw, iu); //A^T
                }
            }
        }

        for (index_t iw=0; iw < numActW; ++iw)
        {
            const int iiw = actives_w(iw);
            for (index_t jf=0; jf < numActF; ++jf) // Build M
            {
                const int jjf = actives_f(jf);
                sysMatrix.coeffRef(iiw+fsz+usz, jjf) += localMatM(iw,jf);
            }
            for (index_t ju=0; ju < numActU; ++ju) // Build A
            {
                //gsInfo << ju << " of " << numActW << std::endl;
                const int jju = actives_u(ju);
                //gsInfo << iif << " -- " << std::endl;
                //gsInfo << " test L2G 02 " << std::endl;

                if ( uMap.is_free_index(jju) )
                {
                    sysMatrix.coeffRef(iiw+fsz+usz, jju+fsz) += localMatA(iw, ju);
                }
                else // Umap.is_boundary_index(jj)
                {
                    const int bb = uMap.global_to_bindex(jju);
                    //gsInfo << "eliminated DoF: " <<  eliminatedDofs(bb,0)<< std::endl;
                    if (math::abs(localMatA(iw, ju)*eliminatedDofs( uMap.global_to_bindex(jju),0 ))> 1e-12 )
                    {
                        //gsInfo <<" In Main visitor, added to the rhs:"
                        //       << localMatA(iw, ju)*eliminatedDofs( uMap.global_to_bindex(jju),0 ) << std::endl;
                    }
                    rhsMatrix.coeffRef(iiw+fsz+usz, 0) -= localMatA(iw, ju)*eliminatedDofs(bb,0);
                }
            }
        }
    }

protected:

    // Right hand side
    T m_alpha;

protected:
    // Basis values
    gsMatrix<T> basisVals_f, basisVals_w;
    std::vector<gsMatrix<T> > basisData_u;
    gsMatrix<T> physGrad_u, physBasisLaplace_u;
    gsMatrix<index_t> actives_f, actives_u, actives_w;

protected:
    // Local matrices
    gsMatrix<T> localMatA;
    gsMatrix<T> localMatM;
    gsMatrix<T> localMatMu;
    gsMatrix<T> localMataM;
    gsMatrix<T> localRhs_f; //Should be empty
    gsMatrix<T> localRhs_d;
    gsMatrix<T> localRhs_w; //Should be empty
    //gsMatrix<T> rhsData;
    //gsMatrix<T> disData;


    // Local values of the right hand side //Do this in a boundary visitor
    // gsMatrix<T> rhsVals;
};


} // namespace gismo

