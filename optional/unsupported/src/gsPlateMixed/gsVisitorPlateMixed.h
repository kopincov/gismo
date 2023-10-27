/** @file gsVisitorPlateMixed.h

    @brief Element visitor (volume integrals) for mixed formulation of Kirchhoff-Love plates

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/

#pragma once

#include <gsAssembler/gsGaussRule.h>
#include <gsPde/gsShellMixedPde.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorPlateMixed
{
public:

    /// Constructor
    gsVisitorPlateMixed(const gsPde<T> & pde, gsPlateMixedAssembler<T>& assembler): assembler(assembler)
    {
        m_pde_ptr= dynamic_cast<const gsShellMixedPde<T>* >(&pde);
        m_rhsf = static_cast<const gsShellMixedPde<T>*>(&pde)->force() ;
    }

    /// Function to initialize the assembly procedure.
    void initialize(const gsBasisRefs<T> & bases,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        bQ = assembler.system().colBasis(0);
        bPsi = assembler.system().colBasis(1);
        basesSize = bases.size();

        basisData.resize(basesSize);
        actives.resize(basesSize);
        physGrad.resize(basesSize);
        numActive.resize(basesSize);

        // Setup Quadrature (use phi basis, because higher degree)
        rule = gsGaussRule<T>(bases[bPsi], options.getReal("quA"), options.getInt("quB"));

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    /// Evaluate on element.
    inline void evaluate(const gsBasisRefs<T> & bases,
                         const gsGeometry<T>  & geo,
                         const gsMatrix<T>    & quNodes)
    {
        md.points = quNodes;

        for (size_t i = 0; i < basesSize; i++)
        {
            // Compute the active basis functions
            // Assumes actives are the same for all quadrature points on the elements
            bases[i].active_into(md.points.col(0), actives[i]);
            numActive[i] = actives[i].rows();

            // Evaluate basis functions on element
            bases[i].evalAllDers_into(md.points, 1, basisData[i]);
        }

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo.computeMap(md);

        // Evaluate right-hand side at the geometry points
        m_rhsf->eval_into(md.values[0], fVals);
        //std::cout<<"fVals"<<std::endl<<std::setprecision(20)<<fVals<<std::endl;

        // Initialize local matrix/rhs
        localMatA00.setZero(numActive[bQ], numActive[bQ]);
        localMatA11.setZero(2 * numActive[bPsi], 2 * numActive[bPsi]);
        localMatA10.setZero(2 * numActive[bPsi], numActive[bQ]);

        localMatB.setZero(numActive[bQ], numActive[bQ]);
        localMatC.setZero(numActive[bQ], numActive[bQ]);

        localRhsf.setZero(numActive[bQ], 1);

        // Initialize auxiliary matrices
        pId.setZero(3, numActive[bQ]);
        symCurlPhi.setZero(3, 2 * numActive[bPsi]);
    }

    inline void assemble(gsDomainIterator<T> & ,
                         const gsVector<T>   & quWeights)
    {
        std::vector<gsMatrix<T> > bVals(basesSize); // ursprünglichen version reference
        for(size_t i=0; i<basesSize; i++)
            bVals[i] = basisData[i][0];

        std::vector<gsMatrix<T> > bGrads(basesSize);
        for(size_t i=0; i<basesSize; i++)
            bGrads[i] = basisData[i][1];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute physical gradients at k as a Dim x NumActive matrix
            for(size_t i=0; i<basesSize; i++)
                transformGradients(md, k, bGrads[i], physGrad[i]);

            // Compute material matrix
            computeMaterialMatrix();

            // Compute mixed formulation quantities
            computeMixedFormulationQuantities(bVals, k);

            //m_C = gsMatrix<T, 3, 3>::Identity();
            //m_C(2,2) = 0.5;

            m_CModInv = m_C.inverse();

            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);

            // A00
            localMatA00.noalias() += weight * pId.transpose() * m_CModInv * pId;

            // A11
            localMatA11.noalias() += weight * symCurlPhi.transpose() * m_CModInv * symCurlPhi;

            // A10
            localMatA10.noalias() += weight * symCurlPhi.transpose() * m_CModInv * pId;

            // B
            localMatB.noalias() += -weight * physGrad[bQ].transpose() * physGrad[bQ];

            // Local rhs vector contribution
            localRhsf.noalias() += -weight *  bVals[bQ].col(k) * fVals.col(k);

            // C
            localMatC.noalias() += -weight * bVals[bQ].col(k) * bVals[bQ].col(k).transpose();

        }

/*

        gsMatrix<T, 2, 2> idMatrix = gsMatrix<T, 2, 2>::Identity();
        gsMatrix<T,2,2> Curl_i;
        gsMatrix<T,2,2> Curl_j;
        gsMatrix<T,2,2> symCurl_i;
        gsMatrix<T,2,2> symCurl_j;

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);

            // Compute physical gradients at k as a Dim x NumActive matrix
            for(size_t i=0; i<basesSize; i++)
                geoEval.transformGradients(k, bGrads[i], physGrad[i]);

            // A00
            for (index_t i = 0; i < numActive[bQ]; i++)
                for (index_t j = 0; j < numActive[bQ]; j++)
                    // Exploit symmetry
                    //for (index_t j = i; j < numActive; j++)
                {
                    localMatA00(i, j) +=  weight * (m_stCoefInv * bVals[bQ](j,k) * idMatrix).trace() * bVals[bQ](i,k);
                }


            // A11
            for (index_t i = 0; i < numActive[bPsi]; i++)
                for (index_t j = 0; j < numActive[bPsi]; j++)
                    // Exploit symmetry of A
                    //for (index_t j = i; j < numActive; j++)
                {
                    for(index_t m = 0; m<2; m++)
                        for(index_t l=0; l<2; l++)
                        {
                            Curl_i = gsMatrix<T, 2, 2>::Zero();
                            Curl_i(m,0)=physGrad[bPsi](1,i);
                            Curl_i(m,1)=-physGrad[bPsi](0,i);
                            Curl_j = gsMatrix<T, 2, 2>::Zero();
                            Curl_j(l,0)=physGrad[bPsi](1,j);
                            Curl_j(l,1)=-physGrad[bPsi](0,j);

                            symCurl_i = 0.5*(Curl_i+Curl_i.transpose());
                            symCurl_j = 0.5*(Curl_j+Curl_j.transpose());

                            localMatA11(m*numActive[bPsi]+i, l*numActive[bPsi]+j) += weight * (symCurl_i.cwiseProduct(m_stCoefInv*symCurl_j)).sum();
                        }
                }


            // A10
            for (index_t i = 0; i < numActive[bPsi]; i++)
                for (index_t j = 0; j < numActive[bQ]; j++)
                    for(index_t m = 0; m<2; m++)
                    {
                        Curl_i = gsMatrix<T, 2, 2>::Zero();
                        Curl_i(m,0)=physGrad[bPsi](1,i);
                        Curl_i(m,1)=-physGrad[bPsi](0,i);
                        symCurl_i = 0.5*(Curl_i+Curl_i.transpose());
                        localMatA10(m*numActive[bPsi]+i, j) += weight * (m_stCoefInv*symCurl_i).trace() * bVals[bQ](j,k);
                    }


            // B
            localMatB.noalias() += weight * (physGrad[bQ].transpose() * physGrad[bQ]);

            // Local rhs vector contribution
            localRhsf.noalias() += weight *  (bVals[bQ].col(k) * fVals.col(k));

            // C
            localMatC.noalias() += -weight * bVals[bQ].col(k) * bVals[bQ].col(k).transpose();

        }
        */

    }

    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {

        // Map patch-local DoFs to global DoFs
        actives_vec.resize(4);

        system.mapColIndices(actives[bQ], patchIndex, actives_vec[0], 0);
        system.mapColIndices(actives[bPsi], patchIndex, actives_vec[1], 1);
        system.mapColIndices(actives[bPsi], patchIndex, actives_vec[2], 2);
        system.mapColIndices(actives[bQ], patchIndex, actives_vec[3], 3);

        // build block information
        gsVector<index_t> p_vec(1);
        gsVector<index_t> phi_vec(2);
        gsVector<index_t> u_vec(1);
        p_vec << 0;
        phi_vec << 1,2;
        u_vec << 3;

        // A00
        system.pushToMatrix(localMatA00, actives_vec[0], eliminatedDofs[0], 0, 0);

        // A11
        system.pushToMatrix(localMatA11, actives_vec, eliminatedDofs, phi_vec, phi_vec);

        // A10, A01
        system.pushToMatrix(localMatA10, actives_vec, eliminatedDofs, phi_vec, p_vec);
        system.pushToMatrix(localMatA10.transpose(), actives_vec, eliminatedDofs, p_vec, phi_vec);

        // B and F, B^T
        system.push(localMatB, localRhsf, actives_vec, eliminatedDofs, u_vec, p_vec);
        system.pushToMatrix(localMatB.transpose(), actives_vec, eliminatedDofs, p_vec, u_vec);

        // C
        //system.pushToMatrix(localMatC, actives_w, eliminatedDofs[3], 3, 3);

    }

    void computeMaterialMatrix()
    {

        //const T E = m_pde_ptr->m_E;
        const T nu = m_pde_ptr->m_nu;

        //const T C_constant = E/(1-nu*nu);
        const T C_constant = 1;

        m_C(0,0) = 1;
        m_C(1,1) = 1;
        m_C(2,2) = 0.5*(1-nu);
        m_C(1,0) = m_C(0,1) = nu;
        m_C(2,0) = m_C(0,2) = 0;
        m_C(2,1) = m_C(1,2) = 0;

        m_C *=C_constant;


    }

    void computeMixedFormulationQuantities(const std::vector<gsMatrix<T> > & bVals, const index_t k)
    {
        pId.setZero(3, numActive[bQ]);
        symCurlPhi.setZero(3, 2*numActive[bPsi]);

        index_t j;

        for (index_t i = 0; i!= numActive[bQ]; ++i) // basis function
        {
            // pId
            pId(0,i) = bVals[bQ](i,k);
            pId(1,i) = bVals[bQ](i,k);
        }


        for (index_t i = 0; i!= numActive[bPsi]; ++i) // basis function
        {
            // symCurl
            j = 0;
            symCurlPhi(0,i+j*numActive[bPsi]) = physGrad[bPsi](1,i);

            j= 1;
            symCurlPhi(1,i+j*numActive[bPsi]) = -physGrad[bPsi](0,i);

            j = 0;
            symCurlPhi(2,i+j*numActive[bPsi]) = -0.5*physGrad[bPsi](0,i);
            j = 1;
            symCurlPhi(2,i+j*numActive[bPsi]) = 0.5*physGrad[bPsi](1,i);

        }

        // Voigt notation 1*M_12

    }


    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    protected:
    // assembler
    gsPlateMixedAssembler<T>& assembler;
    size_t bQ;
    size_t bPsi;
    size_t basesSize;

    // Pde
    const gsShellMixedPde<T>* m_pde_ptr;

    // rhs function
    const gsFunction<T> * m_rhsf;
    // Material matrix
    gsMatrix<T,3,3> m_C;
    gsMatrix<T,3,3> m_CModInv;


    // Local values rhs function
    gsMatrix<T> fVals;


protected:
    // Basis values
    std::vector<std::vector<gsMatrix<T> > > basisData;
    std::vector<gsMatrix<index_t> > actives;
    std::vector<gsMatrix<T> >		 physGrad;
    std::vector<index_t >            numActive;
    std::vector<gsMatrix<index_t> > actives_vec;


protected:

    // Mixed formulation quantities
    gsMatrix<T,3> symCurlPhi;
    gsMatrix<T,3> pId;
    gsMatrix<T> bGrads_k;

    // Local matrices
    gsMatrix<T> localMatA00, localMatA10, localMatA11, localMatB, localMatC;
    gsMatrix<T> localRhsf;

    gsMapData<T> md;
};


} // namespace gismo

