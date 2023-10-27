/** @file gsVisitorShellMixed.h

    @brief Element visitor (volume integrals) for mixed formulation of Kirchhoff-Love shells

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/

#pragma once

#include <gsAssembler/gsGaussRule.h>
#include <gsPde/gsShellMixedPde.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsCore/gsRationalBasis.h>
#include <gsCore/gsBasis.h>
#include <gsPlateMixed/gsShellMixedAssembler.h>

namespace gismo
{

template <class T>
class gsVisitorShellMixed
{
public:

    /// Constructor
    gsVisitorShellMixed(const gsPde<T> & pde, const gsShellMixedAssembler<T>& assembler):m_assembler(assembler)
    {
        m_pde_ptr = dynamic_cast<const gsShellMixedPde<T> * >(&pde);
        m_force = m_pde_ptr->force();
    }

    /// Function to initialize the assembly procedure.
    void initialize(const gsBasisRefs<T> & bases,
                    const index_t          patchIndex,
                    const gsOptionList   & options,
                          gsQuadRule<T>  & rule)
    {
        // Set options
        Mmixed = options.getSwitch("Mmixed");
        Nmixed = options.getSwitch("Nmixed");

        // Store basis information
        pBasis = m_assembler.system().colBasis(0);
        phiBasis = m_assembler.system().colBasis(1);

        N11Basis = m_assembler.system().colBasis(3);
        N22Basis = m_assembler.system().colBasis(4);
        N12Basis = m_assembler.system().colBasis(5);

        uBasis = m_assembler.system().colBasis(6);

        basesSize = bases.size();

        basisData.resize(basesSize);
        actives.resize(basesSize);
        numActive.resize(basesSize);
        bGrads_k.resize(basesSize);


        // Setup Quadrature
        rule = gsGaussRule<T>(bases[0], options.getReal("quA"), options.getInt("quB"));

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_JACOBIAN | NEED_MEASURE | NEED_2ND_DER;
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
            if (!Mmixed)
                bases[i].evalAllDers_into(md.points, 2, basisData[i]);
            else
                bases[i].evalAllDers_into(md.points, 1, basisData[i]);
        }
        numActiveN = numActive[N11Basis] + numActive[N22Basis] + numActive[N12Basis];

        // Compute geometry related quantities
        geo.computeMap(md);

        // Compute 3rd derivatives of geometry
        gComputeDeriv3Old(geo);
        //gComputeDeriv3(geoEval, md.points); // uses implementation of 3rd derivative for NURBS geometry

        // Evaluate right-hand side at the geometry points
        m_force->eval_into(md.values[0], forceVals);

        // Initialize local matrix/rhs
        localMatA00.setZero(numActive[pBasis], numActive[pBasis]);
        localMatA11.setZero(2 * numActive[phiBasis], 2 * numActive[phiBasis]);
        localMatA10.setZero(2 * numActive[phiBasis], numActive[pBasis]);

        localMatB0.setZero(3 * numActive[uBasis], numActive[pBasis]);
        localMatB1.setZero(3 * numActive[uBasis], 2 * numActive[phiBasis]);

        localMatC.setZero(3 * numActive[uBasis], 3 * numActive[uBasis]);
        localMatC00.setZero(numActiveN, numActiveN);
        localMatC10.setZero(3 * numActive[uBasis], numActiveN);


        localRhs.setZero(3 * numActive[uBasis], 1);

        // Initialize auxiliary matrices
        pId.setZero(3, numActive[pBasis]);
        symCurlPhi.setZero(3, 2 * numActive[phiBasis]);

        N.setZero(3, numActiveN);

        mStrain.setZero(3, 3 * numActive[uBasis]);
        bStrainFirstOrder.setZero(3, 3 * numActive[uBasis]);
        bStrain.setZero(3, 3 * numActive[uBasis]);

        // Initialize auxiliary maps (derivatives to column index)
        map_gDerSecond[0][0] = 0;
        map_gDerSecond[1][1] = 1;
        map_gDerSecond[0][1] = 2;
        map_gDerSecond[1][0] = 2;

        map_gDerThird[0][0][0] = 0;
        map_gDerThird[1][1][1] = 1;
        map_gDerThird[0][0][1] = 2;
        map_gDerThird[0][1][0] = 2;
        map_gDerThird[1][0][0] = 2;
        map_gDerThird[0][1][1] = 3;
        map_gDerThird[1][0][1] = 3;
        map_gDerThird[1][1][0] = 3;

        map_contraBDer[0][0] = 0;
        map_contraBDer[1][1] = 1;
        map_contraBDer[0][1] = 2;
        map_contraBDer[1][0] = 3;

        // Map matrix entry (alpha, beta) to component in Voigt notation
        VoigtVtoM.col(0) << 0, 1, 0;
        VoigtVtoM.col(1) << 0, 1, 1;
    }

    inline void assemble(      gsDomainIterator<T> & element,
                         const gsVector<T>         & quWeights)
    {
        std::vector<gsMatrix<T> > bVals(basesSize);
        std::vector<gsMatrix<T> > bGrads(basesSize);
        std::vector<gsMatrix<T> > bGrads2(basesSize);

        for (size_t i = 0; i < basesSize; i++)
        {
            bVals[i] = basisData[i][0];
            bGrads[i] = basisData[i][1];

            if (!Mmixed)
                bGrads2[i] = basisData[i][2];
        }


        const T thickness = m_pde_ptr->m_thickness;

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute geometry quantities
            computeGeometryQuantities(k);

            // Compute material matrix
            computeMaterialMatrix();

            // Compute mixed formulation quantities
            computeMixedFormulationQuantities(bVals, bGrads, k);
            if (Nmixed)
                computeNmixedQuantities(bVals, bGrads, k);


            for (size_t i = 0; i < basesSize; i++)
            {
                bGrads_k[i] = bGrads[i].col(k);
                bGrads_k[i].resize(2, numActive[i]);
            }

            // Compute shell quantities
            computeShellQuantities(bVals[uBasis], bGrads[uBasis], bGrads2[uBasis], k);

            //m_C = gsMatrix<T, 3, 3>::Identity();
            m_CModInv = (thickness * thickness * thickness / 12.0 * md.measure(k) * m_C).inverse();
            //m_CModInv = (1/12.0 * md.measure(k) * m_C).inverse();
            m_CModInvN = (thickness * md.measure(k) * m_C).inverse();
            //m_CModInvN = (md.measure(k) * m_C).inverse();

            // compute contravariant components of forceVals in the covariant basis
            forceVals.col(k) = F.inverse() * forceVals.col(k);

            // Store weight
            const T weight = quWeights[k];

            if (Mmixed)
            {
                // A00
                localMatA00.noalias() += weight * pId.transpose() * m_CModInv * pId;

                // A11
                localMatA11.noalias() += weight * symCurlPhi.transpose() * m_CModInv * symCurlPhi;

                // A10
                localMatA10.noalias() += weight * symCurlPhi.transpose() * m_CModInv * pId;

                // B0
                localMatB0.bottomRows(numActive[uBasis]).noalias() +=
                    weight * bGrads_k[uBasis].transpose() * bGrads_k[pBasis];
                localMatB0.noalias() += -weight * bStrainFirstOrder.transpose() * pId;

                // B1
                localMatB1.noalias() += -weight * bStrainFirstOrder.transpose() * symCurlPhi;
            }
            else
            {
                localMatC.noalias() +=
                    -weight * md.measure(k) * thickness * thickness * thickness / 12.0 * bStrain.transpose() * m_C
                        * bStrain;
            }


            if (Nmixed)
            {
                // C00
                localMatC00.noalias() += weight * N.transpose() * m_CModInvN * N;

                // C10
                localMatC10.noalias() += -weight * mStrain.transpose() * N;

            }
            else
            {
                // C
                localMatC.noalias() += -weight * md.measure(k) * thickness * mStrain.transpose() * m_C * mStrain;
                //localMatC.noalias() += -weight * md.measure(k) * mStrain.transpose()* m_C * mStrain;
            }


            // Local rhs vector contribution
            for (index_t j = 0; j != 3; ++j)
                localRhs.middleRows(j * numActive[uBasis], numActive[uBasis]).noalias() +=
                    -weight * md.measure(k) * forceVals(j, k) * bVals[uBasis].col(k);

        } // end loop quad points

        //gsInfo<<"localMatA00 \n"<<localMatA00<<"\n";
        //gsInfo<<"localMatA11 \n"<<localMatA11<<"\n";
        //gsInfo<<"localMatA10 \n"<<localMatA10<<"\n";
        //gsInfo<<"localMatB0 \n"<<localMatB0<<"\n";
        //gsInfo<<"localMatB1 \n"<<localMatB1<<"\n";
        //gsInfo<<"localMatC \n"<<localMatC<<"\n";
        //gsInfo<<"localRhs \n"<<localRhs<<"\n";

        //gsInfo<<"mStrain \n"<<mStrain<<"\n";
    }

    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        //const T thickness = m_pde_ptr->m_thickness;

        // Map patch-local DoFs to global DoFs
        actives_vec.resize(9);
        if(Mmixed)
        {
            system.mapColIndices(actives[pBasis], patchIndex, actives_vec[0], 0);
            system.mapColIndices(actives[phiBasis], patchIndex, actives_vec[1], 1);
            system.mapColIndices(actives[phiBasis], patchIndex, actives_vec[2], 2);
        }

        if(Nmixed)
        {
            system.mapColIndices(actives[N11Basis], patchIndex, actives_vec[3], 3);
            system.mapColIndices(actives[N22Basis], patchIndex, actives_vec[4], 4);
            system.mapColIndices(actives[N12Basis], patchIndex, actives_vec[5], 5);
        }

        system.mapColIndices(actives[uBasis], patchIndex, actives_vec[6], 6);
        system.mapColIndices(actives[uBasis], patchIndex, actives_vec[7], 7);
        system.mapColIndices(actives[uBasis], patchIndex, actives_vec[8], 8);

        // build block information
        gsVector<index_t> p_vec(1);
        gsVector<index_t> phi_vec(2);
        gsVector<index_t> N_vec(3);
        gsVector<index_t> u_vec(3);
        p_vec << 0;
        phi_vec << 1,2;
        N_vec << 3,4,5;
        u_vec << 6,7,8;


        if(Mmixed)
        {
            // A00
            system.pushToMatrix(localMatA00, actives_vec[0], eliminatedDofs[0], 0, 0);

            // A11
            system.pushToMatrix(localMatA11, actives_vec, eliminatedDofs, phi_vec, phi_vec);

            // A10, A01
            system.pushToMatrix(localMatA10, actives_vec, eliminatedDofs, phi_vec, p_vec);
            system.pushToMatrix(localMatA10.transpose(), actives_vec, eliminatedDofs, p_vec, phi_vec);

            // B0, B0^T
            system.pushToMatrix(localMatB0, actives_vec, eliminatedDofs, u_vec, p_vec);
            system.pushToMatrix(localMatB0.transpose(), actives_vec, eliminatedDofs, p_vec, u_vec);
            //system.pushToMatrix(thickness*thickness*localMatB0, actives_vec, eliminatedDofs, u_vec, p_vec);
            //system.pushToMatrix(thickness*thickness*localMatB0.transpose(), actives_vec, eliminatedDofs, p_vec, u_vec);

            // B1, B1^T
            system.pushToMatrix(localMatB1, actives_vec, eliminatedDofs, u_vec, phi_vec);
            system.pushToMatrix(localMatB1.transpose(), actives_vec, eliminatedDofs, phi_vec, u_vec);
            //system.pushToMatrix(thickness*thickness*localMatB1, actives_vec, eliminatedDofs, u_vec, phi_vec);
            //system.pushToMatrix(thickness*thickness*localMatB1.transpose(), actives_vec, eliminatedDofs, phi_vec, u_vec);
        }

        if(Nmixed)
        {
            // C00
            system.pushToMatrix(localMatC00, actives_vec, eliminatedDofs, N_vec, N_vec);

            // C10, C01
            system.pushToMatrix(localMatC10, actives_vec, eliminatedDofs, u_vec, N_vec);
            system.pushToMatrix(localMatC10.transpose(), actives_vec, eliminatedDofs, N_vec, u_vec);
            //system.pushToMatrix(thickness * localMatC10, actives_vec, eliminatedDofs, u_vec, N_vec);
            //system.pushToMatrix(localMatC10.transpose(), actives_vec, eliminatedDofs, N_vec, u_vec);
        }

        // C and F
        system.push(localMatC, localRhs, actives_vec, eliminatedDofs, u_vec, u_vec);

    }

    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 protected:
    void computeGeometryQuantities(const index_t k)
    {
        // first derivative
        normal(md, k,normalVec);
        normalVec.normalize();
        F.leftCols(2) = md.jacobian(k);
        F.col(2)      = normalVec;
        FInvT = F.inverse().transpose();

        // first ff
        a_co = F.leftCols(2).transpose()*F.leftCols(2);
        a_contra = FInvT.leftCols(2).transpose()*FInvT.leftCols(2);

        // second derivative
        gsMatrix<T> temp = md.deriv2(k);
        temp.resize(3,3);
        gDerSecond = temp.transpose();

        // second ff
        b_co_V = gDerSecond.transpose()*normalVec;
        b_co(0,0)= b_co_V(0,0);
        b_co(1,1)= b_co_V(1,0);
        b_co(1,0) = b_co(0,1) = b_co_V(2,0);

        b_mixed = a_contra * b_co;

        // third ff
        c_co = b_mixed.transpose()*b_co;

        // derivative normal
        gsVector<T> m_normalDer= gDerSecond.col(0).cross(F.col(1)) + F.col(0).cross(gDerSecond.col(2));
        m_normalDer /= md.measure(k);
        normalDer.col(0)= m_normalDer - m_normalDer.dot(normalVec) * normalVec;

        m_normalDer= gDerSecond.col(2).cross(F.col(1)) + F.col(0).cross(gDerSecond.col(1));
        m_normalDer /= md.measure(k);
        normalDer.col(1)= m_normalDer - m_normalDer.dot(normalVec) * normalVec;

        // derivative contravariant basis
        gsVector<T,3> contraBCoefs;
        contraBCoefs(0)= -FInvT.col(0).dot(gDerSecond.col(0));
        contraBCoefs(1)= -FInvT.col(0).dot(gDerSecond.col(2));
        contraBCoefs(2)= -FInvT.col(0).dot(normalDer.col(0));
        contraBDer.col(0) = FInvT * contraBCoefs;

        contraBCoefs(0)= -FInvT.col(1).dot(gDerSecond.col(2));
        contraBCoefs(1)= -FInvT.col(1).dot(gDerSecond.col(1));
        contraBCoefs(2)= -FInvT.col(1).dot(normalDer.col(1));
        contraBDer.col(1) = FInvT * contraBCoefs;

        contraBCoefs(0)= -FInvT.col(1).dot(gDerSecond.col(0));
        contraBCoefs(1)= -FInvT.col(1).dot(gDerSecond.col(2));
        contraBCoefs(2)= -FInvT.col(1).dot(normalDer.col(0));
        contraBDer.col(2) = FInvT * contraBCoefs;

        contraBCoefs(0)= -FInvT.col(0).dot(gDerSecond.col(2));
        contraBCoefs(1)= -FInvT.col(0).dot(gDerSecond.col(1));
        contraBCoefs(2)= -FInvT.col(0).dot(normalDer.col(1));
        contraBDer.col(3) = FInvT * contraBCoefs;

        // third derivative
        gsMatrix<T> temp_deriv3 = gDerThird_quNodes.col(k);
        temp_deriv3.resize(4,3);
        gDerThird = temp_deriv3.transpose();
    }

    void computeMaterialMatrix()
    {

        /*
        const T lambda = m_pde_ptr->m_lambda;
        const T mu = m_pde_ptr->m_mu;

        const T C_constant = 4*lambda*mu/(lambda+2*mu);

        m_C(0,0) = C_constant*a_contra(0,0)*a_contra(0,0) + 2*mu*(2*a_contra(0,0)*a_contra(0,0));
        m_C(1,1) = C_constant*a_contra(1,1)*a_contra(1,1) + 2*mu*(2*a_contra(1,1)*a_contra(1,1));
        m_C(2,2) = C_constant*a_contra(0,1)*a_contra(0,1) + 2*mu*(a_contra(0,0)*a_contra(1,1) + a_contra(0,1)*a_contra(0,1));
        m_C(1,0) = m_C(0,1) = C_constant*a_contra(0,0)*a_contra(1,1) + 2*mu*(2*a_contra(0,1)*a_contra(0,1));
        m_C(2,0) = m_C(0,2) = C_constant*a_contra(0,0)*a_contra(0,1) + 2*mu*(2*a_contra(0,0)*a_contra(0,1));
        m_C(2,1) = m_C(1,2) = C_constant*a_contra(0,1)*a_contra(1,1) + 2*mu*(2*a_contra(0,1)*a_contra(1,1));
*/

        const T E = m_pde_ptr -> m_E;
        const T nu = m_pde_ptr -> m_nu;

        const T C_constant = E/(1-nu*nu);

        m_C(0,0) = (1-nu)/2*(2*a_contra(0,0)*a_contra(0,0)) + nu*(a_contra(0,0)*a_contra(0,0));
        m_C(1,1) = (1-nu)/2*(2*a_contra(1,1)*a_contra(1,1)) + nu*(a_contra(1,1)*a_contra(1,1));
        m_C(2,2) = (1-nu)/2*(a_contra(0,0)*a_contra(1,1) + a_contra(0,1)*a_contra(0,1)) + nu*(a_contra(0,1)*a_contra(0,1));
        m_C(1,0) = m_C(0,1) = (1-nu)/2*(2*a_contra(0,1)*a_contra(0,1)) + nu*(a_contra(0,0)*a_contra(1,1));
        m_C(2,0) = m_C(0,2) = (1-nu)/2*(2*a_contra(0,0)*a_contra(0,1)) + nu*(a_contra(0,0)*a_contra(0,1));
        m_C(2,1) = m_C(1,2) = (1-nu)/2*(2*a_contra(0,1)*a_contra(1,1)) + nu*(a_contra(0,1)*a_contra(1,1));

        m_C*= C_constant;



    }

    void computeMixedFormulationQuantities(const std::vector<gsMatrix<T> > & bVals, const std::vector<gsMatrix<T> > & bGrads, const index_t k)
    {
        pId.setZero(3, numActive[pBasis]);
        symCurlPhi.setZero(3, 2*numActive[phiBasis]);

        for (index_t i = 0; i!= numActive[pBasis]; ++i) // basis function
        {
            // pId
            pId(0,i) = bVals[pBasis](i,k);
            pId(1,i) = bVals[pBasis](i,k);
        }

        index_t j;
        for (index_t i = 0; i!= numActive[phiBasis]; ++i) // basis function
        {
            // symCurl
            j = 0;
            symCurlPhi(0,i+j*numActive[phiBasis]) = bGrads[phiBasis](2*i+1,k);

            j= 1;
            symCurlPhi(1,i+j*numActive[phiBasis]) = -bGrads[phiBasis](2*i,k);

            j = 0;
            symCurlPhi(2,i+j*numActive[phiBasis]) = -0.5*bGrads[phiBasis](2*i,k);
            j = 1;
            symCurlPhi(2,i+j*numActive[phiBasis]) = 0.5*bGrads[phiBasis](2*i+1,k);
        }

        // Voigt notation 1*M_12

    }

    void computeNmixedQuantities(const std::vector<gsMatrix<T> > & bVals, const std::vector<gsMatrix<T> > & bGrads, const index_t k)
    {
        N.setZero(3, numActiveN);

        // N11
        for (index_t i = 0; i!= numActive[N11Basis]; ++i) // basis function N11
            N(0,i) = bVals[N11Basis](i,k);
        // N22
        for (index_t i = 0; i!= numActive[N22Basis]; ++i) // basis function N22
            N(1,i + numActive[N11Basis]) = bVals[N22Basis](i,k);
        // N12
        for (index_t i = 0; i!= numActive[N12Basis]; ++i) // basis function N12
            N(2,i + numActive[N11Basis] + numActive[N22Basis]) = bVals[N12Basis](i,k);


        // Voigt notation 1*N_12

    }

    void computeShellQuantities(const gsMatrix<T> & bVals, const gsMatrix<T> & bGrads, const gsMatrix<T> & bGrads2, const index_t k)
    {
        mStrain.setZero(3, 3*numActive[uBasis]);
        bStrainFirstOrder.setZero(3, 3*numActive[uBasis]);
        bStrain.setZero(3, 3*numActive[uBasis]);

        gsMatrix<T,3,1> temp;
        gsMatrix<T,3,1> temp1;
        index_t j;
        for (index_t i = 0; i!= numActive[uBasis]; ++i) // basis function
        {
            // membrane strain
            // epsilon
            j = 0;
            mStrain(0,i+j*numActive[uBasis]) += bGrads(2*i,k);

            j = 1;
            mStrain(1,i+j*numActive[uBasis]) += bGrads(2*i+1,k);

            j = 0;
            mStrain(2,i+j*numActive[uBasis]) += 0.5*bGrads(2*i+1,k);
            j = 1;
            mStrain(2,i+j*numActive[uBasis]) += 0.5*bGrads(2*i,k);

            // b
            j = 2;

            mStrain(0,i+j*numActive[uBasis]) += -b_co_V(0,0) * bVals(i,k);
            mStrain(1,i+j*numActive[uBasis]) += -b_co_V(1,0) * bVals(i,k);
            mStrain(2,i+j*numActive[uBasis]) += -b_co_V(2,0) * bVals(i,k);


            //---------------------------------------
            // bending strain first order terms

            // v3-terms
            // c
            j = 2;
            bStrainFirstOrder(0,i+j*numActive[uBasis]) += -c_co(0,0) * bVals(i,k);
            bStrainFirstOrder(1,i+j*numActive[uBasis]) += -c_co(1,1) * bVals(i,k);
            bStrainFirstOrder(2,i+j*numActive[uBasis]) += -c_co(1,0) * bVals(i,k);

            // Gamma
            temp.setZero(3,1);
            for(index_t sigma = 0; sigma!=2; ++sigma)
                temp += FInvT.col(sigma)*bGrads(2*i+sigma,k);

            bStrainFirstOrder.col(i+j*numActive[uBasis]) += -gDerSecond.transpose() * temp;


            // loop over tangential components u_1, u_2 of basis functions
            for (index_t tau = 0; tau!= 2; ++tau)
            {
                const index_t s = tau*numActive[uBasis];

                // membrane strain
                // Gamma
                mStrain.col(i+s) += -gDerSecond.transpose() * FInvT.col(tau)*bVals(i,k);


                //---------------------------------------
                // bending strain first order terms

                /*
                // b old version
                temp = b_mixed(0,0) * gDerSecond.col(0) + b_mixed(1,0) * gDerSecond.col(2);
                bStrainFirstOrder(0,i+s) += 2*(b_mixed(tau,0)*bGrads(2*i,k) - bVals(i,k)*FInvT.col(tau).dot(temp));

                temp = b_mixed(0,1) * gDerSecond.col(2) + b_mixed(1,1) * gDerSecond.col(1);
                bStrainFirstOrder(1,i+s) += 2*(b_mixed(tau,1)*bGrads(2*i+1,k) - bVals(i,k)*FInvT.col(tau).dot(temp));

                temp = b_mixed(0,0) * gDerSecond.col(2) + b_mixed(1,0) * gDerSecond.col(1);
                temp1= b_mixed(0,1) * gDerSecond.col(0) + b_mixed(1,1) * gDerSecond.col(2);
                bStrainFirstOrder(2,i+s) += b_mixed(tau,0)*bGrads(2*i+1,k) - bVals(i,k)*FInvT.col(tau).dot(temp) +
                                            b_mixed(tau,1)*bGrads(2*i,k) - bVals(i,k)*FInvT.col(tau).dot(temp1);
                 */


                // loop over components Strain E_11, E_22, E_12
                for(index_t comp =0; comp!=3; ++comp)
                {
                    index_t alpha = VoigtVtoM(comp,0);
                    index_t beta = VoigtVtoM(comp,1);

                    // b
                    temp.setZero(3,1);
                    temp1.setZero(3,1);
                    for(index_t sigma = 0; sigma!=2; ++sigma)
                    {
                        temp += b_mixed(sigma,alpha) * gDerSecond.col(map_gDerSecond[sigma][beta]);
                        temp1 += b_mixed(sigma,beta) * gDerSecond.col(map_gDerSecond[sigma][alpha]);
                    }

                    bStrainFirstOrder(comp,i+s) += b_mixed(tau,alpha) * bGrads(2*i+beta,k) - bVals(i,k)*FInvT.col(tau).dot(temp) +
                            b_mixed(tau,beta) * bGrads(2*i+alpha,k);
                    // -bVals(i,k)*FInvT.col(tau).dot(temp1);
                    // 2.term derivative b:  bVals(i,k)*FInvT.col(tau).dot(temp1)


                    // derivative b
                    // 1.term
                    temp.setZero(3,1);
                    temp1.setZero(3,1);
                    for(index_t sigma = 0; sigma!=2; ++sigma)
                    {
                        temp += b_co(sigma,beta) * FInvT.col(sigma);
                        temp1 += b_co(sigma,beta) * contraBDer.col(map_contraBDer[alpha][sigma]) + normalDer.col(alpha).dot(gDerSecond.col(map_gDerSecond[sigma][beta])) * FInvT.col(sigma)
                                + normalVec.dot(gDerThird.col(map_gDerThird[alpha][sigma][beta])) * FInvT.col(sigma);
                    }

                    bStrainFirstOrder(comp,i+s) += bVals(i,k) * (contraBDer.col(map_contraBDer[alpha][tau]).dot(temp) + FInvT.col(tau).dot(temp1));
                }


                // 2.term
                /*
                bStrainFirstOrder(0,i+s) += bVals(i,k)*FInvT.col(tau).dot(b_mixed(0,0) * gDerSecond.col(0) + b_mixed(1,0) * gDerSecond.col(2));
                bStrainFirstOrder(1,i+s) += bVals(i,k)*FInvT.col(tau).dot(b_mixed(0,1) * gDerSecond.col(2) + b_mixed(1,1) * gDerSecond.col(1));
                bStrainFirstOrder(2,i+s) += bVals(i,k)*FInvT.col(tau).dot(b_mixed(0,1) * gDerSecond.col(0) + b_mixed(1,1) * gDerSecond.col(2));
                */

                // 3.term
                temp.setZero(3,1);
                for(index_t sigma = 0; sigma!=2; ++sigma)
                    temp += b_mixed(tau,sigma)*bVals(i,k)*FInvT.col(sigma);

                bStrainFirstOrder.col(i+s) += -gDerSecond.transpose() * temp;

            } // end loop tau

            // bending strain hessian
            if(!Mmixed)
            {
                j = 2;
                bStrain(0,i+j*numActive[uBasis]) += bGrads2(3*i,k);
                bStrain(1,i+j*numActive[uBasis]) += bGrads2(3*i+1,k);
                bStrain(2,i+j*numActive[uBasis]) += bGrads2(3*i+2,k);
            }

        } // end loop i


        bStrain += bStrainFirstOrder;

        // Voigt notation 2*E_12
        mStrain.row(2)*=2.0;
        bStrainFirstOrder.row(2) *=2;
        bStrain.row(2)*=2.0;


    }

// implementation of 3rd derivative for Bspline geometry (for d=2)
    void gComputeDeriv3Old(const gsGeometry<T> &geometry)
    {
        // special case parameter space dim =2

        //const gsTensorBasis<2, T>& basis = dynamic_cast<const gsTensorBasis<2, T>& >(geometry.basis());
        const gsTensorBasis<2, T>& basis = dynamic_cast<const gsTensorBasis<2, T>& >(geometry.basis().source());
        //const gsRationalBasis<2, T>& basis = dynamic_cast<const gsRationalBasis<2, T>& >(geometry.basis());


        std::vector< gsMatrix<T> >values[2];
        gsVector<index_t, 2> numActive_dir;
        gsMatrix<index_t> actives_dir[2];

        // loop over parameter directions
        for (index_t i = 0; i < 2; ++i)
        {
            // evaluate basis functions/derivatives of 1D basis
            basis.component(i).evalAllDers_into( md.points.row(i), 3, values[i] );

            // number of basis functions (in each direction)
            numActive_dir(i) = values[i].front().rows();

            gsMatrix<T> quNodeComp_i(1,1);
            quNodeComp_i(0,0) = md.points(i,0);
            basis.component(i).active_into(quNodeComp_i, actives_dir[i]);
        }

        gDerThird_quNodes.setZero(3*4,md.points.cols());

        std::vector< gsMatrix<T> > values0 = values[0];
        std::vector< gsMatrix<T> > values1 = values[1];

        for (index_t i0=0; i0<numActive_dir(0); ++i0)
            for (index_t i1=0; i1<numActive_dir(1); ++i1)
            {
                unsigned active0 = actives_dir[0](i0,0);
                unsigned active1 = actives_dir[1](i1,0);

                gDerThird_quNodes.middleRows(0,3) += geometry.coef(basis.index(active0,active1)).transpose() * values0[3].row(i0).cwiseProduct(values1[0].row(i1));
                gDerThird_quNodes.middleRows(3,3) += geometry.coef(basis.index(active0,active1)).transpose() * values0[0].row(i0).cwiseProduct(values1[3].row(i1));
                gDerThird_quNodes.middleRows(6,3) += geometry.coef(basis.index(active0,active1)).transpose() * values0[2].row(i0).cwiseProduct(values1[1].row(i1));
                gDerThird_quNodes.middleRows(9,3) += geometry.coef(basis.index(active0,active1)).transpose() * values0[1].row(i0).cwiseProduct(values1[2].row(i1));
            }

    }

// uses implementation of 3rd derivative for NURBS geometry (for d=2)
/*
    void gComputeDeriv3(const gsGeometryEvaluator<T> &geoEval, const gsMatrix<T>& quNodes)
    {
        //const gsRationalBasis<gsTensorBSplineBasis<2, T> >& basis = dynamic_cast<const gsRationalBasis< gsTensorBSplineBasis<2, T> >& >(geoEval.geometry().basis());
        //const gsTensorBasis<2, T>& basis = dynamic_cast<const gsTensorBasis<2, T>& >(geoEval.geometry().basis());

        const gsBasis<T>& basis = geoEval.geometry().basis();

        gsMatrix<T> B;
        gsMatrix<index_t> actives;

        // compute third derivatives
        basis.deriv3_into(quNodes,B);
        // compute active functions
        basis.active_into(quNodes,actives);

        // compute result as linear combination of
        // "coefs(actives)" and B
        gsBasis<T>::linearCombination_into( geoEval.geometry().coefs(), actives, B, gDerThird_quNodes);

    }
*/






protected:
    // options
    bool Mmixed;
    bool Nmixed;

    // Assembler
    const gsShellMixedAssembler<T>& m_assembler;

    // Bases
    size_t pBasis;
    size_t phiBasis;
    size_t N11Basis;
    size_t N22Basis;
    size_t N12Basis;
    size_t uBasis;
    size_t basesSize;

    // Pde
    const gsShellMixedPde<T>* m_pde_ptr;

    // Pointer to surface force
    const gsFunction<T> * m_force;
    // local values of the surface force
    gsMatrix<T> forceVals;

    // Geometry quantities
    gsVector<T> normalVec; // Normal to the shell centerline
    gsMatrix<T,3,2> normalDer;
    gsMatrix<T,3,3> F; // deformation gradient (a_1, a_2, a_3)
    gsMatrix<T,3,3> FInvT; // (a^1, a^2, a^3)

    gsMatrix<T,3,2> gDerFirst;
    gsMatrix<T,3,3> gDerSecond;
    gsMatrix<T,3,4> gDerThird;
    gsMatrix<T> gDerThird_quNodes;

    gsMatrix<T,2,2> a_co; // first ff covariant components
    gsMatrix<T,2,2> a_contra; // first ff contravariant components

    gsMatrix<T,3,1> b_co_V; // second ff covariant components Voigt notation (b_11, b_22, b_12)^T
    gsMatrix<T,2,2> b_co;
    gsMatrix<T,2,2> b_mixed;

    gsMatrix<T,2,2> c_co; // third ff covariant components

    gsMatrix<T,3,4> contraBDer; // first derivative of the contravariant basis a^1, a^2

    // auxiliary maps
    index_t map_gDerSecond [2][2];
    index_t map_gDerThird [2][2][2];
    index_t map_contraBDer [2][2];

    // Voigt notation
    gsMatrix<index_t, 3, 2> VoigtVtoM;

    // Material matrix
    gsMatrix<T,3,3> m_C;
    gsMatrix<T,3,3> m_CModInv;
    gsMatrix<T,3,3> m_CModInvN;


    // Mixed formulation quantities
    gsMatrix<T,3> symCurlPhi;
    gsMatrix<T,3> pId;
    gsMatrix<T,3> N;
    std::vector<gsMatrix<T> > bGrads_k;

    // Shell quantities
    gsMatrix<T,3> mStrain; // membrane strain epsilon
    gsMatrix<T,3> bStrainFirstOrder; // first order terms bending strain kappa
    gsMatrix<T,3> bStrain; // bending strain kappa


    // Basis values
    std::vector<std::vector<gsMatrix<T> > > basisData;
    std::vector<gsMatrix<index_t> > actives;
    std::vector<index_t >           numActive;
    index_t numActiveN;
    std::vector<gsMatrix<index_t> > actives_vec;

    // Local matrices
    gsMatrix<T> localMatA00, localMatA10, localMatA11;
    gsMatrix<T> localMatB0, localMatB1;
    gsMatrix<T> localMatC, localMatC00, localMatC10;
    gsMatrix<T> localRhs;

    gsMapData<T> md;
};


} // namespace gismo

