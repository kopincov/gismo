/** @file gsVisitorPlateMixedSimplySupp.h

    @brief Boundary element visitor for the contribution of free boundaries

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/

#pragma once
#include <gsPlateMixed/gsPlateMixedAssembler.h>
#include <gsPde/gsShellMixedPde.h>
#include <gsCore/gsFunctionExpr.h>

namespace gismo
{

template <class T>
class gsVisitorPlateMixedSimplySupp
{
public:

    gsVisitorPlateMixedSimplySupp(const gsPde<T> & pde, const boundary_condition<T> & s, gsPlateMixedAssembler<T>& assembler)
        : assembler(assembler), pde(dynamic_cast<const gsShellMixedPde<T>& >(pde)), side(s.side()), data_ptr(s.function().get())
    { }


    void initialize(const gsBasisRefs<T> & bases,
                    const index_t          patchIndex,
                    const gsOptionList   & options,
                          gsQuadRule<T>  & rule)
    {
        // Set options
        phiBcMethod = assembler.options().getInt("phiBcMethodSs");
        wHomogenPsiq = assembler.options().getSwitch("wHomogenPsiq");
        wPenaltyPsiq = assembler.options().getSwitch("wPenaltyPsiq");

        bQ = assembler.system().colBasis(0);
        bPsi = assembler.system().colBasis(1);
        lambdaBasis = assembler.system().colBasis(4);
        basesSize = bases.size();

        basisData.resize(basesSize);
        actives.resize(basesSize);
        physGrad.resize(basesSize);
        numActive.resize(basesSize);

        // Setup Quadrature (harmless slicing occurs)
        rule = gsGaussRule<T>(bases[bPsi], options.getReal("quA"), options.getInt("quB"),
                              side.direction()); // use phi basis, because higher degree

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_OUTER_NORMAL | NEED_GRAD_TRANSFORM;

        // Compute penalty parameter
        const int deg = bases[bPsi].maxDegree();
        penalty = (deg + bases[bPsi].dim()) * (deg + 1) * T(2.5);

        m_patchIndex = patchIndex;
    }


    // Evaluate on element.
    inline void evaluate(const gsBasisRefs<T> & bases, // to do: more unknowns
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

        // Compute function at corners
        side.getContainedCorners(2, corners);
        if (sideOrientation(side) == -1) // clock-wise
        {
            std::vector<boxCorner> corners_temp(corners);
            corners[0] = corners_temp[1];
            corners[1] = corners_temp[0];
        }

        cornerIndex.resize(2);
        cornerIndex[0] = bases[lambdaBasis].functionAtCorner(corners[0]);
        cornerIndex[1] = bases[lambdaBasis].functionAtCorner(corners[1]);

        // Compute geometry related values
        geo.computeMap(md);

        // Evaluate data
        data_ptr->eval_into(md.values[0], dataVals);
        gradVals = dataVals.topRows(2);
        MVals = dataVals.middleRows(2, 3);

        // non-homogeneous BC M
        gsFunctionExpr<T> Mnn_yint, Mnt_yint;
        if (side == boundary::east || side == boundary::west)
        {
            // sol 1
            //Mnn_yint=gsFunctionExpr<T>("Pi^2*Cos(Pi*x)*Sin(2*Pi*y)*1/(2*Pi)",2);
            //Mnt_yint=gsFunctionExpr<T>("2*Pi^2*Sin(Pi*x)*Cos(2*Pi*y)*1/(2*Pi)",2);


            // sol 2
            //Mnn_yint=gsFunctionExpr<T>("-24*x^2*y^5*1/5 + Pi^2*Cos(Pi*x)*Sin(2*Pi*y)*1/(2*Pi)",2);
            //Mnt_yint=gsFunctionExpr<T>("-32*x^3*y^4*1/4 + 2*Pi^2*Sin(Pi*x)*Cos(2*Pi*y)*1/(2*Pi)",2);


        }
        else if (side == boundary::north || side == boundary::south)
        {

            //Mnn_yint=gsFunctionExpr<T>("Pi^2*Cos(Pi*x)*Sin(2*Pi*y)*1/(2*Pi)",2);
            //Mnt_yint=gsFunctionExpr<T>("-2*Pi^2*Sin(Pi*x)*Cos(2*Pi*y)*1/(2*Pi)",2);

            //Mnn_yint=gsFunctionExpr<T>("-24*x^2*y^5*1/5 + Pi^2*Cos(Pi*x)*Sin(2*Pi*y)*1/(2*Pi)",2);
            //Mnt_yint=gsFunctionExpr<T>("+32*x^3*y^4*1/4 - 2*Pi^2*Sin(Pi*x)*Cos(2*Pi*y)*1/(2*Pi)",2);

        }

        Mnn_yint = gsFunctionExpr<T>("0", 2);
        Mnt_yint = gsFunctionExpr<T>("0", 2);

        Mnn_yint.eval_into(md.values[0], Mnn_yintVals);
        Mnt_yint.eval_into(md.values[0], Mnt_yintVals);


        // Initialize local matrix/rhs
        localMatA11Nitsche.setZero(2 * numActive[bPsi], 2 * numActive[bPsi]);
        MatA10Nitsche_psip.setZero(assembler.numActivesFreeBoundary, 2 * numActive[bPsi]);

        MatA00_psiq_psip.setZero(assembler.numActivesFreeBoundary, assembler.numActivesFreeBoundary);
        MatA01_psiq.setZero(assembler.numActivesFreeBoundary, 2 * numActive[bPsi]);

        localMatLn1.setZero(numActive[lambdaBasis], 2 * numActive[bPsi]);

        localRhsfpsi.setZero(2 * numActive[bPsi], 1);
        localRhsfmuN.setZero(numActive[bQ], 1);
    }

    inline void assemble(gsDomainIterator<T> & element,
                   const gsVector<T>         & quWeights)
    {
        std::vector<gsMatrix<T> > bVals(basesSize); // ursprünglichen version reference
        for(size_t i=0; i<basesSize; i++)
            bVals[i] = basisData[i][0];

        std::vector<gsMatrix<T> > bGrads(basesSize);
        for(size_t i=0; i<basesSize; i++)
            bGrads[i] = basisData[i][1];


        gsMatrix<T> psiVals;
        psiVals.setZero(numActive[bQ], 2);

        //gsMatrix<T, 2, 2> idMatrix = gsMatrix<T, 2, 2>::Identity();
        gsMatrix<T,2,2> jac_i;
        gsMatrix<T,2,2> jac_j;
        gsMatrix<T,2,2> Curl_i;
        gsMatrix<T,2,2> Curl_j;
        gsMatrix<T,2,2> symCurl_i;
        gsMatrix<T,2,2> symCurl_j;
        gsMatrix<T,2,1> bValsVec_i;
        gsMatrix<T,2,1> bValsVec_j;

        for (index_t k = 0; k < quWeights.rows(); ++k)
        {
            // compute physical gradients at k as a Dim x NumActive matrix
            for(size_t i=0; i<basesSize; i++)
                transformGradients(md, k, bGrads[i], physGrad[i]);

            // Compute the outer normal vector on the side
            outerNormal(md, k, side, unormal);

            // Multiply quadrature weight by the measure of normal
            const T weight = quWeights[k] * unormal.norm();

            // Get penalty parameter
            T mu = 0;

            switch(phiBcMethod)
            {
            case gsPlateMixedAssembler<T>::nitsche: mu = penalty / element.getCellSize(); break; // Nitsche prescribe value
            case gsPlateMixedAssembler<T>::nitscheDerivative: mu = penalty * element.getCellSize(); break; // Nitsche prescribe derivative
            }

            // Compute unit outer normal and tangent
            unormal.normalize();
            gsVector<T,2> tangent;
            tangent(0) = -unormal(1);
            tangent(1) = unormal(0);

            // Set M_mat at quadrature point k
            gsMatrix<T,2,2> M_mat;
            M_mat(0,0) = MVals(0,k);
            M_mat(1,1) = MVals(1,k);
            M_mat(1,0) = M_mat(0,1) = MVals(2,k);

            psiP_inhomM << Mnn_yintVals(0,k) * unormal + Mnt_yintVals(0,k) * tangent;
            //psiP_inhomM << Mnn_yintVals(0,k) * unormal;
            //psiP_inhomM += weight * (M_mat*unormal).dot(unormal) * unormal;
            //psiP_inhomM += weight * ((M_mat*unormal).dot(unormal) * unormal + (M_mat*unormal).dot(tangent) * tangent);

            // A11
            for (index_t i = 0; i < numActive[bPsi]; i++)
                for (index_t j = 0; j < numActive[bPsi]; j++)
                    for(index_t m = 0; m<2; m++)
                        for(index_t l=0; l<2; l++)
                        {
                            Curl_i = gsMatrix<T, 2, 2>::Zero();
                            Curl_i(m,0)=physGrad[bPsi](1,i);
                            Curl_i(m,1)=-physGrad[bPsi](0,i);
                            symCurl_i = 0.5*(Curl_i+Curl_i.transpose());
                            Curl_j = gsMatrix<T, 2, 2>::Zero();
                            Curl_j(l,0)=physGrad[bPsi](1,j);
                            Curl_j(l,1)=-physGrad[bPsi](0,j);
                            symCurl_j = 0.5*(Curl_j+Curl_j.transpose());

                            jac_i = gsMatrix<T, 2, 2>::Zero();
                            jac_i.row(m)=physGrad[bPsi].col(i);
                            jac_j = gsMatrix<T, 2, 2>::Zero();
                            jac_j.row(l)=physGrad[bPsi].col(j);

                            bValsVec_i = gsMatrix<T, 2, 1>::Zero();
                            bValsVec_i(m)=bVals[bPsi](i,k);
                            bValsVec_j = gsMatrix<T, 2, 1>::Zero();
                            bValsVec_j(l)=bVals[bPsi](j,k);

                            // Nitsche
                            switch(phiBcMethod)
                            {
                            case gsPlateMixedAssembler<T>::nitsche: localMatA11Nitsche(m*numActive[bPsi] + i, l*numActive[bPsi] + j) +=   mu * weight * bValsVec_i.dot(unormal) * bValsVec_j.dot(unormal); break; // p(phi.n, psi.n)
                            case gsPlateMixedAssembler<T>::nitscheDerivative: localMatA11Nitsche(m*numActive[bPsi] + i, l*numActive[bPsi] + j) +=   mu * weight * ((jac_i*tangent).dot(unormal) * (jac_j*tangent).dot(unormal)); break; // p((\grad phi)_tn, (\grad psi)_tn)
                            }

                            localMatA11Nitsche(m*numActive[bPsi] + i, l*numActive[bPsi] + j) += weight * ((symCurl_j * tangent).dot(unormal) * bValsVec_i.dot(unormal) +  (symCurl_i * tangent).dot(unormal) * bValsVec_j.dot(unormal)); // s(phi, psi.n) + s(psi, phi.n)
                        }

            // psiq,p term and psiq,phi term (previous elements)
            // A01
            for ( index_t i =0; i<assembler.numActivesFreeBoundary; i++)
                for (index_t j = 0; j < numActive[bPsi]; j++)
                    for(index_t l=0; l<2; l++)
                    {

                        Curl_j = gsMatrix<T, 2, 2>::Zero();
                        Curl_j(l,0)=physGrad[bPsi](1,j);
                        Curl_j(l,1)=-physGrad[bPsi](0,j);
                        symCurl_j = 0.5*(Curl_j+Curl_j.transpose());

                        jac_j = gsMatrix<T, 2, 2>::Zero();
                        jac_j.row(l)=physGrad[bPsi].col(j);

                        bValsVec_j = gsMatrix<T, 2, 1>::Zero();
                        bValsVec_j(l)=bVals[bPsi](j,k);

                        if ( assembler.system().rowMapper(0).is_free_index(assembler.activesFreeBoundaryMap(i)) )
                        {
                            // Homogenization
                            if(wHomogenPsiq)
                            {
                                MatA01_psiq(i,l*numActive[bPsi] + j) += -weight * assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(unormal) * (symCurl_j * tangent).dot(unormal); // s(phi, psi[q]) (.,g)
                                //MatA01_psiq(i,l*numActive[bPsi] + j) += -weight * assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(symCurl_j * tangent); // s(phi, psi[q]) (.,g)    // modification
                            }

                            if(wPenaltyPsiq)
                            {
                                MatA01_psiq(i,l*numActive[bPsi] + j) += -mu * weight * assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(unormal) * bValsVec_j.dot(unormal); // p(phi, psi[q]) (.,g)
                            }

                            // Nitsche
                            MatA10Nitsche_psip(i,l*numActive[bPsi] + j) += -weight * assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(unormal) * (symCurl_j * tangent).dot(unormal); // s(psi, psi[p].n) undo symmetrization

                            if(phiBcMethod == gsPlateMixedAssembler<T>::nitsche)
                                MatA10Nitsche_psip(i,l*numActive[bPsi] + j) += -mu * weight* assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(unormal) * bValsVec_j.dot(unormal); // p(psi[p].n, psi.n)

                        }
                    }

            // A00
            // Homogenization
            if(wPenaltyPsiq)
            {
                for ( index_t i =0; i<assembler.numActivesFreeBoundary; i++)
                    for ( index_t j =0; j<assembler.numActivesFreeBoundary; j++)
                        if ( assembler.system().rowMapper(0).is_free_index(assembler.activesFreeBoundaryMap(i)) )
                            if ( assembler.system().rowMapper(0).is_free_index(assembler.activesFreeBoundaryMap(j)) )
                                MatA00_psiq_psip(i, j) += mu * weight * assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(unormal) * assembler.psiq.row(assembler.activesFreeBoundaryMap(j)).dot(unormal); // p(psi[p], psi[q]) (g,g)
            }



            // Lagrange multipliers
            real_t lambdaNCoef_Ln;

            for ( index_t i =0; i<numActive[lambdaBasis]; i++)
            {
                // set coefficients for lambdaN, lambdaT at corners

                // setting if not at corner
                lambdaNCoef_Ln = 1;
                //lambdaTCoef_Ln = 0 because BC on ss

                if(actives[lambdaBasis](i,0) == cornerIndex[0]) // start-point edge (counter clock-wise): edge k wrt. corner
                {
                    switch(pde.m_corners[m_patchIndex][corners[0]])
                    {
                    case gsShellMixedPde<T>::sf: lambdaNCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[0]](0,3); break;
                        //lambdaTCoef_Ln = 0 because BC on ss
                    default: break;
                    }
                }
                else if(actives[lambdaBasis](i,0) == cornerIndex[1]) // end-point edge (counter clock-wise): edge k-1 wrt. corner
                {
                    switch(pde.m_corners[m_patchIndex][corners[1]])
                    {
                    case gsShellMixedPde<T>::sf: lambdaNCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[1]](0,2); break;
                        //lambdaTCoef_Ln = 0 because BC on ss
                    default: break;
                    }

                }

                for (index_t j = 0; j < numActive[bPsi]; j++)
                    for(index_t l=0; l<2; l++)
                    {
                        jac_j = gsMatrix<T, 2, 2>::Zero();
                        jac_j.row(l)=physGrad[bPsi].col(j);

                        // Ln1
                        localMatLn1(i,l*numActive[bPsi] + j) += weight * ((jac_j *tangent).dot(unormal) * lambdaNCoef_Ln) * bVals[lambdaBasis](i,k); //lambdaTCoef_Ln = 0 because BC on ss

                    }

                // non-homogeneous BC M
                // fmuN
                localRhsfmuN(i,0) += weight * (M_mat*unormal).dot(unormal) * lambdaNCoef_Ln * bVals[lambdaBasis](i,k);

            }

            // non-homogeneous BC M
            // Nitsche
            // fpsi
            if(phiBcMethod == gsPlateMixedAssembler<T>::nitsche)
            {
                for (index_t i = 0; i < numActive[bPsi]; i++)
                    for(index_t m = 0; m<2; m++)
                    {
                        Curl_i = gsMatrix<T, 2, 2>::Zero();
                        Curl_i(m,0)=physGrad[bPsi](1,i);
                        Curl_i(m,1)=-physGrad[bPsi](0,i);
                        symCurl_i = 0.5*(Curl_i+Curl_i.transpose());

                        bValsVec_i = gsMatrix<T, 2, 1>::Zero();
                        bValsVec_i(m)=bVals[bPsi](i,k);

                        localRhsfpsi(m*numActive[bPsi] + i,0) += weight * psiP_inhomM.dot(unormal) * (symCurl_i * tangent).dot(unormal); // s(psi, psi[p].n) undo symmetrization
                        localRhsfpsi(m*numActive[bPsi] + i,0) += mu * weight* psiP_inhomM.dot(unormal) * bValsVec_i.dot(unormal); // p(psi[p].n, psi.n)
                    }

            }
            
            
            // non-homogeneous BC w
            // fpsi
            for (index_t i = 0; i < numActive[bPsi]; i++)
                for(index_t m = 0; m<2; m++)
                {
                    jac_i = gsMatrix<T, 2, 2>::Zero();
                    jac_i.row(m)=physGrad[bPsi].col(i);

                    localRhsfpsi(m*numActive[bPsi] + i,0) += -weight * ((jac_i *tangent).dot(tangent) * gradVals.col(k).dot(tangent));

                }
            


        }

    }

    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> >   & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        actives_vec.resize(5);

        system.mapColIndices(actives[bQ], patchIndex, actives_vec[0], 0);
        system.mapColIndices(actives[bPsi], patchIndex, actives_vec[1], 1);
        system.mapColIndices(actives[bPsi], patchIndex, actives_vec[2], 2);
        system.mapColIndices(actives[bQ], patchIndex, actives_vec[3], 3);
        if(phiBcMethod == gsPlateMixedAssembler<T>::lagrange)
            system.mapColIndices(actives[lambdaBasis], patchIndex, actives_vec[4], 4);

        actives_phi1 = actives_vec[1];
        actives_phi2 = actives_vec[2];

        // build block information
        gsVector<index_t> phi_vec(2);
        gsVector<index_t> lambdaN_vec(1);
        phi_vec << 1,2;
        lambdaN_vec << 4;

        // A11
        if(phiBcMethod == gsPlateMixedAssembler<T>::nitsche || phiBcMethod == gsPlateMixedAssembler<T>::nitscheDerivative)
        {
            // Nitsche terms
            system.pushToMatrix(localMatA11Nitsche, actives_vec, eliminatedDofs, phi_vec, phi_vec);

            // Nitsche psip terms
            system.pushToMatrix(MatA10Nitsche_psip.block(0,0,assembler.numActivesFreeBoundary,numActive[bPsi]).transpose(), actives_phi1, assembler.activesFreeBoundaryMap, eliminatedDofs[0], 1, 0);
            system.pushToMatrix(MatA10Nitsche_psip.block(0,numActive[bPsi],assembler.numActivesFreeBoundary,numActive[bPsi]).transpose(), actives_phi2, assembler.activesFreeBoundaryMap, eliminatedDofs[0], 2, 0);


            // Homogenization: psi_q terms
            // A00
            system.pushToMatrix(MatA00_psiq_psip, assembler.activesFreeBoundaryMap, assembler.activesFreeBoundaryMap, eliminatedDofs[0], 0, 0);
            // A01
            system.pushToMatrix(MatA01_psiq.block(0,0,assembler.numActivesFreeBoundary,numActive[bPsi]), assembler.activesFreeBoundaryMap, actives_phi1, eliminatedDofs[1], 0, 1);
            system.pushToMatrix(MatA01_psiq.block(0,numActive[bPsi],assembler.numActivesFreeBoundary,numActive[bPsi]), assembler.activesFreeBoundaryMap, actives_phi2, eliminatedDofs[2], 0, 2);
        }

        // Ln1
        if(phiBcMethod == gsPlateMixedAssembler<T>::lagrange)
        {
            system.pushToMatrix(localMatLn1, actives_vec, eliminatedDofs, lambdaN_vec, phi_vec);
            system.pushToMatrix(localMatLn1.transpose(), actives_vec, eliminatedDofs, phi_vec, lambdaN_vec);

            // fmuN
            system.pushToRhs(localRhsfmuN, actives_vec, lambdaN_vec);

        }
        
        // fpsi
        system.pushToRhs(localRhsfpsi, actives_vec, phi_vec);


    }


protected:


    int phiBcMethod;
    bool wHomogenPsiq;
    bool wPenaltyPsiq;

    // Additional data for VisitorPlateMixedNewPsiq
    gsPlateMixedAssembler<T>& assembler;
    const gsShellMixedPde<T>& pde;
    index_t m_patchIndex;
    size_t bQ;
    size_t bPsi;
    size_t lambdaBasis;
    size_t basesSize;

    gsVector<T,2> psiP_inhomM;

    // Neumann function
    boxSide side;
    const gsFunction<T> * data_ptr;
    gsMatrix<T> dataVals;
    gsMatrix<T> gradVals;
    gsMatrix<T> MVals;
    gsMatrix<T> Mnn_yintVals;
    gsMatrix<T> Mnt_yintVals;
    real_t penalty;

    // Basis values
    std::vector<std::vector<gsMatrix<T> > > basisData;
    std::vector<gsMatrix<index_t> > actives;
    std::vector<gsMatrix<T> >		 physGrad;
    std::vector<index_t >            numActive;
    std::vector<boxCorner> corners;
    std::vector<index_t> cornerIndex;

    std::vector<gsMatrix<index_t> > actives_vec;
    gsMatrix<index_t> actives_phi1;
    gsMatrix<index_t> actives_phi2;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> pVals;


    // Local matrix
    //Nitsche
    gsMatrix<T> localMatA11Nitsche;
    gsMatrix<T> MatA10Nitsche_psip;

    // psiq
    gsMatrix<T> MatA00_psiq;
    gsMatrix<T> MatA00_psip;
    gsMatrix<T> MatA00_psiq_psip;
    gsMatrix<T> MatA01_psiq;

    // Lagrange multiplier
    gsMatrix<T> localMatLn1;
    gsMatrix<T> localMatD;
    
    // Homogenization w
    gsMatrix<T> localRhsfq;
    gsMatrix<T> localRhsfpsi;
    gsMatrix<T> localRhsfmuN;

    gsMapData<T> md;
};

} // namespace gismo
