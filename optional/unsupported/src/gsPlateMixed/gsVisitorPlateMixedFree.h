/** @file gsVisitorPlateMixedFree.h

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

namespace gismo
{

template <class T>
class gsVisitorPlateMixedFree
{
public:
    gsVisitorPlateMixedFree(const gsPde<T> & pde, const boundary_condition<T> & s, gsPlateMixedAssembler<T>& assembler)
        : assembler(assembler), pde(dynamic_cast<const gsShellMixedPde<T>& >(pde)), side(s.side()), data_ptr(s.function().get())
    {

    }

    void initialize(const gsBasisRefs<T> & bases,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        // Set options
        phiBcMethod = assembler.options().getInt("phiBcMethodF");
        wHomogenPsiq = assembler.options().getSwitch("wHomogenPsiq");
        wPenaltyPsiq = assembler.options().getSwitch("wPenaltyPsiq");
        pBoundary0 = assembler.options().getSwitch("pBoundary0");
        lambdaTMean0 = assembler.options().getSwitch("lambdaTMean0");

        bQ = assembler.system().colBasis(0);
        bPsi = assembler.system().colBasis(1);
        lambdaBasis = assembler.system().colBasis(4);
        basesSize = bases.size();

        basisData.resize(basesSize);
        basisData_corners.resize(basesSize);
        actives.resize(basesSize);
        physGrad.resize(basesSize);
        numActive.resize(basesSize);

        // Setup Quadrature (harmless slicing occurs)
        rule = gsGaussRule<T>(bases[bPsi],
                              options.getReal("quA"),
                              options.getInt("quB"),
                              side.direction()); // use phi basis, because higher degree

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_OUTER_NORMAL | NEED_GRAD_TRANSFORM;

        // Compute penalty parameter
        const int deg = bases[bPsi].maxDegree();
        penalty = (deg + bases[bPsi].dim()) * (deg + 1) * T(2.5);

        m_patchIndex = patchIndex;

        int_DivM = 0;
        psiP_inhomM << 0, 0;

        firstTime = true;
        geoEval_corners = memory::make_shared(getEvaluator(md.flags, pde.patches()[patchIndex]));
    }

    // Evaluate on element.
    inline void evaluate(const gsBasisRefs<T> & bases,
                         const gsGeometry<T>  & geo,
                               gsMatrix<T>    & quNodes)
    {
        md.points = quNodes;
        orientation = geo.orientation();

        // Flip quadrature nodes if necessary (s.t. counter-clockwise)
        gsMatrix<T> quNodes_tmp(quNodes);
        if((sideOrientation(side) * geo.orientation()) == -1)
        {
            for(index_t i=0; i<md.points.cols(); i++)
            {
                if(side.direction() == 0)
                    md.points(1,i) = 1 - quNodes_tmp(1,md.points.cols()-1-i);
                else
                    md.points(0,i) = 1 - quNodes_tmp(0,md.points.cols()-1-i);
            }
        }

        m_quNodes = md.points;

        if(firstTime==true)
        {
            // Compute function at corners
            side.getContainedCorners(2, corners); // ordered counter clock-wise
            if(sideOrientation(side) == -1) // if clock-wise, swap
            {
                std::vector<boxCorner> corners_temp(corners);
                corners[0] = corners_temp[1];
                corners[1] = corners_temp[0];
            }

            cornerIndex.resize(2);
            cornerIndex[0] = bases[lambdaBasis].functionAtCorner(corners[0]);
            cornerIndex[1] = bases[lambdaBasis].functionAtCorner(corners[1]);

            gsVector<bool> cornerParameters = corners[0].parameters(2);
            cornersMatrix(0,0)= cornerParameters(0);
            cornersMatrix(1,0)= cornerParameters(1);
            cornerParameters = corners[1].parameters(2);
            cornersMatrix(0,1)= cornerParameters(0);
            cornersMatrix(1,1)= cornerParameters(1);

            //std::cout<<"cornersMatrix \n "<<cornersMatrix<<std::endl;

            //geoEval_corners = geoEval.; // make copy

            firstTime=false;
        }


        for(size_t i=0; i<basesSize; i++)
        {
            // Compute the active basis functions
            // Assumes actives are the same for all quadrature points on the elements
            bases[i].active_into(md.points.col(0), actives[i]);
            numActive[i]   = actives[i].rows();

            // Evaluate basis functions on element
            bases[i].evalAllDers_into( md.points, 1, basisData[i]);
            bases[i].evalAllDers_into( cornersMatrix, 1, basisData_corners[i]);
        }


        // Compute geometry related values
        geo.computeMap(md);
        geoEval_corners->evaluateAt(cornersMatrix);

        // Evaluate data
        data_ptr->eval_into(md.values[0], dataVals);
        MVals = dataVals.topRows(3);
        gradMVals = dataVals.middleRows(3,6);


        // non-homogeneous BC M
        gsFunctionExpr<T> Mnn_yint, Mnt_yint;

        if(side == boundary::east || side == boundary::west)
        {

            // sol1
            //Mnn_yint=gsFunctionExpr<T>("Pi^2*Cos(Pi*x)*Sin(2*Pi*y)*1/(2*Pi)",2);
            //Mnt_yint=gsFunctionExpr<T>("2*Pi^2*Sin(Pi*x)*Cos(2*Pi*y)*1/(2*Pi)",2);


            // sol2
            //Mnn_yint=gsFunctionExpr<T>("-24*x^2*y^5*1/5 + Pi^2*Cos(Pi*x)*Sin(2*Pi*y)*1/(2*Pi)",2);
            //Mnt_yint=gsFunctionExpr<T>("-32*x^3*y^4*1/4 + 2*Pi^2*Sin(Pi*x)*Cos(2*Pi*y)*1/(2*Pi)",2);

        }
        else if(side == boundary::north || side == boundary::south)
        {

            //Mnn_yint=gsFunctionExpr<T>("Pi^2*Cos(Pi*x)*Sin(2*Pi*y)*1/(2*Pi)",2);
            //Mnt_yint=gsFunctionExpr<T>("-2*Pi^2*Sin(Pi*x)*Cos(2*Pi*y)*1/(2*Pi)",2);


            //Mnn_yint=gsFunctionExpr<T>("-24*x^2*y^5*1/5 + Pi^2*Cos(Pi*x)*Sin(2*Pi*y)*1/(2*Pi)",2);
            //Mnt_yint=gsFunctionExpr<T>("+32*x^3*y^4*1/4 - 2*Pi^2*Sin(Pi*x)*Cos(2*Pi*y)*1/(2*Pi)",2);
        }

        Mnn_yint=gsFunctionExpr<T>("0",2);
        Mnt_yint=gsFunctionExpr<T>("0",2);

        Mnn_yint.eval_into(md.values[0], Mnn_yintVals);
        Mnt_yint.eval_into(md.values[0], Mnt_yintVals);


        // sol2
        //gsFunctionExpr<T> u=gsFunctionExpr<>("Cos(Pi*x)*Cos(2*Pi*y)+2*x^4*y^4",2);
        //u.eval_into(geoEval_corners->values(), uVals);


        // Initialize local matrix/rhs
        localMatA11Nitsche.setZero(2*numActive[bPsi], 2*numActive[bPsi]);
        localMatA10Nitsche_p.setZero(2*numActive[bPsi], numActive[bQ]);
        localMatA10Nitsche_psip.setZero(2*numActive[bPsi], numActive[bQ]);
        MatA10Nitsche_psip.setZero(assembler.numActivesFreeBoundary, 2*numActive[bPsi]);

        localMatA00_psiq.setZero(numActive[bQ], numActive[bQ]);
        localMatA01_psiq.setZero(numActive[bQ], 2*numActive[bPsi]);

        MatA00_psiq.setZero(assembler.numActivesFreeBoundary, numActive[bQ]);
        MatA00_psip.setZero(numActive[bQ], assembler.numActivesFreeBoundary);
        MatA00_psiq_psip.setZero(assembler.numActivesFreeBoundary, assembler.numActivesFreeBoundary);
        MatA01_psiq.setZero(assembler.numActivesFreeBoundary, 2*numActive[bPsi]);


        localMatwN.setZero(numActive[bQ], numActive[bQ]);

        localMatLn0.setZero(numActive[lambdaBasis], numActive[bQ]);
        localMatLn1.setZero(numActive[lambdaBasis], 2*numActive[bPsi]);
        localMatLt0.setZero(numActive[lambdaBasis], numActive[bQ]);
        localMatLt1.setZero(numActive[lambdaBasis], 2*numActive[bPsi]);
        localMatLnMean.setZero(1, numActive[lambdaBasis]);
        localMatLtMean.setZero(1, numActive[lambdaBasis]);

        localRhsfq.setZero(numActive[bQ], 1);
        localRhsfpsi.setZero(2*numActive[bPsi], 1);
        localRhsfv.setZero(numActive[bQ], 1);
        localRhsfmuN.setZero(numActive[bQ], 1);
        localRhsfmuT.setZero(numActive[bQ], 1);

        Rhsfq.setZero(assembler.numActivesFreeBoundary, 1);

        localMatB1.setZero(numActive[bQ], 2*numActive[bPsi]);
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

        gsMatrix<T, 2, 2> idMatrix = gsMatrix<T, 2, 2>::Identity();
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

            gsMatrix<T,2,1> DivM;
            DivM(0,0) = gradMVals(0,k) + gradMVals(5,k);
            DivM(1,0) = gradMVals(2,k) + gradMVals(4,k);

            gsMatrix<T,2,2> d1M_mat;
            d1M_mat(0,0) = gradMVals(0,k);
            d1M_mat(1,1) = gradMVals(1,k);
            d1M_mat(1,0) = d1M_mat(0,1) = gradMVals(2,k);

            gsMatrix<T,2,2> d2M_mat;
            d2M_mat(0,0) = gradMVals(3,k);
            d2M_mat(1,1) = gradMVals(4,k);
            d2M_mat(1,0) = d2M_mat(0,1) = gradMVals(5,k);

            real_t dtMnt = (d1M_mat*unormal).dot(tangent) * tangent(0) + (d2M_mat*unormal).dot(tangent) * tangent(1);

            psiP_inhomM << Mnn_yintVals(0,k) * unormal;
            //psiP_inhomM << Mnn_yintVals(0,k) * unormal + Mnt_yintVals(0,k) * tangent;
            //psiP_inhomM += weight * ((M_mat*unormal).dot(unormal) * unormal + (M_mat*unormal).dot(tangent) * tangent);


            // Psiq
            //if(phiBcMethod == gsPlateMixedAssembler<T>::nitsche)
            {
                gsVector<T> upperCorner(m_quNodes.col(k));
                psiVals= assembler.applyVisitorPsiq(0, side, upperCorner, element);
            }


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
                            case gsPlateMixedAssembler<T>::nitsche: localMatA11Nitsche(m*numActive[bPsi] + i, l*numActive[bPsi] + j) +=   mu * weight * bValsVec_i.dot(bValsVec_j); break; // p(phi, psi)
                            case gsPlateMixedAssembler<T>::nitscheDerivative: localMatA11Nitsche(m*numActive[bPsi] + i, l*numActive[bPsi] + j) +=  mu * weight * ((jac_i*tangent).dot(unormal) * (jac_j*tangent).dot(unormal)); // p((\grad phi)_tn, (\grad psi)_tn)
                                localMatA11Nitsche(m*numActive[bPsi] + i, l*numActive[bPsi] + j) +=   mu * weight * ((jac_i*tangent).dot(tangent) * (jac_j*tangent).dot(tangent)); break; // // p((\grad phi)_tt, (\grad psi)_tt)
                            }

                            localMatA11Nitsche(m*numActive[bPsi] + i, l*numActive[bPsi] + j) += weight * (bValsVec_i.dot(symCurl_j * tangent) + bValsVec_j.dot(symCurl_i * tangent)); // s(phi, psi) + s(psi, phi)

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

                        jac_i = gsMatrix<T, 2, 2>::Zero();
                        jac_i.row(m)=physGrad[bPsi].col(i);

                        bValsVec_i = gsMatrix<T, 2, 1>::Zero();
                        bValsVec_i(m)=bVals[bPsi](i,k);

                        // Nitsche
                        switch(phiBcMethod)
                        {
                        case gsPlateMixedAssembler<T>::nitsche: localMatA10Nitsche_psip(m*numActive[bPsi] + i, j) += -mu * weight* psiVals.row(j).dot(bValsVec_i); break; // p(psi[p], psi) (l,.)
                        case gsPlateMixedAssembler<T>::nitscheDerivative: localMatA10Nitsche_p(m*numActive[bPsi] + i, j) += mu * weight * bVals[bQ](j,k)*(jac_i*tangent).dot(unormal); break; // p(-p, (\grad psi)_tn)
                        }

                        localMatA10Nitsche_psip(m*numActive[bPsi] + i, j) += -weight * psiVals.row(j).dot(symCurl_i * tangent); // s(psi, psi[p]) undo symmetrization (.,l)

                        localMatA10Nitsche_psip(m*numActive[bPsi] + i, j) += weight * bValsVec_i.dot(bVals[bQ](j,k)*idMatrix*tangent); // c(p,psi)  ibp term p



                    }
            // A01
            for (index_t i = 0; i < numActive[bQ]; i++)
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

                        // Homogenization
                        if(wHomogenPsiq)
                            localMatA01_psiq(i,l*numActive[bPsi] + j) += -weight * (psiVals.row(i)).dot(symCurl_j * tangent); // s(phi, psi[q]) (.,l)

                        if(wPenaltyPsiq)
                            localMatA01_psiq(i,l*numActive[bPsi] + j) += -mu * weight * psiVals.row(i).dot(bValsVec_j); // p(phi, psi[q]) (.,l)
                        //case derivative: MatA01_psiq(i,l*numActive[bPsi] + j) += mu * weight * bVals[bQ](i,k)* (jac_j*tangent).dot(tangent); break; // TODO: prüfen

                    }



            // A00
            for (index_t i = 0; i < numActive[bQ]; i++)
                for (index_t j = 0; j < numActive[bQ]; j++)
                {
                    // Homogenization
                    if(wHomogenPsiq)
                        localMatA00_psiq(i, j) += -weight* ((psiVals.row(i)).dot(bVals[bQ](j,k)* idMatrix * tangent)); // c(p,psi[q]) ibp term p (.,l)

                    if(wPenaltyPsiq)
                        localMatA00_psiq(i, j) += mu* weight* psiVals.row(i).dot(psiVals.row(j)); // p(phi[p], psi[q]) (l,l)
                    //case derivative: localMatA00_psiq(i, j) += mu * weight * bVals[bQ](j,k)* bVals[bQ](i,k); break; // TODO: prüfen

                }



            // psi_q,p term and psi_q,phi term, Nitsche terme (previous elements)
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
                                MatA01_psiq(i,l*numActive[bPsi] + j) += -weight * (assembler.psiq.row(assembler.activesFreeBoundaryMap(i))).dot(symCurl_j * tangent); // s(phi, psi[q]) (.,g)
                            if(wPenaltyPsiq)
                                MatA01_psiq(i,l*numActive[bPsi] + j) += -mu * weight * (assembler.psiq.row(assembler.activesFreeBoundaryMap(i))).dot(bValsVec_j); // p(phi, psi[q]) (.,g)

                            // Nitsche
                            MatA10Nitsche_psip(i,l*numActive[bPsi] + j) += -weight * assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(symCurl_j * tangent); // s(psi, psi[p]) undo symmetrization (.,g)

                            if(phiBcMethod == gsPlateMixedAssembler<T>::nitsche)
                                MatA10Nitsche_psip(i,l*numActive[bPsi] + j) += -mu * weight* assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(bValsVec_j); // p(psi[p], psi) (g,.)


                        }
                    }


            // A00
            for ( index_t i =0; i<assembler.numActivesFreeBoundary; i++)
                for (index_t j = 0; j < numActive[bQ]; j++)
                    if ( assembler.system().rowMapper(0).is_free_index(assembler.activesFreeBoundaryMap(i)) )
                    {
                        // Homogenization
                        if(wHomogenPsiq)
                            MatA00_psiq(i, j) += -weight* assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(bVals[bQ](j,k)* idMatrix * tangent); // c(p,psi[q]) ibp term p (.,g)

                        if(wPenaltyPsiq)
                            MatA00_psiq(i, j) += mu * weight* assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(psiVals.row(j)); // p(phi[p], psi[q]) (l,g)
                    }

            if(wPenaltyPsiq)
            {
                // Homogenization
                for ( index_t j =0; j<assembler.numActivesFreeBoundary; j++)
                    for (index_t i = 0; i < numActive[bQ]; i++)
                        if ( assembler.system().rowMapper(0).is_free_index(assembler.activesFreeBoundaryMap(j)) )
                            MatA00_psip(i, j) += mu * weight* assembler.psiq.row(assembler.activesFreeBoundaryMap(j)).dot(psiVals.row(i)); // p(phi[p], psi[q]) (g,l)


                // Homogenization
                for ( index_t i =0; i<assembler.numActivesFreeBoundary; i++)
                    for ( index_t j =0; j<assembler.numActivesFreeBoundary; j++)
                        if ( assembler.system().rowMapper(0).is_free_index(assembler.activesFreeBoundaryMap(i)) )
                            if ( assembler.system().rowMapper(0).is_free_index(assembler.activesFreeBoundaryMap(j)) )
                                MatA00_psiq_psip(i, j) += mu * weight * assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(assembler.psiq.row(assembler.activesFreeBoundaryMap(j))); // p(phi[p], psi[q]) (g,g)

            }



            // Lagrange multipliers Ln0, Ln1, Lt0, Lt1
            real_t lambdaNCoef_Ln, lambdaTCoef_Ln, lambdaNCoef_Lt, lambdaTCoef_Lt;

            for ( index_t i =0; i<numActive[lambdaBasis]; i++)
            {
                // set coefficients for lambdaN, lambdaT at corners

                // setting interior (not at corner)
                lambdaNCoef_Ln = 1;
                lambdaTCoef_Ln = 0;
                lambdaNCoef_Lt = 0;
                lambdaTCoef_Lt = 1;

                if(actives[lambdaBasis](i,0) == cornerIndex[0]) // start-point edge (counter clock-wise): edge k wrt. corner
                {
                    switch(pde.m_corners[m_patchIndex][corners[0]])
                    {
                    case gsShellMixedPde<T>::sf: lambdaNCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[0]](0,3);
                        lambdaTCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[0]](0,1);
                        lambdaNCoef_Lt = 0; lambdaTCoef_Lt =0; break;
                    case gsShellMixedPde<T>::ff: lambdaNCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[0]](0,3);
                        lambdaTCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[0]](0,1);
                        lambdaNCoef_Lt = pde.m_cornerCouplingCoefs[m_patchIndex][corners[0]](1,3);
                        lambdaTCoef_Lt = pde.m_cornerCouplingCoefs[m_patchIndex][corners[0]](1,1); break;
                    default: break;
                    }
                }
                else if(actives[lambdaBasis](i,0) == cornerIndex[1]) // end-point edge (counter clock-wise): edge k-1 wrt. corner
                {
                    switch(pde.m_corners[m_patchIndex][corners[1]])
                    {
                    case gsShellMixedPde<T>::sf: lambdaNCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[1]](0,2);
                        lambdaTCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[1]](0,0);
                        lambdaNCoef_Lt = 0; lambdaTCoef_Lt =0; break;
                    case gsShellMixedPde<T>::ff: lambdaNCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[1]](0,2);
                        lambdaTCoef_Ln = pde.m_cornerCouplingCoefs[m_patchIndex][corners[1]](0,0);
                        lambdaNCoef_Lt = pde.m_cornerCouplingCoefs[m_patchIndex][corners[1]](1,2);
                        lambdaTCoef_Lt = pde.m_cornerCouplingCoefs[m_patchIndex][corners[1]](1,0); break;
                    default: break;
                    }

                }

                for (index_t j = 0; j < numActive[bPsi]; j++)
                    for(index_t l=0; l<2; l++)
                    {
                        jac_j = gsMatrix<T, 2, 2>::Zero();
                        jac_j.row(l)=physGrad[bPsi].col(j);


                        // Ln1
                        localMatLn1(i,l*numActive[bPsi] + j) += weight * ((jac_j *tangent).dot(unormal) * lambdaNCoef_Ln + (jac_j *tangent).dot(tangent) *lambdaTCoef_Ln ) * bVals[lambdaBasis](i,k);

                        // Lt1
                        localMatLt1(i,l*numActive[bPsi] + j) += weight * ((jac_j *tangent).dot(unormal) * lambdaNCoef_Lt + (jac_j *tangent).dot(tangent) *lambdaTCoef_Lt ) * bVals[lambdaBasis](i,k);
                    }

                for (index_t j = 0; j < numActive[bQ]; j++)
                {
                    // Ln0
                    localMatLn0(i,j) += weight * bVals[bQ](j,k) * lambdaNCoef_Ln * bVals[lambdaBasis](i,k);
                    // Lt0
                    localMatLt0(i,j) += weight * bVals[bQ](j,k) * lambdaNCoef_Lt * bVals[lambdaBasis](i,k); // only non-zero at free coners, lambdaNCoef_Lt = 0 in interior
                }

                // Ltmean
                localMatLtMean(0,i) += weight * lambdaTCoef_Lt * bVals[lambdaBasis](i,k);
                localMatLnMean(0,i) += weight * lambdaTCoef_Ln * bVals[lambdaBasis](i,k);

                // non-homogeneous BC M
                // fmuN rhs Ln1
                //localRhsfmuN(i,0) += weight * ((M_mat * unormal).dot(unormal) * lambdaNCoef_Ln + (M_mat * unormal).dot(tangent) *lambdaTCoef_Ln ) * bVals[lambdaBasis](i,k);

                // fmuT rhs Lt1
                //localRhsfmuT(i,0) += weight * ((M_mat * unormal).dot(unormal) * lambdaNCoef_Lt + (M_mat * unormal).dot(tangent) *lambdaTCoef_Lt ) * bVals[lambdaBasis](i,k);


                // fmuN
                //localRhsfmuN(i,0) += weight * (M_mat*unormal).dot(unormal) * lambdaNCoef_Ln * bVals[bQ](i,k);

                // fmuT
                //int_DivM += weight * DivM.dot(unormal);
                //localRhsfmuT(i,0) += weight * ((M_mat*unormal).dot(tangent) + int_DivM) * lambdaTCoef_Lt * bVals[bQ](i,k);
                //localRhsfmuT(i,0) += weight * (M_mat*unormal).dot(tangent) * lambdaTCoef_Lt * bVals[bQ](i,k);


            }

            // non-homogeneous BC M

            // Nitsche
            // fq

            for (index_t i = 0; i < numActive[bQ]; i++)
            {
                localRhsfq(i, 0) += -mu* weight* psiVals.row(i).dot(psiP_inhomM); // p(phi[p], psi[q]) (l,l)
            }

            for ( index_t i =0; i<assembler.numActivesFreeBoundary; i++)
            {
                if ( assembler.system().rowMapper(0).is_free_index(assembler.activesFreeBoundaryMap(i)) )
                    Rhsfq(i,0) += -mu * weight* assembler.psiq.row(assembler.activesFreeBoundaryMap(i)).dot(psiP_inhomM); // p(phi[p], psi[q]) (l,g)

            }


            // Nitsche
            // fpsi
            for (index_t i = 0; i < numActive[bPsi]; i++)
                for(index_t m = 0; m<2; m++)
                {
                    Curl_i = gsMatrix<T, 2, 2>::Zero();
                    Curl_i(m,0)=physGrad[bPsi](1,i);
                    Curl_i(m,1)=-physGrad[bPsi](0,i);
                    symCurl_i = 0.5*(Curl_i+Curl_i.transpose());

                    bValsVec_i = gsMatrix<T, 2, 1>::Zero();
                    bValsVec_i(m)=bVals[bPsi](i,k);

                    localRhsfpsi(m*numActive[bPsi] + i, 0) += weight * (psiP_inhomM).dot(symCurl_i * tangent); // s(psi, psi[p]) undo symmetrization (.,l)
                    localRhsfpsi(m*numActive[bPsi] + i, 0) += mu * weight* psiP_inhomM.dot(bValsVec_i); // p(psi[p], psi) (l,.)

                }

            // fv
            for ( index_t i =0; i<numActive[bQ]; i++)
            {
                localRhsfv(i,0) += -weight * (dtMnt + DivM.dot(unormal)) * bVals[bQ](i,k);
                //localRhsfv(i,0) += -weight * (DivM.dot(unormal)) * bVals[bQ](i,k);
            }



            if(pBoundary0)
            {
                // p=0 on boundary additional term B1
                for ( index_t i =0; i<numActive[bQ]; i++)
                    for (index_t j = 0; j < numActive[bPsi]; j++)
                        for(index_t l=0; l<2; l++)
                        {
                            jac_j = gsMatrix<T, 2, 2>::Zero();
                            jac_j.row(l)=physGrad[bPsi].col(j);

                            localMatB1(i,l*numActive[bPsi] + j) += weight * (jac_j *tangent).dot(tangent) * physGrad[bQ].col(i).dot(tangent);
                        }
            }

        }

        //std::cout<<"localRhsfv"<<localRhsfv<<std::endl;


        // Flip corners if necessary (s.t. counter-clockwise)
        gsVector<T> element_lowerCorner(element.lowerCorner());
        gsVector<T> element_lowerCorner_tmp(element.lowerCorner());
        gsVector<T> element_upperCorner(element.upperCorner());
        if((sideOrientation(side) * orientation) == -1)
        {
            if(side.direction() == 0)
            {
                element_lowerCorner[1] = 1 - element_upperCorner[1];
                element_upperCorner[1] = 1 - element_lowerCorner_tmp[1];
            }
            else
            {
                element_lowerCorner[0] = 1 - element_upperCorner[0];
                element_upperCorner[0] = 1 - element_lowerCorner_tmp[0];
            }
        }


        // assemble psiq on whole element
        //if(phiBcMethod == gsPlateMixedAssembler<T>::nitsche)
        {
            assembler.applyVisitorPsiq(0, side, element_upperCorner, element);
        }

        // inhomogeneous corner contribution
/*
        for(index_t k= 0; k<2; k++)
        if(element.lowerCorner()(0) <= cornersMatrix(0,k) && element.upperCorner()(0) >= cornersMatrix(0,k) && element.lowerCorner()(1) <= cornersMatrix(1,k) && element.upperCorner()(1) >= cornersMatrix(1,k))
        {
            std::cout<<"at corner:\n"<<cornersMatrix.col(k)<<std::endl;


            // compute physical gradients at k as a Dim x NumActive matrix
            for(size_t i=0; i<basesSize; i++)
                geoEval_corners->transformGradients(k, bGrads[i], physGrad[i]);

            // Compute the outer normal vector on the side
            geoEval_corners->outerNormal(k, side, unormal);

            unormal.normalize();
            gsVector<T,2> tangent;
            tangent(0) = -unormal(1);
            tangent(1) = unormal(0);

            // fpsi
            for (index_t i = 0; i < numActive[bPsi]; i++)
                for(index_t m = 0; m<2; m++)
                {
                    jac_i = gsMatrix<T, 2, 2>::Zero();
                    jac_i.row(m)=physGrad[bPsi].col(i);

                    //if(k == 0)
                    //    localRhsfpsi(m*numActive[bPsi] + i,0) += (jac_i *tangent).dot(tangent) * uVals(0,k);
                    //if(k == 1)
                    //    localRhsfpsi(m*numActive[bPsi] + i,0) += -(jac_i *tangent).dot(tangent) * uVals(0,k);

                }

        }
*/

    }

    inline void localToGlobal(const index_t patchIndex,
                              std::vector<gsMatrix<T> >   & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {

        // Map patch-local DoFs to global DoFs
        actives_vec.resize(7);
        system.mapColIndices(actives[bQ], patchIndex, actives_vec[0], 0);
        system.mapColIndices(actives[bPsi], patchIndex, actives_vec[1], 1);
        system.mapColIndices(actives[bPsi], patchIndex, actives_vec[2], 2);
        system.mapColIndices(actives[bQ], patchIndex, actives_vec[3], 3);

        if(phiBcMethod == gsPlateMixedAssembler<T>::lagrange)
        {
            system.mapColIndices(actives[lambdaBasis], patchIndex, actives_vec[4], 4);
            system.mapColIndices(actives[lambdaBasis], patchIndex, actives_vec[5], 5);
        }

        if(lambdaTMean0)
        {
            actives_vec[6].resize(1,1);
            actives_vec[6] << 0;
            eliminatedDofs[6].resize(1,1);
            eliminatedDofs[6] << 0;
        }

        actives_p = actives_vec[0];
        actives_phi1 = actives_vec[1];
        actives_phi2 = actives_vec[2];

        // build block information
        gsVector<index_t> p_vec(1);
        gsVector<index_t> phi_vec(2);
        gsVector<index_t> w_vec(1);
        gsVector<index_t> lambdaN_vec(1);
        gsVector<index_t> lambdaT_vec(1);
        gsVector<index_t> lambdaMean_vec(1);
        p_vec << 0;
        phi_vec << 1,2;
        w_vec << 3;
        lambdaN_vec << 4;
        lambdaT_vec << 5;
        lambdaMean_vec << 6;


        if(phiBcMethod == gsPlateMixedAssembler<T>::nitsche || phiBcMethod == gsPlateMixedAssembler<T>::nitscheDerivative )
        {
            // A11
            // Nitsche
            system.pushToMatrix(localMatA11Nitsche, actives_vec, eliminatedDofs, phi_vec, phi_vec);

            // A10
            // Nitsche p terms
            system.pushToMatrix(localMatA10Nitsche_p, actives_vec, eliminatedDofs, phi_vec, p_vec);

            // Nitsche psip terms
            system.pushToMatrix(localMatA10Nitsche_psip, actives_vec, eliminatedDofs, phi_vec, p_vec);
            system.pushToMatrix(MatA10Nitsche_psip.block(0,0,assembler.numActivesFreeBoundary,numActive[bPsi]).transpose(), actives_phi1, assembler.activesFreeBoundaryMap, eliminatedDofs[0], 1, 0);
            system.pushToMatrix(MatA10Nitsche_psip.block(0,numActive[bPsi],assembler.numActivesFreeBoundary,numActive[bPsi]).transpose(), actives_phi2, assembler.activesFreeBoundaryMap, eliminatedDofs[0], 2, 0);


            // Homogenization: psi_q terms

            // A00
            system.pushToMatrix(localMatA00_psiq, actives_p, eliminatedDofs[0], 0, 0);
            // A01
            system.pushToMatrix(localMatA01_psiq, actives_vec, eliminatedDofs, p_vec, phi_vec);

            // A00
            system.pushToMatrix(MatA00_psiq, assembler.activesFreeBoundaryMap, actives_p, eliminatedDofs[0], 0, 0);
            system.pushToMatrix(MatA00_psip, actives_p, assembler.activesFreeBoundaryMap, eliminatedDofs[0], 0, 0);
            system.pushToMatrix(MatA00_psiq_psip, assembler.activesFreeBoundaryMap, assembler.activesFreeBoundaryMap, eliminatedDofs[0], 0, 0);

            // A01
            system.pushToMatrix(MatA01_psiq.block(0,0,assembler.numActivesFreeBoundary,numActive[bPsi]), assembler.activesFreeBoundaryMap, actives_phi1, eliminatedDofs[1], 0, 1);
            system.pushToMatrix(MatA01_psiq.block(0,numActive[bPsi],assembler.numActivesFreeBoundary,numActive[bPsi]), assembler.activesFreeBoundaryMap, actives_phi2, eliminatedDofs[2], 0, 2);

            // fq
            system.pushToRhs(localRhsfq, actives_p, 0);
            system.pushToRhs(Rhsfq, assembler.activesFreeBoundaryMap, 0);

            // fpsi
            system.pushToRhs(localRhsfpsi, actives_vec, phi_vec);

        }


        if(phiBcMethod == gsPlateMixedAssembler<T>::lagrange)
        {

            // Ln0
            system.pushToMatrix(localMatLn0, actives_vec, eliminatedDofs, lambdaN_vec, p_vec);
            system.pushToMatrix(localMatLn0.transpose(), actives_vec, eliminatedDofs, p_vec, lambdaN_vec);

            // Ln1
            system.pushToMatrix(localMatLn1, actives_vec, eliminatedDofs, lambdaN_vec, phi_vec);
            system.pushToMatrix(localMatLn1.transpose(), actives_vec, eliminatedDofs, phi_vec, lambdaN_vec);

            // Lt0
            system.pushToMatrix(localMatLt0, actives_vec, eliminatedDofs, lambdaT_vec, p_vec);
            system.pushToMatrix(localMatLt0.transpose(), actives_vec, eliminatedDofs, p_vec, lambdaT_vec);

            // Lt1
            system.pushToMatrix(localMatLt1, actives_vec, eliminatedDofs, lambdaT_vec, phi_vec);
            system.pushToMatrix(localMatLt1.transpose(), actives_vec, eliminatedDofs, phi_vec, lambdaT_vec);

            // fmuN
            system.pushToRhs(localRhsfmuN, actives_vec, lambdaN_vec);

            // fmuT
            system.pushToRhs(localRhsfmuT, actives_vec, lambdaT_vec);

        }

        if(lambdaTMean0)
        {
            // LtMean
            system.pushToMatrix(localMatLtMean, actives_vec, eliminatedDofs, lambdaMean_vec, lambdaT_vec);
            system.pushToMatrix(localMatLtMean.transpose(), actives_vec, eliminatedDofs, lambdaT_vec, lambdaMean_vec);

            // LnMean
            system.pushToMatrix(localMatLnMean, actives_vec, eliminatedDofs, lambdaMean_vec, lambdaN_vec);
            system.pushToMatrix(localMatLnMean.transpose(), actives_vec, eliminatedDofs, lambdaN_vec, lambdaMean_vec);
        }


        // fv
        system.pushToRhs(localRhsfv, actives_vec, w_vec);


        if(pBoundary0)
        {
            // B1, B1^T
            system.pushToMatrix(localMatB1, actives_vec, eliminatedDofs, w_vec, phi_vec);
            system.pushToMatrix(localMatB1.transpose(), actives_vec, eliminatedDofs, phi_vec, w_vec);
        }




    }


protected:

    int phiBcMethod;
    bool wHomogenPsiq;
    bool wPenaltyPsiq;
    bool pBoundary0;
    bool lambdaTMean0;

    bool firstTime;
    memory::shared_ptr<gsGeometryEvaluator<T> > geoEval_corners;
    gsMatrix<T,2,2> cornersMatrix;

    // Additional data for VisitorPlateMixedNewPsiq
    gsPlateMixedAssembler<T>& assembler;
    const gsShellMixedPde<T>& pde;
    index_t m_patchIndex;
    size_t bQ;
    size_t bPsi;
    size_t lambdaBasis;
    size_t basesSize;

    double int_DivM;
    gsVector<T,2> psiP_inhomM;

    // Neumann function
    boxSide side;
    const gsFunction<T> * data_ptr;
    gsMatrix<T> dataVals;
    gsMatrix<T> MVals;
    gsMatrix<T> gradMVals;
    gsMatrix<T> Mnn_yintVals;
    gsMatrix<T> Mnt_yintVals;
    //gsMatrix<T> uVals;
    real_t penalty;

    // Basis values
    std::vector<std::vector<gsMatrix<T> > > basisData;
    std::vector<std::vector<gsMatrix<T> > > basisData_corners;
    std::vector<gsMatrix<index_t> > actives;
    std::vector<gsMatrix<T> >		 physGrad;
    std::vector<index_t >            numActive;
    std::vector<boxCorner> corners;
    std::vector<index_t> cornerIndex;


    std::vector<gsMatrix<index_t> > actives_vec;
    gsMatrix<index_t> actives_p;
    gsMatrix<index_t> actives_phi1;
    gsMatrix<index_t> actives_phi2;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> pVals;
    gsMatrix<T> m_quNodes;

    // Local matrix and rhs
    gsMatrix<T> localMatwN;

    //Nitsche
    gsMatrix<T> localMatA11Nitsche;
    gsMatrix<T> localMatA10Nitsche_p;
    gsMatrix<T> localMatA10Nitsche_psip;
    gsMatrix<T> MatA10Nitsche_psip;

    // psiq
    gsMatrix<T> localMatA00_psiq;
    gsMatrix<T> localMatA01_psiq;
    gsMatrix<T> MatA00_psiq;
    gsMatrix<T> MatA00_psip;
    gsMatrix<T> MatA00_psiq_psip;
    gsMatrix<T> MatA01_psiq;

    gsMatrix<T> localMatLn0;
    gsMatrix<T> localMatLn1;
    gsMatrix<T> localMatLt0;
    gsMatrix<T> localMatLt1;
    gsMatrix<T> localMatLnMean;
    gsMatrix<T> localMatLtMean;

    gsMatrix<T> localMatB1;

    // non-homogeneous BC M
    gsMatrix<T> localRhsfq;
    gsMatrix<T> localRhsfpsi;
    gsMatrix<T> localRhsfv;
    gsMatrix<T> localRhsfmuN;
    gsMatrix<T> localRhsfmuT;
    gsMatrix<T> Rhsfq;

    gsMapData<T> md;
    int orientation;
};


} // namespace gismo
