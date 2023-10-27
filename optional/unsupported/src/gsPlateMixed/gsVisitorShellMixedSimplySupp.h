/** @file gsVisitorShellMixedSimplySupp.h

    @brief Boundary element visitor for the contribution of simplysupp boundaries

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/

#pragma once
#include <gsPlateMixed/gsShellMixedAssembler.h>

namespace gismo
{

template <class T>
class gsVisitorShellMixedSimplySupp
{
public:

    gsVisitorShellMixedSimplySupp(const gsPde<T> & pde, const boundary_condition<T> & s, const gsShellMixedAssembler<T>& assembler)
        : m_assembler(assembler), pde(dynamic_cast<const gsShellMixedPde<T>& >(pde)), side(s.side())
    { }


    void initialize(const gsBasisRefs<T> & bases,
                    const index_t patchIndex,
                    const gsOptionList& options,
                    gsQuadRule<T>    & rule)
    {
        // Store basis information
        pBasis = m_assembler.system().colBasis(0);
        phiBasis = m_assembler.system().colBasis(1);
        uBasis = m_assembler.system().colBasis(6);
        lambdaBasis = m_assembler.system().colBasis(9);
        basesSize = bases.size();

        basisData.resize(basesSize);
        actives.resize(basesSize);
        numActive.resize(basesSize);
        bGrads_k.resize(basesSize);

        // Setup Quadrature (harmless slicing occurs)
        rule = gsGaussRule<T>(bases[0], options.getReal("quA"), options.getInt("quB"), side.direction());

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_OUTER_NORMAL | NEED_GRAD_TRANSFORM;

        m_patchIndex = patchIndex;
    }


    // Evaluate on element.
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

        /*
        std::cout<<"corners"<<std::endl;
        std::cout<<corners[0]<<std::endl;
        std::cout<<corners[1]<<std::endl;
        */

        // Compute geometry related quantities
        geo.computeMap(md);

        // Initialize local matrix/rhs
        localMatLn1.setZero(numActive[lambdaBasis], 2 * numActive[phiBasis]);
    }

    inline void assemble(gsDomainIterator<T> & /*element*/,
                         gsVector<T> const   & quWeights)
    {
        std::vector<gsMatrix<T> > bVals(basesSize);
        std::vector<gsMatrix<T> > bGrads(basesSize);

        for(size_t i=0; i<basesSize; i++)
        {
            bVals[i] = basisData[i][0];
            bGrads[i] = basisData[i][1];
        }

        gsMatrix<T,2,2> jac_j;


        for (index_t k = 0; k < quWeights.rows(); ++k)
        {

            // Compute the outer normal vector on the side // TODO: normal on parameter space
            //outerNormal(md, k, side, unormal);
            //unormal.normalize();

            // normal on parameter space
            gsVector<T,2> tangent;
            tangent.setOnes(2);
            tangent(side.direction()) = 0;
            tangent *=sideOrientation(side);

            unormal.resize(2);
            unormal(0) = tangent(1);
            unormal(1) = -tangent(0);

            //gsInfo<<"side \n"<<side<<"\n";
            //gsInfo<<"unormal \n"<<unormal<<"\n";

            for(size_t i=0; i<basesSize; i++)
            {
                bGrads_k[i] = bGrads[i].col(k);
                bGrads_k[i].resize(2, numActive[i]);
            }

            // Multiply quadrature weight by the measure of normal
            const T weight = quWeights[k];

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

                for (index_t j = 0; j < numActive[phiBasis]; j++)
                    for(index_t l=0; l<2; l++)
                    {
                        jac_j = gsMatrix<T, 2, 2>::Zero();
                        jac_j.row(l)=bGrads_k[phiBasis].col(j);

                        // Ln1
                        localMatLn1(i,l*numActive[phiBasis] + j) += weight * ((jac_j *tangent).dot(unormal) * lambdaNCoef_Ln) * bVals[lambdaBasis](i,k); //lambdaTCoef_Ln = 0 because BC on ss

                    }

            }

            /* old

            // L lagrange multiplier (grad phi)_tn
            for ( index_t i =0; i<numActive; i++)
                for (index_t j = 0; j < numActive; j++)
                    for(index_t l=0; l<2; l++)
                    {
                        jac_j = gsMatrix<T, 2, 2>::Zero();
                        jac_j.row(l)=bGrads_k.col(j);

                        if(actives(i,0) == cornerIndex[0])
                        {
                            if(pde.m_corners[m_patchIndex][corners[0]] == gsShellMixedPde<T>::sf)
                                localMatLn1(i,l*numActive + j) += weight * (pde.m_cornerCouplingCoefs[m_patchIndex][corners[0]](3,0)*(jac_j *tangent).dot(unormal)) * bVals(i,k);
                            if(pde.m_corners[m_patchIndex][corners[0]] == gsShellMixedPde<T>::none)
                                localMatLn1(i,l*numActive + j) += weight * (jac_j *tangent).dot(unormal) * bVals(i,k);

                        }
                        else
                        {
                            if(actives(i,0) == cornerIndex[1])
                            {
                                if(pde.m_corners[m_patchIndex][corners[1]] == gsShellMixedPde<T>::sf)
                                    localMatLn1(i,l*numActive + j) += weight * (pde.m_cornerCouplingCoefs[m_patchIndex][corners[1]](2,0)*(jac_j *tangent).dot(unormal)) * bVals(i,k);
                                if(pde.m_corners[m_patchIndex][corners[1]] == gsShellMixedPde<T>::none)
                                    localMatLn1(i,l*numActive + j) += weight * (jac_j *tangent).dot(unormal) * bVals(i,k);

                            }
                            else
                                localMatLn1(i,l*numActive + j) += weight * (jac_j *tangent).dot(unormal) * bVals(i,k);
                        }

                    }
                    */

        }

    }

    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> >   & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        actives_vec.resize(10);
        system.mapColIndices(actives[phiBasis], patchIndex, actives_vec[1], 1);
        system.mapColIndices(actives[phiBasis], patchIndex, actives_vec[2], 2);

        system.mapColIndices(actives[lambdaBasis], patchIndex, actives_vec[9], 9);

        // build block information
        gsVector<index_t> phi_vec(2);
        gsVector<index_t> lambdaN_vec(1);
        phi_vec << 1,2;
        lambdaN_vec << 9;

        // Ln1
        system.pushToMatrix(localMatLn1, actives_vec, eliminatedDofs, lambdaN_vec, phi_vec);
        system.pushToMatrix(localMatLn1.transpose(), actives_vec, eliminatedDofs, phi_vec, lambdaN_vec);

    }


protected:

    // Assembler
    const gsShellMixedAssembler<T>& m_assembler;

    // Bases
    size_t pBasis;
    size_t phiBasis;
    size_t uBasis;
    size_t lambdaBasis;
    size_t basesSize;

    // Pde
    const gsShellMixedPde<T>& pde;

    index_t m_patchIndex;

    // Neumann function
    const gsFunction<T> * pData_ptr;
    boxSide side;
    real_t penalty;

    // Basis values
    std::vector<std::vector<gsMatrix<T> > > basisData;
    std::vector<gsMatrix<index_t> > actives;
    std::vector<index_t >            numActive;
    std::vector<gsMatrix<T> > bGrads_k;
    std::vector<gsMatrix<index_t> > actives_vec;

    std::vector<boxCorner> corners;
    std::vector<index_t> cornerIndex;

    // Normal and Neumann values
    gsVector<T> unormal;

    // Local matrix
    gsMatrix<T> localMatLn1;

    gsMapData<T> md;
};


} // namespace gismo
