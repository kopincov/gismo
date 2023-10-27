/** @file gsVisitorPlateMixedClamped.h

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
class gsVisitorPlateMixedClamped
{
public:

    gsVisitorPlateMixedClamped(const gsPde<T> & pde, const boundary_condition<T> & s, gsPlateMixedAssembler<T>& assembler)
        : assembler(assembler), pde(dynamic_cast<const gsShellMixedPde<T>& >(pde)), side(s.side()), gradData_ptr(s.function().get())
    { }


    void initialize(const gsBasisRefs<T> & bases,
                    const index_t          patchIndex,
                    const gsOptionList   & options,
                          gsQuadRule<T>  & rule)
    {
        // Set options
        bQ = assembler.system().colBasis(0);
        bPsi = assembler.system().colBasis(1);
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

        // Compute geometry related values
        geo.computeMap(md);

        // Evaluate the Neumann data
        gradData_ptr->eval_into(md.values[0], gradVals);

        // Initialize local matrix/rhs
        localRhsfq.setZero(numActive[bQ], 1);
        localRhsfpsi.setZero(2 * numActive[bPsi], 1);
    }

    inline void assemble(gsDomainIterator<T> & ,
                         const gsVector<T>   & quWeights)
    {
        std::vector<gsMatrix<T> > bVals(basesSize); // ursprünglichen version reference
        for (size_t i = 0; i < basesSize; i++)
            bVals[i] = basisData[i][0];

        std::vector<gsMatrix<T> > bGrads(basesSize);
        for (size_t i = 0; i < basesSize; i++)
            bGrads[i] = basisData[i][1];


        //gsMatrix<T, 2, 2> idMatrix = gsMatrix<T, 2, 2>::Identity();
        gsMatrix<T, 2, 2> jac_i;

        for (index_t k = 0; k < quWeights.rows(); ++k)
        {
            // compute physical gradients at k as a Dim x NumActive matrix
            for (size_t i = 0; i < basesSize; i++)
                transformGradients(md, k, bGrads[i], physGrad[i]);

            // Compute the outer normal vector on the side
            outerNormal(md, k, side, unormal);

            // Multiply quadrature weight by the measure of normal
            const T weight = quWeights[k] * unormal.norm();

            unormal.normalize();
            gsVector<T, 2> tangent;
            tangent(0) = -unormal(1);
            tangent(1) = unormal(0);

            // fq
            /*
            for (index_t i = 0; i < numActive[bQ]; i++)
            {
                localRhsfq(i,0) += -weight * gradVals.col(k).dot(unormal) * bVals[bQ](i,k);
            }
            */

            // fpsi
            for (index_t i = 0; i < numActive[bPsi]; i++)
            {
                for (index_t m = 0; m < 2; m++)
                {
                    jac_i = gsMatrix<T, 2, 2>::Zero();
                    jac_i.row(m) = physGrad[bPsi].col(i);

                    localRhsfpsi(m * numActive[bPsi] + i, 0) += -weight
                        * ((jac_i * tangent).dot(unormal) * gradVals.col(k).dot(unormal)
                            + (jac_i * tangent).dot(tangent) * gradVals.col(k).dot(tangent));
                }
            }
        }
    }

    inline void localToGlobal(const index_t                     patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                                    gsSparseSystem<T>         & system)
    {
        // Map patch-local DoFs to global DoFs
        actives_vec.resize(5);

        system.mapColIndices(actives[bQ], patchIndex, actives_vec[0], 0);
        system.mapColIndices(actives[bPsi], patchIndex, actives_vec[1], 1);
        system.mapColIndices(actives[bPsi], patchIndex, actives_vec[2], 2);
        system.mapColIndices(actives[bQ], patchIndex, actives_vec[3], 3);

        // build block information
        gsVector<index_t> p_vec(1);
        gsVector<index_t> phi_vec(2);
        gsVector<index_t> u_vec(1);
        p_vec << 0;
        phi_vec << 1, 2;
        u_vec << 3;

        // fq
        //system.pushToRhs(localRhsfq, actives_vec, p_vec);

        // fpsi
        system.pushToRhs(localRhsfpsi, actives_vec, phi_vec);
    }

private:
    int phiBcMethod;
    bool wHomogenPsiq;
    bool wPenaltyPsiq;

    // Additional data for VisitorPlateMixedNewPsiq
    gsPlateMixedAssembler<T>& assembler;
    const gsShellMixedPde<T>& pde;
    index_t m_patchIndex;
    size_t bQ;
    size_t bPsi;
    size_t basesSize;

    // Neumann function
    boxSide side;
    const gsFunction<T> * gradData_ptr;
    gsMatrix<T> gradVals;
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
    // Homogenization w
    gsMatrix<T> localRhsfq;
    gsMatrix<T> localRhsfpsi;

    gsMapData<T> md;
};

} // namespace gismo
