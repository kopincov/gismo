/** @file gsVisitorSpaceTimeSlice.h

    @brief A DG interface visitor for the Parabolic problem .

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#pragma once

#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{


template <class T>
class gsVisitorSpaceTimeInterface
{
public:
    gsVisitorSpaceTimeInterface(const gsPde<T> & pde) :
        d(0)
    { }

    void initialize(const gsBasis<T>        & basis1,
                    const gsBasis<T>        & basis2,
                    const boundaryInterface & bi,
                    const gsOptionList      & options,
                          gsQuadRule<T>     & rule)
    {
        d = basis1.dim();
        side1 = bi.first();
        side2 = bi.second();

        bActives1 = basis1.boundary(side1);

        bActives2 = basis1.boundary(side2);

        const T quA = options.getReal("quA") + 1;
        const index_t quB = options.getInt("quB");

        // Setup Quadrature
        rule = gsGaussRule<T>(basis1, quA, quB, side1.direction());


        // Set Geometry evaluation flags
        md1.flags = md2.flags = NEED_VALUE | NEED_MEASURE | NEED_OUTER_NORMAL;
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>    & B1, // to do: more unknowns
                         const gsGeometry<T> & geo1,
                         const gsBasis<T>    & B2, // to do: more unknowns
                         const gsGeometry<T> & geo2,
                               gsMatrix<T>   & quNodes1,
                               gsMatrix<T>   & quNodes2)
    {
        md1.points = quNodes1;
        md2.points = quNodes2;

        changed_ordering = false;
        // Compute the active basis functions
        B1.active_into(md1.points.col(0), actives1);
        B2.active_into(md2.points.col(0), actives2);

        //actives1 = intersectMatrix(bActives1,actives1);
        //  actives2 = intersectMatrix(bActives2,actives2);


        // Evaluate basis functions and their first derivatives
        B1.evalAllDers_into(md1.points, 0, basisData1);
        B2.evalAllDers_into(md2.points, 0, basisData2);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo1.computeMap(md1);
        geo2.computeMap(md2);

        outerNormal(md1, 0, side1, unormal);
        if (unormal.tail(1)(0, 0) < 0)
        {
            actives1.swap(actives2);
            basisData1[0].swap(basisData2[0]);
            changed_ordering = true;
        }

        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        // Initialize local matrices
        E_mixed.setZero(numActive1, numActive2);
        E_top.setZero(numActive2, numActive2);
    }

    // assemble on element
    inline void assemble(gsDomainIterator<T> & element1,
                         gsDomainIterator<T> & element2,
                         gsVector<T>         & quWeights)
    {
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            outerNormal(md1, k, side1, unormal);

            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T weight = quWeights[k] * unormal.norm();


            // Take blocks of values and derivatives of basis functions
            const typename gsMatrix<T>::Block val1 = basisData1[0].block(0, k, numActive1, 1);
            const typename gsMatrix<T>::Block val2 = basisData2[0].block(0, k, numActive2, 1);

            // Compute element matrices
            E_mixed.noalias() += weight * (val1 * val2.transpose());
            E_top.noalias() += weight * (val2 * val2.transpose());
        }
    }

    inline void localToGlobal(const index_t                     patchIndex1,
                              const index_t                     patchIndex2,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                                    gsSparseSystem<T>         & system)
    {
        // Map patch-local DoFs to global DoFs
        if (!changed_ordering)
        {
            system.mapColIndices(actives1, patchIndex1, actives1);
            system.mapColIndices(actives2, patchIndex2, actives2);
        }
        else
        {
            system.mapColIndices(actives1, patchIndex2, actives1);
            system.mapColIndices(actives2, patchIndex1, actives2);
        }

        gsMatrix<T> localRhs1, localRhs2;
        localRhs1.setZero(actives1.rows(), 1);
        //  localRhs2.setZero(actives2.rows(),1);

        system.push(-E_mixed.transpose(), localRhs1, actives2, actives1, eliminatedDofs.front(), 0, 0);
        system.push(E_top, localRhs1, actives2, actives2, eliminatedDofs.front(), 0, 0);
    }

private:
    gsMatrix<unsigned>::uPtr intersectMatrix(gsMatrix<unsigned> & matBig, gsMatrix<unsigned> & matSmall)
    {
        gsMatrix<unsigned> *mat = new gsMatrix<unsigned>(matBig.rows(), 1);
        int n = 0;
        int j_start = 0;
        for (int i = 0; i < matSmall.rows(); ++i)
        {
            for (int j = j_start; j < matBig.rows(); ++j)
            {
                if (matSmall(i, 0) > matBig(j, 0))
                {
                    j_start = j;
                    continue;
                }
                else if (matSmall(i, 0) == matBig(j, 0))
                {
                    (*mat)(n, 0) = matSmall(i, 0);
                    n++;
                    j_start = j + 1;
                    break;
                }
                else if (matSmall(i, 0) < matBig(j, 0))
                    break;

            }
        }
        mat->conservativeResize(n, 1);
        return memory::make_unique(mat);
    }

private:

    // Penalty constant
    T penalty;

    // Side
    boxSide side1, side2;

    // dimension of the problem
    unsigned d;

private:
    gsMatrix<index_t> bActives1, bActives2;
    gsMapData<T> md1, md2;

    // Basis values etc
    std::vector<gsMatrix<T> > basisData1, basisData2;
    gsMatrix<index_t> actives1  , actives2;

    gsVector<T> unormal;

    bool checked_ordering;
    bool changed_ordering;

    // Auxiliary element matrices
    gsMatrix<T> E_mixed, E_top;
};

template <class T>
class gsVisitorSpaceTimeInitial
{
public:
    gsVisitorSpaceTimeInitial(const gsPde<T> & pde, const boundary_condition<T> & s) :
        initialData(s.function().get()), side(s.side()), d(0)
    { }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        d = basis.dim();

        const T quA = options.getReal("quA");
        const index_t quB = options.getInt("quB");
        rule = gsGaussRule<T>(basis, quA, quB, side.direction());

        bActives = basis.boundary(side);

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_OUTER_NORMAL;
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>    & basis, // to do: more unknowns
                         const gsGeometry<T> & geo,
                               gsMatrix<T>   & quNodes)
    {
        md.points = quNodes;
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(md.points.col(0), actives);
        //   actives = intersectMatrix(bActives,actives);

        const index_t numActive = actives.rows();

        // Evaluate basis values on element
        basis.evalAllDers_into(md.points, 0, basisData);

        // Compute geometry related values
        geo.computeMap(md);

        // Evaluate the Dirichlet data
        initialData->eval_into(md.values[0], initData);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, initialData->targetDim());
    }

    // assemble on element
    inline void assemble(gsDomainIterator<T> & element,
                         gsVector<T>         & quWeights)
    {
        //const unsigned d = element1.dim();
        const index_t numActive = actives.rows();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            outerNormal(md, k, side, unormal);

            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T weight = quWeights[k] * unormal.norm();

            // Take blocks of values and derivatives of basis functions
            const typename gsMatrix<T>::Block val = basisData[0].block(0, k, numActive, 1);

            // Compute element matrices
            localRhs.noalias() += weight * (val * initData.col(k).transpose());

            localMat.noalias() += weight * (val * val.transpose());

        }
    }

    inline void localToGlobal(const index_t                     patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>               & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, actives, eliminatedDofs.front(), 0, 0);
    }

private:
    gsMatrix<index_t>::uPtr intersectMatrix(gsMatrix<index_t> & matBig, gsMatrix<index_t> & matSmall )
    {
        gsMatrix<index_t> *mat = new gsMatrix<index_t>(matBig.rows(), 1);
        int n = 0;
        int j_start = 0;
        for (int i = 0; i < matSmall.rows(); ++i)
        {
            for (int j = j_start; j < matBig.rows(); ++j)
            {
                if (matSmall(i, 0) > matBig(j, 0))
                {
                    j_start = j;
                    continue;
                }
                else if (matSmall(i, 0) == matBig(j, 0))
                {
                    (*mat)(n, 0) = matSmall(i, 0);
                    n++;
                    j_start = j + 1;
                    break;
                }
                else if (matSmall(i, 0) < matBig(j, 0))
                    break;

            }
        }
        mat->conservativeResize(n, 1);
        return memory::make_unique(mat);
    }

    const gsFunction<T> * initialData;

    // Side
    boxSide side;

    // dimension of the problem
    unsigned d;

    gsMatrix<index_t> bActives;

    // Basis values etc
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<index_t> actives;

    // Normal and Initial values
    gsVector<T> unormal;
    gsMatrix<T> initData;

    // Auxiliary element matrices
    gsMatrix<T> localMat, localRhs;
    gsMapData<T> md;
};

} // namespace gismo

