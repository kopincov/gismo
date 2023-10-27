/** @file gsVisitorShellMixedLineLoad.h

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
class gsVisitorShellMixedLineLoad
{
public:

    gsVisitorShellMixedLineLoad(const gsPde<T> & pde, const boundary_condition<T> & s, const gsShellMixedAssembler<T>& assembler)
        : m_assembler(assembler), pde(dynamic_cast<const gsShellMixedPde<T>& >(pde)), side(s.side()), lineLoad_ptr(s.function().get())
    { }


    void initialize(const gsBasisRefs<T> & bases,
                    const index_t          patchIndex,
                    const gsOptionList   & options,
                          gsQuadRule<T>  & rule)
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

        // Compute geometry related quantities
        geo.computeMap(md);

        // Evaluate right-hand side at the geometry points
        lineLoad_ptr->eval_into(md.values[0], forceVals);

        // Initialize local matrix/rhs
        localRhsfv.setZero(3 * numActive[uBasis], 1);
    }

    inline void assemble(      gsDomainIterator<T> & element,
                         const gsVector<T>         & quWeights)
    {
        std::vector<gsMatrix<T> > bVals(basesSize);
        bVals[uBasis] = basisData[uBasis][0];

        gsMatrix<T, 3, 3> F;
        gsVector<T> normalVec;


        for (index_t k = 0; k < quWeights.rows(); ++k)
        {

            // Compute the outer normal vector on the side // TODO: normal on parameter space
            //outerNormal(md, k, side, unormal);
            //const T measure = unormal.norm();

            //gsInfo<<"side \n"<<side<<"\n";
            //gsInfo<<"unormal \n"<<unormal<<"\n";

            // Multiply quadrature weight by the measure of normal
            const T weight = quWeights[k];

            // compute contravariant components of forceVals in the covariant basis
            // F
            normal(md, k, normalVec);
            normalVec.normalize();
            F.leftCols(2) = md.jacobian(k);
            F.col(2) = normalVec;

            //gsInfo<<"normalVec \n"<<normalVec<<"\n";

            forceVals.col(k) = F.inverse() * forceVals.col(k);

            // Local rhs vector contribution
            for (index_t j = 0; j != 3; ++j)
                localRhsfv.middleRows(j * numActive[uBasis], numActive[uBasis]).noalias() +=
                    -weight * forceVals(j, k) * bVals[uBasis].col(k);

        }
    }

    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> >   & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        //const T thickness = pde.m_thickness;

        // Map patch-local DoFs to global DoFs
        actives_vec.resize(9);

        system.mapColIndices(actives[uBasis], patchIndex, actives_vec[6], 6);
        system.mapColIndices(actives[uBasis], patchIndex, actives_vec[7], 7);
        system.mapColIndices(actives[uBasis], patchIndex, actives_vec[8], 8);


        // build block information
        gsVector<index_t> u_vec(3);
        u_vec << 6,7,8;

        // Rhsfv
        system.pushToRhs(localRhsfv, actives_vec, u_vec);
        //system.pushToRhs(1/thickness * localRhsfv, actives_vec, u_vec);

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

    // lineLoad function
    boxSide side;
    const gsFunction<T> * lineLoad_ptr;
    gsMatrix<T> forceVals;

    // Basis values
    std::vector<std::vector<gsMatrix<T> > > basisData;
    std::vector<gsMatrix<index_t> > actives;
    std::vector<index_t >            numActive;
    std::vector<gsMatrix<index_t> > actives_vec;

    // Normal and Neumann values
    gsVector<T> unormal;

    // Local matrix
    gsMatrix<T> localRhsfv;
    gsMapData<T> md;
};

} // namespace gismo
