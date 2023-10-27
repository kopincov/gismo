/** @file gsVisitorBiharmonic.h

    @brief Visitor for a simple Biharmonic equation with G1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn & P. Weinmüller
*/

#pragma once

#include <gsMatrix/gsMatrix.h>
#include <gsMatrix/gsVector.h>

namespace gismo
{

/** \brief Visitor for the biharmonic equation.
 *
 * Assembles the bilinear terms
 * \f[ (\Delta u,\Delta v)_\Omega \text{ and } (f,v)_\Omega \f]
 * For \f[ u = g \quad on \quad \partial \Omega \f],
 *
 */

template <class T>
class gsVisitorBiharmonic
{
public:

    gsVisitorBiharmonic(const gsPde<T> & pde)
    { 
        rhs_ptr = static_cast<const gsG1BiharmonicPde<T>&>(pde).rhs() ;

        g1_basis_ptr = static_cast<const gsG1BiharmonicPde<T>&>(pde).G1Basis() ;

    }

    /** \brief Constructor for gsVisitorBiharmonic.
     *
     * \param[in] rhs Given right-hand-side function/source term that, for
     */
    gsVisitorBiharmonic(const gsFunction<T> & rhs) :
        rhs_ptr(&rhs)
    {
        GISMO_ASSERT( rhs.targetDim() == 1 ,"Not yet tested for multiple right-hand-sides");
    }

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T>    & rule)
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// NB!

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
    }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options, 
                    gsQuadRule<T>    & rule)
    {
        m_patchIndex = patchIndex;

        // Setup Quadrature
        rule = gsQuadrature::get(basis, options);

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>       & basis, // to do: more unknowns
                         const gsGeometry<T>    & geo, // patch
                         gsMatrix<T>            & quNodes)
    {
        index_t n = g1_basis_ptr->at(m_patchIndex).nPatches();

        md.points = quNodes;
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);

        // TODO nur wenn es nötig ist, hinzufügen!!
        gsVector<unsigned> linspace(n);
        //linspace << 25,26,27,28,29,30,31,32,33,34;
        linspace.setLinSpaced(n,basis.size(),basis.size()+n); // basis.size() = max number TODO

        actives.conservativeResize(actives.rows()+n,actives.cols());
        actives.bottomRows(n) = linspace;

        numActive = actives.rows();

        //gsInfo << "qunodes 2: " << g1_basis_ptr->at(1).patch(1).deriv(md.points) << "\n";

        // Evaluate basis functions on element
        basis.evalAllDers_into(md.points, 2, basisData);

        // Evaluate g1 basis functions
        gsMatrix<T> eval_temp(n,md.points.cols());
        gsMatrix<T> deriv_temp(2*n,md.points.cols());
        gsMatrix<T> deriv2_temp(3*n,md.points.cols());


        if (m_patchIndex == 0)
            gsWriteParaview(g1_basis_ptr->at(m_patchIndex).patch(6),"Basis_Visitor",5000);
        if (m_patchIndex == 1)
            gsWriteParaview(g1_basis_ptr->at(m_patchIndex).patch(6),"Basis_Visitor2",5000);

        for (index_t i = 0; i < n; i++)
        {

            eval_temp.row(i) = g1_basis_ptr->at(m_patchIndex).patch(i).eval(md.points);
            deriv_temp.block(2*i,0,2,md.points.cols()) = g1_basis_ptr->at(m_patchIndex).patch(i).deriv(md.points);
            deriv2_temp.block(3*i,0,3,md.points.cols()) = g1_basis_ptr->at(m_patchIndex).patch(i).deriv2(md.points);


        }

        basisData_g1.resize(3);
        basisData_g1[0] = eval_temp;
        basisData_g1[1] = deriv_temp;
        basisData_g1[2] = deriv2_temp;

        gsMatrix<T> temp1(basisData.at(0).rows() + eval_temp.rows(),eval_temp.cols());
        temp1.block(0,0,basisData.at(0).rows(),eval_temp.cols()) = basisData[0];
        temp1.block(basisData.at(0).rows(),0,eval_temp.rows(),eval_temp.cols()) = eval_temp;

        gsMatrix<T> temp2(basisData.at(1).rows() + deriv_temp.rows(),deriv_temp.cols());
        temp2.block(0,0,basisData.at(1).rows(),deriv_temp.cols()) = basisData[1];
        temp2.block(basisData.at(1).rows(),0,deriv_temp.rows(),deriv_temp.cols()) = deriv_temp;

        gsMatrix<T> temp3(basisData.at(2).rows() + deriv2_temp.rows(),deriv2_temp.cols());
        temp3.block(0,0,basisData.at(2).rows(),deriv2_temp.cols()) = basisData[2];
        temp3.block(basisData.at(2).rows(),0,deriv2_temp.rows(),deriv2_temp.cols()) = deriv2_temp;

        basisData.clear();
        basisData.resize(3);
        basisData[0] = temp1;
        basisData[1] = temp2;
        basisData[2] = temp3;

        //gsInfo << "md points " << geo. <<"\n";

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo.computeMap(md);

        gsInfo << "basis function: " << m_patchIndex << " mit " << deriv_temp.at(0) << "\n";

        // Evaluate right-hand side at the geometry points
        rhs_ptr->eval_into(md.values[0], rhsVals); // Dim: 1 X NumPts

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
    }


    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData[0];
        gsMatrix<T> & basisGrads = basisData[1];
        gsMatrix<T> & basis2ndDerivs = basisData[2];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);

            gsMatrix<T> temp_matrix;
            transformLaplaceHgrad(md, k, basisData_g1[1], basisData_g1[2], temp_matrix);

            // Compute physical laplacian at k as a 1 x numActive matrix
            transformLaplaceHgrad(md, k, basisGrads, basis2ndDerivs, physBasisLaplace);

            //gsInfo << "physBasisLaplace" << physBasisLaplace << "\n";
            transformGradients(md, k, basisGrads, physGrads);

            // (\Delta u, \Delta v)
            localMat.noalias() += weight * (physBasisLaplace.transpose() * physBasisLaplace);

            localRhs.noalias() += weight * ( basisVals.col(k) * rhsVals.col(k).transpose() ) ;
        }
    }

    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        //gsInfo << "actives local: " << actives << "\n";

        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, eliminatedDofs[0], 0, 0);
    }

    /*
    inline void localToGlobal(const gsDofMapper     & mapper,
                              const gsMatrix<T>     & eliminatedDofs,
                              const index_t patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
        // Local Dofs to global dofs
        mapper.localToGlobal(actives, patchIndex, actives);
        //const int numActive = actives.rows();

        for (index_t i=0; i < numActive; ++i)
        {
            const int ii = actives(i);
            if ( mapper.is_free_index(ii) )
            {
                rhsMatrix.row(ii) += localRhs.row(i);

                for (index_t j=0; j < numActive; ++j)
                {
                    const int jj = actives(j);
                    if ( mapper.is_free_index(jj) )
                    {
                        sysMatrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else
                    {
                        rhsMatrix.row(ii).noalias() -= localMat(i, j) *
                            eliminatedDofs.row( mapper.global_to_bindex(jj) );
                    }
                }
            }
        }
    }
    */


protected:
    // Right hand side
    const gsFunction<T> * rhs_ptr;

    const std::vector< gsMultiPatch<>> * g1_basis_ptr;

protected:
    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T>        physBasisLaplace;
    gsMatrix<T>        physGrads;
    gsMatrix<index_t> actives;
    index_t numActive;

    // For the G1 Basis:
    std::vector<gsMatrix<T> > basisData_g1;

    index_t m_patchIndex;


protected:
    // Local values of the right hand side
    gsMatrix<T> rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    gsMapData<T> md;
};


} // namespace gismo

