/** @file gsStokesIterativeSolver.hpp

    @brief Solver class for the Stokes problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#pragma once
#include "gsStokesIterativeSolver.h"
namespace gismo {


template<class T>
void gsStokesIterativeSolver<T>::initWithStokesAssembler()
{
    GISMO_ASSERT(m_stokesAssembler->fieldDim(0) == m_stokesAssembler->patches().parDim(),
                 "Template parameter d should be the same as parDim og the geometry");

    m_matrix = m_stokesAssembler->systemMatrix();
    m_rhs = m_stokesAssembler->systemRhs();

    m_dofMapper = m_stokesAssembler->dofMapperVector();
    m_idofs = m_matrix.cols();
    m_sz_A = m_dofMapper[0]->freeSize();
    m_sz_B = m_dofMapper[1]->freeSize();
    m_tarDim = m_stokesAssembler->fieldDim(0);
    m_nu = m_stokesAssembler->getnu();
    m_patches = m_stokesAssembler->patches();
    m_bcondition = m_stokesAssembler->boundaryConditions(0);
    m_dirStrategy = m_stokesAssembler->getDirStrategy();
    m_geoTrans = m_stokesAssembler->getGeoTrans();

    for (index_t k = 0; k < m_tarDim ; ++k)
    {
        m_bases_u.push_back(m_stokesAssembler->basis(k));
    }
    m_basis_p = m_stokesAssembler->basis(m_tarDim);

}

template<class T>
template<unsigned d>
gsMultiGridOp<>::Ptr gsStokesIterativeSolver<T>::initMultiGrid(index_t numPreSmo, index_t numPostSmo, index_t coarseLvlMG)
{
    GISMO_ASSERT(m_patches.nPatches() == 1,"Only implemented for single Patch!");

    //Initilase bases for each level (and component for velocity)
    //use to find transfere matricis
    std::vector< std::vector< gsTensorBSplineBasis<d,real_t>  *> >bases_u(m_tarDim);

    //Storage for the removed knots in the coarse basis
    std::vector< std::vector< std::vector< std::vector<real_t> > > > removedKnots_u(m_tarDim);

    //Push back the finist basis
    for ( index_t k = 0; k< m_tarDim; ++k)
        bases_u[k].push_back( memory::convert_ptr<gsTensorBSplineBasis<d,real_t> >(m_bases_u[k].basis(0).clone()).release() );

    //Check if grid should be coarsen or not
    bool coarsen = true;
    for ( index_t k = 0; k< m_tarDim; ++k)
        coarsen = coarsen && (bases_u[k][0]->size() > coarseLvlMG);

    //Coarsen bases until they are less then coarseLvlMG
    index_t lvl = 0;
    while(coarsen)
    {
        for ( index_t k = 0; k< m_tarDim; ++k)
        {
            bases_u[k].push_back( bases_u[k][lvl]->clone().release() );
            bases_u[k][lvl+1]->uniformCoarsen();
        }

        ++lvl;
        //Check if grid should be coarsen again or not
        for ( index_t k = 0; k< m_tarDim; ++k)
            coarsen = coarsen && (bases_u[k][lvl]->size() > coarseLvlMG);
    }

    // Reorder the bases from coarsest to finest
    for ( index_t k = 0; k< m_tarDim; ++k)
        std::reverse( bases_u[k].begin(), bases_u[k].end() );

    const int numTransfer = bases_u[0].size() - 1;
    gsInfo << "The number of matrix transfers is: " << numTransfer <<std::endl;

    //Create the transfer matrix for each level
    std::vector<gsTensorBSplineBasis<d,real_t>* > tmpBases(m_tarDim);
    std::vector<gsSparseMatrix<real_t, RowMajor> >transMat_u(m_tarDim);
    std::vector<index_t> nRows_u(m_tarDim);
    std::vector<index_t> nCols_u(m_tarDim);

    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices   ( numTransfer );
    for (int i = 0; i < numTransfer; ++i)
    {
        for ( index_t k = 0; k< m_tarDim; ++k)
        {
            tmpBases[k] = bases_u[k][i]->clone().release();
            tmpBases[k]->refine_withTransfer( transMat_u[k], removedKnots_u[k][numTransfer - 1 - i] );
            nRows_u[k] = transMat_u[k].rows();
            nCols_u[k] = transMat_u[k].cols();
        }

        freeAll(tmpBases);

        //Take the transfer matrix for each component of the velocity and
        //create a big block diagonal transfer matrix for the whole velocity field
        gsSparseEntries<real_t> entries;
        index_t totalNonZeros = 0;
        for ( index_t k = 0; k< m_tarDim; ++k)
            totalNonZeros += transMat_u[k].nonZeros();
        entries.reserve( totalNonZeros );

        index_t rowShift = 0;
        index_t colShift = 0;
        for ( index_t k = 0; k< m_tarDim; ++k)
        {
            for (int j=0; j < transMat_u[k].outerSize(); ++j)  //Not sure if this is optimal!
                for (gsSparseMatrix<real_t, Eigen::RowMajor>::InnerIterator it(transMat_u[k],j); it; ++it)
                {
                    entries.add(rowShift + it.row(), colShift + it.col(), it.value());
                }
            rowShift += nRows_u[k];
            colShift += nCols_u[k];
        }

        transferMatrices[i] = gsSparseMatrix<real_t>(rowShift,colShift);
        transferMatrices[i].setFromTriplets(entries.begin(), entries.end());
        transferMatrices[i].makeCompressed();
    }

    //Create DoF mapper for each level
    std::vector<gsDofMapper> dofMapper;

    for (int lvl2 = 0; lvl2 < numTransfer; ++lvl2)
    {
        //Make a vector fo multibasis with only velocity components.
        std::vector<const gsMultiBasis<T> *> bases_u_ptr(m_tarDim);
        for (index_t k = 0; k < m_tarDim; ++k)
            bases_u_ptr[k] = new gsMultiBasis<T>(*bases_u[k][lvl2]);

        //Create the correct DoF mapper
        if (m_geoTrans == DIV_CONFORMING && m_dirStrategy == dirichlet::eliminatNormal)
            dofMapper.push_back(*makeDofMapperForVelocityDivConforming(bases_u_ptr, m_bcondition));
        else
            dofMapper.push_back(*makeVectorValuedDofMapper(bases_u_ptr, m_bcondition, m_dirStrategy));
        freeAll(bases_u_ptr);
    }
    //Add from the finest discretization
    const gsDofMapper tmpDoFmapper(*m_dofMapper[0]);
    dofMapper.push_back(tmpDoFmapper);

    // incorporate boundary conditions
    gsDebug << "Setting up boundary conditions... " << std::flush;
    std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatricesWithBC   ( numTransfer );
    for ( index_t i=0; i<numTransfer; ++i )
        gsMultiBasis<>::combineTransferMatrices( std::vector< gsSparseMatrix<real_t, RowMajor> >(1,give(transferMatrices[i])), dofMapper[i], dofMapper[i+1], transferMatricesWithBC[i] );
    transferMatrices.clear();

    gsDebug << "done.\n";

    gsDebug << "Caling gsMultiGrid..." << std::endl;
    gsSparseMatrix<real_t> A_tmp = m_matrix.block(0,0,m_sz_A,m_sz_A);
    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( A_tmp, transferMatricesWithBC);

    for (index_t i=1; i<mg->numLevels(); ++i)
        mg->setSmoother(i,makeGaussSeidelOp(mg->matrix(i)));

    mg->setNumPreSmooth(numPreSmo);
    mg->setNumPostSmooth(numPostSmo);

    //Memory clean up!
    for (size_t k = 0; k < bases_u.size() ; ++k)
        freeAll(bases_u[k]);


    return mg;
}
template<class T>
typename gsMatrixOp<typename gsMatrix<T>::Base >::Ptr gsStokesIterativeSolver<T>::evalMassInv()
{
    //Check if inverse is calculated (if not calculate it)
    if (m_inverseMass.cols() != m_sz_B && m_inverseMass.rows() != m_sz_B)
    {
        //Get mass matrix for the pressure space
        gsGenericAssembler<T> massConst(m_patches, m_basis_p);
        const gsSparseMatrix<T> & massMatrixBtmp = massConst.assembleMass();
        massMatrixBtmp.cols();
        const gsSparseMatrix<T> & massMatrixB = massConst.fullMatrix();
        //Turn it into a dense matrix for invertion
        gsMatrix<T> denseMassB(massMatrixB);
        m_inverseMass = denseMassB.inverse();
        GISMO_ASSERT(m_inverseMass.cols() == m_sz_B && m_inverseMass.rows() == m_sz_B,
                     "Wrong dimention: Something went wrong in inverting the mass matrix");
    }
    return makeMatrixOp(m_inverseMass);
}


}// namespace gismo
