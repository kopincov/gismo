/** @file gsMultiGridAdapter.h

    @brief Provides a unified way to set up multigrid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include<gsSolver/gsMatrixOp.h>
#include<gsSolver/gsSimplePreconditioners.h>
#include<gsMultiGrid/gsGridHierarchy.h>
#include<gsMultiGrid/gsMultiGrid.h>
#include<gsMultiGrid/gsMassSmoother.h>

namespace gismo {

// Setup a multigrid object based on the available data
// Only used by IETI
template< typename T >
typename gsMultiGridOp<T>::uPtr getMultiGrid( gsMultiBasis<T> mb, typename gsSparseMatrix<T>::Ptr systemMatrix,
    const gsBoundaryConditions<T>& bc, const gsOptionList& opt, const gsMultiPatch<T>* mp = NULL )
{
    typedef typename gsSparseMatrix<T, RowMajor>::Ptr TransferMatrixPtr;

    gsGridHierarchy<T> gh = gsGridHierarchy<T>::buildByCoarsening(give(mb),bc,opt);
    std::string smoother = opt.askString("Smoother", "unknown");
    const size_t sz = gh.getTransferMatrices().size();

    std::vector<TransferMatrixPtr> transferMatrixPtrs(sz);
    for (size_t i = 0; i < sz; ++i)
        transferMatrixPtrs[i] = TransferMatrixPtr(new gsSparseMatrix<T, RowMajor>(gh.getTransferMatrices()[i]));

    typename gsMultiGridOp<T>::uPtr mg = gsMultiGridOp<T>::make( give(systemMatrix), give(transferMatrixPtrs) );
    mg->setOptions(opt);

    for (index_t i = 1; i < mg->numLevels(); ++i)
    {
        if (smoother=="Richardson")
        {
            typename gsPreconditionerOp<T>::Ptr result = makeRichardsonOp(mg->matrix(i));
            result->setOptions(opt);
            mg->setSmoother(i, result);
        }
        else if (smoother=="Jacobi")
        {
            typename gsPreconditionerOp<T>::Ptr result = makeJacobiOp(mg->matrix(i));
            result->setOptions(opt);
            mg->setSmoother(i, result);
        }
        else if (smoother=="GaussSeidel")
        {
            typename gsPreconditionerOp<T>::Ptr result = makeGaussSeidelOp(mg->matrix(i));
            result->setOptions(opt);
            mg->setSmoother(i, result);

        }
        else if (smoother=="SubspaceCorrectedMassSmoother")
        {
            if ( mb.nBases() == 1 )
            {
                typename gsPreconditionerOp<T>::Ptr result = gsPreconditionerFromOp<T>::make(
                    mg->underlyingOp(i),
                    makeSubspaceCorrectedMassSmootherOperator(
                        gh.getMultiBases()[i][0],
                        opt.askReal( "Damping", 0.2 ),
                        bc,
                        opt.askReal( "Alpha", 0 )
                    ),
                    opt.askReal( "OuterDamping", .75 )
                );
                mg->setSmoother(i, result);
            }
            else
            {
                typename gsPreconditionerOp<T>::Ptr result = gsAdditiveSmoother::make(
                        mg->underlyingOp(i),
                        setupPiecewisePreconditioner(
                            mg->matrix(i),
                            makeSubspaceCorrectedMassSmootherOperatorsDirichlet(
                                gh.getMultiBases()[i],
                                opt.askReal( "Damping", 0.2 ),
                                opt.askReal( "Alpha", 0 )
                            ),
                            gh.getMultiBases()[i],
                            bc,
                            opt
                        ),
                        opt.askReal( "OuterDamping", .75 )
                    );
                mg->setSmoother(i, result);
            }
        }
        else
        {
            GISMO_ERROR( "The smoother \"" + smoother + "\" is not known.\n" );
        }

    }
    return mg;

}

// Setup a multigrid object based on the available data
// Only used by IETI
template< typename T >
typename gsMultiGridOp<T>::uPtr getMultiGrid( gsMultiBasis<T> mb, gsSparseMatrix<T> systemMatrix,
    const gsBoundaryConditions<T>& bc, const gsOptionList& opt, const gsMultiPatch<T>* mp = NULL )
{ return getMultiGrid<T>( give(mb), systemMatrix.moveToPtr(), bc, opt, mp ); }


} // namespace gismo
